_author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2017-08-23"

from snakemake.exceptions import MissingInputException
import os
from os.path import join
from random import *

rule:
    version: 0.1

localrules:
    all

# variables
home_dir = os.environ['HOME']
home = home_dir

wrapper_dir = os.environ['HOME'] + "/Development/snakemake-wrappers/bio"

REF_GENOME = "hg19"
REF_VERSION = config["references"][REF_GENOME]["version"][1]

# includes
include_prefix= os.environ['HOME'] + "/Development/JCSMR-Tremethick-Lab/Breast/snakemake/rules/"

# input functions
def getAllBAMs(duplicates):
    fn = []
    for i in config["samples"]["ChIP-Seq"]["runID"]:
        for j in config["samples"]["ChIP-Seq"][i]:
            fn.append(join("ChIP-Seq/" + i + "/" + config["processed_dir"] + "/" + REF_VERSION + "/bowtie2/" + duplicates + "/" + j + ".Q" + config["alignment_quality"] + ".sorted.bam"))
    return(fn)

def get_sample_labels(wildcards):
    sl = []
    runIDs = config["samples"][wildcards["assayID"]]["runID"]
    for i in runIDs:
        for k in config["samples"][wildcards["assayID"]][i].keys():
            sl.append(k)
    return(sl)


def macs2OutputFiles():
    fn = []
    for i in config["samples"]["ChIP-Seq"]["ChIP-Input"].keys():
        for j in config["samples"]["ChIP-Seq"]["replicates"][config["samples"]["ChIP-Seq"]["ChIP-Input"][i]["ChIP"]]:
            contrasts = i + "/" + j
            path = "/".join(["ChIP-Seq",
                            "processed_data",
                            REF_VERSION,
                            "macs2",
                            contrasts,
                            "callpeak"])
            fn.append(path)
    return(fn)


def macs2OutputFilesPseudoReps():
    fn = []
    for i in config["samples"]["ChIP-Seq"]["ChIP-Input"].keys():
        for j in config["samples"]["ChIP-Seq"]["replicates"][config["samples"]["ChIP-Seq"]["ChIP-Input"][i]["ChIP"]]:
            pseudorep1 = i + "/" + j + "/pseudo_rep1"
            pseudorep2 = i + "/" + j + "/pseudo_rep2"
            pseudorep1 = "/".join(["ChIP-Seq",
                            "processed_data",
                            REF_VERSION,
                            "macs2",
                            pseudorep1,
                            "callpeak"])
            pseudorep2 = "/".join(["ChIP-Seq",
                            "processed_data",
                            REF_VERSION,
                            "macs2",
                            pseudorep2,
                            "callpeak"])
            fn.append(pseudorep1)
            fn.append(pseudorep2)
    return(fn)


def getMergedInputBAM(wildcards):
    fn = []
    contrast = wildcards["contrast"].split("/")[0]
    inp = config["samples"]["ChIP-Seq"]["ChIP-Input"][contrast]["Input"]
    inp = config["samples"]["ChIP-Seq"]["replicates"][inp]
    for i in config["samples"]["ChIP-Seq"]["runID"]:
        for j in config["samples"]["ChIP-Seq"][i].keys():
            for k in inp:
                if k == j:
                    f = "/".join([wildcards["assayID"],
                                  i,
                                  config["processed_dir"],
                                  REF_VERSION,
                                  "bowtie2",
                                  "duplicates_removed",
                                  j + ".Q" + config["alignment_quality"] + ".sorted.bam"
                                  ])
                    fn.append(f)
    return(fn)


def getChIPBam(wildcards):
    fn = []
    for i in config["samples"]["ChIP-Seq"]["runID"]:
        for j in config["samples"]["ChIP-Seq"][i].keys():
                if j == wildcards["unit"]:
                        f = "/".join([wildcards["assayID"],
                                      i,
                                      config["processed_dir"],
                                      REF_VERSION,
                                      "bowtie2",
                                      "duplicates_removed",
                                      j + ".Q" + config["alignment_quality"] + ".sorted.bam"
                                      ])
                        fn.append(f)
    return(fn)


def getChIPBamPseudoRep(wildcards):
    fn = []
    for i in config["samples"]["ChIP-Seq"]["runID"]:
        for j in config["samples"]["ChIP-Seq"][i].keys():
                if j == wildcards["unit"]:
                        f = "/".join([wildcards["assayID"],
                                      i,
                                      wildcards["processed_dir"],
                                      REF_VERSION,
                                      "bowtie2",
                                      "duplicates_removed",
                                      j + ".Q" + config["alignment_quality"] + "/" + wildcards["pseudo"] + ".sorted.bam"
                                      ])
                        fn.append(f)
    return(fn)

# set targets here

PROCESSED_BAMs_pseudo_reps = expand("{assayID}/{file1}",
                                     assayID = "ChIP-Seq",
                                     file1 = [ i + "/" + config["processed_dir"] + "/" + REF_VERSION + "/bowtie2/duplicates_removed/" + j + ".Q" + config["alignment_quality"] + "/" + "pseudo_rep" + rep + ".sorted.bam.bai" \
                                        for i in config["samples"]["ChIP-Seq"]["runID"]\
                                            for j in config["samples"]["ChIP-Seq"][i]\
                                                for rep in ["1", "2"]])

MACS2_output = macs2OutputFiles()
MACS2_output_pseudo_reps = macs2OutputFilesPseudoReps()

rule make_pseudo_replicates:
    params:
        seed = print(randint(0, 1000)),
        fraction = 0.5
    threads:
        4
    log:
        "{assayID}/{runID}/{outdir}/{reference_version}/bowtie2/duplicates_removed/{unit}.Q{qual}/log.txt"
    input:
        "{assayID}/{runID}/{outdir}/{reference_version}/bowtie2/duplicates_removed/{unit}.Q{qual}.sorted.bam"
    output:
        protected("{assayID}/{runID}/{outdir}/{reference_version}/bowtie2/duplicates_removed/{unit}.Q{qual}/pseudo_{rep}.sorted.bam")
    shell:
        """
            samtools view -b --threads {threads} -s {params.seed}.{params.fraction} -o {output} {input} 1>{log} 2>{log}
        """

rule index_pseudo_replicates:
    threads:
        4
    log:
    input:
        "{assayID}/{runID}/{outdir}/{reference_version}/bowtie2/duplicates_removed/{unit}.Q{qual}/pseudo_{rep}.sorted.bam"
    output:
        "{assayID}/{runID}/{outdir}/{reference_version}/bowtie2/duplicates_removed/{unit}.Q{qual}/pseudo_{rep}.sorted.bam.bai"
    shell:
        """
            samtools index -@ {threads} {input} {output}
        """

rule macs2_callpeak_replicates:
    version:
        0.2
    threads:
        1
    params:
        gsize=config["program_parameters"]["macs2"]["gsize"],
        filetype="BAM",
        verbosity=config["program_parameters"]["macs2"]["verbosity"],
	    macs2_binary=home + config["program_parameters"]["macs2"]["binary"]
    log:
        "{assayID}/{outdir}/{reference_version}/macs2/{contrast}/{unit}/callpeak/callpeak.log"
    input:
        chip=getChIPBam,
        input=getMergedInputBAM
    output:
        "{assayID}/{outdir}/{reference_version}/macs2/{contrast}/{unit}/callpeak"
    shell:
        """
            {params.macs2_binary} callpeak -t {input.chip}\
                           -c {input.input}\
                           --gsize {params.gsize}\
                           -f {params.filetype}\
                           --name {wildcards.contrast}\
                           --outdir {output}\
                           --verbose {params.verbosity}\
                           --bdg\
            1>>{log} 2>>{log}
        """

rule macs2_callpeak_pseudoreplicates:
    version:
        0.1
    threads:
        1
    params:
        gsize=config["program_parameters"]["macs2"]["gsize"],
        filetype="BAM",
        verbosity=config["program_parameters"]["macs2"]["verbosity"],
        name=lambda wildcards: ".".join([wildcards["unit"], wildcards["pseudo"]]),
        macs2_binary=home + config["program_parameters"]["macs2"]["binary"]
    log:
        "{assayID}/{outdir}/{reference_version}/macs2/{contrast}/{unit}/callpeak/{pseudo}/callpeak.log"
    input:
        chip=getChIPBamPseudoRep,
        input=getMergedInputBAM
    output:
        "{assayID}/{outdir}/{reference_version}/macs2/{contrast}/{unit}/callpeak/{pseudo}"
    shell:
        """
            {params.macs2_binary} {wildcard.macs2_command} -t {input.chip}\
                                           -c {input.input}\
                                           --gsize {params.gsize}\
                                           -f {params.filetype}\
                                           --name {params.name}\
                                           --outdir {output}\
                                           --verbose {params.verbosity}\
                                           --bdg\
            1>>{log} 2>>{log}
        """


rule all:
    input:
        PROCESSED_BAMs_pseudo_reps,
        MACS2_output,
	    MACS2_output_pseudo_reps
