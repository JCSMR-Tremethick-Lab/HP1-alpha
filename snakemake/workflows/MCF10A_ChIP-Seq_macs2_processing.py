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


# set targets here

PROCESSED_BAMs_pseudo_reps = expand("{assayID}/{file1}",
                                     assayID = "ChIP-Seq",
                                     file1 = [ i + "/" + config["processed_dir"] + "/" + REF_VERSION + "/bowtie2/duplicates_removed/" + j + ".Q" + config["alignment_quality"] + "/" + "pseudo_rep" + rep + ".sorted.bam.bai" \
                                        for i in config["samples"]["ChIP-Seq"]["runID"]\
                                            for j in config["samples"]["ChIP-Seq"][i]\
                                                for rep in ["1", "2"]])



rule make_pseudo_replicates:
    params:
        seed = print(randint(0,1000)),
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
            samtools view -b --threads {threads} -s {params.seed}.{params.fraction} {input} > {output} 2 > {log}
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

rule all:
    input:
        PROCESSED_BAMs_pseudo_reps
