_author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2017-08-23"

from snakemake.exceptions import MissingInputException
import os
from os.path import join

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

PROCESSED_BAMs_dups_removed = expand("{assayID}/{file1}",
                                     assayID = "ChIP-Seq",
                                     file1 = [ i + "/" + config["processed_dir"] + "/" + REF_VERSION + "/bowtie2/duplicates_removed/" + j + ".Q" + config["alignment_quality"] + ".sorted.bam.bai" \
                                        for i in config["samples"]["ChIP-Seq"]["runID"] \
                                            for j in config["samples"]["ChIP-Seq"][i]]),

PROCESSED_BAMs_dups_marked = expand("{assayID}/{file2}",
                                     assayID = "ChIP-Seq",
                                     file2 = [ i + "/" + config["processed_dir"] + "/" + REF_VERSION + "/bowtie2/duplicates_marked/" + j + ".Q" + config["alignment_quality"] + ".sorted.bam.bai" \
                                     for i in config["samples"]["ChIP-Seq"]["runID"] \
                                        for j in config["samples"]["ChIP-Seq"][i]])



rule bam_rmdup:
    input:
        rules.bam_mark_duplicates.output
    output:
        protected("{assayID}/{runID}/{outdir}/{reference_version}/bowtie2/duplicates_removed/{unit}.Q{qual}.sorted.bam")
    shell:
        "samtools rmdup {input} {output}"

rule bam_rmdup_index:
    params:
        qual = config["alignment_quality"]
    input:
        rules.bam_rmdup.output
    output:
        protected("{assayID}/{runID}/{outdir}/{reference_version}/bowtie2/duplicates_removed/{unit}.Q{qual}.sorted.bam.bai")
    shell:
        "samtools index {input} {output}"

rule all:
    input:
        PROCESSED_BAMs_dups_removed,
        PROCESSED_BAMs_dups_marked
