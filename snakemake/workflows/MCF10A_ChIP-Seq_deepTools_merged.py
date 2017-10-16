from snakemake.exceptions import MissingInputException
import os
from os.path import join

_author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2017-10-16"

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
include_prefix = os.environ['HOME'] + "/Development/JCSMR-Tremethick-Lab/HP1-alpha/snakemake/rules/"

# target files

# functions
def getBAMbyReplicates(wildcards):
    fn=[]
    for i in config["samples"]["ChIP-Seq"]["runID"]:
        for j in config["samples"]["ChIP-Seq"][i]:
            if j in config["samples"]["ChIP-Seq"]["replicates"][wildcards["replicates"]]:
                fn.append(join("ChIP-Seq/" +
                               i +
                               "/" +
                               config["processed_dir"] +
                               "/" +
                               REF_VERSION +
                               "/bowtie2/" +
                               wildcards["duplicates"] +
                               "/" +
                               j +
                               ".Q" +
                               config["alignment_quality"] +
                               ".sorted.bam"))
                print(i, j)
    return(fn)

# rules
rule bamMerge:
    version:
        0.1
    params:
        outputFormat="--output-fmt BAM"  # ToDo: move
    log:
        "logs/{replicates}.bamMerge.log"
    threads:
        16
    input:
        getBAMbyReplicates
    output:
        "{assayID}/{runID}/{outdir}/{reference_version}/bowtie2/{duplicates}/{replicates}.Q20.sorted.bam"
    shell:
        """
            samtools merge -f {output} {input} --threads {threads} {params.outputFormat} 1>>{log} 2>>{log}
        """

rule indexMerged:
    version:
        0.1
    log:
        "logs/{replicates}.indexMerged.log"
    threads:
        16
    input:
        rules.bamMerge.output
    output:
        "{assayID}/{runID}/{outdir}/{reference_version}/bowtie2/{duplicates}/{replicates}.Q20.sorted.bam.bai"
    shell:
        """
            samtools index {input} {output} -@ {threads} 1>>{log} 2>>{log}
        """

# target rules
rule all:
    input:
        expand("{assayID}/{runID}/{outdir}/{reference_version}/bowtie2/duplicates_removed/{replicates}.Q20.sorted.bam.bai",
               assayID="ChIP-Seq",
               runID="merged",
               outdir=config["processed_dir"],
               reference_version=REF_VERSION,
               replicates=['MCF10A_WT_HP1b_ChIP',
                           'MCF10A_shHP1a_HP1b_ChIP',
                           'MCF10A_shHP1b_Input',
                           'MCF10A_WT_Input',
                           'MCF10A_shH2AZ_Input',
                           'MCF10A_shHP1b_HP1a_ChIP',
                           'MCF10A_shH2AZ_HP1b_ChIP',
                           'MCF10A_shH2AZ_HP1a_ChIP',
                           'MCF10A_shHP1a_Input',
                           'MCF10A_WT_HP1a_ChIP',
                           'MCF10A_WT_H2AZ_ChIP'])
