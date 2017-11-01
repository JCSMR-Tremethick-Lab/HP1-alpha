from snakemake.exceptions import MissingInputException
import os
from os.path import join

_author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2017-11-01"

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

# build targets
MergedBIGWIGs = expand("{assayID}/{file}.bw",
                         assayID="ChIP-Seq",
                         file=["merged/" + config["processed_dir"] + "/" + REF_VERSION + "/deepTools/bamCoverage/normal/RPKM/duplicates_removed/"  + j \
                                for j in config["samples"]["ChIP-Seq"]["replicates"].keys()])


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
    return(fn)


def cliParametersBamCoverage(wildcards):
    a = config["program_parameters"][wildcards["application"]][wildcards["tool"]][wildcards["mode"]]
    b = str()
    for (key, val) in a.items():
        if val == " ":
            f = key + " "
            b = b + f
        else:
            f = key + "=" + val + " "
            b = b + f
    if wildcards["mode"] == "MNase":
        b = b + "--MNase"
    return(b.rstrip())

# rules
rule bamMerge:
    version:
        0.1
    params:
        outputFormat = "--output-fmt BAM"  # ToDo: move
    log:
        "logs/{replicates}.bamMerge.log"
    threads:
        16
    input:
        getBAMbyReplicates
    output:
        "{assayID}/merged/{outdir}/{reference_version}/{duplicates}/{replicates}.bam"
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
        "{assayID}/merged/{outdir}/{reference_version}/{duplicates}/{replicates}.bam.bai"
    shell:
        """
            samtools index {input} {output} -@ {threads} 1>>{log} 2>>{log}
        """

# rules section
rule bamCoverageMerged:
    version:
        0.1
    params:
        deepTools_dir = home + config["program_parameters"]["deepTools"]["deepTools_dir"],
        ignore = config["program_parameters"]["deepTools"]["ignoreForNormalization"],
        program_parameters = cliParametersBamCoverage
    threads:
        lambda wildcards: int(str(config["program_parameters"]["deepTools"]["threads"]).strip("['']"))
    input:
        bai = mergeBams("{assayID}/merged/{outdir}/{reference_version}/{duplicates}/{replicates}.bam.bai"),
        bam = mergeBams("{assayID}/merged/{outdir}/{reference_version}/{duplicates}/{replicates}.bam")
    output:
        "{assayID}/merged/{outdir}/{reference_version}/{application}/{tool}/{mode}/{norm}/{duplicates}/{replicates}.bw"
    shell:
        """
            {params.deepTools_dir}/bamCoverage --bam {input.bam} \
                                               --outFileName {output} \
                                               --outFileFormat bigwig \
                                               {params.program_parameters} \
                                               --numberOfProcessors {threads} \
                                               --normalizeUsingRPKM \
                                               --ignoreForNormalization {params.ignore}
        """


# build targets
rule all:
    input:
        MergedBIGWIGs
