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
BIGWIGs = expand("{assayID}/{file}.bw",
                 assayID="ChIP-Seq",
                 file=[i + "/" + config["processed_dir"] + "/" + REF_VERSION + "/deepTools/bamCoverage/normal/RPKM/duplicates_removed/"  + j + ".Q" + config["alignment_quality"] \
                    for i in config["samples"]["ChIP-Seq"]["runID"] \
                        for j in config["samples"]["ChIP-Seq"][i]])

MergedBIGWIGs = expand("{assayID}/{file}.bw",
                        assayID="ChIP-Seq",
                        file=["merged" +
                                "/" +
                                config["processed_dir"] +
                                "/" +
                                REF_VERSION +
                                "/deepTools/bamCoverage/normal/RPKM/duplicates_removed/"  +
                                j +
                                ".Q" +
                                config["alignment_quality"] \
                                for j in ['MCF10A_WT_HP1b_ChIP',
                                          'MCF10A_shHP1a_HP1b_ChIP',
                                          'MCF10A_shHP1b_Input',
                                          'MCF10A_WT_Input',
                                          'MCF10A_shH2AZ_Input',
                                          'MCF10A_shHP1b_HP1a_ChIP',
                                          'MCF10A_shH2AZ_HP1b_ChIP',
                                          'MCF10A_shH2AZ_HP1a_ChIP',
                                          'MCF10A_shHP1a_Input',
                                          'MCF10A_WT_HP1a_ChIP',
                                          'MCF10A_WT_H2AZ_ChIP']])


# input functions
def getAllFASTQ(wildcards):
    fn = []
    for i in config["samples"][wildcards["assayID"]][wildcards["runID"]]:
        for j in config["samples"][wildcards["assayID"]][wildcards["runID"]][i]:
            fn.append("/".join([wildcards["assayID"], wildcards["runID"], config["raw_dir"], j]))
    return(fn)


def getBAMbyCondition(wildcards):
    fn = []
    for i in config["samples"]["ChIP-Seq"]["runID"]:
        for j in config["samples"]["ChIP-Seq"]["conditions"][wildcards["condition"]][i]:
            fn.append(join("ChIP-Seq/" + i + "/" + config["processed_dir"] + "/" + REF_VERSION + "/bowtie2/" + wildcards["duplicates"] + "/" + j + ".Q" + config["alignment_quality"] + ".sorted.bam"))
    return(fn)


def getAllBAMs(wildcards):
    fn = []
    for i in config["samples"]["ChIP-Seq"]["runID"]:
        for j in config["samples"]["ChIP-Seq"][i]:
            fn.append("/".join(["ChIP-Seq",
                                i,
                                config["processed_dir"],
                                REF_VERSION,
                                "bowtie2",
                                wildcards["duplicates"],
                                j + ".Q" + config["alignment_quality"] + ".sorted.bam"]))
    return(fn)


def get_sample_labels(wildcards):
    sl = []
    runIDs = config["samples"][wildcards["assayID"]]["runID"]
    for i in runIDs:
        for k in config["samples"][wildcards["assayID"]][i].keys():
            sl.append(k)
    return(sl)


def getSampleLabelsByCondition(wildcards):
    sl = []
    for i in config["samples"][wildcards["assayID"]]["conditions"][wildcards["condition"]]:
        for j in config["samples"][wildcards["assayID"]]["conditions"][wildcards["condition"]][i]:
            sl.append(j)
    return(sl)


def cli_parameters_computeMatrix(wildcards):
    if wildcards["command"] == "reference-point":
        a = config["program_parameters"][wildcards["application"]][wildcards["tool"]][wildcards["command"]]
        a["--referencePoint"]=wildcards["referencePoint"]
        return(a)
    if wildcards["command"] == "scale-regions":
        a = config["program_parameters"][wildcards["application"]][wildcards["tool"]][wildcards["command"]][wildcards["region"]]
        return(a)


def cli_parameters_normalization(wildcards):
    if wildcards["norm"] == "RPKM":
        a = "--normalizeUsingRPKM"
    elif wildcards["norm"] == "1xcoverage":
        a = " ".join(("--normalizeTo1x", config["references"][REF_GENOME]["effectiveSize"]))
    return(a)


def cli_parameters_bamCoverage(wildcards):
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


def getComputeMatrixInput(wildcards):
    fn = []
    for i in config["samples"]["ChIP-Seq"]["runID"]:
        for j in config["samples"]["ChIP-Seq"][i]:
            fn.append("/".join(["ChIP-Seq",
                                i,
                                config["processed_dir"],
                                REF_VERSION,
                                "deepTools",
                                "bamCoverage",
                                wildcards["mode"],
                                wildcards["norm"],
                                wildcards["duplicates"],
                                j + ".Q" + config["alignment_quality"] + ".bw"]))
    return(fn)


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
        "{assayID}/{outdir}/{reference_version}/merged/{duplicates}/{replicates}.bam.bai"
    shell:
        """
            samtools index {input} {output} -@ {threads} 1>>{log} 2>>{log}
        """

<<<<<<< HEAD
rule bamCoverage:
    version:
        0.2
    params:
        deepTools_dir = home + config["program_parameters"]["deepTools"]["deepTools_dir"],
        ignore = config["program_parameters"]["deepTools"]["ignoreForNormalization"],
        program_parameters = cli_parameters_bamCoverage
    threads:
        lambda wildcards: int(str(config["program_parameters"]["deepTools"]["threads"]).strip("['']"))
    input:
        bam = lambda wildcards: wildcards["assayID"] + "/" + wildcards["runID"] + "/" + wildcards["outdir"] + "/" + wildcards["reference_version"] + "/bowtie2/" + wildcards["duplicates"] + "/" +  wildcards["sample"] + ".sorted.bam"
    output:
        "{assayID}/{runID}/{outdir}/{reference_version}/{application}/{tool}/{mode}/{norm}/{duplicates}/{sample}.bw"
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
=======
# rule bamCoverage:
#     version:
#         0.1
#     params:
#         deepTools_dir = home + config["program_parameters"]["deepTools"]["deepTools_dir"],
#         ignore = config["program_parameters"]["deepTools"]["ignoreForNormalization"],
#         program_parameters = cli_parameters_bamCoverage
#     threads:
#         lambda wildcards: int(str(config["program_parameters"]["deepTools"]["threads"]).strip("['']"))
#     input:
#         bam = lambda wildcards: wildcards["assayID"] + "/" + wildcards["runID"] + "/" + wildcards["outdir"] + "/" + wildcards["reference_version"] + "/bowtie2/" + wildcards["duplicates"] + "/" +  wildcards["sample"] + ".sorted.bam"
#     output:
#         "{assayID}/{runID}/{outdir}/{reference_version}/{application}/{tool}/{mode}/{norm}/{duplicates}/{sample}.bw"
#     shell:
#         """
#             {params.deepTools_dir}/bamCoverage --bam {input.bam} \
#                                                --outFileName {output} \
#                                                --outFileFormat bigwig \
#                                                {params.program_parameters} \
#                                                --numberOfProcessors {threads} \
#                                                --normalizeUsingRPKM \
#                                                --ignoreForNormalization {params.ignore}
#         """
>>>>>>> 6bb07b2d823916ec3957e6d0ddfed2f3ba370808
# target rules
rule all:
    input:
        MergedBIGWIGs,
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
