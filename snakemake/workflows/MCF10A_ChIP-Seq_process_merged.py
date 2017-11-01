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
                 file=[i + "/" + config["processed_dir"] + "/" + REF_VERSION + "/deepTools/bamCoverage/normal/RPKM/duplicates_removed/"  + j + ".Q" + config["alignment_quality"]\
                    for i in config["samples"]["ChIP-Seq"]["runID"] \
                        for j in config["samples"]["ChIP-Seq"][i]])

MergedBIGWIGs = expand("{assayID}/{file}.bw",
                         assayID="ChIP-Seq",
                         file=["merged/" + config["processed_dir"] + "/" + REF_VERSION + "/deepTools/bamCoverage/normal/RPKM/duplicates_removed/"  + j \
                                for j in config["samples"]["ChIP-Seq"]["replicates"].keys()])


# input functions
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
        if wildcards["region"] in ["allGenes", "intergenicRegions"]:
            a = config["program_parameters"][wildcards["application"]][wildcards["tool"]][wildcards["command"]][wildcards["region"]]
            return(a)
        else:
            a = config["program_parameters"][wildcards["application"]][wildcards["tool"]][wildcards["command"]]["default"]
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


def getComputeMatrixInputMerged(wildcards):
    fn = []
    for i in config["samples"]["ChIP-Seq"]["replicates"].keys():
        fn.append("/".join(["ChIP-Seq",
                            "merged",
                            config["processed_dir"],
                            REF_VERSION,
                            "deepTools",
                            "bamCoverage",
                            wildcards["mode"],
                            wildcards["norm"],
                            wildcards["duplicates"],
                            i + ".bw"]))
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
    return(fn)

# rules section
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

rule bamCoverageMerged:
    version:
        0.1
    params:
        deepTools_dir = home + config["program_parameters"]["deepTools"]["deepTools_dir"],
        ignore = config["program_parameters"]["deepTools"]["ignoreForNormalization"],
        program_parameters = cli_parameters_bamCoverage
    threads:
        lambda wildcards: int(str(config["program_parameters"]["deepTools"]["threads"]).strip("['']"))
    input:
        bam = lambda wildcards: wildcards["assayID"] + "/merged/" + wildcards["outdir"] + "/" + wildcards["reference_version"] + "/" + wildcards["duplicates"] + "/" +  wildcards["replicates"] + ".bam"
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

rule computeMatrix:
    version:
        0.3
    params:
        deepTools_dir = home + config["program_parameters"]["deepTools"]["deepTools_dir"],
        program_parameters = lambda wildcards: ' '.join("{!s}={!s}".format(key, val.strip("\\'")) for (key, val) in cli_parameters_computeMatrix(wildcards).items())
    threads:
        32
    input:
        file = getComputeMatrixInputMerged,
        region = lambda wildcards: home + config["program_parameters"]["deepTools"]["regionFiles"][wildcards["reference_version"]][wildcards["region"]]
    output:
        matrix_gz = "{assayID}/{runID}/{outdir}/{reference_version}/{application}/{tool}/{command}/{duplicates}/{referencePoint}/{region}_{mode}_{norm}.matrix.gz"
    shell:
        """
            {params.deepTools_dir}/computeMatrix {wildcards.command} \
                                                 --regionsFileName {input.region} \
                                                 --scoreFileName {input.file} \
                                                 --missingDataAsZero \
                                                 --skipZeros \
                                                 --numberOfProcessors {threads} \
                                                 {params.program_parameters} \
                                                 --outFileName {output.matrix_gz}
        """

rule plotProfileMerged:
    version:
        0.1
    params:
        deepTools_dir = home + config["program_parameters"]["deepTools"]["deepTools_dir"],
    input:
        matrix_gz = "{assayID}/{runID}/{outdir}/{reference_version}/{application}/computeMatrix/{command}/{duplicates}/{referencePoint}/{region}_{mode}_{norm}.matrix.gz"
    output:
        figure = "{assayID}/{runID}/{outdir}/{reference_version}/{application}/{tool}/{command}/{duplicates}/{referencePoint}/{plotType}.{mode}.{norm}.{region}.pdf",
        data = "{assayID}/{runID}/{outdir}/{reference_version}/{application}/{tool}/{command}/{duplicates}/{referencePoint}/{plotType}.{mode}.{norm}.{region}.data",
        regions = "{assayID}/{runID}/{outdir}/{reference_version}/{application}/{tool}/{command}/{duplicates}/{referencePoint}/{plotType}.{mode}.{norm}.{region}.bed"
    shell:
        """
            {params.deepTools_dir}/plotProfile --matrixFile {input.matrix_gz} \
                                               --outFileName {output.figure} \
                                               --outFileNameData {output.data} \
                                               --outFileSortedRegions {output.regions} \
                                               --plotType {wildcards.plotType}
        """

rule bigwigCompareMerged:
    version:
        0.1
    params:
        deepTools_dir = home + config["program_parameters"]["deepTools"]["deepTools_dir"],
    input:
        chip = lambda wildcards: "/".join(wildcards["assayID"],
                                           "merged",
                                           wildcards["outdir"],
                                           wildcards["reference_version"],
                                           wildcards["application"],
                                           "bamCoverage",
                                           wildcards["mode"],
                                           wildcards["norm"],
                                           wildcards["duplicates"],
                                           config["ChIP-Seq"]["ChIP-Input"][wildcards["condition"]]["ChIP"],
                                           ".bw"),
        input = lambda wildcards: "/".join(wildcards["assayID"],
                                           "merged",
                                           wildcards["outdir"],
                                           wildcards["reference_version"],
                                           wildcards["application"],
                                           "bamCoverage",
                                           wildcards["mode"],
                                           wildcards["norm"],
                                           wildcards["duplicates"],
                                           config["ChIP-Seq"]["ChIP-Input"][wildcards["condition"]]["Input"],
                                           ".bw")
    output:
        "{assayID}/merged/{outdir}/{reference_version}/{application}/{tool}/{mode}/{norm}/{duplicates}/{condition}_{contrast}_{ratio}.bw"
    shell:
        """
            {params.deepTools_dir}/bigwigCompare --bigwig1 {input.chip}\
                                                 --bigwig2 {input.input}\
                                                 --ratio {wildcards.ratio}\
                                                 --outFileFormat bigwig\
                                                 --outFileName {output}
        """

# target rules
rule all:
    input:
        expand("{assayID}/{runID}/{outdir}/{reference_version}/{application}/{tool}/{mode}/{norm}/{duplicates}/{condition}_{contrast}_{ratio}.bw",
               assayID="ChIP-Seq",
               runID="merged",
               outdir=config["processed_dir"],
               reference_version=REF_VERSION,
               application=["deepTools"],
               tool=["plotProfile"],
               mode=["normal"],
               norm=["RPKM"],
               duplicates=["duplicates_removed"],
               condition=["MCF10A_WT_HP1a",
                          "MCF10A_WT_HP1b",
                          "MCF10A_WT_H2AZ",
                          "MCF10A_shHP1a_HP1b",
                          "MCF10A_shHP1b_HP1a",
                          "MCF10A_shH2AZ_HP1a",
                          "MCF10A_shH2AZ_HP1b"],
               contrast=["ChIP-Input"],
               ratio="log2"),
        expand("{assayID}/{runID}/{outdir}/{reference_version}/{application}/{tool}/{command}/{duplicates}/{referencePoint}/{plotType}.{mode}.{norm}.{region}.{suffix}",
               assayID="ChIP-Seq",
               runID="merged",
               outdir=config["processed_dir"],
               reference_version=REF_VERSION,
               application="deepTools",
               tool="plotProfile",
               command=["scale-regions"],
               duplicates=["duplicates_removed"],
               referencePoint="TSS",
               plotType="se",
               mode=["normal"],
               norm=["RPKM"],
               region=["allGenes", "intergenicRegions", "conditionMCF10A_shH2AZHP1a", "conditionMCF10A_shHP1ab", "conditionMCF10A_shHP1a", "conditionMCF10A_shHP1b", "conditionMCF10A_WT"],
               suffix=["pdf", "data", "bed"]),
        expand("{assayID}/merged/{outdir}/{reference_version}/duplicates_removed/{replicates}.bam.bai",
               assayID="ChIP-Seq",
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
                           'MCF10A_WT_H2AZ_ChIP']),
        MergedBIGWIGs
