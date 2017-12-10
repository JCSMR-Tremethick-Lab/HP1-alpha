# from snakemake.exceptions import MissingInputException
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

# functions
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

def getComputeMatrixInputMerged(wildcards):
    fn = []
    path = "/".join([wildcards["assayID"],
                     wildcards["runID"],
                     config["processed_dir"],
                     REF_VERSION,
                     wildcards["application"]]) + "/"
    if wildcards["source"] == "bamCoverage":
        for i in config["samples"]["ChIP-Seq"]["replicates"].keys():
            fn.append(path + "/".join([wildcards["source"],
                                wildcards["mode"],
                                wildcards["norm"],
                                wildcards["duplicates"],
                                i + ".bw"]))
    elif wildcards["source"] == "bigwigCompare":
        for i in config["samples"]["ChIP-Seq"][wildcards["contrast"]].keys():
            fn.append(path + "/".join([wildcards["source"],
                                wildcards["mode"],
                                wildcards["norm"],
                                wildcards["duplicates"],
                                wildcards["contrast"],
                                wildcards["ratio"],
                                i + ".bw"]))
    return(fn)


def getRegionFiles(wildcards):
    fn = []
    for i in config["program_parameters"]["deepTools"]["regionFiles"][wildcards["reference_version"]][wildcards["region"]].values():
        fn.append(home + i)
    return(fn)


def debugWildcards(wildcards):
    print(wildcards)

# subworkflows section
subworkflow mergeBams:
    workdir:  home + "/Data/Tremethick/HP1-alpha"
    snakefile: home + "/Development/JCSMR-Tremethick-Lab/HP1-alpha/snakemake/workflows/subworkflows/mergeBam.py"


rule bigwigCompareMerged:
    version:
        0.1
    params:
        deepTools_dir = home + config["program_parameters"]["deepTools"]["deepTools_dir"],
        debug = debugWildcards
    threads:
        32
    input:
        chip = lambda wildcards: "/".join((wildcards["assayID"],
                                           "merged",
                                           wildcards["outdir"],
                                           wildcards["reference_version"],
                                           wildcards["application"],
                                           "bamCoverage",
                                           wildcards["mode"],
                                           wildcards["norm"],
                                           wildcards["duplicates"],
                                           ".".join((config["samples"]["ChIP-Seq"]["ChIP-Input"][wildcards["condition"]]["ChIP"], "bw")))),
        input = lambda wildcards: "/".join((wildcards["assayID"],
                                           "merged",
                                           wildcards["outdir"],
                                           wildcards["reference_version"],
                                           wildcards["application"],
                                           "bamCoverage",
                                           wildcards["mode"],
                                           wildcards["norm"],
                                           wildcards["duplicates"],
                                           ".".join((config["samples"]["ChIP-Seq"]["ChIP-Input"][wildcards["condition"]]["Input"], "bw"))))
    output:
        "{assayID}/{runID}/{outdir}/{reference_version}/{application}/{tool}/{mode}/{norm}/{duplicates}/{contrast}/{ratio}/{condition}.bw"
    shell:
        """
            {params.deepTools_dir}/bigwigCompare --bigwig1 {input.chip}\
                                                 --bigwig2 {input.input}\
                                                 --ratio {wildcards.ratio}\
                                                 --numberOfProcessors {threads} \
                                                 --outFileFormat bigwig\
                                                 --outFileName {output}
        """

rule computeMatrixContrasts:
    version:
        0.3
    params:
        deepTools_dir = home + config["program_parameters"]["deepTools"]["deepTools_dir"],
        program_parameters = lambda wildcards: ' '.join("{!s}={!s}".format(key, val.strip("\\'")) for (key, val) in cli_parameters_computeMatrix(wildcards).items())
    threads:
        32
    input:
        file = getComputeMatrixInputMerged,
        region = getRegionFiles
    output:
        matrix_gz = "{assayID}/{runID}/{outdir}/{reference_version}/{application}/{tool}/{command}/{duplicates}/{referencePoint}/{contrast}/{ratio}/{source}.{region}.{mode}.{norm}.matrix.gz"
    shell:
        """
            {params.deepTools_dir}/computeMatrix {wildcards.command} \
                                                 --regionsFileName {input.region} \
                                                 --scoreFileName {input.file} \
                                                 --missingDataAsZero \
                                                 --skipZeros \
                                                 --binSize 10 \
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
        matrix_gz = "{assayID}/{runID}/{outdir}/{reference_version}/{application}/{tool}/{command}/{duplicates}/{referencePoint}/{region}.{mode}.{norm}.matrix.gz"
    output:
        figure = "{assayID}/{runID}/{outdir}/{reference_version}/{application}/{tool}/{command}/{duplicates}/{referencePoint}/{plotType}.{source}.{region}.{mode}.{norm}.pdf",
        data = "{assayID}/{runID}/{outdir}/{reference_version}/{application}/{tool}/{command}/{duplicates}/{referencePoint}/{plotType}.{source}.{region}.{mode}.{norm}.data",
        regions = "{assayID}/{runID}/{outdir}/{reference_version}/{application}/{tool}/{command}/{duplicates}/{referencePoint}/{plotType}.{source}.{region}.{mode}.{norm}.bed"
    shell:
        """
            {params.deepTools_dir}/plotProfile --matrixFile {input.matrix_gz} \
                                               --outFileName {output.figure} \
                                               --outFileNameData {output.data} \
                                               --outFileSortedRegions {output.regions} \
                                               --plotType {wildcards.plotType}
        """

rule plotProfileContrasts:
    version:
        0.1
    params:
        deepTools_dir = home + config["program_parameters"]["deepTools"]["deepTools_dir"],
    input:
        matrix_gz = "{assayID}/{runID}/{outdir}/{reference_version}/{application}/{tool}/{command}/{duplicates}/{referencePoint}/{contrast}/{ratio}/{source}.{region}.{mode}.{norm}.matrix.gz"
    output:
        figure = "{assayID}/{runID}/{outdir}/{reference_version}/{application}/{tool}/{command}/{duplicates}/{referencePoint}/{contrast}/{ratio}/{plotType}.{source}.{region}.{mode}.{norm}.pdf",
        data = "{assayID}/{runID}/{outdir}/{reference_version}/{application}/{tool}/{command}/{duplicates}/{referencePoint}/{contrast}/{ratio}/{plotType}.{source}.{region}.{mode}.{norm}.data",
        regions = "{assayID}/{runID}/{outdir}/{reference_version}/{application}/{tool}/{command}/{duplicates}/{referencePoint}/{contrast}/{ratio}/{plotType}.{source}.{region}.{mode}.{norm}.bed"
    shell:
        """
            {params.deepTools_dir}/plotProfile --matrixFile {input.matrix_gz} \
                                               --outFileName {output.figure} \
                                               --outFileNameData {output.data} \
                                               --outFileSortedRegions {output.regions} \
                                               --plotType {wildcards.plotType}
        """


# target rules
rule all:
    input:
        expand("{assayID}/{runID}/{outdir}/{reference_version}/{application}/{tool}/{mode}/{norm}/{duplicates}/{contrast}/{ratio}/{condition}.bw",
               assayID="ChIP-Seq",
               runID="merged",
               outdir=config["processed_dir"],
               reference_version=REF_VERSION,
               application=["deepTools"],
               tool=["bigwigCompare"],
               mode=["normal"],
               norm=["RPKM"],
               duplicates=["duplicates_removed"],
               contrast=["ChIP-Input"],
               ratio="log2",
               condition=["MCF10A_WT_HP1a",
                          "MCF10A_WT_HP1b",
                          "MCF10A_WT_H2AZ",
                          "MCF10A_shHP1a_HP1b",
                          "MCF10A_shHP1b_HP1a",
                          "MCF10A_shH2AZ_HP1a",
                          "MCF10A_shH2AZ_HP1b"]),
        expand("{assayID}/{runID}/{outdir}/{reference_version}/{application}/{tool}/{command}/{duplicates}/{referencePoint}/{contrast}/{ratio}/{plotType}.{source}.{region}.{mode}.{norm}.{suffix}",
               assayID="ChIP-Seq",
               runID="merged",
               outdir=config["processed_dir"],
               reference_version=REF_VERSION,
               application=["deepTools"],
               tool=["plotProfile"],
               command=["scale-regions"],
               duplicates=["duplicates_removed"],
               referencePoint="TSS",
               contrast=["ChIP-Input"],
               ratio="log2",
               plotType="se",
               source="bigwigCompare",
               region=["coding", "non-coding"],
               mode=["normal"],
               norm=["RPKM"],
               suffix=["pdf", "data", "bed"]),
        # expand("{assayID}/{runID}/{outdir}/{reference_version}/{application}/{tool}/{command}/{duplicates}/{referencePoint}/{plotType}.{source}.{region}.{mode}.{norm}.{suffix}",
        #        assayID="ChIP-Seq",
        #        runID="merged",
        #        outdir=config["processed_dir"],
        #        reference_version=REF_VERSION,
        #        application="deepTools",
        #        tool="plotProfile",
        #        command=["scale-regions", "reference-point"],
        #        duplicates=["duplicates_removed"],
        #        referencePoint="TSS",
        #        plotType="se",
        #        source = ["bamCoverage", "bigwigCompare"],
        #        region=["coding"],
        #        mode=["normal"],
        #        norm=["RPKM"],
        #        suffix=["pdf", "data", "bed"]),
