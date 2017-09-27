_author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2017-08-29"

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
include_prefix= os.environ['HOME'] + "/Development/JCSMR-Tremethick-Lab/HP1-alpha/snakemake/rules/"

# target files
BIGWIGs = expand("{assayID}/{file}.bw",
                 assayID = "ChIP-Seq",
                 file = [ i + "/" + config["processed_dir"] + "/" + REF_VERSION + "/deepTools/bamCoverage/normal/RPKM/duplicates_removed/"  + j + ".Q" + config["alignment_quality"]\
                    for i in config["samples"]["ChIP-Seq"]["runID"] \
                        for j in config["samples"]["ChIP-Seq"][i]])


# input functions
def getAllFASTQ(wildcards):
    fn =[]
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
        a["--referencePoint"] = wildcards["referencePoint"]
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

# rules
rule bamCoverage:
    version:
        0.1
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

rule computeMatrix:
    version:
        0.2
    params:
        deepTools_dir = home + config["program_parameters"]["deepTools"]["deepTools_dir"],
        program_parameters = lambda wildcards: ' '.join("{!s}={!s}".format(key, val.strip("\\'")) for (key, val) in cli_parameters_computeMatrix(wildcards).items())
    threads:
        lambda wildcards: int(str(config["program_parameters"]["deepTools"]["threads"]).strip("['']"))
    input:
        file = getComputeMatrixInput,
        region = lambda wildcards: home + config["program_parameters"]["deepTools"]["regionFiles"][wildcards["reference_version"]][wildcards["region"]]
    output:
        matrix_gz = "{assayID}/{outdir}/{reference_version}/{application}/{tool}/{command}/{duplicates}/{referencePoint}/{region}_{mode}_{norm}.matrix.gz"
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

rule plotProfile:
    version:
        0.1
    params:
        deepTools_dir = home + config["program_parameters"]["deepTools"]["deepTools_dir"],
    input:
        matrix_gz = "{assayID}/{outdir}/{reference_version}/{application}/computeMatrix/{command}/{duplicates}/{referencePoint}/{region}_{mode}_{norm}.matrix.gz"
    output:
        figure = "{assayID}/{outdir}/{reference_version}/{application}/{tool}/{command}/{duplicates}/{referencePoint}/{plotType}.{mode}.{norm}.{region}.pdf",
        data = "{assayID}/{outdir}/{reference_version}/{application}/{tool}/{command}/{duplicates}/{referencePoint}/{plotType}.{mode}.{norm}.{region}.data",
        regions = "{assayID}/{outdir}/{reference_version}/{application}/{tool}/{command}/{duplicates}/{referencePoint}/{plotType}.{mode}.{norm}.{region}.bed"
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
        BIGWIGs,
        expand("{assayID}/{outdir}/{reference_version}/{application}/{tool}/{command}/{duplicates}/{referencePoint}/{plotType}.{mode}.{norm}.{region}.{suffix}",
               assayID = "ChIP-Seq",
               outdir = config["processed_dir"],
               reference_version = REF_VERSION,
               application = "deepTools",
               tool = "plotProfile",
               command = ["scale-regions"],
               duplicates = ["duplicates_removed"],
               referencePoint = "TSS",
               plotType = "se",
               mode = ["normal"],
               norm = ["RPKM"],
               region = ["allGenes", "intergenicRegions"],
               suffix = ["pdf", "data", "bed"])
