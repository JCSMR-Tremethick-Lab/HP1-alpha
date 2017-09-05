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
include_prefix= os.environ['HOME'] + "/Development/JCSMR-Tremethick-Lab/HP1-alpha/snakemake/rules/"

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
        if i in config["samples"]["ChIP-Seq"]["conditions"][wildcards["condition"]].keys():
            for j in config["samples"]["ChIP-Seq"]["conditions"][wildcards["condition"]][i]:
                fn.append(join("ChIP-Seq/" + i + "/" + config["processed_dir"] + "/" + REF_VERSION + "/bowtie2/" + wildcards["duplicates"] + "/" + j + ".Q" + config["alignment_quality"] + ".sorted.bam"))
    return(fn)

def getAllBAMs(wildcards):
    fn = []
    for i in config["samples"]["ChIP-Seq"]["runID"]:
        for j in config["samples"]["ChIP-Seq"][i]:
            fn.append(join("ChIP-Seq/" + i + "/" + config["processed_dir"] + "/" + REF_VERSION + "/bowtie2/" + wildcards["duplicates"] + "/" + j + ".Q" + config["alignment_quality"] + ".sorted.bam"))
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

NPZ_FILES_marked = expand("{assayID}/{outdir}/{reference_version}/deepTools/multiBamSummary/{duplicates}/results.npz",
                          assayID = "ChIP-Seq",
                          outdir = "processed_data",
                          reference_version = REF_VERSION,
                          duplicates = ["duplicates_removed", "duplicates_marked"])

# set targets here
PROCESSED_BAMs = expand("{assayID}/{file}",
                        assayID = "ChIP-Seq",
                        file = [ i + "/" + config["processed_dir"] + "/" + REF_VERSION + "/bowtie2/duplicates_removed"  + "/" + j + ".Q" + config["alignment_quality"] + ".sorted.bam.bai" \
                            for i in config["samples"]["ChIP-Seq"]["runID"] \
                                for j in config["samples"]["ChIP-Seq"][i]])

BIGWIGs = expand("{assayID}/{file}.bw",
                 assayID = "ChIP-Seq",
                 file = [ i + "/" + config["processed_dir"] + "/" + REF_VERSION + "/deepTools/bamCoverage/normal/RPKM/duplicates_removed/"  + j + ".Q" + config["alignment_quality"]\
                    for i in config["samples"]["ChIP-Seq"]["runID"] \
                        for j in config["samples"]["ChIP-Seq"][i]])

# actual rules
rule multiBamSummary:
    version:
        0.2
    params:
        deepTools_dir = home + config["program_parameters"]["deepTools"]["deepTools_dir"],
        binSize = config["program_parameters"]["deepTools"]["binSize"],
        labels = get_sample_labels
    threads:
        16
    input:
        getAllBAMs
    output:
        npz = "{assayID}/{outdir}/{reference_version}/deepTools/multiBamSummary/{duplicates}/results.npz"
    shell:
        """
            {params.deepTools_dir}/multiBamSummary bins --bamfiles {input} \
                                                        --labels {params.labels} \
                                                        --numberOfProcessors {threads} \
                                                        --centerReads \
                                                        --binSize {params.binSize} \
                                                        --outFileName {output.npz}
        """

rule plotCorrelation_heatmap:
    params:
        deepTools_dir = home + config["program_parameters"]["deepTools"]["deepTools_dir"],
        plotTitle = lambda wildcards: "Correlation heatmap - " + wildcards["duplicates"]
    input:
        npz = "{assayID}/{outdir}/{reference_version}/deepTools/multiBamSummary/{duplicates}/results.npz"
    output:
        png = "{assayID}/{outdir}/{reference_version}/deepTools/plotCorrelation/{duplicates}/heatmap_SpearmanCorr_readCounts.pdf",
        tab = "{assayID}/{outdir}/{reference_version}/deepTools/plotCorrelation/{duplicates}/heatmap_SpearmanCorr_readCounts.tab"
    shell:
        """
            {params.deepTools_dir}/plotCorrelation --corData {input.npz} \
                                                   --corMethod spearman \
                                                   --skipZeros \
                                                   --plotTitle "{params.plotTitle}" \
                                                   --whatToPlot heatmap \
                                                   --colorMap RdYlBu \
                                                   --plotNumbers \
                                                   -o {output.png} \
                                                   --outFileCorMatrix {output.tab}
        """

rule plotPCA:
    params:
        deepTools_dir = home + config["program_parameters"]["deepTools"]["deepTools_dir"],
        plotTitle = lambda wildcards: "PCA - " + wildcards["duplicates"]
    input:
        npz = "{assayID}/{outdir}/{reference_version}/deepTools/multiBamSummary/{duplicates}/results.npz"
    output:
        png = "{assayID}/{outdir}/{reference_version}/deepTools/plotPCA/{duplicates}/PCA_readCounts.{output_format}"
    shell:
        """
            {params.deepTools_dir}/plotPCA --corData {input.npz} \
                                           --plotFile {output.png} \
                                           --plotTitle "{params.plotTitle}"
        """

rule bamPEFragmentSize:
    params:
        deepTools_dir = home + config["program_parameters"]["deepTools"]["deepTools_dir"],
        plotTitle = lambda wildcards: "BAM PE " + wildcards["duplicates"] + " fragment size",
        labels = get_sample_labels,
    threads:
        lambda wildcards: int(str(config["program_parameters"]["deepTools"]["threads"]).strip("['']"))
    input:
        getAllBAMs
    output:
        "{assayID}/{outdir}/{reference_version}/deepTools/bamPEFragmentSize/{duplicates}/histogram_duplicates_marked.{output_format}"
    shell:
        """
            {params.deepTools_dir}/bamPEFragmentSize --bamfiles {input} \
                                                     --samplesLabel {params.labels} \
                                                     --numberOfProcessors {threads} \
                                                     --histogram {output}
        """

rule plotFingerprint:
    params:
        deepTools_dir = home + config["program_parameters"]["deepTools"]["deepTools_dir"],
        plotTitle = lambda wildcards: "BAM PE " + wildcards["duplicates"] + " fingerprint",
        labels = getSampleLabelsByCondition
    threads:
        lambda wildcards: int(str(config["program_parameters"]["deepTools"]["threads"]).strip("['']"))
    input:
        getBAMbyCondition
    output:
        "{assayID}/{outdir}/{reference_version}/deepTools/plotFingerprint/{duplicates}/{condition}_fingerprints.{output_format}"
    shell:
        """
            {params.deepTools_dir}/plotFingerprint --bamfiles {input} \
                                                   --numberOfProcessors {threads} \
                                                   --centerReads \
                                                   --plotTitle "{params.plotTitle}" \
                                                   --labels {params.labels} \
                                                   --skipZeros \
                                                   --plotFile {output}
        """

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


# target rules
rule all:
    input:
        expand(["{assayID}/{outdir}/{reference_version}/deepTools/plotPCA/{duplicates}/PCA_readCounts.{output_format}",
                "{assayID}/{outdir}/{reference_version}/deepTools/plotCorrelation/{duplicates}/heatmap_SpearmanCorr_readCounts.{output_format}",
                "{assayID}/{outdir}/{reference_version}/deepTools/plotCorrelation/{duplicates}/heatmap_SpearmanCorr_readCounts.tab",
                "{assayID}/{outdir}/{reference_version}/deepTools/bamPEFragmentSize/{duplicates}/histogram_duplicates_marked.{output_format}",
                "{assayID}/{outdir}/{reference_version}/deepTools/plotFingerprint/{duplicates}/{condition}_fingerprints.{output_format}"],
               assayID = "ChIP-Seq",
               outdir = config["processed_dir"],
               reference_version = REF_VERSION,
               duplicates = ["duplicates_marked", "duplicates_removed"],
               condition = ["MCF10A_WT", "MCF10A_shHP1b", "MCF10A_shH2AZ", "MCF10A_shHP1a"],
               output_format = ["pdf", "png"]),
        BIGWIGs
