_author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2015-04-22"
__modified__ = "2017-05-11"

from snakemake.exceptions import MissingInputException
import os
from os.path import join


rule:
    version: 0.1

localrules:
    all, run_STAR

# variables
REF_GENOME="hg19"
REF_VERSION="GRCh37_hg19_ensembl75"
HOME=os.environ['HOME']

wrapper_dir = os.environ['HOME'] + "/Development/snakemake-wrappers/bio"

include_prefix= os.environ['HOME'] + "/Development/JCSMR-Tremethick-Lab/HP1-alpha/snakemake/rules/"

rule star_align_full:
    version:
        0.6
    threads:
        16
    params:
        trim_dir = config["trim_dir"],
        genomeLoad = "LoadAndKeep",
        limitBAMsortRAM = 16000000000,
    input:
        read1 = "{assayID}/{runID}/{processed_dir}/trimmed_data/{unit}_R1_001.QT.CA.fastq.gz",
        read2 = "{assayID}/{runID}/{processed_dir}/trimmed_data/{unit}_R2_001.QT.CA.fastq.gz",
        index = lambda wildcards: HOME + config["references"][REF_GENOME]["STAR"][wildcards["reference_version"]]
    output:
        "{assayID}/{runID}/{processed_dir}/{reference_version}/STAR/full/{unit}.bam"
    shell:
        """
            STAR --runMode alignReads \
                 --runThreadN {threads} \
                 --genomeDir {input.index} \
                 --readFilesIn {input.read1} {input.read2} \
                 --readFilesCommand zcat \
                 --outTmpDir /home/sebastian/tmp/{wildcards.unit} \
                 --outSAMmode Full \
                 --outSAMattributes Standard \
                 --outSAMtype BAM SortedByCoordinate \
                 --outStd BAM_SortedByCoordinate \
                 --alignEndsType EndToEnd\
                 > {output}
        """

rule star_align_IRFinder:
    version:
        0.1
    threads:
        16
    params:
        trim_dir = config["trim_dir"],
        genomeLoad = "LoadAndKeep",
        limitBAMsortRAM = 16000000000,
    input:
        read1 = "{assayID}/{runID}/{processed_dir}/trimmed_data/{unit}_R1_001.QT.CA.fastq.gz",
        read2 = "{assayID}/{runID}/{processed_dir}/trimmed_data/{unit}_R2_001.QT.CA.fastq.gz",
        index = lambda wildcards: HOME + config["references"][REF_GENOME]["STAR"][wildcards["reference_version"]]
    output:
        "{assayID}/{runID}/{processed_dir}/{reference_version}/IRFinder/{unit}.bam"
    shell:
        """
            STAR --runMode alignReads \
                 --runThreadN {threads} \
                 --genomeDir {input.index} \
                 --readFilesIn {input.read1} {input.read2} \
                 --readFilesCommand zcat \
                 --outTmpDir /home/sebastian/tmp/{wildcards.unit} \
                 --outSAMmode Full \
                 --outSAMattributes Standard \
                 --outSAMtype BAM SortedByCoordinate \
                 --outStd BAM_SortedByCoordinate \
                 --alignEndsType EndToEnd\
                 > {output}
        """

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


def getBAMbyReplicatesSTAR(wildcards):
    fn=[]
    for i in config["samples"]["RNA-Seq"]["runID"]:
        for j in config["samples"]["RNA-Seq"][i]:
            if j in config["samples"]["RNA-Seq"]["replicates"][wildcards["replicates"]]:
                fn.append(join("RNA-Seq/" +
                               i +
                               "/" +
                               config["processed_dir"] +
                               "/" +
                               REF_VERSION +
                               "/STAR/full/" +
                               j + ".bam"))
    return(fn)


def getComputeMatrixInput(wildcards):
    fn = []
    for j in config["samples"]["RNA-Seq"]["replicates"].keys():
        fn.append("/".join(["RNA-Seq/merged",
                            config["processed_dir"],
                            REF_VERSION,
                            "deepTools",
                            "bamCoverage",
                            wildcards["mode"],
                            wildcards["norm"],
                            j + ".bw"]))
    return(fn)


# rules
rule bam_index_STAR_output:
    version:
        0.2
    input:
        "{assayID}/{runID}/{processed_dir}/{reference_version}/STAR/full/{unit}.bam"
    output:
        "{assayID}/{runID}/{processed_dir}/{reference_version}/STAR/full/{unit}.bam.bai"
    wrapper:
        "file://" + wrapper_dir + "/samtools/index/wrapper.py"

rule bamMerge:
    version:
        0.1
    params:
        outputFormat = "--output-fmt BAM"  # ToDo: move
    log:
        "{assayID}/merged/{outdir}/{reference_version}/{replicates}.bamMerge.log"
    threads:
        16
    input:
        getBAMbyReplicatesSTAR
    output:
        "{assayID}/merged/{outdir}/{reference_version}/{replicates}.bam"
    shell:
        """
            samtools merge -f {output} {input} --threads {threads} {params.outputFormat} 1>>{log} 2>>{log}
        """

rule indexMerged:
    version:
        0.1
    log:
        "{assayID}/merged/{outdir}/{reference_version}/{replicates}.indexMerged.log"
    threads:
        16
    input:
        rules.bamMerge.output
    output:
        "{assayID}/merged/{outdir}/{reference_version}/{replicates}.bam.bai"
    shell:
        """
            samtools index {input} {output} -@ {threads} 1>>{log} 2>>{log}
        """

rule bamCoverageMerged:
    version:
        0.1
    params:
        deepTools_dir = HOME + config["program_parameters"]["deepTools"]["deepTools_dir"],
        ignore = config["program_parameters"]["deepTools"]["ignoreForNormalization"],
        program_parameters = cli_parameters_bamCoverage
    threads:
        lambda wildcards: int(str(config["program_parameters"]["deepTools"]["threads"]).strip("['']"))
    input:
        bam = lambda wildcards: wildcards["assayID"] + "/merged/" + wildcards["outdir"] + "/" + wildcards["reference_version"] + "/" +  wildcards["replicates"] + ".bam"
    output:
        "{assayID}/merged/{outdir}/{reference_version}/{application}/{tool}/{mode}/{norm}/{replicates}.bw"
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
        deepTools_dir = HOME + config["program_parameters"]["deepTools"]["deepTools_dir"],
        program_parameters = lambda wildcards: ' '.join("{!s}={!s}".format(key, val.strip("\\'")) for (key, val) in cli_parameters_computeMatrix(wildcards).items())
    threads:
        lambda wildcards: int(str(config["program_parameters"]["deepTools"]["threads"]).strip("['']"))
    input:
        file = getComputeMatrixInput,
        region = lambda wildcards: HOME + config["program_parameters"]["deepTools"]["regionFiles"][wildcards["reference_version"]][wildcards["region"]]
    output:
        matrix_gz = "{assayID}/{outdir}/{reference_version}/{application}/{tool}/{command}/{referencePoint}/{region}_{mode}_{norm}.matrix.gz"
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
        matrix_gz = "{assayID}/{outdir}/{reference_version}/{application}/computeMatrix/{command}/{referencePoint}/{region}_{mode}_{norm}.matrix.gz"
    output:
        figure = "{assayID}/{outdir}/{reference_version}/{application}/{tool}/{command}/{referencePoint}/{plotType}.{mode}.{norm}.{region}.pdf",
        data = "{assayID}/{outdir}/{reference_version}/{application}/{tool}/{command}/{referencePoint}/{plotType}.{mode}.{norm}.{region}.data",
        regions = "{assayID}/{outdir}/{reference_version}/{application}/{tool}/{command}/{referencePoint}/{plotType}.{mode}.{norm}.{region}.bed"
    shell:
        """
            {params.deepTools_dir}/plotProfile --matrixFile {input.matrix_gz} \
                                               --outFileName {output.figure} \
                                               --outFileNameData {output.data} \
                                               --outFileSortedRegions {output.regions} \
                                               --plotType {wildcards.plotType}
        """

# targets
STARS = expand("{assayID}/{file}",
              assayID = "RNA-Seq",
              file = [ i + "/" + config["processed_dir"] + "/" + REF_VERSION + "/STAR/full/" + j + ".bam.bai"\
                  for i in config["samples"]["RNA-Seq"]["runID"] \
                      for j in config["samples"]["RNA-Seq"][i]])

MERGED_STARS = expand("{assayID}/merged/{outdir}/{reference_version}/{replicates}.bam.bai",
               assayID="RNA-Seq",
               outdir=config["processed_dir"],
               reference_version=REF_VERSION,
               replicates=['MCF10A_WT',
                           'MCF10A_Scramble',
                           'MCF10A_shH2AZHP1a',
                           'MCF10A_shHP1a',
                           'MCF10A_shHP1ab',
                           'MCF10A_shHP1b'])

MergedBIGWIGs = expand("{assayID}/{file}.bw",
                         assayID="RNA-Seq",
                         file=["merged/" + config["processed_dir"] + "/" + REF_VERSION + "/deepTools/bamCoverage/normal/RPKM/"  + j \
                                for j in config["samples"]["RNA-Seq"]["replicates"].keys()])

PLOTs = expand("{assayID}/{outdir}/{reference_version}/{application}/{tool}/{command}/{referencePoint}/{plotType}.{mode}.{norm}.{region}.{suffix}",
               assayID="RNA-Seq",
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
               suffix=["pdf", "data", "bed"],
               region = ["allGenes", "intergenicRegions"])

rule all:
    input:
        MergedBIGWIGs
