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
                           'MCF10A_shHP1b']),
rule all:
    input:
        STARS,
	MERGED_STARS
