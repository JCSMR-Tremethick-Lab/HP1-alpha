_author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2015-04-22"
__modified__ = "2017-05-11"

from snakemake.exceptions import MissingInputException
import os

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
        8
    params:
        trim_dir = config["trim_dir"]
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

rule bam_index_STAR_output:
    version:
        0.2
    input:
        "{assayID}/{runID}/{processed_dir}/{reference_version}/STAR/full/{unit}.bam"
    output:
        "{assayID}/{runID}/{processed_dir}/{reference_version}/STAR/full/{unit}.bam.bai"
    wrapper:
        "file://" + wrapper_dir + "/samtools/index/wrapper.py"

# targets
STARS = expand("{assayID}/{file}",
              assayID = "RNA-Seq",
              file = [ i + "/" + config["processed_dir"] + "/" + REF_VERSION + "/STAR/full/" + j + ".bam.bai"\
                  for i in config["samples"]["RNA-Seq"]["runID"] \
                      for j in config["samples"]["RNA-Seq"][i]])

rule all:
    input:
        STARS
