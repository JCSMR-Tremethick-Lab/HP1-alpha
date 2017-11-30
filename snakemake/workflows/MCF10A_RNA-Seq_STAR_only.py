_author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2015-04-22"
__modified__ = "2017-05-11"

from snakemake.exceptions import MissingInputException
import os

rule:
    version: 0.4

localrules:
    all, run_kallisto, run_STAR, run_htseq, run_cutadapt

home = os.environ['HOME']

wrapper_dir = os.environ['HOME'] + "/Development/snakemake-wrappers/bio"

include_prefix= os.environ['HOME'] + "/Development/JCSMR-Tremethick-Lab/HP1-alpha/snakemake/rules/"

assay_dir = "NB501086_0114_B_Azad_JCSMR_hRNAseq"

include:
     include_prefix + "perform_fastqc.py"
include:
    include_prefix + "perform_cutadapt.py"
include:
    include_prefix + "run_kallisto.py"
include:
    include_prefix + "run_STAR.py"

rule run_fastqc:
    input:
        "RNA-Seq/NB501086_0114_B_Azad_JCSMR_hRNAseq/processed_data/reports/"

rule run_kallisto:
    input:
        expand("{assayID}/NB501086_0114_B_Azad_JCSMR_hRNAseq/{outdir}/{reference_version}/kallisto/{unit}",
               assayID = "RNA-Seq",
               outdir = config["processed_dir"],
               #reference_version = config["references"]["hg38"]["version"],
               reference_version = "GRCh38_RefSeq_NM",
               unit = config["samples"]["RNA-Seq"]["NB501086_0114_B_Azad_JCSMR_hRNAseq"].keys()),
        expand("{assayID}/NB501086_0082_RDomaschenz_JCSMR_mRNAseq/{outdir}/{reference_version}/kallisto/{unit}",
               assayID = "RNA-Seq",
               outdir = config["processed_dir"],
               #reference_version = config["references"]["hg38"]["version"],
               reference_version = "GRCh38_RefSeq_NM",
               unit = config["samples"]["RNA-Seq"]["NB501086_0082_RDomaschenz_JCSMR_mRNAseq"].keys())

rule run_STAR:
    input:
        expand("{assayID}/NB501086_0114_B_Azad_JCSMR_hRNAseq/{outdir}/{reference_version}/STAR/full/{unit}.aligned.bam",
               assayID = "RNA-Seq",
               outdir = config["processed_dir"],
               reference_version = config["references"]["hg38"]["version"],
               unit = config["samples"]["RNA-Seq"]["NB501086_0114_B_Azad_JCSMR_hRNAseq"].keys())

rule run_STAR_untrimmed:
    input:
        expand("{assayID}/{runID}/{outdir}/{reference_version}/STAR/full/untrimmed/{unit}.aligned.bam",
               assayID = "RNA-Seq",
               runID = "NB501086_0114_B_Azad_JCSMR_hRNAseq",
               outdir = config["processed_dir"],
               reference_version = config["references"]["hg38"]["version"],
               unit = config["samples"]["RNA-Seq"]["NB501086_0114_B_Azad_JCSMR_hRNAseq"].keys())

rule run_cutadapt:
    input:
        expand("{assayID}/{runID}/{outdir}/{trim_data}/{unit}_{suffix}.QT.CA.fastq.gz",
               assayID = "RNA-Seq",
               runID = "NB501086_0114_B_Azad_JCSMR_hRNAseq",
               outdir = config["processed_dir"],
               trim_data = config["trim_dir"],
               unit = config["samples"]["RNA-Seq"]["NB501086_0114_B_Azad_JCSMR_hRNAseq"].keys(),
               suffix = ["R1_001", "R2_001"])

rule all:
    input:
        expand("{assayID}/{runID}/{outdir}/{reference_version}/kallisto/{unit}",
               assayID = "RNA-Seq",
               runID = "NB501086_0114_B_Azad_JCSMR_hRNAseq",
               outdir = config["processed_dir"],
               reference_version = config["references"]["hg38"]["version"],
               unit = config["samples"]["RNA-Seq"]["NB501086_0114_B_Azad_JCSMR_hRNAseq"].keys()),
        expand("{assayID}/{runID}/{outdir}/{reference_version}/kallisto/{unit}",
               assayID = "RNA-Seq",
               runID = "NB501086_0082_RDomaschenz_JCSMR_mRNAseq",
               outdir = config["processed_dir"],
               reference_version = config["references"]["hg38"]["version"],
               unit = config["samples"]["RNA-Seq"]["NB501086_0082_RDomaschenz_JCSMR_mRNAseq"].keys()),
        expand("{assayID}/{runID}/{outdir}/{reference_version}/STAR/full/{unit}.aligned.bam",
               assayID = "RNA-Seq",
               runID = "NB501086_0114_B_Azad_JCSMR_hRNAseq",
               outdir = config["processed_dir"],
               reference_version = config["references"]["hg38"]["version"],
               unit = config["samples"]["RNA-Seq"]["NB501086_0114_B_Azad_JCSMR_hRNAseq"].keys()),
