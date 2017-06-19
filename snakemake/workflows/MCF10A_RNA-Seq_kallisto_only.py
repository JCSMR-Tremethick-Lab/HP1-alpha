_author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2015-04-22"

from snakemake.exceptions import MissingInputException
import os

rule:
    version: 0.3

localrules:
    all, run_kallisto, run_STAR, run_htseq, run_cutadapt

wrapper_dir = os.environ['HOME'] + "/Development/snakemake-wrappers/bio"

include_prefix= os.environ['HOME'] + "/Development/JCSMR-Tremethick-Lab/Breast/snakemake/rules/"

include:
    include_prefix + "run_kallisto.py"

rule run_kallisto:
    input:
        expand("{assayID}/{runID}/{outdir}/{reference_version}/kallisto/{unit}",
               assayID = "RNA-Seq",
               runID = ["NB501086_0114_B_Azad_JCSMR_hRNAseq", "NB501086_0082_RDomaschenz_JCSMR_mRNAseq"],
               outdir = config["processed_dir"],
               reference_version = config["references"]["hg38"]["version"],
               unit = config["samples"]["RNA-Seq"])
