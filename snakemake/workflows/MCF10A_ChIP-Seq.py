_author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2017-08-23"

from snakemake.exceptions import MissingInputException
import os

rule:
    version: 0.1

localrules:
    all

wrapper_dir = os.environ['HOME'] + "/Development/snakemake-wrappers/bio"

include_prefix= os.environ['HOME'] + "/Development/JCSMR-Tremethick-Lab/Breast/snakemake/rules/"

include:
    include_prefix + "run_kallisto.py"

rule AdapterRemoval:
    params:
    input:
    output:
    run:

rule all:
    input:
        expand("{assayID}/{runID}/{outdir}/{trim_dir}/{unit}",
               assayID = "ChIP-Seq",
               runID = config["ChIP-Seq"]["runID"],
               outdir = config["processed_dir"],
               trim_dir = config["trim_dir"]
               unit = config["samples"]["ChIP-Seq"][config["ChIP-Seq"]["runID"][0]])
