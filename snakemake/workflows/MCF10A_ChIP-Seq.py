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

wrapper_dir = os.environ['HOME'] + "/Development/snakemake-wrappers/bio"

include_prefix= os.environ['HOME'] + "/Development/JCSMR-Tremethick-Lab/Breast/snakemake/rules/"

include:
    include_prefix + "run_kallisto.py"

# set targets here
TRIMMED_FASTQ1 = expand("ChIP-Seq/{file}",
                        file = [ i + "/" + config["processed_dir"] + "/" + config["trim_dir"] + "/" + j + "_R1.fastq.gz" \
                            for i in config["samples"]["ChIP-Seq"]["runID"] \
                                for j in config["samples"]["ChIP-Seq"][i]])

TRIMMED_FASTQ2 = expand("ChIP-Seq/{file}",
                        file = [ i + "/" + config["processed_dir"] + "/" + config["trim_dir"] + "/" + j + "_R2.fastq.gz" \
                            for i in config["samples"]["ChIP-Seq"]["runID"] \
                                for j in config["samples"]["ChIP-Seq"][i]])

rule AdapterRemoval:
    params:
        threads = config["AdapterRemoval"]["threads"]
    input:
        read1 = config["samples"][wildcards["assayID"]][wildcards["runID"]][wildcards["unit"]][0],
        read2 = config["samples"][wildcards["assayID"]][wildcards["runID"]][wildcards["unit"]][1]
    output:
        read1 = "{assayID}/{runID}/{outdir}/{trim_dir}/{unit}_R1.fastq.gz",
        read2 = "{assayID}/{runID}/{outdir}/{trim_dir}/{unit}_R2.fastq.gz"
    wrapper:
        "0.1.0/bio/AdapterRemoval"

rule run_AdapterRemoval:
    input:
        TRIMMED_FASTQ1,
        TRIMMED_FASTQ2
