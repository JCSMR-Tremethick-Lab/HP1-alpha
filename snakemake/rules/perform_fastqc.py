__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2016-08-08"

from snakemake.exceptions import MissingInputException
import os

def getAllFASTQ(wildcards):
    fn =[]
    for i in config["samples"][wildcards["assayID"]][wildcards["runID"]]:
        for j in config["samples"][wildcards["assayID"]][wildcards["runID"]][i]:
            fn.append("RNA-Seq/NB501086_0114_B_Azad_JCSMR_hRNAseq/fastq/" + j)
    return(fn)

rule dummy:
    input:
        "RNA-Seq/NB501086_0114_B_Azad_JCSMR_hRNAseq/processed_data/reports/"

rule fastqc:
    version:
        0.3
    threads:
        16
    input:
        getAllFASTQ
    output:
        "{assayID}/{runID}/{processed_dir}/{reports_dir}/"
    shell:
        "fastqc {input} --threads {threads} --noextract --outdir  {output}"
