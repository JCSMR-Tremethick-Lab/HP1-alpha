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

wrapper_dir = os.environ['HOME'] + "/Development/snakemake-wrappers/bio"

REF_GENOME = "hg19"
REF_VERSION = config["references"][REF_GENOME]["version"][1]

# includes
include_prefix= os.environ['HOME'] + "/Development/JCSMR-Tremethick-Lab/Breast/snakemake/rules/"

include:
    include_prefix + "run_kallisto.py"

# input functions
def getAllFASTQ(wildcards):
    fn =[]
    for i in config["samples"][wildcards["assayID"]][wildcards["runID"]]:
        for j in config["samples"][wildcards["assayID"]][wildcards["runID"]][i]:
            fn.append("/".join([wildcards["assayID"], wildcards["runID"], config["raw_dir"], j]))
    return(fn)

# set targets here
TRIMMED_FASTQ1 = expand("ChIP-Seq/{file}",
                        file = [ i + "/" + config["processed_dir"] + "/" + config["trim_dir"] + "/" + j + "_R1.fastq.gz" \
                            for i in config["samples"]["ChIP-Seq"]["runID"] \
                                for j in config["samples"]["ChIP-Seq"][i]])

TRIMMED_FASTQ2 = expand("ChIP-Seq/{file}",
                        file = [ i + "/" + config["processed_dir"] + "/" + config["trim_dir"] + "/" + j + "_R2.fastq.gz" \
                            for i in config["samples"]["ChIP-Seq"]["runID"] \
                                for j in config["samples"]["ChIP-Seq"][i]])

FASTQC_OUTPUT = expand("ChIP-Seq/{file}",
                        file = [ i + "/" + config["processed_dir"] + "/" + config["reports_dir"] + "/" + k \
                            for i in config["samples"]["ChIP-Seq"]["runID"] \
                                for j in config["samples"]["ChIP-Seq"][i].values() \
                                    for k in j])

BAMs = expand("{assayID}/{file}",
              assayID = "ChIP-Seq",
              file = [ i + "/" + config["processed_dir"] + "/" + REF_VERSION + "/bowtie2/"  + "/" + j + ".bam" \
                  for i in config["samples"]["ChIP-Seq"]["runID"] \
                      for j in config["samples"]["ChIP-Seq"][i]])

# rule move_fastq:
#     input:
#         "{assayID}/{runID}/{unit}"
#     output:
#         "{assayID}/{runID}/{raw_dir}/{unit}"
#     run:
#        """
#             mv {input} {output}
#        """

rule AdapterRemoval:
    threads:
        lambda wildcards: int(str(config["program_parameters"]["AdapterRemoval"]["threads"]).strip("['']"))
    input:
        read1 = lambda wildcards: wildcards["assayID"] + "/" + wildcards["runID"] + "/" + config["raw_dir"] + "/" + config["samples"][wildcards["assayID"]][wildcards["runID"]][wildcards["unit"]][0],
        read2 = lambda wildcards: wildcards["assayID"] + "/" + wildcards["runID"] + "/" + config["raw_dir"] + "/" + config["samples"][wildcards["assayID"]][wildcards["runID"]][wildcards["unit"]][1]
    output:
        read1 = "{assayID}/{runID}/{outdir}/{trim_dir}/{unit}_R1.fastq.gz",
        read2 = "{assayID}/{runID}/{outdir}/{trim_dir}/{unit}_R2.fastq.gz"
    wrapper:
        "file://" + wrapper_dir + "/AdapterRemoval/wrapper.py"

rule fastqc:
    version:
        0.2
    threads:
        lambda wildcards: int(str(config["program_parameters"]["fastqc"]["threads"]).strip("['']"))
    input:
        getAllFASTQ
    output:
        "{assayID}/{runID}/{processed_dir}/{reports_dir}/"
    shell:
        home_dir + "/bin/fastqc {input} --noextract --threads {threads} --outdir {output}"

rule bowtie2_pe:
    version:
        "0.2"
    params:
        max_in = config["program_parameters"]["bowtie2"]["max_insert"],
        sample = lambda wildcards: wildcards["unit"][:-1],
        bt2_index = home_dir + config["references"][REF_GENOME]["bowtie2"][REF_VERSION]
    threads:
        lambda wildcards: int(str(config["program_parameters"]["bowtie2"]["threads"]).strip("['']"))
    input:
        read1="{assayID}/{runID}/{outdir}/trimmed_data/{unit}_R1.fastq.gz",
        read2="{assayID}/{runID}/{outdir}/trimmed_data/{unit}_R2.fastq.gz",
    output:
        protected("{assayID}/{runID}/{outdir}/" + REF_VERSION + "/bowtie2/{unit}.bam")
    shell:
        """
            bowtie2 \
            -x {params.bt2_index}\
            --no-mixed \
            --no-discordant \
            --maxins {params.max_in} \
            --threads {threads}\
            --rg-id '{wildcards.unit}' \
            --rg 'LB:{wildcards.unit}' \
            --rg 'SM:{params.sample}' \
            --rg 'PL:Illumina' \
            --rg 'PU:NA' \
            -1 {input.read1} \
            -2 {input.read2} \
            | samtools view -Sb - > {output}
        """

rule run_AdapterRemoval:
    input:
        TRIMMED_FASTQ1,
        TRIMMED_FASTQ2

rule run_fastqc:
    input:
        expand("ChIP-Seq/{runID}/processed_data/reports/",
               runID = config["samples"]["ChIP-Seq"]["runID"])

rule all:
    input:
        TRIMMED_FASTQ1,
        TRIMMED_FASTQ2,
        BAMs
