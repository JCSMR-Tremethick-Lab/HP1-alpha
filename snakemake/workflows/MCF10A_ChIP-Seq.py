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
include_prefix= os.environ['HOME'] + "/Development/JCSMR-Tremethick-Lab/Breast/snakemake/rules/"

# input functions
def getAllFASTQ(wildcards):
    fn =[]
    for i in config["samples"][wildcards["assayID"]][wildcards["runID"]]:
        for j in config["samples"][wildcards["assayID"]][wildcards["runID"]][i]:
            fn.append("/".join([wildcards["assayID"], wildcards["runID"], config["raw_dir"], j]))
    return(fn)

def getAllBAMs(duplicates):
    fn = []
    for i in config["samples"]["ChIP-Seq"]["runID"]:
        for j in config["samples"]["ChIP-Seq"][i]:
            fn.append(join("ChIP-Seq/" + i + "/" + config["processed_dir"] + "/" + REF_VERSION + "/bowtie2/" + duplicates + "/" + j + ".Q" + config["alignment_quality"] + ".sorted.bam"))
    return(fn)

def get_sample_labels(wildcards):
    sl = []
    runIDs = config["samples"][wildcards["assayID"]]["runID"]
    for i in runIDs:
        for k in config["samples"][wildcards["assayID"]][i].keys():
            sl.append(k)
    return(sl)


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
              file = [ i + "/" + config["processed_dir"] + "/" + REF_VERSION + "/bowtie2"  + "/" + j + ".bam" \
                  for i in config["samples"]["ChIP-Seq"]["runID"] \
                      for j in config["samples"]["ChIP-Seq"][i]])

PROCESSED_BAMs_dups_removed = expand("{assayID}/{file1}",
                                     assayID = "ChIP-Seq",
                                     file1 = [ i + "/" + config["processed_dir"] + "/" + REF_VERSION + "/bowtie2/duplicates_removed/" + j + ".Q" + config["alignment_quality"] + ".sorted.bam.bai" \
                                        for i in config["samples"]["ChIP-Seq"]["runID"] \
                                            for j in config["samples"]["ChIP-Seq"][i]]),

PROCESSED_BAMs_dups_marked = expand("{assayID}/{file2}",
                                     assayID = "ChIP-Seq",
                                     file2 = [ i + "/" + config["processed_dir"] + "/" + REF_VERSION + "/bowtie2/duplicates_marked/" + j + ".Q" + config["alignment_quality"] + ".sorted.bam.bai" \
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

# rule AdapterRemoval:
#     threads:
#         lambda wildcards: int(str(config["program_parameters"]["AdapterRemoval"]["threads"]).strip("['']"))
#     input:
#         read1 = lambda wildcards: wildcards["assayID"] + "/" + wildcards["runID"] + "/" + config["raw_dir"] + "/" + config["samples"][wildcards["assayID"]][wildcards["runID"]][wildcards["unit"]][0],
#         read2 = lambda wildcards: wildcards["assayID"] + "/" + wildcards["runID"] + "/" + config["raw_dir"] + "/" + config["samples"][wildcards["assayID"]][wildcards["runID"]][wildcards["unit"]][1]
#     output:
#         read1 = "{assayID}/{runID}/{outdir}/{trim_dir}/{unit}_R1.fastq.gz",
#         read2 = "{assayID}/{runID}/{outdir}/{trim_dir}/{unit}_R2.fastq.gz"
#     wrapper:
#         "file://" + wrapper_dir + "/AdapterRemoval/wrapper.py"
#
# rule fastqc:
#     version:
#         0.2
#     threads:
#         lambda wildcards: int(str(config["program_parameters"]["fastqc"]["threads"]).strip("['']"))
#     input:
#         getAllFASTQ
#     output:
#         "{assayID}/{runID}/{processed_dir}/{reports_dir}/"
#     shell:
#         home_dir + "/bin/fastqc {input} --noextract --threads {threads} --outdir {output}"

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

rule bam_quality_filter:
    # params:
    #     qual = config["alignment_quality"]
    input:
        "{assayID}/{runID}/{outdir}/" + REF_VERSION + "/bowtie2/{unit}.bam"
    output:
        temp("{assayID}/{runID}/{outdir}/{reference_version}/bowtie2/quality_filtered/{unit}.Q{qual}.bam")
    shell:
        "samtools view -b -h -q {params.qual} {input} > {output}"

rule bam_sort:
    params:
        qual = config["alignment_quality"]
    threads:
        4
    input:
        rules.bam_quality_filter.output
    output:
        "{assayID}/{runID}/{outdir}/{reference_version}/bowtie2/sorted/{unit}.Q{qual}.sorted.bam"
    shell:
        "samtools sort -@ {threads} {input} -T {wildcards.unit}.Q{params.qual}.sorted -o {output}"

rule bam_mark_duplicates:
    params:
        qual = config["alignment_quality"],
        picard = home + config["program_parameters"]["picard_tools"]["jar"],
        temp = home + config["temp_dir"]
    threads:
        4
    input:
        rules.bam_sort.output
    output:
        protected("{assayID}/{runID}/{outdir}/{reference_version}/bowtie2/duplicates_marked/{unit}.Q{qual}.sorted.bam")
    shell:
        """
            java -Djava.io.tmpdir={params.temp} \
            -Xmx24G \
            -jar {params.picard} MarkDuplicates \
            INPUT={input}\
            OUTPUT={output}\
            ASSUME_SORTED=true\
            METRICS_FILE={output}.metrics.txt
        """

rule bam_index:
    params:
        qual = config["alignment_quality"]
    input:
        rules.bam_mark_duplicates.output
    output:
        protected("{assayID}/{runID}/{outdir}/{reference_version}/bowtie2/duplicates_marked/{unit}.Q{qual}.sorted.bam.bai")
    shell:
        "samtools index {input} {output}"

rule bam_rmdup:
    input:
        rules.bam_mark_duplicates.output
    output:
        protected("{assayID}/{runID}/{outdir}/{reference_version}/bowtie2/duplicates_removed/{unit}.Q{qual}.sorted.bam")
    shell:
        "samtools rmdup {input} {output}"

rule bam_rmdup_index:
    params:
        qual = config["alignment_quality"]
    input:
        rules.bam_rmdup.output
    output:
        protected("{assayID}/{runID}/{outdir}/{reference_version}/bowtie2/duplicates_removed/{unit}.Q{qual}.sorted.bam.bai")
    shell:
        "samtools index {input} {output}"


# target rules
# rule run_AdapterRemoval:
#     input:
#         TRIMMED_FASTQ1,
#         TRIMMED_FASTQ2
#
# rule run_fastqc:
#     input:
#         expand("ChIP-Seq/{runID}/processed_data/reports/",
#                runID = config["samples"]["ChIP-Seq"]["runID"])

rule all:
    input:
        PROCESSED_BAMs_dups_removed,
        PROCESSED_BAMs_dups_marked
