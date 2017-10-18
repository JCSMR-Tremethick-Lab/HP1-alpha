__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2016-04-10"

from snakemake.exceptions import MissingInputException

#wrapper_dir = "/home/sebastian/Development/snakemake-wrappers/bio"

rule kallisto_quant:
    params:
        bootstraps = config["program_parameters"]["kallisto"]["bootstraps"],
        trim_dir = config["trim_dir"]
    threads:
        lambda wildcards: int(str(config["program_parameters"]["kallisto"]["threads"]).strip("['']"))
    input:
        read1 = "{assayID}/{runID}/{processed_dir}/trimmed_data/{unit}_R1_001.QT.CA.fastq.gz",
        read2 = "{assayID}/{runID}/{processed_dir}/trimmed_data/{unit}_R2_001.QT.CA.fastq.gz",
        ki = lambda wildcards: HOME + config["references"][REF_GENOME]["kallisto"][wildcards["reference_version"]]
    output:
        protected("{assayID}/{runID}/{processed_dir}/{reference_version}/kallisto/{unit}")
    shell:
        """
            kallisto quant --index={input.ki} \
                           --output-dir={output} \
                           --threads={threads} \
                           --bootstrap-samples={params.bootstraps} \
                           {input.read1} {input.read2}
        """
