
# rules
rule bamMerge:
    version:
        1
    conda:
        "./envs/bam.yaml"
    params:
        outputFormat = "--output-fmt BAM"  # ToDo: move
    log:
        "logs/{replicates}.bamMerge.log"
    threads:
        16
    input:
        getBAMbyReplicates
    output:
        "{assayID}/merged/{outdir}/{reference_version}/{duplicates}/{replicates}.bam"
    shell:
        """
            samtools merge -f {output} {input} --threads {threads} {params.outputFormat} 1>>{log} 2>>{log}
        """

rule indexMerged:
    version:
        0.1
    log:
        "logs/{replicates}.indexMerged.log"
    threads:
        16
    input:
        rules.bamMerge.output
    output:
        "{assayID}/merged/{outdir}/{reference_version}/{duplicates}/{replicates}.bam.bai"
    shell:
        """
            samtools index {input} {output} -@ {threads} 1>>{log} 2>>{log}
        """

# rules section 
rule bamCoverageMerged:
    version:
        0.1
    params:
        deepTools_dir = home + config["program_parameters"]["deepTools"]["deepTools_dir"],
        ignore = config["program_parameters"]["deepTools"]["ignoreForNormalization"],
        program_parameters = cliParametersBamCoverage
    threads:
        lambda wildcards: int(str(config["program_parameters"]["deepTools"]["threads"]).strip("['']"))
    input:
        bai = "{assayID}/merged/{outdir}/{reference_version}/{duplicates}/{replicates}.bam.bai",
        bam = "{assayID}/merged/{outdir}/{reference_version}/{duplicates}/{replicates}.bam"
    output:
        "{assayID}/merged/{outdir}/{reference_version}/{application}/{tool}/{mode}/{norm}/{duplicates}/{replicates}.bw"
    shell:
        """
            {params.deepTools_dir}/bamCoverage --bam {input.bam} \
                                               --outFileName {output} \
                                               --outFileFormat bigwig \
                                               {params.program_parameters} \
                                               --numberOfProcessors {threads} \
                                               --normalizeUsingRPKM \
                                               --ignoreForNormalization {params.ignore}
        """
