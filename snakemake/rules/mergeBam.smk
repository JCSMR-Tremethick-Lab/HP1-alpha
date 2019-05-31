rule bamMerge:
    version:
        1
    conda:
        "../envs/bam.yaml"
    params:
        outputFormat = "--output-fmt BAM"  # ToDo: move
    threads:
        16
    input:
        getBAMbyReplicates
    output:
        "{assayID}/merged/{outdir}/{reference_version}/{duplicates}/{replicates}.bam"
    shell:
        """
            samtools merge -f {output} {input} --threads {threads} {params.outputFormat}
        """

rule indexMerged:
    version:
        1
    conda:
        "../envs/bam.yaml"
    threads:
        16
    input:
        rules.bamMerge.output
    output:
        "{assayID}/merged/{outdir}/{reference_version}/{duplicates}/{replicates}.bam.bai"
    shell:
        """
            samtools index {input} {output} -@ {threads}
        """

# rules section 
rule bamCoverageMerged:
    version:
        1
    conda:
        "../envs/deeptools.yaml"
    params:
        ignore = config["program_parameters"]["deepTools"]["ignoreForNormalization"],
        program_parameters = cliParametersBamCoverage
    threads:
        16
    input:
        bai = "{assayID}/merged/{outdir}/{reference_version}/{duplicates}/{replicates}.bam.bai",
        bam = "{assayID}/merged/{outdir}/{reference_version}/{duplicates}/{replicates}.bam"
    output:
        "{assayID}/merged/{outdir}/{reference_version}/{application}/{tool}/{mode}/{norm}/{duplicates}/{replicates}.bw"
    shell:
        """
            bamCoverage --bam {input.bam} \
                        --outFileName {output} \
                        --outFileFormat bigwig \
                        {params.program_parameters} \
                        --numberOfProcessors {threads} \
                        --normalizeUsing RPKM \
                        --ignoreForNormalization {params.ignore}
        """

rule bigwigCompareMerged:
    version:
        1
    conda:
        "../envs/deeptools.yaml"
    threads:
        32
    input:
        chip = lambda wildcards: "/".join((wildcards["assayID"],
                                           "merged",
                                           wildcards["outdir"],
                                           wildcards["reference_version"],
                                           wildcards["application"],
                                           "bamCoverage",
                                           wildcards["mode"],
                                           wildcards["norm"],
                                           wildcards["duplicates"],
                                           ".".join((config["samples"]["ChIP-Seq"]["ChIP-Input"][wildcards["condition"]]["ChIP"], "bw")))),
        input = lambda wildcards: "/".join((wildcards["assayID"],
                                           "merged",
                                           wildcards["outdir"],
                                           wildcards["reference_version"],
                                           wildcards["application"],
                                           "bamCoverage",
                                           wildcards["mode"],
                                           wildcards["norm"],
                                           wildcards["duplicates"],
                                           ".".join((config["samples"]["ChIP-Seq"]["ChIP-Input"][wildcards["condition"]]["Input"], "bw"))))
    output:
        "{assayID}/{runID}/{outdir}/{reference_version}/{application}/{tool}/{mode}/{norm}/{duplicates}/{contrast}/{ratio}/{condition}.bw"
    shell:
        """
            bigwigCompare --bigwig1 {input.chip}\
                          --bigwig2 {input.input}\
                          --ratio {wildcards.ratio}\
                          --numberOfProcessors {threads} \
                          --outFileFormat bigwig\
                          --outFileName {output}
        """
