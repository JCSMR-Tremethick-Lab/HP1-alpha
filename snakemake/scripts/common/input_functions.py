def getBAMbyReplicates(wildcards):
    """input function for collecting all BAM files per replicate"""
    fn=[]
    for i in config["samples"]["ChIP-Seq"]["runID"]:
        for j in config["samples"]["ChIP-Seq"][i]:
            if j in config["samples"]["ChIP-Seq"]["replicates"][wildcards["replicates"]]:
                fn.append(join("ChIP-Seq/" +
                               i +
                               "/" +
                               config["processed_dir"] +
                               "/" +
                               REF_VERSION +
                               "/bowtie2/" +
                               wildcards["duplicates"] +
                               "/" +
                               j +
                               ".Q" +
                               config["alignment_quality"] +
                               ".sorted.bam"))
    return(fn)
