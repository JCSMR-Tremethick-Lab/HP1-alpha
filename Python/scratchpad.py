import json
import os
from pprint import pprint
from os.path import join
from snakemake.io import expand, glob_wildcards

with open("config_ChIP-Seq.json") as data_file:
    config = json.load(data_file)

REF_GENOME = "hg19"
REF_VERSION = config["references"][REF_GENOME]["version"][1]


glob_wildcards("ChIP-Seq/merged/{outdir}/{reference_version}/{application}/bigwigCompare/{mode}/{norm}/{duplicates}/{contrast}/{ratio}/{condition}.bw")
wildcards = dict()
wildcards = {"assayID" : "RNA-Seq", "processed_dir" : "processed_data", "runID": "NB501086_0114_B_Azad_JCSMR_hRNAseq", "reports_dir" : "reports", "reference_version" : "GRCh38_ensembl84_ERCC"}

wildcards = {"assayID" : "ChIP-Seq",
             "processed_dir" : "processed_data",
             "runID": "merged",
             "reports_dir" : "reports",
             "reference_version" : "GRCh37_hg19_ensembl75",
             "unit" : "MCF10A_shHP1b_Input_1b",
             "condition" : "MCF10A_shH2AZ",
             "duplicates" : "duplicates_removed",
             "mode" : "normal",
             "norm" : "RPKM",
             "application" : "deepTools",
             "source" : "bigwigCompare",
             "contrast" : "ChIP-Input",
             "ratio" : "log2",
             "macs2_command" : "callpeak"
             }

def getRegionFiles(wildcards):
    fn = []
    for i in config["program_parameters"]["deepTools"]["regionFiles"][wildcards["reference_version"]][wildcards["region"]].values():
        fn.append(home + i)
    return(fn)

def getComputeMatrixInputMerged(wildcards):
    fn = []
    path = "/".join([wildcards["assayID"],
                     wildcards["runID"],
                     config["processed_dir"],
                     REF_VERSION,
                     wildcards["application"]]) + "/"
    if wildcards["source"] == "bamCoverage":
        for i in config["samples"]["ChIP-Seq"]["replicates"].keys():
            fn.append(path + "/".join([wildcards["source"],
                                wildcards["mode"],
                                wildcards["norm"],
                                wildcards["duplicates"],
                                i + ".bw"]))
    elif wildcards["source"] == "bigwigCompare":
        for i in config["samples"]["ChIP-Seq"][wildcards["contrast"]].keys():
            fn.append(path + "/".join([wildcards["source"],
                                wildcards["mode"],
                                wildcards["norm"],
                                wildcards["duplicates"],
                                wildcards["contrast"],
                                wildcards["ratio"],
                                i + ".bw"]))
    return(fn)


def macs2OutputFiles():
    fn = []
    for i in config["samples"]["ChIP-Seq"]["ChIP-Input"].keys():
        for j in config["samples"]["ChIP-Seq"]["replicates"][config["samples"]["ChIP-Seq"]["ChIP-Input"][i]["ChIP"]]:
            contrasts = i + "/" + j
            path = "/".join(["ChIP-Seq",
                            "processed_data",
                            REF_VERSION,
                            "macs2",
                            contrasts,
                            "callpeak"])
            pseudorep1 = i + "/" + j + "_pseudo_rep1"
            pseudorep2 = i + "/" + j + "_pseudo_rep2"
            pseudorep1 = "/".join(["ChIP-Seq",
                            "processed_data",
                            REF_VERSION,
                            "macs2",
                            pseudorep1,
                            "callpeak"])
            pseudorep2 = "/".join(["ChIP-Seq",
                            "processed_data",
                            REF_VERSION,
                            "macs2",
                            pseudorep2,
                            "callpeak"])
            fn.append(path)
            fn.append(pseudorep1)
            fn.append(pseudorep2)
    return(fn)


MACS2_output = macs2OutputFiles(wildcards)

def getMergedInputBAM(wildcards):
    fn = []
    inp = config["samples"]["ChIP-Seq"]["ChIP-Input"][wildcards["contrast"]]["Input"]
    inp = config["samples"]["ChIP-Seq"]["replicates"][inp]
    for i in config["samples"]["ChIP-Seq"]["runID"]:
        for j in config["samples"]["ChIP-Seq"][i].keys():
            for k in inp:
                if k == j:
                    f = "/".join([wildcards["assayID"],
                                  i,
                                  config["processed_dir"],
                                  REF_VERSION,
                                  "bowtie2",
                                  "duplicates_removed",
                                  j + ".Q" + config["alignment_quality"] + ".sorted.bam"
                                  ])
                    fn.append(f)
    return(fn)

def getChIPBam(wildcards):
    fn = []
    for i in config["samples"]["ChIP-Seq"]["runID"]:
        for j in config["samples"]["ChIP-Seq"][i].keys():
                if j == wildcards["unit"]:
                        f = "/".join([wildcards["assayID"],
                                      i,
                                      config["processed_dir"],
                                      REF_VERSION,
                                      "bowtie2",
                                      "duplicates_removed",
                                      j + ".Q" + config["alignment_quality"] + ".sorted.bam"
                                      ])
                        fn.append(f)
    return(fn)

def getBAMbyCondition(wildcards):
    fn = []
    for i in config["samples"]["ChIP-Seq"]["runID"]:
        if i in config["samples"]["ChIP-Seq"]["conditions"][wildcards["condition"]].keys():
            for j in config["samples"]["ChIP-Seq"]["conditions"][wildcards["condition"]][i]:
                fn.append(join("ChIP-Seq/" + i + "/" + config["processed_dir"] + "/" + REF_VERSION + "/bowtie2/" + wildcards["duplicates"] + "/" + j + ".Q" + config["alignment_quality"] + ".sorted.bam"))
    return(fn)


def getAllFASTQ(wildcards):
    fn =[]
    for i in config["samples"][wildcards["assayID"]][wildcards["runID"]]:
        for j in config["samples"][wildcards["assayID"]][wildcards["runID"]][i]:
            fn.append("/".join([wildcards["assayID"], wildcards["runID"], config["raw_dir"], j]))
    return(fn)


for i in config["samples"]["ChIP-Seq"]["runID"]:
    for j in config["samples"]["ChIP-Seq"][i]:
        join("ChIP-Seq/" + i + "/" + config["processed_dir"] + "/" + config["trim_dir"] + "/" + j + "_R1.fastq.gz")
        join("ChIP-Seq/" + i + "/" + config["processed_dir"] + "/" + config["trim_dir"] + "/" + j + "_R2.fastq.gz")

for i in config["samples"]["ChIP-Seq"]["runID"]:
    for j in config["samples"]["ChIP-Seq"][i]:
        join("ChIP-Seq/" + i + "/" + config["processed_dir"] + "/" + REF_VERSION + "/bowtie2/duplicates_removed/" + j + ".Q" + config["alignment_quality"] + ".sorted.bam")

for i in config["samples"]["ChIP-Seq"]["runID"]:
    for j in config["samples"]["ChIP-Seq"][i].values():
        for k in j:
            print(k)


config["references"]["hg38"]["STAR"][wildcards["reference_version"]]
