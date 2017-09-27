import json
from pprint import pprint
from os.path import join

with open("config_ChIP-Seq.json") as data_file:
    config = json.load(data_file)

REF_GENOME = "hg19"
REF_VERSION = config["references"][REF_GENOME]["version"][1]

wildcards = dict()
wildcards = {"assayID" : "RNA-Seq", "processed_dir" : "processed_data", "runID": "NB501086_0114_B_Azad_JCSMR_hRNAseq", "reports_dir" : "reports", "reference_version" : "GRCh38_ensembl84_ERCC"}

wildcards = {"assayID" : "ChIP-Seq",
             "processed_dir" : "processed_data",
             "runID": "NB501086_0138_TSoboleva_JCSMR_Mouse_ChIPseq",
             "reports_dir" : "reports",
             "reference_version" : "GRCh37_hg19_ensembl75",
             "unit" : "MCF10A_shHP1b_Input_1b",
             "condition" : "MCF10A_shH2AZ",
             "duplicates" : "duplicates_removed",
             "mode" : "normal",
             "norm" : "RPKM",
             "application" : "deepTools"}

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
