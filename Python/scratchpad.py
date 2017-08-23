import json
from pprint import pprint
from os.path import join

with open("config_ChIP-Seq.json") as data_file:
    config = json.load(data_file)

wildcards = dict()
wildcards = {"assayID" : "RNA-Seq", "processed_dir" : "processed_data", "runID": "NB501086_0114_B_Azad_JCSMR_hRNAseq", "reports_dir" : "reports", "reference_version" : "GRCh38_ensembl84_ERCC"}

wildcards = {"assayID" : "ChIP-Seq",
             "processed_dir" : "processed_data",
             "runID": "NB501086_0136_TSoboleva_JCSMR_Mouse_ChIPseq",
             "reports_dir" : "reports",
             "reference_version" : "GRCh37_hg19_ensembl75",
             "unit" : "MCF10A_shHP1b_Input_1b"}

def getAllFASTQ(wildcards):
    fn =[]
    for i in config["samples"][wildcards["assayID"]][wildcards["runID"]]:
        for j in config["samples"][wildcards["assayID"]][wildcards["runID"]][i]:
            fn.append("ChIP-Seq/NB501086_0136_TSoboleva_JCSMR_Mouse_ChIPseq/fastq/" + j)
    return(fn)


for i in config["samples"]["ChIP-Seq"]["runID"]:
    for j in config["samples"]["ChIP-Seq"][i]:
        join("ChIP-Seq/" + i + "/" + config["processed_dir"] + "/" + config["trim_dir"] + "/" + j + "_R1.fastq.gz")
        join("ChIP-Seq/" + i + "/" + config["processed_dir"] + "/" + config["trim_dir"] + "/" + j + "_R2.fastq.gz")


config["references"]["hg38"]["STAR"][wildcards["reference_version"]]
