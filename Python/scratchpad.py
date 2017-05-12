import json
from pprint import pprint

with open("config.json") as data_file:
    config = json.load(data_file)

wildcards = dict()
wildcards = {"assayID" : "RNA-Seq", "processed_dir" : "processed_data", "runID": "NB501086_0114_B_Azad_JCSMR_hRNAseq", "reports_dir" : "reports", "reference_version" : "GRCh38_ensembl84_ERCC"}

def getAllFASTQ(wildcards):
    fn =[]
    for i in config["samples"][wildcards["assayID"]][wildcards["runID"]]:
        for j in config["samples"][wildcards["assayID"]][wildcards["runID"]][i]:
            fn.append("RNA-Seq/NB501086_0114_B_Azad_JCSMR_hRNAseq/fastq/" + j)
    return(fn)

config["references"]["hg38"]["STAR"][wildcards["reference_version"]]
