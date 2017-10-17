from snakemake.exceptions import MissingInputException
import os

_author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2015-04-22"

rule:
    version: 0.3

localrules:
    all, run_kallisto, run_STAR, run_htseq, run_cutadapt

wrapper_dir = os.environ['HOME'] + "/Development/snakemake-wrappers/bio"

include_prefix= os.environ['HOME'] + "/Development/JCSMR-Tremethick-Lab/Breast/snakemake/rules/"

include:
    include_prefix + "run_kallisto.py"

# variables
REF_GENOME="hg19"
REF_VERSION="GRCh37_hg19_ensembl75"

# targets
BEARS = expand("{assayID}/{file}",
              assayID = "RNA-Seq",
              file = [ i + "/" + config["processed_dir"] + "/" + REF_VERSION + "/kallisto/" + j \
                  for i in config["samples"]["RNA-Seq"]["runID"] \
                      for j in config["samples"]["RNA-Seq"][i]])

rule run_kallisto:
    input:
        BEARS
