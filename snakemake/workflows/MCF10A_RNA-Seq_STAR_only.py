_author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2015-04-22"
__modified__ = "2017-05-11"

from snakemake.exceptions import MissingInputException
import os

rule:
    version: 0.1

localrules:
    all, run_STAR

# variables
REF_GENOME="hg19"
REF_VERSION="GRCh37_hg19_ensembl75"
HOME=os.environ['HOME']

wrapper_dir = os.environ['HOME'] + "/Development/snakemake-wrappers/bio"

include_prefix= os.environ['HOME'] + "/Development/JCSMR-Tremethick-Lab/HP1-alpha/snakemake/rules/"

assay_dir = "NB501086_0114_B_Azad_JCSMR_hRNAseq"

include:
    include_prefix + "run_STAR.py"

# targets
STARS = expand("{assayID}/{file}",
              assayID = "RNA-Seq",
              file = [ i + "/" + config["processed_dir"] + "/" + REF_VERSION + "/STAR/full/" + j + ".bam.bai"\
                  for i in config["samples"]["RNA-Seq"]["runID"] \
                      for j in config["samples"]["RNA-Seq"][i]])

rule all:
    input:
        STARS
