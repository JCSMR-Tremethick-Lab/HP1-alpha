#!/bin/bash

# bigwigCompare commands
source activate py27

export source_dir="/home/sebastian/Data/Tremethick/HP1-alpha/ChIP-Seq/merged/processed_data/GRCh37_hg19_ensembl75/deepTools/bigwigCompare/normal/RPKM/duplicates_removed/ChIP-Input/log2"
bigwigCompare --bigwig1 $source_dir/MCF10A_shHP1b_HP1a.bw\
              --bigwig2 $source_dir/MCF10A_WT_HP1a.bw\
              --numberOfProcessors 32\
              --region chr10\
              --binSize 10\
              --skipNonCoveredRegions\
              --ratio subtract\
              --outFileName MCF10A_HP1a_WT-shHP1b.bw

bigwigCompare --bigwig1 $source_dir/MCF10A_shHP1a_HP1b.bw\
              --bigwig2 $source_dir/MCF10A_WT_HP1b.bw\
              --numberOfProcessors 32\
              --region chr10\
              --binSize 10\
              --skipNonCoveredRegions\
              --ratio log2\
              --outFileName MCF10A_HP1b_WT-shHP1a.bw

bigwigCompare --bigwig1 $source_dir/MCF10A_shHP1b_HP1a.bw\
              --bigwig2 $source_dir/MCF10A_WT_HP1a.bw\
              --numberOfProcessors 32\
              --region chr10\
              --binSize 10\
              --skipNonCoveredRegions\
              --ratio log2\
              --outFileName MCF10A_HP1a_WT-shHP1b.bw
