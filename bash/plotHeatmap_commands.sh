/home/sebastian/miniconda3/envs/py27/bin/plotProfile --matrixFile ChIP-Seq/merged/processed_data/GRCh37_hg19_ensembl75/deepTools/computeMatrix/scale-regions/duplicates_removed/TSS/coding_normal_RPKM.matrix.gz\
                                                     --outFileName ChIP-Seq/merged/processed_data/GRCh37_hg19_ensembl75/deepTools/plotProfile/scale-regions/duplicates_removed/TSS/test.pdf\
                                                     --plotType se\
                                                     --perGroup\
                                                     --numPlotsPerRow 3


/home/sebastian/miniconda3/envs/py27/bin/plotProfile --matrixFile ChIP-Seq/merged/processed_data/GRCh37_hg19_ensembl75/deepTools/computeMatrix/scale-regions/duplicates_removed/TSS/ChIP-Input/log2/bigwigCompare.coding.normal.RPKM.matrix.gz\
                                                     --outFileName ChIP-Seq/merged/processed_data/GRCh37_hg19_ensembl75/deepTools/plotProfile/scale-regions/duplicates_removed/TSS/ChIP-Input/log2/se.bigwigCompare.coding.normal.RPKM.pdf\
                                                     --plotType se\
                                                     --perGroup\
                                                     --numPlotsPerRow 3


/home/sebastian/miniconda3/envs/py27/bin/computeMatrix scale-regions --regionsFileName /home/sebastian/Data/References/Annotations/Homo_sapiens/hg19/Ensembl/RepeatMasker.bed\
                                                                     --scoreFileName ChIP-Seq/merged/processed_data/GRCh37_hg19_ensembl75/deepTools/bigwigCompare/normal/RPKM/duplicates_removed/ChIP-Input/log2/MCF10A_WT_HP1a.bw\
                                                                                     ChIP-Seq/merged/processed_data/GRCh37_hg19_ensembl75/deepTools/bigwigCompare/normal/RPKM/duplicates_removed/ChIP-Input/log2/MCF10A_WT_HP1b.bw\
                                                                                     ChIP-Seq/merged/processed_data/GRCh37_hg19_ensembl75/deepTools/bigwigCompare/normal/RPKM/duplicates_removed/ChIP-Input/log2/MCF10A_WT_H2AZ.bw\
                                                                                     ChIP-Seq/merged/processed_data/GRCh37_hg19_ensembl75/deepTools/bigwigCompare/normal/RPKM/duplicates_removed/ChIP-Input/log2/MCF10A_shHP1a_HP1b.bw\
                                                                                     ChIP-Seq/merged/processed_data/GRCh37_hg19_ensembl75/deepTools/bigwigCompare/normal/RPKM/duplicates_removed/ChIP-Input/log2/MCF10A_shHP1b_HP1a.bw\
                                                                                     ChIP-Seq/merged/processed_data/GRCh37_hg19_ensembl75/deepTools/bigwigCompare/normal/RPKM/duplicates_removed/ChIP-Input/log2/MCF10A_shH2AZ_HP1a.bw\
                                                                                     ChIP-Seq/merged/processed_data/GRCh37_hg19_ensembl75/deepTools/bigwigCompare/normal/RPKM/duplicates_removed/ChIP-Input/log2/MCF10A_shH2AZ_HP1b.bw\
                                                                                     --missingDataAsZero\
                                                                                     --skipZeros\
                                                                                     --binSize 10\
                                                                                     --numberOfProcessors 32\
                                                                                     --upstream=100 --downstream=100 --regionBodyLength=200 --unscaled5prime=0 --unscaled3prime=0\
                                                                                     --outFileName ChIP-Seq/merged/processed_data/GRCh37_hg19_ensembl75/deepTools/computeMatrix/scale-regions/duplicates_removed/TSS/ChIP-Input/log2/bigwigCompare.repeats.normal.RPKM.matrix.gz

/home/sebastian/miniconda3/envs/py27/bin/plotProfile --matrixFile ChIP-Seq/merged/processed_data/GRCh37_hg19_ensembl75/deepTools/computeMatrix/scale-regions/duplicates_removed/TSS/ChIP-Input/log2/bigwigCompare.repeats.normal.RPKM.matrix.gz\
                                                     --outFileName ChIP-Seq/merged/processed_data/GRCh37_hg19_ensembl75/deepTools/plotProfile/scale-regions/duplicates_removed/TSS/ChIP-Input/log2/se.bigwigCompare.repeats.normal.RPKM.pdf\
                                                     --outFileNameData ChIP-Seq/merged/processed_data/GRCh37_hg19_ensembl75/deepTools/plotProfile/scale-regions/duplicates_removed/TSS/ChIP-Input/log2/se.bigwigCompare.repeats.normal.RPKM.data\
                                                     --outFileSortedRegions ChIP-Seq/merged/processed_data/GRCh37_hg19_ensembl75/deepTools/plotProfile/scale-regions/duplicates_removed/TSS/ChIP-Input/log2/se.bigwigCompare.repeats.normal.RPKM.bed\
                                                     --plotType se

/home/sebastian/miniconda3/envs/py27/bin/plotHeatmap --matrixFile ChIP-Seq/merged/processed_data/GRCh37_hg19_ensembl75/deepTools/computeMatrix/scale-regions/duplicates_removed/TSS/ChIP-Input/log2/bigwigCompare.repeats.normal.RPKM.matrix.gz\
                                                     --outFileName ChIP-Seq/merged/processed_data/GRCh37_hg19_ensembl75/deepTools/plotHeatmap/scale-regions/duplicates_removed/TSS/ChIP-Input/log2/bigwigCompare.repeats.normal.RPKM.pdf\
                                                     --outFileNameMatrix ChIP-Seq/merged/processed_data/GRCh37_hg19_ensembl75/deepTools/plotHeatmap/scale-regions/duplicates_removed/TSS/ChIP-Input/log2/bigwigCompare.repeats.normal.RPKM.tab\
                                                     
/home/sebastian/miniconda3/envs/py27/bin/plotHeatmap --matrixFile ChIP-Seq/merged/processed_data/GRCh37_hg19_ensembl75/deepTools/computeMatrix/scale-regions/duplicates_removed/TSS/ChIP-Input/log2/bigwigCompare.coding.normal.RPKM.matrix.gz\
                                                     --outFileName ChIP-Seq/merged/processed_data/GRCh37_hg19_ensembl75/deepTools/plotHeatmap/scale-regions/duplicates_removed/TSS/ChIP-Input/log2/bigwigCompare.coding.normal.RPKM.pdf\
                                                     --outFileNameMatrix ChIP-Seq/merged/processed_data/GRCh37_hg19_ensembl75/deepTools/plotHeatmap/scale-regions/duplicates_removed/TSS/ChIP-Input/log2/bigwigCompare.coding.normal.RPKM.tab\
                                                     

# only plot regions of conditionMCF10A_shHP1a
/home/sebastian/miniconda3/envs/py27/bin/computeMatrix scale-regions --regionsFileName /home/sebastian/Data/References/Annotations/Homo_sapiens/hg19/Ensembl/RepeatMasker.bed\
                                                                     --scoreFileName ChIP-Seq/merged/processed_data/GRCh37_hg19_ensembl75/deepTools/bigwigCompare/normal/RPKM/duplicates_removed/ChIP-Input/log2/MCF10A_WT_HP1a.bw\
                                                                                     ChIP-Seq/merged/processed_data/GRCh37_hg19_ensembl75/deepTools/bigwigCompare/normal/RPKM/duplicates_removed/ChIP-Input/log2/MCF10A_WT_HP1b.bw\
                                                                                     ChIP-Seq/merged/processed_data/GRCh37_hg19_ensembl75/deepTools/bigwigCompare/normal/RPKM/duplicates_removed/ChIP-Input/log2/MCF10A_WT_H2AZ.bw\
                                                                                     ChIP-Seq/merged/processed_data/GRCh37_hg19_ensembl75/deepTools/bigwigCompare/normal/RPKM/duplicates_removed/ChIP-Input/log2/MCF10A_shHP1a_HP1b.bw\
                                                                                     ChIP-Seq/merged/processed_data/GRCh37_hg19_ensembl75/deepTools/bigwigCompare/normal/RPKM/duplicates_removed/ChIP-Input/log2/MCF10A_shHP1b_HP1a.bw\
                                                                                     ChIP-Seq/merged/processed_data/GRCh37_hg19_ensembl75/deepTools/bigwigCompare/normal/RPKM/duplicates_removed/ChIP-Input/log2/MCF10A_shH2AZ_HP1a.bw\
                                                                                     ChIP-Seq/merged/processed_data/GRCh37_hg19_ensembl75/deepTools/bigwigCompare/normal/RPKM/duplicates_removed/ChIP-Input/log2/MCF10A_shH2AZ_HP1b.bw\
                                                                                     --missingDataAsZero\
                                                                                     --skipZeros\
                                                                                     --binSize 10\
                                                                                     --numberOfProcessors 32\
                                                                                     --upstream=100 --downstream=100 --regionBodyLength=200 --unscaled5prime=0 --unscaled3prime=0\
                                                                                     --outFileName ChIP-Seq/merged/processed_data/GRCh37_hg19_ensembl75/deepTools/computeMatrix/scale-regions/duplicates_removed/TSS/ChIP-Input/log2/bigwigCompare.repeats.normal.RPKM.matrix.gz
                                                                                                          