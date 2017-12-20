cd /home/sebastian/Data/Tremethick/HP1-alpha/Figures
export deepTools_dir="/home/sebastian/miniconda3/envs/py27/bin"
export bigWig_dir="/home/sebastian/Data/Tremethick/HP1-alpha/ChIP-Seq/merged/processed_data/GRCh37_hg19_ensembl75/deepTools/bigwigCompare/normal/RPKM/duplicates_removed/ChIP-Input/log2"
export bed_dir="/home/sebastian/Data/Tremethick/HP1-alpha/AnnotationData/GRCh37_hg19_ensembl75"
${deepTools_dir}/computeMatrix scale-regions --regionsFileName ${bed_dir}allGenes.bed\
                                             --scoreFileName ${bigWig_dir}/MCF10A_WT_HP1a.bw\
                                                             ${bigWig_dir}/MCF10A_WT_HP1b.bw\
                                                             ${bigWig_dir}/MCF10A_WT_H2AZ.bw\
                                                             ${bigWig_dir}/MCF10A_shHP1a_HP1b.bw\
                                                             ${bigWig_dir}/MCF10A_shHP1b_HP1a.bw\
                                                             ${bigWig_dir}/MCF10A_shH2AZ_HP1a.bw\
                                                             ${bigWig_dir}/MCF10A_shH2AZ_HP1b.bw\
                                                             --missingDataAsZero\
                                                             --skipZeros\
                                                             --binSize 10\
                                                             --numberOfProcessors 32\
                                                             --upstream=3000 --downstream=1500 --regionBodyLength=5000 --unscaled5prime=200 --unscaled3prime=200\
                                                             --outFileName bigwigCompare.allGenes.normal.RPKM.matrix.gz 1>>scale-regions.log 2>>scale-regions.log &

${deepTools_dir}/computeMatrix scale-regions --regionsFileName ${bed_dir}hg19_LTR_repeats.bed\
                                             --scoreFileName ${bigWig_dir}/MCF10A_WT_HP1a.bw\
                                                             ${bigWig_dir}/MCF10A_WT_HP1b.bw\
                                                             ${bigWig_dir}/MCF10A_WT_H2AZ.bw\
                                                             ${bigWig_dir}/MCF10A_shHP1a_HP1b.bw\
                                                             ${bigWig_dir}/MCF10A_shHP1b_HP1a.bw\
                                                             ${bigWig_dir}/MCF10A_shH2AZ_HP1a.bw\
                                                             ${bigWig_dir}/MCF10A_shH2AZ_HP1b.bw\
                                                             --missingDataAsZero\
                                                             --skipZeros\
                                                             --binSize 10\
                                                             --numberOfProcessors 16\
                                                             --upstream=500 --downstream=500 --regionBodyLength=280 --unscaled5prime=30 --unscaled3prime=30\
                                                             --outFileName bigwigCompare.LTR_repeats.normal.RPKM.matrix.gz 1>>scale-regions.log 2>>scale-regions.log &

${deepTools_dir}/computeMatrix scale-regions --regionsFileName ${bed_dir}hg19_GSAT_repeats.bed\
                                             --scoreFileName ${bigWig_dir}/MCF10A_WT_HP1a.bw\
                                                             ${bigWig_dir}/MCF10A_WT_HP1b.bw\
                                                             ${bigWig_dir}/MCF10A_WT_H2AZ.bw\
                                                             ${bigWig_dir}/MCF10A_shHP1a_HP1b.bw\
                                                             ${bigWig_dir}/MCF10A_shHP1b_HP1a.bw\
                                                             ${bigWig_dir}/MCF10A_shH2AZ_HP1a.bw\
                                                             ${bigWig_dir}/MCF10A_shH2AZ_HP1b.bw\
                                                             --missingDataAsZero\
                                                             --skipZeros\
                                                             --binSize 10\
                                                             --numberOfProcessors 16\
                                                             --upstream=500 --downstream=500 --regionBodyLength=900 --unscaled5prime=90 --unscaled3prime=90\
                                                             --outFileName bigwigCompare.GSAT_repeats.normal.RPKM.matrix.gz 1>>scale-regions.log 2>>scale-regions.log &

for i in hg19_SINE_repeats.bed hg19_LINE_repeats.bed hg19_Low_complexity_repeats.bed;
do
  if [ $i = "hg19_Low_complexity_repeats.bed" ]
  then
    export rbl=100
    export unscaled=10
  else
    export rbl=300
    export unscaled=30
  fi
  export name=$(echo $i|cut -f 1 -d ".")
  ${deepTools_dir}/computeMatrix scale-regions --regionsFileName ${bed_dir}/${i}\
                                               --scoreFileName ${bigWig_dir}/MCF10A_WT_HP1a.bw\
                                                               ${bigWig_dir}/MCF10A_WT_HP1b.bw\
                                                               ${bigWig_dir}/MCF10A_WT_H2AZ.bw\
                                                               ${bigWig_dir}/MCF10A_shHP1a_HP1b.bw\
                                                               ${bigWig_dir}/MCF10A_shHP1b_HP1a.bw\
                                                               ${bigWig_dir}/MCF10A_shH2AZ_HP1a.bw\
                                                               ${bigWig_dir}/MCF10A_shH2AZ_HP1b.bw\
                                                               --missingDataAsZero\
                                                               --skipZeros\
                                                               --binSize 10\
                                                               --numberOfProcessors 32\
                                                               --upstream=500 --downstream=500 --regionBodyLength=${rbl} --unscaled5prime=${unscaled} --unscaled3prime=${unscaled}\
                                                               --outFileName bigwigCompare.${name}.normal.RPKM.matrix.gz 1>>scale-regions.log 2>>scale-regions.log

done


${deepTools_dir}/plotHeatmap --matrixFile bigwigCompare.GSAT_repeats.normal.RPKM.matrix.gz\
                             --outFileName bigwigCompare.GSAT_repeats.normal.RPKM.pdf\
                             --outFileNameMatrix bigwigCompare.GSAT_repeats.normal.RPKM.tab\

${deepTools_dir}/plotHeatmap --matrixFile bigwigCompare.allGenes.normal.RPKM.matrix.gz\
                             --outFileName bigwigCompare.allGenes.normal.RPKM.pdf\
                             --outFileNameMatrix bigwigCompare.allGenes.normal.RPKM.tab\

${deepTools_dir}/plotHeatmap --matrixFile bigwigCompare.LTR_repeats.normal.RPKM.matrix.gz\
                             --outFileName bigwigCompare.LTR_repeats.normal.RPKM.pdf\
                             --outFileNameMatrix bigwigCompare.LTR_repeats.normal.RPKM.tab\
