#!/bin/bash

cd /home/sebastian/Data/Tremethick/HP1-alpha/Figures
export deepTools_dir="/home/sebastian/miniconda3/envs/py27/bin"
export bigWig_dir="/home/sebastian/Data/Tremethick/HP1-alpha/ChIP-Seq/merged/processed_data/GRCh37_hg19_ensembl75/deepTools/bigwigCompare/normal/RPKM/duplicates_removed/ChIP-Input/log2"
export bed_dir="/home/sebastian/Data/Tremethick/HP1-alpha/AnnotationData/GRCh37_hg19_ensembl75/"

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
                                                             --upstream=3000 --downstream=3000 --regionBodyLength=5000 --unscaled5prime=200 --unscaled3prime=200\
                                                             --outFileName bigwigCompare.allGenes.normal.RPKM.matrix.gz 1>>scale-regions.log 2>>scale-regions.log &

${deepTools_dir}/computeMatrix scale-regions --regionsFileName ${bed_dir}/allGenesRNASeq.bed\
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
                                                             --upstream=3000 --downstream=3000 --regionBodyLength=5000 --unscaled5prime=0 --unscaled3prime=0\
                                                             --outFileName bigwigCompare.allGenesRNASeq.bed.normal.RPKM.matrix.gz 1>>scale-regions.log 2>>scale-regions.log &


#-------------------------------------------------------------------------------
# exons
#-------------------------------------------------------------------------------

#------------------scaled exons
${deepTools_dir}/computeMatrix scale-regions --regionsFileName ${bed_dir}/allExons.bed\
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
                                                             --upstream=1000 --downstream=1000 --regionBodyLength=300 --unscaled5prime=200 --unscaled3prime=200\
                                                             --outFileName bigwigCompare.allExons.normal.RPKM.matrix.gz 1>>scale-regions.log 2>>scale-regions.log &

${deepTools_dir}/plotHeatmap --matrixFile bigwigCompare.allExons.normal.RPKM.matrix.gz\
                             --outFileName bigwigCompare.allExons.normal.RPKM.pdf\
                             --outFileNameMatrix bigwigCompare.allExons.normal.RPKM.tab


${deepTools_dir}/computeMatrix scale-regions --regionsFileName ${bed_dir}/expressedExons.bed ${bed_dir}/silentExons.bed\
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
                                                             --upstream=500 --downstream=500 --regionBodyLength=300 --unscaled5prime=50 --unscaled3prime=50\
                                                             --outFileName bigwigCompare.exons.scaled.RPKM.matrix.gz 1>>scale-regions.log 2>>scale-regions.log

${deepTools_dir}/plotHeatmap --matrixFile bigwigCompare.exons.scaled.RPKM.matrix.gz\
                             --outFileName bigwigCompare.exons.scaled.RPKM.pdf\
                             --outFileNameMatrix bigwigCompare.exons.scaled.RPKM.tab

#------------------reference point
${deepTools_dir}/computeMatrix reference-point --regionsFileName ${bed_dir}/allExons.bed\
                                               --referencePoint TSS\
                                               --beforeRegionStartLength 1000\
                                               --afterRegionStartLength 1000\
                                               --scoreFileName ${bigWig_dir}/MCF10A_WT_HP1a.bw\
                                                               ${bigWig_dir}/MCF10A_WT_HP1b.bw\
                                                               ${bigWig_dir}/MCF10A_WT_H2AZ.bw\
                                                               ${bigWig_dir}/MCF10A_shHP1a_HP1b.bw\
                                                               ${bigWig_dir}/MCF10A_shHP1b_HP1a.bw\
                                                               ${bigWig_dir}/MCF10A_shH2AZ_HP1a.bw\
                                                               ${bigWig_dir}/MCF10A_shH2AZ_HP1b.bw\
                                               --samplesLabel "WT - HP1a ChIP"\
                                                              "WT - HP1b ChIP"\
                                                              "WT - H2A.Z ChIP"\
                                                              "HP1a KD - HP1b ChIP"\
                                                              "HP1b KD - HP1a ChIP"\
                                                              "H2A.Z KD - HP1a ChIP"\
                                                              "H2A.Z KD - HP1b ChIP"\
                                               --binSize 10\
                                               --numberOfProcessors 32\
                                               --outFileName bigwigCompare.intronExon.refPoint.RPKM.matrix.gz 1>>scale-regions.log 2>>scale-regions.log &

${deepTools_dir}/plotHeatmap --matrixFile bigwigCompare.intronExon.refPoint.RPKM.matrix.gz\
                             --outFileName bigwigCompare.intronExon.refPoint.RPKM.pdf\
                             --outFileNameMatrix bigwigCompare.intronExon.refPoint.normal.RPKM.tab\
                             --refPointLabel "intron/exon"

${deepTools_dir}/computeMatrix reference-point --regionsFileName ${bed_dir}/expressedExons.bed ${bed_dir}/silentExons.bed\
                                               --referencePoint TSS\
                                               --beforeRegionStartLength 1000\
                                               --afterRegionStartLength 1000\
                                               --scoreFileName ${bigWig_dir}/MCF10A_WT_HP1a.bw\
                                                               ${bigWig_dir}/MCF10A_WT_HP1b.bw\
                                                               ${bigWig_dir}/MCF10A_WT_H2AZ.bw\
                                                               ${bigWig_dir}/MCF10A_shHP1a_HP1b.bw\
                                                               ${bigWig_dir}/MCF10A_shHP1b_HP1a.bw\
                                                               ${bigWig_dir}/MCF10A_shH2AZ_HP1a.bw\
                                                               ${bigWig_dir}/MCF10A_shH2AZ_HP1b.bw\
                                               --samplesLabel "WT - HP1a ChIP"\
                                                              "WT - HP1b ChIP"\
                                                              "WT - H2A.Z ChIP"\
                                                              "HP1a KD - HP1b ChIP"\
                                                              "HP1b KD - HP1a ChIP"\
                                                              "H2A.Z KD - HP1a ChIP"\
                                                              "H2A.Z KD - HP1b ChIP"\
                                               --binSize 10\
                                               --numberOfProcessors 32\
                                               --outFileName bigwigCompare.intronExon.twoGroups.refPoint.RPKM.matrix.gz 1>>scale-regions.log 2>>scale-regions.log

${deepTools_dir}/plotHeatmap --matrixFile bigwigCompare.intronExon.twoGroups.refPoint.RPKM.matrix.gz\
                             --outFileName bigwigCompare.intronExon.twoGroups.refPoint.RPKM.pdf\
                             --outFileNameMatrix bigwigCompare.intronExon.twoGroups.refPoint.RPKM.tab\
                             --refPointLabel "intron/exon"


${deepTools_dir}/computeMatrix reference-point --regionsFileName ${bed_dir}/allExons.bed\
                                               --referencePoint TES\
                                               --beforeRegionStartLength 1000\
                                               --afterRegionStartLength 1000\
                                               --nanAfterEnd\
                                               --scoreFileName ${bigWig_dir}/MCF10A_WT_HP1a.bw\
                                                               ${bigWig_dir}/MCF10A_WT_HP1b.bw\
                                                               ${bigWig_dir}/MCF10A_WT_H2AZ.bw\
                                                               ${bigWig_dir}/MCF10A_shHP1a_HP1b.bw\
                                                               ${bigWig_dir}/MCF10A_shHP1b_HP1a.bw\
                                                               ${bigWig_dir}/MCF10A_shH2AZ_HP1a.bw\
                                                               ${bigWig_dir}/MCF10A_shH2AZ_HP1b.bw\
                                               --samplesLabel "WT - HP1a ChIP"\
                                                              "WT - HP1b ChIP"\
                                                              "WT - H2A.Z ChIP"\
                                                              "HP1a KD - HP1b ChIP"\
                                                              "HP1b KD - HP1a ChIP"\
                                                              "H2A.Z KD - HP1a ChIP"\
                                                              "H2A.Z KD - HP1b ChIP"\
                                               --missingDataAsZero\
                                               --skipZeros\
                                               --binSize 10\
                                               --numberOfProcessors 32\
                                               --outFileName bigwigCompare.exonIntron.refPoint.RPKM.matrix.gz 1>>scale-regions.log 2>>scale-regions.log &

${deepTools_dir}/plotHeatmap --matrixFile bigwigCompare.exonIntron.refPoint.RPKM.matrix.gz\
                             --outFileName bigwigCompare.exonIntron.refPoint.RPKM.pdf\
                             --outFileNameMatrix bigwigCompare.exonIntron.refPoint.normal.RPKM.tab\
                             --refPointLabel "exon/intron"

#-------------------------------------------------------------------------------
# repeats
#-------------------------------------------------------------------------------
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

${deepTools_dir}/plotHeatmap --matrixFile bigwigCompare.GSAT_repeats.normal.RPKM.matrix.gz\
                             --outFileName bigwigCompare.GSAT_repeats.normal.RPKM.pdf\
                             --outFileNameMatrix bigwigCompare.GSAT_repeats.normal.RPKM.tab\

${deepTools_dir}/computeMatrix scale-regions --regionsFileName ${bed_dir}hg19_SATs_repeats.bed\
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
                                                             --numberOfProcessors 30\
                                                             --upstream=500 --downstream=500 --regionBodyLength=900 --unscaled5prime=90 --unscaled3prime=90\
                                                             --outFileName bigwigCompare.SATs.normal.RPKM.matrix.gz 1>>scale-regions.log 2>>scale-regions.log

${deepTools_dir}/plotHeatmap --matrixFile bigwigCompare.SATs.normal.RPKM.matrix.gz \
                             --outFileName bigwigCompare.SATs.normal.RPKM.pdf\
                             --outFileNameMatrix bigwigCompare.SATs.normal.RPKM.tab\

${deepTools_dir}/plotHeatmap --kmeans 4\
                             --matrixFile bigwigCompare.SATs.normal.RPKM.matrix.gz \
                             --outFileName bigwigCompare.SATs.kmeans4.RPKM.pdf\
                             --outFileNameMatrix bigwigCompare.SATs.kmeans4.RPKM.tab\

${deepTools_dir}/plotProfile --kmeans 4\
                             --plotType se\
                             --matrixFile bigwigCompare.SATs.normal.RPKM.matrix.gz \
                             --outFileName bigwigCompare.SATs.kmeans4.profile.RPKM.pdf


# plot coverage across genes that change in expression quartiles
${deepTools_dir}/computeMatrix scale-regions --regionsFileName ${bed_dir}/Q2_to_Q3.bed\
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
                                                             --upstream=3000 --downstream=3000 --regionBodyLength=5000 --unscaled5prime=200 --unscaled3prime=200\
                                                             --outFileName bigwigCompare.Q2_to_Q3.normal.RPKM.matrix.gz 1>>scale-regions.log 2>>scale-regions.log &

${deepTools_dir}/plotHeatmap --matrixFile bigwigCompare.Q2_to_Q3.normal.RPKM.matrix.gz\
                             --outFileName bigwigCompare.Q2_to_Q3.normal.RPKM.pdf\
                             --outFileNameMatrix bigwigCompare.Q2_to_Q3.normal.RPKM.tab

${deepTools_dir}/computeMatrix scale-regions --regionsFileName ${bed_dir}/Q3_to_Q4.bed\
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
                                                             --upstream=3000 --downstream=3000 --regionBodyLength=5000 --unscaled5prime=200 --unscaled3prime=200\
                                                             --outFileName bigwigCompare.Q3_to_Q4.normal.RPKM.matrix.gz 1>>scale-regions.log 2>>scale-regions.log &

${deepTools_dir}/plotHeatmap --matrixFile bigwigCompare.Q3_to_Q4.normal.RPKM.matrix.gz\
                             --outFileName bigwigCompare.Q3_to_Q4.normal.RPKM.pdf\
                             --outFileNameMatrix bigwigCompare.Q3_to_Q4.normal.RPKM.tab

${deepTools_dir}/computeMatrix scale-regions --regionsFileName ${bed_dir}/Q2_to_Q4.bed\
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
                                                             --upstream=3000 --downstream=3000 --regionBodyLength=5000 --unscaled5prime=200 --unscaled3prime=200\
                                                             --outFileName bigwigCompare.Q2_to_Q4.normal.RPKM.matrix.gz 1>>scale-regions.log 2>>scale-regions.log &

${deepTools_dir}/plotHeatmap --matrixFile bigwigCompare.Q2_to_Q4.normal.RPKM.matrix.gz\
                             --outFileName bigwigCompare.Q2_to_Q4.normal.RPKM.pdf\
                             --outFileNameMatrix bigwigCompare.Q2_to_Q4.normal.RPKM.tab



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

${deepTools_dir}/plotHeatmap --matrixFile bigwigCompare.allGenes.normal.RPKM.matrix.gz\
                             --outFileName bigwigCompare.allGenes.normal.RPKM.pdf\
                             --outFileNameMatrix bigwigCompare.allGenes.normal.RPKM.tab\

${deepTools_dir}/plotHeatmap --matrixFile bigwigCompare.LTR_repeats.normal.RPKM.matrix.gz\
                             --outFileName bigwigCompare.LTR_repeats.normal.RPKM.pdf\
                             --outFileNameMatrix bigwigCompare.LTR_repeats.normal.RPKM.tab\

${deepTools_dir}/plotHeatmap --matrixFile bigwigCompare.allGenes.normal.RPKM.matrix.gz\
                             --outFileName bigwigCompare.allGenes.normal.RPKM.kmeans6.pdf\
                             --outFileNameMatrix bigwigCompare.allGenes.normal.RPKM.kmeans6.tab\
                              --outFileSortedRegions bigwigCompare.allGenes.normal.RPKM.kmeans6.bed\
                             --kmeans 6\

${deepTools_dir}/plotHeatmap --matrixFile bigwigCompare.allGenesRNASeq.bed.normal.RPKM.matrix.gz\
                             --outFileName bigwigCompare.allGenesRNASeq.normal.RPKM.kmeans6.pdf\
                             --outFileNameMatrix bigwigCompare.allGenesRNASeq.normal.RPKM.kmeans6.tab\
                              --outFileSortedRegions bigwigCompare.allGenesRNASeq.normal.RPKM.kmeans6.bed\
                             --kmeans 6


#--------------
# make individual compute matrices
${deepTools_dir}/computeMatrix scale-regions --regionsFileName ${bed_dir}allGenes.bed\
                                             --scoreFileName ${bigWig_dir}/MCF10A_WT_HP1a.bw\
                                                             --missingDataAsZero\
                                                             --skipZeros\
                                                             --binSize 10\
                                                             --numberOfProcessors 8\
                                                             --upstream=3000 --downstream=1500 --regionBodyLength=5000 --unscaled5prime=200 --unscaled3prime=200\
                                                             --outFileName MCF10A_WT_HP1a.bigwigCompare.allGenes.normal.RPKM.matrix.gz 1>>scale-regions.log 2>>scale-regions.log &

${deepTools_dir}/plotHeatmap --matrixFile MCF10A_WT_HP1a.bigwigCompare.allGenes.normal.RPKM.matrix.gz\
                             --outFileName MCF10A_WT_HP1a.bigwigCompare.allGenes.normal.RPKM.kmeans6.pdf\
                             --outFileNameMatrix MCF10A_WT_HP1a.bigwigCompare.allGenes.normal.RPKM.kmeans6.tab\
                             --outFileSortedRegions MCF10A_WT_HP1a.bigwigCompare.allGenes.normal.RPKM.kmeans6.bed\
                             --kmeans 6


${deepTools_dir}/computeMatrix scale-regions --regionsFileName ${bed_dir}allGenes.bed\
                                             --scoreFileName ${bigWig_dir}/MCF10A_WT_HP1b.bw\
                                                             --missingDataAsZero\
                                                             --skipZeros\
                                                             --binSize 10\
                                                             --numberOfProcessors 8\
                                                             --upstream=3000 --downstream=1500 --regionBodyLength=5000 --unscaled5prime=200 --unscaled3prime=200\
                                                             --outFileName MCF10A_WT_HP1b.bigwigCompare.allGenes.normal.RPKM.matrix.gz 1>>scale-regions.log 2>>scale-regions.log &
