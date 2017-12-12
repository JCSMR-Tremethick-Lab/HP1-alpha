# testing peak model parameters on H2AZ WT data
dataDir="/home/sebastian/Data/Tremethick/HP1-alpha/ChIP-Seq/merged/processed_data/GRCh37_hg19_ensembl75/duplicates_removed"
export dataDir

macs2 predictd -i $dataDir/MCF10A_WT_H2AZ_ChIP.bam\
               --outdir predict_output\
               --rfile H2AZ_peaks.R\
               --gsize hs

macs2 callpeak -t $dataDir/MCF10A_WT_H2AZ_ChIP.bam\
               -c $dataDir/MCF10A_WT_Input.bam\
               --gsize hs\
               --format BAMPE\
               --outdir callpeak_BAMPE\
               -n MCF10A_WT_H2AZ_pooled 1>>callpeak.log 2>>callpeak.log &
