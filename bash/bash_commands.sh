mv MCF10A_shH2AZ_HP1a_1_S13_R1_001.fastq.gz MCF10A_shH2AZHP1a_1_S13_R1_001.fastq.gz
mv MCF10A_shH2AZ_HP1a_1_S13_R2_001.fastq.gz MCF10A_shH2AZHP1a_1_S13_R2_001.fastq.gz
mv MCF10A_shH2AZ_HP1a_2_S14_R1_001.fastq.gz MCF10A_shH2AZHP1a_2_S14_R1_001.fastq.gz
mv MCF10A_shH2AZ_HP1a_2_S14_R2_001.fastq.gz MCF10A_shH2AZHP1a_2_S14_R2_001.fastq.gz
mv MCF10A_shH2AZ_HP1a_3_S15_R1_001.fastq.gz MCF10A_shH2AZHP1a_3_S15_R1_001.fastq.gz
mv MCF10A_shH2AZ_HP1a_3_S15_R2_001.fastq.gz MCF10A_shH2AZHP1a_3_S15_R2_001.fastq.gz
mv MCF10A_shHP1a_b_1_S10_R1_001.fastq.gz MCF10A_shHP1ab_1_S10_R1_001.fastq.gz
mv MCF10A_shHP1a_b_1_S10_R2_001.fastq.gz MCF10A_shHP1ab_1_S10_R2_001.fastq.gz
mv MCF10A_shHP1a_b_2_S11_R1_001.fastq.gz MCF10A_shHP1ab_2_S11_R1_001.fastq.gz
mv MCF10A_shHP1a_b_2_S11_R2_001.fastq.gz MCF10A_shHP1ab_2_S11_R2_001.fastq.gz
mv MCF10A_shHP1a_b_3_S12_R1_001.fastq.gz MCF10A_shHP1ab_3_S12_R1_001.fastq.gz
mv MCF10A_shHP1a_b_3_S12_R2_001.fastq.gz MCF10A_shHP1ab_3_S12_R2_001.fastq.gz

for i in $(ls | cut -f 1-3 -d "_" | uniq); do cont1=$(echo $(ls ${i}*)|cut -f 1 -d " "); cont2=$(echo $(ls ${i}*)|cut -f 2 -d " "); echo "\""$i"\" : [\""$cont1"\", \""${cont2}"\"]," ; done
