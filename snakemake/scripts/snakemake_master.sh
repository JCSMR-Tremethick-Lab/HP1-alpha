#!/bin/bash
#PBS -P kv78
#PBS -l walltime=24:00:00
#PBS -l wd
#PBS -e /home/150/sxk150/qsub_error
#PBS -o /home/150/sxk150/qsub_out
#PBS -l ncpus=1
#PBS -l mem=4G
#PBS -M skurscheid@gmail.com
#PBS -m abe
#PBS -N snakemake_master
#PBS -q express

source ~/.bashrc

/short/rl2/miniconda3/envs/snakemake/bin/snakemake -s /home/sxk150/HP1-alpha all -n\
        --configfile /home/sxk150/HP1-alpha/configs/config_ChIP-Seq.json\
        --use-conda\
        --jobs 32\
        -d /short/kv78/HP1-alpha\
        --local-cores 1\
        -pr\
        --cluster "qsub -P {cluster.P}\
                    -l ncpus={cluster.ncpus} \
                    -q {cluster.queue} \
                    -l mem={cluster.mem} \
                    -l wd\
                    -l walltime={cluster.walltime}\
                    -e {cluster.error_out_dir} \
                    -o {cluster.std1_out_dir}" \
        --cluster-config /home/sxk150/HP1-alpha/configs/cluster.json\
        --rerun-incomplete


