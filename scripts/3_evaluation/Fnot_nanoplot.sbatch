#!/bin/bash -l
#SBATCH -D /pylon5/mc5phkp/ljcohen/kfish_ONT_qc
#SBATCH -J Fnot-nanoplot
#SBATCH -o /pylon5/mc5phkp/ljcohen/kfish_ONT_qc/nanoplot-Fnot%j.o
#SBATCH -e /pylon5/mc5phkp/ljcohen/kfish_ONT_qc/nanoplot-Fnot%j.o
#SBATCH -t 8:00:00
#SBATCH -p RM-small
#SBATCH --ntasks-per-node 4
#SBATCH --cpus-per-task 1

source ~/.bashrc
source activate ONT

cd /pylon5/mc5phkp/ljcohen/kfish_ONT_qc

NanoPlot --verbose --N50 --loglength -o Fnot_nanoplot_fixed --maxlength 25 -p Fnot --raw --store --fastq F_nottii_trimmed.Q5.fq
