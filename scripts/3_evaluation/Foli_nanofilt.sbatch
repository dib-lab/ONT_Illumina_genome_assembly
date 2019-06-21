#!/bin/bash -l
#SBATCH -D /pylon5/mc5phkp/ljcohen/kfish_ONT_qc
#SBATCH -J Foli-nanoplot
#SBATCH -o /pylon5/mc5phkp/ljcohen/kfish_ONT_qc/nanoplot-Foli%j.o
#SBATCH -e /pylon5/mc5phkp/ljcohen/kfish_ONT_qc/nanoplot-Foli%j.o
#SBATCH -t 8:00:00
#SBATCH -p RM-small
#SBATCH --ntasks-per-node 4
#SBATCH --cpus-per-task 1

source ~/.bashrc
source activate ONT

cd /pylon5/mc5phkp/ljcohen/kfish_ONT_qc

NanoPlot --verbose --N50 --loglength -o Foli_nanoplot_fixed --maxlength 25 -p Foli --raw --store --fastq F_olivaceus_trimmed.Q5.fq
