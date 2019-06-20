#!/bin/bash -l
#SBATCH -D /pylon5/mc5phkp/ljcohen/kfish_canu
#SBATCH -J Foli-canu
#SBATCH -o /pylon5/mc5phkp/ljcohen/kfish_canu/canu-Foli%j.o
#SBATCH -e /pylon5/mc5phkp/ljcohen/kfish_canu/canu-Foli%j.o
#SBATCH -t 120:00:00
#SBATCH -p LM
#SBATCH --mem=500GB
#SBATCH --ntasks-per-node 14
#SBATCH --cpus-per-task 1

source ~/.bashrc
source activate ONT

cd /pylon5/mc5phkp/ljcohen/kfish_canu

canu useGrid=false -p Foli -d Foli-canu genomeSize=1.1m gnuplotTested=true -nanopore-raw F_olivaceus.fastq
