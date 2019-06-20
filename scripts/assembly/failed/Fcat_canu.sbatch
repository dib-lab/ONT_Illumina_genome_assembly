#!/bin/bash -l
#SBATCH -D /pylon5/mc5phkp/ljcohen/kfish_canu
#SBATCH -J Fcat-canu
#SBATCH -o /pylon5/mc5phkp/ljcohen/kfish_canu/canu-Fcat%j.o
#SBATCH -e /pylon5/mc5phkp/ljcohen/kfish_canu/canu-Fcat%j.o
#SBATCH -t 120:00:00
#SBATCH -p LM
#SBATCH --mem=1000GB

source ~/.bashrc
source activate ONT

cd /pylon5/mc5phkp/ljcohen/kfish_canu

canu useGrid=false stopOnReadQuality=false -p Fcat -d Fcat-canu genomeSize=1.1m gnuplotTested=true -nanopore-raw F_catenatus.fastq
