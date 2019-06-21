#!/bin/bash -l
#SBATCH -D /pylon5/mc5phkp/ljcohen/kfish_pilon/Foli/
#SBATCH -J Foli_busco
#SBATCH -o /pylon5/mc5phkp/ljcohen/kfish_pilon/Foli/busco-Foli%j.o
#SBATCH -e /pylon5/mc5phkp/ljcohen/kfish_pilon/Foli/busco-Foli%j.o
#SBATCH -t 6:00:00
#SBATCH -p RM
#SBATCH --ntasks-per-node 8
#SBATCH --cpus-per-task 1
source activate dammit_conda
export AUGUSTUS_CONFIG_PATH=/pylon5/mc5phkp/ljcohen/config
run_BUSCO.py --in Foli_pilon.fasta -m genome --out Foli --lineage_path /pylon5/bi5fpmp/ljcohen/dammit/busco2db/eukaryota_odb9
