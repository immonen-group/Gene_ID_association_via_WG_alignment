#!/bin/bash -l
#SBATCH -A naiss2024-5-135
#SBATCH -p core
#SBATCH -n 20
#SBATCH -t 10:00:00
#SBATCH -J nucmer_maculatus_ENA_vs_HiC
#SBATCH -o nucmer_maculatus_ENA_vs_HiC.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user milena.trabert@ebc.uu.se

module load bioinfo-tools MUMmer/4.0.0rc1
#nucmer -t 20 ref:path qry:path+ 

ASSEMBLIES_DIR=/proj/naiss2023-6-65/Milena/martyna_SNP_stuff/gene_association

# reference (later called Cmac-M)
CMAC_ENA=${ASSEMBLIES_DIR}/Cmac_from_ENA_GCA_951848785.1.fasta.masked

# query (later called Cmac-G)
CMAC_HiC=${ASSEMBLIES_DIR}/Cmac_HiC_zooeco.fasta



nucmer -t 20 $CMAC_ENA $CMAC_HiC

