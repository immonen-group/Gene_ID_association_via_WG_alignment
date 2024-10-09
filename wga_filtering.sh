#!/bin/bash -l
#SBATCH -A naiss2024-5-135
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 2-00:00:00
#SBATCH -J delta_filter
#SBATCH -o delta_filter.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user milena.trabert@ebc.uu.se

module load bioinfo-tools MUMmer/4.0.0rc1
#nucmer -t 20 ref:path qry:path+

ASSEMBLIES_DIR=/proj/naiss2023-6-65/Milena/martyna_SNP_stuff/gene_association

# query
CMAC_ENA=${ASSEMBLIES_DIR}/Cmac_from_ENA_GCA_951848785.1.fasta.masked

# reference(s)
CMAC_HiC=${ASSEMBLIES_DIR}/Cmac_HiC_zooeco.fasta

# nucmer -t 20 $CMAC_ENA $CMAC_HiC


echo "start filter"
# the filtering step takes really long, much longer than for alignments of different bruchid species
# maybe because the alignment of two bruchid samples identifies many more matches? therefore making the file much bigger and filtering takes a long time?
delta-filter -m out.delta > out.delta_filtered
echo "done filter"
echo "start show-coords"
show-coords -rTH -I 90 out.delta_filtered > out.delta_filtered.coords