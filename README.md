# Associate gene positions between two annotated assemblies via WGA
Assoicate gene positions of two separate annotated assemblies of the same species via whole genome alignments. In principle this is similar to how annotation liftover works where you associate sequences via an alignment, but I used two whole genome alignments (each assembly is used as the query once) instead of just aligning protein sequences. This is because the C. maculatus genomes used here have many paralogs and I wanted to be certain that the paralogs are associated with each other correctly by also taking information from the surrounding sequences into account.

# input data

* **HiC assembly**: also referred to as **Cmac-G** from [Arnqvist et al. 2024](https://academic.oup.com/g3journal/article/14/2/jkad266/7471854)
* **ENA assembly**: also referred to as **Cmac-M** from [Kaufmann et al. 2023](https://academic.oup.com/mbe/article/40/8/msad167/7227908) (Y_s haplotype)

## basic workflow

### alignment
make an alignment with nucmer using the script `cmac_wga.sh`, and then filter the alignent with `wga_filtering.sh`. This results then in an alignent file that shows the positions of the query regions and where they aligned to the reference. There are two alignment files, one for each input assembly where it was used as the query.

### gene position transfer with `gene_association.py`
For more technical details see the comments in the python script directly.
1. Find annotation positions (contig name and coordinates) of a list of all feature IDs of interest using the annotation files. This is done with the class implemented in `parse_gff.py`, where the IDs of all features of the class `Gene` are extracted.
2. Read the alignment file into the `ContigAlignments` class, that contains the contig name of the query sequence as `name` and then a list of `Region` class instances that represent all aligned regions in `name`.
3. Calculate the inferred coordinates of the genes from 1. using the alignment coordinates in 2. and get a list of inferred gene coordinates

### get ids of transferred gene positions
Use bedtools to get overlaps between above inferred gene positions and gene positions in the original annotation. bedtools accepts *any overlap at all*, at least 1bp. If filtering by percent or length of overlap is desired then the output file can be filtered afterwards. The filenames here come from `gene_association.py`

Use bedtools to get the overlapping positions. `mRNA` overlaps transcripts, but change for any feature of your choice. awk is used to get rid of useless information from the gff file.
```
grep "\tmRNA"  Cmac-G_HiC.gtf | awk -F'\t' '$4 ~ /^[0-9]+$/' |  bedtools intersect -wa -wb -a stdin -b Cmac-G_gene_coordinates_ENA_ref_complete.tsv Cmac-G_gene_coordinates_HiC_ref_complete.tsv -filenames > Cmac-G_gene_matches_both_alns_complete.tsv
```
parse the output file
```
awk -F "\t" '{print $1, $4, $5, $9, $10, $11, $12, $13, $14}' Cmac-G_gene_matches_both_alns_complete.tsv | awk -F ";" '{print $1, $NF}' > Cmac-G_gene_matches_both_alns_simplified_complete.tsv 
```

### outfile
The output file at the end has these columns:

1-4 : original HiC annotation (scaffold, coordinates, gene ID)

5 : I used information from two alignments, one with ENA as reference and one with HiC as reference, and this tells which one was used to infer this position

6-8 : HiC gene positions of ENA annotated genes (HiC scaffold and coordinates inferred via alignment as explained above).

9: original ENA gene ID

## Things to note about the output file

1. There are some duplicate ENA gene IDs. This happens when the gene found a HiC partner for both ENA and HiC reference alignments, so it appears twice.
2. If you check for overlap with transcripts and not genes, the for example cmac_gs1_g9969 matches CALMACT00000001113 and CALMACT00000001114, which are transcripts of the same gene CALMACG00000000693.
3. Sometimes an ENA-inferred position overlaps with two or more separate HiC genes (for example cmac_gs1_g2162), I just included all of these. It’s possible to filter them for largest overlap probably, but I don’t think that’s
4. It says [picked_reason "no_overlap"] behind a HiC gene ID sometimes, and that is part of the functional annotation from the original HiC annotation file, it doesn’t have anything to do with gene overlap from my script. I don’t know what it means.
5. If a gene ID from either annotation is missing from my file, then it does not overlap with any gene from the other annotation in the alignments.

# infographic
<img src="https://github.com/milena-t/Gene_position_association_via_WG_alignment/blob/main/ENA_HiC_gene_coordinate_association.png" alt="drawing" width="450"/>
