## find gene positions of genes of interest in the Cmac-M annotation on which Martyna did her analysis 
## in the Cmac-G assembly and annotation that Göran used for the balancing selection stuff 

## I will just name them Cmac-M and Cmac-G for Martyna and Göran because that is the least confusing way I can think of.
## Cmac-M = our non-chromosome level assembly and annotation that Martyna did her analysis on, genes of interest are in this annotation
## Cmac-G = Göran's chromosome level assembly and annotation that has the data on the data on balancing selection per gene

## the alignment file is tab separated and has 9 columns (some columns listed in https://github.com/mummer4/mummer are not applicable):
#     [S1]    Start of the alignment region in the reference sequence.
#     [E1]    End of the alignment region in the reference sequence.
#     [S2]    Start of the alignment region in the query sequence.
#     [E2]    End of the alignment region in the query sequence.
#     [LEN 1] Length of the alignment region in the reference sequence, measured in nucleotides.
#     [LEN 2] Length of the alignment region in the query sequence, measured in nucleotides.
#     [% IDY] Percent identity of the alignment, calculated as the (number of exact matches) / ([LEN 1] + insertions in the query).
#     [REF]   The reference FastA ID (Cmac-M)
#     [QUE]   The query FastA ID. (Cmac-G)


## contig name formatting in Cmac-M (bane of my existence)

# the scaffolds in the assembly are named like this: (except utg000312c_1) which has a "c" instead of an "l" for magical computer reasons
# >ENA|CATOUR010000001|CATOUR010000001.1 Callosobruchus maculatus genome assembly, contig: utg000001l_1
# >ENA|CATOUR010000002|CATOUR010000002.1 Callosobruchus maculatus genome assembly, contig: utg000002l_1
# >ENA|CATOUR010000003|CATOUR010000003.1 Callosobruchus maculatus genome assembly, contig: utg000003l_1
# >ENA|CATOUR010000004|CATOUR010000004.1 Callosobruchus maculatus genome assembly, contig: utg000005l_1 #!!!# note that the utg-names skip nmbers while the ENA-names do not! 
# the alignment coordinate file truncates the name after the first space to this 
# >ENA|CATOUR010000001|CATOUR010000001.1
# >ENA|CATOUR010000002|CATOUR010000002.1
# >ENA|CATOUR010000003|CATOUR010000003.1
# >ENA|CATOUR010000004|CATOUR010000004.1
# the annotation file has the utg names at the end
# utg000001l_1
# utg000002l_1
# utg000003l_1
# utg000005l_1

# because the numbers between ETA and utg contig names do not correspond I need to create a lookup-table (as a dictionary)


import sys

# import parse gff from a different folder
sys.path.insert(1, '/Users/miltr339/Box Sync/code/annotation_pipeline/evaluate_annotation') # caution: path[0] is reserved for script path (or '' in REPL)
import parse_gff as gff

import typing as tp
from dataclasses import dataclass, field
# dataclass generates stuff like __init__() for classes automatically 
from enum import Enum
import time


def find_annotation_position(annotation_filepath, gene_ID_list = [], annotation_names_dict = {}, verbose = False):
    # if no gene id list is given extract the information for all genes in the annotation
    # there is no failsafe if the annotation_names_dict keys does not match the contig names in the annotation!
    gff_dict = gff.parse_gff3_general(annotation_filepath, verbose = verbose)
    if len(gene_ID_list)>0 and len(annotation_names_dict) == 0:
        gff_dict_focal_genes = { key : [value.contig, int(value.start), int(value.stop)] for key, value in gff_dict.items() if key in gene_ID_list}
    elif len(gene_ID_list) == 0 and len(annotation_names_dict)>0:
        gff_dict_focal_genes = { key : [annotation_names_dict[value.contig], int(value.start), int(value.stop)] for key, value in gff_dict.items() if value.category == gff.FeatureCategory.Gene}
    elif len(gene_ID_list)>0 and len(annotation_names_dict)>0:
        gff_dict_focal_genes = { key : [annotation_names_dict[value.contig], int(value.start), int(value.stop)] for key, value in gff_dict.items() if key in gene_ID_list}
    
    else:
        gff_dict_focal_genes = { key : [value.contig, int(value.start), int(value.stop)] for key, value in gff_dict.items() if value.category == gff.FeatureCategory.Gene}

    return(gff_dict_focal_genes)



def create_M_names_dict(filepath):
    # read the contig names in the Cmac-M assembly and parse them into a dictionary that looks like { annotation_contig : alignment_coords_contig }
    names_dict = {}
    with open(filepath, "r") as file:
        for line in file:
            line = line.strip().split(" ")
            names_dict[line[-1]] = line[0].split(">")[-1]
    return(names_dict)


# read the alignment coordinates file in a general way, so that it does not matter which assembly was the query and which the reference. 
# make a dictionary where the top-level keys are all the unique contig names that appear for both reference and query, and the values
# are themselves dictionaries that further parse all the lines that contain the key contig names. 
# 
# Example (note that cmac-M_contig1 appears in both entries, this method purposefully duplicates data like this to make the searching more general)
# dict =  {
#   cmac-M_contig1 : [ 
#       [ cmac-M_start1, cmac-M_end1, Cmac-G_contig1, cmac-G_start1 , cmac-G_end1 ] , 
#       [ cmac-M_start2, cmac-M_end2, Cmac-G_contig2, cmac-G_start2 , cmac-G_end2 ] , 
#   ]
#   cmac-G_contig1 : [ 
#       [ cmac-G_start1, cmac-G_end1, Cmac-M_contig1, cmac-M_start1 , cmac-M_end1 ] , 
#       [ cmac-G_start2, cmac-G_end2, Cmac-M_contig2, cmac-M_start2 , cmac-M_end2 ] , 
#   ]
# }
# 
# Use this to return the corresponding match for the search query (first find contig name, and then return the target list who's interval contains the gene coordinates)
# 
# Or I could try to write a class?

@dataclass
class Region:
    contig_name:str
    start:int 
    end:int 
    partner_name:str
    partner_start:int
    partner_end:int
    seq_ident:float
    # target_region:tp.List[self.partner_name, partner_start, partner_end] = field(default_factory=list) # initialize with empty dictionary

    @property
    def show(self):
        print(f"\tfrom {self.start} to {self.end} aligns to {self.partner_name} from {self.partner_start} to {self.partner_end} with {self.seq_ident}% identity")

@dataclass
class ContigAlignments:
    name:str
    alignment_regions:tp.List[Region] = field(default_factory=list)

    def add_region(self, newregion:Region):
        self.alignment_regions.append(newregion)

    @property
    def show(self):
        print(f"there are {len(self.alignment_regions)} alignments on contig: {self.name}")
        if len(self.alignment_regions) < 20:
            for partner in self.alignment_regions:
                partner.show


# parse alignment file into the above defined dataclasses
def parse_alignment_coordinates(coords_filepath, verbose = True):
    alignments = {} 
    if verbose:
        print(f"coordinates file that is being parsed: {coords_filepath}")
        start_time = time.perf_counter()
    
    with open(coords_filepath, "r") as coords_file:
        for line in coords_file:
            start1, end1, start2, end2, len1, len2, seq_id, name1, name2 =line.strip().split()
            # read lines into Region classes from both directions
            region1 = Region(contig_name=name1, start=int(start1), end=int(end1), partner_name=name2, partner_start=int(start2), partner_end=int(end2), seq_ident=float(seq_id))
            region2 = Region(contig_name=name2, start=int(start2), end=int(end2), partner_name=name1, partner_start=int(start1), partner_end=int(end1), seq_ident=float(seq_id))
            
            # continue to populate alignments dictionary
            read_contigs = alignments.keys()

            # if the contig has not been part of a previous alignment, add it to the dictionary
            if region1.contig_name not in read_contigs:
                alignments[region1.contig_name] = ContigAlignments(region1.contig_name, [region1])
            # if the contig is already present in the list from other alignment hits
            elif region1.contig_name in read_contigs:
                alignments[region1.contig_name].add_region(region1)

            # (same with region2)
            if region2.contig_name not in read_contigs:
                alignments[region2.contig_name] = ContigAlignments(region2.contig_name, [region2])
            elif region2.contig_name in read_contigs:
                alignments[region2.contig_name].add_region(region2)
    
    if verbose:
        end_time = time.perf_counter()
        execution_time = end_time - start_time
        print(f"\tparsing time: {execution_time:.2f} seconds")

    return(alignments)       


# Find corresponding alignment regions in partner alignment
def get_gene_aligned_regions(gene_positions_dict, parsed_alignment_file, verbose = True):
    # aligned_gene_regions = {gene_ID : [target_contig, target_start, target_end]}
    aligned_gene_regions = {}
    genes_in_aligned_regions = 0

    # if verbose:
    #     print(f"finding aligned regions of {len(gene_positions_dict)} genes.")

    for gene_ID, gene_position in gene_positions_dict.items():

        source_contig, gene_start, gene_end = gene_position
        num_overlaps = 0 # if a region with a gene is aligned more than once
        
        if verbose:
            print(f"{gene_ID} at: \t{gene_position}")
        
        try:
            current_contig_aligned_regions = parsed_alignment_file[source_contig]
        except:
            # if the parsed alignment file doesn't have any aligned regions in source_contig, then it won't appear in the parsed_alignment_file
            continue
            # continue to next iteration and skip below for-loop for the current gene_ID

        for alignment in current_contig_aligned_regions.alignment_regions:
            len_overlap = 0
            overlaps_list = []
            target_list = ["", 0, 0] # [target_contig, target_start, target_end]
            # test if the gene_ID overlaps with any of the aligned regions
            # (no need to be a complete overlap, only start or end needs to be inside)
            start_in_aln = gene_start > alignment.start and gene_end < alignment.end
            end_in_aln = gene_end > alignment.start and gene_end < alignment.end
            gene_len = gene_end - gene_start

            if start_in_aln or end_in_aln: 
                genes_in_aligned_regions += 1
                num_overlaps += 1
                overlaps_list.append(alignment)
                target_list[0] = alignment.partner_name

                if start_in_aln and end_in_aln:
                    # calculate target coordinates
                    dist_start = gene_start - alignment.start
                    target_list[1] = alignment.partner_start + dist_start
                    target_list[2] = target_list[1] + gene_len

                # calculate length of overlap
                if start_in_aln and not end_in_aln: 
                    dist_start = abs(gene_start-alignment.start)
                    dist_end = abs(gene_start-alignment.end)
                    if dist_start > dist_end:
                        len_overlap = dist_start
                    else:
                        len_overlap = dist_end
                    # calculate target coordinates
                    target_list[1] = alignment.partner_end - dist_end
                    target_list[2] = alignment.partner_end + gene_len - len_overlap

                if end_in_aln and not start_in_aln:
                    dist_start = abs(gene_end-alignment.start)
                    dist_end = abs(gene_end-alignment.end)
                    if dist_start > dist_end:
                        len_overlap = dist_start
                    else:
                        len_overlap = dist_end
                    # calculate target coordinates
                    target_list[2] = alignment.partner_start - dist_end
                    target_list[1] = alignment.partner_start - gene_len + len_overlap
                
                if verbose:
                    print(f"finding aligned regions of {gene_ID}")
                    print(f"\t target coordinates: {target_list}")
                    print(f"\t length of the overlap: {len_overlap}, length of the gene: {gene_len}")


                # make sure target_list is sorted correctly with target_list[1] < target_list[2]
                if target_list[2]<target_list[1]:
                    target_list = [target_list[0], target_list[2], target_list[1]]
                # filter genes that have negative coordinates calculated
                if target_list[2] < 0 or target_list[1] < 0:
                    continue

                # 
                aligned_gene_regions[gene_ID] = target_list ## TODO fix assignments? somehow it only finds one
                # break 


            if verbose:
                if num_overlaps>1:
                    print("more than one overlap")
        
    if verbose:
        number_of_genes = len(gene_positions_dict)
        print(f"of {number_of_genes} genes, {genes_in_aligned_regions} are in aligned regions")
        
    return(aligned_gene_regions)
        #break


def write_aln_regions_to_file(aln_coords_dict, outfile_name):
    with open(outfile_name, "w") as outfile:
        for gene_ID, aln_region in aln_coords_dict.items():
            aln_region = "\t".join([str(item) for item in aln_region])
            line = f"{aln_region}\t{gene_ID}\n"
            outfile.write(line)
    print(f"outfile written to {outfile_name}")

    
    

if __name__ == "__main__":

    # alignment_coords_filepath = "/proj/naiss2023-6-65/Milena/martyna_SNP_stuff/gene_association/out.delta_filtered.coords" # on uppmax

    alignment_coords_filepath_ENA = "/Users/miltr339/work/martyna_SNP_stuff/gene_associaton_with_HiC_assembly/ENA_ref_out.delta_filtered.coords"
    alignment_coords_filepath_HiC = "/Users/miltr339/work/martyna_SNP_stuff/gene_associaton_with_HiC_assembly/HiC_ref_out.delta_filtered.coords"

    Cmac_M_filepath = "/Users/miltr339/work/martyna_SNP_stuff/gene_associaton_with_HiC_assembly/Cmac-M.gff"
    Cmac_G_filepath = "/Users/miltr339/work/martyna_SNP_stuff/gene_associaton_with_HiC_assembly/Cmac-G_HiC.gtf"

    # gff hell mitigation
    lookup_table_filepath = "/Users/miltr339/work/martyna_SNP_stuff/gene_associaton_with_HiC_assembly/Cmac-M_complete_contig_names.txt"
    cmac_M_contig_names = create_M_names_dict(lookup_table_filepath)

    martyna_gene_IDs = [
            "cmac_gs1_g13105",
            "cmac_gs1_g13827",
            "cmac_gs1_g15237",
            "cmac_gs1_g15239",
            "cmac_gs1_g15555",
            "cmac_gs1_g15562",
            "cmac_gs1_g15595",
            "cmac_gs1_g15673",
            "cmac_gs1_g17192",
            "cmac_gs1_g1753",
            "cmac_gs1_g19520",
            "cmac_gs1_g19603",
            "cmac_gs1_g20671",
            "cmac_gs1_g20680",
            "cmac_gs1_g20737",
            "cmac_gs1_g20822",
            "cmac_gs1_g20828",
            "cmac_gs1_g20840",
            "cmac_gs1_g21194",
            "cmac_gs1_g2125",
            "cmac_gs1_g2162",
            "cmac_gs1_g2290",
            "cmac_gs1_g23644",
            "cmac_gs1_g24004",
            "cmac_gs1_g2416",
            "cmac_gs1_g25062",
            "cmac_gs1_g25619",
            "cmac_gs1_g27769",
            "cmac_gs1_g27837",
            "cmac_gs1_g27890",
            "cmac_gs1_g28503",
            "cmac_gs1_g29056",
            "cmac_gs1_g30457",
            "cmac_gs1_g3589",
            "cmac_gs1_g3622",
            "cmac_gs1_g3654",
            "cmac_gs1_g3736",
            "cmac_gs1_g3747",
            "cmac_gs1_g3773",
            "cmac_gs1_g459",
            "cmac_gs1_g470",
            "cmac_gs1_g613",
            "cmac_gs1_g6346",
            "cmac_gs1_g6472",
            "cmac_gs1_g7360",
            "cmac_gs1_g7472",
            "cmac_gs1_g7474",
            "cmac_gs1_g7475",
            "cmac_gs1_g7476",
            "cmac_gs1_g7480",
            "cmac_gs1_g7526",
            "cmac_gs1_g7654",
            "cmac_gs1_g784",
            "cmac_gs1_g8508",
            "cmac_gs1_g8509",
            "cmac_gs1_g9470",
            "cmac_gs1_g9667",
            "cmac_gs1_g9726",
            "cmac_gs1_g9856",
    ]


    #####################################################
    ###### make gene association table for Martyna ######
    #####################################################
    if True:
        cmac_M_parsed = gff.parse_gff3_general(Cmac_M_filepath)
        cmac_M_gene_IDs = [key for key, value in cmac_M_parsed.items() if value.category == gff.FeatureCategory.Gene]
        print(f"{len(cmac_M_gene_IDs)} genes in the Cmac-M annotation")
        # print(cmac_M_contig_names["utg000312c_1"])

        # make dictionary of all genes of interest: { gene_ID : [contig_name,start,end]}
        cmac_M_gene_pos = find_annotation_position(Cmac_M_filepath, gene_ID_list= cmac_M_gene_IDs, annotation_names_dict=cmac_M_contig_names)

        # parse alignment file into a dict sorted by contigs 
        alignment_pairs_dict_ENA = parse_alignment_coordinates(coords_filepath=alignment_coords_filepath_ENA)
        alignment_pairs_dict_HiC = parse_alignment_coordinates(coords_filepath=alignment_coords_filepath_HiC)

        aligned_genes_ENA = get_gene_aligned_regions(cmac_M_gene_pos, alignment_pairs_dict_ENA, verbose=False)
        print(f" --> done ENA alignemnt")
        aligned_genes_HiC = get_gene_aligned_regions(cmac_M_gene_pos, alignment_pairs_dict_HiC, verbose=False)
        print(f" --> done HiC alignemnt")

        print()
        write_aln_regions_to_file(aligned_genes_ENA, "/Users/miltr339/work/martyna_SNP_stuff/gene_associaton_with_HiC_assembly/Cmac-G_gene_coordinates_ENA_ref_complete.tsv")
        write_aln_regions_to_file(aligned_genes_HiC, "/Users/miltr339/work/martyna_SNP_stuff/gene_associaton_with_HiC_assembly/Cmac-G_gene_coordinates_HiC_ref_complete.tsv")

        ### use bedtools to get the genes from the gff that intersect with these gene coordinates
        #
        # grep "\tgene"  Cmac-G_HiC.gtf | awk -F'\t' '$4 ~ /^[0-9]+$/' |  bedtools intersect -wa -wb -a stdin -b Cmac-G_gene_coordinates_ENA_ref_complete.tsv Cmac-G_gene_coordinates_HiC_ref_complete.tsv -filenames > Cmac-G_gene_matches_both_alns_complete.tsv
        # awk -F "\t" '{print $1, $4, $5, $9, $10, $11, $12, $13, $14}' Cmac-G_gene_matches_both_alns_complete.tsv | awk -F ";" '{print $1, $NF}' > Cmac-G_gene_matches_both_alns_simplified_complete.tsv 
        # 



    #####################################################
    ###### make gene association table for Martyna ######
    #####################################################
    if False:

        # gene IDs of the genes associated with significant SNPs
        # extracted with 
        #   cut -f2 background_OGs.tsv | sort | uniq

        # make dictionary of all genes of interest: { gene_ID : [contig_name,start,end]}
        cmac_M_gene_pos = find_annotation_position(Cmac_M_filepath, gene_ID_list= martyna_gene_IDs, annotation_names_dict=cmac_M_contig_names)

        # parse alignment file into a dict sorted by contigs 
        alignment_pairs_dict_ENA = parse_alignment_coordinates(coords_filepath=alignment_coords_filepath_ENA)
        alignment_pairs_dict_HiC = parse_alignment_coordinates(coords_filepath=alignment_coords_filepath_HiC)

        aligned_genes_ENA = get_gene_aligned_regions(cmac_M_gene_pos, alignment_pairs_dict_ENA, verbose=False)
        aligned_genes_HiC = get_gene_aligned_regions(cmac_M_gene_pos, alignment_pairs_dict_HiC, verbose=False)

        write_aln_regions_to_file(aligned_genes_ENA, "/Users/miltr339/work/martyna_SNP_stuff/gene_associaton_with_HiC_assembly/Cmac-G_gene_coordinates_ENA_ref.tsv")
        write_aln_regions_to_file(aligned_genes_HiC, "/Users/miltr339/work/martyna_SNP_stuff/gene_associaton_with_HiC_assembly/Cmac-G_gene_coordinates_HiC_ref.tsv")
        

        ### use bedtools to get the genes from the gff that intersect with these gene coordinates
        #
        # grep "\tgene"  Cmac-G_HiC.gtf | awk -F'\t' '$4 ~ /^[0-9]+$/' |  bedtools intersect -wa -wb -a stdin -b Cmac-G_gene_coordinates_* -filenames > Cmac-G_gene_matches_both_alns.tsv
        # awk -F "\t" '{print $1, $4, $5, $9, $10, $11, $12, $13, $14}' Cmac-G_gene_matches_both_alns.tsv | awk -F ";" '{print $1, $NF}' > Cmac-G_gene_matches_both_alns_simplified.tsv 
        # 
