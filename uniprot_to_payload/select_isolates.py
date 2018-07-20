import sys
from uniprot_viruses import *
import ncbi_taxonomy
import pprint
import cPickle as pickle
import os.path

# this will generate a report of the isolates found within a proteomes directory and their corresponding species

def main():
    pickle_file = "data_out/parsed-proteomes"

    proteomes_list = "data_in/proteomes-all.tab"
    mapfile_proteomes = "data_in/proteomes.entryids"
    dirname_proteomes = "data_in/proteomes"

    ncbi_tax_file = "../taxonomy/names.dmp"
    taxtree_filename = "../taxonomy/nodes.dmp"

    pogs_file = "data_in/pVOGs"

    out_isolates = "data_out/isolates_amb"
    out_isolates_all = "data_out/isolates_all"

    out_fasta_dir = "data_out/fastafiles"
    out_pfam_file = "data_out/pfam"
    out_go_file = "data_out/goterms"
    out_functions_file = "data_out/functions"
    out_pdb_file = "data_out/pdb"
    out_uniprot_file = "data_out/uniprot"
    out_host_file = "data_out/hosts"
    out_species_file = "data_out/string_v11.species.tsv"
    out_synonyms_file = "data_out/synonyms"
    out_best_synonyms_file = "data_out/synonyms_best"
    out_functions_file = "data_out/functions"
    out_pogs_file = "data_out/pogs"

    if os.path.isfile(pickle_file):
        chosen_proteins = pickle.load( open( pickle_file, "rb" ) )
    else:
        # read proteome to entry mapping
        proteomes = parse_proteome_mapping(proteomes_list, mapfile_proteomes)

        # generate taxmap based on all entries in ncbi taxonomy so that we can map isolate -> species
        get_levels = True
        taxtree, taxlevel = ncbi_taxonomy.parse_taxtree(taxtree_filename, get_levels)
        unused, taxnames = ncbi_taxonomy.parse_ncbi_taxonomy_names(ncbi_tax_file, True)

        taxmap = {}
        for taxid in taxtree:
            species = ncbi_taxonomy.climb_tax_tree_to_level(taxid, taxtree, taxlevel, "species")
            if species == -1:
                next # above species
            taxmap[taxid] = species

        # read the tax and uniprot data
        chosen_proteins = parse_uniprot_promeome_dir(proteomes, dirname_proteomes, taxmap, True) # isolate_is_species
        sys.stderr.write("parsed uniprot xml\n")

        chosen_proteins = cleave_polyproteins(chosen_proteins)

        chosen_proteins = dedup(chosen_proteins)
        print "have " + str(len(chosen_proteins)) + " left"
        chosen_proteins = remove_too_small_genomes(chosen_proteins)

        pickle.dump(chosen_proteins, open( pickle_file, "wb" ) )


    if False:
        # run this when getting names for textmining
        # takes a really long time because it blasts swiss-prot against (some of) tremble to find more synonyms
        filename_uniprot_xml_all = "data_in/uniprot_viruses.xml"
        allspecies = parse_uniprot_xml(filename_uniprot_xml_all, taxmap)
        allspecies = cleave_polyproteins(allspecies)
    else:
        allspecies = []
    
    # print simple derivatives of the uniprot data
    print_fasta(out_fasta_dir, chosen_proteins, False, True)    # False, True -> no name details, don't print polyproteins

    print_pfam(out_pfam_file, chosen_proteins)
    print_go(out_go_file, chosen_proteins)
    print_pdb(out_pdb_file, chosen_proteins)
    print_uniprot(out_uniprot_file, chosen_proteins)
    print_hosts(out_host_file, chosen_proteins)
    print_species(out_species_file, chosen_proteins)
    print_synonyms(out_synonyms_file, out_best_synonyms_file, chosen_proteins, allspecies)
    print_functions(out_functions_file, chosen_proteins)
    map_pogs(pogs_file, out_pogs_file, chosen_proteins)

    # print output that is needed for text mining
    #print_synonyms(out_synonyms_file, out_best_synonyms_file, chosen_proteins, chosen_proteins)

    sys.stderr.write("done\n")

if __name__ == "__main__": main()
