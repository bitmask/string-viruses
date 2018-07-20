
from uniprot_viruses import *
import ncbi_taxonomy

def main():
    filename_uniprot = "data_in/uniprotkb-viruses-reviewed-date-keywords-chain.txt"
    taxfilename = "../taxonomy/viruses.categories.txt"
    taxfilename_edit = "../taxonomy/dengue-edit"
    ncbi_tax_file = "../taxonomy/names.dmp"
    taxtree_filename = "../taxonomy/nodes.dmp"
    simap_filename = "../SIMAP/virus_vs_all/homologs"


    out_clades_file = "data_out/initial_clade_assignments"
    
    # TODO read in the object output by uniprot_to_payload
    # but for now, we do it this way

    # read the tax and uniprot data
    taxmap = parse_taxids(taxfilename, taxfilename_edit)  # mapping isolates to species
    uniprot = parse_uniprot(filename_uniprot, taxmap)
    
    # choose the isolates that will be in STRING
    chosen = choose_isolates(uniprot)
    # then get all the proteins that correspond to these isolates
    chosen_proteins = get_chosen_proteins(uniprot, chosen)
    # cleave any polyproteins based on the chain entries in uniprot
    chosen_proteins = cleave_polyproteins(chosen_proteins)

    # convert host name to taxid
    ncbi_tax, taxid_to_name = ncbi_taxonomy.parse_ncbi_taxonomy_names(ncbi_tax_file, True) # mapping name to taxids
    chosen_proteins = convert_host_name_to_taxids(chosen_proteins, ncbi_tax)

    get_levels = True
    taxtree, taxlevel = ncbi_taxonomy.parse_taxtree(taxtree_filename, get_levels)

    # TODO all the above should be encapsulated

    clades = generate_clades(chosen_proteins, taxtree, taxlevel, "Family")
    print_clades(out_clades_file, clades)

if __name__ == "__main__": main()
