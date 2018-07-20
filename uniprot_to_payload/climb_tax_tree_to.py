
from uniprot_viruses import *
import ncbi_taxonomy
import sys

def main():
    taxtree_filename = "../taxonomy/nodes.dmp"
    get_levels = True
    taxtree, taxlevel = ncbi_taxonomy.parse_taxtree(taxtree_filename, get_levels)

    species_file = sys.argv[1]
    with open(species_file, 'r') as f:
        for l in f:
            taxid = l.rstrip('\n')
            lineage = ncbi_taxonomy.climb_tax_tree_to_level(taxid, taxtree, taxlevel, "species")
            print str(taxid) + "\t" + str(lineage)
        

if __name__ == "__main__": main()
