
from uniprot_viruses import *
import ncbi_taxonomy
import sys

def main():
    taxtree_filename = "../taxonomy/nodes.dmp"
    get_levels = False
    taxtree = ncbi_taxonomy.parse_taxtree(taxtree_filename, get_levels)

    species = sys.argv[1:]
    for taxid in species:
        lineage = ncbi_taxonomy.climb_tax_tree(taxid, taxtree)
        print str(taxid) + "\t" + str(lineage)
        

if __name__ == "__main__": main()
