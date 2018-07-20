
from uniprot_viruses import *
import ncbi_taxonomy
import sys

def main():
    taxtree_filename = "../taxonomy/nodes.dmp"
    get_levels = True
    taxtree, taxlevel = ncbi_taxonomy.parse_taxtree(taxtree_filename, get_levels)
    stoplevel = "Family"
    stopcriteria = [
                    944644, # marseille
                    549779, # mimi
                    10486, # irido
                    10501, # phycodna
                    43682, # asco
                    10240, # pox
                    137992, # asfar
                   ]

    species = sys.argv[1:]
    for taxid in species:
        family = ncbi_taxonomy.climb_tax_tree_to_level(taxid, taxtree, taxlevel, stoplevel)
        print str(taxid) + "\t" + str(family)
        #if family in stopcriteria:
            #print taxid
        

if __name__ == "__main__": main()
