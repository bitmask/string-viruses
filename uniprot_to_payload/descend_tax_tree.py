
import ncbi_taxonomy
import sys

def main():
    taxtree_filename = "../taxonomy/nodes.dmp"
    taxtree = ncbi_taxonomy.parse_taxtree_down(taxtree_filename)
    level = "superkingdom"
    andbelow = True # include nodes below the level
    excl_level = False # exclude the nodes at the level

    print "Warning: edit this code to know what it's doing"
    print "or better, add some command line switches"

    ncbi_tax_file = "../taxonomy/names.dmp"
    reverse = True
    taxnames, taxid_to_name = ncbi_taxonomy.parse_ncbi_taxonomy_names(ncbi_tax_file, reverse)

    parents = sys.argv[1:]
    for taxid in parents:
        children = ncbi_taxonomy.descend_tax_tree(taxid, taxtree, taxid_to_name, level, andbelow, excl_level, reached=True)
        for child in children:
            print str(child)

if __name__ == "__main__": main()
