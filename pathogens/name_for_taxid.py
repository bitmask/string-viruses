
import ncbi_taxonomy
import sys

def main():
    ncbi_tax_file = "../taxonomy/names.dmp"
    get_levels = False
    taxtree, taxnames = ncbi_taxonomy.parse_ncbi_taxonomy_names(ncbi_tax_file, True)

    host_pathogen_list = sys.argv[1]
    with open(host_pathogen_list, 'r') as f:
        for l in f:
            taxid = l.rstrip('\n')
            taxid = int(taxid)
            if taxid in taxnames:
                print str(taxid) + "\t" + taxnames[taxid]
            else:
                print taxid
        

if __name__ == "__main__": main()
