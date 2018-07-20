
import ncbi_taxonomy
import sys

def main():
    taxtree_filename = "../taxonomy/nodes.dmp"
    get_levels = False
    taxtree = ncbi_taxonomy.parse_taxtree(taxtree_filename, get_levels)

    host_pathogen_list = sys.argv[1]
    with open(host_pathogen_list, 'r') as f:
        for l in f:
            host, pathogen = l.rstrip('\n').split('\t')
        
            lineage = ncbi_taxonomy.climb_tax_tree(host, taxtree)
            host_domain = ''
            if 2 in lineage:
                host_domain = "bacteria"
            elif 2157 in lineage:
                host_domain = "archaea"
            elif 9606 in lineage:
                host_domain = "human"
            elif 7742 in lineage:
                host_domain = "vertebrates"
            elif 33208 in lineage:
                host_domain = "animals"
            elif 3193 in lineage:
                host_domain = "plants"
            elif 6960 in lineage:
                host_domain = "insects"
            elif 4751 in lineage:
                host_domain = "fungi"
            elif 2759 in lineage:
                host_domain = "eukaryota"
            else:
                host_domain = "none"

            lineage = ncbi_taxonomy.climb_tax_tree(pathogen, taxtree)
            pathogen_domain = ''
            if 2 in lineage:
                pathogen_domain = "bacteria"
            elif 2157 in lineage:
                pathogen_domain = "archaea"
            elif 5794 in lineage:
                pathogen_domain = "apicomplexa (protists)"
            elif 5653 in lineage:
                pathogen_domain = "kinetoplastida (protists)"
            elif 554915 in lineage:
                pathogen_domain = "amoebozoa"
            elif 6340 in lineage:
                pathogen_domain = "annelida (worms)"
            elif 6157 in lineage:
                pathogen_domain = "platyhelminthes (worms)"
            elif 6231 in lineage:
                pathogen_domain = "nematoda (worms)"
            elif 10232 in lineage:
                pathogen_domain = "acanthocephala (worms)"
            elif 33340 in lineage:
                pathogen_domain = "neoptera (fleas, louse)"
            elif 4751 in lineage:
                pathogen_domain = "fungi"
            elif 7742 in lineage:
                pathogen_domain = "vertebrates"
            elif 2759 in lineage:
                pathogen_domain = "eukaryota"
            elif 10239 in lineage:
                pathogen_domain = "viruses"
            else:
                pathogen_domain = "none"

            print host + "\t" + host_domain + "\t" + pathogen + "\t" + pathogen_domain
        

if __name__ == "__main__": main()
