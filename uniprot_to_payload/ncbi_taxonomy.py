def parse_taxtree(filename, returnlevel=False):
    with open(filename) as f:
        # read in the ncbi tax nodes.dmp file
        taxtree = {}
        taxlevel = {}
        for line in f:
            array = line.rstrip("\n").split("|")
            taxid = int(array[0].strip("\t"))
            parent = int(array[1].strip("\t"))
            level = array[2].strip("\t")
            taxtree[taxid] = parent
            taxlevel[taxid] = level
        if returnlevel:
            return (taxtree, taxlevel)
        else:
            return taxtree

def parse_taxtree_down(filename):
    with open(filename) as f:
        # read in the ncbi tax nodes.dmp file and store descendants 
        taxtree = {}
        taxlevel = {}
        for line in f:
            array = line.rstrip("\n").split("|")
            taxid = int(array[0].strip("\t"))
            parent = int(array[1].strip("\t"))
            level = array[2].strip("\t")
            pair = (taxid, level)
            if parent in taxtree:
                taxtree[parent].append(pair)
            else:
                taxtree[parent] = [pair]
        return taxtree

def climb_tax_tree_to_level(taxid, taxtree, taxlevel, stoplevel):
    ''' stoplevel is the name of the taxonomic level we want to stop at eg 'Family'
    '''
    taxid = int(taxid)
    if taxid in taxtree and taxid in taxlevel:
        parent = taxid
        parentlevel = taxlevel[taxid]
    else:
        # TODO this is a problem
        return -1
    #print "climbing"
    while parent != 1 and parentlevel.lower() != stoplevel.lower():
        #print "\t parent is " + str(parent)
        if parent in taxtree and parent in taxlevel:
            parent = taxtree[parent]
            parentlevel = taxlevel[parent]
            #print "\t \tis now " + str(parent) + " at " + str(parentlevel)
        else:
            #print "parent " + str(parent) + " not found in taxtree or taxlevel.  this is very bad"
            return -1
    #print "found " + str(parent)
    return parent

def climb_tax_tree(taxid, taxtree):
    # climb the tree to the root and return the lineage
    taxid = int(taxid)

    if taxid in taxtree:
        parent = taxtree[taxid]
    else:
        # TODO this is a problem
        #print "id " + str(taxid) + " not found in taxtree.  this is bad"
        return []

    lineage = [taxid]
    while parent != 1:
        lineage.insert(0, parent)
        if parent in taxtree:
            parent = taxtree[parent]
        else:
            print "parent " + str(parent) + " not found in taxtree.  this is very bad"
            return -1
    return lineage

def flatten(container):
    # thanks https://stackoverflow.com/questions/10823877/what-is-the-fastest-way-to-flatten-arbitrarily-nested-lists-in-python
    for i in container:
        if isinstance(i, list) or isinstance(i, tuple):
            for j in flatten(i):
                yield j
        else:
            yield i

def desc_recursive(taxid, taxtree, taxid_to_name, level, andbelow, excl_level, reached, depth, override):
    matching_children = []
    #print str(depth) + " desc_recursive taxid " + str(taxid) + " level " + str(level) + " reached " + str(reached) + " override " + str(override)
    if taxid in taxtree:
        children = taxtree[taxid]
        for child in children:
            child_taxid, child_level = child
            #print str(depth) + " child name " + taxid_to_name[child_taxid]
            reached = reached
            local_override = override # don't transfer this value to the siblings

            # if we're at a level with shitty data, ignore it and things below it
            if "environmental samples" in taxid_to_name[child_taxid] or "unclassified" in taxid_to_name[child_taxid] or "uncultured" in taxid_to_name[child_taxid]:
                #print "shitty"
                local_override = True

            #print str(depth) + " found child " + str(child_taxid) + " level " + str(child_level)
            if not local_override and (child_level == level or (andbelow == True and reached == True)):
                excl_this_level = False # maybe we want everything under the given level, and not the actual level
                if child_level == level:
                    reached = True # have reached the desired level
                    if excl_level:
                        excl_this_level = True
                #print str(depth) + " added child "
                if not excl_this_level:
                    matching_children.append(child_taxid)
                    #matching_children.append(taxid_to_name[child_taxid])

            new_children = desc_recursive(child_taxid, taxtree, taxid_to_name, level, andbelow, excl_level, reached, depth+1, local_override)

            if new_children:
                matching_children.append(new_children)
        #print str(depth) + " return " + str(matching_children)
        return matching_children
    else:
        # reached the end of the recursion 
        # (unless we are in a case where a parent is refered to but isn't in the tax tree -- haven't seen this happen yet)
        #print str(depth) + " return nothing " + str(matching_children)
        return matching_children

def descend_tax_tree(taxid, taxtree, taxid_to_name, level, andbelow, excl_level, reached=False):
    # return all the nodes from a given branch
    #reached  # have we reached the desired level yet?

    taxid = int(taxid)
    # everything at the specified level (and below) for the given taxid

    depth = 0
    override = False # used to not print env samples and unclassified crap

    #flatten the list
    return list(flatten(desc_recursive(taxid, taxtree, taxid_to_name, level, andbelow, excl_level, reached, depth, override)))

def parse_ncbi_taxonomy_names(ncbi_tax_file, reverse):
    taxnames = {}
    taxid_to_name = {}
    with open(ncbi_tax_file) as f:
        for line in f:
            taxid, unused1, name, unused2, unused3, unused4, nametype, unused5 = line.rstrip("\n").split("\t")
            taxnames[name] = taxid
            if nametype == "scientific name":
                taxid_to_name[int(taxid)] = name
    if reverse:
        return taxnames, taxid_to_name
    else:
        return taxnames

def host_kingdom(virus, chosen_proteins, taxtree):
    eukaryotes = 2759
    bacteria = 2
    archaea = 2157

    if virus in chosen_proteins:
        array = chosen_proteins[virus]['hosts_taxids'] # [0]
        if len(array) > 0:
            if len(array[0]) > 0:
                lineage = climb_tax_tree(array[0], taxtree) # assume that all hosts will share the same kingdom
                if archaea in lineage:
                    return "archaea"
                if bacteria in lineage:
                    return "bacteria"
                if eukaryotes in lineage:
                    return "eukaryotes"
    return ''


def host_family(virus, chosen_proteins, taxtree, taxlevel):
    if virus in chosen_proteins:
        array = chosen_proteins[virus][0]['hosts_taxids']
        if len(array) > 0:
            if len(array[0]) > 0:
                taxid = climb_tax_tree_to_level(array[0], taxtree, taxlevel, "family") # TODO need family for all hosts
                # TODO get name
                return taxid
    return ''


def which_family(virus, taxtree, taxlevel):
    taxid = climb_tax_tree_to_level(virus, taxtree, taxlevel, "family")
    if taxid > 0:
        return taxid
    else:
        #print "could not find a family for " + str(virus) 
        return -1

def which_order(virus, taxtree, taxlevel):
    taxid = climb_tax_tree_to_level(virus, taxtree, taxlevel, "order")
    if taxid > 0:
        return taxid
    else:
        #print "could not find an order for " + str(virus) 
        return -1


def which_level(virus, taxtree):
    bacteria = 2
    archaea = 2157
    human = 9606
    vertebrates = 7742
    metazoa = 33208
    eukaryotes = 2759
    viruses = 10239
    

    lineage = climb_tax_tree(virus, taxtree)
    if archaea in lineage:
        return "archaea"
    if bacteria in lineage:
        return "bacteria"
    if human in lineage:
        return "human"
    if vertebrates in lineage:
        return "vertebrates"
    if metazoa in lineage:
        return "metazoa"
    if eukaryotes in lineage:
        return "eukaryotes"
    if viruses in lineage:
        return "viruses"

    #TODO
    #print "could not identify level from lineage " + str(lineage)
    return "unknown"


def which_viral_clade(virus, taxtree):

    viruses_dsDNA = 35237
    viruses_dsRNA = 35325
    viruses_ssDNA = 29258
    viruses_ssRNA_pos = 35278
    viruses_ssRNA_neg = 35301
    viruses_rt = 35268
    viruses_rtDNA = 11632
    viruses_rtRNA = [186534, 10404]

    lineage = climb_tax_tree(virus, taxtree)
    if viruses_dsDNA in lineage:
        return "dsDNA"
    if viruses_dsRNA in lineage:
        return "dsRNA"
    if viruses_ssDNA in lineage:
        return "ssDNA"
    if viruses_ssRNA_pos in lineage:
        return "+ssRNA"
    if viruses_ssRNA_neg in lineage:
        return "-ssRNA"
    if viruses_rt in lineage:
        if viruses_rtDNA in lineage:
            return "rtDNA"
        else:
            if viruses_rtRNA[0] in lineage or viruses_rtRNA[1] in lineage:
                return "rtRNA"
            else:
                return "rt" # unclassified
    return "viruses" # unclassified

def split_string(stringid):
    allitems = stringid.split(".")
    taxid = allitems[0]
    proteinid = ".".join(allitems[1:])
    return taxid, proteinid

