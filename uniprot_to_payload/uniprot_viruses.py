import sys
import re
import os
import shutil
import json
import gzip
import ncbi_taxonomy as ncbi
import pprint
import collections
import shlex
import pprint
from subprocess import Popen, PIPE


def read_taxfile(taxmap, filename):
    fh = open(filename, 'r')
    for line in fh:
        v, species, isolate = line.rstrip("\n").split("\t")
        taxmap[int(isolate)] = int(species)
    fh.close()
    return taxmap


def parse_taxids(filename, taxfilename_edit):
    try:
        taxmap = {}
        taxmap = read_taxfile(taxmap, filename)          # read in the standard file
        taxmap = read_taxfile(taxmap, taxfilename_edit) # then read in the exceptions for Dengue and others
        return taxmap
    except:
        raise
            
def filter_uniprot_xml(filename):
    '''
    Filter out all non-virus entries from the xml
    '''

    import xml.etree.cElementTree as ET
    try:
        tree = ET.parse(filename)
    except:
        pass

    root = tree.getroot()

    ns = {'uniprot': 'http://uniprot.org/uniprot'}

    for entry in root.findall('uniprot:entry', ns):
        organism = entry.find('uniprot:organism', ns)
        if entry.find('uniprot:organism', ns).find('uniprot:lineage', ns).find('uniprot:taxon', ns).text != "Viruses":
            root.remove(entry)

    tree.write('trimmed.xml')

def parse_proteome_mapping(proteomes_file, entries_file):
    proteomes = {}
    with open(proteomes_file, 'r') as f:
        next(f) # skip header
        for l in f:
            c = l.rstrip("\n").split("\t")
            proteomes[c[0]] = {'isolate_taxid': c[2], 'size': c[3], 'proteins': []}
    with open(entries_file, 'r') as f:
        for l in f:
            c = l.rstrip("\n").split("\t")
            prots = c[7].split(", ")
            ref_prots = []
            for p in prots:
                p = re.sub(r': .*', '', p)
                if p in proteomes.keys():
                    proteomes[p]['proteins'].append(c[0])
                    ref_prots.append(p)
            if (len(ref_prots) > 0):
                if (len(ref_prots) > 1):
                    print "multiple reference proteomes for " + l
            else:
                print "could not find proteome in " + str(prots)
    # sanity checks
    for p in proteomes.keys():
        if int(proteomes[p]['size']) != len(proteomes[p]['proteins']):
            print "proteome size mismatch for " + p + ": declared size " + str(proteomes[p]['size']) + " ne " + str(len(proteomes[p]['proteins']))
    return proteomes
            

def parse_uniprot_promeome_dir(proteomes, dirname, taxmap, isolate_is_species):
    # if isolate_is_species is true: each isolate will be a "first class" species in string, ie isolate_id.protein_id is the string identifier
    allspecies = {}
    for proteome in proteomes.keys():
        for protein in (proteomes[proteome]['proteins']):
            xml = dirname + "/" + protein + ".xml"
            this_protein = parse_uniprot_single_xml(xml, taxmap, isolate_is_species)
            for species, stuff in this_protein.iteritems():
                if species in allspecies:
                    allspecies[species] += stuff
                else:
                    allspecies[species] = stuff
    return allspecies
                

    #import glob
    #for xml in glob.glob(dirname + "/*.xml"):
    #    this_protein = parse_uniprot_single_xml(xml, taxmap)
    #    for species, stuff in this_protein.iteritems():
    #        if species in allspecies:
    #            allspecies[species] += stuff
    #        else:
    #            allspecies[species] = stuff
    #return allspecies

def parse_uniprot_single_xml(filename, taxmap, isolate_is_species):
    '''
    Read one inidividual xml files (representing one protein) that has been downloaded to a dir, eg as part of a proteome

    XML format has changed between eggnog viruses version 1 and 2
    '''

    #sys.stderr.write("parsing " + filename + "\n")
    
    import xml.etree.cElementTree as ET
    try:
        tree = ET.parse(filename)
    except:
        raise

    root = tree.getroot()
    allspecies = {}

    source = "" # trEMBL or swiss prot
    entry_name = ""
    short_name = ""
    entry_ids = []
    isolate_name = ""
    isolate_name_details = ""
    isolate_abbrev = ""
    isolate_taxid = 0
    protein_name = ""
    aliases = []
    component_aliases = {}
    component_short_names = {}
    gene_name = ""
    date = ""
    sequence = ""
    hosts = [] # unused
    hosts_taxids = []
    chains = []
    pfam = {}
    go = []
    complete_prot = False
    reference_prot = False
    proteomes = [] # what proteomes is this protein a member of?
    polyprotein = False
    species = 0
    string_id = ""
    embl_ids = []
    gene_ids = []
    function = ""
    pdb = []
    isoforms = [] # TODO

    for child in root:

        #if "entry" in entry.tag:
            #date = entry.attrib['created']
            #source = entry.attrib['dataset'] 
        #else:
            #continue
        #for child in entry:
        if "name" in child.tag:
            entry_name = child.text
        if "accession" in child.tag:
            entry_ids.append(child.text)
        if "organism" in child.tag and "Host" not in child.tag:
            for taxentry in child:
                if "name" in taxentry.tag:
                    if taxentry.attrib['type'] == "scientific":
                        isolate_name = taxentry.text
                    if taxentry.attrib['type'] == "common":
                        isolate_name_details = taxentry.text
                    if taxentry.attrib['type'] == "synonym":
                        isolate_abbrev = taxentry.text
                if "dbReference" in taxentry.tag:
                    if taxentry.attrib['type'] == "NCBI Taxonomy":
                        isolate_taxid = int(taxentry.attrib['id'])
        if "protein" in child.tag:
            for rname in child:
                for name in rname:
                    if "fullName" in name.tag:
                        if 'recommendedName' in rname.tag:
                            protein_name = name.text
                            if re.search('olyprotein', protein_name):
                                polyprotein = True
                        else:
                            aliases.append(name.text)
                            pass
                    if 'shortName' in name.tag:
                        aliases.append(name.text)
                        if not short_name:
                            short_name = name.text
                        pass

                # get the recommened and alternate names for the chains
                if 'component' in rname.tag:
                    rec = ""
                    for r in rname:
                        for n in r:
                            if 'recommendedName' in r.tag:  # this logic is right, it is complicated
                                if 'fullName' in n.tag:
                                    rec = n.text
                                if 'shortName' in n.tag:
                                    if rec in component_aliases:
                                        component_aliases[rec].append(n.text)
                                    else:
                                        component_aliases[rec] = [n.text]
                                    component_short_names[rec] = n.text
                            if 'alternativeName' in r.tag:
                                if 'Name' in n.tag:
                                    if rec in component_aliases:
                                        component_aliases[rec].append(n.text)
                                    else:
                                        component_aliases[rec] = [n.text]
        if "gene" in child.tag:
            for name in child:
                if not gene_name:
                    gene_name = name.text
                else:
                    aliases.append(name.text) # add multiple gene names as aliases

        if "sequence" in child.tag:
            sequence = child.text.replace('\n', '')
        if "organismHost" in child.tag:
            for ref in child:
                if "dbReference" in ref.tag:
                    hosts_taxids.append(ref.attrib['id'])
        if "feature" in child.tag:
            if child.attrib['type'] == "chain":
                chain = {}
                chain['id'] = child.attrib['id']
                if 'status' not in child.attrib: # some chains have status='potential' and these seem to be crap
                    chain['name'] = child.attrib['description']
                    if chain['name'] in component_aliases:
                        chain['aliases'] = list(set(component_aliases[chain['name']]))
                    else:
                        chain['aliases'] = []
                    if chain['name'] in component_short_names:
                        chain['short_name'] = component_short_names[chain['name']]
                    complete = True
                    for loc in child:
                        for pos in loc:
                            if "begin" in pos.tag:
                                if "status" in pos.attrib and pos.attrib['status'] == "unknown":
                                    complete = False
                                else:
                                    chain['start'] = pos.attrib['position']
                            if "end" in pos.tag:
                                if 'status' in pos.attrib and pos.attrib['status'] == 'unknown':
                                    complete = False
                                else:
                                    chain['end'] = pos.attrib['position']
                    if complete:
                        chains.append(chain)
        if "dbReference" in child.tag:
            if child.attrib['type'] == "Proteomes":
                proteomes.append(child.attrib['id'])
            if child.attrib['type'] == "Pfam":
                pfam[child.attrib['id']] = {} 
            if child.attrib['type'] == "GO":
                go.append(child.attrib['id'])
            if child.attrib['type'] == "EMBL":
                for p in child:
                    if p.attrib['type'] == 'protein sequence ID':
                        embl_ids.append(p.attrib['value'])
            if child.attrib['type'] == "RefSeq":
                embl_ids.append(child.attrib['id'])
                for p in child:
                    if 'type' in p.attrib:
                        if p.attrib['type'] == 'nucleotide sequence ID':
                            embl_ids.append(p.attrib['value'])
            if child.attrib['type'] == "GeneID":
                gene_ids.append(child.attrib['id'])
            if child.attrib['type'] == "PDB":
                pdb.append(child.attrib['id'])


        if "keyword" in child.tag:
            if child.attrib['id'] == 'KW-0181':
                complete_prot = True
            if child.attrib['id'] == 'KW-1185':
                reference_prot = True

        if "comment" in child.tag:
            if child.attrib['type'] == "function":
                for c in child:
                    if "text" in c.tag:
                        function = c.text

            for c in child:
                if "isoform" in c.tag:
                    for p in c:
                        isoform_id = ""
                        add = None
                        if "id" in p.tag:
                            isoform_id = p.text
                            add = True
                        if "text" in p.tag:
                            if "frameshift" in p.text:
                                add = False # these will be very different proteins than those produced via readthrough, alt start, etc
                        if add:
                            isoform_id = isoform_id.split("-")[0]
                            isoforms.append(isoform_id)
    isoforms = list(set(isoforms))  # no duplicates
    if len(isoforms) == 1:  # if there is only one element, then don't remove it
        isoforms = []

    # viral species comes from the NCBI Taxonomy file -- roll isolate up to species level
    if isolate_is_species:
        species = int(isolate_taxid)
        string_id = str(species) + "." + str(entry_name)
    else:
        if int(isolate_taxid) in taxmap:
            species = taxmap[int(isolate_taxid)]
            string_id = str(species) + "." + str(entry_name)
        else:
            sys.stderr.write("skipping " + str(isolate_taxid) + " because it doesn't have a species (maybe update viruses.categories.txt)" + "\n")
            next

    if isolate_name_details != "":
        isolate_name += " (" + isolate_name_details + ")"
    if isolate_abbrev == "":
        isolate_abbrev = isolate_name
    else:
        isolate_name += " (" + isolate_abbrev + ")"

    this_protein = {'string_id': string_id, 'entry_id': entry_ids, 'entry_name': entry_name, 'isolate_name': isolate_name, 'isolate_abbrev': isolate_abbrev, 'gene_name': gene_name, 'protein_name': protein_name, 'short_name': short_name, 'species': species, 'isolate_taxid': isolate_taxid, 'reference_proteome': reference_prot, 'complete_proteome': complete_prot, 'proteomes': proteomes, 'polyprotein': polyprotein, 'date': date, 'hosts_taxids': list(set(hosts_taxids)), 'sequence': sequence, 'chain': chains, 'pfam': pfam, 'go': go, 'aliases': list(set(aliases)), 'embl_ids': embl_ids, 'gene_ids': gene_ids, 'source': source, 'function': function, 'pdb': pdb, 'isoforms': isoforms}

    #if not reference_prot:
        #sys.stderr.write("**** not a reference proteome " + string_id + "\n")
    
    if int(species) > 0:
        if isolate_taxid in allspecies:
            # then we are appending this_protein to the array
            if isolate_is_species:
                allspecies[isolate_taxid].append(this_protein)
        else:
            allspecies[isolate_taxid] = [this_protein]

    return allspecies

def parse_uniprot_xml(filename, taxmap):
    '''
    Generates the same datastructure as parse_uniprot, but reads an xml file instead.  this contains more information than the tsv dump, so use it instead
    '''

    return parse_uniprot_single_xml(filename, taxmap)


def parse_uniprot(filename, taxmap):
    try:
        allspecies = {}   # species: [ {}, {}, ... ]
        fh = open(filename, 'r')
        count = 0
        for line in fh:
            if count < 1:
                count += 1
                continue
            else:
                #(e, entry_name, s, protein_names, gene_name, gnp, gns, gnol, gnorf, organism_name, 
                #organism_id, virus_hosts, l, g, ap, acts, bs, dnab, mb, nb, 
                #func, taxid, cr_biogrid, cr_dip, cr_intact, cr_mint, cr_string, cr_ensemble, cr_geneid, iw, 
                #taxlin, cr_pfam, cr_tirgrfam, cr_supfam, cr,interpro, cr_prosite, fr, sequence, creation_date, mod_date, 
                #seqmod_date, keywords) = 
                cols = line.rstrip('\n').split('\t')

                entry_id = cols[0]
                entry_name = cols[1]
                isolate_name = cols[9]
                isolate_taxid = int(cols[10])
                protein_name = cols[3]
                gene_name = cols[4] 
                date = cols[38]
                sequence = cols[37]
                hosts = cols[11]
                chain = cols[43]
                pfam = cols[31]

                # viral species comes from the NCBI Taxonomy file -- roll isolate up to species level
                if isolate_taxid in taxmap:
                    species = taxmap[isolate_taxid]
                else:
                    print "skipping " + str(isolate_taxid) + " because it doesn't look like a virus (maybe update viruses.categories.txt)"
                    next


                string_id = str(species) + "." + str(entry_name)


                # is this protein from a genome that is the reference or is complete?
                keywords = cols[41]
                ref = False
                comp = False
                if re.search('Reference proteome', keywords):
                    ref = True
                if re.search('Complete proteome', keywords):
                    comp = True

                polyprotein = False
                if re.search('olyprotein', cols[3]):
                    polyprotein = True

                this_protein = {'string_id': string_id, 'entry_id': entry_id, 'entry_name': entry_name, 'isolate_name': isolate_name, 'gene_name': gene_name, 'protein_name': protein_name, 'species': species, 'isolate_taxid': isolate_taxid, 'reference_proteome': ref, 'complete_proteome': comp, 'polyprotein': polyprotein, 'date': date, 'hosts': hosts, 'sequence': sequence, 'chain': chain, 'pfam': pfam}

                if species in allspecies:
                    # then we are appending this_protein to the array
                    allspecies[species].append(this_protein)
                else:
                    allspecies[species] = [this_protein]

            count += 1
        return allspecies
    except:
        raise

def parse_human_proteins(filename):
    human = {} # uniprot: string identifiers
    fh = open(filename, 'r')
    for line in fh:
        species, string, uniprot, name = line.rstrip("\n").split("\t")
        human[uniprot] = string
    return human
        


def print_all_isolate(outfile, uniprot):
    with open(outfile, 'w') as f:
        for species, isolate in uniprot.iteritems():
            sortedspecies = sorted(isolate, key=lambda x: (x['polyprotein'], -x['reference_proteome'], -x['complete_proteome'], x['date']))
            for ss in sortedspecies:
                f.write(str(species) + "\t" + str(ss['isolate_taxid']) + "\t" + str(ss['entry_name']) + "\t" + str(ss['polyprotein']) + "\t" + str(ss['reference_proteome']) + "\t" + str(ss['complete_proteome']) + "\t" + str(ss['date']) + "\n")


def print_isolate_species(outfile, chosen_proteins):
    with open(outfile, 'w') as f:
        for species in chosen_proteins:
            for protein in chosen_proteins[species]:
                f.write(str(species) + "\t" + str(protein['isolate_taxid']) + "\t" + str(protein['entry_name']) + "\t" + protein['source'] + "\n")
                #f.write(str(protein['isolate_taxid']) + "\t" + str(species) + "\t" + str(protein['entry_name']) + "\t" + str(protein['polyprotein']) + "\t" + str(protein['reference_proteome']) + "\t" + str(protein['complete_proteome']) + "\t" + str(protein['date']) + "\n")

def isvirus(speciestax, chosen_proteins):
    if int(speciestax) in chosen_proteins:
        return True
    else:
        return False

def ishost(virus, nonvirus, chosen_proteins):
    virus = int(virus)
    if virus in chosen_proteins:
        #array = chosen_proteins[virus][0]['hosts_taxids'].split(",")
        array = chosen_proteins[virus][0]['hosts_taxids']
        return str(nonvirus) in array
    else:
        return False

def parse_pfam(pfam_file):
    domains = {}
    try:
        with open(pfam_file, 'r') as f:
            for l in f.readlines():
                if l.startswith('#') or l == "\n":
                    continue
                string_id, align_start, align_end, env_start, env_end, hmm_accession, hmm_name, hmm_type, hmm_start, hmm_end, hmm_length, bit_score, e_value, significance, clan = l.rstrip("\n").split()
                species, entry_name = split_string(string_id)

                this = {'string_id': string_id,
                        'species': species,
                        'entry_name': entry_name,
                        'align_start': align_start,
                        'align_end': align_end,
                        'env_start': env_start,
                        'env_end': env_end,
                        'hmm_accession': hmm_accession,
                        'hmm_name': hmm_name,
                        'hmm_type': hmm_type,
                        'hmm_start': hmm_end,
                        'hmm_length': hmm_length,
                        'bit_score': bit_score,
                        'e_value': e_value,
                        'significance': significance,
                        'clan': clan,
                        }

                pfam_id = hmm_accession.split('.')[0]
                if string_id in domains:
                    domains[string_id][pfam_id] = this
                else:
                    domains[string_id] = {pfam_id: this}

    except:
        pass
    return domains


def update_pfam(uniprot, pfam_dir):
    for species, stuff in uniprot.iteritems():
        pfam_file = os.path.join(pfam_dir, str(species) + ".fa.pfam")
        domains = parse_pfam(pfam_file)
        if len(domains.keys()) > 0:
            acc = 0;
            for isolate in stuff:
                string_id = str(isolate['species']) + "." + isolate['entry_name']
                if string_id in domains:
                    uniprot[species][acc]['pfam'] = domains[string_id]
                acc+=1

    return uniprot


def split_string(string_id):
    allitems = string_id.split(".")
    taxid = allitems[0]
    proteinid = ".".join(allitems[1:])
    return taxid, proteinid

def parse_simap(filename, chosen_proteins, taxtree, taxlevel, taxid_to_name, make_symmetric):
    with gzip.open(filename) as f:
        simap = {}
        print_counter = 100000
        acc = 0
        for line in f:
            acc += 1
            if acc == print_counter:
                acc = 0
                print "processed 100000"

            stringA, stringB, bit_score, perc_ident, perc_sim, startA, endA, startB, endB = line.rstrip("\n").split("\t")

            taxA, proteinidA = split_string(stringA)
            taxB, proteinidB = split_string(stringB)

            taxA = int(taxA)
            taxB = int(taxB)

            if taxA in chosen_proteins:
                nameA = chosen_proteins[taxA][0]['isolate_name']

                # use abbreviated isolate names
                m = re.match("(.*?) \(", nameA)
                if m is None:
                    abrev = nameA
                else:
                    abrev = m.groups()[0]
                nameA = abrev

                genenameA = chosen_proteins[taxA][0]['gene_name']
                proteinnameA = "" #chosen_proteins[taxA][0]['protein_name']
            else:
                if taxA in taxid_to_name:
                    nameA = taxid_to_name[taxA]
                else:
                    nameA = str(taxA)
                genenameA = str(stringA)
                proteinnameB = str(stringA)

            if taxB in chosen_proteins:
                nameB = chosen_proteins[taxB][0]['isolate_name']

                #use abbreviated isolate names // now have 'isolate_abbrev' change to use this if this code is actually still used
                m = re.match("(.*?) \(", nameB)
                if m is None:
                    abrev = nameB
                else:
                    abrev = m.groups()[0]
                nameB = abrev
                
                genenameB = chosen_proteins[taxB][0]['gene_name']
                proteinnameB = "" #chosen_proteins[taxB][0]['protein_name']
            else:
                if taxB in taxid_to_name:
                    nameB = taxid_to_name[taxB]
                else:
                    nameB = str(taxB)
                genenameB = str(stringB)
                proteinnameB = str(stringB)

            isvirusA = isvirus(taxA, chosen_proteins)
            isvirusB = isvirus(taxB, chosen_proteins)

            ishostB = ishost(taxA, taxB, chosen_proteins)  # is B the host of A?
            ishostA = ishost(taxB, taxA, chosen_proteins)

            isself = False
            if taxA == taxB:
                isself = True 

            levelA = ncbi.which_level(taxA, taxtree)
            levelB = ncbi.which_level(taxB, taxtree)

            familyA = ncbi.which_family(taxA, taxtree, taxlevel)
            familyB = ncbi.which_family(taxB, taxtree, taxlevel)
            if familyA in taxid_to_name:
                familynameA = taxid_to_name[familyA]
            else:
                familynameA = str(familyA)
            if familyB in taxid_to_name:
                familynameB = taxid_to_name[familyB]
            else:
                familynameB = str(familyB)

            orderA = ncbi.which_order(taxA, taxtree, taxlevel)
            orderB = ncbi.which_order(taxB, taxtree, taxlevel)
            if orderA in taxid_to_name:
                ordernameA = taxid_to_name[orderA]
            else:
                ordernameA = str(orderA)
            if orderB in taxid_to_name:
                ordernameB = taxid_to_name[orderB]
            else:
                ordernameB = str(orderB)

            viralcladeA = ''
            viralcladeB = ''
            hostkingdomA = ''
            hostkingdomB = ''
            hostfamilyA = ''
            hostfamilyB = ''
            if isvirusA:
                viralcladeA = ncbi.which_viral_clade(taxA, taxtree)
                hostkingdomA = ncbi.host_kingdom(taxA, chosen_proteins, taxtree)
                hostfamilyA = ncbi.host_family(taxA, chosen_proteins, taxtree, taxlevel)
            if isvirusB:
                viralcladeB = ncbi.which_viral_clade(taxB, taxtree)
                hostkingdomB = ncbi.host_kingdom(taxB, chosen_proteins, taxtree)
                hostfamilyB = ncbi.host_family(taxB, chosen_proteins, taxtree, taxlevel)


            # create this scruture
            #virus <: virus
            #virus : host
            #virus : other

            forward =                     {'stringA': stringA, 
                                           'stringB': stringB,
                                           'taxA': taxA,
                                           'taxB': taxB,
                                           'proteinidA': proteinidA,
                                           'proteinidB': proteinidB,
                                           'isvirusA': isvirusA,
                                           'isvirusB': isvirusB,
                                           'nameA': nameA,
                                           'nameB': nameB,
                                           'genenameA': genenameA,
                                           'genenameB': genenameB,
                                           'proteinnameA': proteinnameA,
                                           'proteinnameB': proteinnameB,
                                           'familyA': familyA,
                                           'familyB': familyB,
                                           'familynameA': familynameA,
                                           'familynameB': familynameB,
                                           'orderA': orderA,
                                           'orderB': orderB,
                                           'ordernameA': ordernameA,
                                           'ordernameB': ordernameB,
                                           'levelA': levelA,
                                           'levelB': levelB,
                                           'viralcladeA': viralcladeA,
                                           'viralcladeB': viralcladeB,
                                           'hostfamilyA': hostfamilyA,
                                           'hostfamilyB': hostfamilyB,
                                           'hostkingdomA': hostkingdomA,
                                           'hostkingdomB': hostkingdomB,
                                           'ishostA': ishostA,
                                           'ishostB': ishostB,
                                           'isself': isself,
                                           'bit_score': bit_score,
                                           'perc_ident': perc_ident,
                                           'perc_sim': perc_sim,
                                           'startA': startA,
                                           'endA': endA,
                                           'startB': startB,
                                           'endB': endB,
                                          }

            backward =                    {'stringA': stringB, 
                                           'stringB': stringA,
                                           'taxA': taxB,
                                           'taxB': taxA,
                                           'proteinidA': proteinidB,
                                           'proteinidB': proteinidA,
                                           'isvirusA': isvirusB,
                                           'isvirusB': isvirusA,
                                           'nameA': nameB,
                                           'nameB': nameA,
                                           'genenameA': genenameB,
                                           'genenameB': genenameA,
                                           'proteinnameA': proteinnameB,
                                           'proteinnameB': proteinnameA,
                                           'familyA': familyB,
                                           'familyB': familyA,
                                           'familynameA': familynameB,
                                           'familynameB': familynameA,
                                           'orderA': orderB,
                                           'orderB': orderA,
                                           'ordernameA': ordernameB,
                                           'ordernameB': ordernameA,
                                           'levelA': levelB,
                                           'levelB': levelA,
                                           'viralcladeA': viralcladeB,
                                           'viralcladeB': viralcladeA,
                                           'hostfamilyA': hostfamilyB,
                                           'hostfamilyB': hostfamilyA,
                                           'hostkingdomA': hostkingdomB,
                                           'hostkingdomB': hostkingdomA,
                                           'ishostA': ishostB,
                                           'ishostB': ishostA,
                                           'isself': isself,
                                           'bit_score': bit_score,
                                           'perc_ident': perc_ident,
                                           'perc_sim': perc_sim,
                                           'startA': startB,
                                           'endA': endB,
                                           'startB': startA,
                                           'endB': endA,

                                          }

            string1 = stringA
            string2 = stringB


            if isvirusA:
                if isvirusB:
                    if make_symmetric:
                        thing = forward
                        symmetric = backward
                    else:
                        if taxA < taxB:
                            thing = forward
                        else:
                            thing = backward
                            string1 = stringB
                            string2 = stringA
                else:
                    thing = forward
            else:
                thing = backward
                string1 = stringB
                string2 = stringA
                    

            if string1 in simap:
                simap[string1][string2] = thing
            else:
                simap[string1] = {}
                simap[string1][string2] = thing

            if make_symmetric and isvirusA and isvirusB:
                if string2 in simap:
                    simap[string2][string1] = symmetric
                else:
                    simap[string2] = {}
                    simap[string2][string1] = symmetric



    print "done parsing simap"
    return simap
                    
def normalize_bitscore(simap):
    missingself = 999999999

    for stringA in simap:
        for stringB in simap[stringA]:

            num = float(simap[stringA][stringB]['bit_score']) 

            if stringA in simap[stringA]:
                den1 = float(simap[stringA][stringA]['bit_score'])
            else:
                print "missing self " + stringA
                den1 = missingself

            if stringB in simap:
                if stringB in simap[stringB]:
                    den2 = float(simap[stringB][stringB]['bit_score'])
                else:
                    print "missing self den2 " + stringB
                    den2 = missingself
            else:
                # non-virus
                den2 = missingself

            norm = num / min(den1, den2)
            simap[stringA][stringB]['bit_score_A'] = den1
            simap[stringA][stringB]['bit_score_B'] = den2
            simap[stringA][stringB]['bit_score_norm'] = norm
    return simap
                                                               
def print_data(out_filename, simap):
    with open(out_filename, "w") as f:
        f.write("stringA\tstringB\ttaxA\ttaxB\tproteinidA\tproteinidB\tisvirusA\tisvirusB\tnameA\tnameB\tgenenameA\tgenenameB\tproteinnameA\tproteinnameB\tfamilyA\tfamilyB\tfamilynameA\tfamilynameB\torderA\torderB\tordernameA\tordernameB\tlevelA\tlevelB\tviralcladeA\tviralcladeB\thostfamilyA\thostfamilyB\thostkingdomA\thostkingdomB\tishostA\tishostB\tisself\tbit_score\tbit_score_A\tbit_score_b\tbit_score_norm\tperc_ident\tperc_sim\n")
        for stringA in simap:
            for stringB in simap[stringA]:
                f.write(str(simap[stringA][stringB]['stringA']) + "\t")
                f.write(str(simap[stringA][stringB]['stringB']) + "\t")
                f.write(str(simap[stringA][stringB]['taxA']) + "\t")
                f.write(str(simap[stringA][stringB]['taxB']) + "\t")
                f.write(str(simap[stringA][stringB]['proteinidA']) + "\t")
                f.write(str(simap[stringA][stringB]['proteinidB']) + "\t")
                f.write(str(simap[stringA][stringB]['isvirusA']) + "\t")
                f.write(str(simap[stringA][stringB]['isvirusB']) + "\t")
                f.write(str(simap[stringA][stringB]['nameA']) + "\t")
                f.write(str(simap[stringA][stringB]['nameB']) + "\t")
                f.write(str(simap[stringA][stringB]['genenameA']) + "\t")
                f.write(str(simap[stringA][stringB]['genenameB']) + "\t")
                f.write(str(simap[stringA][stringB]['proteinnameA']) + "\t")
                f.write(str(simap[stringA][stringB]['proteinnameB']) + "\t")
                f.write(str(simap[stringA][stringB]['familyA']) + "\t")
                f.write(str(simap[stringA][stringB]['familyB']) + "\t")
                f.write(str(simap[stringA][stringB]['familynameA']) + "\t")
                f.write(str(simap[stringA][stringB]['familynameB']) + "\t")
                f.write(str(simap[stringA][stringB]['orderA']) + "\t")
                f.write(str(simap[stringA][stringB]['orderB']) + "\t")
                f.write(str(simap[stringA][stringB]['ordernameA']) + "\t")
                f.write(str(simap[stringA][stringB]['ordernameB']) + "\t")
                f.write(str(simap[stringA][stringB]['levelA']) + "\t")
                f.write(str(simap[stringA][stringB]['levelB']) + "\t")
                f.write(str(simap[stringA][stringB]['viralcladeA']) + "\t")
                f.write(str(simap[stringA][stringB]['viralcladeB']) + "\t")
                f.write(str(simap[stringA][stringB]['hostfamilyA']) + "\t")
                f.write(str(simap[stringA][stringB]['hostfamilyB']) + "\t")
                f.write(str(simap[stringA][stringB]['hostkingdomA']) + "\t")
                f.write(str(simap[stringA][stringB]['hostkingdomB']) + "\t")
                f.write(str(simap[stringA][stringB]['ishostA']) + "\t")
                f.write(str(simap[stringA][stringB]['ishostB']) + "\t")
                f.write(str(simap[stringA][stringB]['isself']) + "\t")
                f.write(str(simap[stringA][stringB]['bit_score']) + "\t")
                f.write(str(simap[stringA][stringB]['bit_score_A']) + "\t")
                f.write(str(simap[stringA][stringB]['bit_score_B']) + "\t")
                f.write(str(simap[stringA][stringB]['bit_score_norm']) + "\t")
                f.write(str(simap[stringA][stringB]['perc_ident']) + "\t")
                f.write(str(simap[stringA][stringB]['perc_sim']) + "\t")
                f.write("\n")

def remove_too_small_genomes(chosen_proteins):
    '''
    Remove any genomes that have less than two proteins
    '''

    to_delete = []

    for species, proteins in chosen_proteins.iteritems():
        if len(proteins) < 2:
            to_delete.append(species)

    # can't delete while iterating
    for species in to_delete: 
        del(chosen_proteins[species])

    return chosen_proteins


#def print_duplicates(chosen_proteins):
#    for species, isolates in chosen_proteins.iteritems():
#        for i in isolates:
#
#


def remove_duplicate_protein(chosen_proteins, species, key, value):
    print "removing protein for species " + str(species)
    edit = chosen_proteins[species]
    acc = 0
    possible = {} # acc: positive attributes
    for i in edit:
        if i[key] == value:
            if i['source'] == "Swiss-Prot":
                possible[acc] = 0
            else:
                possible[acc] = -1
        acc+=1

    found = False
    for acc, pos in possible.iteritems():
        if pos < 0:
            edit[acc] = {}
            found = True
            # if it's all crap, just remove it all
    if not found: # if everything is good, then just remove the first few. meh
        if int(species) == 211044:
            sys.stderr.write("not found " + str(acc) + " " + str(pos) + "\n")
        for i in range(0, len(possible) - 1): # don't remove the last element
            goner = sorted(possible.keys())[i] # hopefully this will pick the same protein if we need to do this multiple times
            edit[goner] = {}
            #TODO add the removed protein to the kept protein so that synonyms will come from both

    # need to keep the indicies constant between removing things.. there must be a better way to do this
    while {} in edit:
        edit.remove({})

    chosen_proteins[species] = edit
    return chosen_proteins

def dedup(chosen_proteins):
    found = 0
    to_remove = []
    blasted = 0
    samegenename = 0
    for species, isolate in chosen_proteins.iteritems():
        sequences = []
        gene_names = []

        for i in isolate:
            seq = i['sequence']
            gn = i['gene_name']
            entry_name = i['entry_name']

            if 'from_polyprotein' in i:
                # polyproteins are carefully deduped already, so skip them
                continue

            for isoform in i['isoforms']:
                # only one isoform should be kept, keep the one with the longest sequence
                seq_isoform = get_protein_by_entry_name(chosen_proteins, species, isoform)
                if seq_isoform:
                    if len(seq) > len(seq_isoform['seq']):
                        to_remove.append({'species': species, 'key': 'entry_name', 'value': i})

            # if sequences match, they are the same thing
            for s in sequences:
                if seq == s['seq']:
                    to_remove.append({'species': species, 'key': 'sequence', 'value': seq})
                    if int(species) == 211044:
                        sys.stderr.write("removing because sequence " + seq + "\n")
            
            # blank gene name is not useful
            if gn == '':
                continue

            # have same names, then blast sequences
            if gn in gene_names:
                for s in sequences:
                    if s['gn'] == gn:
                        samegenename += 1
                        if blast(s['seq'], seq, i['entry_name'], s['entry_name']):
                            blasted += 1
                            to_remove.append({'species': species, 'key': 'sequence', 'value': seq})
                            if int(species) == 211044:
                                sys.stderr.write("removing because blast " + gn + "\n\t" + seq + "\n\t" + s['seq'] + "\n")

            sequences.append({'seq': seq, 'gn': gn, 'entry_name': entry_name})
            gene_names.append(gn)

    for x in to_remove:
        chosen_proteins = remove_duplicate_protein(chosen_proteins, x['species'], x['key'], x['value'])

    sys.stderr.write("dups " + str(found) + "\n")
    sys.stderr.write("blast operations " + str(blasted) + "\n")
    return chosen_proteins



def choose_isolates(allspecies, lookup_file):
    chosen = {}

    for species, isolate in allspecies.iteritems(): # species is an int
        # Ideally we'll take the earliest reference proteome that is not a polyprotein
        # Prefer non-polyprotein, reference proteome, complete_proteome, earliest, alphabetical by entry name
        sortedspecies = sorted(isolate, key=lambda x: (x['polyprotein'], -x['reference_proteome'], -x['complete_proteome'], x['date']))
        chosen[species] = sortedspecies[0]

        #print "species\tisolate\tentry_name\tpolyprotein\treference\tcomplete\tdate"
        #for s in sortedspecies:
        #    print str(s['species']) + "\t" + str(s['isolate_taxid']) + "\t" + str(s['entry_name']) + "\t" + str(s['polyprotein']) + "\t" + str(s['reference_proteome']) + "\t" + str(s['complete_proteome']) + "\t" + str(s['date'])


    # if any species have specified isolates in the lookup_file, then use these isolates instead
    if lookup_file:
        with open(lookup_file, 'r') as f:
            for l in f.readlines():
                species, isolate, entry_name = l.rstrip("\n").split("\t")
                for test_isolate in allspecies[int(species)]:
                    if test_isolate['isolate_taxid'] == int(isolate) and test_isolate['entry_name'] == entry_name:
                        chosen[int(species)] = test_isolate
                        break

    return chosen


def get_chosen_proteins(allspecies, chosen):
    chosen_proteins = {}
    for species, isolate in allspecies.iteritems():
        for i in isolate: 
            if species in chosen:
                if i['isolate_taxid'] == chosen[species]['isolate_taxid']:
                    if species in chosen_proteins:
                        chosen_proteins[species].append(i)
                    else:
                        chosen_proteins[species] = [i]
                else:
                    continue
            else:
                print "can't find " + str(species) + " in chosen"
    return chosen_proteins

def statistics_all(uniprot):
    for species, isolate in uniprot.iteritems():
        taxids = []
        proteins = {}
        n = 0
        for i in isolate:
            n += 1
            taxids.append(i['isolate_taxid'])
            if i['gene_name'] in proteins:
                proteins[i['gene_name']] += 1
            else:
                proteins[i['gene_name']] = 1
        print
        print "Species: " + str(species)
        print "\tIsolates: " + str(n)
        print "\tUnique isolates: " + str(len(set(taxids))) + " " + str(set(taxids))
        print "\tProteins"
        for prot in proteins:
            print "\t\t" + str(proteins[prot]) + "\t--\t" + prot

def statistics_chosen(chosen):
    n_polyproteins = 0
    n = 0
    n_neither = 0
    for species, isolate in chosen.iteritems():
        n += 1
        if isolate['polyprotein'] == True:
            n_polyproteins += 1
            print str(isolate['species']) + " " + isolate['isolate_name']
        if not(isolate['complete_proteome'] or isolate['reference_proteome']):
            n_neither += 1
    print
    print "Statistics ----------------------"
    print "Total species: " + str(n)
    print "Polyproteins: " +  str(n_polyproteins)
    print "Not Ref or Complete: " + str(n_neither)
    print "---------------------------------"

def chains_longer(seen_chains, nc):
    for k, sc in seen_chains.iteritems():
        if int(sc['start']) <= int(nc['start']) and int(sc['end']) >= int(nc['end']):
            return sc['name']
    return False

def chains_shorter(seen_chains, nc):
    for k, sc in seen_chains.iteritems():
        if int(sc['start']) >= int(nc['start']) and int(sc['end']) <= int(nc['end']):
            return True
        if int(sc['start']) <= int(nc['start']) and int(sc['end']) == int(nc['end']):
            return True

def chains_overlap(seen_chains, nc):
    for k, sc in seen_chains.iteritems():
        if int(sc['start']) < int(nc['start']) and int(sc['end']) < int(nc['end']) and int(sc['end']) > int(nc['start']):
            return True

        # TODO remove all debugging crap

def cleave_polyproteins(chosen_proteins):
    ''' cleave polyproteins into functional units based on their chain entries
    '''
    # some viruses are just polyproteins; others are some proteins and a polyprotein; so have to check all proteins
    reported_dup_chain = []
    for species, proteins in chosen_proteins.iteritems():
        for i in proteins:
            if i['polyprotein']:
                #chains = i['chain'].split("; ")
                chains = i['chain']

                # build a map of which residues are split into subchains
                # use this to control which chains are output -- don't want any that contain duplicate proteins

                seen_chains = {}
                sc = {}

                for chain in chains:

                    # the old way of doing this, from the tsv is commented

                    # match >34 and ?98 but not ? as position
                    #m = re.match("CHAIN <?\??([0-9]+) >?\??([0-9]+) (.*?)\..*FTId=(.*).", chain)
                    #print chain
                    #if m is None:
                    #    continue
                    #else:
                        # add new proteins with correct sequence and id
                        #start = int(m.groups()[0])
                        #end = int(m.groups()[1])
                        #name = m.groups()[2]
                        #ftid = m.groups()[3]



                    start = int(chain['start'])
                    end = int(chain['end'])
                    name = chain['name']
                    ftid = chain['id']
                    aliases = chain['aliases']

                    short_name = ''
                    if 'short_name' in chain:
                        short_name = chain['short_name']
                        #print "------ got short name " + chain['short_name']

                    polyprotein = False
                    if re.search('olyprotein', name):
                        polyprotein = True
                    if chain['name'] == i['protein_name']:
                        polyprotein = True

                    this = {'start': start,
                            'end': end,
                            'name': name,
                            'ftid': ftid,
                            'aliases': aliases,
                            'polyprotein': polyprotein,
                            'short_name': short_name,
                            }

                    # select only fragments of polyproteins that contain unique genes
                    # eg if there is polyproteinABC that is cleaved into polyproteinAB and ppC, and ppAB is then cleaved to ppA and ppB, then just keep chains for ppC, ppA and ppB

                    # if there is something longer than this in seen_chains
                    co = chains_longer(seen_chains, this)
                    if co:
                        seen_chains.pop(co)
                        # delete the longer
                        seen_chains[name] = this
                    elif chains_shorter(seen_chains, this):
                        pass # a shorter chain that overlaps this one is already in there
                    elif chains_overlap(seen_chains, this):
                        print "OMG chains overlap, giving in " + i['isolate_name'] + " " + i['protein_name'] + " up on these" 
                        pass
                    else:
                        # just store this in seen_chains 
                        seen_chains[name] = this

                # then store the seen chains
                for k, chain in seen_chains.iteritems():

                    if not polyprotein:
                        for x in range(int(chain['start']), int(chain['end']) + 1):
                            if x in sc:
                                sc[x].append(chain['name'])
                                if i['entry_name'] not in reported_dup_chain:
                                    reported_dup_chain.append(i['entry_name'])
                            else:
                                sc[x] = [ chain['name'] ]
                    dups = 0
                    for k, v in sc.iteritems():
                        if len(v) > 1:
                            dups += 1
                    if dups > 5:
                        # allow a handful of aa overlaps
                        print " ** duplicate chain in " + str(species) + " " + i['entry_name'] + " " + str(dups)
                        print " ** " + str(sc)

                    this_protein = {}
                    # copy over all the existing keys
                    for key in i:
                        if key != 'short_name':
                            # only set the short name if it is in the chain (below)
                            this_protein[key] = i[key]
                    # then update the few that need to be changed
                    this_protein['polyprotein_entry_name'] = i['entry_name']
                    this_protein['entry_name'] = chain['ftid']
                    this_protein['string_id'] = str(this_protein['species']) + "." + str(this_protein['entry_name'])
                    this_protein['protein_name'] = chain['name']
                    this_protein['polyprotein'] = chain['polyprotein']
                    this_protein['sequence'] = i['sequence'][chain['start'] - 1 : chain['end']]
                    this_protein['chain'] = ""
                    this_protein['aliases'] = chain['aliases']
                    this_protein['from_polyprotein'] = True
                    if 'short_name' in chain and chain['short_name'] != '':
                        this_protein['short_name'] = chain['short_name']
                        #print "set chain short name " + this_protein['entry_name'] + " to " + chain['short_name']


                    chosen_proteins[species].append(this_protein)
                
            
            # for HIV-1 and HIV-2 gp-160 that need to be cleaved to gp120 and gp41
            elif i['entry_name'] == 'ENV_HV1H2' or i['entry_name'] == 'ENV_HV2BE':
                i['polyprotein'] = True
                for chain in i['chain']:
                    start = int(chain['start'])
                    end = int(chain['end'])
                    name = chain['name']
                    ftid = chain['id']
                    aliases = chain['aliases']

                    polyprotein = False

                    this_protein = {}
                    # copy over all the existing keys
                    for key in i:
                        this_protein[key] = i[key]
                    # then update the few that need to be changed
                    this_protein['polyprotein_entry_name'] = i['entry_name']
                    this_protein['entry_name'] = ftid
                    this_protein['string_id'] = str(this_protein['species']) + "." + str(this_protein['entry_name'])
                    this_protein['protein_name'] = name
                    this_protein['polyprotein'] = polyprotein
                    this_protein['sequence'] = i['sequence'][start - 1 : end]
                    this_protein['chain'] = ""
                    this_protein['aliases'] = aliases
                    this_protein['from_polyprotein'] = True

                    #print "******* Added " + i['entry_name'] + " " + ftid + " " + str(this_protein)

                    chosen_proteins[species].append(this_protein)


                
    return chosen_proteins

def map_pogs(pogs_file, out_mapped_pogs, chosen_proteins):
    # create dictionary mapping embl ids to my protein data structures
    embl_to_data = {}
    for species, proteins in chosen_proteins.iteritems():
        for p in proteins:
            for pid in p['embl_ids']:
                if pid in embl_to_data:
                    embl_to_data[pid].append(p['string_id'])
                else:
                    embl_to_data[pid] = [p['string_id']]

    with open(pogs_file, "r") as f:
        with open(out_mapped_pogs, "w") as out:
            for l in f.readlines():
                pog, pid = l.rstrip("\n").split("\t")
                if pid in embl_to_data:
                    #chosen_string_id = ''
                    #sortedspecies = sorted(embl_to_data[pid], key=lambda x: (x['polyprotein'], -x['reference_proteome'], -x['complete_proteome'], x['date']))
                    #for s in sortedspecies:
                    #    chosen_string_id = find_closest_chosen_protein(s, chosen_proteins)
                    #    if chosen_string_id != '':
                    #        continue

                    #if chosen_string_id == '':
                    #    print "can't find a match for protein id " + str(pid) + " when mapping pogs"
                    for string_id in embl_to_data[pid]:
                        out.write(pog + "\t" + str(string_id) + "\n")

def generate_clades(chosen_proteins, taxtree, taxlevel, stoplevel):
    ''' generate the initial clades that will feed into eggNOG
    '''
                #this_protein = {'string_id': string_id, 'entry_id': entry_id, 'entry_name': entry_name, 'isolate_name': isolate_name, 'gene_name': gene_name, 'protein_name': protein_name, 'species': species, 'isolate_taxid': isolate_taxid, 'reference_proteome': ref, 'complete_proteome': comp, 'polyprotein': polyprotein, 'date': date, 'hosts': hosts, 'sequence': sequence, 'chain': chain, 'pfam': pfam}

    clades = {}
    for species, proteins in chosen_proteins.iteritems():
        for i in proteins:
            clade = ncbi.climb_tax_tree_to_level(i['species'], taxtree, taxlevel, stoplevel)
            if int(clade) == 1: # if the virus does not have a family, then just put the species
                print "family is undef for " + str(i['species'])
                clade = i['species']
            species_taxid = split_string(i['string_id'])[0]

            if clade in clades:
                if species_taxid in clades[clade]:
                    clades[clade][species_taxid].append({'string_id': i['string_id'], 'sequence': i['sequence']})
                else:
                    clades[clade][species_taxid] = [{'string_id': i['string_id'], 'sequence': i['sequence']}]
            else:
                clades[clade] = {}
                clades[clade][species_taxid] = [{'string_id': i['string_id'], 'sequence': i['sequence']}]

    return clades

def print_clades(out_clades_file, clades):
    ''' print clades
    '''
    with open(out_clades_file, "w") as f:
        acc = 1
        for family, taxdata in clades.iteritems():
            clname = "CL" + '%(number)03d' % {"number": acc}
            taxids = " ".join(taxdata.keys())
            f.write(str(clname) + " " + str(taxids) + "\n")
            acc+=1

def print_clades_fasta(out_clades_dir, clades):
    ''' print clades
    '''
    for name, fam in clades.iteritems():
        with open(out_clades_dir + "/" + str(name) + ".fa", "w") as f:
            for data in fam.values():
                for entry in data:
                    f.write(">" + entry['string_id'] + "\n" + entry['sequence'] + "\n")

def parse_protein_sizes(protein_sizes_filename):
    sizes = {}
    with open(protein_sizes_filename, "r") as f:
        for line in f:
            string_id, size = line.rstrip("\n").split("\t")
            sizes[string_id] = size
    return sizes

def print_host_kingdom(out_host_kingdom_file, chosen_proteins, taxtree):
    with open(out_host_kingdom_file, "w") as f:
        for species, stuff in chosen_proteins.iteritems():
            f.write(str(species) + "\t" + ncbi.host_kingdom(species, chosen_proteins, taxtree) + "\n")

def print_homology(out_homology_filename, simap, sizes):
    with open(out_homology_filename, "w") as f:
        for sA in simap:
            for sB in simap[sA]:
                s = simap[sA][sB]
                string_id = s['stringB']
                if s['isvirusA'] and s['isvirusB']:
                    if string_id in sizes: # otherwise, it is one of the hits against non viral proteins
                        f.write(s['stringA'] + "\t" + s['stringB'] + "\t" + str(s['taxB']) + "\t" + str(s['bit_score']) + "\t" + str(s['startA']) + "\t" + str(s['endA']) + "\t" + str(s['startB']) + "\t" + str(s['endB']) + "\t" + str(sizes[s['stringB']]) + "\n")

def print_nodes(out_nodes_dir, chosen_proteins, version, no_polyproteins):
    ''' print out the nodes to use as a payload
    '''
    for species, isolate in chosen_proteins.iteritems():
        with open(out_nodes_dir + "/" + str(species) + "_v_" + str(version) + ".node", "w") as f:
            for protein in isolate:
                if no_polyproteins and protein['polyprotein']:
                    continue
                else:
                    f.write(str(protein['string_id']) + "\n")

def print_fasta(out_fasta_dir, chosen_proteins, name_details, no_polyproteins):
    ''' open a file handle named on the species and append to this file
    '''
    print "out_fasta_dir " + out_fasta_dir
    for species in chosen_proteins:
        with open(out_fasta_dir + "/" + str(species) + ".fa", "w") as f:
            for protein in sorted(chosen_proteins[species], key=lambda x: (x['entry_name'])):
                #TODO pre-filter for polyproteins
                if no_polyproteins and protein['polyprotein']:
                    continue
                else: 
                    if name_details:
                        f.write(">" + str(protein['string_id']) + "." + str(protein['gene_name']) + "." + str(protein['protein_name']).replace(" ", "_") + "\n")
                    else:
                        f.write(">" + str(protein['string_id']) + "\n")

                    f.write(protein['sequence'] + "\n")

def print_uniprot(out_uniprot_file, chosen):
    with open(out_uniprot_file, "w") as f:
        for species, isolate in chosen.iteritems():
            for protein in chosen[species]:
                if 'polyprotein_entry_name' in protein:
                    f.write(protein['string_id'] + "\t" + protein['polyprotein_entry_name'] + "\n")
                else:
                    f.write(protein['string_id'] + "\t" + protein['entry_name'] + "\t" + "\n")

def filter_basename(basename):
    name = basename
    m = re.match("(.*) serotype", basename) 
    if m is not None:
        name = m.groups()[0]
    m = re.match("(.*) subtype", basename) 
    if m is not None:
        name = m.groups()[0]
    m = re.match("(.*) group", basename)
    if m is not None:
        name = m.groups()[0]
    m = re.match("(.*) genotype", basename)
    if m is not None:
        name = m.groups()[0]
    return name

def print_species(out_species_file, chosen):
    with open(out_species_file, "w") as f:
        for species, isolate in chosen.iteritems():
            isolate = isolate[0]
            org = isolate['isolate_name']
            m = re.match("(.*?) \(", org)
            if m is None:
                basename = filter_basename(org)
                abrev = org
            else:
                basename = filter_basename(m.groups()[0])
                m = re.match(".*\((.*)\)$", org)
                if m is None:
                    abrev = org
                else:
                    if (re.match("\)", m.groups()[0]) or re.match("isolate", m.groups()[0]) or re.match("strain", m.groups()[0])):
                        abrev = ""
                    else:
                        abrev = m.groups()[0]

            if len(org) > 100:
                m = re.search("\(.*\)", org)
                #if m is None:
                org = org[0:99]
                #else:
                    #org = re.sub(" \([^)]*\)$", "", org)

            f.write(str(isolate['species']) + "\t" + org + "\t" + basename + "\t" + abrev + "\t" + "Uniprot" + "\t" + "Core" + "\t" + "unused" + "\n")

def add_synonym(species, protein, name, source, synonyms):
    # do not allow only numbers
    disallowed_regex = re.compile('^\d+$')
    if disallowed_regex.match(name):
        return synonyms

    disallowed_single_words = ("Uncharacterized", "Probable", "Non-structural", "Capsid", "Envelope", "Tegument", "Movement", "Nucleoprotein", "Replicase", "transmembrane", "Membrane", "Matrix", "egress", "Transcriptase", "membrane")
    if name in disallowed_single_words:
        return synonyms

    disallowed_longer = ("Uncharacterized protein", "Putative uncharacterized protein", "Uncharacterised protein", "Putative uncharacterised protein", "protein Uncharacterized", "protein Uncharacterised")
    if name in disallowed_longer:
        return synonyms

    # has passed all filters, so add it
    synonyms[species][protein][name].append(source)
    return synonyms

def extract_synonyms(chosen, species, protein, synonyms, postfix = ''):
    # pull all synonyms out of the uniprot entry

    synonyms = add_synonym(species, protein, chosen['entry_name'], 'UniProtKB-EN' + postfix, synonyms)

    for eid in chosen['entry_id']:
        synonyms = add_synonym(species, protein, eid, 'UniProtKB-EI' + postfix, synonyms)

    # write the gene name if it isn't empty
    if chosen['gene_name']:
        if not 'from_polyprotein' in chosen: # XXX 
            synonyms = add_synonym(species, protein, chosen['gene_name'], 'UniProtKB-GN' + postfix, synonyms)
            if len(chosen['gene_name']) < 3:
                synonyms = add_synonym(species, protein, 'protein ' + chosen['gene_name'], 'UniProtKB-GN' + postfix, synonyms)
                synonyms = add_synonym(species, protein, chosen['gene_name'] + ' protein', 'UniProtKB-GN' + postfix, synonyms)

    if 'short_name' in chosen and chosen['short_name']:
        synonyms = add_synonym(species, protein, chosen['short_name'], 'UniProtKB-SN' + postfix, synonyms)

    for alias in chosen['aliases']:
        synonyms = add_synonym(species, protein, alias, 'UniProtKB-A' + postfix, synonyms)
    
        # uniprot entry name for HCV NS3 is incorrect
        if re.match("NS3P", alias):
            synonyms = add_synonym(species, protein, 'NS3', 'UniProtKB-ARGH' + postfix, synonyms)
            synonyms = add_synonym(species, protein, 'Hepatitis C virus NS3', 'UniProtKB-ARGH' + postfix, synonyms)
    
    # for both proteins that come from polyproteins and those that don't
    if len(str(chosen['protein_name'])) > 0:
        # whole original name
        synonyms = add_synonym(species, protein, chosen['protein_name'], 'UniProtKB-PN' + postfix, synonyms)

    #if postfix:
        #if we are adding names from other proteins, bail out early to not add so much crap
    #    return synonyms

    # then try to clean up the name
    name = chosen['protein_name']
    origname = name

    # delete everything in [] because this contains the chain info and we are doing that better above
    name = re.sub(r'\[[^\]]*?\]', '', name)

    # remove "virion" and "viral" in the name -- of course it is?!
    name = re.sub(r' virion', '', name)
    name = re.sub(r'Virion ', '', name)
    name = re.sub(r' viral', '', name)
    name = re.sub(r'Viral ', '', name)
    if len(name) > 0 and origname != name:
        synonyms = add_synonym(species, protein, name, 'UniProtKB-PV' + postfix, synonyms)
    origname = name

    # remove "putative" -- research confirming this prediction isn't going to use this in the protein name
    name = re.sub(r'Putative ', '', name)
    if len(name) > 0 and origname != name:
        synonyms = add_synonym(species, protein, name, 'UniProtKB-PP' + postfix, synonyms)
    origname = name

    # find all the synonyms in parens
    parens = re.findall(r'\((.*?)\)', name)
    for p in parens:
        # items in brackets
        synonyms = add_synonym(species, protein, p, 'UniProtKB-PNB' + postfix, synonyms)
        if len(p) < 3:
            synonyms = add_synonym(species, protein, "protein " + p, 'UniProtKB-PNB' + postfix, synonyms)
            synonyms = add_synonym(species, protein, p + " protein", 'UniProtKB-PNB' + postfix, synonyms)

    # then delete (greedy) all the things in () because we included them just above
    name = re.sub(r'\([^\)]*?\)', '', name)
    name = name.strip()
    if len(name) > 0 and origname != name:
        # cleaned
        synonyms = add_synonym(species, protein, name, 'UniProtKB-PNC' + postfix, synonyms)
        if len(name) <= 4:
            # items in brackets
            synonyms = add_synonym(species, protein, "protein " + name, 'UniProtKB-PNC' + postfix, synonyms)
            synonyms = add_synonym(species, protein, name + " protein", 'UniProtKB-PNC' + postfix, synonyms)

    # protein A <-> A protein
    pr_match = re.findall(r'^(.* )?protein( .*)?$', name, re.IGNORECASE)
    for m in list(sum(pr_match, ())): # convert tuples to a list
        m = m.strip()
        if len(m) > 0 and m != origname:
            synonyms = add_synonym(species, protein, m, 'UniProtKB-PM' + postfix, synonyms)
            synonyms = add_synonym(species, protein, m + " protein", 'UniProtKB-PM' + postfix, synonyms)
            synonyms = add_synonym(species, protein, "protein " + m, 'UniProtKB-PM' + postfix, synonyms)

    # glycoprotein G1 can also be written as G1 glycoprotein
    gp_match = re.findall(r'(.* )?glycoprotein( .*)?', name, re.IGNORECASE)
    for m in list(sum(gp_match, ())):
        m = m.strip()
        if len(m) > 0 and m != origname:
            synonyms = add_synonym(species, protein, m, 'UniProtKB-GP' + postfix, synonyms)
            synonyms = add_synonym(species, protein, m + " glycoprotein", 'UniProtKB-GP' + postfix, synonyms)
            synonyms = add_synonym(species, protein, "glycoprotein " + m, 'UniProtKB-GP' + postfix, synonyms)
        
    # if the name contains a slash, then include both before the slash and after 
    # although sometimes the format is "bar something/otherthing foo" and we should write "bar something foo" and "bar otherthing foo", I can't distinguish where the boundaries are programatically
    slashes = name.split("/")
    for s in slashes:
        s = s.strip()
        if len(s) > 0 and s != origname:
            synonyms = add_synonym(species, protein, s, 'UniProtKB-PS' + postfix, synonyms)

    ssu = re.sub(r' small subunit', '', name)
    if len(ssu) > 0 and ssu != name:
        synonyms = add_synonym(species, protein, ssu, 'UniProtKB-PSSU' + postfix, synonyms)

    lsu = re.sub(r' large subunit', '', name)
    if len(lsu) > 0 and lsu != name:
        synonyms = add_synonym(species, protein, lsu, 'UniProtKB-PLSU' + postfix, synonyms)

    su = re.sub(r' subunit.*', '', name)
    if len(su) > 0 and su != name:
        synonyms = add_synonym(species, protein, su, 'UniProtKB-PSU' + postfix, synonyms)

    pc = re.sub(r' precursor.*', '', name)
    if len(pc) > 0 and pc != name:
        synonyms = add_synonym(species, protein, pc, 'UniProtKB-PPC' + postfix, synonyms)

    return synonyms

def simplify_name(best_synonym):
    regex = re.compile('\sprotein\s', re.IGNORECASE)
    best_synonym = regex.sub(' ', best_synonym)
    regex = re.compile('^protein\s', re.IGNORECASE)
    best_synonym = regex.sub('', best_synonym)
    regex = re.compile('\sprotein$', re.IGNORECASE)
    best_synonym = regex.sub('', best_synonym)

    regex = re.compile('Uncharacterized\s', re.IGNORECASE)
    best_synonym = regex.sub('', best_synonym)

    regex = re.compile('Putative\s', re.IGNORECASE)
    best_synonym = regex.sub('', best_synonym)

    return best_synonym

def print_synonyms(out_synonyms_file, out_best_synonyms_file, chosen_proteins, allspecies=[]):
    # synonyms of the species -- already covered in textmining code
    synonyms = collections.defaultdict(lambda : collections.defaultdict(lambda : collections.defaultdict(list))) # species -> protein -> synonym -> [sources, ...]

    species_list = []

    for species, isolate in chosen_proteins.iteritems():
        species_list.append(species)
        for chosen in isolate:
            # aliases come from:
            # protein name
            # entry name
            # gene name
            # protein aliases
            # entry names of all isolates


            if (chosen['polyprotein'] == True and 'from_polyprotein' not in chosen):
                # skip the parents of polyproteins
                continue

            protein = chosen['entry_name']
            synonyms = extract_synonyms(chosen, species, protein, synonyms)

    if allspecies:
        # write out fasta file of all proteins
        out_allspecies_fasta = "data_out/blast/allspecies.fa"
        with open(out_allspecies_fasta, "w") as f:
            for species in species_list:
                if species in allspecies:
                    for allprotein in allspecies[species]:
                        f.write(">" + str(allprotein['species']) + "." + allprotein['entry_name'] + "\n" + allprotein['sequence'] + "\n")


        results_file = "data_out/blast/allspecies.results"
        #blast_all(out_allspecies_fasta, results_file) # done offline
        matches = parse_blast_results(results_file)
        for target, data in matches.items():
            target_species, target_entry_name = target.split(".")
            target_species = int(target_species)
            target = {}
            for t in chosen_proteins[target_species]:
                if t['entry_name'] == target_entry_name:
                    target = t
            for query in data:
                query_species, query_entry_name = query.split(".")
                query_species = int(query_species)
                if query_species in allspecies:
                    for allprotein in allspecies[query_species]:
                        # go through all the proteins associated with this to find the right one
                        if allprotein['entry_name'] == query_entry_name:
                            if allprotein['isolate_taxid'] != t['isolate_taxid']:
                                # do not transfer synonyms between proteins in the same isolate
                                postfix = "-T"
                                synonyms = extract_synonyms(allprotein, target_species, target_entry_name, synonyms, postfix + "_" + query_entry_name)


          
    # go through the datastructure and write it
    with open(out_synonyms_file, "w") as f:
        with open(out_best_synonyms_file, "w") as b:
            for species, sphash in synonyms.items():
                all_proteins_for_species = {}
                for protein, prhash in sphash.items():
                    preferred_names = {'UniProtKB-GN': '', 'UniProtKB-PN': '', 'UniProtKB-EN': ''}
                    for synonym, sylist in prhash.items():
                        # store the preferred names
                        for source in sylist:
                            if source == 'UniProtKB-SN':
                                preferred_names['UniProtKB-SN'] = synonym
                            if source == 'UniProtKB-GN':
                                preferred_names['UniProtKB-GN'] = synonym
                            if source == 'UniProtKB-PN':
                                preferred_names['UniProtKB-PN'] = synonym
                            if source == 'UniProtKB-EN':
                                preferred_names['UniProtKB-EN'] = synonym

                        # write all the synonyms
                        sourcelist = " ".join(set(sylist))
                        f.write(str(species) + "\t" + protein + "\t" + synonym + "\t" + sourcelist + "\n")

                    # write the best synonym
                    best_synonym = ''
                    best_source = ''

                    found = False
                    while not found:
                        if 'UniProtKB-SN' in preferred_names and preferred_names['UniProtKB-SN'] != '':
                            best_synonym = preferred_names['UniProtKB-SN']
                            best_source = 'UniProtKB-SN'
                            if best_synonym not in all_proteins_for_species:
                                found = True
                            else:
                                best_synonym = ''

                        if best_synonym == '' and 'UniProtKB-GN' in preferred_names and preferred_names['UniProtKB-GN'] != '':
                            best_synonym = preferred_names['UniProtKB-GN']
                            best_source = 'UniProtKB-GN'
                            if best_synonym not in all_proteins_for_species:
                                found = True
                            else:
                                best_synonym = ''

                        if best_synonym == '' and 'UniProtKB-PN' in preferred_names and preferred_names['UniProtKB-PN'] != '':
                            best_synonym = preferred_names['UniProtKB-PN']
                            best_source = 'UniProtKB-PN'
                            if best_synonym not in all_proteins_for_species:
                                found = True
                            else:
                                best_synonym = ''

                        if best_synonym == '':
                            best_synonym = preferred_names['UniProtKB-EN'] # this always exists
                            best_source = 'UniProtKB-EN'
                            if best_synonym not in all_proteins_for_species:
                                found = True
                            else:
                                print "best sysnonym " + best_synonym + " already taken"
                                best_synonym = best_synonym + "_2"
                                found = True
                                

                    if not found:
                        print "best not found"

                    if best_synonym == '':
                        print "preferred names:"
                        print str(preferred_names)
                        print "** Best synonym empty ** this should not happen ** " + str(species) + " " + str(protein)

                    # if protein is in the name, we don't need it 
                    best_synonym = simplify_name(best_synonym)

                    # all of that doesn't always get it right
                    if protein == "PA_I34A1":
                        best_synonym = "PA"
                    if protein == "PAX_I34A1":
                        best_synonym = "PAX"

                    write = str(species) + "\t" + protein + "\t" + best_synonym + "\t" + best_source + "\n"
                    all_proteins_for_species[best_synonym] = write

                for prot in all_proteins_for_species.values():
                    b.write(prot)
                    f.write(prot) # add the best names into the whole synonyms file too

def print_functions(out_functions_file, chosen_proteins):
    '''
    print the functions of each protein for import to the database
    '''
    with open(out_functions_file, "w") as f:
        for species, proteins in chosen_proteins.iteritems():
            for p in proteins:
                f.write(str(species) + "\t" + p['entry_name'] + "\t" + p['function'] + "\t" + "UniProt_DE" + "\n")

def print_pdb(out_pdb_file, chosen_proteins):
    '''
    print pdb structures for each protein
    '''
    with open(out_pdb_file, "w") as f:
        for species, proteins in chosen_proteins.iteritems():
            for p in proteins:
                for structure in p['pdb']:
                    f.write(str(species) + "\t" + p['entry_name'] + "\t" + structure + "\n")

def parse_blast_results(results_file):
    matches = collections.defaultdict(list)
    with open(results_file) as f:
        for line in f.readlines():
            qseqid, sseqid, pident, qlen, slen, length, mismatch, gapopen, evalue, bitscore = line.split("\t")
            q_tax, q_entry_name = split_string(qseqid)
            s_tax, s_entry_name = split_string(sseqid)

            if q_tax == s_tax: 
                # only count the match if these roll up to the same to the same species
                if q_entry_name != s_entry_name:
                    # but don't count the entries we have already included
                    if pident > 95:
                        if (int(slen) - int(length)) < 0.05 * int(slen):
                            if (int(qlen) - int(length)) < 0.05 * int(qlen):
                                matches[sseqid].append(qseqid)
    return matches
                
def get_protein_by_entry_name(chosen_proteins, species, entry_name):
    if species in chosen_proteins:
        for possible in chosen_proteins[species]:
            if possible['entry_name'] == entry_name:
                return possible
    return None


def find_closest_chosen_protein_by_name(species, name, chosen_proteins):
    if species in chosen_proteins:
        for possible in chosen_proteins[species]:
            if possible['entry_name'] == name:
                return possible['string_id']
            for eid in possible['entry_id']:
                if eid == name:
                    return possible['string_id']
            if possible['protein_name'] == name:
                return possible['string_id']
    return ''

def find_closest_chosen_protein(protein_data, chosen_proteins):
    '''
    for a given protein that may not be from the chosen isolate, return the closest corresponding protein from the chosen isolate
    '''
    species = protein_data['species']
    #print "finding closest protein for " + str(protein_data)
    for possible in chosen_proteins[species]:
        if possible['sequence'] == protein_data['sequence']:
            #print "closest protein by sequence for protein_data is " + possible['string_id']
            #print possible
            return possible['string_id']
        if possible['protein_name'] == protein_data['protein_name']:
            #print "closest protein by name for protein_data is " + possible['string_id']
            #print possible
            return possible['string_id']
        #if possible['protein_name'] == protein_data['aliases'
        # TODO check a bunch of aliases

    #print "could not find a match for species " + str(species)
    #print protein_data
    return ''



def blast_all(allspecies_fasta):
    # blastp -query allspecies.fa -db viruses -outfmt '6 qseqid sseqid pident qlen slen length mismatch gapopen evalue bitscore' > allspecies.results
    #command1 = "blastp -query " + allspecies_fasta + " -db data_out/blast/viruses " + "-outfmt '6 qseqid sseqid pident qlen slen length mismatch gapopen evalue bitscore'"
    #process1 = Popen(shlex.split(command1), stdout=PIPE)
    #out, err = process1.communicate()
    #exit_code = process1.wait()
    return True

def blast(seq1, seq2, en1, en2):
    # run blast on these sequences and if they have above 95% identity, then consider them to be the same
    tmpdir = "/tmp"

    # this is not threadsafe, but we're only using one thread anyway
    tmp1 = tmpdir + "/tmp1.fna"
    tmp2 = tmpdir + "/tmp2.fna"
    with open(tmp1, 'w') as f1:
        f1.write(seq1 + "\n")
    with open(tmp2, 'w') as f2:
        f2.write(seq2 + "\n")

    command1 = "blastp -query " + tmp1 + " -subject " + tmp2 + " -outfmt '6 pident qlen qseq slen sseq'"
    process1 = Popen(shlex.split(command1), stdout=PIPE)
    command2 = "head -1"
    process2 = Popen(shlex.split(command2), stdin=process1.stdout, stdout=PIPE)
    out, err = process2.communicate()
    exit_code = process2.wait()

    if out:
        pident, qlen, qseq, slen, sseq = out.strip().split()
        #print "\t".join(("blastresults", pident, qlen, qseq, slen, sseq, en1, en2))
        if float(pident) > 95:  # 95% identity over 95% of one sequence
            if int(qlen)*.75 <= len(qseq) and int(slen)*.75 <= len(sseq):
                return True
        return False
    return False



def print_allisolates(out_allisolates_file, uniprot):
    with open(out_allisolates_file, "w") as f:
        for species in uniprot:
            for isolate in uniprot[species]:
                f.write(str(isolate['isolate_taxid']) + "\n")

               
def print_edges(out_edges_dir, edges, version):
    for species in edges:
        with open(out_edges_dir + "/" + str(species) + "_v_" + str(version) + ".edges", "w") as f:
            for edge in edges[species]:
                f.write(str(edge['stringA']) + "\t" + str(edge['stringB']) + "\t" + edge['evidence_type'] + "\t" + str(edge['interaction_score']) + "\t\t\n")

    
def print_payload(json_filename, species, version):

    nodefilename = "http://jensenlab.org/helen/hvppi/nodes/" + str(species) + "_v_" + str(version) + ".node"
    edgefilename = "http://jensenlab.org/helen/hvppi/edges/" + str(species) + "_v_" + str(version) + ".edges"

    json_data = { "extension_nodes_file": nodefilename,
                  "extension_edges_file": edgefilename,
                  "name": "Human-Virus PPI",
                }

    with open(json_filename, "w") as f:
        f.write(json.dumps(json_data) + "\n")

def print_indexhtml(html_filename, chosen_proteins, edges, version):
    with open(html_filename, "w") as f:
        # TODO use an html templating language?
        f.write("<html>\n")
        f.write("<body>\n")
        f.write("<script type=\"text/javascript\">\n")
        f.write("<!--\n")
        f.write("var flVer = GetFlashMajorVersion();\n")
        f.write("if(flVer>=10){ \n")

        f.write("$(\"a\").attr('href', function(i, h) {\n")
        f.write("return h + (h.indexOf('?') != -1 ? \"&flash=flVer\" : \"?flash=flVer\");\n")
        f.write("});\n")
        f.write("});\n")

        f.write(" -->\n")
        f.write("</script>\n")

        for species in sorted(chosen_proteins.keys()):
            if species in edges:
                string_id = chosen_proteins[species]['string_id']
                virus = chosen_proteins[species]['isolate_name']
                virus_protein = chosen_proteins[species]['entry_name']

                #print species
                #if species in edges:
                #    print edges[species]
                #print


                # we don't need to find the human node anymore
                #humanid = ""
                #for edge in edges[species]:
                    #if edge['sp_taxA'] == 9606:
                        #humanid = edge['stringA']
                        #break
                    #if edge['sp_taxB'] == 9606:
                        #humanid = edge['stringB']
                        #break

                #if humanid != "":
                f.write("<div id='stringlink'>\n")
                f.write("<a href=\"http://string-db.org/version_9_1/newstring_cgi/show_network_section.pl?targetmode=proteins&caller_identity=hvppi&network_flavor=evidence&additional_network_nodes=10&identifier=" + str(string_id) + "&input_query_species=9606&external_payload=http://jensenlab.org/helen/hvppi/hvppi_" + str(species) + "_v_" + str(version) + ".json\">")

                # TODO -- put the name directly in the structure
                    
                # to generate the image:
                # http://string-db.org/version_9_0/api/image/network?targetmode=proteins&caller_identity=hvppi&network_flavor=evidence&identifier=9606.ENSP00000328173
 
                f.write(str(species) + " " + virus + " " + virus_protein)
                f.write("</a>")
                f.write("</div>\n")



        f.write("</body>\n")
        f.write("</html>\n")

def print_pfam(out_pfam_file, chosen_proteins):
    no_polyproteins = True
    with open(out_pfam_file, "w") as f:
        for species, isolate in chosen_proteins.iteritems():
            for protein in isolate:
                if no_polyproteins and protein['polyprotein']:
                    continue
                else:
                    for domain in protein['pfam'].keys():
                        f.write(str(protein['string_id']) + "\t" + protein['entry_id'][0] + "\t" + str(domain) + "\n")

def print_go(out_go_file, chosen_proteins):
    no_polyproteins = True
    with open(out_go_file, "w") as f:
        for species, isolate in chosen_proteins.iteritems():
            for protein in isolate:
                if no_polyproteins and protein['polyprotein']:
                    continue
                else:
                    for go in protein['go']:
                        f.write(str(protein['string_id']) + "\t" + protein['entry_id'][0] + "\t" + str(go) + "\n")
                   
#def convert_host_name_to_taxids(chosen_proteins, taxnames):
#    for species in sorted(chosen_proteins.keys()):
#        hosts_names = chosen_proteins[species][0]['hosts'].split("; ")
#        hosts_taxid = []
#
#        for host in hosts_names:
#            if host in taxnames:
#                hosts_taxid.append(str(taxnames[host]))
#            else:
#                m = re.match(" ?([^()]*) [(].*", host)
#                if m:
#                    host = m.group(1)
#                    if host in taxnames:
#                        hosts_taxid.append(str(taxnames[host]))
#        chosen_proteins[species][0]['hosts_taxids'] = ",".join(hosts_taxid)
#    return chosen_proteins
        
def print_hosts(host_filename, chosen_proteins):
    with open(host_filename, "w") as f:
        for species in sorted(chosen_proteins.keys()):
            for host_taxid in chosen_proteins[species][0]['hosts_taxids']:
                f.write(str(species) + "\t" + str(host_taxid) + "\n")

 

#def print_hosts_old(host_filename, chosen_proteins, taxnames):
#    with open(host_filename, "w") as f:
#        for species in sorted(chosen_proteins.keys()):
#            hosts_names = chosen_proteins[species][0]['hosts'].split("; ") # just take the first protein
#            hosts_taxid = []
#
#            for host in hosts_names:
#                if host in taxnames:
#                    hosts_taxid.append(str(taxnames[host]))
#                    f.write(str(species) + "\t" + str(taxnames[host]) + "\n")
#                else:
#                    m = re.match(" ?([^()]*) [(].*", host)
#                    if m:
#                        host = m.group(1)
#                        if host in taxnames:
#                            hosts_taxid.append(str(taxnames[host]))
#                            f.write(str(species) + "\t" + str(taxnames[host]) + "\n")
#                        # just write the known hosts, don't bother with the ones we can't find (the're rare)
#                        #else:
#                            #hosts_taxid.append("")
#                    #else:
#                        #hosts_taxid.append("")
#            # write in long format instead of wide
#            #hosts_string = "\t".join(hosts_taxid)
#            #f.write(str(species) + "\t" + hosts_string + "\n")


def convert_edge_ids(edges):
    # turn uniprot into string ids for viruses
    # look up human mapping in the 9606 file
    # convert isolate to species taxids from intact
    return False

def clean_up(out_fasta_dir, out_species_file, out_synonyms_file, out_best_synonyms_file, out_host_file):
    if os.path.exists(out_species_file):
        os.remove(out_species_file)
    if os.path.exists(out_synonyms_file):
        os.remove(out_synonyms_file)
    if os.path.exists(out_best_synonyms_file):
        os.remove(out_best_synonyms_file)
    if os.path.exists(out_fasta_dir):
        shutil.rmtree(out_fasta_dir)
    if os.path.exists(out_host_file):
        os.remove(out_host_file)
    os.mkdir(out_fasta_dir)

def clean_up_dir(out_clades_dir):
    if os.path.exists(out_clades_dir):
        shutil.rmtree(out_clades_dir)
    os.mkdir(out_clades_dir)


def write_version(filename, version):
    with open(filename, "w") as f:
        f.write(str(version) + "\n")

def read_version(filename):
    fh = open(filename, 'r')
    old = 0
    for line in fh:
        old = int(line.rstrip("\n"))
    fh.close()
    return old

def increment_version(filename):
    if os.path.exists(filename):
        old = read_version(filename)
        new = old + 1
    else:
        new = 1

    write_version(filename, new)
    return new

def detect_fusion_proteins(simap):
    # requires that simap be made symmetric for viruses
    # compare each protein to each other protein
    fused = {}
    acc = 0
    for stringA in simap.keys():
        for stringB in simap.keys():
            acc += 1
            if acc == 100000:
                print "looked at 100000"
                acc = 0
            if not stringB in simap[stringA].keys():
                # stringA and stringB are not similar
                if simap[stringA].itervalues().next()['isvirusA']:
                    if simap[stringB].itervalues().next()['isvirusA']:
                        # and they are both virus proteins, so add them to the list
                        candidates = list(set(simap[stringA].keys()) & set(simap[stringB].keys()))
                        if stringA in fused:
                            fused[stringA][stringB] = candidates
                        else:
                            fused[stringA] = {}
                            fused[stringA][stringB] = candidates
    return fused

def print_fp(filename, fused):
    print "writing fusion proteins"
    with open(filename, "w") as f:
        f.write(str(fused))
    print "done writing fusion proteins"

