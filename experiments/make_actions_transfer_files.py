#!/usr/bin/env python
# encoding: utf-8

from __future__ import print_function

import sys
import os
import re

import glob
from collections import defaultdict
import compute_scores

base_dir = "/g/scb/bork/mering/string/data/primary_v10_5/evidence_transfer/"
import_chemicals = False

limited_transfer_channels = """BCRT
ECOC
NCIN
RCTM
go_textmining_cooccur
go_textmining_nlp
kegg
""".strip().split()


unknown_actions = defaultdict(int)

interaction_types = ("activation", "inhibition", "binding", "catalysis", "ptmod", "reaction", "expression", "phenotype")


translation = {
                        # mode, action (inhibition / activation)
    "complex"           : ("binding", ""),
    "complexassembly"   : ("binding", ""),
    "target"            : ("binding", ""),
    "direct"            : ("binding", ""),
    "direct_irrelevant" : ("binding", ""),
    "bind"              : ("binding", ""),
    "binding"           : ("binding", ""),
    "bndb"              : ("binding", ""),
    "gld"               : ("binding", ""),
    "kidb"              : ("binding", ""),
    "enzyme_subunit"    : ("binding", ""),
    "group_subunit"     : ("binding", ""),
    "modulation"        : ("binding", ""),

    "metab"             : ("catalysis", ""),
    "metabolizing"      : ("catalysis", ""),
    "catalysis"         : ("catalysis", ""),

    "biochemicalreaction" : ("reaction", ""),
    "reaction"          : ("reaction", ""),
    "transportwithbiochemicalreaction" : ("reaction", ""),

    "inh"               : ("inhibition", "inhibition"),
    "inhibit"           : ("inhibition", "inhibition"),
    "inhibition"        : ("inhibition", "inhibition"),
    "inhibition-allosteric" : ("inhibition", "inhibition"),
    "inhibition-competitive" : ("inhibition", "inhibition"),
    "inhibition-noncompetitive" : ("inhibition", "inhibition"),

    "act"               : ("activation", "activation"),
    "activate"          : ("activation", "activation"),
    "activation"        : ("activation", "activation"),
    "activation-allosteric" : ("activation", "activation"),

    "expr"              : ("expression", ""),
    "expression"        : ("expression", ""),

    "nacetyl"           : ("ptmod", ""),
    "acetyl"            : ("ptmod", ""),
    "deacetyl"          : ("ptmod", ""),
    "demethy"           : ("ptmod", ""),
    "dephos"            : ("ptmod", ""),
    "demethylation"     : ("ptmod", ""),
    "dephosphorylation" : ("ptmod", ""),
    "glyc"              : ("ptmod", ""),
    "methy"             : ("ptmod", ""),
    "methylation"       : ("ptmod", ""),
    "phos"              : ("ptmod", ""),
    "phosphorylation"   : ("ptmod", ""),
    "ubiq"              : ("ptmod", ""),
    "ubiquitination"    : ("ptmod", ""),

    "samemoa"           : ("phenotype", ""),
    "phenotype"         : ("phenotype", ""),
}

for (k, (mode, action)) in translation.items():
    assert mode in interaction_types, "unknown mode: %s for %s" % (mode, k)
    assert not action or action in interaction_types[:2], "unknown action: %s for %s" % (action, k)


def identify_interaction_type(action_description):
    """
    >>> identify_interaction_type("biochemicalReaction_lr")
    (('reaction', ''), 'lr')
    """

    action_description = action_description.lower()

    direction = ""
    specified_action = None
    interaction_type = None

    for part in action_description.split("_"):

        pl = part.lower()

        if len(part) == 2 and part in ("rx", "xr", "ll", "rr", "rl", "lr"):
            direction = part

        elif pl.startswith("act"):
            specified_action = "activation"

        elif pl.startswith("inh") or pl == "rep":
            specified_action = "inhibition"

        if interaction_type is None:
            interaction_type = translation.get(part)

    if interaction_type is None:
        interaction_type = translation.get(action_description)

    if specified_action is not None:
        (mode, _) = interaction_type

        if mode != "catalysis":
            interaction_type = (mode, specified_action)

    if interaction_type is None:
        unknown_actions[action_description] += 1

    return (interaction_type, direction)


def store_score(score, source, pmids, k, scores, use_max):

    if k in scores:

        (prev_score, prev_sources, prev_pmids) = scores[k]

        if use_max:

            if prev_score < score:
                scores[k] = (score, [source], pmids)

        elif all(pmid not in prev_pmids for pmid in pmids):

            if import_chemicals:
                score = compute_scores.combine_two_scores_protein_chemical(score, prev_score)
            else:
                score = compute_scores.combine_two_scores_protein_protein(score, prev_score)

            if source not in prev_sources:
                prev_sources.append(source)
            scores[k] = (score, prev_sources, prev_pmids + pmids)

    else:
        scores[k] = (score, [source], pmids)


def main():

    scores = {}
    action_scores = defaultdict(dict)

    fhs = {}

    for filename in ["../tm-string-virus-final/nlp/Medline_nlp_interact.tsv"]:
    #for filename in glob.glob(base_dir + "/channel_with_metadata/*interact.tsv"):

        print (filename, file=sys.stderr)

        #channel = re.search(r"([^/]+)_interact.tsv", filename).group(1)
        channel = "textmining_nlp"

        # no action from co-occurrence
        if "cooccur" in channel:
            continue

        if channel in ("textmining_cooccur", "go_textmining_cooccur"):
            use_max = False
        elif channel in ("experiments", "kegg", "RCTM", "BCRT", "ECOC", "NCIN", "databases", "predictions", "textmining_nlp", "go_textmining_nlp"):
            use_max = True
        else:
            sys.exit("Don't know channel: " + channel)

        #for line in os.popen("egrep -iv '^(pathwayneighbors|outer_metabolic_pathway|inner_metabolic_pathway)\t' "+filename):
        with open(filename, "r") as f:
            for line in f:

                fields = line.strip("\n").split("\t")

                if len(fields) == 6:
                    (action, species, item1, item2, score, source) = fields
                    pmid = ""
                elif len(fields) >= 9:
                    (action, species, item1, item2, score, pmid) = fields[:6]
                    source = channel
                elif len(fields) == 5:
                    (action, species, item1, item2, score) = fields
                    pmid = ""
                    source = channel
                else:
                    (action, species, item1, item2, score, pmid, source) = fields

                score = float(score)
                if score < 0.1:
                    continue

                (interaction_type, direction) = identify_interaction_type(action)

                if interaction_type is None:
                    continue

                pmids = []
                if pmid:
                    for s in re.split("[, ]+", pmid):
                        if s[0].isdigit():
                            s = "PMID%09d" % int(s)
                        pmids.append(s)

                if import_chemicals:
                    if not item2.startswith("CID"):
                        continue

                k = (species, item1, item2)

                path = channel if channel in limited_transfer_channels else ""

                store_score(score, source, pmids, k, action_scores[path, interaction_type, direction], use_max)


    for (path, (mode, specified_action), direction), scores in action_scores.items():

        path = "./actions/" + path
        filename = path + "/" + "_".join((mode, specified_action, direction)) + ".tsv"

        if not os.path.exists(path):
            os.mkdir(path)

        fh_out = open(filename, "w")

        for (species, item1, item2), (score, sources, pmids) in scores.iteritems():
            if pmids:
                print(species, item1, item2, score, " ".join(pmids), " ".join(sources), sep="\t", file=fh_out)
            else:
                print(species, item1, item2, score, " ".join(sources), sep="\t", file=fh_out)


    for action_description, count in unknown_actions.items():
        print("Unknown: ", action_description, count, file=sys.stderr)








if __name__ == '__main__':
    main()
