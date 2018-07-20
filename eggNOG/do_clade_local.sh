#! /bin/sh

# Call this script with a parameters for clade (example fungi)
# ./do_clade.sh 4751

echo $1

if [ $# -lt 1 ]
then
    echo "This script needs at least 1 arguments of the clade to compute (e.g. 4751)"
    exit 0;
fi

## ok, first we need to set some environment variables:
## paths to some data, including the protein sequences and the name of a postgres
## database ready to hold the Smith-Waterman results in tables.
## That database must already exist and have the necessary permissions.

#DERIVED_DIR=/bork/emil1/mering/STRING_derived_v8_1; export DERIVED_DIR
#FREEZE_DIR=/bork/emil1/mering/STRING_freeze_v8_1; export FREEZE_DIR

## database settings
#STRING_DATABASE=string_8_1; export STRING_DATABASE
#STRING_DATABASE_HOST=130.60.93.58; export STRING_DATABASE_HOST
#STRING_DATABASE_PORT=5432; export STRING_DATABASE_PORT
#TABLE_NAME_PROTEIN_SIMILARITIES=homology.blast_data; export TABLE_NAME_PROTEIN_SIMILARITIES

STRING_DATABASE=string_9_1
STRING_DATABASE_HOST=127.0.0.1
STRING_DATABASE_PORT=5432
TABLE_NAME_PROTEIN_SIMILARITIES=homology.blast_data
export STRING_DATABASE
export STRING_DATABASE_HOST
export STRING_DATABASE_PORT
export TABLE_NAME_PROTEIN_SIMILARITIES


## setteing for running script
#WORK_DIR=/local/eris/alexande/tmp/eggNOG/$1; export WORK_DIR     # the clade to be computed
#DATA_DIR=/local/eris/alexande/tmp/eggNOG/DATA; export DATA_DIR   # data common to all
#APPS_DIR=/local/eris/alexande/eggNOG; export APPS_DIR            # directory containing scripts
#PERL5LIB=/local/eris/alexande/eggNOG; export PERL5LIB            # directory containing perl modules

WORK_DIR=/home/stitch1/STITCH/eggNOG-pipeline/$1             # the clade to be computed
DATA_DIR=/home/stitch1/STITCH/eggNOG-pipeline/DATA            # data common to all
APPS_DIR=/home/stitch1/STITCH/eggNOG-pipeline                           # directory containing scripts
PERL5LIB=/home/stitch1/STITCH/eggNOG-pipeline                           # directory containing perl modules

export WORK_DIR
export DATA_DIR
export APPS_DIR
export PERL5LIB


CLADE=$1; export CLADE

if [ -f $WORK_DIR ]
then
    echo "$WORK_DIR already exists"
    #exit 0;
else
    mkdir $WORK_DIR
fi

$APPS_DIR/parse_NCBI_taxonomy_and_identify_groupings.pl \
    -tdf $DATA_DIR/nodes.dmp \
    -coif $DATA_DIR/clades_of_interest.txt \
    > $WORK_DIR/dependencies_clades_species.tsv;

cat $WORK_DIR/dependencies_clades_species.tsv | perl -ane 'print "$F[1]\n" if $F[0] == $ENV{CLADE}' \
    > $WORK_DIR/species.local.txt;

cat $WORK_DIR/species.local.txt | $APPS_DIR/prune_species_clades_file.pl \
    > $WORK_DIR/species_and_clades.txt;

$APPS_DIR/prune_species_lookup_file.pl \
    > $WORK_DIR/species.lookup;






$APPS_DIR/get_best_hits.pl \
    -k \
    -d $STRING_DATABASE \
    > $WORK_DIR/best_hits.txt;
















$APPS_DIR/find_inparalogs.pl \
    -d $STRING_DATABASE \
    -b $WORK_DIR/best_hits.txt \
    -h 550 \
    -l 50 \
    -s 100 \
    -m 50 \
    > $WORK_DIR/inparalogs.txt;

cat $WORK_DIR/best_hits.txt | $APPS_DIR/reduce_best_hits_via_inparalogs.pl $WORK_DIR/inparalogs.txt \
    > $WORK_DIR/best_hits_inparalogs.txt;


$APPS_DIR/find_reciprocal_best_hits.pl \
    -i $WORK_DIR/best_hits_inparalogs.txt \
    > $WORK_DIR/reciprocal_best_hits.txt; 

$APPS_DIR/find_and_extend_triangles.pl \
    -r $WORK_DIR/reciprocal_best_hits.txt \
    -sdf $DATA_DIR/species.taxonomic_domain.txt \
    -e off \
    -h 280 \
    -l 0 \
    -s 40 \
    > $WORK_DIR/og1_triangl.txt;

cat $WORK_DIR/og1_triangl.txt | $APPS_DIR/find_and_extend_tuples.pl \
    -r $WORK_DIR/reciprocal_best_hits.txt \
    -sdf $DATA_DIR/species.taxonomic_domain.txt \
    -e off \
    -h 260 \
    -l 60 \
    -s 40 \
    > $WORK_DIR/og1_triangl_tupl.txt; 

cat $WORK_DIR/og1_triangl_tupl.txt | $APPS_DIR/expand_with_inparalogs.pl \
    > $WORK_DIR/og1_triangl_tupl_expanded.txt; 

$APPS_DIR/find_reciprocal_best_hits.pl \
    -i $WORK_DIR/best_hits.txt \
    -c 300000 \
    > $WORK_DIR/reciprocal_best_hits_expanded.txt; 

$APPS_DIR/consolidate_groups.pl \
    -i $WORK_DIR/og1_triangl_tupl_expanded.txt \
    -d $STRING_DATABASE \
    > $WORK_DIR/og2.txt;

cat $WORK_DIR/og2.txt | $APPS_DIR/find_and_extend_triangles.pl \
    -r $WORK_DIR/reciprocal_best_hits_expanded.txt \
    -sdf $DATA_DIR/species.taxonomic_domain.txt \
    -e on \
    -d $STRING_DATABASE \
    -h 280 \
    -l 0 \
    -s 40 \
    > $WORK_DIR/og2_triangl.txt;

cat $WORK_DIR/og2_triangl.txt | $APPS_DIR/find_and_extend_tuples.pl \
    -r $WORK_DIR/reciprocal_best_hits_expanded.txt \
    -sdf $DATA_DIR/species.taxonomic_domain.txt \
    -e on \
    -d $STRING_DATABASE \
    -h 260 \
    -l 60 \
    -s 40 \
    > $WORK_DIR/og2_triangl_tupl.txt;

cat $WORK_DIR/og2_triangl_tupl.txt | $APPS_DIR/merge_orthgroups.pl \
    -r $WORK_DIR/reciprocal_best_hits_expanded.txt \
    -c 250 \
    -w \
    -t2 \
    > $WORK_DIR/og2_triangl_tupl_merged.txt;

$APPS_DIR/consolidate_groups.pl \
    -i $WORK_DIR/og2_triangl_tupl_merged.txt \
    -d $STRING_DATABASE \
    -bc -bk -bt \
    > $WORK_DIR/og3.txt;

cat $WORK_DIR/og3.txt | $APPS_DIR/find_and_extend_triangles.pl \
    -r $WORK_DIR/reciprocal_best_hits_expanded.txt \
    -sdf $DATA_DIR/species.taxonomic_domain.txt \
    -e on \
    -d $STRING_DATABASE \
    -h 280 \
    -l 0 \
    -s 40 \
    > $WORK_DIR/og3_triangl.txt;

cat $WORK_DIR/og3_triangl.txt | $APPS_DIR/find_and_extend_tuples.pl \
    -r $WORK_DIR/reciprocal_best_hits_expanded.txt \
    -sdf $DATA_DIR/species.taxonomic_domain.txt \
    -e on \
    -d $STRING_DATABASE \
    -h 260 \
    -l 60 \
    -s 40 \
    > $WORK_DIR/og3_triangl_tupl.txt;

$APPS_DIR/summarize_orthgroups.pl \
    -i $WORK_DIR/og3_triangl_tupl.txt \
    -d $STRING_DATABASE \
    -s 200 \
    > $WORK_DIR/og3_triangl_tupl_alignments.txt;

$APPS_DIR/assign_orthgroup_protein_positions.pl \
    -d $STRING_DATABASE \
    -o $WORK_DIR/og3_triangl_tupl.txt \
    > $WORK_DIR/og3_triangl_tupl_orthgroups.txt;

$APPS_DIR/get_best_hits.pl \
    -o $WORK_DIR/og3_triangl_tupl.txt \
    -d $STRING_DATABASE \
    > $WORK_DIR/best_hits_outside_orthgroup.txt;

$APPS_DIR/find_reciprocal_best_hits.pl \
    -i $WORK_DIR/best_hits_outside_orthgroup.txt \
    -c 300000 \
    > $WORK_DIR/recipr_best_hits_outside_orthgroups.txt;

$APPS_DIR/find_fusion_candidates.pl \
    -o $WORK_DIR/og3_triangl_tupl.txt \
    -r $WORK_DIR/recipr_best_hits_outside_orthgroups.txt \
    -c 80 \
    -t1 \
    > $WORK_DIR/node_fusion_candidates.txt;

$APPS_DIR/execute_fusion_proteins.pl \
    -o $WORK_DIR/og3_triangl_tupl_orthgroups.txt \
    -c $WORK_DIR/node_fusion_candidates.txt \
    -d $STRING_DATABASE \
    -nmc 35 \
    -nmf 0.5 \
    -mgs 10 \
    > $WORK_DIR/og3_triangl_tupl_fused.txt;

$APPS_DIR/consolidate_groups.pl \
    -i $WORK_DIR/og3_triangl_tupl_fused.txt \
    -d $STRING_DATABASE \
    -bc -bk -bt \
    > $WORK_DIR/og4.txt;

$APPS_DIR/assign_orphans_via_Cognitor.pl \
    -o $WORK_DIR/og4.txt \
    -d $STRING_DATABASE \
    -n 2 \
    -s 65 \
    -c \
    > $WORK_DIR/og6.txt;

cat $WORK_DIR/og6.txt | $APPS_DIR/orthgroup_counter.pl \
    > $WORK_DIR/final_count.txt

cat $WORK_DIR/og4.txt | $APPS_DIR/orthgroup_counter.pl \
    > $WORK_DIR/final_count.stringent.txt;

$APPS_DIR/summarize_orthgroups.pl \
    -i $WORK_DIR/og6.txt \
    -d $STRING_DATABASE \
    -s 100 \
    > $WORK_DIR/final_alignments.txt;n

$APPS_DIR/summarize_orthgroups.pl \
    -i $WORK_DIR/og4.txt \
    -d $STRING_DATABASE \
    -s 100 \
    > $WORK_DIR/final_alignments.stringent.txt;

$APPS_DIR/assign_orthgroup_protein_positions.pl \
    -d $STRING_DATABASE \
    -o $WORK_DIR/og6.txt \
    > $WORK_DIR/new_orthgroups.txt;

$APPS_DIR/assign_orthgroup_protein_positions.pl \
    -d $STRING_DATABASE \
    -o $WORK_DIR/og4.txt \
    > $WORK_DIR/new_orthgroups.stringent.txt;

