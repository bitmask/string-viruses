#! /bin/sh

## ok, first we need to set some environment variables:
## paths to some data, including the protein sequences and the name of a postgres 
## database ready to hold the Smith-Waterman results in tables.
## That database must already exist and have the necessary permissions.


# do these set up steps manually, then run the remainder 

# Inputs to the eggNOG pipeline
# 

# copy string/derived/ to eggNOG-pipeline/DATA
#scp species.names.txt proteins.checksums.txt protein_sizes.txt species.taxonomic_domain.txt protein.shorthands.txt clades_and_names.txt species.lookup helen@stitch:/home/stitch1/STITCH/eggNOG-pipeline/DATA
# make sure that in species.taxonomic_domain.txt col 2 is viruses
#ln -s species.lookup species.lookup.local

# copy stuff from data_out to eggNOG-pipeline/DATA
#scp -r fastafiles helen@stitch:/home/stitch1/STITCH/eggNOG-pipeline/DATA
#scp string_v10.species.tsv helen@stitch:/home/stitch1/STITCH/eggNOG-pipeline/DATA

# generate correct species_and_clades -- treat each virus as it's own clade
# for v5, do not do this, string steps have generated a reasonable clades file, so use it
#cat species.names.txt|grep -v \# |awk ' {printf("CL%03i",NR); print "\t" $1}' > species_and_clades.txt

# get nodes.dmp from ncbi taxonomy


#the following database tables need to be populated:

#  homology.blast_data
#the following can be populated, but it seems to not matter for eggNOG
#  items.species
#  items.proteins
#  items.proteins_sequences
#populating the database is done in the STRING doall_local.bash script

# pretend to have a few more files
#touch cogs+kogs.tsv
#cat data_out/pogs |uniq|awk '{print $2 "\t" $3 "\tstart\tend\t" $1 "\tannot"}' > koonin_ogs.tsv
#ln -s species.lookup species.lookup.local
#ln -s species_and_clades.txt species_and_clades.original.txt
#ln -s species.lookup species.lookup.original



# run this ./do_all_local.sh
# then, generate the sub groups  with ./do_clade_local.sh
# when this is done, take og3_triangl_tupl_orthgroups.txt from each of the levels in eggNOG/clades-indiv-xml
#for i in `cat ../clades_of_interest.txt|grep -v "#"`; do cp ~/mnt/eggnog/$i/og3_triangl_tupl_orthgroups.txt ./$i.out ; done


# start database
#/home/stitch1/STITCH/20150123-install/bin/postgres -D /home/stitch1/STITCH/20150123-install/pgsql/data >/home/stitch1/STITCH/20150123-install/log/postgres.log 2>&1 &
# connect to the database with the following
# /home/stitch1/STITCH/20150123-install/bin/psql -h 127.0.0.1 -p 5432 -U helen -d string_9_1

#DERIVED_DIR=/bork/emil1/mering/STRING_derived_v8_1
#FREEZE_DIR=/bork/emil1/mering/STRING_freeze_v8_1
STRING_DATABASE=string_11
STRING_DATABASE_HOST=127.0.0.1
STRING_DATABASE_PORT=5432
TABLE_NAME_PROTEIN_SIMILARITIES=homology.blast_data
WORK_DIR=/home/stitch1/STITCH/eggNOG-pipeline/all_xml             # the clade to be computed
DATA_DIR=/home/stitch1/STITCH/eggNOG-pipeline/DATA            # data common to all
APPS_DIR=/home/stitch1/STITCH/eggNOG-pipeline                           # directory containing scripts
PERL5LIB=/home/stitch1/STITCH/eggNOG-pipeline                           # directory containing perl modules

export WORK_DIR=/home/people/hecook/projects/eggNOG-pipeline_5/all
export DATA_DIR=/home/people/hecook/projects/eggNOG-pipeline_5/DATA
export APPS_DIR=/home/people/hecook/projects/eggNOG-pipeline_5/
export PERL5LIB=/home/people/hecook/projects/eggNOG-pipeline_5/

export WORK_DIR=/home/stitch1/STITCH/eggNOG-pipeline/all_xml             # the clade to be computed
export DATA_DIR=/home/stitch1/STITCH/eggNOG-pipeline/DATA            # data common to all
export APPS_DIR=/home/stitch1/STITCH/eggNOG-pipeline                           # directory containing scripts
export PERL5LIB=/home/stitch1/STITCH/eggNOG-pipeline                           # directory containing perl modules
export STRING_DATABASE=string_9_1
export STRING_DATABASE_HOST=127.0.0.1
export STRING_DATABASE_PORT=5432
export TABLE_NAME_PROTEIN_SIMILARITIES=homology.blast_data

#export DERIVED_DIR
#export FREEZE_DIR
export STRING_DATABASE
export STRING_DATABASE_HOST
export STRING_DATABASE_PORT
export TABLE_NAME_PROTEIN_SIMILARITIES
export WORK_DIR
export DATA_DIR
export APPS_DIR
export PERL5LIB

## first, we compute the best hits of each protein.


$APPS_DIR/get_best_hits.pl \
    -k \
    -d $STRING_DATABASE \
    > $WORK_DIR/best_hits_all_proteins.txt


# run this manually
# generate the initial clade assignments; writes initial_clade_assignments; rename this species_and_clades.txt
# python generate_clades_for_eggNOG.py

## then, we exclude everything that Lars has already assigned previously to COG or KOGs 
## and do the first round of NOG-assignments.

$APPS_DIR/get_COG_assignments_by_lars.pl \
    > $WORK_DIR/COG_groups.txt
$APPS_DIR/make_COG_exclusion_file.pl \
    -i $WORK_DIR/COG_groups.txt \
    > $WORK_DIR/COG_assigned_proteins.txt


#cat $DATA_DIR/species.lookup | $APPS_DIR/reduce_species_lookup_file.pl $WORK_DIR/COG_assigned_proteins.txt \
#    > $WORK_DIR/species.lookup.local

$APPS_DIR/find_inparalogs.pl \
    -d $STRING_DATABASE \
    -b $WORK_DIR/best_hits_all_proteins.txt \
    -h 550 \
    -l 50 \
    -s 100 \
    -m 50 \
    > $WORK_DIR/inparalogs.txt

    
# dependencies -- database input for get_best_hits.pl

$APPS_DIR/get_best_hits.pl \
    -k \
    -d $STRING_DATABASE \
    > $WORK_DIR/best_hits.txt

cat $WORK_DIR/best_hits.txt | $APPS_DIR/reduce_best_hits_via_inparalogs.pl $WORK_DIR/inparalogs.txt \
    > $WORK_DIR/best_hits_inparalogs.txt

$APPS_DIR/find_reciprocal_best_hits.pl \
    -i $WORK_DIR/best_hits_inparalogs.txt \
    > $WORK_DIR/reciprocal_best_hits.txt

# ok to here, but the next file is empty

# $APPS_DIR/find_and_extend_triangles.pl -r $WORK_DIR/reciprocal_best_hits.txt -e off -h 280 -l 0 -s 40 -sdf $DATA_DIR/species.taxonomic_domain.txt > $WORK_DIR/og1_triangl.txt
echo '#' | $APPS_DIR/find_and_extend_triangles.pl \
    -r $WORK_DIR/reciprocal_best_hits.txt \
    -sdf $DATA_DIR/species.taxonomic_domain.txt \
    -e off \
    -h 280 \
    -l 0 \
    -s 40 \
    > $WORK_DIR/og1_triangl.txt

cat $WORK_DIR/og1_triangl.txt | $APPS_DIR/find_and_extend_tuples.pl \
    -r $WORK_DIR/reciprocal_best_hits.txt \
    -sdf $DATA_DIR/species.taxonomic_domain.txt \
    -e off \
    -h 260 \
    -l 60 \
    -s 40 \
    > $WORK_DIR/og1_triangl_tupl.txt

cat $WORK_DIR/og1_triangl_tupl.txt | $APPS_DIR/expand_with_inparalogs.pl \
    > $WORK_DIR/og1_triangl_tupl_expanded.txt

$APPS_DIR/consolidate_groups.pl \
    -i $WORK_DIR/og1_triangl_tupl_expanded.txt \
    -d $STRING_DATABASE \
    > $WORK_DIR/og2.txt

## now, we join both assignments, and continue from here on with both in the pot.

cat $WORK_DIR/COG_groups.txt $WORK_DIR/og2.txt \
    > $WORK_DIR/og2_joined.txt

## now, move the local species.lookup out of the way (this way, the official will get used again from now on.

mv $WORK_DIR/species.lookup.local $WORK_DIR/species.lookup.first_steps_only

$APPS_DIR/find_reciprocal_best_hits.pl \
    -i $WORK_DIR/best_hits_all_proteins.txt \
    -c 1000000 \
    > $WORK_DIR/reciprocal_best_hits_expanded.txt

echo "*****************************************"
echo "og2_joined > find_and_extend_triangles > og2_triangl.txt"
echo "*****************************************"

cat $WORK_DIR/og2_joined.txt | $APPS_DIR/find_and_extend_triangles.pl \
    -r $WORK_DIR/reciprocal_best_hits_expanded.txt \
    -sdf $DATA_DIR/species.taxonomic_domain.txt \
    -e on \
    -d $STRING_DATABASE \
    -h 280 \
    -l 0 \
    -s 40 \
    > $WORK_DIR/og2_triangl.txt

echo "*****************************************"
echo "og2_triangl > find_and_extend_tuples > og2_triangl_tupl.txt"
echo "*****************************************"

cat $WORK_DIR/og2_triangl.txt | $APPS_DIR/find_and_extend_tuples.pl \
    -r $WORK_DIR/reciprocal_best_hits_expanded.txt \
    -sdf $DATA_DIR/species.taxonomic_domain.txt \
    -e on \
    -d $STRING_DATABASE \
    -h 260 \
    -l 60 \
    -s 40 \
    > $WORK_DIR/og2_triangl_tupl.txt

echo "*****************************************"
echo "merge_orthgroups > og2_triangl_tupl_merged.txt"
echo "*****************************************"

cat $WORK_DIR/og2_triangl_tupl.txt | $APPS_DIR/merge_orthgroups.pl \
    -r $WORK_DIR/reciprocal_best_hits_expanded.txt \
    -c 200 \
    -w \
    -t2 \
    > $WORK_DIR/og2_triangl_tupl_merged.txt

echo "*****************************************"
echo "consolidate_groups > og3.txt"
echo "*****************************************"

$APPS_DIR/consolidate_groups.pl \
    -i $WORK_DIR/og2_triangl_tupl_merged.txt \
    -d $STRING_DATABASE \
    -bc -bk -bt \
    > $WORK_DIR/og3.txt

echo "*****************************************"
echo "og3 > find_and_extend_triangles > og3_triangl.txt"
echo "*****************************************"

cat $WORK_DIR/og3.txt | $APPS_DIR/find_and_extend_triangles.pl \
    -r $WORK_DIR/reciprocal_best_hits_expanded.txt \
    -sdf $DATA_DIR/species.taxonomic_domain.txt \
    -e on \
    -d $STRING_DATABASE \
    -h 280 \
    -l 0 \
    -s 40 \
    -dngc \
    -dngk \
    > $WORK_DIR/og3_triangl.txt

echo "*****************************************"
echo "og3_triangle > find_and_extend_tuples > og3_triangl_tupl.txt"
echo "*****************************************"

cat $WORK_DIR/og3_triangl.txt | $APPS_DIR/find_and_extend_tuples.pl \
    -r $WORK_DIR/reciprocal_best_hits_expanded.txt \
    -sdf $DATA_DIR/species.taxonomic_domain.txt \
    -e on \
    -d $STRING_DATABASE \
    -h 260 \
    -l 60 \
    -s 40 \
    -mpc 5 \
    -dngc \
    -dngk \
    > $WORK_DIR/og3_triangl_tupl.txt

echo "*****************************************"
echo "summarize_orthgroups > og3_triangl_tupl_alignments.txt"
echo "*****************************************"
# commented out in latest version
$APPS_DIR/summarize_orthgroups.pl \
    -i $WORK_DIR/og3_triangl_tupl.txt \
    -d $STRING_DATABASE \
    -s 200 \
    > $WORK_DIR/og3_triangl_tupl_alignments.txt
# /home/stitch1/STITCH/eggNOG-pipeline/summarize_orthgroups.pl -i /home/stitch1/STITCH/eggNOG-pipeline/all_new/og3_triangl_tupl.txt -d string_9_1 -s 200 > $WORK_DIR/og3_triangl_tupl_alignments.txt

## from here on it is about fusion ...

echo "*****************************************"
echo "assign_orthgroup_protein_positions.pl"
echo "*****************************************"
# commented out in latest version
$APPS_DIR/assign_orthgroup_protein_positions.pl \
    -d $STRING_DATABASE \
    -o $WORK_DIR/og3_triangl_tupl.txt \
    > $WORK_DIR/og3_triangl_tupl_orthgroups.txt

# $APPS_DIR/assign_orthgroup_protein_positions.pl -d $STRING_DATABASE -o $WORK_DIR/og3_triangl_tupl.txt > $WORK_DIR/og3_triangl_tupl_orthgroups.txt



echo "*****************************************"
echo "get_best_hits > best_hits_outside_orthgroup.txt"
echo "*****************************************"

$APPS_DIR/get_best_hits.pl \
    -o $WORK_DIR/og3_triangl_tupl.txt \
    -d $STRING_DATABASE \
    > $WORK_DIR/best_hits_outside_orthgroup.txt

echo "*****************************************"
echo "find_reciprocal_best_hits > recipr_best_hits_outside_orthgroups.txt"
echo "*****************************************"

#was: -c 900000 \
$APPS_DIR/find_reciprocal_best_hits.pl \
    -i $WORK_DIR/best_hits_outside_orthgroup.txt \
    -c 100000000000 \
    > $WORK_DIR/recipr_best_hits_outside_orthgroups.txt

echo "*****************************************"
echo "find_fusion_candidates > node_fusion_candidates.txt"
echo "*****************************************"

$APPS_DIR/find_fusion_candidates.pl \
    -o $WORK_DIR/og3_triangl_tupl.txt \
    -r $WORK_DIR/recipr_best_hits_outside_orthgroups.txt \
    -c 80 \
    -t1 \
    -ltpc 10 \
    > $WORK_DIR/node_fusion_candidates.txt

# from here are empty

echo "*****************************************"
echo "execute_fusion_proteins > og3_trangl_tupl_fused.txt"
echo "*****************************************"

$APPS_DIR/execute_fusion_proteins.pl \
    -o $WORK_DIR/og3_triangl_tupl_orthgroups.txt \
    -c $WORK_DIR/node_fusion_candidates.txt \
    -d $STRING_DATABASE \
    -nmc 35 \
    -nmf 0.5 \
    -mgs 10 \
    > $WORK_DIR/og3_triangl_tupl_fused.txt

echo "*****************************************"
echo "consolidate_groups > og4.txt"
echo "*****************************************"

$APPS_DIR/consolidate_groups.pl \
    -i $WORK_DIR/og3_triangl_tupl_fused.txt \
    -d $STRING_DATABASE \
    -bc -bk -bt \
    > $WORK_DIR/og4.txt

echo "*****************************************"
echo "assign_orphans_via_Cognitor > og5.txt"
echo "*****************************************"

$APPS_DIR/assign_orphans_via_Cognitor.pl \
    -o $WORK_DIR/og4.txt \
    -d $STRING_DATABASE \
    -n 2 \
    -s 70 \
    -c \
    > $WORK_DIR/og5.txt

# dependencies 
# v7_orthgroups.txt
# v7_checksums.txt
# v8_checksums.txt


#$APPS_DIR/NOGs_update.pl \
#    -ioo $WORK_DIR/v7_orthgroups.txt \
#    -ion $WORK_DIR/og5.txt \
#    -nco $WORK_DIR/v7_checksums.txt \
#    -ncn $WORK_DIR/v8_checksums.txt \
#    > $WORK_DIR/og6.txt
#
# dependencies 
# KOG_COG_mappings.manually_edited.txt 

#$APPS_DIR/implement_external_KOG_COG_mappings.pl \
#    -if $WORK_DIR/og6.txt \
#    -mf $WORK_DIR/KOG_COG_mappings.manually_edited.txt \
#    -sdf $WORK_DIR/species.taxonomic_domain.txt \
#    > $WORK_DIR/og7.txt
#
#cat $WORK_DIR/og7.txt | $APPS_DIR/orthgroup_counter.pl \
#    > $WORK_DIR/final_count.txt
#
#$APPS_DIR/summarize_orthgroups.pl \
#    -i $WORK_DIR/og6.txt \
#    -d $STRING_DATABASE \
#    -s 200 \
#    > $WORK_DIR/final_alignments.txt
#
#$APPS_DIR/assign_orthgroup_protein_positions.pl \
#    -d $STRING_DATABASE \
#    -o $WORK_DIR/og7.txt \
#    > $WORK_DIR/new_orthgroups_inclCOGKOGmapping.txt
#
#$APPS_DIR/assign_orthgroup_protein_positions.pl \
#    -d $STRING_DATABASE \
#    -o $WORK_DIR/og6.txt \
#    > $WORK_DIR/new_orthgroups_normal.txt
#
#$APPS_DIR/make_average_species_sequence_distance_matrix.pl \
#    -d swdata_v6 \
#    -o $WORK_DIR/og7.txt \
#    -l 50 \
#    -s 10 \
#    > $WORK_DIR/average_species_similarity.txt
#
#
#$APPS_DIR/make_COG_per_COG_evolution_speed_correction.pl \
#    -d swdata_v6 \
#    -o $WORK_DIR/og7.txt \
#    -sd $WORK_DIR/average_species_similarity.txt \
#    -l 50 \
#    -s 10 \
#    -i 500 \
#    > $WORK_DIR/COG_per_COG_correction.fits.txt
