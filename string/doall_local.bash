#! /bin/bash


# 

## Allright, this script here is the main routine for making a STRING database update 
## (at least with respect to getting the data into the database tables).
##
## Use it by uncommenting a couple of commands at a time, running them, and then
## going to the next commands. Do not run all commands at the same time,
## as some commands require actions external to this script which in turn
## depend on earlier commands from within this script. Am I making any sense at all ... ?    
##  
## Also, please make sure that you do not overwrite older files ... get them
## out of the way before running this script. You might be interested in older
## files for comparison/reference. Essentially, make backup-copies of the 
## 'primary' and 'derived' directories, and then empty them before you start.


#  # set up postgresql on a new machine
#  sudo yum install postgresql
#  sudo yum install postgresql-server
#  sudo postgresql-setup initdb
#  sudo chkconfig postgresql on
#  sudo service postgresql start
#  # create database
#  sudo -u postgres psql
#  create database string_11;
#  sudo -u postgres psql -d string_11
#  create user readonly with password 'aCg00aisWP7GVWNp' nologin;
#  create user stringrw with password 'l4i9YYoV5b7qthJ9';
#  create user mering with password 'l4i9YYoV5b7qthJ9';
#  create user helen with password 'l4i9YYoV5b7qthJ9';
#
#  set up /var/lib/pgsql/data/pg_ident.conf
#  set up /var/lib/pgsql/data/pg_hba.conf
# 
#  # role public?
#  sudo -u postgres psql -d string_11 < string/maintenance/sql/queries/create_item_schema_tables.sql
#  sudo -u postgres psql -d string_11  < string/maintenance/sql/queries/create_evidence_schema_tables.sql
#  sudo -u postgres psql -d string_11  < string/maintenance/sql/queries/create_homology_schema_tables.sql
#  sudo -u postgres psql -d string_11  < string/maintenance/sql/queries/create_network_schema_tables.sql
#  sudo -u postgres psql -d string_11  < string/maintenance/sql/queries/create_externalpayload_schema_tables.sql
#  sudo -u postgres psql -d string_11  < string/maintenance/sql/queries/create_loginuser_schema_tables.sql
#  
#  # do I need to do this?
#  #grant select on all tables in schema items, network, evidence, homology to stringro;
#  #grant select on all tables in schema items to stringrw;
#  
#
#  string_11 is viruses first
#  string_10_5 imports the string database dumps, then adds viruses (and is the correct version)
#
#  
# to change location of data directory
# edit /etc/systemd/system/multi-user.target.wants/postgresql.service
# Environment=PGDATA=/data/pgsql/data
#
# sudo vi /var/lib/pgsql/data/postgresql.conf
# data_directory = /data/pgsql/data
#
# 
# then  (was /var/lib/pgsql/data)
# sudo vi /data/pgsql/data/pg_ident.conf
# sudo vi /data/pgsql/data/pg_hba.conf
#  
#  # end of database setup

grant all privileges on database string_11 to helen;

alter database string_11 owner to helen;
alter schema items owner to helen;
alter schema network owner to helen;
alter schema evidence owner to helen;
alter schema homology owner to helen;
alter schema external_payload owner to helen;
alter schema loginuser owner to helen;
alter schema public owner to helen;
psql -tc "select 'alter table ' || schemaname || '.' || tablename || ' owner to helen;' from pg_tables where schemaname not in ('pg_catalog', 'information_schema');" string_11 | psql -a string_11

# these don't work!! command appears to be successful, but no permissions show in \dp
grant all privileges on all tables in schema items to helen;
grant all privileges on all tables in schema network to helen;
grant all privileges on all tables in schema evidence to helen;
grant all privileges on all tables in schema homology to helen;
grant all privileges on all tables in schema external_payload to helen;
grant all privileges on all tables in schema loginuser to helen;
grant all privileges on all tables in schema public to helen;






# start of the script

export DATABASE_NAME=real_10_5
export DATABASE_HOST=localhost
export DATABASE_PORT=5432
export DATABASE_USER=helen

# TODO what is the difference between loginuser and frontend user? -- frontend user can write to tables?? and what does the mering user do?

#STRING_LOGINUSER_DATABASE_NAME=string_11
#STRING_LOGINUSER_DATABASE_USER=stringrw
#STRING_LOGINUSER_DATABASE_HOST=localhost
#STRING_LOGINUSER_DATABASE_PORT=5432

## the mechanism below allows to retrieve most data from database, but homology data from a separate database.
## this is copied from the webserver-code, note that "NEWSTRING_db.pm" and "NEWSTRING_globals.pm"
## are identical here and in the server code (CvM actually uses symlinks).

#export STRING_FRONTEND_DATABASE_VERSION=10_5
export STRING_FRONTEND_DATABASE_VERSION=11
#export STRING_FRONTEND_DATABASE_NAME=string_10_5
export STRING_FRONTEND_DATABASE_NAME=string_11
export STRING_FRONTEND_DATABASE_HOST=localhost
export STRING_FRONTEND_DATABASE_PORT=5432
export STRING_FRONTEND_DATABASE_USER=helen

export STRING_FRONTEND_HOMOLOGY_DATABASE_VERSION=10_5
export STRING_FRONTEND_HOMOLOGY_DATABASE_NAME=string_10_5
export STRING_FRONTEND_HOMOLOGY_DATABASE_HOST=localhost
export STRING_FRONTEND_HOMOLOGY_DATABASE_USER=helen
export STRING_FRONTEND_HOMOLOGY_DATABASE_PORT=5432

#export WORKING_DIR=/g/bork/mering/string/data/primary_v10/
export WORKING_DIR=/home/helen/derived_11
cd $WORKING_DIR

# on helen's laptop
#export FREEZE_DIR=/Users/bitmask/projects/phd/string/freeze
#export DERIVED_DIR=/Users/bitmask/projects/phd/string/derived_proteomes
#export OUTPUT_DIR=/Users/bitmask/projects/phd/string/derived_proteomes/output

# on blue
#export FREEZE_DIR=/home/helen/freeze_10_5
#export DERIVED_DIR=/home/helen/derived_10_5
#export OUTPUT_DIR=/home/helen/derived_10_5
export FREEZE_DIR=/home/helen/freeze_11
export DERIVED_DIR=/home/helen/derived_11
export OUTPUT_DIR=/home/helen/derived_11

# set up for viruses eggnog version 2 -- run on computerome


## the first action is to check and consolidate the names of the organisms in question. 
## this is first done by a script (see next command), but the result of that script
## needs to be checked for consistency, and edited manually if needed (some species may have 
## to be removed, because they are incomplete, or redundant, and some names may 
## have to be shortened to fit the displays).


# helen
# generate fastafiles without the species taxid with 
#for i in *.fa; do cat $i|perl -nale 'if ($_ =~ ">") { @ar = split /\./;  print ">" . $ar[1] } else {print $_}'> ../fastafiles_nospecies/$i; done
# copy data_out into the derived directory
# sshmount blue
# cd ~/mnt/blue/string_build/derived
# cp ~/Google\ Drive/Human-Virus-PPI/data/uniprot_to_payload/data_out/string_v11.species.tsv ./string_v11.species.preliminary.tsv
# cp -r ~/Google\ Drive/Human-Virus-PPI/data/uniprot_to_payload/data_out/fastafiles .
# cp ~/projects/phd/string/derived_proteomes/all.sorted.tsv ~/mnt/blue/string_build/working

# OK 11
../maintenance/make_species_name_file_and_check_protein_counts.pl > species.names.raw.txt
mv species.names.raw.txt species.names.txt

## edit and verify the file just produced (the final file should be called 'species.names.txt').

## next, let's create the species.lookup file, the protein-shorthands file, the protein-sizes file, and 
## the species_and_clades file. The latter groups species that are very close (for the orthology scripts), 
## and may have to be adjusted manually after being produced.

# OK 11
cat species.names.txt | ../maintenance/make_species_lookup_file.pl > species.lookup


# OK 11
cat species.lookup | ../maintenance/make_protein_shorthands.pl > protein.shorthands.txt
# OK 11
../maintenance/get_protein_sizes.pl > protein_sizes.txt
# OK 11
../maintenance/make_species_and_clades_file.pl > species_and_clades.txt
# edit get_crc64<tab> to not prepend species to name
../maintenance/get_crc64_checksum_entries.pl > proteins.checksums.txt
# edit assign_tax<tab> to add viruses as domain
cat species.names.txt | ../maintenance/assign_taxonomic_domain_to_species.pl > species.taxonomic_domain.txt

# TODO don't have gene_position.tsv required for this
## ../../maintenance/computations/neighborhood_blue_mode/runall_create_runs.pl > runs.genes.txt

# TODO don't have /kegg/kegg_benchmarking.CONN_maps_in.tsv required for this
## ../maintenance/compute_prior.py > priors.txt

# OK
## check the homology data (currently produced by SIMAP's Thomas Rattai), and create the corresponding database tables.
../maintenance/check_simap_completeness.pl >& simap.completeness.report





# XXX block after this replaces this block because my data is small

##
## this next line is run externally, using qsub: 

## find /g/bork/mering/STRING_freeze_v10/simap/fileshare.csb.univie.ac.at/f8f962bbc097af/ -name '2013*.gz' | perl -ne '($gh) = /(2013\d+.*gz)/; print "$gh\n";' | /g/bork/mering/string/maintenance/setup_simap_pruning_jobs.pl > submitall.sh

## next, sort those separately ("sort -k 1g -k 2g"), and then merge-sort them while piping them through 
## 'uniq' (sort -k 1g -k 2g --merge sorted_files.*.tsv | uniq all.sorted.tsv)
## NOTE: sorting cannot be done on alpha in HD, will crash with out of memory. Do the sorting on atlas in Zurich, where it works just fine, using 'setup_simap_sorting_jobs.pl' to run three sort jobs in parallel. The 'sort' command on atlas is multithreaded, even. 
## NOTE 2: on atlas, everything can be done in one command: sort --temporary-directory=/home/mering/string/data/primary_v10/ -k 1g -k 2g parsed.all.txt | uniq | gzip -c > sorted.all.gz

# ../maintenance/create_best_hit_matrix.pl | gzip -3 > $OUTPUT_DIR/fill.best_hit_per_species.table.sql.gz




# XXX start here
# my data is small
../maintenance/generate_homology_table_data.pl ../SIMAP/homologs.gz > SIMAP_homologs.parsed_unsorted.tsv
sort -k 1g -k 2g SIMAP_homologs.parsed_unsorted.tsv | uniq > all.sorted.tsv
#open this in vi and switch \t for ,


# the next one depends on data already being in the database, see below
../maintenance/create_best_hit_matrix.pl | gzip -3 > $OUTPUT_DIR/fill.best_hit_per_species.table.sql.gz



## the next block creates a set of files to fill most of the central 'parts-list' tables in STRING.


# OK 11
mkdir output
mkdir aliases_and_descriptions
# pretend to have some files
touch COGS.txt
touch NOGS.txt
touch final_orthgroups.txt
touch aliases_and_descriptions/alias_best.tsv
touch aliases_and_descriptions/alias_all.tsv
touch aliases_and_descriptions/text_all.tsv
touch aliases_and_descriptions/text_best.tsv
touch gene_position.tsv
touch runs.genes.txt
../maintenance/parse_primary_data.pl > $OUTPUT_DIR/populate_items_schema_1.sql

# TODO script convert_to_utf8.pl not in repository
## ../maintenance/convert_to_utf8.pl $OUTPUT_DIR/populate_items_schema_1.sql > $OUTPUT_DIR/populate_items_schema_1.utf8.sql

# edit this script to refer to fastafiles not fastafiles_db
../maintenance/parse_protein_sequences_from_fasta_files.pl > $OUTPUT_DIR/populate_items_schema_2.sql

# TODO still needed to edit names by hand because they are too long
sudo -u postgres psql -d string_11 < ../output/populate_items_schema_1.sql
sudo -u postgres psql -d string_11 < ../output/populate_items_schema_2.sql

# only the homology needs to be in the database to run eggNOG

# then write this stuff into the database with
#cat string_v10.species.tsv|awk -F"\t" {'print $1 "," $2 "," $3 ",virus," $6'} > load_db_items.species
# edit this file by hand to remove fields that are too big
#/home/stitch1/STITCH/20150123-install/bin/psql -h 127.0.0.1 -p 5432 -U helen -d string_9_1
#copy items.species from '/home/people/helen/load_db_items.species' (FORMAT CSV, DELIMITER(',')); 
# cut the species stuff out of populate_items_schema_1.sql by hand 
#/home/stitch1/STITCH/20150123-install/bin/psql -h 127.0.0.1 -p 5432 -U helen -d string_9_1 < ~/populate_items_schema_1.sql
#/home/stitch1/STITCH/20150123-install/bin/psql -h 127.0.0.1 -p 5432 -U helen -d string_9_1 < ~/populate_items_schema_2.sql

# these tables seem to be optional for the rest of the pipeline to work
#../maintenance/parse_protein_sequences_from_fasta_files.pl > $OUTPUT_DIR/populate_items_schema_2.sql
#../maintenance/convert_to_utf8.pl $OUTPUT_DIR/populate_items_schema_1.sql > $OUTPUT_DIR/populate_items_schema_1.utf8.sql
#../maintenance/parse_primary_data.pl > $OUTPUT_DIR/populate_items_schema_1.sql

# write the homology data into the database

psql -d string_10_5 < fill.best_hit_per_species.table.sql
psql -d string_10_5 < <( cat <(echo "Copy homology.similarity_data from stdin; " ) all.sorted.tsv)


# then, can run the eggnog pipeline on stitch!

# XXX XXX XXX done to here

##../../maintenance/ncbi_clades.py
## cat $DERIVED_DIR/ncbi_taxonomy/ncbi_taxonomy.sql | gzip -3 > $OUTPUT_DIR/populate_items_schema_4.sql.gz

## play in the above two files, then create the fulltext-search field and the index 
## (this is done externally by playing in 'create_annotation_fulltext_index.sql')
## lastly, create all the other indices and run 'vacuum analyze'.



# TODO requires chromosome positions -- and KEGG pathway info and ... 
## next, compute the neighborhood scores ... last time, machine 'dag' was used for this: needs xmgrace, GD, ...
## ../../maintenance/computations/neighborhood_blue_mode/runall_compute_scores.pl

## external: execute the 'fastcoc' script to precompute the mutual information.
## NOTE: There are two parameters which need checking in the fastcoc script: Max species and when to ignore after collapse!!
## NOTE: There are also parameters in the first script below ... that is the names of eukaryotic species.
## note also - comments above are slightly out of date. There have been various changes in the fastcoc setup, see there ...

###########################################################################################
## Protein-mode scores and evidence transfer
##
##
###########################################################################################

# ran
## ../maintenance/create_self_hits_dump_file.pl > $OUTPUT_DIR/self_hits.blast_data.tsv

# ran
## ../maintenance/make_compact_score_files_for_transfer.pl 

# OK
## NOTE: run the next step on the 'atlas' machine in Zurich, to handle the substantial memory load.
 ../../maintenance/create_transfer_input_files.pl                              

## ../../maintenance/check_symmetry_of_transfer_input_files.pl > evidence_transfer/input_files_symmetry_report.txt
## ../../maintenance/check_counts_per_evidence_type.pl evidence_transfer/input/*.gz > evidence_transfer/input.counts.checking.txt

## then, do the actual transfer - this is external (run on a cluster using qsub; currently the cluster 'alpha' at EMBL is used)

## external command is: mering@alpha:transfer/transfer_evidence.py submit

## check that transfer was successful
## ../../maintenance/transfer/transfer_evidence.py check

## ok, now combine direct scores and transferred scores, and play that in.

../maintenance/precompute_protein_links.pl > $OUTPUT_DIR/populate_network_schema_6.sql

## then, do the transfer reporting. for this, first create compact files describing what is worth reporting:

## ../../maintenance/create_compact_combinedscore_dumps.pl

## then, for the transfer reporting (and the parsing of the results), this is external again.

## mering@alpha[lately mkuhn]:
## /g/bork/mering/string/data/primary/evidence_transfer/transfer_evidence.py submit report

## now, remember that the next two python scripts will require a perl-program to run in the background and to make 
## available some functions to python, via remote-procedure-calls on some port above 8000.
## If the cluster is to be used, this perl script has to be running on each node, started manually by bypassing
##  the submission system ... 
## mering@alpha: 
## /g/bork/mering/string/maintenance/launch_compute_server.bash

## @alpha:
## parse_evidence_transfer_reports.py submit

# cat evidence_transfer/parsed_reports/*.sql.gz > ../derived/populate_evidence_schema_8.sql.gz

## @sigma:
## import_actions.py submit

## now, concatenate the information on actions:

## cat evidence_transfer/sql_actions/*.sql | gzip -c > ../derived/populate_network_schema_8.sql.gz
## cat evidence_transfer/sql_sources/*.sql | gzip -c > ../derived/populate_evidence_schema_9.sql.gz

##################################################################################
## COG-mode scores
##
##
##################################################################################
 
## first, we collect non genomic-context scores from the protein mode.

## ../../maintenance/precompute_COG_evidence_channels.pl > ../derived/COG_evidence_associations.txt

## next, compute an initial version of the full scores, using the genomic-context sigmoids of the last release 
## (at this point, only the raw scores of the genomic-context channels are needed).
##
## for minor releases, these steps can be skipped, since the sigmoids will not change ...

## ../../maintenance/precompute_orthgroup_links.pl preparation > ../derived/populate_network_schema_7.sql

## play in the above file, and compute new sigmoids using command below:

# ../../maintenance/compute_orthgroup_performance_curves.pl -nbn 1000 -nbf 50 -nbc 1000 -ipegf -ipegc \
# > ../derived/COG_mode_genomic_context_sigmoids.txt

## (note, to get the actual sigmoid curve parameters, run 'xmgrace' on machine 'string32' ...)

## and re-run the orthgroup precomputation with the correct sigmoids (these have to be updated in the code!)
## we now also use the final combined-score cutoff (for example 120).

## ../../maintenance/precompute_orthgroup_links.pl final > ../derived/populate_network_schema_7.sql

## ok, now read in the accessory stuff needed to trace and present evidence in the experiments and databases channels.

##################################################################################
## fill evidence tables 
##
##
##################################################################################

# EVIDENCE FOR VIRUSES
# we are adding on to the database, so don't need to import everything that the main scripts import

# from ~/derived_10_5   OK
#../maintenance/create_experiment_database_annotation_tables.pl > populate_evidence_schema_1.sql


## ../../maintenance/create_experiment_database_annotation_tables.pl > ../derived/populate_evidence_schema_1.sql

## next, read in accessory data regarding the textmining.
## please note: this will also import the abstracts referred to by the experimental data created in the previous step,
## --> so, make sure that has been played in and is up-to-date.

## this next steps needs to be run on machine 'atlas' in Zurich, consumes too much memory for delta.
#  OK
../maintenance/parse_textmining_data.pl > ../derived/populate_evidence_schema_2.sql

# OK
## ../../maintenance/parse_nlp_sentences.pl > ../derived/populate_evidence_schema_3.sql

## fusion reporting: the actual reporting is done by Milan ... 

## two alternatives for the sorting (I used the second option on 'delta.embl' last time):

# sort -T fusion_reports +0n +1n +2n +3n fusion_reports/fusion.reports.txt > fusion_reports/fusion.reports.sorted.txt
# sort -T fusion_reports -k 1,1n -k 2,2n -k 3,3n -k 4,4n fusion_reports/fusion.reports.txt > fusion_reports/fusion.reports.sorted.txt

# cat fusion_reports/fusion.reports.sorted.txt | ../../maintenance/parse_fusion_reports.pl \
# > ../derived/populate_evidence_schema_4.sql

## ../../maintenance/connect_orthgroups_to_abstracts_etal.pl > ../derived/populate_evidence_schema_5.sql

## ../../maintenance/make_blast_database_files.pl 

###################################################################################
## ivica (smart) link-out and link-in:
##
## note: "link-in" means from a URL in SMART that links to STRING, and vice versa
###################################################################################

## next command needs to be run on the production machine - or else the network layouts may not be fully compatible ...

## /srv/servers/string/cgi_bin/generate_link_in_data_for_smartdb.pl > smart_link_in_data_v8_1.txt
## NOTE: the above is now done by Milan - so the call is obsolete !!

## from here on, continue as before.

## cat ../derived/smart_link_ins/string_linkout_info_v6_2.txt__ | ../../maintenance/reduce_smart_link_in_data_to_model_organisms.pl > ../derived/smart_link_ins/string_linkout_info_v6_2.txt

## cat smartdb_link_outs/smart_links_for_string | ../../maintenance/generate_smartdb_link_out_data.pl \
## > ../derived/populate_items_schema_3.sql

##################################################################################
## sitemaps 
##
##################################################################################

# ../../maintenance/generate_sitemaps.py

##################################################################################
## download files: this is external to this script. 
##
## go to the sub-directory 'download_files' and check the script 
## 'make_download_files.bash'
##################################################################################

##################################################################################
## updates ... these should eventually move more upstream, into the COPY statements
##################################################################################

## ../../maintenance/linkout/update_proteins_orthgroups.pl \
## > ../derived/update_items_schema_8.sql

##################################################################################
## things not to forget !!!!
##
##
##################################################################################

## A)   Alex has a script that will try to make as many useful back-links 
##      for our evidence sets as possible.
##
## B)   Alex also has a script to provide correct link-outs to the proteins
##      themselves.
