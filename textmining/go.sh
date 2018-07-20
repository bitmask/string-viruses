#!/bin/zsh


################################################################################
# Textmining
################################################################################


# whitelist for stopwords built in project/tm-string

#cat virus_names_all.tsv | awk '{print $0 "\t" 5}' > virus_names_all_score.tsv
#cat virus_names_all_score.tsv string_names.tsv > host_and_virus_names.tsv
#cat viruses.species.tsv|awk '{print $1 "\t" $1}' > virus_self_pairs.tsv
#cat virus_type_pairs.tsv virus_self_pairs.tsv > virus_pairs.tsv

zcat ../string_10_5_corpus/all.tsv.gz | ../tagger/tagcorpus --types=types-autodetect --entities=string_entities.tsv --names=host_and_virus_names.tsv --stopwords=all_global_pluswhitelist --type-pairs=virus_pairs.tsv --threads=16 --groups=string_groups.tsv --autodetect --out-pairs=Medline.outpairs > Medline.output

#./string_cooccur_pairs.pl --entities_file=string_entities.tsv --pairs_file=Medline.outpairs Medline

# get kegg_benchmark.tsv
# get viruses.species.tsv
#zcat Medline_cooccur_pairs.tsv.gz|awk -F'\t' '$1==9606 {print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5}'|awk -F'\t' '$3 == 9606 {print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5}' > human
#./string_cooccur_interact_viruses.pl Medline

#zcat string_virus_cooccur_interact.tsv.gz|perl -nale '@c = split "\t"; @t = split "-", @c[1]; print $t[0] . "." . $c[2] . "\t" . $t[1] . "." . $c[3] . "\t" . @c[4]*1000 ;' > string_virus_cooccur_interact_tmp
#perl ../experiments/extract-virus.pl --scores=string_virus_cooccur_interact.tsv.gz --viruses=viruses.species.tsv --output=textmining_interact.tsv --mode=textmining


#./string_complex.pl Medline
#./string_info.pl Medline

################################################################################
#NLP
################################################################################

#mkdir nlp
#copy nlp directories with templates from textmining repo
#
#cd nlp
#./nlp_split.pl Medline 28
#cd ..

# submit the nlp jobs
#for i in `seq 0 27`; do qsub -W group_list=jensenlab -A jensenlab -l nodes=1:ppn=1,mem=8gb,walltime=100000 -d /home/people/hecook/projects/tm-string-virus-final -v idx=$i nlp.sh; sleep 1; done

# when those are done, run:
#zcat Medline*_interact.tsv.gz | gzip -c > Medline_nlp_interact.tsv.gz 

################################################################################
# MERGE interactions
################################################################################

perl merge_scores.pl --output=lars --nlp=nlp/Medline_nlp_interact.tsv --textmining=textmining_interact.tsv --experiments=../experiments/experiments_interact.tsv --shorthands=protein.shorthands.all > virus.protein.links.detailed.10.5.txt
perl merge_scores.pl --output=db --nlp=nlp/Medline_nlp_interact.tsv --textmining=textmining_interact.tsv --experiments=../experiments/experiments_interact.tsv --shorthands=protein.shorthands.all > db_import_virus_host_network

