

# generate textmining groups files

for i in *.members.tsv; do cat $i |
sed 1d|perl -nale '
%h = (ssRNA => 439488, Retrotranscribing => 35268, dsDNA => 35237,Viruses => 10239, Herpesvirales => 548681, ssDNA => 29258, ssRNA_positive => 35278, Retroviridae => 11632, Ligamenvirales => 1511857, 'Caudovirales' => 28883, Mononegavirales => 11157,Tymovirales => 675063, Nidovirales => 76804,Picornavirales => 464095,dsRNA => 35325, ssRNA_negative => 35301);
@cols = split /\t/;
@prot = split /,/, $cols[4];
foreach $p (@prot) {
    @s = split /\./, $p;
    $taxid = $h{$cols[0]};
    print $cols[1] . "\t" . $taxid . "\t" . $s[0] . "\t" . $s[1];
}' > $i.list; done


# generate entities for the groups
for i in *.list; do cat $i|tab 2 1|sort|uniq > $i.ent; done
COUNT=600000000
cat *.ent|awk -v count=$COUNT -F"\t" '{print count "\t" $0; count+=1}'|sort -k 3 > groups.entities

# join on the protein entities, and generate the groups file -- child, parent
cat *.list|sort -k 4 > groups.proteins 
# add the protein serialno to the eggnog id
join -1 4 -2 3 -o 2.1,1.1 groups.proteins ../../../textmining/tm-vir-host/entities.sort3 | sort -k 2 > groups.temp
join -1 2 -2 3 -o 1.1,2.1 groups.temp groups.entities > groups




cat data_out/string_v10.species.tsv|cut -f1|sort|uniq | xargs python climb_tax_tree.py > data_out/taxtree_climbed
for i in `cat ../eggNOG/clades_of_interest.txt|grep -v "#"`; do grep $i data_out/taxtree_climbed|awk -v var=$i '{print var "\t" $0}'; done > data_out/taxtree_climbed_annot


library(dplyr)
options(max.print=100000)
options(dplyr.width = Inf)

d <- read.table("~/Google Drive/Human-Virus-PPI/data/eggNOG/eggNOG-from-site/site/Picornavirales.members.tsv.list", header=F, sep="\t")
colnames(d) <- c("nog", "level", "taxid", "prot")
d <- tbl_df(d)

length(unique(d$prot))

# get names from the textmining data
e <- read.table("~/Google Drive/Human-Virus-PPI/data/uniprot_to_payload/data_out/tm_sy_entities.uniq", header=F, sep="\t")
colnames(e) <- c("serial", "taxid", "prot")
e <- tbl_df(e)


n <- read.table("~/Google Drive/Human-Virus-PPI/data/uniprot_to_payload/data_out/tm_sy_names.uniq", header=F, sep="\t", quote="")
colnames(n) <- c("serial", "name")
n <- tbl_df(n)

d <- d %>% arrange(desc(prot))

j <- inner_join(d, e, by=c("prot" = "prot")) %>% 
    inner_join(n, by=c("serial" = "serial")) %>% 
    filter(!grepl('protein', name, ignore.case=T) & prot != name) %>% 
    group_by(prot) %>% 
    summarize(nog=first(nog), taxid=first(taxid.x), names=paste(name, collapse=" "))

k <- j %>% group_by(nog) %>% summarize(prots=paste(prot, collapse=" "), taxids=paste(taxid, collapse=" "), names=paste(names, collapse="; "))

head(k, n=1) %>% print(n=1)

k %>% filter(nog == "ENOG411EPWY")
