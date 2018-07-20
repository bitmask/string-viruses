



library(ggplot2)
library(dplyr)
scores <- read.table("~/Google Drive/Human-Virus-PPI/data/virus_interactions/species_scores33")
colnames(scores) <- c("tax1", "tax2", "channel", "score")
scores <- tbl_df(scores) %>% filter(score != 0)
levels(scores$channel) <- c("experiments", "transfer (ex)", "text mining", "transfer (tm)")
cols <- c("#AE6962", "#DBAE9B", "#6A8E8C", "#99B4AF")

viruses <- read.table("~/Google Drive/Human-Virus-PPI/data/uniprot_to_payload/data_in/string_v11.species.manual.tsv", sep="\t", fill=T, quote="")
colnames(viruses) <- c("tax1", "longname", "shortname", "abbrev", "uniprot", "core")

scores <- left_join(scores, viruses)
scores <- left_join(scores, viruses, c("tax2"="tax1"))
k <- scores

t <- k %>% mutate(shortname = ifelse(is.na(shortname.x), as.character(shortname.y), as.character(shortname.x)))
t <- t %>% mutate(abbrev = ifelse(is.na(abbrev.x), as.character(abbrev.y), as.character(abbrev.x)))
t <- t %>% mutate(taxid = ifelse(tax1 %in% viruses$tax1, as.character(tax1), as.character(tax2)))
t$shortname.x <- NULL
t$longname.x <- NULL
t$core.x <- NULL
t$uniprot.x <- NULL
t$uniprot.y <- NULL
t$core.y <- NULL
t$longname.y <- NULL
t$shortname.y <- NULL
t$abbrev.x <- NULL
t$abbrev.y <- NULL
t$shortname <- factor(t$shortname)
t[t$abbrev == "Murid cytomegalovirus",]$abbrev <- "MCMV"
t[t$abbrev == "Human cytomegalovirus",]$abbrev <- "HCMV"
t$abbrev <- factor(t$abbrev)
t$taxid <- factor(t$taxid)
summary(t)
t <- t %>% mutate(intra = ifelse(as.character(tax1) == as.character(tax2), "intra", "inter"))
t$intra = factor(t$intra)
g <- t %>% group_by(taxid, shortname, abbrev) %>% summarize(hosts=length(unique(c(unique(tax2), unique(tax1)))), total_n=n())

l <- left_join(t, g)

ffi <- l %>% group_by(taxid, shortname, abbrev, intra, channel, total_n ) %>% summarize(n=n())
ff <- l %>% group_by(taxid, shortname, abbrev, channel, total_n ) %>% summarize(n=n())


hosts <- read.table("~/Google Drive/Human-Virus-PPI/data/uniprot_to_payload/data_out/hosts", header=F, sep="\t")
colnames(hosts) <- c("taxid", "host")
hosts$taxid <- factor(hosts$taxid)

ho <- left_join(t, hosts) 

ho %>% group_by(channel, intra) %>% summarize(n=n())



# bar chart of interaction type for the top viruses
ggplot(ff %>% filter(total_n>2500), aes(x=reorder(shortname, total_n), y=n/2, fill=channel)) +  geom_bar(alpha=1, stat="identity") + theme_bw() + theme(axis.text.x=element_text(angle=0)) + scale_fill_manual(values = cols) + scale_color_manual(values=cols) + theme(axis.title.x=element_blank(), axis.title.y=element_blank()) + coord_flip() + theme(legend.position = c(0.9, 0.2))

# number of interactions for each channel
ggplot(ffi , aes(x=channel, y=n, fill=channel)) +  geom_histogram(stat="identity") + theme_bw() + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + scale_fill_manual(values = cols) + scale_color_manual(values=cols) + theme(axis.title.x=element_blank(), axis.title.y=element_blank())

# box plot of score by channel
ggplot(l %>% filter(total_n>0), aes(x=channel, y=score, fill=channel)) +  geom_boxplot() + theme_bw() + theme() + scale_fill_manual(values = cols) + scale_color_manual(values=cols) + theme(axis.title.x=element_blank(), axis.title.y=element_blank()) + theme(legend.position="none")

# box plot of intra/inter interaction counts by channel
ggplot(ffi , aes(x=intra, y=n, fill=channel)) +  geom_histogram(stat="identity") + theme_bw() + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + scale_fill_manual(values = cols) + scale_color_manual(values=cols) + theme(axis.title.x=element_blank(), axis.title.y=element_blank())
# the other way around
ggplot(ffi , aes(x=channel, y=n, fill=intra)) +  geom_histogram(stat="identity") + theme_bw() + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + scale_fill_manual(values = cols) + scale_color_manual(values=cols) + theme(axis.title.x=element_blank(), axis.title.y=element_blank())

# plot all interactions for top viruses
ggplot(l %>% filter(total_n>1000), aes(x=reorder(abbrev, total_n), y=score, col=channel)) +  geom_point(alpha=0.3) + theme_bw() + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + scale_fill_manual(values = cols) + scale_color_manual(values=cols) + theme(axis.title.x=element_blank(), axis.title.y=element_blank()) + coord_flip()



counts <- read.table("~/Google Drive/Human-Virus-PPI/data/virus_interactions/protein_counts", header=F)
colnames(counts) <- c("taxid", "proteins")
counts$taxid <- factor(counts$taxid)

vc <- left_join(ff, counts)
vc_simple <- vc %>% group_by(taxid,shortname,proteins) %>% summarize(total_n=max(total_n)) %>% filter(total_n>2500)
ggplot(vc_simple, aes(x=reorder(shortname, total_n/proteins), y=total_n/(2*proteins))) +  geom_bar(alpha=1, stat="identity") + theme_bw() + theme(axis.text.x=element_text(angle=0)) + scale_fill_manual(values = cols) + scale_color_manual(values=cols) + theme(axis.title.x=element_blank(), axis.title.y=element_blank()) + coord_flip() + theme(legend.position = c(0.9, 0.2))
ggsave("~/Google Drive/Human-Virus-PPI/data/virus_interactions/viruses_interactiontypes_normalized.pdf")

library(cowplot)

un <- ungroup(ff) %>% arrange(-total_n) %>% select(taxid) 
un_taxid <- unique(un$taxid)[1:20]

n <- ungroup(vc_simple) %>% arrange(-total_n/proteins) %>% select(taxid)
n_taxid <- unique(n$taxid)[1:20]


unnorm <- ggplot(ff %>% filter(taxid %in% un_taxid), aes(x=reorder(shortname, total_n), y=n/2, fill=channel)) +
    geom_bar(alpha=1, stat="identity") +
    scale_fill_manual(values = cols) +
    scale_color_manual(values=cols) +
    coord_flip() + 
    theme_bw() +
    theme(axis.text.x=element_text(angle=0), axis.title.y=element_blank(), legend.position = c(0.9, 0.3)) +
    ylab("Evidences")

norm <- ggplot(vc_simple %>% filter(taxid %in% n_taxid), aes(x=reorder(shortname, total_n/proteins), y=total_n/(2*proteins))) + 
    geom_bar(alpha=1, stat="identity") + 
    scale_fill_manual(values = cols) + 
    scale_color_manual(values=cols) + 
    coord_flip() + 
    theme_bw() + 
    theme(axis.text.x=element_text(angle=0), axis.title.y=element_blank(), legend.position = "none") + 
    ylab("Evidence/Protein")

pdf(file="combined.pdf", height=7, width=9)
plot_grid(unnorm, norm, ncol=1, align="v", labels=c('A', 'B'))
dev.off()


#hiv, siv
ggplot(scores %>% filter(tax1==11676 | tax1==11723), aes(x=score, fill=channel)) + geom_histogram(position="dodge") + theme_bw() + scale_fill_manual(values=cols, drop=FALSE) + facet_grid(abbrev ~ ., scales="free_y")
#influenza a, b
ggplot(scores %>% filter((tax1==11320 & (tax2==9606 | tax2==11320)) | (tax1==11520 & (tax2==9606 | tax2==11520))), aes(x=score, fill=channel)) + geom_histogram(position="dodge") + theme_bw() + scale_fill_manual(values=cols, drop=FALSE) + facet_grid(abbrev ~ ., scales="free_y")
ggsave("~/Google Drive/Human-Virus-PPI/data/virus_interactions/influenza_scores.pdf")
#hpvs
ggplot(scores %>% filter(tax1 %in% c(333924, 333926, 333930, 334207, 37118, 10602, 334203, 334204, 334205, 333928, 333929, 337039, 337041, 337042, 337043, 337044, 337048, 337049, 337050, 337051, 333767, 333754, 333757, 333766)), aes(x=score, fill=channel)) + geom_histogram(position="dodge") + theme_bw() + scale_fill_manual(values=cols, drop=FALSE) + facet_grid(tax1 ~ ., scales="free_y")
#herpes
ggplot(scores %>% filter(tax1 %in% c(10298, 10310, 10335, 10359, 10372, 10376, 37296, 32603, 32604)), aes(x=score, fill=channel)) + geom_histogram(position="dodge") + theme_bw() + scale_fill_manual(values=cols, drop=FALSE) + facet_grid(tax1 ~ ., scales="free_y")
#*cmv
ggplot(scores %>% filter((tax1==10359 & (tax2==9606 | tax2==10359)) | (tax1==10366 & (tax2==10090 | tax2==10366))), aes(x=score, fill=channel)) + geom_histogram(position="dodge") + theme_bw() + scale_fill_manual(values=cols, drop=FALSE) + facet_grid(abbrev ~ ., scales="free_y")
ggsave("~/Google Drive/Human-Virus-PPI/data/virus_interactions/cytomegalovirus_scores.pdf")


# generate plots for presentation
ggplot(scores, aes(x=channel, y=score, col=channel)) + geom_boxplot() + theme_bw() + scale_color_manual(values=cols, drop=FALSE) + theme(legend.position="none")
ggsave("~/Google Drive/Human-Virus-PPI/data/virus_interactions/overview.pdf")

ggplot(scores, aes(x=channel, fill=channel)) + geom_histogram(stat="count") + theme_bw() + scale_fill_manual(values=cols, drop=FALSE) + theme(legend.position="none")
ggsave("~/Google Drive/Human-Virus-PPI/data/virus_interactions/overview-count.pdf")



# plot the inter and intra virus interactions 

library(tidyr)
p <- ffi %>% group_by(taxid, intra, total_n) %>% summarize(n=sum(n)) %>% mutate(p=n/total_n) 
q <- spread(p %>% select(taxid, intra, n, total_n), "intra", "n")
q[is.na(q$inter),]$inter <- 0
q[is.na(q$intra),]$intra <- 0
ggplot(q, aes(x=inter, y=intra, col=log(total_n), label=taxid)) + geom_point() + theme_bw() + scale_x_log10() + scale_y_log10() + coord_fixed() + geom_abline(slope=1,intercept=0) + geom_text(size=2, vjust=-1)
ggsave("~/Google Drive/Human-Virus-PPI/data/virus_interactions/intra-inter.pdf")




# how many interactions have evidence from both tm and exp interactions

int <- read.table("~/Google Drive/Human-Virus-PPI/data/virus_interactions/lars_db_with_transfer_sym_all_noself")
colnames(int) <- c("string1", "string2", "u", "un", "unu", "unus", "unuse", "ex", "unused", "tm", "comb")
int <- tbl_df(int) %>% select("string1", "string2", "tm", "ex", "comb")
int <- int %>% separate("string1", c("tax1", "prot1"), "\\.", remove=T, extra="merge") %>% separate("string2", c("tax2", "prot2"), "\\.", remove=T, extra="merge")
int$tax1 <- as.numeric(int$tax1)
int$tax2 <- as.numeric(int$tax2)
int_v <- left_join(int, viruses, by=c("tax1"="tax1"))
int_v <- left_join(int_v, viruses, by=c("tax2"="tax1"))

t <- int_v %>% mutate(shortname = ifelse(is.na(shortname.x), as.character(shortname.y), as.character(shortname.x)))
t <- t %>% mutate(abbrev = ifelse(is.na(abbrev.x), as.character(abbrev.y), as.character(abbrev.x)))
t <- t %>% mutate(taxid = ifelse(tax1 %in% viruses$tax1, as.character(tax1), as.character(tax2)))
t <- t %>% mutate(host_taxid = ifelse(tax1 %in% viruses$tax1, as.character(tax2), as.character(tax1)))
t$shortname.x <- NULL
t$longname.x <- NULL
t$core.x <- NULL
t$uniprot.x <- NULL
t$uniprot.y <- NULL
t$core.y <- NULL
t$longname.y <- NULL
t$shortname.y <- NULL
t$abbrev.x <- NULL
t$abbrev.y <- NULL
t$shortname <- factor(t$shortname)
t[t$abbrev == "Murid cytomegalovirus",]$abbrev <- "MCMV"
t[t$abbrev == "Human cytomegalovirus",]$abbrev <- "HCMV"
t$abbrev <- factor(t$abbrev)
t$taxid <- factor(t$taxid)

ggplot(t, aes(x=tm, y=ex)) + geom_point()

t[t$tm > 0,]$tm <- 1
t[t$ex > 0,]$ex <- 1


int_t <- t %>% group_by(taxid) %>% summarize(tmex=sum(tm>0 & ex>0), tm=sum(tm>0 & ex==0), ex=sum(tm==0 & ex>0))

int_t %>% filter(ex==0 & tmex==0) # 154 viruses have only tm evidence
int_t %>% filter(ex==0 & tmex>0) # 77 viruses have all ex evidence confirmed by tm, and extra tm
int_t %>% filter(tm==0 & tmex==0) # 3 viruses have only ex evidence
int_t %>% filter(tm==0 & tmex>0) # 2 viruses have all tm evidence confirmed by ex, and extra ex

#154+77+3+2
#236 viruses total


int_v <- t %>% group_by(taxid, host_taxid) %>% summarize(tmex=sum(tm>0 & ex>0), tm=sum(tm>0 & ex==0), ex=sum(tm==0 & ex>0))
s <- gather(int_v, key="key", value="value", c(3:5))
ggplot(s %>% filter(host_taxid == 9606, taxid==11676), aes(x=taxid, y=host_taxid, col=key, size=value)) + geom_point(alpha=0.3) + theme_bw() + facet_grid(key ~ .)


int_v %>% filter
exconf <- int_v %>% filter(ex==0, tmex>0)
exconf %>% ungroup() %>% summarize(sum=sum(tmex))

overlap <- int_v %>% mutate(peroverlap = tmex / (tmex+tm+ex))
ggplot(overlap, aes(x=taxid, y=host_taxid, size=peroverlap)) + geom_point(alpha=0.5)


overlap %>% filter(tm+ex+tmex > 100, host_taxid==9606) %>% print(n=100)

sum <- overlap %>% group_by(host_taxid) %>% summarize(mean=mean(peroverlap))
ggplot(sum %>% filter(mean > 0.3), aes(x=host_taxid, y=mean)) + geom_point()













# old crap follows



library(ggplot2)
library(dplyr)
scores <- read.table("~/Google Drive/Human-Virus-PPI/data/virus_interactions/db_import_virus_host_network_transfertest.r2")
colnames(scores) <- c("tax1", "tax2", "channel", "score")
scores <- tbl_df(scores) %>% filter(score != 0)
levels(scores$channel) <- c("experiments", "transfer (ex)", "text mining", "transfer (tm)")
cols <- c("#AE6962", "#DBAE9B", "#6A8E8C", "#99B4AF")

viruses <- read.table("~/Google Drive/Human-Virus-PPI/data/uniprot_to_payload/data_in/string_v11.species.manual.tsv", sep="\t", fill=T, quote="")
colnames(viruses) <- c("tax1", "longname", "shortname", "abbrev", "uniprot", "core")

scores <- left_join(scores, viruses)

#hiv, siv
ggplot(scores %>% filter(tax1==11676 | tax1==11723), aes(x=score, fill=channel)) + geom_histogram(position="dodge") + theme_bw() + scale_fill_manual(values=cols, drop=FALSE) + facet_grid(abbrev ~ ., scales="free_y")
#influenza a, b
ggplot(scores %>% filter((tax1==11320 & (tax2==9606 | tax2==11320)) | (tax1==11520 & (tax2==9606 | tax2==11520))), aes(x=score, fill=channel)) + geom_histogram(position="dodge") + theme_bw() + scale_fill_manual(values=cols, drop=FALSE) + facet_grid(abbrev ~ ., scales="free_y")
ggsave("~/Google Drive/Human-Virus-PPI/data/virus_interactions/influenza_scores.pdf")
#hpvs
ggplot(scores %>% filter(tax1 %in% c(333924, 333926, 333930, 334207, 37118, 10602, 334203, 334204, 334205, 333928, 333929, 337039, 337041, 337042, 337043, 337044, 337048, 337049, 337050, 337051, 333767, 333754, 333757, 333766)), aes(x=score, fill=channel)) + geom_histogram(position="dodge") + theme_bw() + scale_fill_manual(values=cols, drop=FALSE) + facet_grid(tax1 ~ ., scales="free_y")
#herpes
ggplot(scores %>% filter(tax1 %in% c(10298, 10310, 10335, 10359, 10372, 10376, 37296, 32603, 32604)), aes(x=score, fill=channel)) + geom_histogram(position="dodge") + theme_bw() + scale_fill_manual(values=cols, drop=FALSE) + facet_grid(tax1 ~ ., scales="free_y")
#*cmv
ggplot(scores %>% filter((tax1==10359 & (tax2==9606 | tax2==10359)) | (tax1==10366 & (tax2==10090 | tax2==10366))), aes(x=score, fill=channel)) + geom_histogram(position="dodge") + theme_bw() + scale_fill_manual(values=cols, drop=FALSE) + facet_grid(abbrev ~ ., scales="free_y")
ggsave("~/Google Drive/Human-Virus-PPI/data/virus_interactions/cytomegalovirus_scores.pdf")


# generate plots for presentation
ggplot(scores, aes(x=channel, y=score, col=channel)) + geom_boxplot() + theme_bw() + scale_color_manual(values=cols, drop=FALSE) + theme(legend.position="none")
ggsave("~/Google Drive/Human-Virus-PPI/data/virus_interactions/overview.pdf")

ggplot(scores, aes(x=channel, fill=channel)) + geom_histogram(stat="count") + theme_bw() + scale_fill_manual(values=cols, drop=FALSE) + theme(legend.position="none")
ggsave("~/Google Drive/Human-Virus-PPI/data/virus_interactions/overview-count.pdf")





#:r

install.packages("dplyr")
install.packages("ggplot2")
library(dplyr)
library(ggplot2)
d <- read.table("~/Google Drive/Human-Virus-PPI/data/virus_interactions/species_scores33", sep="\t")
d <- tbl_df(d)
colnames(d) <- c("tax1", "tax2", "channel", "score")
d <- d %>% filter(tax1 >= tax2)
d$taxa <- factor(d$tax1)
d$taxb <- factor(d$tax2)
d <- d %>% filter(tax1 != 1, tax2 != 1, score > .150)
d
summary(d)

viruses <- read.table("~/Google Drive/Human-Virus-PPI/data/uniprot_to_payload/data_in/string_v11.species.manual.tsv", sep="\t", fill=T, quote="")
colnames(viruses) <- c("taxid", "longname", "shortname", "abbrev", "uniprot", "core")
viruses$longname <- NULL
viruses$uniprot <- NULL
viruses$core <- NULL
viruses$taxid <- factor(viruses$taxid)

# join virus data
j <- left_join(d, viruses, c("taxa"="taxid")) 
j$taxa <- factor(j$taxa)
summary(j)
j

# join virus data for viruses with taxid > host
k <- left_join(j, viruses, c("taxb"="taxid"))
k$taxb <- factor(k$taxb) 
summary(k)
k

#scores <- left_join(scores, viruses, c("tax2"="taxid"))
#k <- scores

t <- k %>% mutate(shortname = ifelse(is.na(shortname.x), as.character(shortname.y), as.character(shortname.x)))
t <- t %>% mutate(abbrev = ifelse(is.na(abbrev.x), as.character(abbrev.y), as.character(abbrev.x)))
t <- t %>% mutate(taxid = ifelse(tax1 %in% viruses$tax1, as.character(tax1), as.character(tax2)))
t$shortname.x <- NULL
t$longname.x <- NULL
t$core.x <- NULL
t$uniprot.x <- NULL
t$uniprot.y <- NULL
t$core.y <- NULL
t$longname.y <- NULL
t$shortname.y <- NULL
t$abbrev.x <- NULL
t$abbrev.y <- NULL
t$shortname <- factor(t$shortname)
t[t$abbrev == "Murid cytomegalovirus",]$abbrev <- "MCMV"
t[t$abbrev == "Human cytomegalovirus",]$abbrev <- "HCMV"
t$abbrev <- factor(t$abbrev)
t$taxid <- factor(t$taxid)
summary(t)
t <- t %>% mutate(intra = ifelse(as.character(taxa) == as.character(tax2), "intra", "inter"))
t$intra = factor(t$intra)
g <- t %>% group_by(taxid, shortname, abbrev) %>% summarize(hosts=length(unique(c(unique(tax2), unique(tax1)))), total_n=n())

l <- left_join(t, g)

ffi <- l %>% group_by(taxid, shortname, abbrev, intra, channel, total_n ) %>% summarize(n=n())
ff <- l %>% group_by(taxid, shortname, abbrev, channel, total_n ) %>% summarize(n=n())


g <- t %>% group_by(tax1, shortname, abbrev) %>% summarize(hosts=length(unique(c(unique(tax2), unique(tax1)))), total_n=n())
g
summary(g)

h <- left_join(t, g)
h <- h %>% mutate(intratype = ifelse(channel=="txt", ifelse(as.character(taxa) == as.character(taxb), "txtintra", "txtinter"), "exp"))
h <- h %>% mutate(intratype = ifelse(intratype=="exp", ifelse(as.character(tax1) == as.character(tax2), "expintra", "expinter"), intratype))
h <- h %>% mutate(intra = ifelse(as.character(tax1) == as.character(tax2), "intra", "inter"))
h$intratype = factor(h$intratype)
h$intra = factor(h$intra)
h
summary(h)


# for i in *.fa; do BN=$(basename $i .fa); echo -n $BN; cat $i|grep ">"|wc -l ; done > protein_counts
counts <- read.table("~/Google Drive/Human-Virus-PPI/data/uniprot_to_payload/data_out/fastafiles/protein_counts", header=F)
colnames(counts) <- c("taxid", "proteins")
counts$taxid <- factor(counts$taxid)
h <- left_join(h, counts, by.x="taxa", by.y="taxid")
h$taxid <- factor(h$taxid)

f <- h %>% group_by(tax1, shortname, abbrev, channel, proteins) %>% summarize(n=n(), total_n=max(total_n))
f$taxid <- factor(f$taxid)

h <- h %>% mutate(type = ifelse(channel == "exp", "experimental", "text mining"))
ff <- h %>% group_by(tax1, shortname, abbrev, intra, channel ) %>% summarize(n=n(), total_n=max(total_n))
ff$taxid <- factor(ff$taxid)
ff$tax1 <- factor(ff$tax1)
ff$tax2 <- factor(ff$tax2)

hs <- h %>% mutate(conf = ifelse(score >= 900, ">=0.9", ifelse(score >= 700, ">=0.7", ifelse(score >= 400, ">=0.4", ">=0.15"))))
hs$conf = factor(hs$conf, levels=c('>=0.9','>=0.7','>=0.4','>=0.15'))

fs <- hs %>% group_by(taxid, shortname, abbrev, type, proteins, conf) %>% summarize(n=n(), total_n=max(total_n))
fis <- hs %>% group_by(taxid, shortname, abbrev, intratype, type, intra, proteins, conf) %>% summarize(n=n(), total_n=max(total_n))

v <- h %>% group_by(taxid, shortname, abbrev, total_n, proteins) %>% summarize()


# how many hosts 
hosts <- read.table("~/Google Drive/Human-Virus-PPI/data/uniprot_to_payload/data_out/hosts", header=F, sep="\t")
colnames(hosts) <- c("taxid", "host")
hosts$taxid <- factor(hosts$taxid)

ho <- left_join(v, hosts)
ho$taxid <- factor(ho$taxid)
# how many human viruses?
ho %>% filter(host == 9606) %>% group_by(shortname) %>% summarize()

ho <- ho %>% mutate(human = ifelse(host==9606, "human", ""))
ho <- ho %>% group_by(taxid, shortname, abbrev, total_n, proteins) %>% summarize(hosts = n(), human=max(human))


# how many non-human hosts?
length(setdiff( unique(as.vector(rbind(as.character(h$taxa), as.character(h$taxb)))), hosts$taxid)) - 1

# all rolled up
ggplot(v %>% filter(total_n>100), aes(y=reorder(abbrev, total_n), x=total_n)) + geom_point()

# BETTER PLOTS FOR PAPER
ggplot(ff %>% filter(total_n>2500), aes(x=reorder(shortname, total_n), y=n/2, fill=channel)) +  geom_bar(alpha=1, stat="identity") + theme_bw() + theme(axis.text.x=element_text(angle=0)) + scale_fill_manual(values = cols) + scale_color_manual(values=cols) + theme(axis.title.x=element_blank(), axis.title.y=element_blank()) + coord_flip() + theme(legend.position = c(0.9, 0.2))
ggsave("~/Google Drive/Human-Virus-PPI/data/virus_interactions/viruses_interactiontypes.pdf")

ggplot(l %>% filter(total_n>1000), aes(x=reorder(abbrev, total_n), y=score, col=channel)) +  geom_point(alpha=0.3) + theme_bw() + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + scale_fill_manual(values = cols) + scale_color_manual(values=cols) + theme(axis.title.x=element_blank(), axis.title.y=element_blank()) + coord_flip()

ggplot(l %>% filter(total_n>0), aes(x=channel, y=score, fill=channel)) +  geom_boxplot() + theme_bw() + theme() + scale_fill_manual(values = cols) + scale_color_manual(values=cols) + theme(axis.title.x=element_blank(), axis.title.y=element_blank()) + theme(legend.position="none")

ggplot(ffi , aes(x=intra, y=n/2, fill=channel)) +  geom_histogram(stat="identity") + theme_bw() + theme(axis.text.x=element_text(angle=0,hjust=0.5,vjust=0.5)) + scale_fill_manual(values = cols) + scale_color_manual(values=cols) + theme(axis.title.x=element_blank(), axis.title.y=element_blank())

ggplot(ffi , aes(x=channel, y=n/2, fill=channel)) +  geom_histogram(stat="identity") + theme_bw() + theme(axis.text.x=element_text(angle=0,hjust=0.5,vjust=0.5)) + scale_fill_manual(values = cols) + scale_color_manual(values=cols) + theme(axis.title.x=element_blank(), axis.title.y=element_blank())



# NUMBERS FOR THE TABLE
summar <- h %>% group_by(intratype) %>% summarize(n=n(), total_n=sum(total_n))

#average number of proteins in a virus


# with facet by confidence
cols <- c("#AE6962", "#99B4AF")
ggplot(fs %>% filter(total_n>100), aes(x=reorder(abbrev, total_n/proteins), y=n, col=type)) +  geom_point(alpha=0.8) + theme_bw() + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + scale_y_log10() + scale_fill_manual(values = cols) + scale_color_manual(values=cols) + facet_grid( . ~ conf )
ggsave("~/Google Drive/Human-Virus-PPI/data/virus_interactions/confidence.pdf")

# with facet by intra and confidence
cols <- c("#AE6962", "#99B4AF")
ggplot(fis %>% filter(total_n>30), aes(x=reorder(abbrev, total_n), y=n, col=intratype)) +  geom_point(alpha=0.8) + theme_bw() + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + scale_y_log10() + scale_fill_manual(values = cols) + scale_color_manual(values=cols) + facet_grid( conf ~ intra )
ggsave("~/Google Drive/Human-Virus-PPI/data/virus_interactions/intractions-intra-confidence.pdf")

# number of proteins per virus
ggplot(f %>% filter(total_n>100), aes(x=reorder(abbrev, total_n/proteins), y=1, label=proteins)) +  geom_text(size=3) + theme_bw() + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + scale_y_log10() + scale_fill_manual(values = cols) + scale_color_manual(values=cols)
ggsave("~/Google Drive/Human-Virus-PPI/data/virus_interactions/proteins.pdf")

# number of hosts
ggplot(ho %>% filter(total_n>100), aes(x=reorder(abbrev, total_n/proteins), y=1, label=hosts)) +  geom_text(size=3) + theme_bw() + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + scale_y_log10() + scale_fill_manual(values = cols) + scale_color_manual(values=cols) + facet_grid(human ~ .)
ggsave("~/Google Drive/Human-Virus-PPI/data/virus_interactions/hosts.pdf")

# total interactions
nrow(h)

