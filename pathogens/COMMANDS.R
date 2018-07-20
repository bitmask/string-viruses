library(dplyr)
library(ggplot2)

# prepare the list of host:pathogen pairs

d <- read.table("~/Google Drive/Human-Virus-PPI/data/pathogen-host/host_pathogen_allsources_domains", quote="", sep="\t")
d <- tbl_df(d)
colnames(d) <- c('host', 'host_domain', 'pathogen', 'pathogen_domain')
d %>% distinct() %>% group_by(host_domain, pathogen_domain) %>% summarize(n=n()) %>% print(n=200)

d %>% distinct() %>% filter(host == 9606) %>% group_by(host_domain, pathogen_domain) %>% summarize(n=n())

# get names
pathogens <- read.table("~/Google Drive/Human-Virus-PPI/data/pathogen-host/pathogens.names", header=FALSE, sep="\t", quote="", fill=TRUE)
colnames(pathogens) = c("taxid", "pathogen_name")
pathogens <- tbl_df(pathogens)
hosts <- read.table("~/Google Drive/Human-Virus-PPI/data/pathogen-host/hosts.names", header=FALSE, sep="\t", quote="")
colnames(hosts) = c("taxid", "host_name")
hosts <- tbl_df(hosts)


# filter by tax domain
f <- d %>% filter(host_domain != "none", pathogen_domain != "none") %>%  # remove crap
   filter(pathogen_domain != "viruses") %>%  # remove viruses, because their host-virus list comes from uniprot
   filter(host_domain != "bacteria") %>% # remove bacteria as a host because these are most likely just wrong
   filter(host_domain != "vertebrates" & pathogen_domain!="vertebrates")  # remove bats as pathogen of human, and human as pathogen of chickens

# filter specific model organisms out (from list of model organisms on wikipedia)
f <- f %>% filter(
            pathogen!= 511145, # e coli mg1655
            pathogen!= 83333, # e coli k-12
            pathogen!= 1423, # b subtilis
            pathogen!= 155892, # Caulobacter crescentus
            pathogen!= 2097, # Mycoplasma genitalium
            pathogen!= 668, #Aliivibrio fischeri
            pathogen!= 1142, #Synechocystis --- need to climb tree to here
            pathogen!= 294, #Pseudomonas fluorescens
            pathogen!= 354, #Azotobacter vinelandii
            pathogen!= 1902, #Streptomyces coelicolor

            pathogen!= 3055, #Chlamydomonas reinhardtii
            pathogen!= 5963, #Stentor coeruleus
            pathogen!= 44689, #Dictyostelium discoideum
            pathogen!= 5911, #Tetrahymena thermophila
            pathogen!= 2903, #Emiliania huxleyi
            pathogen!= 35128, #Thalassiosira pseudonana

            pathogen!= 33169, #Eremothecium gossypii
            pathogen!= 162425, #Aspergillus nidulans
            pathogen!= 5346, #Coprinopsis cinerea
            pathogen!= 5207, #Cryptococcus neoformans
            pathogen!= 5141, #Neurospora crassa
            pathogen!= 4932, #Saccharomyces cerevisiae
            pathogen!= 559292, #Saccharomyces cerevisiae S288C
            pathogen!= 5334, #Schizophyllum commune
            pathogen!= 4896, #Schizosaccharomyces pombe
            pathogen!= 284812, #Schizosaccharomyces pombe 972h-
            pathogen!= 5270 #Ustilago maydis
            ) 

f <- f %>% left_join(hosts, by=c("host"="taxid"))
f <- f %>% left_join(pathogens, by=c("pathogen"="taxid"))
f %>% select(host, host_name, host_domain, pathogen, pathogen_name, pathogen_domain) %>% print(n=1000)

f %>% select(host, pathogen) %>% write.table(file="host_pathogen_trimmed", sep="\t", col.names=F, row.names=F)

string <- read.table("~/Google Drive/Human-Virus-PPI/data/pathogen-host/string_v11_species_tospecies", quote="", sep="\t")
colnames(string) <- c("instring", "species")
instring <- string$instring
species <- string[string$species != 1,]$species
allstring <- union(instring, species) # contains any taxid in string, or it rolled up to species level


# you need to run climb_tax_tree on output above to generate these files

hpt1 <- read.table("~/Google Drive/Human-Virus-PPI/data/pathogen-host/hpt1_species", sep="\t")
colnames(hpt1) <- c("host", "hostspecies")
hpt1 <- tbl_df(hpt1) %>% filter(hostspecies != 1)
hpt2 <- read.table("~/Google Drive/Human-Virus-PPI/data/pathogen-host/hpt2_species", sep="\t")
colnames(hpt2) <- c("pathogen", "pathogenspecies")
hpt2 <- tbl_df(hpt2) %>% filter(pathogenspecies != 1)

# 
f1 <- left_join(f, hpt1)
f2 <- left_join(f1, hpt2)

a <- f2 %>% filter(host %in% allstring & pathogen %in% allstring) %>% select(host, pathogen)
b <- f2 %>% filter(host %in% allstring & pathogenspecies %in% allstring) %>% select(host, pathogenspecies)
colnames(b) <- c("host", "pathogen")
c <- f2 %>% filter(hostspecies %in% allstring & pathogen %in% allstring) %>% select(hostspecies, pathogen)
colnames(c) <- c("host", "pathogen")
d <- f2 %>% filter(hostspecies %in% allstring & pathogenspecies %in% allstring) %>% select(hostspecies, pathogenspecies)
colnames(d) <- c("host", "pathogen")

rbind(a, b, c, d) %>% group_by(host, pathogen) %>% summarize() %>% select(pathogen, host) %>% write.table(file="species_selection_host_pathogen", sep="\t", col.names=F, row.names=F)

# this needs to be combined with the viruses file, and it's done.







sum <- f %>% group_by(host_domain, pathogen_domain) %>% summarize(n=n()) 
sum %>% print(n=200)

sum$host_domain <- factor(sum$host_domain, levels=c("human", "animals", "plants", "fungi", "eukaryota", "bacteria", "archaea", "none"))
sub$pathogen_domain <- factor(sum$pathogen_domain, levels=c("bacteria", "archaea", "fungi", "apicomplexa (protists)", "kinetoplastida (protists)", "amoebozoa", "annelida (worms)", "nematoda (worms)", "platyhelminthes (worms)", "neoptera (fleas, louse)", "eukaryota", "viruses", "none"))



ggplot(sum, aes(x=host_domain, y=pathogen_domain, label=n)) + geom_text()




# textmining results

tax <- read.table("~/Google Drive/Human-Virus-PPI/data/pathogen-host/host_pathogen_allsources_domains", header=FALSE, sep="\t", quote="")
hosttax <- tax[,1:2] %>% tbl_df() %>% distinct()
colnames(hosttax) <- c("taxid", "host_domain")
pathtax <- tax[,3:4] %>% tbl_df() %>% distinct()
colnames(pathtax) <- c("taxid", "pathogen_domain")


pathogens <- left_join(pathogens, pathtax)
hosts <- left_join(hosts, hosttax)



tm <- read.table("~/Google Drive/Human-Virus-PPI/data/pathogen-host/textmining/pathogen_host_network_hpidbonly.r", header=FALSE, sep="\t")
colnames(tm) <- c("tax1", "tax2", "channel", "score")
tm <- tbl_df(tm) %>% filter(score > 0, tax1!=tax2)
summary(tm)

tmn <- left_join(tm, pathogens, by=c("tax1"="taxid"))
tmn <- left_join(tmn, pathogens, by=c("tax2"="taxid"))
tmn <- tmn %>% mutate(pathogen_name = ifelse(is.na(pathogen_name.x), as.character(pathogen_name.y), as.character(pathogen_name.x)))
tmn$pathogen_name.x <- NULL
tmn$pathogen_name.y <- NULL
tmn <- tmn %>% mutate(pathogen_domain = ifelse(is.na(pathogen_domain.x), as.character(pathogen_domain.y), as.character(pathogen_domain.x)))
tmn$pathogen_domain.x <- NULL
tmn$pathogen_domain.y <- NULL

tmn <- left_join(tmn, hosts, by=c("tax1"="taxid"))
tmn <- left_join(tmn, hosts, by=c("tax2"="taxid"))
tmn <- tmn %>% mutate(host_name = ifelse(is.na(host_name.x), as.character(host_name.y), as.character(host_name.x)))
tmn$host_name.x <- NULL
tmn$host_name.y <- NULL
tmn <- tmn %>% mutate(host_domain = ifelse(is.na(host_domain.x), as.character(host_domain.y), as.character(host_domain.x)))
tmn$host_domain.x <- NULL
tmn$host_domain.y <- NULL

tmn %>% group_by(pathogen_name, pathogen_domain, host_name, host_domain) %>% summarize(n=n()) %>% filter(!grepl("virus", pathogen_name)) %>% arrange(-n) %>% print(n=100)

dom <- tmn %>% group_by(pathogen_domain, host_domain) %>% summarize(n=n()) %>%  arrange(-n) %>% print(n=100)
ggplot(dom %>% filter(pathogen_domain != "viruses"), aes(x=host_domain, y=pathogen_domain, label=n, fill=log(n))) + geom_tile(aes(fill=log(n))) + scale_fill_gradient(low="white", high="steelblue") + geom_text() + theme_minimal() 


ggplot(tm, aes(x=score)) + geom_histogram(binwidth=.01)
