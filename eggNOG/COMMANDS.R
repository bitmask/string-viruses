#install.packages("dplyr")
library(dplyr)
library(ggplot2)
setwd("~/Google Drive/Human-Virus-PPI/data/eggNOG/clades-indiv-xml")
setwd("~/Google Drive/Human-Virus-PPI/data/eggNOG/clades-indiv-proteomes")


#nog = orthology group from eggNOG
#mog = merged nogs

# statistics for merging
# "how many nogs are merged into mogs?"

details <- read.table("10239.out.eggnog.final.details", header=FALSE, sep="\t")
#details <- read.table("10239.out.details", header=FALSE, sep="\t")
summary(details)
colnames(details) <- c("mog", "nog", "stringid", "descr", "pfam", "clan", "tax")
details$id <- as.integer(gsub("MOG", "", details$mog))
nog_details <- tbl_df(details) %>% group_by(nog) %>% summarize(mog=first(mog), proteinspernog=n(), id=first(id)) %>% group_by(mog) %>% summarize(nogspermog = n(), proteinspermog = sum(proteinspernog), id=first(id))
arrange(nog_details, -nogspermog)
ggplot(nog_details, aes(x=nogspermog)) + geom_histogram() + theme_bw() + xlab("NOGs per MOG")
ggsave("nogs_per_mog.pdf")
ggplot(nog_details, aes(x=proteinspermog)) + geom_histogram(binwidth=1) + theme_bw() + xlab("Proteins per MOG")
filter(nog_details, nogspermog >1) %>% arrange(-nogspermog)



# statistics for length information

# for NOGS

length_nog_stats <- read.table("10239.out.nog.length", header=TRUE)
length_stats <- length_nog_stats

length_mog_stats <- read.table("10239.out.length", header=TRUE)
length_stats <- length_mog_stats

length_df <- tbl_df(length_stats)
summary(length_df)
ggplot(length_df, aes(x=reorder(nog, length), y=length)) + geom_boxplot() + xlab("NOG")
ggsave("NOG_length_boxplot.pdf")
ggplot(length_df, aes(x=reorder(nog, length), y=length)) + geom_boxplot() + xlab("MOG")
ggsave("MOG_length_boxplot.pdf")


length_df <- mutate(length_df, scaled = length / mean(length))


length_summary <- group_by(length_df, nog) %>% summarize( smin=min(length/mean(length)), smean=mean(length/mean(length)), smax=max(length/mean(length)), ssd=sd(length/mean(length)), min=min(length), mean=mean(length), max=max(length), sd=sd(length), size=n() )
#length_summary <- group_by(length_df, nog) %>% summarize(min=min(length), max=max(length), mean=mean(length), sd=sd(length), size=n()) %>% mutate(z = sd/mean) 
summary(length_summary)
# distribution of range of lengths



ggplot(length_summary, aes(x=mean, y=sd)) + geom_point(alpha=0.2) + theme_bw() + scale_x_log10() + scale_y_log10() + xlab("Mean of NOG protein lengths") + ylab("Standard deviation of NOG protein lengths")
ggsave("NOG_length_log.pdf")
ggplot(length_summary, aes(x=mean, y=sd)) + geom_point(alpha=0.2) + theme_bw() + scale_x_log10() + scale_y_log10() + xlab("Mean of MOG protein lengths") + ylab("Standard deviation of MOG protein lengths")
ggsave("MOG_length_log.pdf")


ggplot(length_summary, aes(x=size, y=ssd)) + geom_point(alpha=0.2) + scale_y_log10() + scale_x_log10()
ggplot(length_summary, aes(x=size, y=max/min)) + geom_point(alpha=0.2) + scale_y_log10() + scale_x_log10()

# size of each mog in proteins
#ggplot(length_summary, aes(x=reorder(nog, size), y=size)) + geom_point(alpha=0.2) + scale_y_log10()
#arrange(length_summary, -size)

filter(length_summary, ssd > .75)
filter(length_df, nog=="NOG06292")

ggplot(length_summary, aes(x=reorder(nog, ssd), y=ssd)) + geom_point(alpha=0.2) + scale_y_log10()

arrange(length_summary, -ssd)

length_summary_high <- length_summary %>% filter( z > 0.15)
length_summary_low <- length_summary %>% filter( z <= 0.15)
ggplot(filter(length_df, nog %in% length_summary_high$nog), aes(x=reorder(nog, length), y=length)) + geom_boxplot()

arrange(length_summary, -sd)






# statistics for architecture information

arch_nog_stats <- read.table("10239.out.nog.architectures", header=T)
arch_stats <- arch_nog_stats
arch_df <- tbl_df(arch_stats)
summary(arch_df)
ggplot(arch_df, aes(x=unique_clan)) + geom_histogram() + theme_bw() + xlab("Unique clans per NOG")
ggsave("unique_clans_per_NOG.pdf")

arch_nog_stats <- read.table("10239.out.architectures", header=T)
arch_stats <- arch_nog_stats
arch_df <- tbl_df(arch_stats)
summary(arch_df)
ggplot(arch_df, aes(x=unique_clan)) + geom_histogram() + theme_bw() + xlab("Unique clans per MOG")
ggsave("unique_clans_per_MOG.pdf")
ggplot(arch_df, aes(y=unique_clan/total_entries, x=reorder(name, unique_clan/total_entries))) + geom_point(alpha=0.2)

# slope shows the number of unique clans that are within each group; would like all groups to be on x=y, if they are much lower then there are too many clans for the size of the group
ggplot(arch_df, aes(y=total_entries/unique_clan, x=total_entries)) + geom_point(alpha=0.2)

ggplot(arch_df, aes(x=reorder(name, total_entries), y=unique_clan)) + geom_point(alpha=0.2)

arrange(arch_df, -unique_arch)
filter(arch_df, unique_clan / total_entries > 0.5 & total_entries > 2)







# statistics for depth information

depth_nog_stats <- read.table("10239.out.nog.depth", header=F, sep="\t")
depth_stats <- depth_nog_stats

depth_stats <- read.table("10239.out.depth", header=F, sep="\t")


colnames(depth_stats) <- c("nog", "depth", "details")
depth_df <- tbl_df(depth_stats)
arrange(depth_df, depth)
arrange(depth_df, -depth)
group_by(depth_df, depth) %>% summarize(count=n())
ggplot(depth_df, aes(x=depth)) + geom_histogram()
ggsave("NOG_depth.pdf")
ggsave("MOG_depth.pdf")







# statistics for consistency 

parent <- read.table("10239.out.details", header=F, sep="\t")[,c(1,2,3)]
colnames(parent) <- c("mog", "nog", "stringid")
child <- read.table("29258.out.details", header=F, sep="\t")[,c(1,2,3)]
colnames(child) <- c("mog", "nog", "stringid")

parent_df <- tbl_df(parent)
group_by(parent_df, mog) %>% summarize(stringid=concatenate(stringid))




# OLD


# statistics for clan information

stats <- read.table("10239.out.architectures", header=TRUE)
summary(stats)
all <- tbl_df(stats) %>% mutate(perc = unique_clan / total_entries)

# select the MOGs that have the highest percentage of unique architectures on clan basis
suspect <- filter(all, (perc > 0.50 & total_entries > 2) | (unique_clan > 3))
# number of them as percentage of all MOGs
100 * length(rownames(suspect)) / length(rownames(all))


hist(stats$unique_arch, breaks=60)
hist(stats$unique_clan, breaks=7, xlab="Unique clan architectures per orthology group")
plot(stats$unique_clan, stats$total_entries)
hist(stats$total_entries, breaks=70, xlab="Proteins per merged orthology group")



# other

hist.data = hist(stats$unique, plot=F)
hist.data$counts = log10(hist.data$counts)

dev.new(width=4, height=4)
hist(x)

dev.new(width=4, height=4)
plot(hist.data, ylab='log10(Frequency)')

a <- read.table("~/Desktop/
