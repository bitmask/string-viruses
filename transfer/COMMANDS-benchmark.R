# get 9606 output from transfer and calibrated data
# sagittarius: /mnt/mnemo4/damian/STRING_v10_5/STRING_derived_v10_5/transfer/transfer_output/
# join the pre and post benchmark files together on id
for i in actions exper text; do cat $i|awk '{print $1 "-" $2 "-" $3 "\t" $6}' > $i.out; done
for i in actions exper text; do join $i.out ../orig/$i.out |awk '{print $2 "\t" $3}'|cut -c 1-11|sort| uniq> $i.benchmark; done

# actions needs to be joined on all fields
cat actions|awk '{print $1 "-" $2 "-" $3 "-" $4 "-" $5 "\t" $6}' > actions.out


# plot curves for each 
library(ggplot2)
actions <- read.table("~/Google Drive/Human-Virus-PPI/data/transfer/benchmark.actions")[c(1,7,2,3,5,14,6)]
colnames(actions) <- c('tax1', 'tax2', 'prot1', 'prot2', 'type', 'score', 'tscore')
ggplot(actions, aes(x=score, y=tscore*1000, col=type)) + geom_point() + geom_abline(intercept=0, slope=1)

experiments <- read.table("~/Google Drive/Human-Virus-PPI/data/transfer/benchmark.experiments")[c(1,7,2,3,5,14,6)]
colnames(experiments) <- c('tax1', 'tax2', 'prot1', 'prot2', 'type', 'score', 'tscore')
ggplot(rbind(actions, experiments), aes(x=score, y=tscore*1000, col=type)) + geom_point() + geom_abline(intercept=0, slope=1)

# get cellular benchmark from zurich, orig scores output from transfer and calibrated scores
actions.bench <- read.table("~/Google Drive/Human-Virus-PPI/data/transfer/9606.cellular.benchmark/actions.benchmark", sep="\t", header=FALSE)
colnames(actions.bench) <- c("calib", "orig")
ggplot(actions.bench, aes(x=orig, y=calib)) + geom_point(alpha=0.1) + xlim(c(0,1)) + ylim(c(0,1)) + geom_abline(intercept=.1, slope=0.125) + geom_abline(intercept=-0.45, slope=1.1) + geom_abline(intercept=-11, slope=12)
ggsave("~/Google Drive/Human-Virus-PPI/data/transfer/9606.cellular.benchmark/actions.pdf")

exper.bench <- read.table("~/Google Drive/Human-Virus-PPI/data/transfer/9606.cellular.benchmark/exper.benchmark", sep="\t", header=FALSE)
colnames(exper.bench) <- c("calib", "orig")
ggplot(exper.bench, aes(x=orig, y=calib)) + geom_point() + xlim(c(0,1)) + ylim(c(0,1)) + geom_abline(intercept=.02, slope=0.25) + geom_abline(intercept=-0.95, slope=2)
ggsave("~/Google Drive/Human-Virus-PPI/data/transfer/9606.cellular.benchmark/exper.pdf")

text.bench <- read.table("~/Google Drive/Human-Virus-PPI/data/transfer/9606.cellular.benchmark/text.benchmark", sep="\t", header=FALSE)
colnames(text.bench) <- c("calib", "orig")
ggplot(text.bench, aes(x=orig, y=calib)) + geom_point() + xlim(c(0,1)) + ylim(c(0,1)) + geom_abline(intercept=0, slope=0.25) + geom_abline(intercept=-0.54, slope=1.3)
ggsave("~/Google Drive/Human-Virus-PPI/data/transfer/9606.cellular.benchmark/text.pdf")
