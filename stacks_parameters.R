setwd("~/Downloads/new_subset_parameters/")
library(ggplot2)
library(RColorBrewer)
library(reshape2)

############################################################################
  # Stacks parameter check
############################################################################

# i) total number of assembled loci, simply by summing all values of the second column; 
# ii) the number of polymorphic loci, by summing all all values of the second column, 
    # except the first one (which tells you how many loci are monomorphic); 
# iii) total number of SNPs, by multiplying the second column with the first column.


count <- 1
files <- list.files(pattern="*_snp_distribution.txt", full.names = T)

for (i in files[1:9]){
  table <- read.delim(i, skip=1, header=T)
  table$n_loci_percent<- table$n_loci/sum(table$n_loci)
  table$m<- count
  table$total_snps <- table$n_snps*table$n_loci
  write.table(table, "distributions.tsv", append=T, row.names=F, col.names = F)
  snp_count <- data.frame("m"= count, "n_snps"=sum(table$n_loci))
  write.table(snp_count, "total_count.tsv", append=T, row.names=F, col.names = F)

  numrow <- nrow(table)
  table_sum <- data.frame(m=count ,total_loci = sum(table$n_loci), total_poly_loci=sum(table$n_loci[2:numrow]))
  write.table(table_sum, "summary_distributions.tsv", append=T, row.names=F, col.names = F)
  count <- count + 1
  
}

library(ggplot2)

snp_count<-read.delim("total_count.tsv", sep=" ", header=F)
names(snp_count)<-c("M", "n_snps")
snp_count$M<-as.factor(snp_count$M)
#snp_count[,2] <- c(6318, 7868, 8270, 8323, 8281, 8220, 8151, 8077, 8034)
ggplot(snp_count, aes(x=M, y=n_snps))+geom_point()

#p<-ggplot(data=snp_count, aes(x=M, y=n_snps)) +
#  geom_point() + theme_classic() + scale_y_continuous(limits = c(0, 3000))
#p


loci_counts<-read.delim("summary_distributions.tsv", sep=" ", header=F)
names(loci_counts)<-c("M", "total_loci", "total_polymorphic_loci")
loci_counts <- melt(loci_counts,  id.vars = 'M', variable.name = 'loci')
ggplot(loci_counts, aes(M, value)) +
  geom_line(aes(colour = loci))+scale_x_continuous(labels=as.character(loci_counts$M),breaks=loci_counts$M)


snp_table<-read.delim("distributions.tsv", sep=" ", header=F)
names(snp_table)<- c("n_snps","n_loci", "n_loci_percent", "m") 
snp_table$n_loci_percent<-snp_table$n_loci_percent*100
snp_table$n_snps<-ifelse(snp_table$n_snps < 9, snp_table$n_snps, "9 +")
snp_table$n_snps<-as.factor(snp_table$n_snps)
snp_table$m<-as.factor(snp_table$m)

q<-ggplot(data = snp_table) + 
  geom_col(aes(x=n_snps, y=n_loci_percent, fill=m), position="dodge") + scale_fill_brewer(palette = "PuBu")+ theme_grey()
q

############################################################################
# PCA plots
############################################################################





