#install.packages("Pstat")
install.packages("caret")
devtools::install_github('spflanagan/gwscaR')
library(Pstat)
library(gwscaR)
library(ade4)
library(caret)

# data(test)
# Pst(test)
# TracePst(test,va=0,ci=1,boot=10,pe=0.95,Fst=-1,,xm=2,pts=30)

#myPCA <- gm.prcomp(lent.gpa$coords)
setwd("~/Documents/Desktop_Items/Research/")

#PST for the Top View

# Read in the data
PCvals <- read.csv("/Users/QuinnThomas/Documents/Desktop_Items/Research/PCvals.csv", header=TRUE)
beakVals <- read.csv("/Users/QuinnThomas/Documents/Desktop_Items/Research/beakVals.csv", header=TRUE)
#Read in metadata
indNames <- read.csv("/Users/QuinnThomas/Documents/Desktop_Items/Research/sitePlants.csv", header=TRUE)
varNames <- read.csv("/Users/QuinnThomas/Documents/Desktop_Items/Research/varLabels.csv", header=TRUE)

#condense beak opening and beak depth data by plant
beakOpening <- aggregate( beakOpening ~ id, beakVals, mean )
beakDepth <- aggregate( beakDepth ~ id, beakVals, mean )

# checkpoint
length(indNames)
length(varNames)
length(PCvals[,1])

# format the individuals vecctor to extract sites
inds <- sapply(strsplit(indNames[,2], "-"), '[', 2)
sites2 <- gsub('.{1}$', '', inds)
sites <- paste0("S_", sites2)

# The commented out code was already run and saved as a csv file which I read in the uncommented code

# read in K pop structure
# K_map <- read.delim("sampleNames_popmap2.2.csv", header=T, sep=",")
# K_inds <- sapply(strsplit(K_map[,1], "-"), '[', 1)
# K_map[,1] <- K_inds
# K_map <- K_map[,c(1,9)]
# K_csv <- left_join(as.data.frame(inds), as.data.frame(K_map), by=c("inds"="Sample_ID"))
# write.csv(K_csv, file="PCA_K.csv")
# I made some adjustments based on assumptions
K_pops <- read.csv(file="PCA_K_inf.csv", header=T)
inf_K <- K_pops[,2]

# ?left_join
# ?strsplit

# make a dataframe with one column of sample names, one column of population identification 
# and columns for each relevant PC (maybe just the first two or three)

#I propose we calculate Pst for a few different "populations"

# 1. proposed genetic populations according to structure
  # Avg FST = 0.084
# 2. variety populations
  # Avg FST = 0.127
# 3. site populations
  # Avg FST = 0.269

#Just a way to check the breakdown of data
summary_var <- table(varNames[,2]) # min 3, max fremontii
summary_site <- table(sites) # min 1, max 21
summary_K <- table(inf_K) # min 4, max 120

Var.df <- data.frame(sample=indNames[,2], pop=varNames[,2], site=sites, K_group=inf_K, PC1=PCvals[,1], PC2=PCvals[,2], PC3=PCvals[,3], beakOpening=beakOpening[,2], beakDepth=beakDepth[,2])
write.csv(Var.df, file="~/Documents/Desktop_Items/Research/PST_df.csv")
#Var.df$site <- as.character(Var.df$site)
#Var.df$K_group <- as.character(Var.df$K_group)

var_frame <- Var.df[,c(2,5:9)]
site_frame <- Var.df[,c(3,5:9)]
K_frame <- Var.df[,c(4:9)]
length(site_frame[,1])
length(site_frame[,2])
length(var_frame[,1])
length(var_frame[,2])
length(K_frame[,1])
length(K_frame[,2])

levels(as.factor(site_frame[,1]))
?pairwise.pst

var_pst_pc1 <- pairwise.pst(var_frame[1:2], levels(as.factor(var_frame[,1])))
var_pst_pc2 <- pairwise.pst(var_frame[,c(1,3)], levels(as.factor(var_frame[,1])))
var_pst_pc3 <- pairwise.pst(var_frame[,c(1,4)], levels(as.factor(var_frame[,1])))
var_pst_bo <- pairwise.pst(var_frame[,c(1,5)], levels(as.factor(var_frame[,1])))
var_pst_bd <- pairwise.pst(var_frame[,c(1,6)], levels(as.factor(var_frame[,1])))

mean(var_pst_pc1)

table(site_frame[,1])
#site_frame2 <- site_frame[-c(24,187,180,278),]
site_pst_pc1 <- pairwise.pst(site_frame[,c(1,2)], levels(as.factor(site_frame[,1])))
site_pst_pc2 <- pairwise.pst(site_frame[,c(1,3)], levels(as.factor(site_frame[,1])))
site_pst_pc3 <- pairwise.pst(site_frame[,c(1,4)], levels(as.factor(site_frame[,1])))
site_pst_bo <- pairwise.pst(site_frame[,c(1,5)], levels(as.factor(site_frame[,1])))
site_pst_bd <- pairwise.pst(site_frame[,c(1,6)], levels(as.factor(site_frame[,1])))

sitepst <- BootPst(site_frame[,c(1,2)], va = 1,  csh=1, Rp=low_count_site)
data(test)

table(K_frame[,1])
K_pst_pc1 <- pairwise.pst(K_frame[,c(1,2)], levels(as.factor(K_frame[,1])))
K_pst_pc2 <- pairwise.pst(K_frame[,c(1,2)], levels(as.factor(K_frame[,1])))
K_pst_pc3 <- pairwise.pst(K_frame[,c(1,2)], levels(as.factor(K_frame[,1])))
K_pst_bo <- pairwise.pst(K_frame[,c(1,2)], levels(as.factor(K_frame[,1])))
K_pst_bd <- pairwise.pst(K_frame[,c(1,2)], levels(as.factor(K_frame[,1])))

# remove K inds with NA
K_frame_f <- K_frame[-c(24, 203:206),]
  
#considering all populations
#do this for c/h2 0.25

low_count_var <- c("araneosus","australis","floribundus") #removing populations with <10 individuals
low_count_site <- c("S_15","S_54","S_59", "S_77", "S_66", "S_61", "S_58", "S_18") #removing populations with <2 individuals
low_count_K <- c("2") #removing populations with <10 individuals

variables <- c("PC1", "PC2", "PC3", "beakOpening", "beakDepth")

?TracePst
#practice <- Pst(var_frame, ci=1, csh=1, Pw=c("maricopae", "wilsonii"))

#PST according to variety
TracePst(var_frame, ci=1, boot=10, Fst= 0.258, Rp=low_count_var)
#PST according to K
TracePst(K_frame_f, ci=1, boot=10, Fst=0.194, Rp=low_count_K)
#PST according to variety
TracePst(site_frame, ci=1, boot=10, Fst=0.315, Rp=low_count_site)


?mantel.randtest()

t_var_pst_pc1 <- t(var_pst_pc1)
t_site_pst_pc1 <- t(site_pst_pc1)
t_K_pst_pc1 <- t(K_pst_pc1)

plot(hist(as.dist(t_var_pst_pc1)))
plot(hist(as.dist(t_site_pst_pc1)))
plot(hist(as.dist(t_K_pst_pc1)))

pstAdj <- function(x) {
  return(x/(1-x))
}

# process <- preProcess(t_var_pst_pc1, method=c("range"))
# norm_scale <- predict(process, t_var_pst_pc1)
# 
scale_data <- as.data.frame(scale(t_var_pst_pc1))
scale_data_site <- as.data.frame(scale(t_site_pst_pc1))
# 
tadj_var_pst_pc1 <-log(t_var_pst_pc1)
tadj_site_pst_pc1 <-log(t_site_pst_pc1)

understand <- as.dist(ordered_fst_Hudson)
understand2 <- as.dist(t_site_pst_pc1)
row.names(ordered_fst_Hudson)
row.names(t_site_pst_pc1)

mantel.randtest(as.dist(ordered_fst_WC), as.dist(t_site_pst_pc1))

plot(hist(as.dist(norm_scale)))
plot(hist(as.dist(scale_data)))
plot(hist(as.dist(scale_data_site)))


# For a site comparison we have to remove everything the two matrices don't have in common
remove_fst <- c("1", "17", "25", "29", "34", "35", "36", "37", "39", "47", "51", 
                "55", "56", "57", "62", "63", "8", "9")
remove_pst <- c("S_15","S_2", "S_20", "S_21", "S_4", "S_67", "S_76", "S_77")

new_fst_Hud <- ordered_fst_Hudson[!row.names(ordered_fst_Hudson) %in% remove_fst, !row.names(ordered_fst_Hudson) %in% remove_fst]
new_pst_site <- t_site_pst_pc1[!row.names(t_site_pst_pc1) %in% remove_pst, !row.names(t_site_pst_pc1) %in% remove_pst]


newV.df <- geomorph.data.frame(lent.gpa, variety = beak_meta2[,3], sites= beak_meta2[,2], sample=beak_meta[,1],
                               beakWidth = beakOpenings, beakDepth = beakDepths)   


mantel.randtest(as.dist(new_fst_Hud), as.dist(t_site_pst_pc1))
gl.ibd(Dgen = as.dist(new_fst_Hud), Dgeo = as.dist(new_pst_site))

View(cbind(row.names(new_fst_Hud), row.names(new_pst_site)))

getwd()
mean_pst_fst <- read.csv("site_pst.csv")
mean_pst_fst[,1] <- factor(mean_pst_fst[,1], levels=mean_pst_fst[,1])

dev.off()
ggplot(mean_pst_fst, aes(Variable, Mean)) +        # ggplot2 plot with confidence intervals
  geom_point() +
  geom_errorbar(aes(ymin = LC, ymax = UC))













#PST for the Top View

# Read in the data
side_PCvals <- read.csv("/Users/QuinnThomas/Documents/Desktop_Items/Research/Side_PCvals.csv", header=TRUE)
#Read in metadata
side_indNames <- read.csv("/Users/QuinnThomas/Documents/Desktop_Items/Research/indNamesSide.csv", header=TRUE)
side_varNames <- read.csv("/Users/QuinnThomas/Documents/Desktop_Items/Research/varNamesSide.csv", header=TRUE)

side_varNames[c(27:39,44:53),2] <- "fremontii"
side_varNames[40:43,2] <- "borreganus"
side_sitesp <- sapply(strsplit(side_indNames[241:267,2], "-"), '[', 2)
side_sitesp2 <- c(side_indNames[1:240,2],side_sitesp)
side_sitesp2[243] <- "76E"
side_sites <- gsub('.{1}$', '', side_sitesp2)

# The commented out code was already run and saved as a csv file which I read in the uncommented code

# read in K pop structure
# side_K_map <- read.delim("popmap_K5.txt", header=F, sep="\t")
# side_K_inds <- sapply(strsplit(side_K_map[,1], "-"), '[', 1)
# side_K_map[,1] <- side_K_inds
# side_K_csv <- left_join(as.data.frame(side_indNames[,2]), as.data.frame(side_K_map), by=c("side_indNames[, 2]"="V1"))
# write.csv(side_K_csv, file="side_PCA_K.csv")
# I made some adjustments based on assumptions
side_K_pops <- read.csv(file="side_PCA_K_inf.csv", header=T)
side_inf_K <- side_K_pops[,2]

summary_side_var <- table(side_varNames[,2]) # min 3, max fremontii
summary_side_site <- table(side_sites) # min 1, max 21
summary_side_K <- table(side_inf_K) # min 4, max 120

Side.df <- data.frame(sample=side_indNames[,2], pop=side_varNames[,2], site=side_sites, K_group=side_inf_K, PC1=side_PCvals[,1], PC2=side_PCvals[,2], PC3=side_PCvals[,3])
side_var_frame <- Side.df[,c(2,5:7)]
side_site_frame <- Side.df[,c(3,5:7)]
side_K_frame <- Side.df[,c(4:7)]

side_site_pst_pc1 <- pairwise.pst(side_site_frame[,c(1,2)], levels(as.factor(side_site_frame[,1])))
side_site_pst_pc2 <- pairwise.pst(side_site_frame[,c(1,3)], levels(as.factor(side_site_frame[,1])))
side_site_pst_pc3 <- pairwise.pst(side_site_frame[,c(1,4)], levels(as.factor(side_site_frame[,1])))

Pst(side_site_frame[,c(1,3)], ci=1, csh=0.2)

s_mean_pst_fst <- read.csv("side_site_pst.csv")
s_mean_pst_fst[,1] <- factor(s_mean_pst_fst[,1], levels=s_mean_pst_fst[,1])

dev.off()
ggplot(s_mean_pst_fst, aes(Variable, Mean)) +        # ggplot2 plot with confidence intervals
  geom_point() +
  geom_errorbar(aes(ymin = LC, ymax = UC))


# remove K inds with NA
side_K_frame_f <- side_K_frame[-c(24, 191:194),]

side_low_count_var <- c("araneosus","australis","floribundus", "kennedyi") #removing populations with <10 individuals
side_low_count_site <- c("15","54","58", "77", "66", "61", "69", "18", "71") #removing populations with <2 individuals
side_low_count_K <- c("2") #removing populations with <10 individuals

#PST according to variety
TracePst(side_var_frame, ci=1, boot=10, Fst=0.127, Rp=side_low_count_var)
#PST according to K
TracePst(side_K_frame_f, ci=1, boot=10, Fst=0.084, Rp=side_low_count_K)
#PST according to variety
TracePst(side_site_frame, ci=1, boot=10, Fst=0.315, Rp=side_low_count_site)

