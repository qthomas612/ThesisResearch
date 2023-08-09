library(dplyr)
library(data.table)
library(adegenet)
library(ape)
library(pegas)
library("seqinr")
library("ggplot2")
library(devtools)
library(RColorBrewer)
library(dartR)
#install.packages("geosphere")                # Install geosphere package
library("geosphere") 
# install.packages("poppr")
library(poppr)
library(ade4)
#install.packages("graph4lg")
library(graph4lg)

#install_github("jgx65/hierfstat")
library("hierfstat")


# Here we are going to: read in our Master info sheet,
# Subset it by plant, site, variety and GPS decimal degrees
# join our population map with this information

setwd("~/Documents/Desktop_Items/Research/")

#read in all sample info
masterList <- read.csv("MasterPlant.csv")
#subset by sample, site, var and coordinates
subsetCols <- data.frame(masterList[,c(3,12,13)])
names(subsetCols)[names(subsetCols) == 'Sample_ID'] <- 'Sample'

setwd("~/Downloads/2023_lent/May_lent/")

#read in our population map of samples we have SNP data for
popmap <- read.csv("May_K2.csv", header=T)

#now we want to do a left join with popmap
all_meta_data <- left_join(popmap, subsetCols, by=c("Sample_ID"="Sample"))
#all_meta_data2 <- all_meta_data[order(all_meta_data$Sample_ID),]
#all_meta_data[32,1] <- "20E-2"
#all_meta_data[29,1] <- "1E-1"
#all_meta_data <- all_meta_data[-124,] #remove 43L
all_meta_data$varFactor <- as.numeric(as.factor(all_meta_data[,2]))

all_meta_data2 <- all_meta_data[order(all_meta_data$Sample_ID),]

#write.table(all_meta_data[,c(1,7)], file="Dec15_popmap1_gp_var.txt", quote=F, sep="\t", row.names = F)
#write.table(all_meta_data[,c(1,3)], file="Dec15_popmap1_gp_site.txt", quote=F, sep="\t", row.names = F)

all_site_coords <- as.data.table(all_meta_data[,c(3,7,8)])
all_site_coords2 <- data.frame(all_site_coords)
all_var_coords <- as.data.table(all_meta_data[,c(2,7,8)])
# insert nigricalycis: (30E) 35.10335, -119.418467 and (33D) 34.93875, -118.933783
all_var_coords[16,2] <- 35.10335
all_var_coords[16,3] <- -119.418467
all_var_coords[17,2] <- 34.93875
all_var_coords[17,3] <- -118.933783
all_ind_coords <- as.data.table(all_meta_data[,c(1,7,8)])


latvec <-c()
longvec<-c()
for (x in all_site_coords2$Site){
  lat =all_site_coords2[all_site_coords2$Site == x, "Lat"][1]
  lenlat = length(all_site_coords2[all_site_coords2$Site == x, "Lat"])
  long = all_site_coords2[all_site_coords2$Site == x, "Long"][1]
  lenlong = length(all_site_coords2[all_site_coords2$Site == x, "Long"])
  #tmplat <- rep(lat, lenlat)
  latvec <- append(latvec, lat)
  #tmplong <- rep(long, lenlong)
  longvec <- append(longvec, long)
}
site_coords_ind <- data.frame(all_site_coords2$Site, latvec, longvec)
# List of replacements needed: 30[16], 33[17], 35[18:19], 37[20], 44[34:35], 45[36:38], 48[44:47], 
# 63[76:79], 65[81:83], 77[116:119]
# 17, 25, 30, 33, 37, 44 32.0258333, -111.559117

site_coords_ind[5:6,2] <- 32.0258333 #replace the NAs for site 10
site_coords_ind[5:6,3] <- -111.559117 #replace the NAs for site 10
site_coords_ind[16,2] <- 35.10335 #replace the NAs for site 30
site_coords_ind[16,3] <- -119.4185 #replace the NAs for site 30
site_coords_ind[17,2] <- 34.93877 #replace the NAs for site 33
site_coords_ind[17,3] <- -118.9337#replace the NAs for site 33
site_coords_ind[18:19,2] <- 34.49958 #replace the NAs for site 35
site_coords_ind[18:19,3] <- -117.4106 #replace the NAs for site 35
site_coords_ind[20,2] <- 34.33623 #replace the NAs for site 37
site_coords_ind[20,3] <- -117.6037 #replace the NAs for site 37
site_coords_ind[34:35,2] <- 37.05157 #replace the NAs for site 44
site_coords_ind[34:35,3] <- -113.2733 #replace the NAs for site 44
site_coords_ind[36:38,2] <- 36.87802 #replace the NAs for site 45
site_coords_ind[36:38,3] <- -112.6480 #replace the NAs for site 45
site_coords_ind[44:47,2] <- 34.71380 #replace the NAs for site 48
site_coords_ind[44:47,3] <- -111.7793 #replace the NAs for site 48
site_coords_ind[76:79,2] <- 38.79937 #replace the NAs for site 63
site_coords_ind[76:79,3] <- -112.8255 #replace the NAs for site 63
site_coords_ind[81:83,2] <- 36.52380 #replace the NAs for site 65
site_coords_ind[81:83,3] <- -114.9415 #replace the NAs for site 65
site_coords_ind[116:119,2] <- 40.6104167 #replace the NAs for site 77
site_coords_ind[116:119,3] <- -119.7117 #replace the NAs for site 77

#our Q sites don't have coordinates so they will be removed with this line of code
site_coords <- aggregate( . ~ all_site_coords2.Site, data=site_coords_ind, FUN=mean)

var_coords <- aggregate( . ~ Pop_ID, data=all_var_coords, FUN=mean)
#ind_coords <- all_ind_coords[-c(67,71,210,213:215),]
#ind_coords2 <- all_ind_coords

##############################################################################
# PCA
##############################################################################

# data("H3N2")
# pop(H3N2) <- factor(H3N2$other$epid)
# View(H3N2@tab)
# da_pc <- dapc(H3N2, var.contrib=FALSE, scale=FALSE, n.pca=150, n.da=5)

gp_may.1 <- read.genepop("populations.snps.gen")
gp_may <- read.genepop("populations.snps.gen")
row.names(gp_may.1@tab)
row.names(gp_may.1@tab[-c(20,22,90,158,159,160,182),])
gp_may@tab <- gp_may.1@tab[-c(20,22,90,158,159,160,182),]
gp_may@ploidy <- gp_may.1@ploidy[-c(20,22,90,158,159,160,182)]
#out_gp <- genind_to_genepop(gp_M2, output = "galex_input.txt")
gp_may2 <- genind2genpop(gp_may)

ind_order <- as.data.frame(row.names(gp_may$tab))
ind_order$Sample_ID <-ind_order$`row.names(gp_may$tab)`

pop_info <- left_join(ind_order, all_meta_data)
pop_info[]
gp_may@pop <- as.factor(pop_info$Pop_ID)


gp_tab_i<- tab(gp_may, freq=TRUE, NA.method="mean")
gp_tab_b <- tab(gp_may, freq=F, NA.method="mean")
?tab()
#ri <- c("34E", "69I", "45F", "69G")
#gp_tab <- gp_tab_i[!row.names(gp_tab_i) %in% ri,]

pca.ind <- dudi.pca(gp_tab_i, center=TRUE, scale=FALSE, scannf = FALSE, nf=2)

grp <- find.clusters(gp_may, max.n.clust=50) #I did 150 PCs and 11 clusters
grp_ids <- as.numeric(as.vector(grp$grp))
grp_labs <- row.names(gp_may@tab)
dapc_clusters <- cbind(grp_ids, grp_labs)
write.csv(dapc_clusters, "~/Downloads/2023_lent/May_lent/dapc_K11_assignments.csv")

da_pc <- dapc(gp_may, grp$grp) # 120 PCs, 3 LDs
da_pc_var <- dapc(gp_may) # 120 PCs, 3 LDs
scatter(da_pc_var)
scatter(da_pc_var, scree.da=FALSE, scree.pca=F, bg="white", cstar=0,col=q_colors, solid=.7,cex=4,clab=0, leg=TRUE)
scatter(da_pc, scree.da=FALSE, scree.pca=F, bg="white", cstar=0,col=q_colors, solid=.7,cex=3,clab=0, leg=TRUE)



ind_names <- row.names(gp_may@tab)

id_table <- cbind(ind_names,grp_ids)

q_colors <- c("dodgerblue2", "#E31A1C", # red
                           "green4",
                           "#6A3D9A", # purple
                           "#FF7F00", # orange
                           "gold1",
                           "skyblue2", "#FB9A99", # lt pink
                           "palegreen2",
                           "#CAB2D6", # lt purple
                           "#FDBF6F", # lt orange
                           "gray70",
                           "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
                           "darkturquoise", "green1", "yellow4", "yellow3",
                           "darkorange4", "brown"
)


##############################################################################
# Visualize PCA
##############################################################################

# Create a data.frame containing individual coordinates
ind_coords = as.data.frame(pca.ind$li)

# Rename columns of dataframe
colnames(ind_coords) = c("Axis1","Axis2")

# Add a column containing individuals
ind_coords$Ind = row.names(gp_tab_i)

#ind_coords <- subset(ind_coords, !(ind_coords$Ind %in% c("Q-63", "63.2", "63.3")))
nrow(ind_coords)

# Add a column with the site IDs
sub_meta = subset(all_meta_data2, all_meta_data2$Sample_ID %in% ind_coords$Ind)
#sub_meta = sub_meta[-93,]
#subset(sub_meta, !(sub_meta$Sample_ID %in% c("Q-63", "63.2", "63.3")))[,3]
ind_coords$Site <- sub_meta[,3]
ind_coords$Variety = sub_meta[,2]

# Calculate centroid (average) position for each population
centroid = aggregate(cbind(Axis1, Axis2) ~ Site, data = ind_coords, FUN = mean)
centroid2 = aggregate(cbind(Axis1, Axis2) ~ Variety, data = ind_coords, FUN = mean)
centroid2[4,2:3] <- c(-8,-1)

# Analyse how much percent of genetic variance is explained by each axis
percent = pca.ind$eig/sum(pca.ind$eig)*100

# Add centroid coordinates to ind_coords dataframe
ind_coords = left_join(ind_coords, centroid, by = "Site", suffix = c("",".cen"))
ind_coords = left_join(ind_coords, centroid2, by = "Variety", suffix = c("",".cen"))


ind_coords$araneosus <- as.factor(ind_coords$Variety == "araneosus")
ind_coords$fremontii <- as.factor(ind_coords$Variety == "fremontii")
ind_coords$variablis <- as.factor(ind_coords$Variety == "variablis")
ind_coords$antonius <- as.factor(ind_coords$Variety == "antonius")

# Define colour palette
cols = rainbow(70)
cols2 = rainbow(18)

q_colors2 <- c("darkturquoise", "dodgerblue2", "#E31A1C", # red
                              "green4",
                              "#6A3D9A", # purple
                              "#FF7F00", # orange
                              "gold1", "green1",
                              "skyblue2", "#FB9A99", # lt pink
                              "palegreen2",
                              "#CAB2D6", "blue1",# lt purple
                              "gray70",
                              "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
                              "yellow4", "yellow3",
                              "brown")
                              


# Custom x and y labels
xlab = paste("Axis 1 (", format(round(percent[1], 1), nsmall=1)," %)", sep="")
ylab = paste("Axis 2 (", format(round(percent[2], 1), nsmall=1)," %)", sep="")

# Custom theme for ggplot2
ggtheme = theme(axis.text.y = element_text(colour="black", size=12),
                axis.text.x = element_text(colour="black", size=12),
                axis.title = element_text(colour="black", size=12),
                panel.border = element_rect(colour="black", fill=NA, size=1),
                panel.background = element_blank(),
                plot.title = element_text(hjust=0.5, size=15) 
)

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector2 = as.vector(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))))
col_vector <- c(col_vector2, "#E495A5", "#ABB065", "#39BEB1")

# Scatter plot axis 1 vs. 2
ggplot(data = ind_coords, aes(x = Axis1, y = Axis2))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  # spider segments
  #geom_segment(aes(xend = Axis1.cen, yend = Axis2.cen, colour = Site), show.legend = FALSE)+
  # points
  geom_point(aes(fill = Site), shape = 21, size = 3, show.legend = T)+
  # centroids
  #geom_label(data = centroid, aes(label = Site, fill = Site), size = 4, show.legend = FALSE)+
  # colouring
  scale_fill_manual(values = col_vector)+
  scale_colour_manual(values = col_vector)+
  # custom labels
  labs(x = xlab, y = ylab)+
  ggtitle("A. lentiginosus Site PCA")+
  # custom theme
  ggtheme

araneosus_colors <- c("gray","dodgerblue2")
fremontii_colors <- c("gray", "#FF7F00")
variablis_colors <- c("gray", "blue1")
antonius_colors <- c("gray","darkturquoise")
factorVar <- as.factor(ind_coords$Variety)

# Scatter plot axis 1 vs. 2
ggplot(data = ind_coords, aes(x = Axis1, y = Axis2, group=Variety, fill=Variety))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  # spider segments
  #geom_segment(aes(xend = Axis1.cen, yend = Axis2.cen, colour = Variety), show.legend = FALSE)+
  # points
  #stat_ellipse()+
  geom_point(aes(fill = Variety), shape = 21, size = 3, show.legend = T)+
  # centroids
  #geom_label(data = centroid2, aes(label = Variety, fill = Variety), size = 3, show.legend = FALSE)+
  #geom_text(data=subset(ind_coords, Variety == "antonius"), aes(label=Site))+
  # colouring
  scale_fill_manual(values = q_colors)+
  scale_colour_manual(values = q_colors)+
  # custom labels
  labs(x = xlab, y = ylab)+
  ggtitle("A. lentiginosus Variety PCA")+
  # custom theme
  ggtheme




##############################################################################
# Geo dists - by site
##############################################################################

# Convert long-lat to km distances
swap_site_coords <- cbind(site_coords[,3], site_coords[,2])

dist_mat_m <- distm(swap_site_coords, fun = distGeo)
dist_mat_km <- dist_mat_m %/% 1000

row.names(dist_mat_m) <- site_coords[,1]
colnames(dist_mat_m) <- site_coords[,1]

row.names(dist_mat_km) <- site_coords[,1]
colnames(dist_mat_km) <- site_coords[,1]

# remove rows and cols where there's only one individual per site
rem_inds <- c("1", "25", "30", "33", "37", "39", "42", "56", "63S", "64", "85", "Q", "V")
dist_mat_m_df <- as.data.frame(dist_mat_m)
dist_mat_km_df <- as.data.frame(dist_mat_km)

write.csv(dist_mat_km_df, file="km_btw_sites.csv")

#Syntax to drop columns using select()
mat_m_df <- select(dist_mat_m_df, -c("1", "39", "42", "56", "64", "85"))
mat_km_df <- select(dist_mat_km_df, -c("1", "39", "42", "56", "64", "85"))

dist_mat_df_m<-mat_m_df[!(row.names(mat_m_df) %in% rem_inds),]
dist_mat_df_km<-mat_km_df[!(row.names(mat_km_df) %in% rem_inds),]


row.names(gp_may@tab)
##############################################################################
# formatting
##############################################################################


#create new gp objects
gp_may_var <- gp_may
gp_may_var@pop <- as.factor(pop_info$Pop_ID)

#create new gp objects
gp_may_site <- gp_may
gp_may_site@pop <- as.factor(pop_info$Site)

#create new gp objects
gp_may_ind <- gp_may 
gp_may_ind@pop <- as.factor(row.names(gp_may@tab))


#convert gp to to gi
gi_may_var <- genind2genpop(gp_may_var)
gi_may_site <- genind2genpop(gp_may_site)
gi_may_ind <- genind2genpop(gp_may_ind)
class(gi_may_var)

#export to genalex format
out_ga1 <- genind2genalex(gp_may_var, filename = "galex_may_var_input.csv")
out_ga2 <- genind2genalex(gp_may_site, filename = "galex_may_site_input.csv")



# #IBD by var
gene_dist_var <- dist.genpop(gi_may_var, method=2)
# we have to remove, palmeri, sierrae, and micans
gdv_temp <- as.matrix(gene_dist_var)
gdv_temp <- gdv_temp[-c(11,14,17),-c(11,14,17)]
gene_dist_var <- as.dist(gdv_temp)
geo_dist_var <- dist(var_coords)
dim(as.matrix(geo_dist_var))
dim(as.matrix(gene_dist_var))
ibd_var <- mantel.randtest(gene_dist_var, geo_dist_var)
ibd_var


# #IBD by site
gene_dist_site <- dist.genpop(gi_may_site, method=2)
gds_temp <- as.matrix(gene_dist_site)
gds_temp <- gds_temp[-c(3,6,35,38,63:64),-c(3,6,35,38,63:64)] # remove inds to get data frames to match
gene_dist_site <- as.dist(gds_temp)
#remove 63S, 65S, Q and V
geo_dist_site <- dist(site_coords)
#missing values for sites: 17, 25, 30, 33, 37, 44 
dim(as.matrix(geo_dist_site))
dim(as.matrix(gene_dist_site))

ibd_site <- mantel.randtest(gene_dist_site, geo_dist_site)
ibd_site

# 
# #IBD by site and variety
ibd_s <- gl.ibd(Dgen = as.matrix(gene_dist_site), Dgeo= as.matrix(geo_dist_site), permutations=999)
ibd_v <- gl.ibd(Dgen = as.matrix(gene_dist_var), Dgeo= as.matrix(geo_dist_var), permutations=999)


##############################################################################
# Nei and WC FST
##############################################################################

#for this we need to use genind objects so make sure that all of our matrices are converted
# we are not computing fst for inds because you can't have a single individual representing a population


# FST  var according to Nei + Weir and Cockerham genetic distances
fst_var_WC <- pairwise.WCfst(gp_may_var)
fst_var_Nei <- pairwise.neifst(gp_may_var)

# FST site according to Nei's + Weir and Cockerham genetic distances
fst_site_WC <- pairwise.WCfst(gp_may_site)
fst_site_Nei <- pairwise.neifst(gp_may_site)


write.csv(fst_site_WC, file="~/Downloads/2023_lent/May_lent/site_fst_WC.csv")
write.csv(fst_site_Nei, file="~/Downloads/2023_lent/May_lent/site_fst_Nei.csv")

write.csv(fst_var_WC, file="~/Downloads/2023_lent/May_lent/var_fst_WC.csv")
write.csv(fst_var_Nei, file="~/Downloads/2023_lent/May_lent/var_fst_Nei.csv")


# FST for all sites with only one individuals will be 0  
# lets go ahead and remove those columns and rows from gen and geo matrices

#Check for Q
#inds_to_remove <- c("2A", "4D", "21C", "20E-2", "75F")

fstAdj <- function(x) {
  return(x/(1-x))
}

plot(hist(as.dist(fst_site_WC)))
plot(hist(as.dist(fst_site_Nei)))

fst_site_adj_WC <- apply(fst_site_WC, 2, fstAdj)
fst_site_adj_Nei <- apply(fst_site_Nei, 2, fstAdj)

fst_var_adj_WC <- apply(fst_var_WC, 2, fstAdj)
fst_var_adj_Nei <- apply(fst_var_Nei, 2, fstAdj)

plot(hist(as.dist(fst_site_adj_WC)))
plot(hist(as.dist(fst_site_adj_Nei)))
# fst_site_adj_WC <- fst_site_adj_WC[-70, -70]
# fst_site_adj_Nei <- fst_site_adj_Nei[-70, -70]

ordered_fst_WC_s2 <- fst_site_adj_WC[order(colnames(fst_site_adj_WC)), order(colnames(fst_site_adj_WC))]
ordered_fst_Nei_s2 <- fst_site_adj_Nei[order(colnames(fst_site_adj_Nei)), order(colnames(fst_site_adj_Nei))]


# ordering <- c(7,23, 41, 49, 56:58, 2:6, 10:14, 16:22, 24:32, 34:40, 42:48, 50:55)
# ordering2 <- c(7, 1:6, 10:14, 23, 16:22, 24:32, 41, 34:40, 49, 42:48, 56, 50:54, 57:58)

#ordered_fst_WC <- fst_site_adj_WC[ordering2,ordering2]
#ordered_fst_Nei <- fst_site_adj_Nei[ordering2,ordering2]


ordered_fst_WC_s <- ordered_fst_WC_s2[-c(1,3,6,8,9,11,12,16,18,28,35,36,38,60,63:64),-c(1,3,6,8,9,11,12,16,18,28,35,36,38,60,63:64)]
ordered_fst_Nei_s <- ordered_fst_Nei_s2[-c(1,3,6,8,9,11,12,16,18,28,35,36,38,60,63:64),-c(1,3,6,8,9,11,12,16,18,28,35,36,38,60,63:64)]


row.names(ordered_fst_WC_s) == row.names(dist_mat_df_km)
colnames(ordered_fst_WC_s) <- row.names(dist_mat_df_km)

row.names(ordered_fst_Nei_s) == row.names(dist_mat_df_km)
colnames(ordered_fst_Nei_s) <- row.names(dist_mat_df_km)

dim(fst_site_WC)
dim(dist_mat_km)

##############################################################################
# Hudson's FST
##############################################################################
library(KRIS)


# Hudson pairwise FST for site
site_idxs <- list()
num=1
for (sitepop in levels(gp_may_site@pop)){
  site_idx <- which(gp_may_site@pop==sitepop)
  site_idxs[[num]] <- c(site_idx)
  num = num+1
}

length(levels(gp_may_site@pop))

col_site <- list()
for (num in 1:length(levels(gp_may_site@pop))){
  siteV <- rep(NA, num)
  for (num2 in (num+1):63){
    pairwise_site <- fst.hudson(gp_may_site@tab, site_idxs[[num]], site_idxs[[num2]])
    siteV <- append(siteV, pairwise_site)
  }
  col_site[[num]] <- siteV
}
fst_site_hudson <- do.call("cbind", col_site)
row.names(fst_site_hudson) <- as.vector(levels(gp_may_site@pop))


# Hudson pairwise FST for var
var_idxs <- list()
num=1
for (varpop in levels(gp_may_var@pop)){
  var_idx <- which(gp_may_var@pop==varpop)
  var_idxs[[num]] <- c(var_idx)
  num = num+1
}

length(levels(gp_may_var@pop))

col_var <- list()
for (num in 1:length(levels(gp_may_var@pop))){
  varV <- rep(NA, num)
  for (num2 in (num+1):21){
    pairwise_var <- fst.hudson(gp_may_var@tab, var_idxs[[num]], var_idxs[[num2]])
    varV <- append(varV, pairwise_var)
  }
  col_var[[num]] <- varV
}
fst_var_hudson <- do.call("cbind", col_var)
row.names(fst_var_hudson) <- as.vector(levels(gp_may_var@pop))

write.csv(fst_var_hudson, file="var_fst_Hudson.csv")
write.csv(fst_site_hudson, file="site_fst_Hudson.csv")

# After this I formatted the Hudson fst csv file so that it matched the format of the other two csv files

hud_vec <- as.vector(fst_site_hudson) #convert matrix to vector
hud_vec<-hud_vec[!is.na(hud_vec)] # remove NAs

hud.mean <- mean(hud_vec) # This is the value we will use for FST in our PST plots
sample.n <- length(hud_vec)
sample.sd <- sd(hud_vec)
sample.se <- sample.sd/sqrt(sample.n)
alpha = 0.05
degrees.freedom = sample.n - 1
t.score = qt(p=alpha/2, df=degrees.freedom,lower.tail=F)
margin.error <- t.score * sample.se

lower.bound <- hud.mean - margin.error
upper.bound <- hud.mean + margin.error


##############################################################################
# FST adjustment
##############################################################################
rfst_site_Hud <- read.csv("site_fst_Hudson_formatted.csv")
snames <- rfst_site_Hud[,1]
rfst_site_Hud <- rfst_site_Hud[,-1]
row.names(rfst_site_Hud) <- snames
colnames(rfst_site_Hud) <- snames

plot(hist(as.dist(rfst_site_Hud)))

fst_site_adj_Hud <- apply(rfst_site_Hud, 2, fstAdj)
dim(fst_site_adj_Hud)
dim(dist_mat_df_km)
row.names(fst_site_adj_Hud) 
row.names(dist_mat_df_km)

#ordered_fst_Hudson <- fst_site_adj_Hud[ordering2,ordering2]
#ordered_fst_Hud <- rfst_site_Hud[ordering2,ordering2]
ordered_fst_Hudson <- rfst_site_Hud[-c(1,3,6,8,9,11,12,16,18,28,35,36,38,60,63:64),-c(1,3,6,8,9,11,12,16,18,28,35,36,38,60,63:64)]
dim(ordered_fst_Hudson)
row.names(ordered_fst_Hudson) == row.names(dist_mat_df_km)
row.names(ordered_fst_WC_s) == row.names(dist_mat_df_km)


##############################################################################
# IBD
##############################################################################
ibd_site_WC_m <- mantel.randtest(as.dist(ordered_fst_WC_s), as.dist(dist_mat_df_m))
ibd_site_Nei_m <- mantel.randtest(as.dist(ordered_fst_Nei_s), as.dist(dist_mat_df_m))
ibd_site_Hudson_m <- mantel.randtest(as.dist(ordered_fst_Hudson), as.dist(dist_mat_df_m))

ibd_site_WC_km <- mantel.randtest(as.dist(ordered_fst_WC), as.dist(dist_mat_df_km))
ibd_site_Nei_km <- mantel.randtest(as.dist(ordered_fst_Nei), as.dist(dist_mat_df_km))
ibd_site_Hudson_km <- mantel.randtest(as.dist(ordered_fst_Hudson), as.dist(dist_mat_df_km))

ibd_SITE_m <- gl.ibd(Dgen = as.dist(ordered_fst_WC_s), Dgeo= as.dist(dist_mat_df_m), permutations=999, plot.out = T)
ibd_SITE_km <- gl.ibd(Dgen = as.dist(ordered_fst_WC_s), Dgeo= as.dist(dist_mat_df_km), permutations=999, plot.out = T)
Hud_ibd_SITE_km <- gl.ibd(Dgen = as.dist(ordered_fst_Hudson), Dgeo= as.dist(dist_mat_df_km), permutations=999, plot.out = T)
Nei_ibd_SITE_km <- gl.ibd(Dgen = as.dist(ordered_fst_Nei_s), Dgeo= as.dist(dist_mat_df_km), permutations=999, plot.out = T)

##############################################################################
# AMOVA
##############################################################################

###  Tutorial practice
# data(Aeut)
# strata(Aeut) <- other(Aeut)$population_hierarchy[-1]
# agc <- as.genclone(Aeut)
# agc
# amova.result <- poppr.amova(agc, ~Pop/Subpop)

gp_vs <- gp_may_var

strata_df <- as.data.frame(row.names(gp_vs$tab))
strata_df$Variety <- gp_vs$pop
strata_df$Site <-gp_may_site$pop
# strata_df[32,2] <- "20"
# strata_df[98:100,2] <- "44"
# strata_df[177,2] <- "Q"

strata(gp_vs) <- strata_df[-1]
gc_vs <- as.genclone(gp_vs)

amova_var_site <- poppr.amova(gc_vs, ~Variety/Site)
amova_var <- poppr.amova(gc_vs, ~Variety)
amova_site <- poppr.amova(gc_vs, ~Site)

amova_vs_signif <- randtest(amova_var_site, nrepet = 999)
amova_v_signif <- randtest(amova_var, nrepet = 999)
amova_s_signif <- randtest(amova_site, nrepet = 999)

plot(amova_vs_signif)
plot(amova_v_signif)
plot(amova_s_signif)

##############################################################################
# Heterozygosity estimates
##############################################################################

gp_may_all <- gp_may
gp_may_all$pop <- as.factor(rep(1, nrow(gp_may@tab)))

# adegenet heterozygosity from a genpop object
hs_var <- Hs(gp_may_var)
hs_site <- Hs(gp_may_site)
hs_ind <- Hs(gp_may_ind)
hs_all <- Hs(gp_may_all)

# adegenet hardy weinberg equilibrium tests
hw_var <- hw.test(gp_may_var)

# adegenet inbreeding
## estimate inbreeding - return proba density functions
Fdens <- inbreeding(gp_may_var, res.type="function")
## estimate inbreeding - return maximum likelihood estimates
Fest <- inbreeding(gp_may_var, res.type = "estimate")
mostInbred <- which.max(Fest)
plot(Fdens[[mostInbred]], ylab = "Density", xlab = "F",
     main = paste("Probability density of F values\nfor", names(mostInbred)))
abline(v = Fest[mostInbred], col = "red", lty = 2)
legend("topright", legend = "MLE", col = "red", lty = 2)

##############################################################################
# NJT - look like shit
##############################################################################

gene_dist_ind <- dist.genpop(gi_may_ind, method=3)

#Individuals
Ind_tre <- nj(as.dist(gene_dist_ind))
class(Ind_tre)
length(Ind_tre)
Ind_tre <- ladderize(Ind_tre)
plot(Ind_tre, cex=.6)
title("Neighbor Joining tree: Individual")

#Site
Site_tre <- nj(gene_dist_site)
class(Site_tre)
length(Site_tre)
Site_tre <- ladderize(Site_tre)
plot(Site_tre, cex=.6)
title("Neighbor Joining Tree: Site")

#Variety
Var_tre <- nj(gene_dist_var)
class(Var_tre)
length(Var_tre)
Var_tre <- ladderize(Var_tre)
plot(Var_tre, cex=.6)
title("Neighbor Joining Tree: Variety")

