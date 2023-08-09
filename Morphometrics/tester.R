#Install all the required packages for the analysis
#install.packages("devtools")
library(devtools)
#devtools::install_github('geomorphR/geomorph',ref="Stable")
library(geomorph)
library(data.table)
library(dplyr)
library(RRPP)

#Set my working directory
setwd("/Users/QuinnThomas/Documents/Desktop_Items/landmarkFiles/CSV_files/")
wd = "/Users/QuinnThomas/Documents/Desktop_Items/landmarkFiles/CSV_files/"
list_csv <- list.files(wd, pattern = "csv")
#Read all csv files listed from my working directory and format it to include relevant information 
DT = do.call(rbind, lapply(list_csv, fread))
#DT = rbindlist(lapply(list_csv, fread), use.names=FALSE)
masterList <- read.csv("/Users/QuinnThomas/Documents/Desktop_Items/landmarkFiles/lentiginosusMaster.csv")
formatMaster <- masterList[c(1:13)]
class(formatMaster)
#pull plant and variety information
idVariety2 <- masterList[c(2,3)]
nrow(idVariety)
hmmm <- read.csv("~/Documents/Desktop_Items/Research/Var_K_map.csv")
idVariety <- left_join( idVariety2, hmmm)

#Remove missing data
newdata <- na.omit(DT)
numrows = nrow(newdata)

#Format data into a 3d array
Sliced <- aperm(`dim<-`(t(newdata), c(6, 30, (numrows/30))), c(2,1,3))


#Save one individual fruit so we can set our sliding coordinates
# FirstSample <- array(read.csv(file.choose(),
#                nrow=30,
#                header=TRUE))
# XYlandmarks<- FirstSample[,5:6]

dev.off()
#create a matrix to indicate which points slide with respect to their anchors
#aerial_curvers <- define.sliders(XYlandmarks, nsliders=26, write.file = T)
#aerial_curvers2 <- define.sliders(XYlandmarks, nsliders=22, write.file = T)

ac <- read.csv("curveslide.txt", header = T)
ac2 <- read.csv("curveslide2.txt", header = T)

aerial_curvers <- ac
aerial_curvers2 <- ac2

#Verify that our data does not contain any missing values
no_Missing <- na.omit(newdata)
any(is.na(no_Missing))
#pull only landmark information
lentiginosusLand <- as.matrix(no_Missing[,c(5:6)])
#reformat to 3d array
A <- arrayspecs(lentiginosusLand, 30, 2)
any(is.na(A))
#if any missing tell me which
which(is.na(A))
#check for infinite?
which(!is.finite(A))
all(sapply(A, is.finite))


#geomorph procrustes of all fruits
lent.gpa<- gpagen(
  A = A,
  curves = aerial_curvers,
  ProcD = FALSE)

#geomorph procrustes of all fruits
lent.gpaC2<- gpagen(
  A = A,
  curves = aerial_curvers2,
  ProcD = FALSE)

#plot(procrustes)
plot(lent.gpa)
plot(lent.gpaC2)

#Principal component analysis
myPCA <- gm.prcomp(lent.gpaC2$coords)

#plantID for all fruits
q_labels <- as.factor(no_Missing[["sample_ID"]])
plot(myPCA, col=rainbow(39))
sites.col <- levels(as.factor(sites3))
legend("topright", legend=sites.col, pch=19, title='Plant', col=rainbow(39))



#plantID for every individual landmark
vec <- as.vector(no_Missing[[1]])
#plantID for every individual fruit
labels <- vec[seq(1, length(vec), 30)]
#confirm it is what you expect
class(labels)
levels(as.factor(labels))
#length(varLabels)
labels <- recode(labels, "3E" = "19-3E", "3J" = "19-3J", "4E" = "19-4E", "4H" = "19-4H", 
                 "4I" = "19-4I", "4J" = "19-4J", "4K" = "19-4K", "4L" = "19-4L", "5H" = "19-5H", 
                 "5J" = "19-5J", "5L" = "19-5L", "6C" = "19-6C", "6D" = "19-6D", "6F" = "19-6F", 
                 "6G" = "19-6G", "6H" = "19-6H", "6I" = "19-6I")

dddf <- as.data.frame(labels)
mergeVar2 <- merge(dddf, idVariety, by.x = "labels", by.y = "Sample", no.dups=FALSE,  all.x=TRUE, all.y=FALSE)
which(is.na(mergeVar2[1]))
#REMOVE: rows 2E
mergeVar2 <- mergeVar2[-c(621:630),]

#formatting
rowNum = nrow(A)
my_arr <- array(A, dim = c(30,2,length(labels)))
length(my_arr)
length(labels)
#subset by plant
subset.coords <- coords.subset(A=my_arr, group=labels)
#Mean shape for a given plant
ind_means <- as.array(lapply(subset.coords, mshape))

#confirm tht you have what you expect
nrow(ind_means)
nrow(varLabels)

#format matrix for all mean shapes of plants
combined <- do.call(rbind,ind_means)
B <- arrayspecs(combined, 30,2)

#Get variety information for each data point
sitePlants <- names(ind_means)
#write.csv(sitePlants, file="~/Documents/Desktop_Items/Research/sitePlants.csv")
labelFrame <- as.data.frame(sitePlants)
mergeVar <- merge(labelFrame, idVariety, by.x = "sitePlants", by.y = "Sample", no.dups=TRUE,  all.x=TRUE, )
#mergeVar <- mergeVar[46,]
varLabels <- as.vector(mergeVar[["Variety"]])
#Remove palans 19-43L duplicate
varLabels <- varLabels[-c(132)]

K_labels <- as.vector(mergeVar$K.pop)
K_labels <- K_labels[-c(132)]
write.csv(varLabels, file="~/Documents/Desktop_Items/Research/varLabels.csv")


#remove maricopae and wilsonii bc they are skewing our dataset
marIndex <- which(varLabels == 'maricopae')
wilIndex <- which(varLabels == 'wilsonii')
varIndex <- c(marIndex,wilIndex)
sitePlants2 <- sitePlants[-varIndex]
varLabels2 <- varLabels[-varIndex]

#format
sites <-gsub('.{1}$', '', sitePlants)
sites2 <-gsub('.{1}$', '', sitePlants2)
sites3 <-gsub('.{1}$', '', labels)

levels(as.factor(sites))

#Subset data so that maricopae and wilsonii are excluded.
subset.mw <- coords.subset(A=B, group=varLabels)
subset.mw <- subset.mw[-c(6,11)]

#formatting data without maricopae and wilsonii.
combinedMW = matrix(, ncol=2)
for(val in subset.mw){
  #r=val
  y <- aperm(val, c(1, 3, 2))
  dim(y) <- c(prod(dim(val)[-2]), dim(val)[2])
  combinedMW <- rbind(combinedMW, y)
  #combinedMW <- rbind(combinedMW, subset.mw[[val]])
}
combinedMW<- combinedMW[-1,]
C <- arrayspecs(combinedMW, 30,2)


#check to make sure everything is as we expect
length(labelFrame)
length(sites)
length(varLabels)
length(K_labels)
dim(B)
    

#Procrustes of plants
lent2.gpa<- gpagen(
  A = B,
  curves = aerial_curvers,
  ProcD = FALSE)
lent2.1gpa<- gpagen(
  A = B,
  curves = aerial_curvers2,
  ProcD = FALSE)

#Procrustes of plants without maricopae and wilsonii
lent3.gpa<- gpagen(
  A = C,
  curves = aerial_curvers2,
  ProcD = FALSE)


subset.var <- coords.subset(A=lent2.gpa$coords, group=varLabels)
subset.var2 <- coords.subset(A=lent.gpa$coords, group=as.vector(mergeVar2[,2]))

plot(lent2.gpa)
plot(lent2.1gpa)
plot(lent3.gpa)
pca2 <- gm.prcomp(lent2.gpa$coords)
pca2.1 <- gm.prcomp(lent2.1gpa$coords)
summary(pca2)

#export Principal components so we can use in Pstat program
write.csv(pca2$x, "/Users/QuinnThomas/Documents/Desktop_Items/Research/PCvals.csv", row.names = FALSE)


#legend and plot design
col.rainbow <- rainbow(14)
col.11 <- col.rainbow[-c(5, 9)]
mycols <- colors()[c(153, 31, 8, 27, 47, 91, 34, 41, 139, 143, 642, 117, 220, 590)]

q_colors <- c("dodgerblue2", "#E31A1C", # red
                           "green4",
                           "#6A3D9A", # purple
                           "#FF7F00", # orange
                           "gold1",
                           "skyblue2", "#FB9A99", # lt pink
                           "palegreen2",
                           "#CAB2D6", "blue1", # lt purple
                           "gray70",
                           "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
                           "darkturquoise", "green1", "yellow4", "yellow3",
                           "darkorange4", "brown"
)

dev.off()

#plotting PCA for individual plants
facVar <- as.factor(varLabels)
facK <- as.factor(K_labels)
levels(facVar)
plot(pca2, col=q_colors[facK], pch=19)
legend("topright", legend=levels(as.factor(varLabels)), pch=19, title='Varieties', col=q_colors)
text(pca2$x, labels=sites)

plot(pca2.1, col=q_colors[facVar], pch=19)
legend("topright", legend=levels(facVar), pch=19, title='Varieties', col=q_colors)
text(pca2.1$x, labels=sites)

#blank plot
plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
legend("topleft", legend=levels(facVar), pch=19, title='Varieties', col=q_colors)

#Here are the min and max shapes for each principal component
plot(pca2$shapes$shapes.comp1$min)
plot(pca2$shapes$shapes.comp1$max)
plot(pca2$shapes$shapes.comp2$min)
plot(pca2$shapes$shapes.comp2$max)

plot(pca2.1$shapes$shapes.comp1$min)
plot(pca2.1$shapes$shapes.comp1$max)
plot(pca2.1$shapes$shapes.comp2$min)
plot(pca2.1$shapes$shapes.comp2$max)
#blank plot 
plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')

#plotting PCA for individual plants (w/o maricopae and wilsonii)
plot(lent3.gpa)
pca3 <- gm.prcomp(lent3.gpa$coords)
length(sites2)
length(varLabels2)
#PC1 and PC2
q_colors_mw <- q_colors[-c(7,13)]
plot(pca3, col=q_colors_mw[as.factor(varLabels2)], pch=19)
legend("topright", legend=levels(as.factor(varLabels2)), pch=19, title='Varieties', col=q_colors_mw)
text(pca3$x, labels=sites2)
#PC2 and PC3
plot(x=pca3$x[,2], y=pca3$x[,3], col=q_colors[as.factor(varLabels2)], pch=19)
#plot(pca3, col=col.10[as.factor(sites2)], pch=19)
text(pca3$x, labels=sites2)
legend("topright", legend=levels(as.factor(varLabels2)), pch=19, title='Varieties', col=q_colors)
#legend("topright", legend=levels(as.factor(sites2)), pch=c(1), title='Sites', col=unique(as.factor(sites2)))



###############################################################################

'''
The function quantifies the relative amount of shape variation attributable to one or more factors in
a linear model and estimates the probability of this variation ("significance") for a null model, via
distributions generated from resampling permutations. Data input is specified by a formula (e.g.,
y~X), where ’y’ specifies the response variables (Procrustes shape variables), and ’X’ contains one
or more independent variables (discrete or continuous)
'''

gdf <- geomorph.data.frame(lent2.gpa,
  site = sites,
  variety = varLabels) # geomorph data frame
gdf.mat <- as.matrix(gdf)

gdf2 <- geomorph.data.frame(lent3.gpa,
        site = sites2, variety = varLabels2) # geomorph data frame
gdf3 <- geomorph.data.frame(lent.gpa, site = sites3,
        variety = as.vector(mergeVar2[,2]), individual = labels) # geomorph data frame
gdf4 <- geomorph.data.frame(lent.gpaC2, site = sites3,
                            variety = as.vector(mergeVar2[,2]), individual = labels) # geomorph data frame

length(sites3)
length(varLabels)


gmat <- two.d.array(lent2.gpa$coords)
gmat2 <- two.d.array(lent.gpa$coords)
rrpp.df<-rrpp.data.frame(shape=gmat,CS=lent2.gpa$Csize,variety=varLabels, site=sites)
rrpp.df2<-rrpp.data.frame(shape=gmat2,CS=lent.gpa$Csize,variety=as.vector(mergeVar2[,2]), site=sites3, individual=labels)
rrpp.df3<-rrpp.data.frame(shape=gmat2,CS=lent.gpaC2$Csize,variety=as.vector(mergeVar2[,2]), site=sites3, individual=labels)


#Nested ANOVA
#rr.nest<-lm.rrpp(shape~variety/site,data=rrpp.df,iter=999,print.progress=F)
#anova(rr.nest)

#Two way Nested ANOVA... I think this might be redundant?
rr.nest2<-lm.rrpp(shape~variety/site+site,data=rrpp.df,iter=999,print.progress=F)
anova(rr.nest2)
levels(as.factor(rrpp.df$site))

#rr.nest3.1<-lm.rrpp(shape~CS * variety/site,data=rrpp.df,iter=999,print.progress=F)
#anova(rr.nest3.1)
#rr.nest3.2<-lm.rrpp(shape~(variety/site),data=rrpp.df,iter=999,print.progress=F)
#anova(rr.nest3.2)

#install.packages("olsrr")
library(olsrr)

rr.nest4 <- lm.rrpp(shape~(variety/site)+site, data=rrpp.df2,iter=999,print.progress=F)
anova(rr.nest4)
plot(rr.nest4)
mdis<-morphol.disparity(f1 = rr.nest4, groups = ~ variety, data = rrpp.df2,
                         iter = 999, print.progress = FALSE)


coef(rr.nest4, test = TRUE)
rr.nest5 <- lm.rrpp(shape~(site/individual)+individual, data=rrpp.df2,iter=999,print.progress=F)
anova(rr.nest5)
coef(rr.nest5, test = TRUE)

#I ran this for almost 48 hours and it never made it past the first step
#rr.ind <- lm.rrpp(shape~(variety/site/individual), data=rrpp.df2,iter=999,print.progress=T)
#anova(rr.ind)
#coef(rr.ind, test = TRUE)


PWT <- pairwise(fit, groups = interaction(Pupfish$Sex, Pupfish$Pop))
summary(PWT, confidence = 0.95)

#Here are my ANOVAs! Some information for myself
#Use a one-way ANOVA when you have collected data about one categorical independent variable and one quantitative dependent variable. The independent variable should have at least three levels (i.e. at least three different groups or categories).
#A two-way ANOVA is used to estimate how the mean of a quantitative variable changes according to the levels of two categorical variables. Use a two-way ANOVA when you want to know how two independent variables, in combination, affect a dependent variable.
#To test whether two variables have an interaction effect in ANOVA, simply use an asterisk instead of a plus-sign in the model:
#The p-value is low (p < 0.001), so it appears that x real impact on y.


#One way ANOVA on entire dataset for variety and centroid size
fit1.var <- procD.lm(coords ~ variety,
                 data = gdf, iter = 999, turbo = TRUE,
                 RRPP = TRUE, print.progress = FALSE) # randomize residuals

fit1.csize <- procD.lm(coords ~ Csize,
                 data = gdf, iter = 999, turbo = TRUE,
                 RRPP = TRUE, print.progress = FALSE) # randomize residuals

#One way ANOVA on partial dataset for variety and centroid size
fit4.var <- procD.lm(coords ~ variety,
                 data = gdf2, iter = 999, turbo = TRUE,
                 RRPP = TRUE, print.progress = FALSE) # randomize residuals
fit4.csize <- procD.lm(coords ~ Csize,
                 data = gdf2, iter = 999, turbo = TRUE,
                 RRPP = TRUE, print.progress = FALSE) # randomize residuals

PWT <- pairwise(fit1.var, groups = fit1.var$data$variety)
summary(PWT)
?pairwise

#Two way ANOVA on where variety INTERACTS with site (both= x) and coords =y
fit2 <- procD.lm(coords ~ variety+log(Csize),
                 data = gdf, iter = 999, turbo = TRUE,
                 RRPP = TRUE, print.progress = FALSE) # randomize raw values
fit3 <- procD.lm(coords ~ site * Csize,
                 data = gdf, iter = 999, turbo = TRUE,
                 RRPP = TRUE, print.progress = FALSE) # randomize raw values
fitAllom <- procD.lm(coords ~ log(Csize), data=gdf, print.progress = FALSE)
fitAllom2 <- procD.lm(coords ~ variety*log(Csize), data=gdf, iter=999, 
                     turbo = TRUE,RRPP = TRUE, print.progress = FALSE)


summary(fit1.var)
summary(fit2)
summary(fitAllom)
plot(fit1)


###############################################################################

#T-tests

install.packages("dispRity")
library("dispRity")

#example code
geomorph_df <- geomorph.data.frame(procrustes, species = plethodon$species)
geomorph.ordination(geomorph_df)

## Calculating disparity from dispRity or geomorph::morphol.disparity
geomorph_disparity <- geomorph::morphol.disparity(coords ~ 1,
                                                  groups= ~ species, data = geomorph_df)
dispRity_disparity <- dispRity(geomorph.ordination(geomorph_df),
                               metric = function(X) return(sum(X^2)/nrow(X)))

## Extracting the raw disparity values
geomorph_val <- round(as.numeric(geomorph_disparity$Procrustes.var), 15)
dispRity_val <- as.vector(summary(dispRity_disparity, digits = 15)$obs)
## Comparing the values (to the 15th decimal!)
geomorph_val == dispRity_val # all TRUE


## Measuring disparity as a distribution
disparity_var <- dispRity(bootstrapped_data, metric = variances)
## Differences between the concatenated bootstrapped values of the subsets
test.dispRity(disparity_var, test = t.test, comparisons = "pairwise",
              concatenate = TRUE, correction = "bonferroni")


###############################################################################


# Examples use geometric morphometric data on pupfishes
# See the package, geomorph, for details about obtaining such data
# Body Shape Analysis (Multivariate) --------------
data("Pupfish")
# Note:
dim(Pupfish$coords) # highly multivariate!
Pupfish$logSize <- log(Pupfish$CS)
# Note: one should use all dimensions of the data but with this
# example, there are many. Thus, only three principal components
# will be used for demonstration purposes.
Pupfish$Y <- ordinate(Pupfish$coords)$x[, 1:3]
## Pairwise comparisons of LS means
# Note: one should increase RRPP iterations but a
# smaller number is used here for demonstration
# efficiency. Generally, iter = 999 will take less
# than 1s for these examples with a modern computer.
fit1 <- lm.rrpp(Y ~ logSize + Sex * Pop, SS.type = "I",
                data = Pupfish, print.progress = FALSE, iter = 999)
summary(fit1, formula = FALSE)
anova(fit1)
pup.group <- interaction(Pupfish$Sex, Pupfish$Pop)
pup.group
PW1 <- pairwise(fit1, groups = pup.group)
PW1
# distances between means
summary(PW1, confidence = 0.95, test.type = "dist")
summary(PW1, confidence = 0.95, test.type = "dist", stat.table = FALSE)
# absolute difference between mean vector lengths
summary(PW1, confidence = 0.95, test.type = "DL")
# correlation between mean vectors (angles in degrees)
summary(PW1, confidence = 0.95, test.type = "VC",
        angle.type = "deg")
# Can also compare the dispersion around means
summary(PW1, confidence = 0.95, test.type = "var")
## Pairwise comparisons of slopes
fit2 <- lm.rrpp(Y ~ logSize * Sex * Pop, SS.type = "I",
                data = Pupfish, print.progress = FALSE, iter = 199)
summary(fit2, formula = FALSE)
anova(fit1, fit2)
# Using a null fit that excludes all factor-covariate
# interactions, not just the last one
PW2 <- pairwise(fit2, fit.null = fit1, groups = pup.group,
                covariate = Pupfish$logSize, print.progress = FALSE)
PW2
# distances between slope vectors (end-points)
summary(PW2, confidence = 0.95, test.type = "dist")
summary(PW2, confidence = 0.95, test.type = "dist", stat.table = FALSE)
# absolute difference between slope vector lengths
summary(PW2, confidence = 0.95, test.type = "DL")
# correlation between slope vectors (and angles)
summary(PW2, confidence = 0.95, test.type = "VC",
        angle.type = "deg")
# Can also compare the dispersion around group slopes
summary(PW2, confidence = 0.95, test.type = "var")


#MY DATA
twoD.lent <- two.d.array(lent2.gpa$coords)
dim(twoD.lent)

# Note: one should use all dimensions of the data but with this
# example, there are many. Thus, only three principal components
# will be used for demonstration purposes.

# Here we look at all varieties and all sites

gdf$Y <- ordinate(twoD.lent)$x[, 1:3]

#myDF<-rrpp.data.frame(shape=gmat,CS=lent2.gpa$Csize,variety=varLabels, site=sites)

fit1 <- lm.rrpp(Y ~(variety/site), SS.type = "I",
               data = gdf, print.progress = FALSE, iter = 999)
fit2 <- lm.rrpp(Y ~(site), SS.type = "I",
                data = gdf, print.progress = FALSE, iter = 999)
summary(fit1)
anova(fit1)
coef(fit1, test = TRUE)
anova(fit2)
coef(fit2, test = TRUE)

lent.group <- interaction(gdf$site, gdf$variety)

PW1 <- pairwise(fit1, groups = gdf$variety)
PW1.2 <- pairwise(fit2, groups = gdf$site)
pairwiseSite1<-summary(PW1, confidence = 0.95, test.type = "dist")
pairwiseSite1.2<-summary(PW1.2, confidence = 0.95, test.type = "dist")
View(pairwiseSite1[["pairwise.tables"]][["P"]])
View(pairwiseSite1.2[["pairwise.tables"]][["P"]])


CV1 <- looCV(fit1)
summary(CV1)
group <- interaction(Pupfish$Pop, Pupfish$Sex)
plot(CV1, flip = 1, pch = 19, col = as.factor(gdf$variety))


fit5 <- lm.rrpp(coords ~ Pop*Sex, data = Pupfish, iter = 0)
CV5 <- looCV(fit5)
summary(CV5)
group <- interaction(Pupfish$Pop, Pupfish$Sex)
plot(CV5, flip = 1, pch = 19, col = group)

n <- NROW(Pupfish$coords)
p <- NCOL(Pupfish$coords)
set.seed(1001)
Yr <- matrix(rnorm(n * p), n, p) # random noise
fit6 <-lm.rrpp(Yr ~ Pop*Sex, data = Pupfish, iter = 0)
CV6 <- looCV(fit6)
plot(CV6, pch = 19, col = group)

PWY <- pairwise(fit2, groups = lent.group)
pairwiseSite<-summary(PWY, confidence = 0.95, test.type = "dist")
View(pairwiseSite[["x"]][["means.dist"]][["obs"]])

summary(PW1, confidence = 0.95, test.type = "dist")
summary(PW1, confidence = 0.95, test.type = "dist", stat.table = FALSE)
# absolute difference between mean vector lengths
summary(PW1, confidence = 0.95, test.type = "DL")
# correlation between mean vectors (angles in degrees)
summary(PW1, confidence = 0.95, test.type = "VC",
        angle.type = "deg")

PW2 <- pairwise(fit1, groups = gdf$site)
summary(PW2, confidence = 0.95, test.type = "dist")
summary(PW2, confidence = 0.95, test.type = "dist", stat.table = FALSE)
# absolute difference between mean vector lengths
summary(PW2, confidence = 0.95, test.type = "DL")
# correlation between mean vectors (angles in degrees)
summary(PW2, confidence = 0.95, test.type = "VC",
        angle.type = "deg")

PW3 <- pairwise(fit1, groups = gdf$site)
summary(PW3, confidence = 0.95, test.type = "dist")
summary(PW3, confidence = 0.95, test.type = "dist", stat.table = FALSE)
# absolute difference between mean vector lengths
summary(PW3, confidence = 0.95, test.type = "DL")
# correlation between mean vectors (angles in degrees)
summary(PW3, confidence = 0.95, test.type = "VC",
        angle.type = "deg")

# Let's break down by variety

borreganusLabels <- c()
fremontiiLabels <- c()
kennedyiLabels <- c()
palansLabels <- c()
salinusLabels <- c()
vitreusLabels <- c()
wilsoniiLabels <- c()
yuccanusLabels <- c()

for(z in rep(1:nrow(mergeVar))){
  if(mergeVar[z,2] == "borreganus"){
    borreganusLabels <- c(borreganusLabels,mergeVar[z,1])
  }
  if(mergeVar[z,2] == "fremontii"){
    fremontiiLabels <- c(fremontiiLabels,mergeVar[z,1])
  }
  if(mergeVar[z,2] == "kennedyi"){
    kennedyiLabels <- c(kennedyiLabels,mergeVar[z,1])
  }
  if(mergeVar[z,2] == "palans"){
    palansLabels <- c(palansLabels,mergeVar[z,1])
  }
  if(mergeVar[z,2] == "salinus"){
    salinusLabels <- c(salinusLabels,mergeVar[z,1])
  }
  if(mergeVar[z,2] == "vitreus"){
    vitreusLabels <- c(vitreusLabels,mergeVar[z,1])
  }
  if(mergeVar[z,2] == "wilsonii"){
    wilsoniiLabels <- c(wilsoniiLabels,mergeVar[z,1])
  }
  if(mergeVar[z,2] == "yuccanus"){
    yuccanusLabels <- c(yuccanusLabels,mergeVar[z,1])
  }
}


borreganus <- geomorph.data.frame(subset.var$borreganus,
                    plant = borreganusLabels,
                    site = gsub('.{1}$', '', borreganusLabels)) # geomorph data frame
twoD.bor <- two.d.array(subset.var$borreganus)
borreganus$Y <- ordinate(twoD.bor)$x[, 1:4]

fit.bor <- lm.rrpp(Y ~ site * plant, SS.type = "I",
                data = borreganus, print.progress = FALSE, iter = 999)

PWsite.bor <- pairwise(fit.bor, groups = borreganus$site)
summary(PWsite.bor, confidence = 0.95, test.type = "dist")
summary(PWsite.bor, confidence = 0.95, test.type = "DL")
summary(PWsite.bor, confidence = 0.95, test.type = "VC", angle.type = "deg")

bor.group <- interaction(borreganus$plant, borreganus$site)
PWplant.bor <- pairwise(fit.bor, groups = bor.group)
summary(PWplant.bor, confidence = 0.95, test.type = "dist")
summary(PWplant.bor, confidence = 0.95, test.type = "DL")
summary(PWplant.bor, confidence = 0.95, test.type = "VC", angle.type = "deg")

fremontii <- geomorph.data.frame(subset.var$fremontii,
                                  plant = fremontiiLabels,
                                  site = gsub('.{1}$', '', fremontiiLabels)) # geomorph data frame
twoD.fre <- two.d.array(subset.var$fremontii)
fremontii$Y <- ordinate(twoD.fre)$x[, 1:4]

fit.fre <- lm.rrpp(Y ~ site * plant, SS.type = "I",
                   data = fremontii, print.progress = FALSE, iter = 999)

PWsite.fre <- pairwise(fit.fre, groups = fremontii$site)
summary(PWsite.fre, confidence = 0.95, test.type = "dist")
summary(PWsite.fre, confidence = 0.95, test.type = "DL")
summary(PWsite.fre, confidence = 0.95, test.type = "VC", angle.type = "deg")

fre.group <- interaction(fremontii$plant, fremontii$site)
PWplant.fre <- pairwise(fit.fre, groups = fre.group)
summary(PWplant.fre, confidence = 0.95, test.type = "dist")
summary(PWplant.fre, confidence = 0.95, test.type = "DL")
summary(PWplant.fre, confidence = 0.95, test.type = "VC", angle.type = "deg")

kennedyi <- geomorph.data.frame(subset.var$kennedyi,
                                 plant = kennedyiLabels,
                                 site = gsub('.{1}$', '', kennedyiLabels)) # geomorph data frame
twoD.ken <- two.d.array(subset.var$kennedyi)
kennedyi$Y <- ordinate(twoD.ken)$x[, 1:4]

fit.ken <- lm.rrpp(Y ~ site * plant, SS.type = "I",
                   data = kennedyi, print.progress = FALSE, iter = 999)

PWsite.ken <- pairwise(fit.ken, groups = kennedyi$site)
summary(PWsite.ken, confidence = 0.95, test.type = "dist")
summary(PWsite.ken, confidence = 0.95, test.type = "DL")
summary(PWsite.ken, confidence = 0.95, test.type = "VC", angle.type = "deg")

palans <- geomorph.data.frame(subset.var$palans,
                                plant = palansLabels,
                                site = gsub('.{1}$', '', palansLabels)) # geomorph data frame
twoD.pal <- two.d.array(subset.var$palans)
palans$Y <- ordinate(twoD.pal)$x[, 1:4]

fit.pal <- lm.rrpp(Y ~ site * plant, SS.type = "I",
                   data = palans, print.progress = FALSE, iter = 999)

PWsite.pal <- pairwise(fit.pal, groups = palans$site)
summary(PWsite.pal, confidence = 0.95, test.type = "dist")
summary(PWsite.pal, confidence = 0.95, test.type = "DL")
summary(PWsite.pal, confidence = 0.95, test.type = "VC", angle.type = "deg")

salinus <- geomorph.data.frame(subset.var$salinus,
                                plant = salinusLabels,
                                site = gsub('.{1}$', '', salinusLabels)) # geomorph data frame
twoD.sal <- two.d.array(subset.var$salinus)
salinus$Y <- ordinate(twoD.sal)$x[, 1:4]

fit.sal <- lm.rrpp(Y ~ site * plant, SS.type = "I",
                   data = salinus, print.progress = FALSE, iter = 999)

PWsite.sal <- pairwise(fit.sal, groups = kennedyi$site)
summary(PWsite.sal, confidence = 0.95, test.type = "dist")
summary(PWsite.sal, confidence = 0.95, test.type = "DL")
summary(PWsite.sal, confidence = 0.95, test.type = "VC", angle.type = "deg")
###############################################################################

#morphol.disparity?
#still a bit confused on interpreting this.
#The function estimates morphological disparity and performs pairwise comparisons to identify differences among groups. 
#Morphological disparity is estimated as the Procrustes variance, overall or for groups, using residuals of a linear model fit. 
#Procrustes variance is the same sum of the diagonal elements of the group covariance matrix divided 
#by the number of observations in the group (e.g., Zelditch et al. 2012).

# Morphological disparity for entire data set
morphol.disparity(coords ~ 1, groups = NULL, data = gdf,
                  iter = 999, print.progress = FALSE)
# Morphological disparity for entire data set, accounting for allometry
morphol.disparity(coords ~ Csize, groups= NULL, data = gdf,
                  iter = 999, print.progress = FALSE)

# Morphological partial disparities for overall mean
good1 <- morphol.disparity(coords ~ 1, groups= ~ variety, partial = TRUE,
                  data = gdf, iter = 999, print.progress = FALSE)
# Morphological NOT partial disparities for overal mean
good2 <- morphol.disparity(coords ~ 1, groups= ~ variety, partial = FALSE,
                           data = gdf, iter = 999, print.progress = FALSE)
# Morphological disparity without covariates, using group means
good3<- morphol.disparity(coords ~ variety, groups= ~variety,
                  data = gdf, iter = 999, print.progress = FALSE)
# With a previously defined ANOVA

good4<-morphol.disparity(f1 = fit2, groups = ~ variety, data = gdf,
                  iter = 999, print.progress = FALSE)



###############################################################################

#plotReftoTarget and plot mean shape per variety

#first take the mean shape of all coordinates.
msh <- mshape(lent2.gpa$coords)
plot(msh)

#here we take the mean shape of the two varieties that are separating from the rest 
meanaraneosus <- mshape(subset.var$araneosus)
meanaustralis <- mshape(subset.var$australis)
meanborreganus <- mshape(subset.var$borreganus)
meanfloribundus <- mshape(subset.var$floribundus)
meanfremontii <- mshape(subset.var$fremontii)
meankennedyii <- mshape(subset.var$kennedyi)
meanmaricopae <- mshape(subset.var$maricopae)
meannigricalycis <- mshape(subset.var$nigricalycis)
meanpalans <- mshape(subset.var$palans)
meansalinus <- mshape(subset.var$salinus)
meanvariablis <- mshape(subset.var$variabilis)
meanvitreus <- mshape(subset.var$vitreus)
meanwilsonii <- mshape(subset.var$wilsonii)
meanyuccanus <- mshape(subset.var$yuccanus)


plot(meanaraneosus)
plot(meanaustralis)
plot(meanborreganus)
plot(meanfloribundus)
plot(meanfremontii)
plot(meankennedyii)
plot(meanmaricopae)
plot(meannigricalycis)
plot(meanpalans)
plot(meansalinus)
plot(meanvariablis)
plot(meanvitreus)
plot(meanwilsonii)
plot(meanyuccanus)


gpa.araneosus <- gpagen(A = subset.var$araneosus,
       curves = aerial_curvers2,
       ProcD = FALSE)
gpa.australis <- gpagen(A = subset.var$australis,
       curves = aerial_curvers2,
       ProcD = FALSE)
gpa.borreganus <- gpagen(A = subset.var$borreganus,
       curves = aerial_curvers2,
       ProcD = FALSE)
gpa.floribundus <- gpagen(A = subset.var$floribundus,
       curves = aerial_curvers2,
       ProcD = FALSE)
gpa.fremontii <- gpagen(A = subset.var$fremontii,
       curves = aerial_curvers2,
       ProcD = FALSE)
gpa.kennedyii <- gpagen(A = subset.var$kennedyi,
       curves = aerial_curvers2,
       ProcD = FALSE)
gpa.maricopae <- gpagen(A = subset.var$maricopae,
       curves = aerial_curvers2,
       ProcD = FALSE)
gpa.nigricalycis <- gpagen(A = subset.var$nigricalycis,
       curves = aerial_curvers2,
       ProcD = FALSE)
gpa.palans <- gpagen(A = subset.var$palans,
       curves = aerial_curvers2,
       ProcD = FALSE)
gpa.salinus <- gpagen(A = subset.var$salinus,
       curves = aerial_curvers2,
       ProcD = FALSE)
gpa.variablis <- gpagen(A = subset.var$variabilis,
       curves = aerial_curvers2,
       ProcD = FALSE)
gpa.vitreus <- gpagen(A = subset.var2$vitreus,
       curves = aerial_curvers2,
       ProcD = FALSE)
gpa.wilsonii <- gpagen(A = subset.var$wilsonii,
       curves = aerial_curvers2,
       ProcD = FALSE)
gpa.wilsonii2 <- gpagen(A = subset.var2$wilsonii,
                       ProcD = FALSE)
gpa.yuccanus <- gpagen(A = subset.var$yuccanus,
       curves = aerial_curvers2,
       ProcD = FALSE)


plot(gpa.araneosus)
plot(gpa.australis)
plot(gpa.borreganus)
plot(gpa.floribundus)
plot(gpa.fremontii)
plot(gpa.kennedyii)
plot(gpa.maricopae)
plot(gpa.nigricalycis)
plot(gpa.palans)
plot(gpa.salinus)
plot(gpa.variablis)
plot(gpa.vitreus)
plot(gpa.wilsonii)
plot(gpa.wilsonii2)
plot(gpa.yuccanus)

#and then we plot how they differ to the mean and to each other.
plotRefToTarget(meanwilsonii, msh)
plotRefToTarget(meanmaricopae, msh)
plotRefToTarget(meanwilsonii, meanmaricopae)

###############################################################################

#Calculate Centroid size
#from what I understand, the Csize data is from before the procrustes transformation https://anatomypubs.onlinelibrary.wiley.com/doi/full/10.1002/ar.23065

#blank plot
plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')

#principal component of size-shape analysis
pc.plot <- pc
plotAllometry(fit2, size = gdf$Csize, logsz = TRUE, 
                         method = "PredLine", 
                         pch = 19, col = q_colors[as.factor(gdf$variety)])
legend("topright", legend=levels(interaction(gdf$variety)), pch=19, title='Varieties', col=q_colors)

plotAllometry(fitAllom, size = gdf$Csize, logsz = TRUE, 
              method = "PredLine", 
              pch = 19, col = q_colors[as.factor(gdf$variety)])
legend("topright", legend=levels(interaction(gdf$variety)), pch=19, title='Varieties', col=q_colors)


###############################################################################
#At the individual level:

#Plant 48A or 49G - wilsonii
p48C.gpa <- gpagen(A=subset.coords$`19-48C`, curves = aerial_curvers2, ProcD = FALSE)
p49G.gpa <- gpagen(A=subset.coords$`19-49G`, curves = aerial_curvers2, ProcD = FALSE)
plot(p48C.gpa)
plot(p49G.gpa)

#Plant 7B or 21B - borreganus
#WOAH check these out.....
p7B.gpa <- gpagen(A=subset.coords$`19-7B`, curves = aerial_curvers2, ProcD = FALSE)
p21C.gpa <- gpagen(A=subset.coords$`19-21C`, curves = aerial_curvers2, ProcD = FALSE)
plot(p7B.gpa)
plot(p21C.gpa)

#Troubleshooting
Read21C <- array(read.csv(file.choose(), header=TRUE))
XY21Clandmarks<- Read21C[,5:6]
plot(XY21Clandmarks[271:300,])

#Plant 72A or 66B or 23A - fremontii
p72A.gpa <- gpagen(A=subset.coords$`19-72A`, curves = aerial_curvers2, ProcD = FALSE)
p66B.gpa <- gpagen(A=subset.coords$`19-66B`, curves = aerial_curvers2, ProcD = FALSE)
p23A.gpa <- gpagen(A=subset.coords$`19-23A`, curves = aerial_curvers2, ProcD = FALSE)
plot(p72A.gpa)
plot(p66B.gpa)
plot(p23A.gpa)

#Plant 44A or 46C - vitreus
p44A.gpa <- gpagen(A=subset.coords$`19-44A`, curves = aerial_curvers2, ProcD = FALSE)
p46C.gpa <- gpagen(A=subset.coords$`19-46C`, curves = aerial_curvers2, ProcD = FALSE)
plot(p44A.gpa)
plot(p46C.gpa)

Read44A <- array(read.csv(file.choose(), header=TRUE))
XY44Alandmarks<- Read44A[,5:6]
plot(XY44Alandmarks[91:120,])

#Plant 40A or 2A - yuccanus
p40A.gpa <- gpagen(A=subset.coords$`19-40A`, curves = aerial_curvers2, ProcD = FALSE)
p2A.gpa <- gpagen(A=subset.coords$`19-2A`, curves = aerial_curvers2, ProcD = FALSE)
plot(p40A.gpa)
plot(p2A.gpa)


###############################################################################
#At the site level:
dim(my_arr)
length(sites3)
sites3 <-gsub('.{1}$', '', labels)
subset.sites <- coords.subset(A=my_arr, group=sites3)

#Plant 48A or 49G - wilsonii
s48.gpa <- gpagen(A=subset.sites$`19-48`, curves = aerial_curvers2, ProcD = FALSE)
s49.gpa <- gpagen(A=subset.sites$`19-49`, curves = aerial_curvers2, ProcD = FALSE)
plot(s48.gpa)
plot(s49.gpa)

#Plant 7B or 21B - borreganus
s7.gpa <- gpagen(A=subset.sites$`19-7`, curves = aerial_curvers2, ProcD = FALSE)
s21.gpa <- gpagen(A=subset.sites$`19-21`, curves = aerial_curvers2, ProcD = FALSE)
plot(s7.gpa)
plot(s21.gpa)

#Plant 72A or 66B or 23A - fremontii
s72.gpa <- gpagen(A=subset.sites$`19-72`, curves = aerial_curvers2, ProcD = FALSE)
s66.gpa <- gpagen(A=subset.sites$`19-66`, curves = aerial_curvers2, ProcD = FALSE)
s23.gpa <- gpagen(A=subset.sites$`19-23`, curves = aerial_curvers2, ProcD = FALSE)
plot(s72.gpa)
plot(s66.gpa)
plot(s23.gpa)

#Plant 44A or 46C - vitreus
s44.gpa <- gpagen(A=subset.sites$`19-44`, curves = aerial_curvers2, ProcD = FALSE)
s46.gpa <- gpagen(A=subset.sites$`19-46`, curves = aerial_curvers2, ProcD = FALSE)
plot(s44.gpa)
plot(s46.gpa)

#Plant 40A or 2A - yuccanus
s40.gpa <- gpagen(A=subset.sites$`19-40`, curves = aerial_curvers2, ProcD = FALSE)
s2.gpa <- gpagen(A=subset.sites$`19-2`, curves = aerial_curvers2, ProcD = FALSE)
plot(s40.gpa)
plot(s2.gpa)


###############################################################################

#Integration is the cohesion among traits that results from interactions of the biological 
#processes producing the phenotypic structures under study. Modularity refers to the relative 
#degrees of connectivity in systems—a module is a unit that is tightly integrated internally 
#but relatively independent from other such modules. In other words, modularity is about 
#differences in the degree of integration of parts within and between sets of traits.


integration.test(
  A,
  A2 = NULL,
  partition.gp = NULL,
  iter = 999,
  seed = 5,
  print.progress = TRUE
)

data(plethodon) 
Y.gpa<-gpagen(plethodon$land)    #GPA-alignment    
#landmarks on the skull and mandible assigned to partitions
land.gps<-c("A","A","A","A","A","B","B","B","B","B","B","B") 
IT <- integration.test(Y.gpa$coords, partition.gp=land.gps, iter=999)
summary(IT) # Test summary
plot(IT) # PLS plot

nrow(Y.gpa$coords)

###############################################################################

#Mantel test?







###############################################################################
#Beak characteristics

beaks <- newdata[,4:6]
sampleID <- as.vector(newdata[[1]])
div_sampleID <- sampleID[seq(1, length(sampleID), 30)]
beakSites <- gsub('.{1}$', '', div_sampleID)
beak_meta <- cbind(div_sampleID, beakSites)
beak_meta2 <- left_join(data.frame(beak_meta), data.frame(idVariety), by=c("div_sampleID"="Sample"))
#again I don't know what the issue with 43L is but there are always join issues
beak_meta2 <- beak_meta2[-c(1089:1098),]

nrow(beaks)/30
length(div_sampleID)
nrow(beak_meta2)

df.15 = beaks[seq(15, nrow(beaks), 30), ]
df.16 = beaks[seq(16, nrow(beaks), 30), ]
df.17 = beaks[seq(17, nrow(beaks), 30), ]


library(gdata)
library(multcompView)
mid.df <- interleave(df.16, df.17)

#create an id column for both dfs
#this id will dictate the order of rows
df.15$id <- 1:nrow(df.15)
mid.df$id <- rep(1:(nrow(mid.df)/2),each=2)
#rbind the two data.frames
fin.df <- rbind(df.15, mid.df)
#and now just order based on their ids
fin.df <- fin.df[order(fin.df$id), ]
dim(fin.df)

beak.df <- aperm(`dim<-`(t(fin.df), c(4, 3, (nrow(fin.df)/3))), c(2,1,3))
beak.df[1,2,3]
beak.df[1:3,,3]

#EXAMPLE
points <- fin.df[1:3]
x3 = x2 - x1
x <- fin.df[3,2] - fin.df[1,2]
y <- fin.df[3,3] - fin.df[1,3]
beakOpen <- sqrt((x^2)+(y^2))

beakOpenings <- c()
for (val in rep(1:2515)){
  #one fruit at a time
  mini <- beak.df[,,val]
  #hypotenuse
  x <-mini[3,2] - mini[1,2]
  y <- mini[3,3] - mini[1,3]
  beakOpen <- sqrt((x^2)+(y^2))
  beakOpenings <- c(beakOpenings,beakOpen)
}

write.csv(beakOpenings, file="beakOpenings.csv")

beakDepths <- c()
for (val in rep(1:2515)){
  #one fruit at a time
  mini <- beak.df[,,val]
  #calculate the midpoint
  x <-(mini[3,2]+mini[1,2])/2
  y <- (mini[3,3]+mini[1,3])/2
  #hypotenuse 
  xdiff <- x - mini[2,2]
  ydiff <- y - mini[2,3]
  beakDepth <- sqrt((xdiff^2)+(ydiff^2))
  #add to vector
  beakDepths <- c(beakDepths,beakDepth)
}

write.csv(beakDepths, file="beakDepths.csv")

#Combine everything back together
length(beakDepths)
length(beakOpenings)
length(labels)

length(mergeVar2[,2])

newVars.df <- data.frame(id=labels, beakOpening=beakOpenings, beakDepth=beakDepths)
newVars.df <- geomorph.data.frame(lent.gpa, site = sites3,
            plant = labels, variety= mergeVar2[,2], beakWidth = beakOpenings, beakDepth = beakDepths) # geomorph data frame
newVars.mat <- cbind(beakOpenings, beakDepths, lent.gpa$Csize, mergeVar2[,2], sites3)

newV.df <- geomorph.data.frame(lent.gpa, variety = beak_meta2[,3], sites= beak_meta2[,2], sample=beak_meta[,1],
                               beakWidth = beakOpenings, beakDepth = beakDepths)   

length(beak_meta$Variety)
length(beak_meta$beakSites)
length(beak_meta$div_sampleID)
length(beakOpenings)

write.csv(newVars.df, "/Users/QuinnThomas/Documents/Desktop_Items/Research/beakVals.csv", row.names = FALSE)

#Morphological disparity, the measure of morphological variation among species and higher taxa, has been at the core of an important research program in paleobiology over the last 25 years.
hist(newVars.df$beakWidth)
hist(newVars.df$beakDepth)

cor(beakDepths, beakOpenings)

plot(beakDepths, col=q_colors[as.factor(beak_meta$Variety)], pch=19)
plot(beakOpenings, col=q_colors[as.factor(beak_meta$Variety)], pch=19)
legend("topright", legend=levels(as.factor(beak_meta$Variety)), pch=19, title='Varieties', col=q_colors)


plot(c(0,0))

boxplot(newV.df$beakWidth ~ newV.df$variety,col = q_colors)
par(cex.axis=0.75)
boxplot(newV.df$beakDepth ~ newV.df$variety,col = q_colors)
summary(newV.df)
mat <- matrix(as.numeric(newVars.mat[,1:3]), ncol = 3)
measurePCA <- prcomp(mat, center = TRUE, scale. = TRUE)

model=lm( newV.df$beakWidth ~ newV.df$variety )
ANOVA=aov(model)
TUKEY<-TukeyHSD(ANOVA, 'newV.df$variety', conf.level=0.95)
plot(TUKEY , las=1 , col=col.rainbow[as.factor(newVars.df$variety)])


model2=lm( newV.df$beakDepth ~ newV.df$variety )
ANOVA2=aov(model2)
TUKEY2<-TukeyHSD(ANOVA2, 'newV.df$variety', conf.level=0.95)

generate_label_df <- function(TUKEY, variable){
  
  # Extract labels and factor levels from Tukey post-hoc 
  Tukey.levels <- TUKEY[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  
  #I need to put the labels in the same order as in the boxplot :
  Tukey.labels$treatment=rownames(Tukey.labels)
  Tukey.labels=Tukey.labels[order(Tukey.labels$treatment) , ]
  return(Tukey.labels)
}

# Apply the function on my dataset
LABELS <- generate_label_df(TUKEY, "newV.df$variety")
LABELS2 <- generate_label_df(TUKEY2, "newV.df$variety")

levels(as.factor(LABELS[,1]))
levels(as.factor(LABELS2[,1]))


Blues8 <- colorRampPalette(brewer.pal(9,"Blues"))(8)
Blues6 <- colorRampPalette(brewer.pal(9,"Blues"))(6)


beakBox <- k
boxplot(newV.df$beakWidth ~ newV.df$variety , ylim=c(min(newV.df$beakWidth) , 1.1*max(newV.df$beakWidth)) , col=Blues8[as.factor(LABELS[,1])] , ylab="Beak Opening Width", xlab="Variety", main="Beak Opening Width by variety")
legend("topleft", legend=levels(as.factor(LABELS[,1])), pch=19, title='Levels', col=Blues8)
beakBox2 <- k
boxplot(newV.df$beakDepth ~ newV.df$variety , ylim=c(min(newV.df$beakDepth) , 1.1*max(newV.df$beakDepth)) , col=Blues6[as.factor(LABELS2[,1])] , ylab="Beak Depth", xlab="Variety", main="Beak opening Depth by Variety")
legend("topleft", legend=levels(as.factor(LABELS2[,1])), pch=19, title='Levels', col=Blues6)


plot(measurePCA$x, col=col.rainbow[as.factor(newVars.df$variety)], pch=19)
legend("topright", legend=levels(as.factor(newVars.df$variety)), pch=19, title='Varieties', col=col.rainbow)
summary(measurePCA)

library(devtools)
install_github("vqv/ggbiplot")
library(ggbiplot)
ggbiplot(measurePCA, obs.scale = 1, var.scale = 1,
         groups = newVars.df$variety, ellipse = TRUE, circle = TRUE,ellipse.prob = 0.68) +
  scale_color_discrete(name = '') +
  theme(legend.direction = 'horizontal', legend.position = 'top')
print(g)
install.packages("psych")
library(psych)

