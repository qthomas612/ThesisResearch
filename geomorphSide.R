install.packages("devtools")

library(devtools)

devtools::install_github('geomorphR/geomorph',ref="Stable")

library(geomorph)
library(data.table)
library(dplyr)

#For practice
data("plethspecies")
View(plethspecies$phy)

#Pull files from working directory
setwd("/Users/QuinnThomas/Documents/Desktop_Items/landmarkFiles/Side_CSV_files/")
wd = "/Users/QuinnThomas/Documents/Desktop_Items/landmarkFiles/Side_CSV_files/"
side_list_csv <- list.files(wd, pattern = "csv")
DT_side = do.call(rbind, lapply(side_list_csv, fread))
masterList <- read.csv("/Users/QuinnThomas/Documents/Desktop_Items/landmarkFiles/lentiginosusMaster.csv")
formatMaster <- masterList[c(1:13)]
idVariety <- masterList[c(2,3)]

#Remove missing data
sideData <- na.omit(DT_side) 
#confirm data type
class(sideData)

numSrows = nrow(sideData)

#I kind of forget what this is
sideSliced <- aperm(`dim<-`(t(sideData), c(6, 30, (numSrows/30))), c(2,1,3))


#land <- Sliced[,5:6,]
#tail(land)

#Save one individual fruit so we can set our sliding coordinates
FirstSample <- array(read.csv(file.choose(),
                              nrow=30,
                              header=TRUE))
XYSidelandmarks<- FirstSample[,5:6]


#This is conserved for all of our matrices from the same view AKA only do it once
#input an example xy matrix(hiihihihi) and interactively determine which
#landmarks are fixed and which ones slide(26 sliders)
side_curvers <- define.sliders(XYSidelandmarks, nsliders=29, write.file = FALSE)

#write.csv(side_curvers, file="side_curvers.csv")
#numericSide <- as.numeric(sideSliced[,5:6,])
#sideArr <- array(numericSide, dim = c(30,2,numSrows))
#check for missing values
#any(is.na(sideArr))
#noMissing <- estimate.missing(myarr,method="Reg")

#
no_MissingS <- na.omit(sideData)
any(is.na(no_MissingS))
sidelentiginosusLand <- as.matrix(no_MissingS[,c(5:6)])
A_side <- arrayspecs(sidelentiginosusLand, 30, 2)

#double check
any(is.na(A_side))
which(is.na(A_side))

#geomorph procrustes
side.gpa<- gpagen(
  A = A_side,
  curves = side_curvers,
  ProcD = FALSE)

#plot(procrustes)
plot(side.gpa)

sidePCA <- gm.prcomp(side.gpa$coords)

side_labels <- as.factor(no_MissingS[["sample_ID"]])
plot(sidePCA, col=side_labels)
levels(side_labels)

#EXAMPLE for extracting every nth element
#extracted_vec <- vec[seq(1, length(vec), 4)]
side_vec <- as.vector(no_MissingS[[1]])
side_labels.21 <- side_vec[seq(1, length(side_vec), 30)]
length(levels(as.factor(side_labels.21)))
levels(as.factor(side_labels.21))
side_labels <- recode(side_labels.21, "3E" = "19-3E", "3J" = "19-3J", "4E" = "19-4E", "4H" = "19-4H", 
                 "4I" = "19-4I", "4J" = "19-4J", "4K" = "19-4K", "4L" = "19-4L", "5H" = "19-5H",
                 "5J" = "19-5J", "5L" = "19-5L", "6C" = "19-6C", "6D" = "19-6D", "6F" = "19-6F", 
                 "6G" = "19-6G", "6H" = "19-6H", "6I" = "19-6I", "76A"="19-76A", "76B"="19-76B", "76C" ="19-76C", 
                 "76D"="19-76D", "76E"="19-76DE", "76F"="19-76F", "76G"="19-76G", "76H" ="19-76H", "761"="19-76J", "76K"="19-76K")
levels(as.factor(side_labels.2))

dddfS <- as.data.frame(side_labels)
mergeVarS2 <- merge(dddfS, idVariety, by.x = "side_labels", by.y = "Sample", no.dups=FALSE,  all.x=TRUE, all.y=FALSE)

#REMOVE: rows the extra 43L rows
mergeVarS2 <- mergeVarS2[-c(1008:1017),]

rowSNum = nrow(A_side)
side_arr <- array(A_side, dim = c(30,2,length(side_labels)))
length(side_arr)
subset.side.coords <- coords.subset(A=side_arr, group=side_labels)
ind_S_means <- as.array(lapply(subset.side.coords, mshape))

nrow(ind_S_means)

side_combined <- do.call(rbind,ind_S_means)
B_side <- arrayspecs(side_combined, 30,2)

siteSPlants <- names(ind_S_means)
pre_labelSFrame <- as.data.frame(siteSPlants)
sub_labelsSFrame <- sub("^","19-",pre_labelSFrame[241:267,1])
labelSFrame <- as.data.frame(c(pre_labelSFrame[1:240,1],sub_labelsSFrame))
colnames(pre_labelSFrame)<- "siteSplants"
mergeSVar <- merge(pre_labelSFrame, idVariety, by.x="siteSplants", by.y= "Sample", no.dups=TRUE,  all.x=TRUE, )
mergeSVar[240,2] <- "salinus" 
varSLabels <- as.vector(mergeSVar[["Variety"]])
#Remove palans 19-43L duplicate
varSLabels <- varSLabels[-c(130)]
write.csv(varSLabels, file="/Users/QuinnThomas/Documents/Desktop_Items/Research/varNamesSide.csv")

? merge
formattingNames <- sapply(strsplit(siteSPlants[1:240], "-"), '[', 2)
formatted_siteSPlants <- c(formattingNames, siteSPlants[241:267])
write.csv(formatted_siteSPlants, file="/Users/QuinnThomas/Documents/Desktop_Items/Research/indNamesSide.csv")

#class(subset.var)
S_sites <-gsub('.{1}$', '', siteSPlants)
S_sites2 <-gsub('.{1}$', '', side_labels)

length(siteSPlants)
nrow(labelSFrame)
nrow(mergeSVar)
length(S_sites)
length(varSLabels)

side2.gpa<- gpagen(
  A = B_side,
  curves = side_curvers,
  ProcD = FALSE)


subset.S.var <- coords.subset(A=side2.gpa$coords, group=varSLabels)

plot(side2.gpa)
pcaS2 <- gm.prcomp(side2.gpa$coords)

write.csv(pcaS2$x, "/Users/QuinnThomas/Documents/Desktop_Items/Research/Side_PCvals.csv", row.names = FALSE)

plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
plot(pcaS2$shapes$shapes.comp1$min)
plot(pcaS2$shapes$shapes.comp1$max)
plot(pcaS2$shapes$shapes.comp2$min)
plot(pcaS2$shapes$shapes.comp2$max)

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


plot(pcaS2, col=q_colors[as.factor(varSLabels)], pch=19)
facSVar <- as.factor(varSLabels)
legend("topright", legend=levels(as.factor(varSLabels)), pch=19, title='Varieties', col=q_colors)

##########################################################################################

gmatS <- two.d.array(side2.gpa$coords)
gmatS2 <- two.d.array(side.gpa$coords)

dim(gmatS)
length(varSLabels)
length(S_sites)
dim(gmatS2)
length(as.vector(mergeVarS2[,2]))
length(S_sites2)

#plant df
rrppS.df1<-rrpp.data.frame(shape=gmatS,CS=side2.gpa$Csize,variety=varSLabels, site=S_sites)
#fruit df
rrppS.df2<-rrpp.data.frame(shape=gmatS2,CS=side.gpa$Csize,variety=as.vector(mergeVarS2[,2]), site=S_sites2, individual=side_labels)

#Nested ANOVAs  by plant
rrS.nest1 <- lm.rrpp(shape~variety/site+site, data=rrppS.df1,iter=999,print.progress=F)
anova(rrS.nest1)
coef(rrS.nest1, test = TRUE)

#Nested ANOVAs  by fruit
rrS.nest4 <- lm.rrpp(shape~(variety/site), data=rrppS.df2,iter=999,print.progress=F)
anova(rrS.nest4)
coef(rrS.nest4, test = TRUE)
rrS.nest5 <- lm.rrpp(shape~(site/individual)+individual, data=rrppS.df2,iter=999,print.progress=F)
anova(rrS.nest5)
coef(rrS.nest5, test = TRUE)

##########################################################################################

#Pairwise comparisons
gdfS <- geomorph.data.frame(side2.gpa,
                           site = S_sites,
                           variety = varSLabels) # geomorph data frame
twoDS.lent <- two.d.array(side2.gpa$coords)
gdfS$Y <- ordinate(twoDS.lent)$x[, 1:3]

fitS1 <- lm.rrpp(Y ~(variety/site), SS.type = "I",
                data = gdfS, print.progress = FALSE, iter = 999)
fitS2 <- lm.rrpp(Y ~(site), SS.type = "I",
                data = gdfS, print.progress = FALSE, iter = 999)

sfit1.var <- procD.lm(coords ~ variety,
                     data = gdfS, iter = 999, turbo = TRUE,
                     RRPP = TRUE, print.progress = FALSE) # randomize residuals


PW1S <- pairwise(sfit1.var, groups = sfit1.var$data$variety)
summary(PW1S)
PW1.2S <- pairwise(fitS2, groups = gdfS$site)
pairwiseSite1S<-summary(PW1S, confidence = 0.95, test.type = "dist")
pairwiseSite1.2S<-summary(PW1.2S, confidence = 0.95, test.type = "dist")
View(pairwiseSite1S[["pairwise.tables"]][["P"]])
View(pairwiseSite1.2S[["pairwise.tables"]][["P"]])

##########################################################################################

#Morphological Disparity

# With a previously defined ANOVA
side_dispar<-morphol.disparity(f1 = sfit1.var, groups = ~ variety, data = gdfS,
                         iter = 999, print.progress = FALSE)
write.csv(as.matrix(side_dispar[["PV.dist"]]),file = "disp_side.csv")
write.csv(as.matrix(side_dispar[["PV.dist.Pval"]]),file = "disp_side_p.csv")
##########################################################################################

#plotReftoTarget and plot mean shape per variety

#first take the mean shape of all coordinates.
mshS <- mshape(side2.gpa$coords)
plot(mshS)

#here we take the mean shape of the two varieties that are separating from the rest 
meanSaraneosus <- mshape(subset.S.var$araneosus)
meanSaustralis <- mshape(subset.S.var$australis)
meanSborreganus <- mshape(subset.S.var$borreganus)
meanSfloribundus <- mshape(subset.S.var$floribundus)
meanSfremontii <- mshape(subset.S.var$fremontii)
meanSkennedyii <- mshape(subset.S.var$kennedyi)
meanSmaricopae <- mshape(subset.S.var$maricopae)
meanSnigricalycis <- mshape(subset.S.var$nigricalycis)
meanSpalans <- mshape(subset.S.var$palans)
meanSsalinus <- mshape(subset.S.var$salinus)
meanSvariablis <- mshape(subset.S.var$variabilis)
meanSvitreus <- mshape(subset.S.var$vitreus)
meanSwilsonii <- mshape(subset.S.var$wilsonii)
meanSyuccanus <- mshape(subset.S.var$yuccanus)

plot(meanSaraneosus)
plot(meanSaustralis)
plot(meanSborreganus)
plot(meanSfloribundus)
plot(meanSfremontii)
plot(meanSkennedyii)
plot(meanSmaricopae)
plot(meanSnigricalycis)
plot(meanSpalans)
plot(meanSsalinus)
plot(meanSvariablis)
plot(meanSvitreus)
plot(meanSwilsonii)
plot(meanSyuccanus)

gpaS.araneosus <- gpagen(A = subset.S.var$araneosus,
                        curves = side_curvers,
                        ProcD = FALSE)
gpaS.australis <- gpagen(A = subset.S.var$australis,
                         curves = side_curvers,
                         ProcD = FALSE)
gpaS.borreganus <- gpagen(A = subset.S.var$borreganus,
                         curves = side_curvers,
                         ProcD = FALSE)
gpaS.floribundus <- gpagen(A = subset.S.var$floribundus,
                          curves = side_curvers,
                          ProcD = FALSE)
gpaS.fremontii <- gpagen(A = subset.S.var$fremontii,
                        curves = side_curvers,
                        ProcD = FALSE)
gpaS.kennedyii <- gpagen(A = subset.S.var$kennedyi,
                        curves = side_curvers,
                        ProcD = FALSE)
gpaS.maricopae <- gpagen(A = subset.S.var$maricopae,
                        curves = side_curvers,
                        ProcD = FALSE)
gpaS.nigricalycis <- gpagen(A = subset.S.var$nigricalycis,
                         curves = side_curvers,
                         ProcD = FALSE)
gpaS.palans <- gpagen(A = subset.S.var$palans,
                     curves = side_curvers,
                     ProcD = FALSE)
gpaS.salinus <- gpagen(A = subset.S.var$salinus,
                      curves = side_curvers,
                      ProcD = FALSE)
gpaS.variablis <- gpagen(A = subset.S.var$variabilis,
                        curves = side_curvers,
                        ProcD = FALSE)
gpaS.vitreus <- gpagen(A = subset.S.var$vitreus,
                      curves = side_curvers,
                      ProcD = FALSE)
gpaS.wilsonii <- gpagen(A = subset.S.var$wilsonii,
                       curves = side_curvers,
                       ProcD = FALSE)
gpaS.wilsonii2 <- gpagen(A = subset.S.var$wilsonii,
                        ProcD = FALSE)
gpaS.yuccanus <- gpagen(A = subset.S.var$yuccanus,
                       curves = side_curvers,
                       ProcD = FALSE)

plot(gpaS.araneosus)
plot(gpaS.australis)
plot(gpaS.borreganus)
plot(gpaS.floribundus)
plot(gpaS.fremontii)
plot(gpaS.kennedyii)
plot(gpaS.maricopae)
plot(gpaS.nigricalycis)
plot(gpaS.palans)
plot(gpaS.salinus)
plot(gpaS.variablis)
plot(gpaS.vitreus)
plot(gpaS.wilsonii)
plot(gpaS.wilsonii2)
plot(gpaS.yuccanus)

##########################################################################################

"Redundancy analysis (RDA) was performed using the package vegan (Oksanen et al., 2019) as an alternate method to
testing the association between genetic distance and environmental and morphological
variation. The RDA used principal coordinate analysis axes from microsatellite genotypes
as the predicted variables and environmental and morphological variables as predictor
variables."

###############################################################################

#Calculate Centroid size
#from what I understand, the Csize data is from before the procrustes transformation https://anatomypubs.onlinelibrary.wiley.com/doi/full/10.1002/ar.23065

#blank plot
plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')

#principal component of size-shape analysis
fitAllom2 <- procD.lm(coords ~ log(Csize), data=gdfS, print.progress = FALSE)
summary(fitAllom2)
pcS.plot <- plotAllometry(sfit1.var, size = log(gdfS$Csize), logsz = TRUE, 
                         method = "PredLine", 
                         pch = 19, col = mycols[as.factor(gdfS$variety)])
legend("topright", legend=levels(interaction(gdfS$variety)), pch=19, title='Varieties', col=mycols)











