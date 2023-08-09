setwd("~/Documents/Desktop_Items/Research/")
setwd("~/Downloads/2023_lent/May_lent/")

install.packages(c("ggplot2","gridExtra","label.switching","tidyr","remotes"),repos="https://cloud.r-project.org")

# install pophelper package from GitHub
remotes::install_github('royfrancis/pophelper')

library(pophelper)
library(rlang)
library(dplyr)

setwd("~/Downloads/2023_lent/combined_lent/")
setwd("~/Downloads/2023_lent/May_lent/")
setwd("~/Documents/Desktop_Items/Research/final_structure/")

q_colors <- c("dodgerblue2", "#E31A1C", # red
                           "green4",
                           "#6A3D9A", # purple
                           "#FF7F00", # orange
                           "gold1",
                           "skyblue2", "#FB9A99", # lt pink
                           "palegreen2",
                           "#CAB2D6", # lt purple
                           #"#FDBF6F", # lt orange
                           "gray70",
                           "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
                           "darkturquoise", "green1", "yellow4", "yellow3",
                           "darkorange4", "brown"
)

#list all structure files in the folder
com_files <- list.files(("~/Downloads/2023_lent/combined_lent/filtered_r70/current_outfiles"), pattern = "*_f", full.names=TRUE)
MayK4_files <- list.files(("~/Downloads/2023_lent/May_lent"), pattern = "*_f", full.names=TRUE)
Jan_files <- list.files("~/Documents/Desktop_Items/Research/final_structure", pattern = "*_f", full.names=TRUE)
o_sfiles <- sfiles[c(2:3,1)]
K4_files <- sfiles[c(10:12)]


#makes the above command simpler
allfiles <-readQ(com_files,indlabfromfile=TRUE)
May_files <-readQ(MayK4_files[1:3],indlabfromfile=TRUE)
#subM4 <-readQ("combined_subMay/CM_out_4_f",indlabfromfile=TRUE)
#subJ4 <-readQ("combined_subJan/CJ_out_4_f",indlabfromfile=TRUE)
Jan_files <- readQ(Jan_files,indlabfromfile=TRUE)

#tab <- tabulateQ(qlist=allfiles)
#tabulateQ(qlist = M3NPr80k4)

class(allfiles)

#allfiles$`2_pop2.2_out_f`

#Practice plotting
plotQ(qlist=allfiles[8], exportplot = FALSE, returnplot = TRUE ,exportpath=NULL)
plotQ(qlist=allfiles[8],sortind="all",exportplot = FALSE, returnplot = TRUE ,exportpath=NULL,
      showindlab=TRUE, useindlab=T)

#Read in group information including variety and site
'''
grpinfoNPn <- as.data.frame(read.csv(file="../sampleNames_popmap2.2.csv", header = T))
#formatting
grpinfoNPn[29,1] <- "1E-1"
grpinfoNPn[32,1] <- "20E-2"

row.names(grpinfoNPn) <- grpinfoNPn[,1]
grpinfoNP <- grpinfoNPn[,2:3]
colnames(grpinfoNP) <- c("variety", "site")
class(grpinfoNP)
'''


#grpinfoNPn <- as.data.frame(read.csv(file="~/Downloads/2023_lent/May_lent/struct_MayK2.csv", header = T))
May_grpInfo <- as.data.frame(read.csv(file = "~/Downloads/2023_lent/combined_lent/old/MayK4.csv", header=T))
combined_grpInfo <- as.data.frame(read.csv(file = "~/Downloads/2023_lent/combined_lent/old/r70_output/r70_K4.csv", header=T))
#subM_grpInfo <- as.data.frame(read.csv(file = "combined_subMay/CM_struct4.csv", header=T))
#subJ_grpInfo <- as.data.frame(read.csv(file = "combined_subJan/CJ_struct4.csv", header=T))
Jan_grpInfo <- as.data.frame(read.csv(file = "~/Documents/Desktop_Items/Research/sampleNames_Jan30.csv", header=T))

com_ginfo<-  combined_grpInfo[,1,drop=FALSE]
com_ginfo2 <- combined_grpInfo[,1:2]
may_ginfo2 <- May_grpInfo[,1:2]
#subm_ginfo2 <- subM_grpInfo[,1:2]
#subj_ginfo2 <- subJ_grpInfo[,1:2]
jan_ginfo2 <- Jan_grpInfo[,1:2]

Q_order <- as.data.frame(row.names(allfiles$out_3.1_f))
colnames(Q_order) <- "Sample_ID"
com_ordered_groups <- inner_join(Q_order, com_ginfo2)
com_o_variety <- as.data.frame(com_ordered_groups[,2])
#o_site <- as.data.frame(ordered_groups[,3])

May_order <- as.data.frame(row.names(MK4_files$May_out_4.1_f))
colnames(May_order) <- "Sample_ID"
May_ordered_groups <- inner_join(May_order, may_ginfo2)
May_o_variety <- as.data.frame(May_ordered_groups[,2])

Jan_order <- as.data.frame(row.names(Jan_files$`4_pop2.2_out_f`))
colnames(Jan_order) <- "Sample_ID"
Jan_ordered_groups <- inner_join(Jan_order, jan_ginfo2)
Jan_o_variety <- as.data.frame(Jan_ordered_groups[,2])

# subM_order <- as.data.frame(row.names(subM4$CM_out_4_f))
# colnames(subM_order) <- "Sample_ID"
# subM_ordered_groups <- inner_join(subM_order, subm_ginfo2)
# subM_o_variety <- as.data.frame(subM_ordered_groups[,2])
# 
# subJ_order <- as.data.frame(row.names(subJ4$CJ_out_4_f))
# colnames(subJ_order) <- "Sample_ID"
# subJ_ordered_groups <- inner_join(subJ_order, subj_ginfo2)
# subJ_o_variety <- as.data.frame(subJ_ordered_groups[,2])

out3 <- allfiles[2]
out2 <- allfiles[1]
out3$outfile3.1_f <- out3$outfile3.1_f[1:180,]
out2$outfile2.1_f <- out2$outfile2.1_f[1:180,]

#OK here is the real shit
qgraph <- plotQ(qlist=allfiles[8],sortind="label",exportplot = FALSE, returnplot = TRUE,
      showindlab=TRUE, useindlab=T, grplab = com_o_variety,
      ordergrp = T,returndata=T)
plot(qgraph$plot[[1]])
# qgraphSJ <- plotQ(qlist=subJ4[1],sortind="label",exportplot = FALSE, returnplot = TRUE,
#                 showindlab=TRUE, useindlab=T, grplab = subJ_o_variety,
#                 ordergrp = T,returndata=T)
# qgraphSM <- plotQ(qlist=subM4[1],sortind="label",exportplot = FALSE, returnplot = TRUE,
#                   showindlab=TRUE, useindlab=T, grplab = subM_o_variety,
#                   ordergrp = T,returndata=T)
# plot(qgraphSJ$plot[[1]])
# plot(qgraphSM$plot[[1]])


qgraphM <- plotQ(qlist=MK4_files[2],sortind="label",exportplot = FALSE, returnplot = TRUE,
                showindlab=TRUE, useindlab=T, grplab = May_o_variety,
                ordergrp = T,returndata=T)
qgraphJ <- plotQ(qlist=Jan_files[3],sortind="label",exportplot = FALSE, returnplot = TRUE,
                 showindlab=TRUE, useindlab=T, grplab = Jan_o_variety,
                 ordergrp = T,returndata=T)
plot(qgraphM$plot[[1]])
plot(qgraphJ$plot[[1]])


#Jan30
qgraph <- plotQ(qlist=allfiles[8],sortind="label",exportplot = FALSE, returnplot = TRUE,
                showindlab=TRUE, useindlab=T, grplab = com_o_variety,
                ordergrp = T,returndata=T)
plot(qgraph$plot[[1]])
qgraphk <- plotQ(qlist=allfiles[8],sortind="all",exportplot = FALSE, returnplot = TRUE,
                 showindlab=TRUE, useindlab=T,returndata=T,
                 showsp=F, showlegend = T, basesize=4,
                 clustercol=q_colors[1:4],
                 showtitle = T, titlelab = "Combined Structure K = 3", titlehjust = 0.5,)
plot(qgraphk$plot[[1]])

qgraphk <- plotQ(qlist=May_files[2],sortind="all",exportplot = FALSE, returnplot = TRUE,
                 showindlab=TRUE, useindlab=T,returndata=T,
                 showsp=F, showlegend = T, basesize=4,
                 clustercol=q_colors[1:4],
                 showtitle = T, titlelab = "May Structure K = 3", titlehjust = 0.5,)
plot(qgraphk$plot[[1]])



# #sort by k
# qgraphk <- plotQ(qlist=allfiles_M2[3],sortind="all",exportplot = FALSE, returnplot = TRUE,
#                 showindlab=TRUE, useindlab=T,returndata=T,
#                 showsp=F, showlegend = T, basesize=4,
#                 clustercol=q_colors[1:5],
#                 showtitle = T, titlelab = "Structure K = 5", titlehjust = 0.5,)
# plot(qgraphk$plot[[1]])

#Looking at the above plot we can visualize which varieties primarily belong to which groups so let's reorder them

#reds, variablis -dark blue, araneosus - green
grpsubred <- c("araneosus", "floribundus", "fremontii", "kennedyi", "lentiginosus", "maricopae", "palans", "salinus", "vitreus", "wilsonii", "variablis")
#light blue floribundus -dark blue, ineptus -dark blue, palans -green, salinus -green
grpsublb <- c("antonius", "nigricalycis")
#dark blue, variablis -red, floribundus -lightblue, ineptus -lightblue
grpsubdb <- c("borreganus")
#green, borreganus -red, palans -light blue, salinus -lightblue
grpsubgreen <- c("australis", "yuccanus")


qsorder <- c("araneosus", "vitreus", "palans", "floribundus", "kennedyi", "fremontii", "lentiginosus", "salinus", "maricopae","wilsonii", "variablis",
             "antonius", "nigricalycis",
             "borreganus",
             "australis", "yuccanus")
qsorder2 <- c("araneosus", "vitreus", "maricopae","wilsonii",
             "palans", "floribundus", "kennedyi", "lentiginosus", "salinus", "fremontii", "variablis",
             "antonius", "nigricalycis",
             "borreganus",
             "australis", "yuccanus")

# qsorder2 <- c("antonius", "yuccanus", "borreganus", "australis", "variablis",
#              "nigricalycis", "semotus", "floribundus", "ineptus", 
#              "fremontii", "kennedyi","salinus", "lentiginosus", "palans", 
#              "maricopae", "wilsonii", "vitreus", "araneosus")


qgraph2 <- plotQ(qlist=as.qlist(qgraph$data$qlist),sortind="Cluster2",exportplot = FALSE, returnplot = TRUE,
                showindlab=TRUE, useindlab=T, grplab = data.frame(qgraph$data$grplab),
                ordergrp = F,returndata=T, subsetgrp = qsorder2, 
                showsp=F, showlegend = T, basesize=5, grplabsize=2,
                showtitle = T, titlelab = "Structure K = 4", titlehjust = 0.5,
                clustercol=q_colors[1:5],
                grplabspacer = 0, grplabpos = 0.60, grplabheight = NA)
plot(qgraph2$plot[[1]])

# qgraphsubJ2 <- plotQ(qlist=as.qlist(qgraphSJ$data$qlist),sortind="Cluster2",exportplot = FALSE, returnplot = TRUE,
#                  showindlab=TRUE, useindlab=T, grplab = data.frame(qgraphSJ$data$grplab),
#                  ordergrp = F,returndata=T, subsetgrp = qsorder, 
#                  showsp=F, showlegend = T, basesize=5, grplabsize=2,
#                  showtitle = T, titlelab = "Structure K = 4", titlehjust = 0.5,
#                  clustercol=q_colors[1:5],
#                  grplabspacer = 0, grplabpos = 0.60, grplabheight = NA)
# plot(qgraphsubJ2$plot[[1]])
# 
# qgraphsubM2 <- plotQ(qlist=as.qlist(qgraphSM$data$qlist),sortind="Cluster2",exportplot = FALSE, returnplot = TRUE,
#                  showindlab=TRUE, useindlab=T, grplab = data.frame(qgraphSM$data$grplab),
#                  ordergrp = F,returndata=T, subsetgrp = qsorder, 
#                  showsp=F, showlegend = T, basesize=5, grplabsize=2,
#                  showtitle = T, titlelab = "Structure K = 4", titlehjust = 0.5,
#                  clustercol=q_colors[1:5],
#                  grplabspacer = 0, grplabpos = 0.60, grplabheight = NA)
# plot(qgraphsubM2$plot[[1]])

qgraphM2 <- plotQ(qlist=as.qlist(qgraphM$data$qlist),sortind="Cluster2",exportplot = FALSE, returnplot = TRUE,
                 showindlab=TRUE, useindlab=T, grplab = data.frame(qgraphM$data$grplab),
                 ordergrp = F,returndata=T, subsetgrp = qsorder, 
                 showsp=F, showlegend = T, basesize=5, grplabsize=2,
                 showtitle = T, titlelab = "Structure K = 3", titlehjust = 0.5,
                 clustercol=q_colors[1:5],
                 grplabspacer = 0, grplabpos = 0.60, grplabheight = NA)
plot(qgraphM2$plot[[1]])

qgraphJ2 <- plotQ(qlist=as.qlist(qgraphJ$data$qlist),sortind="Cluster2",exportplot = FALSE, returnplot = TRUE,
                  showindlab=TRUE, useindlab=T, grplab = data.frame(qgraphJ$data$grplab),
                  ordergrp = F,returndata=T, subsetgrp = qsorder, 
                  showsp=F, showlegend = T, basesize=5, grplabsize=2,
                  showtitle = T, titlelab = "Structure K = 4", titlehjust = 0.5,
                  clustercol=q_colors[1:5],
                  grplabspacer = 0, grplabpos = 0.60, grplabheight = NA)
plot(qgraphJ2$plot[[1]])

d1<-data.frame(qgraph$data$grplab)
d1$Sample_ID <-row.names.data.frame(d1)
mergedorder <- left_join(d1, combined_grpInfo[,c(1,3)])
glaborder <- mergedorder[,-c(2)]
#glaborder[203,2] <- "2"
colnames(glaborder) <- c("variety","site")

# qsorder <- c("antonius", "yuccanus", "borreganus", "australis", "variablis",
# "nigricalycis", "semotus", "floribundus", "ineptus", 
# "fremontii", "kennedyi","salinus", "palans", 
# "maricopae", "wilsonii", "vitreus", "lentiginosus", "araneosus")


# this is to plot individual varieties
qgraphvar1 <- plotQ(qlist=as.qlist(qgraph$data$qlist),sortind="label", returnplot = TRUE,
                   showindlab=TRUE, useindlab=T, grplab = glaborder,
                   ordergrp = T,returndata=T, subsetgrp = "fremontii", 
                   basesize=11,
                   showsp=F, showlegend = T, showgrplab = F,
                   showtitle = T, titlelab = "Structure K = 4", titlehjust = 0.5,
                   grplabspacer = 0, grplabpos = 0.25, grplabheight = NA,
                   exportplot = F, exportpath = getwd())
#}

#plot(qgraphvar1$plot[[1]])
glab2 <- as.data.frame(qgraphvar1[["data"]][["grplab"]][[1]][["site"]])
colnames(glab2) <- "site"
qgraphvar2 <- plotQ(qlist=as.qlist(qgraphvar1$data$qlist),sortind="Cluster4", returnplot = TRUE,
                    showindlab=TRUE, useindlab=T, grplab = glab2,
                    ordergrp = T,returndata=T, 
                    basesize=11,
                    showsp=F, showlegend = T, showgrplab = T, showyaxis = F,
                    showtitle = T, titlelab = "yuccanus K = 4", titlehjust = 0.5,
                    clustercol=q_colors[1:5],
                    grplabspacer = 0, grplabpos = 0.25, grplabheight = NA,grplabsize=3.5,
                    exportplot = F)
plot(qgraphvar2$plot[[1]])



##### Create a joined plot
align_allfiles <- alignK(allfiles)


M2_join <- plotQ(qlist=align_allfiles,grplab=o_variety,imgoutput="join",
      sharedindlab=F,exportplot=F, returnplot = T, showlegend = T,
      showindlab = T, useindlab = T, ordergrp = T, subsetgrp = qsorder,
      grplabspacer = 0, grplabpos = 0.25, grplabheight = NA,grplabsize=1.5,
      basesize=4,
      splab=c("K=2", "K=4", "K=16"),
      grplabangle = 45, clustercol=q_colors[1:16])
plot(M2_join$plot[[1]])

