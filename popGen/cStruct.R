#remove.packages("rlang")
#install.packages("rlang", repos='http://cran.us.r-project.org')

#install.packages("data.table", repos='http://cran.us.r-project.org')
#install.packages("dplyr", repos='http://cran.us.r-project.org')
#install.packages("Matrix", repos='http://cran.us.r-project.org')
#install.packages("devtools", repos='http://cran.us.r-project.org')


library(data.table)
library(dplyr)
library(Matrix)

library(devtools)
#install_github("gbradburd/conStruct", ref="covariance_fix")
library("conStruct")

a_freqs <- as.matrix(read.csv("freq.csv", row.names=1, header=T))
crds <- as.matrix(read.csv("coords.csv", row.names=1))
gdists <- as.matrix(read.csv("gDist.csv", row.names=1))


prac2 <- conStruct(
  spatial = T,
  K = 2,
  freqs = a_freqs,
  coords=crds,
  geoDist = gdists,
  prefix = "",
  n.chains = ,
  n.iter = 50000,
  make.figs = T,
  save.files = T,
  control = setNames(list(0.99),"adapt_delta")
)






