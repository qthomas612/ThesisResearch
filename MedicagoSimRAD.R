#medicago --> truncatula , 419Mb, 32.6%GC
#Chromosome sizes in order: 56, 51, 58, 64, 45, 43, 56, 50
library(SimRAD)

#Medicago genome fna file
Medgen.sq <- ref.DNAseq("~/Downloads/Mtrunc.fa")

#Medicago chromosomes
MedCh1.sq <- ref.DNAseq("~/Downloads/MedtruncCh1.fasta")
MedCh2.sq <- ref.DNAseq("~/Downloads/MedtruncCh2.fasta")
MedCh3.sq <- ref.DNAseq("~/Downloads/MedtruncCh3.fasta")
MedCh4.sq <- ref.DNAseq("~/Downloads/MedtruncCh4.fasta")
MedCh5.sq <- ref.DNAseq("~/Downloads/MedtruncCh5.fasta")
MedCh6.sq <- ref.DNAseq("~/Downloads/MedtruncCh6.fasta")
MedCh7.sq <- ref.DNAseq("~/Downloads/MedtruncCh7.fasta")
MedCh8.sq <- ref.DNAseq("~/Downloads/MedtruncCh8.fasta")

#Restriction Enzyme 1
#EcoRI
cs_5p1 <- "G"
cs_3p1 <- "AATTC"
#Restriction Enzyme 2
#MseI
cs_5p2 <- "T"
cs_3p2 <- "TAA"

#Exp. reads per lane / (number of samples * number of loci)
#http://seqanswers.com/forums/archive/index.php/t-17121.html

Gen.dig <- insilico.digest(Medgen.sq, cs_5p1, cs_3p1, cs_5p2, cs_3p2, verbose=TRUE)

#digestion with EcoRI and MseI
MedCh1.dig <- insilico.digest(MedCh1.sq, cs_5p1, cs_3p1, cs_5p2, cs_3p2, verbose=TRUE)
MedCh2.dig <- insilico.digest(MedCh2.sq, cs_5p1, cs_3p1, cs_5p2, cs_3p2, verbose=TRUE)
MedCh3.dig <- insilico.digest(MedCh3.sq, cs_5p1, cs_3p1, cs_5p2, cs_3p2, verbose=TRUE)
MedCh4.dig <- insilico.digest(MedCh4.sq, cs_5p1, cs_3p1, cs_5p2, cs_3p2, verbose=TRUE)
MedCh5.dig <- insilico.digest(MedCh5.sq, cs_5p1, cs_3p1, cs_5p2, cs_3p2, verbose=TRUE)
MedCh6.dig <- insilico.digest(MedCh6.sq, cs_5p1, cs_3p1, cs_5p2, cs_3p2, verbose=TRUE)
MedCh7.dig <- insilico.digest(MedCh7.sq, cs_5p1, cs_3p1, cs_5p2, cs_3p2, verbose=TRUE)
MedCh8.dig <- insilico.digest(MedCh8.sq, cs_5p1, cs_3p1, cs_5p2, cs_3p2, verbose=TRUE)

#fragment type selection (with both cutsites)
MedCh1.sel <- adapt.select(MedCh1.dig, type="AB+BA", cs_5p1, cs_3p1, cs_5p2, cs_3p2)
MedCh2.sel <- adapt.select(MedCh2.dig, type="AB+BA", cs_5p1, cs_3p1, cs_5p2, cs_3p2)
MedCh3.sel <- adapt.select(MedCh3.dig, type="AB+BA", cs_5p1, cs_3p1, cs_5p2, cs_3p2)
MedCh4.sel <- adapt.select(MedCh4.dig, type="AB+BA", cs_5p1, cs_3p1, cs_5p2, cs_3p2)
MedCh5.sel <- adapt.select(MedCh5.dig, type="AB+BA", cs_5p1, cs_3p1, cs_5p2, cs_3p2)
MedCh6.sel <- adapt.select(MedCh6.dig, type="AB+BA", cs_5p1, cs_3p1, cs_5p2, cs_3p2)
MedCh7.sel <- adapt.select(MedCh7.dig, type="AB+BA", cs_5p1, cs_3p1, cs_5p2, cs_3p2)
MedCh8.sel <- adapt.select(MedCh8.dig, type="AB+BA", cs_5p1, cs_3p1, cs_5p2, cs_3p2)

#wide size selection (200-400)
wid.Ch1 <- size.select(MedCh1.sel, min.size = 200, max.size = 400, graph=TRUE, verbose=TRUE)
wid.Ch2 <- size.select(MedCh2.sel, min.size = 200, max.size = 400, graph=TRUE, verbose=TRUE)
wid.Ch3 <- size.select(MedCh3.sel, min.size = 200, max.size = 400, graph=TRUE, verbose=TRUE)
wid.Ch4 <- size.select(MedCh4.sel, min.size = 200, max.size = 400, graph=TRUE, verbose=TRUE)
wid.Ch5 <- size.select(MedCh5.sel, min.size = 200, max.size = 400, graph=TRUE, verbose=TRUE)
wid.Ch6 <- size.select(MedCh6.sel, min.size = 200, max.size = 400, graph=TRUE, verbose=TRUE)
wid.Ch7 <- size.select(MedCh7.sel, min.size = 200, max.size = 400, graph=TRUE, verbose=TRUE)
wid.Ch8 <- size.select(MedCh8.sel, min.size = 200, max.size = 400, graph=TRUE, verbose=TRUE)

#narrow size selection (250-350)
nar.MCh1 <- size.select(MedCh1.sel, min.size = 250, max.size = 350, graph=TRUE, verbose=TRUE)
nar.MCh2 <- size.select(MedCh2.sel, min.size = 250, max.size = 350, graph=TRUE, verbose=TRUE)
nar.MCh3 <- size.select(MedCh3.sel, min.size = 250, max.size = 350, graph=TRUE, verbose=TRUE)
nar.MCh4 <- size.select(MedCh4.sel, min.size = 250, max.size = 350, graph=TRUE, verbose=TRUE)
nar.MCh5 <- size.select(MedCh5.sel, min.size = 250, max.size = 350, graph=TRUE, verbose=TRUE)
nar.MCh6 <- size.select(MedCh6.sel, min.size = 250, max.size = 350, graph=TRUE, verbose=TRUE)
nar.MCh7 <- size.select(MedCh7.sel, min.size = 250, max.size = 350, graph=TRUE, verbose=TRUE)
nar.MCh8 <- size.select(MedCh8.sel, min.size = 250, max.size = 350, graph=TRUE, verbose=TRUE)

#narrow size selection (350-450)
narn.MCh1 <- size.select(MedCh1.sel, min.size = 350, max.size = 450, graph=TRUE, verbose=TRUE)
narn.MCh2 <- size.select(MedCh2.sel, min.size = 350, max.size = 450, graph=TRUE, verbose=TRUE)
narn.MCh3 <- size.select(MedCh3.sel, min.size = 350, max.size = 450, graph=TRUE, verbose=TRUE)
narn.MCh4 <- size.select(MedCh4.sel, min.size = 350, max.size = 450, graph=TRUE, verbose=TRUE)
narn.MCh5 <- size.select(MedCh5.sel, min.size = 350, max.size = 450, graph=TRUE, verbose=TRUE)
narn.MCh6 <- size.select(MedCh6.sel, min.size = 350, max.size = 450, graph=TRUE, verbose=TRUE)
narn.MCh7 <- size.select(MedCh7.sel, min.size = 350, max.size = 450, graph=TRUE, verbose=TRUE)
narn.MCh8 <- size.select(MedCh8.sel, min.size = 350, max.size = 450, graph=TRUE, verbose=TRUE)


#Exp. reads per lane / (number of samples * number of loci)

Nova <- 1300000000
MiSeq <- 22000000

MiCov <- MiSeq/(40*17000)
NovaCov <- Nova/(40*17000)




