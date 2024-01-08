############################################################################################################################################################
# R script to filter data after process_radtags in Stacks and run analyses found in Dorsey et al 2023 Conservation genomics of Dioon holmgrenii (Zamiaceae) reveals
# a history of range expansion, fragmentation, and isolation
# of populations. Conservation Genetics
#
# Only D. holmgrenii populations used here. See script DholmPlusOGDiversityAnalysis_pub.R for analyses with D. planifolium as outgroup and SNP data.
#
############################################################################################################################################################


setwd("/Volumes/CatDisk/Dioon_RADseq_analysis/Dholm/Rworking/sixPops")

library(ade4)
library(adegenet)
library(pegas)
library(graph4lg)
library(ggplot2)
library(ggrepel)
library(plyr)
library(gridExtra)
library(hierfstat)
library(poppr)
library(diveRsity)
library(parallel)
library(boot)
library(data.table)
library(PopGenReport)
library(geosphere)
library(vegan)
library(ape)
library(radiator)

#  From https://rdrr.io/github/romunov/zvau/src/R/writeGenPop.R
writeGenPop <- function(gi, file.name, comment) {
  
  if (is.list(gi)) {
    # do all genind objects have the same number of loci?
    if (length(unique(sapply(gi, nLoc))) != 1) stop("Number of loci per individual genind object in a list is not equal for all.")
    gi.char <- gi[[1]]
    loc.names <- locNames(gi[[1]])
  } else {
    gi.char <- gi
    loc.names <- locNames(gi)
  }
  
  # Calculate the length of two alleles.
  lng <- as.character(na.omit(genind2df(gi.char)[, locNames(gi.char)[1]]))
  lng <- unique(nchar(lng))
  
  stopifnot(length(lng) == 1)
  
  cat(paste(comment, "\n"), file = file.name)
  cat(paste(paste(loc.names, collapse = ", "), "\n"), file = file.name, append = TRUE)
  
  if (is.list(gi)) {
    pop.names <- seq_len(length(gi))
  } else {
    pop.names <- popNames(gi)
  }
  
  for (i in pop.names) {
    cat("pop\n", file = file.name, append = TRUE)
    if (is.list(gi)) {
      intm <- gi[[i]]
      loc.names <- locNames(gi[[i]])
    } else {
      intm <- gi[pop(gi) == i, drop = FALSE]
    }
    ind.names <- indNames(intm)
    intm <- genind2df(intm, sep = "")
    intm[is.na(intm)] <- paste(rep("0", lng), collapse = "")
    out <- cbind(names = paste(ind.names, ",", sep = ""), intm[, loc.names])
    write.table(out, file = file.name, row.names = FALSE, col.names = FALSE, append = TRUE, quote = FALSE)
  }
  
  return(NULL)
}


##############
# Saved data after filtering below
# No outgroup

haps.7070 <- read.genepop(file="/Volumes/CatDisk/Dioon_RADseq_analysis/Dholm/Rworking/sixPops/haps.7070.gen")
levels(haps.7070@pop) <- c("Ixtayutla", "Jamiltepec", "Juchatengo", "Loxicha", "Rancho Limon", "Textitlan")

haps.7070.pop <- genind2genpop(haps.7070)
haps.7070.hier <- genind2hierfstat(haps.7070, pop=haps.7070$pop)

# Teita as OG
# haps.7070.og <- read.genepop(file="/Volumes/HD3/DebianShare/stacks/dholm/sixPlusOG/m3r5p7WL/populations.haps.genepop.gen")
# levels(haps.7070.og@pop) <- c("Ixtayutla", "Jamiltepec", "Juchatengo", "Loxicha", "Rancho Limon", "Textitlan", "Teita")

##############

######################################################################################################
# Filtering Data #####################################################################################

# Import data from Stacks

# haplotypes in genepop format
haps <- read.genepop(file="/Volumes/HD3/DebianShare/stacks/dholm/sixPop/m3r5p6BL/populations.haps.genepop.gen")
# haps@tab[1:10,1:10]

levels(haps@pop) <- c("Ixtayutla", "Jamiltepec", "Juchatengo", "Loxicha", "Rancho Limon", "Textitlan")
pops <- levels(haps@pop)


# haps$pop <- as.factor(sub("-.*","",haps$pop))
# haps$pop <- as.factor(sub("Dho","",haps$pop))

# raw haplotypes from stacks for fineRAD
haps.raw <- read.delim("/Volumes/HD3/DebianShare/stacks/dholm/sixPop/m3r5p6BL/populations.haplotypes.tsv", check.names = F)
haps.raw2 <- haps.raw[which(haps.raw[,1] %in% locNames(haps),arr.ind=T),c(1,3:length(haps.raw))] # remove filtered haps and count column

# check
# identical(as.character(haps.raw2[,1]),locNames(haps))
# identical(indNames(haps),as.character(names(haps.raw2)))

###########################
## use file to run bayescan then return to script

## Import bayescan R code

source("/Volumes/HD4/BayeScan2.1/R_functions/plot_R.r")

## plot to get outliers

# get indices of outliers
myOuts <- plot_bayescan("/Volumes/HD4/Dioon_RADseq_analysis/Dholm/bayescan/sixpops/sixpop.bayes_input_fst.txt",FDR = 0.1, pos=0.1)
#str(myOuts)
# myOuts$outliers
## get list of snps cut from vcf

snpList <- read.table("/Volumes/HD4/Dioon_RADseq_analysis/Dholm/bayescan/sixpops/snpList.txt", col.names=c("chr","pos"))

## get outlier indices

outHaps <- unique(snpList[myOuts$outliers,1]) # get locus names for outliers
# class(outHaps)

# num of loci under selection
length(outHaps)

## Filter loci under selection
## outHaps from BayeScan

remove.haps <- which(locNames(haps) %in% outHaps,arr.ind=T) # get indices in genind of outlier names

haps.neut <- haps[loc=-remove.haps] #remove outlier loci

neut.loci <- locNames(haps.neut)

#######
# raw for finRAD
# remove.raw <- which(haps.raw2[,1] %in% outHaps) # get indices in raw2 for outliers
raw.neut.loci <- which(haps.raw2[,1] %in% neut.loci)
haps.raw.neut <- haps.raw2[raw.neut.loci,]

# check
# identical(locNames(haps.neut),as.character(haps.raw.neut[,1]))


## Filtering missing data

# get proportion of missing data 
# across individuals

haps.typedInd <- propTyped(haps.neut, by="ind")

# snps.typedInd <- propTyped(snps.neut, by="ind")
# across loci

haps.typedLoc <- propTyped(haps.neut, by="loc")

# snps.typedLoc <- propTyped(snps.neut, by="loc")

# Filter by loci, keep with > 0.7 coverage

high.cover.haps <- which(haps.typedLoc >= 0.7)

haps.7.filt <- haps.neut[,loc=high.cover.haps]

# raw for finRAD
haps.raw.7.filt <- haps.raw.neut[high.cover.haps,]

#check
# identical(locNames(haps.7.filt),as.character(haps.raw.7.filt[,1]))

# Filter out Individuals with < 70% coverage, after filtering loci above

haps.7070 <- haps.7.filt[which(haps.typedInd > 0.7),,drop=T]

# drop monomorphic loci
polyLoc <- which(haps.7070@loc.n.all != 1)
# 123184 525169 
# 2279   2311 

haps.7070 <- haps.7070[, loc = polyLoc]

summary(haps.7070)$NA.perc
# 7.257294

writeGenPop(haps.7070, "haps.7070.gen", "Dholm filtered")

# object conversions

haps.7070.mat <- tab(haps.7070, freq=T, NA.method="mean")

# convert to genpop
haps.7070.pop <- genind2genpop(haps.7070)

# convert genpop tab to freq matrix
haps.7070.popmat <- tab(haps.7070.pop, freq=T, NA.method="mean")

# convert to hierfstat
haps.7070.hier <- genind2hierfstat(haps.7070)

# raw for finRAD
hi.cover.inds <- which(haps.typedInd > 0.7)

# Change first column name
names(haps.raw.7.filt)[1] <- "Locus"

keep.raw.inds <- c("Locus",names(hi.cover.inds)) # get names of individuals to keep - plus Locus column

haps.raw.7070 <- haps.raw.7.filt[,keep.raw.inds]
# dim(haps.raw.7070)
# haps.raw.7070[1:5,1:3]

haps.raw.7070.t <- transpose(haps.raw.7070[-1])
# dim(haps.raw.7070.t)
# haps.raw.7070.t[1:5,1:5]
setnames(haps.raw.7070.t, rownames(haps.raw.7070))
rownames(haps.raw.7070.t) <- names(haps.raw.7070)[-1]

haps.raw.7070.t <- data.frame(population=haps.7070@pop, haps.raw.7070.t)

######################################################################################################

######################################################################################################
##### Test for HWE per locus ####

#across all individuals

# exact test

haps.hwe.exact <- hw.test(haps.7070, B=100)

# summary(haps.hwe.exact)
num.out.HW.overall <- length(which(haps.hwe.exact[,4]< 0.05)) # count number of loci rejecting HWE

perc.out.overall <- num.out.HW.overall/2314

hist(haps.hwe.exact[,4], breaks = 20, main = "Dispribution of p-values from HW exact test across all individuals")

# split into separate pops
by.pop <- seppop(haps.7070, drop=T)

# HWE test for each population
by.pop.hwe.exact <- mclapply(by.pop, hw.test, B=1000, mc.cores = 20)

# any loci missing?
# lapply(names(by.pop.hwe.exact), function(x) length(which(is.na(by.pop.hwe.exact$x[,4])))) # no NaNs b/c all loci are present in each population

# count loci out of HW
count.outHW <- function(pop){
  length(which(by.pop.hwe.exact[[pop]][,4]<0.05))
}

num.loc.outHW <- lapply(names(by.pop.hwe.exact), count.outHW)
names(num.loc.outHW) <- names(by.pop.hwe.exact)

num.loc.outHW

# Percentage out per pop
perc.loc.outHW <- lapply(names(num.loc.outHW), function(x){num.loc.outHW[[x]]/nLoc(by.pop[[x]])})

names(perc.loc.outHW) <- popNames(haps.7070)

perc.loc.outHW

# test for significance per Waples 2014
# get critical values for 0.05 and 0.95
qbinom(c(0.05, 0.95), 2316, 0.05)
# 99 133

# plot hist of p-values - should be flat
par(mfrow=c(2,3))
lapply(names(by.pop.hwe.exact), function (x) hist(by.pop.hwe.exact[[x]][,4], main=x, xlab = "P-values"))


######################################################################################################
# Run DAPC ####

dat <- tab(haps.7070, freq=T, NA.method="mean")

# find clusters
grp <- find.clusters(dat, n.start=100, n.pca=150, scale=F, n.clust=6)  # retain 150 PCs
# 6 clusters chosen

#xvalidation to determine proper # pc to retain
crossvalTest <- xvalDapc(dat, grp$grp, n.pca.max = 150, center = T, scale = F, n.pca = NULL, n.rep = 30, xval.plot = T)

crossvalTest[6]
# $`Number of PCs Achieving Lowest MSE`
# [1] "90"

# dapc
haps.da <- dapc(dat, n.pca=90,  grp=grp$grp, n.da=4, var.loadings=T)

summary(haps.da)

compoplot(haps.da)

assignplot(haps.da)

haps.da$var

# equal contribution loading test
numAlleles <- length(haps.da$var.contr[,1])
thresh <- 1/numAlleles
# 0.0001207584

threshAlleles <- apply(haps.da$var.contr, 2, function (x) length(x[which(x>thresh)])) #  length(haps.da$var.contr[which(haps.da$var.contr[,x]>thresh),1])
contAlleles <- threshAlleles/numAlleles

# LD1       LD2       LD3       LD4 
# 0.2476754 0.2441734 0.2403091 0.2271465 

contrib <- loadingplot(haps.da$var.contr, axis=2, threshold = thresh, at = order(haps.da$var.contr[,2]))

par(mfrow = c(2,2))
barplot(sort(haps.da$var.contr[,3]))
abline(h=thresh)
######################################################################################################
# DAPC plotting ####

co <- as.data.frame(haps.da$ind.coord)
co$pop <- haps.7070$pop

cents <- as.data.frame(ddply(co, ~pop, summarize, mean1=mean(LD1), mean2=mean(LD2), mean3=mean(LD3), mean4=mean(LD4)))
#cents.grp <- as.data.frame(ddply(co, ~pop, summarize, mean1=mean(LD1), mean2=mean(LD2), mean3=mean(LD3)))

grp.coords <- as.data.frame(haps.da$grp.coord)

ld12 <- ggplot(co) + geom_point(aes(LD1, LD2, color=co$pop)) + geom_text_repel(data=cents, x=cents$mean1, y=cents$mean2, label = cents$pop) + theme(legend.position = "none") #+ geom_text_repel(data=grp.coords,  x=grp.coords$LD1, y=grp.coords$LD2, label=rownames(grp.coords)

ld34 <- ggplot(co) + geom_point(aes(LD3, LD4, color=co$pop)) + geom_text_repel(data=cents, x=cents$mean3, y=cents$mean4, label = cents$pop) + theme(legend.position = "none") #+ geom_text_repel(data=grp.coords,  x=grp.coords$LD3, y=grp.coords$LD4, label=rownames(grp.coords)

pdf("Dholm.dapc.20221004.pdf",8,12)

arrplot <- grid.arrange(ld12,ld34)

dev.off()

######################################################################################################
# Within pop diversity
######################################################################################################
# allelic richness - Using Metapop2 instead
##################
# arHier <- allelic.richness(haps.7070.hier)
# 
# arHier16 <- allelic.richness(haps.7070.hier, 16)
# summary(arHier16$Ar)
# 
# arHier$min.all
# 
# summary(arHier$Ar)
# 
# dev.off()
# par(mfrow=c(2,3))
# 
# lapply(1:6,function(x) {
#   hist(arHier$Ar[,x], main = "", xlab = "")
#   title(main = levels(haps.7070@pop)[x], xlab = "", cex = 3)
# })
# 
# pdf("ARBoxplots.pdf", 8,6)
# par(mfrow=c(2,3))
# lapply(1:6,function(x) {
#   boxplot(arHier$Ar[,x], main = levels(haps.7070@pop)[x], xlab = "Allelic Richness")
# })
# dev.off()
##################
# Private alleles with poppr
##################
priv.alls.dose <- poppr::private_alleles(haps.7070)

#total count of private allele copies per population
private.count <- rowSums(priv.alls.dose)

#count of alleles that are private to each population
priv.alls <- poppr::private_alleles(haps.7070,count.alleles = F)
priv.alls[1:6,1:10]
priv.alls.sums <- rowSums(priv.alls)

# number of distinct alleles present in each pop
distinctallelesInpops <- apply(haps.7070.pop@tab,1,function(l){length(l[l != 0])})

# percentage of alleles present in each pop that are private
perc.priv <- round((priv.alls.sums/distinctallelesInpops*100),2)

# number of individuals per population
Indcounts <- summary(haps.7070@pop)

# possible number of genotypes per population
all.poss <- Indcounts*nLoc(haps.7070)

# number of missing genotypes per individual
n.na.gt <- apply(is.na(haps.7070.hier),1,sum)

# number of missing genotypes per pop
pop.na.gt <- as.numeric(t(rowsum(n.na.gt,haps.7070.hier$pop)))

# number of present genotypes
actual.genos <- all.poss-pop.na.gt

# actual number of alleles (copies) per pop
total.alls <- actual.genos*2

# percentage of copies that are private alleles
perc.private.copies <- as.numeric(round((private.count/total.alls)*100,2))
# as proportion (freq)
freq.private.copies <- as.numeric(round((private.count/total.alls),4))
names(perc.private.copies) <- names(perc.priv)

# put into table
private.allele.df <- data.frame(names(private.count),Indcounts,actual.genos,total.alls,distinctallelesInpops,private.count,perc.private.copies,priv.alls.sums,perc.priv)

pa.df.sort <- arrange(private.allele.df,perc.priv)
##################
######################################################################################################
# hier.fstat for diversity
##################
basics <- basic.stats(haps.7070, diploid = T, digits = 3)

basics$overall
#    Ho    Hs    Ht   Dst   Htp  Dstp   Fst  Fstp   Fis  Dest 
# 0.283 0.293 0.332 0.039 0.340 0.047 0.118 0.139 0.034 0.067

apply(basics$Fis,2,median, na.rm=T)
#  Ixt1  Jamil    Juc    Lxa    Rli    Tex 
# 0.000 -0.024 -0.026  0.000 -0.030 -0.006 

apply(basics$Fis,2,mean, na.rm=T)

apply(basics$Ho,2,median, na.rm=T)
# Ixtayutla   Jamiltepec   Juchatengo      Loxicha Rancho Limon    Textitlan 
# 0.250        0.265        0.286        0.231        0.278        0.269 

apply(basics$Hs,2,median, na.rm=T)
# Ixtayutla   Jamiltepec   Juchatengo      Loxicha Rancho Limon    Textitlan 
# 0.2810       0.2830       0.3135       0.2540       0.2880       0.3060 

##################
# bootstraps
##################

medianFunc <- function(data, indices) {
  data <- as.matrix(data)
  dt <- data[indices,]
  med <- apply(dt,2,median, na.rm = T)
  return(med)
}


FisBootHier <- boot(basics$Fis, medianFunc, R=10000, parallel = "multicore")
HoBootHier <- boot(basics$Ho, medianFunc, R=10000, parallel = "multicore")
HsBootHier <- boot(basics$Hs, medianFunc, R=10000, parallel = "multicore")
ArBootHier <- boot(arHier$Ar, medianFunc, R=10000, parallel = "multicore")

FhierCIlist <- list()
for (i in 1:6) {
  FhierCIlist[[i]] <- boot.ci(FisBootHier, index = i)
}
allFisCI <- do.call(rbind, (lapply(1:6, function (x) FhierCIlist[[x]]$normal)))

HohierCIlist <- list()
for (i in 1:6) {
  HohierCIlist[[i]] <- boot.ci(HoBootHier, index = i)
}
allHoCI <- do.call(rbind,(lapply(1:6, function (x) HohierCIlist[[x]]$normal)))

HshierCIlist <- list()
for (i in 1:6) {
  HshierCIlist[[i]] <- boot.ci(HsBootHier, index = i)
}
allHsCI <- do.call(rbind,(lapply(1:6, function (x) HshierCIlist[[x]]$normal)))

ArCIlist <- list()
for (i in 1:6) {
  ArCIlist[[i]] <- boot.ci(ArBootHier, index = i)
}
allArCI <- do.call(rbind, (lapply(1:6, function (x) ArCIlist[[x]]$normal)))




# allFisCI <- lapply(1:6, function (x) medBootCI_List[[x]]$normal)
# allFisCI <- do.call(rbind,allFisCI)


######################################################################################################
# Diversity
##################
# !! Change working directory for diveRsity!!
getwd()
setwd("./diveRsity/")

# diveRsity for Fstats
# population differentiation

# bs across loci
haps.diff.bsLoc <- diffCalc("/Volumes/HD4/Dioon_RADseq_analysis/Dholm/Rworking/sixPops/haps.7070.gen", outfile = "haps.diffCalc.20220422", fst=T, pairwise=T, bs_pairwise=T, boots=1000, para= TRUE, ci_type = "loci", bs_locus = TRUE)

stdStats <- read.delim("haps.diffCalc.20220422-[diffCalc]/std_stats.txt")

haps.diff.bsLoc$pairwise$Gst
haps.diff.bsLoc$pairwise$D
haps.diff.bsLoc$std_stats[2315,]
######################################################################################################

setwd("/Volumes/CatDisk/Dioon_RADseq_analysis/Dholm/Rworking/sixPops")

######################################################################################################
#AMOVA
##################
strata(haps.7070) <- data.frame(haps.7070@pop)

nameStrata(haps.7070) <- ~population

amova.results.loci.ade <- poppr.amova(haps.7070, hier = ~population, missing = "loci", within=T)

amova.test <- randtest(amova.results.loci.ade, nrepet = 999)

plot(amova.test)

print(amova.results.loci.ade)

print(amova.test)

######################################################################################################
# Test for Isolation by Distance among populations of D. holmgrenii
##################

# Use Cavalli-Sforza and Edwards Chord distance using heirfstat  (Takezaki & Nei, 1996)
Dgen.chord <- genet.dist(haps.7070)

# calculate great circle distances between all pop pairs
# population coordinates
latlong <- read.delim("/Volumes/HD4/Dioon_RADseq_analysis/Dholm/Rworking/sixPops/latlong.txt", header=T)
# switch lat/long
longlat <- latlong[,2:1]

# add to genpop object
haps.7070.pop$other$xy <- longlat

#get distances
geodist <- distm(longlat, fun = distHaversine)

# set row/col names again
row.names(geodist) <- row.names(longlat)
colnames(geodist) <- rownames(longlat)
# matrix
geodist

#convert to dist object
geoDistObj <- as.dist(geodist) 

# Dgeo <- dist(haps.7070.pop$other$xy)
# # ibd <- mantel.randtest(Dgen,Dgeo,nrepet = 9999)
# # ibd
# 

plot(geoDistObj,Dgen.chord)
abline(lm(Dgen.chord~geoDistObj))

summary(lm(Dgen.chord~geoDistObj))

ibd.chord <- mantel.randtest(Dgen.chord, geoDistObj, nrepet = 9999)

ibd.chord

plot(ibd.chord, main="Distribution of Permutation Results and Observed Mantel Test Value")

#######################################################################################################

######################################################################################################
# plot loci distributions

# Ho
pdf("Ho_freq.pdf",8,6)

par(mfrow=c(2,3), oma = c(2,2,0,0))

lapply(1:6,function(x) {
  hist(basics$Ho[,x], main = levels(haps.7070$pop)[x], cex.main = 2, xlab = "", ylab = "", breaks = (0:20)/20)
  #abline(v=mean(basics$Ho[,x], na.rm = T))
})
mtext("Ho", side = 1, outer = T, cex = 1.3)
mtext("Frequency", side = 2, outer = T, cex = 1.3)

dev.off()

#Hs
pdf("Hs_freq.pdf",8,6)
par(mfrow=c(2,3), oma = c(2,2,0,0))

lapply(1:6,function(x) {
  hist(basics$Hs[,x], main = levels(haps.7070$pop)[x], cex.main = 2, xlab = "", ylab = "", breaks = (0:20)/20)
})
mtext("Hs", side = 1, outer = T, cex = 1.3)
mtext("Frequency", side = 2, outer = T, cex = 1.3)

dev.off()

#F
pdf("F_freq.pdf",8,6)

par(mfrow=c(2,3), oma = c(2,2,0,0))
lapply(1:6,function(x) {
  hist(basics$Fis[,x], main = levels(haps.7070$pop)[x], cex.main = 2, xlab = "", ylab = "")
  # abline(v=mean(basics$Fis[,x], na.rm = T))
})
mtext("F", side = 1, outer = T, cex = 1.3)
mtext("Frequency", side = 2, outer = T, cex = 1.3)

dev.off()


#AR

pdf("AR.pdf",8,6)
par(mfrow=c(2,3), oma = c(2,2,0,0))
lapply(1:6,function(x) {
  hist(arHier$Ar[,x], main = levels(haps.7070$pop)[x], cex.main = 2, xlab = "", ylab = "")
})
mtext("Allelic Richness", side = 1, outer = T, cex = 1.3)
mtext("Frequency", side = 2, outer = T, cex = 1.3)

dev.off()


# Globals

# stdStats <- head(stdStats, -1)
# names(stdStats)
# range(stdStats$Fst, na.rm = T)

pdf("globalFstatDists.pdf",8,3)

par(mfrow=c(1,3), oma = c(0,2,0,0))
fstHist <- hist(stdStats[,7], main = "", xlab = "Fst", ylab = "", cex.lab = 1.5)
fisHist <- hist(stdStats[,6], main = "", xlab = "Fis", ylab = "", cex.lab = 1.5)
fitHist <- hist(stdStats[,8], main = "", xlab = "Fit", ylab = "", cex.lab = 1.5)

mtext("Frequency", side = 2, outer = T, cex = 1)

dev.off()










