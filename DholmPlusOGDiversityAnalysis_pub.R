############################################################################################################################################################
# R script to filter data after process_radtags in Stacks and run analyses found in Dorsey et al 2023 Conservation genomics of Dioon holmgrenii (Zamiaceae) reveals
# a history of range expansion, fragmentation, and isolation
# of populations. Conservation Genetics
#
# SNP data from D. holmgrenii populations and D. planifolium outgroup population used here. See script Dholm_genDiversityAnalysis_pub.R for analyses 
# with of haplotype data from D. holmgrenii only.
#
############################################################################################################################################################



setwd("/Volumes/CatDisk/Dioon_RADseq_analysis/Dholm/Rworking/sixPlusOG")


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
library(HWxtest)
library(PopGenReport)
library(geosphere)
library(vegan)
library(ape)
library(radiator)
library(rangeExpansion)
library(sp)
library(snpStats)
library(rworldmap)


######


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


# Read in Stacks file with Teita as OG
snps <- read.genepop(file="/Volumes/HD3/DebianShare/stacks/dholm/sixPlusOG/m3r5p7WL/singleSNP/populations.snps.genepop.gen")
pops <- c("Ixtayutla", "Jamiltepec", "Juchatengo", "Loxicha", "Rancho Limon", "Textitlan", "Teita")
levels(snps@pop) <- pops

# Convert to VCF
genomic_converter(snps, output="vcf", filename = "snps.vcf")

###########################
# Identify potential loci under selection

## use vcf from stacks file to run bayescan then return here

## Import bayescan R code

source("/Volumes/HD4/BayeScan2.1/R_functions/plot_R.r")

## plot to get outliers

# get indices of outliers
myOuts <- plot_bayescan("/Volumes/HD4/Dioon_RADseq_analysis/Dholm/bayescan/sixPlusOG/run1/sixOG_bayes_input_fst.txt",FDR = 0.1, pos=0.1)

snpList <- read.table("/Volumes/HD4/Dioon_RADseq_analysis/Dholm/bayescan/sixPlusOG/run1/snpkey.txt", col.names=c("index","pos"))


# num of loci under selection
myOuts$nb_outliers

## Filter loci under selection
## from BayeScan

snps.neut <- snps[,loc=-myOuts$outliers] #remove outlier loci

neut.loci <- locNames(snps.neut)

###########################
## Filtering missing data

# get proportion of missing data 
# across individuals
snps.typedInd <- propTyped(snps.neut, by="ind")

# across loci
snps.typedLoc <- propTyped(snps.neut, by="loc")

# Filter by loci, keep with > 0.7 coverage
high.cover.snps <- which(snps.typedLoc >= 0.7)

snps.7.filt <- snps.neut[,loc=high.cover.snps]

# Filter out Individuals with < 70% coverage, after filtering loci above
snps.7070 <- snps.7.filt[which(snps.typedInd > 0.7),,drop=T]

# drop monomorphic loci
polyLoc <- which(snps.7070@loc.n.all != 1)

length(polyLoc)
polyLoc[1:10]

snps.7070.og <- snps.7070[, loc = polyLoc]

summary(snps.7070.og)$NA.perc

writeGenPop(snps.7070.og, "snps.7070.og.genpop", "D holm plus OG snps")
genomic_converter(snps.7070.og, output="vcf", filename = "snps.7070.og.vcf")


##################################
# Randge expansion test

# set files and parameters
snp.file <- "expansion/snps.7070.og.snapp"
coord.file <- "expansion/coordsFile.csv"
ploidy <- 2
region <- list(NULL)

# load data
raw.data <- load.data.snapp(snp.file, coord.file, sep=',', ploidy = ploidy)

# calculate pop level data and pairwise statistics
ex.pop <- make.pop(raw.data, ploidy)

system.time( {psi <- get.all.psi(ex.pop)} )

rownames(psi) <- pops[1:6] #levels(snps.7070.og$pop)[1:6]
colnames(psi) <- pops[1:6] #levels(snps.7070.og$pop)[1:6]

psi

ex.results <- run.regions(region=region, pop=ex.pop, psi=psi)

sum.ex.res <- summary(ex.results)

plot.new()

plot(ex.results)

######################################
# Permutation analysis
######################################

# make an array of randomly assigned pop arrangements at each locus
# used in permute.psi functions below
create.perm.array <- function(nperm, popsObj){
  n.snps <- ncol(popsObj$data)
  n.pops <- popsObj$n
  
  p.array <- array(0, dim = c(n.pops,n.snps,nperm)) # creat an array to hold permutations CHANGE THE DIMS!!!!

  for (p in 1:nperm) {                #for each permutation randomly sample pop order
    for (a in 1:n.snps){
      p.array[,a,p] <- sample(n.pops)
    }
  }
  return(p.array)
}


######################
# parallel

permute.psi.par <- function(nperm,popsObj) {
  n.snps <- ncol(popsObj$data)
  n.pops <- popsObj$n

  index <- p.array[,,nperm] # get one permutation as matrix of pop/locus indices
  d <- lapply(1:n.snps, function(x){popsObj$data[,x][index[,x]]})
  s <- lapply(1:n.snps, function(x){popsObj$ss[,x][index[,x]]})
  p <- list(data=as.data.frame(d), ss=as.data.frame(s), coords=popsObj$coors, n=popsObj$n)
  res <- get.all.psi(p)

  return(res)

}

# Parallel permutations
perms <- seq(1,10000)

numCores <- 20

p.array <- create.perm.array(length(perms), ex.pop)

system.time({
  para.results <- mclapply(perms, permute.psi.par, mc.cores = numCores, popsObj = ex.pop)
})

res.array <- array(unlist(para.results), dim = c(6,6,length(perms)))

#############################
# calculate p-values for psi estimates from permutation distribution of psi values
# uses ecdf function to calculate cumulative distribution function
get.psi.pvals <- function(permObj, psi.table) {
  n.pops <- dim(psi.table)[1]
  psi.pvals <- matrix(NA,nrow = n.pops, ncol = n.pops)
  
  for (r in 1:(n.pops)){
    for (c in 1:n.pops){
      t <- psi.table[r,c]
      P <- ecdf(permObj[r,c,])
      p <- P(t)
      if(t > 0){
        p <- 1-p
      }
      psi.pvals[r,c] <- p
      # p.res <- rbind(p.res, c(r,c,p))
    }
  }
  return(psi.pvals)
}

# plot permutation disctibutions, actual psi values, and associated p-values
plot.psi.pvals <- function(pvalTable, permObj, psi.table){
  n.pops <- dim(psi.table)[1]
  par(mfrow = c(3,2))
  
  for (r in 1:(n.pops-1)){
    for (c in (r+1):n.pops){
      plot(density(permObj[r,c,]), main= paste(row.names(psi.table)[r],"vs", row.names(psi.table)[c]), xlab = "Psi")
      abline(v=psi.table[r,c])
      p <- pvalTable[r,c]
      text(0, 10, p)
    }
  }  
}

pvals.7070 <- get.psi.pvals(res.array, psi)
rownames(pvals.7070) <- pops[1:6] #levels(snps.7070.og$pop)[1:6]
colnames(pvals.7070) <- pops[1:6] #levels(snps.7070.og$pop)[1:6]

pvals.7070

plot.psi.pvals(pvals.7070, res.array, psi)

# show significant psi values
sig.psi <- psi

sig.psi[which(pvals.7070 > 0.05)] <- NA
sig.psi

###############################################
# order expansion of populations
psi.sums <- c(rep(0,6))
names(psi.sums) <- pops[1:6]

for (r in 1:dim(psi)[1]){
  psi.sums[r] <- sum(psi[r,which(pvals.7070[r,]<0.05)])
}

psi.sums
psi.sums.ord <- psi.sums[order(psi.sums)]
psi.sums.ord
#########################################3





