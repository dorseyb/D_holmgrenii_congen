############################################################################################################################################################
# R script to run analyses on morphological data found in Dorsey et al 2023 Conservation genomics of Dioon holmgrenii (Zamiaceae) reveals
# a history of range expansion, fragmentation, and isolation
# of populations. Conservation Genetics
# 
############################################################################################################################################################


setwd("/Users/bdorsey/Documents/Rworking/DholMorph")

library(adegenet)
library(vegan)
library(pairwiseAdonis)
library(pca3d)
library(car)
library(ggplot2)
library(ggrepel)
library(plyr)
library(dplyr)
library(gridExtra)
library(rgl)
library(lattice)
library(MASS)
library(geosphere)
library(usedist)

##################################
# Get data
data <- as.data.frame(read.csv("/Users/bdorsey/Documents/Rworking/DholMorph/Dhol_data_compiled.csv"))
row.names(data) <- data$specimen
levels(data$pop) <- c("Ixtayutla","Jamiltapec","Juchatengo","Rancho Limon","Loxicha","Textitlan", "Zaragosa")

data <- data[order(data$specimen),]

data6 <- data[data$pop %in% c("Ixtayutla","Jamiltapec","Juchatengo","Rancho Limon","Loxicha","Textitlan"),]
data6$pop <- droplevels(data6$pop)

# summary(data6)
# names(data6)

data6num <- data6[5:17]
data6num <- data6num[-c(8)]
# length(data6num)
# names(data6num)

# look for outliers

# mahalanobis distances
m_dist <- mahalanobis(data6num, colMeans(data6num), cov(data6num))

# pdf("outlier.boxplot.pdf",8,10)
# boxplot(m_dist)
# title("Mahalanobis Distances of Raw Data")
# dev.off()

# get critical value from chi-squared dist
chi.crit <- qchisq(c(0.95),12)
# 21.02607

data6.noOuts <- data6num[which(m_dist < chi.crit),]
noOuts.pops <- data6$pop[which(m_dist < chi.crit)]

# of outliers
# length(m_dist[which(m_dist > 21.026)])

# summary(data6.noOuts)



##################################

# DAPC


###############
# infer groups#
###############

crossvalTest <- xvalDapc(data6.noOuts, noOuts.pops, n.pca.max = 12, center = T, scale = T, n.pca = NULL, n.rep = 100, xval.plot = T)

crossvalTest[6]
# $`Number of PCs Achieving Lowest MSE`
# [1] "11"

grp_no <- find.clusters(data6.noOuts, method="kmeans", stat='BIC', max.n.clust=20, center = T, scale=T, n.start = 100, n.pca = 11, n.clust = 8)  #choose 8

da_no <- dapc(data6.noOuts, grp = grp_no$grp, center = T, scale = T, n.pca = 11, n.da = 2) #choose 2 df

#####
#multiple runs show different groupings
par(mfrow=c(2,5), oma=c(2,2,1,1))
for (k in 1:10) {
  grp_no <- find.clusters(data6.noOuts, method="kmeans", stat='BIC', max.n.clust=20, center = T, scale=T, n.start = 100, n.pca = 11, n.clust = 8)  #choose 8
  
  da_no <- dapc(data6.noOuts, grp = grp_no$grp, center = T, scale = T, n.pca = 11, n.da = 2) #choose 2 df
  
  
  scatter(da_no, 1,2, scree.da = F, scree.pca = F, bg ="white", pch=15:20, cstar=F, cellipse = F, col = safe_colorblind_palette_point)
  
  mtext("LD1", side = 1, outer = T)
  mtext("LD2", side = 2, outer = T)
  
  for (i in 1:9) {
    p <- levels(grp_no$grp)[i]
    f <- which(grp_no$grp==p)
    Plot_ConvexHull(da_no$ind.coord[f,1],da_no$ind.coord[f,2], lcolor = safe_colorblind_palette_hull[i])
  }
  
}

#####

numTraits <- length(da_no$var.contr[,1])
thresh <- 1/numTraits
threshTraits <- apply(da_no$var.contr, 2, function (x) x[which(x>thresh)])
threshTraitsCount <- apply(da_no$var.contr, 2, function (x) length(x[which(x>thresh)])) #  length(haps.da$var.contr[which(haps.da$var.contr[,x]>thresh),1])
contTraits <- threshTraits/numTraits


##############
# using pops #
##############

da_no_pops <- dapc(data6.noOuts, grp = noOuts.pops, center = T, scale = T, n.pca = 11, n.da = 2)

summary(da_no_pops)

table(da_no_pops$assign, noOuts.pops)

da_no_pops$var.contr


popthreshTraits <- apply(da_no_pops$var.contr, 2, function (x) x[which(x>thresh)])
popthreshTraitsCount <- apply(da_no_pops$var.contr, 2, function (x) length(x[which(x>thresh)]))


#################
# PLot DAPC

# variable names for legends

morphVars <- c("Trunk Diameter", "Leaflet LxW", "Leaflet Width", "Leaflet Length", "Leaflet Base", "Num. Leaflets", "Petiole Width", "Petiole Length", "Leaf Width", "Leaf Length", "Trunk Height", "Num. Leaves")

morphVars <- rev(morphVars) #fix wrong order

#Function
Plot_ConvexHull<-function(xcoord, ycoord, lcolor){
  hpts <- chull(x = xcoord, y = ycoord)
  # hpts <- c(hpts, hpts[1])
  # lines(xcoord[hpts], ycoord[hpts], col = lcolor)
  polygon(xcoord[hpts], ycoord[hpts], col = lcolor)
}

# set colors
rgbVal <- t(col2rgb(c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")))

safe_colorblind_palette_hull <- rgb(rgbVal, alpha = 50, maxColorValue = 255)

safe_colorblind_palette_point <- rgb(rgbVal, maxColorValue = 255)

# plot k-means DAPC results

pdf("kmeans.morph.dapc1-2.pdf",10,10)
par(oma=c(2,2,1,1), mar = c(0,0,0,0))

scatplot_grp <- scatter(da_no, 1,2, scree.da = F, scree.pca = F, bg ="white", pch=15:20, cstar=F, cellipse = F, col = safe_colorblind_palette_point)

mtext("LD1", side = 1, outer = T, cex = 1.5, line = 0.5)
mtext("LD2", side = 2, outer = T, cex = 1.5)

for (i in 1:9) {
  p <- levels(grp_no$grp)[i]
  f <- which(grp_no$grp==p)
  Plot_ConvexHull(da_no$ind.coord[f,1],da_no$ind.coord[f,2], lcolor = safe_colorblind_palette_hull[i])
}
dev.off()


# plot population DAPC

pdf("pop.morph.dapc.pdf",10,10)
par(oma=c(2,2,1,1), mar = c(0,0,0,0))

scatplot <- scatter(da_no_pops,1,2, scree.da=F, scree.pca = F, bg="white", pch=15:20, cstar=F, cellipse = F, legend = F, clab=F,  txt.leg=c("Ixtayutla","Jamiltapec","Juchatengo","Rancho Limon","Loxicha","Textitlan"), xlim = range(da_no_pops$ind.coord[,1])+c(-1,1), col = safe_colorblind_palette_point)#range(da_no_pops$ind.coord[,2])+c(-0.25,0.25)
legend("topright", legend = c("Ixtayutla","Jamiltapec","Juchatengo","Rancho Limon","Loxicha","Textitlan"), pch=15:20, col = safe_colorblind_palette_point, cex = 1.5, pt.cex = 2)
mtext("LD1", side = 1, outer = T, cex = 1.5, line = 0.5)
mtext("LD2", side = 2, outer = T, cex = 1.5)

for (i in 1:6) {
  p <- levels(noOuts.pops)[i]
  f <- which(noOuts.pops==p)
  Plot_ConvexHull(da_no_pops$ind.coord[f,1],da_no_pops$ind.coord[f,2], lcolor = safe_colorblind_palette_hull[i])
}

points(da_no_pops$grp.coord, pch=15:20, col=safe_colorblind_palette_point, cex=3)

dev.off()

# plot assignment PP

View(assignplot)

# Use image.plot from library(fields) instead of image, adds legend scale
newAssignplot <- function (x, only.grp = NULL, subset = NULL, new.pred = NULL, 
          cex.lab = 0.75, pch = 3) 
{
  require(fields)
  if (!inherits(x, "dapc")) 
    stop("x is not a dapc object")
  if (!is.null(new.pred)) {
    n.new <- length(new.pred$assign)
    x$grp <- c(as.character(x$grp), rep("unknown", n.new))
    x$assign <- c(as.character(x$assign), as.character(new.pred$assign))
    x$posterior <- rbind(x$posterior, new.pred$posterior)
  }
  if (!is.null(only.grp)) {
    only.grp <- as.character(only.grp)
    ori.grp <- as.character(x$grp)
    x$grp <- x$grp[only.grp == ori.grp]
    x$assign <- x$assign[only.grp == ori.grp]
    x$posterior <- x$posterior[only.grp == ori.grp, , drop = FALSE]
  }
  else if (!is.null(subset)) {
    x$grp <- x$grp[subset]
    x$assign <- x$assign[subset]
    x$posterior <- x$posterior[subset, , drop = FALSE]
  }
  n.grp <- ncol(x$posterior)
  n.ind <- nrow(x$posterior)
  Z <- t(x$posterior)
  Z <- Z[, ncol(Z):1, drop = FALSE]
  image.plot(x = 1:n.grp, y = seq(0.5, by = 1, le = n.ind), Z, 
        col = rev(heat.colors(100)), yaxt = "n", ylab = "", 
        xaxt = "n", xlab = "Clusters")
  axis(side = 1, at = 1:n.grp, tick = FALSE, labels = colnames(x$posterior))
  # axis(side = 2, at = seq(0.5, by = 1, le = n.ind), labels = rev(rownames(x$posterior)), 
       # las = 1, cex.axis = cex.lab)
  abline(h = 1:n.ind, col = "lightgrey")
  abline(v = seq(0.5, by = 1, le = n.grp))
  box()
  newGrp <- colnames(x$posterior)
  x.real.coord <- rev(match(x$grp, newGrp))
  y.real.coord <- seq(0.5, by = 1, le = n.ind)
  points(x.real.coord, y.real.coord, col = "deepskyblue2", 
         pch = pch)
  return(invisible(match.call()))
}

# get group midpoints on y-axis for axis labels
popNums <- rev(table(noOuts.pops))
popSums <- cumsum(popNums)
popMids <- popSums[1:length(popSums)] - popNums[1:length(popNums)]/2

# k-means groups
pdf("morph.kmean.assignPlot.pdf",8,8)
newAssignplot(da_no)
abline(h = popSums, col = "black", lwd=3)
axis(2, at = popMids, rev(levels(noOuts.pops)))
dev.off()

# a priori pops
pdf("morph.pops.assignPlot.pdf",10,10)
newAssignplot(da_no_pops)
abline(h = popSums, col = "black", lwd=3)
axis(2, at = popMids, rev(levels(noOuts.pops)))
dev.off()

#############################

# Test for sig diff between pops in morpho space

# Distance in original morpho space
morph_dist <- vegdist(data6.noOuts[5:12], "euclidean")

# diff in disperson can mimic diff in means, test for this
anova(betadisper(morph_dist, noOuts.pops))

# permanova
morph_adonis <- adonis(morph_dist ~ noOuts.pops, permutations = 9999)

# pairwise permanova
pw_morph_adonis <- pairwise.adonis(morph_dist, noOuts.pops, perm = 9999, p.adjust.m = "holm")

