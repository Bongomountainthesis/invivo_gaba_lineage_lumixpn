#!/usr/local/bin/Rscript

library(beadarray)
library(gplots)

load("data/BSData.quantile.RData")

limma11_14 <- read.csv("results/e11_e14.csv")
limma14_p1 <- read.csv("results/e14_p1.csv")

E <- exprs(BSData.quantile)

inds1 <- which(abs(limma11_14[,"logFC"])>1)
inds2 <-which(limma11_14[,"adj.P.Val"]<0.05)

inds1_2 <-intersect(inds1, inds2)

inds3 <- which(abs(limma14_p1[,"logFC"])>1)
inds4 <-which(limma14_p1[,"adj.P.Val"]<0.05)

inds3_4 <-intersect(inds3, inds4)

inds <-intersect(inds1_2, inds3_4)

ids <- limma11_14[inds,"ID"]

ids <- as.character(ids)

filteredE <- E[ids,]

symbol <- limma11_14[inds,"symbol"]

names <- c("E11.5a","E11.5b","E14.5a","E14.5b","E14.5c","E14.5d","P1a","P1b")

postscript(file="results/heatmap.ps", horizontal=FALSE)
heatmap.2(filteredE[,],
		Colv=NA,
		col=greenred(75), 
		scale="none",
		key=TRUE,
		keysize=0.75,
		symkey=FALSE,
		density.info="none",
		trace="none", 
		labRow=NA,
		labCol=NA,
		cexRow=0.75,
	)
dev.off()

####Mfuzz#####
#try hopach to define numbers of clusters
library(hopach)

#compute the distance matrix first
gene.dist <- distancematrix(filteredE, "cosangle")

#now run hopach
gene.hopach <- hopach(filteredE, dmat=gene.dist)

#plot distance matrix
dplot(gene.dist, 
	gene.hopach, 
	ord = "cluster", 
	main = "FACS Mice - Gene Distance Matrix", 
	showclusters = FALSE)

#HOPACH clustering of arrays
array.hopach <- hopach(t(filteredE), d="euclid")
array.hopach$clust$k

#run Mfuzz to plot clusters
library(Mfuzz)
tmp_expr = new('ExpressionSet', exprs=filteredE)
cl = mfuzz(tmp_expr, c=10,m=3)
postscript(file="results/clusters.ps", horizontal=FALSE)
par(las=2)
#mfuzz.plot
matt.plot(tmp_expr,cl=cl,
	   mfrow=c(2,4),
	   new.window = FALSE,
	   time.labels=names,
	   min.mem=0.2,
           ymin = 4,
           ymax = 18,
           xlab = ""
                )
dev.off()

cluster=1
cl[[4]][,cluster]





matt.plot<-function (eset, cl, mfrow = c(1, 1), colo, min.mem = 0, time.labels, new.window = TRUE, ymin=-999, ymax=-999, xlab="Time", ylab="Expression Chan ges")
{
    clusterindex <- cl[[3]]
    memship <- cl[[4]]
    memship[memship < min.mem] <- -1
    colorindex <- integer(dim(exprs(eset))[[1]])
    if (missing(colo)) {
        colo <- c("#FF8F00", "#FFA700", "#FFBF00", "#FFD700",
            "#FFEF00", "#F7FF00", "#DFFF00", "#C7FF00", "#AFFF00",
            "#97FF00", "#80FF00", "#68FF00", "#50FF00", "#38FF00",
            "#20FF00", "#08FF00", "#00FF10", "#00FF28", "#00FF40",
            "#00FF58", "#00FF70", "#00FF87", "#00FF9F", "#00FFB7",
            "#00FFCF", "#00FFE7", "#00FFFF", "#00E7FF", "#00CFFF",
            "#00B7FF", "#009FFF", "#0087FF", "#0070FF", "#0058FF",
            "#0040FF", "#0028FF", "#0010FF", "#0800FF", "#2000FF",
            "#3800FF", "#5000FF", "#6800FF", "#8000FF", "#9700FF",
            "#AF00FF", "#C700FF", "#DF00FF", "#F700FF", "#FF00EF",
            "#FF00D7", "#FF00BF", "#FF00A7", "#FF008F", "#FF0078",
            "#FF0060", "#FF0048", "#FF0030", "#FF0018")
    }
    colorseq <- seq(0, 1, length = length(colo))
    for (j in 1:max(clusterindex)) {
        tmp <- exprs(eset)[clusterindex == j, ]
        tmpmem <- memship[clusterindex == j, j]
        if (((j - 1)%%(mfrow[1] * mfrow[2])) == 0) {
            if (new.window)
                X11()
            par(mfrow = mfrow)
 
            if (sum(clusterindex == j) == 0) {
                if(ymin == -999){ymin <- -1}
                if(ymax == -999) {ymax <- +1}
            }
            else {
                if(ymin == -999) {ymin <- min(tmp)}
                if(ymax == -999) {ymax <- max(tmp)}
            }
            plot.default(x = NA, xlim = c(1, dim(exprs(eset))[[2]]),
                ylim = c(ymin, ymax), xlab = xlab, ylab = ylab,
                main = paste("Cluster", j), axes = FALSE)
            if (missing(time.labels)) {
                axis(1, 1:dim(exprs(eset))[[2]], c(1:dim(exprs(eset))[[2]]))
                axis(2)
            }
           else {
                axis(1, 1:dim(exprs(eset))[[2]], time.labels)
                axis(2)
            }
        }
        else {
            if (sum(clusterindex == j) == 0) {
                if(ymin == -999){ymin <- -1}
                if(ymax == -999) {ymax <- +1}
            }
            else {
                if(ymin == -999) {ymin <- min(tmp)}
                if(ymax == -999) {ymax <- max(tmp)}
            }
            plot.default(x = NA, xlim = c(1, dim(exprs(eset))[[2]]),
                ylim = c(ymin, ymax), xlab = xlab, ylab = ylab,
                main = paste("Cluster", j), axes = FALSE)
            if (missing(time.labels)) {
                axis(1, 1:dim(exprs(eset))[[2]], c(1:dim(exprs(eset))[[2]]))
                axis(2)
            }
            else {
                axis(1, 1:dim(exprs(eset))[[2]], time.labels)
                axis(2)
            }
        }
        if (!(sum(clusterindex == j) == 0)) {
            for (jj in 1:(length(colorseq) - 1)) {
                tmpcol <- (tmpmem >= colorseq[jj] & tmpmem <=
                  colorseq[jj + 1])
                if (sum(tmpcol) > 0) {
                  tmpind <- which(tmpcol)
                  for (k in 1:length(tmpind)) {
                    lines(tmp[tmpind[k], ], col = colo[jj])
                  }
                }
            }
        }
    }
}



