#!/usr/local/bin/Rscript

library(beadarray)
library(gplots)

load("data/BSData.quantile.RData")

dat <- read.csv("results/f_test.csv")

E <- exprs(BSData.quantile)

inds1 <- which(abs(dat[,"E11.5vE14.5"])>0.5)
inds2 <-which(dat[,"adj.P.Val"]<0.05)

inds1_2 <-intersect(inds1, inds2)

inds3 <- which(abs(dat[,"E14.5vP1"])>0.5)
inds4 <-which(dat[,"adj.P.Val"]<0.05)

inds3_4 <-intersect(inds3, inds4)

inds <-intersect(inds1_2, inds3_4)

ids <- dat[inds,"ID"]

ids <- as.character(ids)

filteredE <- E[ids,]

test<-apply(filteredE,1,function(x){any(x>10)})
filteredE2<-filteredE[test,]



symbol <- dat[inds,"symbol"]

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
gene.dist <- distancematrix(filteredE2, "euclid")
pretty<-distancematrix(filteredE, "cor")
postscript(file="results/pretty.ps", horizontal=FALSE)
dplot(gene.dist, 
	gene.hopach, 
	ord = "cluster", 
	main = "FACS Mice - Gene Distance Matrix", 
	showclusters = FALSE)
dev.off()



#now run hopach. K score relates to the number of permitted children
gene.hopach <- hopach(filteredE2, dmat=gene.dist, d="euclid", K=1)

#plot distance matrix
#postscript(file="results/distancematrix.ps", horizontal=FALSE)
dplot(gene.dist, 
	gene.hopach, 
	ord = "cluster", 
	main = "FACS Mice - Gene Distance Matrix", 
	showclusters = FALSE)
#dev.off()

gene.hopach$clust$k



#HOPACH clustering of arrays
array.hopach <- hopach(t(filteredE), d="euclid")
array.hopach$clust$k




#run Mfuzz to plot clusters
library(Mfuzz)
tmp_expr = new('ExpressionSet', exprs=filteredE)
cl = mfuzz(tmp_expr, c=5,m=3)
postscript(file="results/clusters.ps", horizontal=FALSE)
par(las=2)
#mfuzz.plot
matt.plot(tmp_expr,cl=cl,
	   mfrow=c(2,4),
	   new.window = FALSE,
	   time.labels=names,
	   min.mem=0.4,
           ymin = 6,
           ymax = 15,
           xlab = ""
                )
dev.off()

#extract genes and get membership scores
core <- acore(tmp_expr, cl, min.acore=0.3)
#core is a list which seems to be an arse
core1 <- core[[1]]
core2 <- core[[2]]
core3 <- core[[3]]
core4 <- core[[4]]
core5 <- core[[5]]

#annotate them somehow
library(biomaRt)
library(illuminaMousev2BeadID.db)
library(illuminaMousev2.db) 

ids = rownames(core5)
membership =core5[,2]

symbol <- mget(ids, illuminaMousev2SYMBOL, ifnotfound = NA)
sum(sapply(symbol, length) > 1)
symbol <- as.character(symbol)
length(which(symbol=="NA"))

ensembl = mget(ids, illuminaMousev2ENSEMBL, ifnotfound = NA)
length(which(is.na(ensembl)))

sum(sapply(ensembl, length) > 1)
crosshyb <- which(( sapply(ensembl, length) ) > 1)
length(crosshyb)
ensembl[crosshyb] <- NA
ensembl <- as.character(ensembl)
ensembl[ensembl=="NA"] = NA
length(which(is.na(ensembl)))

ensmart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")


filters <- "ensembl_gene_id"
values <- ensembl[!is.na(ensembl)]
attributes <- c("ensembl_gene_id", "chromosome_name","start_position", "end_position", "strand", "description")
ens.anno <- getBM(filters=filters, values=values, attributes=attributes, mart=ensmart)
rownames(ens.anno)<-ens.anno[,1]

anno <- data.frame(
                   ID = as.character(ids),
                   EnsemblID = ensembl,
                   symbol=symbol,
                   ens.anno[ensembl,],
                   stringsAsFactors=F
              )
rownames(anno) <- anno[,"ID"]

#add membership column

anno <- data.frame(anno,core5)

write.csv(anno, "results/Mfuzz/cluster5.csv")



#compute overlap and draw plots
O <- overlap(cl)
postscript(file="
Ptmp <- overlap.plot(cl,over= O,thres = 0.05)

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



