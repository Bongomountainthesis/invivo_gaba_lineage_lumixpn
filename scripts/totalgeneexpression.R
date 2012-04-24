
library(beadarray)
library(limma)
library(biomaRt)
library(illuminaMousev2BeadID.db)



BSData <- get(load("results/BSData.quantile.RData"))
E <- exprs(BSData)

#E.ordered<-E[order(abs(E[,2]), decreasing=TRUE),]
#E.nodup<-E.ordered[!duplicated(E.ordered[,1]),]

ids = rownames(E)
symbol <- mget(ids, illuminaMousev2BeadIDSYMBOL, ifnotfound = NA)
sum(sapply(symbol, length) > 1)
symbol <- as.character(symbol)

mean.ns<-rowMeans(E[,1:5])
mean.neurons<-rowMeans(E[,6:8])

dat<-cbind(rownames=symbol, mean.ns, mean.neurons)

dat.ordered<-dat[order(dat[,2], decreasing=TRUE),]

dat.nodup<-dat.ordered[!duplicated(dat.ordered[,1]),]

genesexpressedinneurons<-dat.nodup[dat.nodup[,3]>=1,]
