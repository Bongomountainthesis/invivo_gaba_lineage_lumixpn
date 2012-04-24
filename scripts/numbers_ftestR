#!/usr/local/bin/Rscript

e1114 <- read.csv("results/e11_e14.csv")
e14p1 <- read.csv("results/e14_p1.csv")

##temporarily, order them by probe ID
e1114.s <- e1114[order(e1114[,"ID"],decreasing=FALSE),]
e14p1.s <- e14p1[order(e14p1[,"ID"],decreasing=FALSE),]

#find rownames of genes that change in both
inds11_14 <- which(e1114.s[,"logFC"]<=-1)
inds14_p1 <- which(abs(e14p1.s[,"logFC"])<=1)

inds <- intersect(inds11_14, inds14_p1)

#find genes at each age
e1114.g <- e1114.s[inds,]
e14p1.g <- e14p1.s[inds,]

#order by FC and then remove duplicates
e1114.og <- e1114.g[order(abs(e1114.g[,"logFC"]),decreasing=TRUE),]
e14p1.og <- e14p1.g[order(abs(e14p1.g[,"logFC"]),decreasing=TRUE),]
e1114.dog <- e1114.og[!duplicated(e1114.og[,"symbol"]),]
e14p1.dog <- e14p1.og[!duplicated(e14p1.og[,"symbol"]),]

#reorder by symbol
e1114.dog<-e1114.dog[order(e1114.dog[,"symbol"],decreasing=TRUE),]
e14p1.dog<-e14p1.dog[order(e14p1.dog[,"symbol"],decreasing=TRUE),]

dat <- data.frame(e1114.dog, e14p1.dog)

#remove all the spare columns
#rm.cols<-c("ID","ID.1","EnsemblID.1")
#keep.cols<-!(colnames(dat) %in% rm.cols)

dat <- dat[,colnames(dat)!=("ID")]
dat <- dat[,colnames(dat)!=("ID.1")]
dat <- dat[,colnames(dat)!=("EnsemblID.1")]
dat <- dat[,colnames(dat)!=("ensembl_gene_id")]
dat <- dat[,colnames(dat)!=("ensembl_gene_id.1")]
dat <- dat[,colnames(dat)!=("symbol.1")]
dat <- dat[,colnames(dat)!=("chromosome_name.1")]
dat <- dat[,colnames(dat)!=("start_position.1")]
dat <- dat[,colnames(dat)!=("end_position.1")]
dat <- dat[,colnames(dat)!=("strand.1")]
dat <- dat[,colnames(dat)!=("t.1")]
dat <- dat[,colnames(dat)!=("description.1")]
dat <- dat[,colnames(dat)!=("B")]
dat <- dat[,colnames(dat)!=("B.1")]
dat <- dat[,colnames(dat)!=("t")]


#rename them to something more sensible
colnames(dat)[8] <- "logFC.e1114"
colnames(dat)[9] <- "AveExpr.e1114"
colnames(dat)[10] <- "P.Value.e1114"
colnames(dat)[11] <- "adj.P.Val.e1114"
colnames(dat)[12] <- "logFC.e14p1"
colnames(dat)[13] <- "AveExpr.e14p1"
colnames(dat)[14] <- "P.Value.e14p1"
colnames(dat)[15] <- "adj.P.Val.e14p1"

#order by fold change in e1114

dat<- dat[order(dat[,"logFC.e1114"], decreasing=TRUE),]

#write.csv
write.csv(dat, "results/siggenes/downsame.csv")









####################old attempts at doing the same thing
#order by FC and then remove duplicates
#e1114.op <- e1114.p[order(abs(e1114.p[,"logFC"]),decreasing=TRUE),]
#e14p1.op <- e14p1.p[order(abs(e14p1.p[,"logFC"]),decreasing=TRUE),]
#e1114.dop <- e1114.op[!duplicated(e1114.op[,"symbol"]),]
#e14p1.dop <- e14p1.op[!duplicated(e14p1.op[,"symbol"]),]

#reorder by symbol so that the data.frame reads it ok
#e1114.dop<-e1114.dop[order(e1114.dop[,"symbol"],decreasing=FALSE),]
#e14p1.dop<-e14p1.dop[order(e14p1.dop[,"symbol"],decreasing=FALSE),]

#find genes at each age that change expression
#e11_14 <- e1114.dop[e1114.dop[,"logFC"]>=1,]
#e14_p1 <- e14p1.dop[e14p1.dop[,"logFC"]<=-1,]

#Find gene symbols that are in both
#inds <- intersect(e11_14[,"symbol"], e14_p1[,"symbol"])


#inds11_14 <- which(e1114.dop[,"logFC"]>1)
#inds14_p1 <- which(e14p1.dop[,"logFC"]>1)

#which<-which(e11_14[,"symbol"] %in% e14p1[,"symbol"])

#Find gene symbols that are in both
#inds <- intersect(e11_14[,"symbol"], e14_p1[,"symbol"])

#inds11_14 <- 

#find the right columns for each age

#e1114.symbols<-e1114.dop[inds,"symbol"]
#e1114.FC<-e1114.dop[inds,"logFC"]
#e1114.pval<-e1114.dop[inds,"adj.P.Val"]
#e14p1.symbols<-e14p1.dop[inds,"symbol"]
#e14p1.FC<-e14p1.dop[inds,"logFC"]
#e14p1.pval<-e14p1.dop[inds,"adj.P.Val"]


#dat.1114<-data.frame(rownames=e1114.symbols,e1114.FC,e1114.pval) 
#dat.14p1<-data.frame(rownames=e14p1.symbols,e14p1.FC,e14p1.pval)





#write.csv(e1114.dop[inds,],"results/siggenes/upup.csv")

