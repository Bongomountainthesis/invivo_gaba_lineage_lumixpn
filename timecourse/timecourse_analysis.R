############ Load, QC and normalise ###############

library(beadarray)

dataFile = "data/Raw_Data_UnNorm_Rewritten.txt"

BSData <- readBeadSummaryData(dataFile=dataFile,
                              skip=0,
                                ProbeID="ProbeID",
                              columns = list(exprs = "AVG_Signal",
                                se.exprs="BEAD_STDERR",
                                NoBeads = "Avg_NBEADS",
                                Detection="Detection Pval"
                                             ),
                              )
#> colnames(exprs(BSData))
# [1] "P4 RNA1"  "Lhx6_2"   "P10 RNA4" "Lhx6_4"   "Lhx6_12"  "Lhx6_1"  
# [7] "Lhx6_11"  "P10 RNA3" "Lhx6_7"   "P10 RNA2" "Lhx6_6"   "Lhx6_3"  
#[13] "P21 RNA2" "P21 RNA1" "Lhx6_8"   "P4 RNA2"  "Lhx6_9"   "P21 RNA3"
#[19] "P10 RNA1" "P4 RNA3"  "P21 RNA4" "Lhx6_10"  "Lhx6_5"   "E13.5_22"
#[25] "E13.5_21" "E12.5_20" "E12.5_12" "E13.5_20" "10"       "E13.5_25"
#[31] "8"        "4"        "E12.5_13" "E13.5_19" "6"        "7"       
#[37] "E12.5_4" 

# assume we're interested in the Lhx6 data at this stage
# 3 timepoints, 4 replicates at each:
# E11.4 = Lhx6_1..4
# E14.5 = Lhx6_5..8
# P1    = Lhx6_9..11
lhx <- paste("Lhx6_", 1:12, sep="")

BSData<-BSData[,lhx]
save(BSData, file="data/BSData.RData")


# some QC plots:

postscript(file="results/Boxplotprenorm.ps", horizontal=FALSE)
boxplot(as.data.frame(log2(exprs(BSData))),  las = 2, outline = FALSE, ylab = "log2(intensity)", col=c(2,2,2,2,4,4,4,4,6,6,6,6), main="Log2 Expression")
dev.off()

postscript(file="results/plotMAXYprenorm.ps", horizontal=FALSE)
plotMAXY(exprs(BSData), arrays = 1:12,  pch = 16)
dev.off()

postscript(file="results/plotDENSITYprenorm.ps", horizontal=FALSE)
E <- log2(exprs(BSData))
plot(density(E[,1]))
for(i in 1:12){
  lines(density(E[,i]),col=i)
}
dev.off()

#From the QC plots, clearly some of these haven't worked:
fail <- c("Lhx6_1", "Lhx6_2", "Lhx6_9", "Lhx6_12")
BSData<-BSData[, !colnames(exprs(BSData)) %in% fail ]

# Quantile Normalise
BSData.quantile = normaliseIllumina(BSData, method = "quantile", transform = "log2")


postscript(file="results/Boxplotpostnorm.ps", horizontal=FALSE)
boxplot(as.data.frame((exprs(BSData.quantile))),  las = 2, outline = FALSE, ylab = "log2(intensity)", col=c(2,2,4,4,4,4,6,6), main="Log2 Expression")
dev.off()

postscript(file="results/plotMAXYpostnorm.ps", horizontal=FALSE)
plotMAXY(exprs(BSData.quantile), arrays = 1:8, log = FALSE, pch = 16)
dev.off()

save(BSData.quantile, file="data/BSData.quantile.RData")

postscript(file="results/plotDENSITYpostnorm.ps", horizontal=FALSE)
E <- exprs(BSData.quantile)
plot(density(E[,1]))
for(i in 1:8){
  lines(density(E[,i]),col=i)
}
dev.off()


################ Limma #################

library(limma)
library(biomaRt)
library(beadarray)

# Bioconductor packages are based on Mark Dunning's ReMOAT reannotations.
# The BeadID package uses numeric ArrayAddressIDs as keys and is for data 
# that has been summarised from bead level - so this is the file you need 
# if you are using the beadarray library to read in the raw bead data.
library(illuminaMousev2BeadID.db) 

# This package has much the same annotation, again based on the ReMOAT process
# but with ProbeIDs (ILMN_*) as keys. This is for data that has already been 
# collapsed into a single value per probe. Which is what we have in this case.
#
# Note that crap probes which don't map to anything or cross map to loads of
# stuff have been removed before annotation.
library(illuminaMousev2.db) 

BSData <- get(load("BSData.quantile.e11MGE.RData"))
Em <- exprs(BSData.quantile)
BSData <- get(load("BSData.quantile.e11e14p1.RData"))
Et <- exprs(BSData.quantile)

# We only have 3 timepoints, which isn't *really* enough to treat the
# data as a proper timecourse, so just run Limma on all three possible
# pairwise comparisons.
design<-matrix(0,nrow=(ncol(E)), ncol=3)
colnames(design) <- c("E11.5","E14.5","P1")
rownames(design) <- colnames(E)
design[1:2,1] <- 1
design[3:6,2] <- 1
design[7:8,3] <- 1

cont.matrix<-makeContrasts(E11.5vE14.5=E14.5-E11.5, 
			   E11.5vP1=P1-E11.5,
                           E14.5vP1=P1-E14.5, 
                           levels=design)


fit<-lmFit(E, design)
fit<-contrasts.fit(fit, cont.matrix)
ebFit<-eBayes(fit)

ids = rownames(E)

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

ebFit$genes = anno


# We're really doing a 1-way ANOVA here, with a post-hoc test between all possible
# comparisons.  The Ftest will tell you which differences are significant, so these are the
# p-values and adjusted p-values you are interested in.
# The contrasts can then tell you where the significant difference *is*. 
f.test<- topTable(ebFit, number=nrow(E))
rownames(f.test)<-f.test$ID
f.test<-f.test[order(f.test[,"P.Value"]),]
write.csv(f.test,"results/f_test.csv",row.names=F)

e11.e14<-topTable(ebFit, coef=1, adjust="BH", number=nrow(E))
rownames(e11.e14)<-e11.e14$ID
e11.e14<-e11.e14[order(e11.e14[,"P.Value"]),]
write.csv(e11.e14,"results/e11_e14.csv",row.names=F)

e11.p1<-topTable(ebFit, coef=2, adjust="BH", number=nrow(E))
rownames(e11.p1)<-e11.p1$ID
e11.p1<-e11.p1[order(e11.p1[,"P.Value"]),]
write.csv(e11.p1,"results/e11_p1.csv",row.names=F)

e14.p1<-topTable(ebFit, coef=3, adjust="BH", number=nrow(E))
rownames(e14.p1)<-e14.p1$ID
e14.p1<-e14.p1[order(e14.p1[,"P.Value"]),]
write.csv(e14.p1,"results/e14_p1.csv",row.names=F)













