#!/usr/local/bin/Rscript

library(beadarray)
library(illuminaMousev2BeadID.db)
library(illuminaMousev2.db)

BSData <- get(load("data/BSData.quantile.RData"))

##for testing
symbol <- as.matrix(c("Dazap2","Bcl6","Kdm4b","Jazf1","Trim39"))

#args <- commandArgs(trailingOnly=TRUE)
#symbol <- as.matrix(c(args[1]))

annot <- mget(ls(illuminaMousev2SYMBOL), illuminaMousev2SYMBOL, ifnotfound = NA)

res <- matrix(nrow = nrow(symbol), ncol = 2)
for(i in 1:nrow(symbol)){
	res[i,1] <- symbol[i,]
	res[i,2] <- annot[which(annot == symbol[i,]),]
	}

res[i,2 <- annot[which(annot == symbol[i,]),]

###meh doesnt work....
symbol <- as.matrix(c("Dazap2","Bcl6","Kdm4b","Jazf1","Trim39"))

range <- range(E)

res1 <- ls(annot[grep(symbol[1,], annot, ignore.case=T)])
res2 <- ls(annot[grep(symbol[2,], annot, ignore.case=T)])
res3 <- ls(annot[grep(symbol[3,], annot, ignore.case=T)])
res4 <- ls(annot[grep(symbol[4,], annot, ignore.case=T)])
res5 <- ls(annot[grep(symbol[5,], annot, ignore.case=T)])

res <- c(res1,res2,res3,res4,res5)
E<- exprs(BSData)

E_res <- E[which(rownames(E) %in% res),]

# bcl6 has two, kdm4b has three

rownames(E_res) <- c("Bcl6","Bcl6b","Dazap2","Jazf1","Kdm4b_2","Kdm4b_3","Kdm4b_1","Trim39")
#rownames(E_res) <- c("Dazap2","Bcl6","Bcl6b","Kdm4b_1","Kdm4b_2","Kdm4b_3","Jazf1","Trim39")

#remove bcl6b

E_res <- E_res[which(!(rownames(E_res) == "Bcl6b")),]

##average up
E_avg <- matrix(nrow=nrow(E_res),ncol=3)

rownames(E_avg) <- rownames(E_res)
colnames(E_avg) <- c("E11.5","E14.5","P1")

for(i in 1:nrow(E_res)){
	E_avg[i,1] <- mean(E_res[i,1:2])
	E_avg[i,2] <- mean(E_res[i,3:6])
	E_avg[i,3] <- mean(E_res[i,7:8])
	}

E_avg <- E_avg[c(1,2,3,5,7),]























