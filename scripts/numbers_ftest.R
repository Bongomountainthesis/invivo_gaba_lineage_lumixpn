#!/usr/local/bin/Rscript
dat <- read.csv("results/f_test.csv")

#find rownames of genes that change in both
inds11_14 <- which(dat[,"E11.5vE14.5"]>=1)
inds14_p1 <- which(dat[,"E14.5vP1"]<=-1)

inds <- intersect(inds11_14, inds14_p1)

#find genes at each age
dat.g <- dat[inds,]

#order by FC and then remove duplicates
dat.og <- dat.g[order(abs(dat.g[,"E11.5vE14.5"]),decreasing=TRUE),]
dat.dog <- dat.og[!duplicated(dat.og[,"symbol"]),]
dat.og <- dat.g[order(abs(dat.g[,"E14.5vP1"]),decreasing=TRUE),]
dat.dog <- dat.og[!duplicated(dat.og[,"symbol"]),]

#remove nonsig genes <0.05
dat.sdog <- dat.dog[dat.dog[,"adj.P.Val"]<=0.05,]
#order by logFC in E11.5vE14.5
dat.osdog<- dat.sdog[order(abs(dat.sdog[,"E11.5vE14.5"]),decreasing=TRUE),]

#write to csv file
write.csv(dat.sdog,"results/siggenes/with_ftest/updown.csv")
