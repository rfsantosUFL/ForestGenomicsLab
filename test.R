#LDna
#### Assign R CMD BATCH Arguments ####
args = commandArgs(TRUE)
CORE=args[1]
MKR.ID=args[2]

#CORE=4
#MKR.ID="cap_fs2"

####LOAD r2 DATA####
library("stringr")
library("data.table")
SCAFFOLDS = 11
LD.r2<-data.frame(NULL)
for (i in seq(1:SCAFFOLDS)) {
  LD.r2<-rbind(LD.r2,(fread(sprintf("LD.%s.scaffold_%s_PROBE.r2v.csv", MKR.ID, i), stringsAsFactor=F,h=T, sep=",")))
}

LD.r2$pos.loc1<-as.numeric(unlist(str_extract(LD.r2$loc1,"[0-9]{3,}")))
LD.r2$pos.loc2<-as.numeric(unlist(str_extract(LD.r2$loc2,"[0-9]{3,}")))
LD.r2$dist.bp<-LD.r2$pos.loc2 - LD.r2$pos.loc1

#Converting LDcorSV LD table to matrix
LD.mat<- as.matrix(as.dist(xtabs(r2 ~ loc2 + loc1, data=LD.r2),diag = T))

library("LDna")
ldna <- LDnaRaw(LD.mat, mc.cores = CORE)

#png(sprintf("LD.clusters.%s.%d.PROBE.png", MKR.ID, 1:2), width = 800, height = 800)
png("LDna.clusters.%d.png", width = 800, height = 800)
clusters <- extractClusters(ldna, min.edges = 20, phi = 2, lambda.lim = NULL,
                            rm.COCs = TRUE, extract = TRUE, plot.tree = TRUE, plot.graph = TRUE)
dev.off()

summary <- summaryLDna(ldna, clusters, LD.mat)

#png("LD.plot.network.%d.png", width = 800, height = 800)
#plotLDnetwork(ldna = ldna, LDmat = LD.mat, option=2, clusters=clusters, summary=summary, full.network=FALSE, include.parent=FALSE, after.merger=FALSE)
#dev.off()
library(parallel)
fun <- function(x){
  setEPS()
  postscript(paste("LDna.plot.network",x,"eps" ,sep="."))
  plotLDnetwork(ldna, LD.mat, option=2, clusters=clusters[x], summary=summary[x,], full.network=F, include.parent=T)
  dev.off()
}
mclapply(1:length(clusters), fun, mc.cores=CORE, mc.preschedule=TRUE)
