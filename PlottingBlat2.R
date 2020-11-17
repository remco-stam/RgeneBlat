#### REMCO STAM 2020
#### Simple script used to visualise BLAT hits
## Step one, turns individual blocks in .psl files into individual lines of a bed file. 
## Step two plost the bed file using the sushi package
## In a third step an annotated reference is added. 


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#BiocManager::install("Sushi")
library("Sushi")


dat <- read.table("/path/to/Rpiamrs_vs_bac_Score93Edit.csv", header = T)
#pslfile <- dat[2,]
#pslfile2 <- dat[3,]

SplitIntoDataframe <- function(pslfile, n){
  starts <- as.vector(pslfile$tStarts)
  len <- as.vector(pslfile$blockSizes)
  
  startsplit <- unlist(strsplit(as.character(starts), ","))
  lensplit <- unlist(strsplit(as.character(len), ","))
  endsplit <- as.numeric(startsplit) + as.numeric(lensplit)
  
  chrom <- rep(pslfile$Tname, length(startsplit))
  start <- as.numeric(startsplit)
  stop <- as.numeric(endsplit)
  name <- rep(pslfile$Qname, length(startsplit))
  score <- as.numeric(rep(pslfile$match, length(startsplit)))
  strand <- rep(pslfile$strand, length(startsplit))
  row <- as.numeric(rep(n, length(startsplit)))
  
  df <- data.frame(
    chrom, start, stop, name, score, strand, row)

  return(df)
}  
  
testdata <- rbind(SplitIntoDataframe(dat[1,],1), SplitIntoDataframe(dat[2,],2), SplitIntoDataframe(dat[3,],3), SplitIntoDataframe(dat[4,],4),
                  SplitIntoDataframe(dat[5,],5), SplitIntoDataframe(dat[6,],6), SplitIntoDataframe(dat[7,],7), SplitIntoDataframe(dat[8,],8),
                  SplitIntoDataframe(dat[9,],9), SplitIntoDataframe(dat[10,],10), SplitIntoDataframe(dat[11,],11), SplitIntoDataframe(dat[12,],12),
                  SplitIntoDataframe(dat[13,],13), SplitIntoDataframe(dat[14,],14), SplitIntoDataframe(dat[15,],15), SplitIntoDataframe(dat[16,],16),
                  SplitIntoDataframe(dat[17,],17))
#Plot BLAT hits
plotBed(testdata,chrom = "BAC_5G_(reversed)",
        chromstart = 0 ,chromend = 80000,
        type = "region", row = "given", rownumber = testdata$row,
        plotbg="grey95", rowlabels=unique(testdata$name),
        rowlabelcex=0.75)
labelgenome("BAC_5G_(reversed)",0,80000,n=15,scale="Kb")
mtext("Rpi-amr mapped to reference",side=3, adj=-0.065,line=0.5,font=2)

#Read annotated NLRs
dat2 <- read.table("/path/to/annotatedNLRs.bed", header = F, skip = 1)
testdata2 <- dat2[,c(1:6)]
testdata2$row <- 0
colnames(testdata2) <-  c("chrom", "start", "stop", "name", "score", "strand", "row") 
testdata2$chrom <- "BAC_5G_(reversed)"
testdata2$name <- "annotated-nlrs"

#Merge annotatedNLRs and BLAT results
testdata3 <- rbind(testdata, testdata2)

#Plot all
plotBed(testdata3,chrom = "BAC_5G_(reversed)",
        chromstart = 0 ,chromend = 80000,
        type = "region", row = "given", rownumber = testdata3$row,
        plotbg="grey95", rowlabels=unique(testdata3$name),
        rowlabelcex=0.75)
labelgenome("BAC_5G_(reversed)",0,80000,n=15,scale="Kb")
mtext("Rpi-amr mapped to reference",side=3, adj=-0.065,line=0.5,font=2)

#File edited for publication in Inkscape