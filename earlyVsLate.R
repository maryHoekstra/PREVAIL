# look at differential expression between early vs. late samples in group A and group B
# early = day 1, late = days 14, 21, 28

earlyVsLate <- function(groupCode) {
  
  # create expression set with desired group code
  oneGroup.eset <- filtered.eset[ , filtered.eset$code==groupCode]
  
  # create expression set with only early and late samples 
  earlyLate.eset <- oneGroup.eset[ , oneGroup.eset$Day == "day 1" | oneGroup.eset$Day == 'day 14' | oneGroup.eset$Day == 'day 21' | oneGroup.eset$Day == 'day 28']
  early <- as.integer(earlyLate.eset$Day == 'day 1')
  late <- as.integer(earlyLate.eset$Day == 'day 14' | earlyLate.eset$Day == 'day 21' | earlyLate.eset$Day == 'day 28')
  
  # create design matrix with the two groups we are interested in 
  design <- cbind(early, late)
  colnames(design) <- c('early', 'late')
  EvsL <- "early-late"
  
  diffs <- runLimma(exprs(earlyLate.eset),design,EvsL)
  
  # get entrez gene IDs for the differentially expressed probes
  diffs$ID <- as.factor(rownames(diffs))
  diffsWithAnno <- merge(diffs, dplyr::select(annoData, ID, Entrez.Gene))
  
  return(diffsWithAnno)
}

splitGeneList <-function(x) {
  return(strsplit(as.character(x), '///')[[1]])
}

getGeneList <- function(diffsWithAnno) {
  geneEntrez <- dplyr::select(diffsWithAnno, Entrez.Gene)
  geneList <- unlist(sapply(geneEntrez[,1], splitGeneList))
  geneList <- geneList[geneList != "---"]
  geneList <- gsub(' ', '', geneList)
  return (geneList)
}


diffsA <- earlyVsLate('A')
geneListA <- getGeneList(diffsA)
write(geneListA,file="/Users/maryhoekstra/Desktop/A.txt")

diffsB <- earlyVsLate('B')
geneListB <- getGeneList(diffsB)
write(geneListB,file="/Users/maryhoekstra/Desktop/B.txt")

commonGenes <- intersect(geneListA,geneListB)
write(commonGenes,file="/Users/maryhoekstra/Desktop/common.txt")

uniqueToA <- setdiff(geneListA,geneListB)
write(uniqueToA,file="/Users/maryhoekstra/Desktop/uniqueA.txt")

uniqueToB <- setdiff(geneListB,geneListA)
write(uniqueToB,file="/Users/maryhoekstra/Desktop/uniqueB.txt")

allGenes <- union(geneListA,geneListB)
write(allGenes,file="/Users/maryhoekstra/Desktop/union.txt")
