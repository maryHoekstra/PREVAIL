# look at differential expression between early vs. late samples in group A and group B

# earlyVsLate takes an Expression Set, a group code, and vectors specifying which days are considered early and which days are considered late
getESet <- function(filteredESet,groupCode,earlyDays,lateDays) {
  
  # create expression set with desired group code
  oneGroup.eset <- filtered.eset[ , filtered.eset$code==groupCode]
  
  # create expression set with only early and late samples 
  earlyLateDays = c(earlyDays,lateDays)
  earlyLateIndices <- logical(ncol(oneGroup.eset))
  for (i in 1:length(earlyLateDays)) {
    earlyLateIndices <- earlyLateIndices | oneGroup.eset$Day==earlyLateDays[i]
  }
  
  earlyLate.eset <- oneGroup.eset[,earlyLateIndices]
  return (earlyLate.eset)
}

createDesign <- function(earlyLateESet) {  
  earlyIndices <- logical(ncol(earlyLateESet))
  for (i in 1:length(earlyDays)) {
    earlyIndices <- earlyIndices | earlyLateESet$Day==earlyDays[i]
  }
  
  lateIndices <- logical(ncol(earlyLateESet))
  for (i in 1:length(lateDays)) {
    lateIndices <- lateIndices | earlyLateESet$Day==lateDays[i]
  }
  
  # create design matrix with the two groups we are interested in 
  design <- cbind(as.numeric(earlyIndices),as.numeric(lateIndices))
  colnames(design) <- c('early', 'late')
  return (design)
}
 
getDiffs <- function(earlyLateESet,design) {
  EvsL <- "early-late"
  diffs <- runLimma(exprs(earlyLateESet),design,EvsL)
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


earlyDays <- c("day 1")
lateDays <- c("day 14", "day 21", "day 28")

earlyVsLateESet_A <- getESet(filtered.eset,'A',earlyDays,lateDays)
designA <- createDesign(earlyVsLateESet_A)
diffsA <- getDiffs(earlyVsLateESet_A,designA)
geneListA <- getGeneList(diffsA)
write(geneListA,file="/Users/maryhoekstra/Desktop/A_day1vs3.txt")

diffsB <- earlyVsLate(filtered.eset,'B',earlyDays,lateDays)
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
