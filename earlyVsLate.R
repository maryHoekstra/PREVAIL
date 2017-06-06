# look at differential expression between early vs. late samples in group A and group B

# subsets an expression set based on group code and samples from specific days
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

# creates design matrix using smaller expression set
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
 
# gets a table of top differentially expressed probes from limma analysis
getDiffs <- function(earlyLateESet,design) {
  EvsL <- "early-late"
  diffs <- runLimma(exprs(earlyLateESet),design,EvsL)
  # get entrez gene IDs for the differentially expressed probes
  diffs$ID <- as.factor(rownames(diffs))
  diffsWithAnno <- merge(diffs, dplyr::select(annoData, ID, Entrez.Gene))
  return(diffsWithAnno)
}

# splits entries containing slashes
splitGeneList <-function(x) {
  return(strsplit(as.character(x), '///')[[1]])
}

# gets a list of entrez gene IDs from the diffs table
getGeneList <- function(diffsWithAnno) {
  geneEntrez <- dplyr::select(diffsWithAnno, Entrez.Gene)
  geneList <- unlist(sapply(geneEntrez[,1], splitGeneList))
  geneList <- geneList[geneList != "---"]
  geneList <- gsub(' ', '', geneList)
  return (geneList)
}


earlyDays <- c("day 1")
lateDays <- c("day 7")

earlyVsLateESet_A <- getESet(filtered.eset,'A',earlyDays,lateDays)
designA <- createDesign(earlyVsLateESet_A)
diffsA <- getDiffs(earlyVsLateESet_A,designA)
geneListA <- getGeneList(diffsA)
write(geneListA,file="/Users/maryhoekstra/Desktop/A_day1vs3.txt")

earlyVsLateESet_B <- getESet(filtered.eset,'B',earlyDays,lateDays)
designB <- createDesign(earlyVsLateESet_B)
diffsB <- getDiffs(earlyVsLateESet_B,designB)
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
