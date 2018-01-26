# EARLY VS LATE
# look at differential expression between early vs. late samples in group A and group B

# last modified: January 15th, 2018

# subsets an expression set based on samples from specific days
getESet <- function(filteredESet,earlyDays,lateDays) {
  
  # create expression set with only early and late samples 
  earlyLateDays = c(earlyDays,lateDays)
  earlyLateIndices <- logical(ncol(filteredESet))
  for (i in 1:length(earlyLateDays)) {
    earlyLateIndices <- earlyLateIndices | filteredESet$Day==earlyLateDays[i]
  }
  earlyLate.eset <- filteredESet[,earlyLateIndices]
  return (earlyLate.eset)
}

# creates design matrix using smaller expression set
createDesign <- function(earlyLateESet,earlyDays,lateDays) {  
  earlyLateVec <- character(ncol(earlyLateESet))
  earlyIndices <- logical(ncol(earlyLateESet))
  for (i in 1:length(earlyDays)) {
    earlyIndices <- earlyIndices | earlyLateESet$Day==earlyDays[i]
  }
  earlyLateVec[earlyIndices]="early"
  lateIndices <- logical(ncol(earlyLateESet))
  for (i in 1:length(lateDays)) {
    lateIndices <- lateIndices | earlyLateESet$Day==lateDays[i]
  }
  earlyLateVec[lateIndices]="late"
  # create design matrix with the four groups we are interested in (A.early, A.late, etc) 
  designVec <- factor(paste(earlyLateESet$code,earlyLateVec,sep="."))
  designMatrix <- model.matrix(~0+designVec)
  colnames(designMatrix) <- levels(designVec)
  return (designMatrix)
}

# gets a table of top differentially expressed probes for a given group from limma analysis
getDiffs <- function(earlyLateESet,groupCode,design) {
  contrastsString <- paste(groupCode,".early-",groupCode,".late",sep="")
  diffs <- runLimma(earlyLateESet,design,contrastsString)
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
lateDays <- c("day 14","day 21","day 28")

earlyVsLateESet <- getESet(filtered.eset,earlyDays,lateDays)
design <- createDesign(earlyVsLateESet,earlyDays,lateDays)
diffsA <- getDiffs(earlyVsLateESet,groupCode="A",design)
geneListA <- getGeneList(diffsA)
write(geneListA,file="/Users/maryhoekstra/Desktop/A_day1vs7.txt")

diffsB <- getDiffs(earlyVsLateESet,groupCode="B",design)
geneListB <- getGeneList(diffsB)
write(geneListB,file="/Users/maryhoekstra/Desktop/B_day1vs7.txt")

commonGenes <- intersect(geneListA,geneListB)
write(commonGenes,file="/Users/maryhoekstra/Desktop/common.txt")

uniqueToA <- setdiff(geneListA,geneListB)
write(uniqueToA,file="/Users/maryhoekstra/Desktop/uniqueA.txt")

uniqueToB <- setdiff(geneListB,geneListA)
write(uniqueToB,file="/Users/maryhoekstra/Desktop/uniqueB.txt")

allGenes <- union(geneListA,geneListB)
write(allGenes,file="/Users/maryhoekstra/Desktop/union.txt")
