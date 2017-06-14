# collapse probes into genes
library(WGCNA)

collapseProbes <- function(filteredExprs,annoData) {
  
  # find indices of our chosen probes in the original annotation matrix
  probeIndices <- match(rownames(filteredExprs),annoData$ID)
  
  #find indices of our chosen probes in the expanded annotation matrix 
  allIDS <- as.character(annoData$Entrez.Gene)
  # get the gene symbols which correspond to the chosen probes
  IDsubset <- allIDS[probeIndices]
  
  # choose probe with highest mean expression when multiple probes map to a gene
  # use only the symbols from the filtered expression set
  collapsedData <- collapseRows(filteredExprs,IDsubset,rownames(filteredExprs),method="MaxMean")
  collapsedExprs <- collapsedData$datETcollapsed
  
  # remove row with "---" as gene symbol
  collapsedExprs <- collapsedExprs[-1,]
  
  return(collapsedExprs) 
}

collapsedExprs <- collapseProbes(filteredExprs,annoData)

# change rownames of collapsed expression matrix

# for every row of the collapsed expression matrix, if Gene Symbol column contains "///", split string into tokens 
entrezIDs <- rownames(collapsedExprs)
for (i in 1:length(entrezIDs)) {
  rowname <- entrezIDs[i]
  if (grepl("///",rowname)) {
    tokens <- strsplit(rowname,"///")
    tokenArray <- tokens[[1]]
    # append first symbol
    entrezIDs[i] <- trimws(tokenArray[1])
  }
}
rownames(collapsedExprs) <- entrezIDs
