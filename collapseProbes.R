# collapse probes into genes
library(WGCNA)

collapseProbes <- function(filteredExprs,annoData) {
  
  # find indices of our chosen probes in the original annotation matrix
  probeIndices <- match(rownames(filteredExprs),annoData$ID)
  
  #find indices of our chosen probes in the expanded annotation matrix 
  allSymbols <- as.character(annoData$Gene.Symbol)
  # get the gene symbols which correspond to the chosen probes
  symbolSubset <- allSymbols[probeIndices]
  
  # choose probe with highest mean expression when multiple probes map to a gene
  # use only the symbols from the filtered expression set
  collapsedData <- collapseRows(filteredExprs,symbolSubset,rownames(filteredExprs),method="MaxMean")
  collapsedExprs <- collapsedData$datETcollapsed
  
  # remove row with "---" as gene symbol
  collapsedExprs <- collapsedExprs[-1,]
  
  return(collapsedExprs) 
}

# change rownames of collapsed expression matrix

# for every row of the collapsed expression matrix, if Gene Symbol column contains "///", split string into tokens 
geneSymbols <- rownames(collapsedExprs)
for (i in 1:length(geneSymbols)) {
  rowname <- geneSymbols[i]
  if (grepl("///",rowname)) {
    tokens <- strsplit(rowname,"///")
    tokenArray <- tokens[[1]]
    # append first symbol
    geneSymbols[i] <- trimws(tokenArray[1])
  }
}
rownames(collapsedExprs) <- geneSymbols
