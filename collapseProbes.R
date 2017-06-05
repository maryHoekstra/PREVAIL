# collapse probes into genes

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
