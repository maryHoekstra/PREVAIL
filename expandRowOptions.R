# option 1 for solving one-to-many probe problem

# expand rows where one probe maps to several genes
# create new data frame to hold probe ID and Gene Symbol
expandedAnnotation <- data.frame(matrix(ncol=2))
colnames(expandedAnnotation) <- c("ID","Gene.Symbol")
# for every row of the annotation dataframe, if Gene Symbol column contains "///", split string into tokens 
for (i in 1:50) {
  row <- PVannotation[i,]
  probeID <- as.character(row$ID)
  geneSymbol <- as.character(row$Gene.Symbol)
  if (grepl("///",geneSymbol)) {
    tokens <- strsplit(geneSymbol,"///")
    tokenArray <- lapply(tokens[[1]],FUN=trimws)
    # create row in new dataframe for every token
    for (j in 1:length(tokenArray)) {
      token <- tokenArray[j]
      expandedAnnotation <- rbind(expandedAnnotation,c(probeID,token))
    }
  }
  else {
    # just append the single ID and gene symbol
    expandedAnnotation <- rbind(expandedAnnotation,c(probeID,geneSymbol))
  }
}

expandedAnnotation <- expandedAnnotation[-1,]


# option 2 for solving one-to-many probe problem

# make new dataframe with columns for probe ID and Gene Symbol
refinedAnnotation <- data.frame(matrix(ncol=2))
colnames(refinedAnnotation) <- c("ID","Gene.Symbol")

# for every row of the annotation dataframe, if Gene Symbol column contains "///", split string into tokens 
for (i in 1:50) {
  row <- PVannotation[i,]
  probeID <- as.character(row$ID)
  geneSymbol <- as.character(row$Gene.Symbol)
  if (grepl("///",geneSymbol)) {
    tokens <- strsplit(geneSymbol,"///")
    tokenArray <- tokens[[1]]
    refinedAnnotation <- rbind(refinedAnnotation,c(probeID,trimws(tokenArray[1])))
  }
  else {
    refinedAnnotation <- rbind(refinedAnnotation,c(probeID,geneSymbol))
  }
}

# remove line containg NA
refinedAnnotation <- refinedAnnotation[-1,]

###