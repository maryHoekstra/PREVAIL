allGenes <- colnames(datExprs)

getModuleGenes <- function(geneList,geneColours,moduleColour) {
  moduleGenes <- (geneColours==moduleColour)
  moduleIDs <- geneList[moduleGenes]
  return (moduleIDs)
}

writeFile <- function(IDs,fileName) {
  write(IDs,file=paste("/Users/maryhoekstra/Desktop/",fileName,".txt",sep=""))  
}

nonpreservedModules <- c("grey60","tan","lightgreen","purple","greenyellow")
for (i in 1:length(nonpreservedModules)) {
  moduleIDs <- getModuleGenes(geneList=allGenes,geneColours = refColours,moduleColour = nonpreservedModules[i])
  writeFile(moduleIDs,paste(nonpreservedModules[i],"_LvsCons_A",sep = ""))
}
