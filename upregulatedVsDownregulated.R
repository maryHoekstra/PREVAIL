# determine which genes are upregulated or downregulated for a given diffs table 

# takes the expanded diffs table and a string indicating the direction of regulation, where "up" means upregulation and "down" means downregulation
subsetByRegulation <- function(diffsWithAnno,directionString) {
  if (directionString=="up") {
    diffs <- diffsWithAnno[diffsWithAnno$logFC>0,]
  }
  else {
    diffs <- diffsWithAnno[diffsWithAnno$logFC<0,]
  }
  geneEntrez <- dplyr::select(diffs, Entrez.Gene)
  geneList <- unlist(sapply(geneEntrez[,1], splitGeneList))
  geneList <- geneList[geneList != "---"]
  geneList <- gsub(' ', '', geneList)
  
  return(geneList)
}

upregulatedGeneListA <- subsetByRegulation(diffsA,"up")
upregulatedGeneListB <- subsetByRegulation(diffsB,"up")
downregulatedGeneListA <- subsetByRegulation(diffsA,"down")
downregulatedGeneListB <-subsetByRegulation(diffsB,"down")

downregulatedGeneListB_unique <- setdiff(downregulatedGeneListB,downregulatedGeneListA)
upregulatedGeneListB_unique <- setdiff(upregulatedGeneListB,upregulatedGeneListA)

write(downregulatedGeneListA,file="/Users/maryhoekstra/Desktop/downregulatedGeneListA.txt")
write(downregulatedGeneListB_unique,file=paste("/Users/maryhoekstra/Desktop/downregulatedGeneListB_unique.txt"))
write(downregulatedGeneListB,file=paste("/Users/maryhoekstra/Desktop/downregulatedGeneListB.txt"))
write(upregulatedGeneListB_unique,file=paste("/Users/maryhoekstra/Desktop/upregulatedGeneListB_unique.txt"))
write(upregulatedGeneListA,file=paste("/Users/maryhoekstra/Desktop/upregulatedGeneListA.txt"))
write(upregulatedGeneListB,file=paste("/Users/maryhoekstra/Desktop/upregulatedGeneListB.txt"))
