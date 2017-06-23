# compare sets of modules from two networks
# two modules are similar if they share 75% of the same genes

# compare each module in network one with each module in network two 
allGenes <- colnames(earlyExprs_A)
modIndices_A = !(modColors_A %in% c("grey","gold"))
modColors_A = modColors_A[modIndices_A]
modIndices_B = !(modColors_B %in% c("grey","gold"))
modColors_B = modColors_B[modIndices_B]
for (i in 1:length(modColors_A)) {
  refModuleIDs <- getModuleGenes(geneList=allGenes,geneColours = geneColors_A,moduleColour = modColors_A[i])
  for (j in 1:length(modColors_B)) {
    print(paste(modColors_A[i],"vs.",modColors_B[j]))
    testModuleIDs <- getModuleGenes(geneList = allGenes, geneColours = geneColors_B, moduleColour = modColors_B[j])
    commonGenes <- intersect(refModuleIDs,testModuleIDs)
    print(paste("num common genes:",length(commonGenes)))
    fractionShared <- length(commonGenes) / max(length(refModuleIDs),length(testModuleIDs))
    if (fractionShared >= 0.75) {
      print(paste(modColors_A[i],"in A and",modColors_B[j],"in B are comparable: f =",fractionShared))
    }
  }
}
