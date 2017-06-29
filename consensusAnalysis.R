# consensus analysis

# we work with two sets:
nSets = 2
# for easier labeling of plots, create a vector holding descriptive names of the two sets
setLabels = c("Group A", "Group B")
shortLabels = c("A","B")
# form multi-set expression data
multiExpr = vector(mode = "list", length = nSets)
multiExpr[[1]] = list(data = lateExprs_A)
multiExpr[[2]] = list(data = lateExprs_B)

# check that the data has the correct format for many functions operating on multiple sets:
exprSize = checkSets(multiExpr)

# find modules for reference set 
refNet <- getModules(datExprs = lateExprs_A,sfPower = 26)
refGeneLabels <- refNet$colors
refColours <- labels2colors(refGeneLabels)
refMEs <- moduleEigengenes(lateExprs_A,refColours)$eigengenes
refMEs = orderMEs(refMEs, greyName = "ME0",greyLast = TRUE)

# find consesus modules
consNet = blockwiseConsensusModules(
  multiExpr, maxBlockSize = 6074, power = c(26,16), minModuleSize = 30,
  deepSplit = 2, networkType= "signed", TOMType = "signed",
  corType = "bicor",
  pamRespectsDendro = FALSE,
  mergeCutHeight = 0.25, numericLabels = TRUE,
  minKMEtoStay = 0,
  verbose = 5)

consLabels <- consNet$colors
consColours <- labels2colors(consLabels)
consMEs <- consNet$multiMEs
# recalculate consMEs to give them color names
consMEsC = multiSetMEs(multiExpr, universalColors = consColours)
# order MEs the same way as the reference group
consMEs <- orderMEs(consMEsC,greyName = "ME0",greyLast = TRUE) 

# compare consensus modules to reference modules
# isolate the module labels in the order they appear in ordered module eigengenes
refLabels = substring(names(refMEs), 3)
consLabels = substring(names(consMEs[[1]]$data), 3)
# numbers of reference and consensus modules
nRefMods = length(refLabels)
nConsMods = length(consLabels)

# calculate the overlaps of each pair of reference-consensus modules
# use the Fisher’s exact test (also known as hypergeometric test) to assign a p-value to each of the pairwise overlaps
# initialize tables of p-values and of the corresponding counts
pTable = matrix(0, nrow = nRefMods, ncol = nConsMods)
CountTbl = matrix(0, nrow = nRefMods, ncol = nConsMods)
# execute all pairwaise comparisons
for (rmod in 1:nRefMods)
  for (cmod in 1:nConsMods)
  {
    refMembers = (refColours == refLabels[rmod]);
    consMembers = (consColours == consLabels[cmod]);
    pTable[rmod, cmod] = -log10(fisher.test(refMembers, consMembers, alternative = "greater")$p.value);
    CountTbl[rmod, cmod] = sum(refColours == refLabels[rmod] & consColours == consLabels[cmod])
  }

# truncate p values smaller than 10^{-50} to 10^{-50}
pTable[is.infinite(pTable)] = 1.3*max(pTable[is.finite(pTable)]);
pTable[pTable>50 ] = 50 ;
# marginal counts (really module sizes)
refModTotals = apply(CountTbl, 1, sum)
consModTotals = apply(CountTbl, 2, sum)
# actual plotting
sizeGrWindow(10,7 )
par(mfrow=c(1,1))
par(cex = 1.0)
par(mar=c(8, 10.4, 2.7, 1)+0.3)
# use function labeledHeatmap to produce the color-coded table with all the trimmings
# numbers in the table indicate gene counts in the intersection of the corresponding modules
# coloring of the table encodes −log(p), with p being the Fisher’s exact test p-value for the overlap of the two modules
# the stronger the red color, the more significant the overlap is
labeledHeatmap(Matrix = pTable,
               xLabels = paste(" ", consLabels),
               yLabels = paste(" ", refLabels),
               colorLabels = TRUE,
               xSymbols = paste("Cons ", consLabels, ": ", consModTotals, sep=""),
               ySymbols = paste("A ", refLabels, ": ", refModTotals, sep=""),
               textMatrix = CountTbl,
               colors = greenWhiteRed(100)[50:100],
               main = "Correspondence of A-specific and A-B consensus modules (Days 14, 21, 28)",
               cex.text = 1.0, cex.lab = 1.0, setStdMargins = FALSE)

# compare consensus eigenegene networks
# are there two modules that are highly related in one set but not the other?

# order eigengenes by consesus hierarchical clustering:
MET = consensusOrderMEs(consMEsC)

# perform differential analysis and plot
sizeGrWindow(8,10)
par(cex = 0.9)
plotEigengeneNetworks(MET, setLabels, marDendro = c(0,2,2,1), marHeatmap = c(3,3,2,1),
                      zlimPreservation = c(0.5, 1), xLabelsAngle = 90)
