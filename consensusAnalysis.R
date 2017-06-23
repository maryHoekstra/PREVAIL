# consensus analysis

# we work with two sets:
nSets = 2
# for easier labeling of plots, create a vector holding descriptive names of the two sets
setLabels = c("Group A", "Group B")
shortLabels = c("A", "B")
# form multi-set expression data
multiExpr = vector(mode = "list", length = nSets)
multiExpr[[1]] = list(data = lateExprs_A)
multiExpr[[2]] = list(data = lateExprs_B)

# check that the data has the correct format for many functions operating on multiple sets:
exprSize = checkSets(multiExpr)

# find modules for reference set 
refNet <- getModules(datExprs = lateExprs_A,sfPower = 12)
refLabels <- refNet$colors
refColours <- labels2colors(refLabels)
refMEs <- moduleEigengenes(lateExprs_A,refColours)$eigengenes
refMEs = orderMEs(refMEs, greyName = "ME0")

# find consesus modules
consNet = blockwiseConsensusModules(
  multiExpr, maxBlockSize = 6074, power = c(26,12), minModuleSize = 30,
  deepSplit = 2, networkType= "signed", TOMType = "signed",
  corType = "bicor",
  pamRespectsDendro = FALSE,
  mergeCutHeight = 0.25, numericLabels = TRUE,
  minKMEtoStay = 0,
  saveTOMs = TRUE, verbose = 5)

consLabels <- consNet$colors
consColours <- labels2colors(consLabels)
# recalculate consMEs to give them color names
consMEsC = multiSetMEs(multiExpr, universalColors = consColours)
# order eigengenes by consesus hierarchical clustering:
MET = consensusOrderMEs(consMEsC)

# perform differential analysis and plot
sizeGrWindow(8,10)
par(cex = 0.9)
plotEigengeneNetworks(MET, setLabels, marDendro = c(0,2,2,1), marHeatmap = c(3,3,2,1),
                      zlimPreservation = c(0.5, 1), xLabelsAngle = 90)

# compare consensus modules to reference modules
# isolate the module labels in the order they appear in ordered module eigengenes
refLabels = substring(names(refMEs), 3)
consLabels = substring(names(consMEs[[1]]$data), 3)
# convert the numeric module labels to color labels
refModules = labels2colors(as.numeric(refLabels))
consModules = labels2colors(as.numeric(consLabels))
# numbers of reference and consensus modules
nRefMods = length(refModules)
nConsMods = length(consModules)
# initialize tables of p-values and of the corresponding counts
pTable = matrix(0, nrow = nRefMods, ncol = nConsMods)
CountTbl = matrix(0, nrow = nRefMods, ncol = nConsMods)
# execute all pairwaise comparisons
for (rmod in 1:nRefMods)
  for (cmod in 1:nConsMods)
  {
    refMembers = (refColours == refModules[rmod]);
    consMembers = (moduleColours == consModules[cmod]);
    pTable[rmod, cmod] = -log10(fisher.test(refMembers, consMembers, alternative = "greater")$p.value);
    CountTbl[rmod, cmod] = sum(refColours == refModules[rmod] & moduleColours == consModules[cmod])
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
labeledHeatmap(Matrix = pTable,
               xLabels = paste(" ", consModules),
               yLabels = paste(" ", refModules),
               colorLabels = TRUE,
               xSymbols = paste("Cons ", consModules, ": ", consModTotals, sep=""),
               ySymbols = paste("Ref ", refModules, ": ", refModTotals, sep=""),
               textMatrix = CountTbl,
               colors = greenWhiteRed(100)[50:100],
               main = "Correspondence of reference set-specific and reference-test consensus modules",
               cex.text = 1.0, cex.lab = 1.0, setStdMargins = FALSE)
    
  