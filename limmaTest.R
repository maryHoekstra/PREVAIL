# perform a limma (Empirical Bayes) moderated t-test
library(limma)
library(statmod)

# given an expression set and contrasts, performs differential expression analysis on expression data
# contrastsDesign is a design matrix with the two groups we are interested in comparing
# contrastsString specifies how we wish to compare the two groups ("early-late","late-early",etc)
runLimma <- function(exprsSet,designMatrix,contrastsString) {
  
  # account for correlations between samples of the same patient
  patientIDs <- exprsSet$patient_id
  corfit <- duplicateCorrelation(exprsSet,designMatrix,block=patientIDs)
  # consensus represents the correlation between measurements made on the same patient 
  print(corfit$consensus)
  
  # fit a linear model using the design
  fit <- lmFit(exprsSet, designMatrix, block=patientIDs,correlation=corfit$consensus)
  
  # create a contrast matrix 
  # a positive logFC for a gene means that it's upregulated in the early group
  cont.matrix <- makeContrasts(contrasts=contrastsString, levels=designMatrix)
  
  # fit a linear model with specified contrasts
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2)
  # adjust p-values based on Holm's method, which controls family-wide error rate
  diffs <- topTable(fit2, p.value=0.05, adjust='holm', number=5000)

  return(diffs)
}


