# perform a limma (Empirical Bayes) moderated t-test
library(limma)

runLimma <- function(exprsMatrix,designMatrix,contrastsString) {
  
  # fit a linear model using the design
  fit <- lmFit(exprsMatrix, designMatrix)
  
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

