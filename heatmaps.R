library(gplots)

createHeatmap <- function(diffs,earlyLateESet,design,titleString) {

  # take max of 500 top differentially expressed probes
  if (nrow(diffs) < 500) 
    numProbes=nrow(diffs)
  else
    numProbes=500
  
  # need the subsetted expression matrix, the diffs table, and the design matrix
  topProbes <- as.character(diffs[1:numProbes,]$ID)
  exprsData <- exprs(earlyLateESet[top500,])
  colColours <- sub(1,"orange",design[,1]) # assign all early samples to gray
  colColours <- sub(0,"blue",colColours) # assign all late samples to purple 

  heatmap.2(exprsData,
            main = titleString, # heat map title
            notecol="black",      # change font color of cell labels to black
            density.info="none",  # turns off density plot inside color legend
            trace="none",         # turns off trace lines inside the heat map
            margins =c(12,9),    # widens margins around plot
            col=greenred(10),
            ColSideColors = colColours,
            lhei = c(1,3)
            #Colv = FALSE
            #breaks=col_breaks,    # enable color transition at specified limits
  )
  par(lend = 1)           # square line ends for the color legend
  legend("topright",      # location of the legend on the heatmap plot
         legend = c("Early","Late"), # category labels
         col = c("orange", "blue"),  # color key
         lty= 1,             # line style
         lwd = 10            # line width
  )
}


titleString <- "Early Vs. Late Expression in Group B"
createHeatmap(diffsB,earlyLate.eset,design,titleString)

titleString <- "Early Vs. Late Expression in Group A"
createHeatmap(diffsA,earlyVsLateESet_A,designA,titleString)
