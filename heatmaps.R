library(gplots)

# creates heatmap using the top differentially expressed probes and a design matrix
createHeatmap <- function(diffs,earlyLateESet,design,titleString) {

  # take max of 500 top differentially expressed probes
  if (nrow(diffs) < 500) {
    numProbes=nrow(diffs)
  }
  else {
    numProbes=500
  }
  
  # need the subsetted expression matrix, the diffs table, and the design matrix
  topProbes <- as.character(diffs[1:numProbes,]$ID)
  exprsData <- exprs(earlyLateESet[topProbes,])
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
  )
  par(lend = 1)           # square line ends for the color legend
  legend("topright",      # location of the legend on the heatmap plot
         legend = c("A","B"), # category labels
         col = c("orange", "blue"),  # color key
         lty= 1,             # line style
         lwd = 10            # line width
  )
}


titleString <- "Early Vs. Late Expression in Group B"
createHeatmap(diffsB,earlyLate.eset,design,titleString)

titleString <- "Early Vs. Late Expression in Group A"
createHeatmap(diffsA,earlyVsLateESet_A,designA,titleString)

titleString <- "Day 1 vs. Day 7 in Group B"
createHeatmap(diffsB,earlyVsLateESet_B,designB,titleString)

# use superheat package

# order samples by timepoint
# given an expression set and a vector of ordered day strings, return the expression matrix with samples in order
groupB.eset = filtered.eset[,filtered.eset$code=="B"]
groupA.eset = filtered.eset[,filtered.eset$code=="A"]

#days <- as.character(group.eset$Day)
exprsData <- exprs(groupB.eset[topProbes,])

orderedExprs <- exprs(groupA.eset[,groupA.eset$Day=="day 1"])
days = c("day 14","day 21","day 28")
dayAssgns <- rep("day 1",times=ncol(orderedExprs))
for (i in 1:length(days)) {
  dayExprs <- exprs(groupA.eset[,groupA.eset$Day==days[i]])
  orderedExprs <- cbind(orderedExprs,dayExprs)
  dayAssgns <- c(dayAssgns,rep(days[i],times=ncol(dayExprs)))
}

# use top 500 differentially expressed genes from top table result
topProbes <- as.character(diffsB[1:500,]$ID)
exprsData <- orderedExprs[topProbes,]

library(superheat)
superheat(exprsData,
          pretty.order.rows = TRUE,
          pretty.order.cols = TRUE,
          left.label.size = 0.05,
          bottom.label.size = 0.05,
          col.dendrogram = TRUE,
          clustering.method = 'hierarchical',
          row.dendrogram = TRUE,
          row.title = "Probes",
          column.title = "Samples",
          left.label.text.size = 2,
          bottom.label.text.size = 2,
          bottom.label.text.angle = 90,
          legend.breaks = c(5,7.5, 10,12.5, 15),
          title = "Early vs. Late in Group B",
          #membership.cols = dayAssgns,
          grid.vline = FALSE)
