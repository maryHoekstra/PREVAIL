library(gplots)
library(RColorBrewer)

earlyDays <- c("day 1")
lateDays <- c("day 14")
earlyVsLateESet <- getESet(filtered.eset,earlyDays,lateDays)
design <- createDesign(earlyVsLateESet)
diffsB <- getDiffs(earlyVsLateESet,groupCode="B",design)
groupB.eset = earlyVsLateESet[,earlyVsLateESet$code=="B"]
# use top 500 differentially expressed genes from top table result
topProbes <- as.character(diffsB[1:500,]$ID)
exprsData <- exprs(groupB.eset[topProbes,])
titleString <- "Day 1 vs. Day 14 in Group B"
dayLabels <- sub("day 14","firebrick1",groupB.eset$Day)
dayLabels <- sub("day 1","lightslateblue",dayLabels)

heatmap.2(exprsData,
          main = titleString, # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),    # widens margins around plot
          col = brewer.pal(9,"YlGnBu"),
          ColSideColors = dayLabels,
          #Colv = NA,
          #key = FALSE,
          keysize = 0.9,
          labRow = ""
)
coords <- locator(1)
par(lend = 1)           # square line ends for the color legend
legend(x = 0.815,y = 1.21, xpd = TRUE,      # location of the legend on the heatmap plot
       legend = c("Day 1","Day 14"), # category labels
       col = c("lightslateblue", "firebrick1"),  # color key
       lty= 1,             # line style
       lwd = 10            # line width
)


# use superheat package

# order samples by timepoint
orderedExprs <- exprs(groupB.eset[,groupB.eset$Day=="day 1"])
days = c("day 14","day 21","day 28")
dayAssgns <- rep("purple",times=ncol(orderedExprs))
for (i in 1:length(days)) {
  dayExprs <- exprs(groupB.eset[,groupB.eset$Day==days[i]])
  orderedExprs <- cbind(orderedExprs,dayExprs)
  dayAssgns <- c(dayAssgns,rep("gray",times=ncol(dayExprs)))
}

# create vector of colours indicated which day each sample belows to
dayLabels <- sub("day 14","salmon",groupB.eset$Day)
dayLabels <- sub("day 1","cornflowerblue",dayLabels)

# scale rows before passing to superheat
scaledData <- scale(t(exprsData),center = TRUE,scale = TRUE)
exprsData <- t(scaledData)

library(superheat)
superheat(exprsData,
          clustering.method = "hierarchical",
          pretty.order.rows = TRUE,
          pretty.order.cols = TRUE,
          left.label.size = 0.05,
          bottom.label.size = 0.1,
          col.dendrogram = TRUE,
          row.dendrogram = TRUE,
          row.title = "Probes",
          column.title = "Samples",
          left.label.text.size = 2,
          bottom.label.text.size = 2,
          bottom.label.text.angle = 90,
          legend.breaks = c(5,7.5, 10,12.5, 15),
          title = "Day 1 vs. Day 14 in Group B",
          bottom.label.col = dayLabels,
          #membership.cols = dayAssgns,
          grid.vline = FALSE)
