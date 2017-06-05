library(VennDiagram)
library(gdata)

GO_BP_groupA = read.xls ("/Users/maryhoekstra/Desktop/GO_BP_EarlyVsLate.xlsx", sheet = 1, header = TRUE)
GO_BP_groupB = read.xls ("/Users/maryhoekstra/Desktop/GO_BP_EarlyVsLate.xlsx", sheet = 2, header = TRUE)
GO_BP_groupA_down = read.xls ("/Users/maryhoekstra/Desktop/GO_BP_EarlyVsLate_regulated.xlsx", sheet = 1, header = TRUE)
GO_BP_groupB_down = read.xls ("/Users/maryhoekstra/Desktop/GO_BP_EarlyVsLate_regulated.xlsx", sheet = 2, header = TRUE)
GO_BP_groupA_up = read.xls ("/Users/maryhoekstra/Desktop/GO_BP_EarlyVsLate_regulated.xlsx", sheet = 3, header = TRUE)
GO_BP_groupB_up = read.xls ("/Users/maryhoekstra/Desktop/GO_BP_EarlyVsLate_regulated.xlsx", sheet = 4, header = TRUE)

GO_BP_groupB_incDay7 = read.xls("/Users/maryhoekstra/Desktop/GO_BP_EarlyVsLate.xlsx", sheet = 3, header = TRUE)
GO_BP_groupA_day1vs7 = read.xls("/Users/maryhoekstra/Desktop/GO_BP_EarlyVsLate.xlsx", sheet = 4, header = TRUE)
GO_BP_groupB_day1vs7 = read.xls("/Users/maryhoekstra/Desktop/GO_BP_EarlyVsLate.xlsx", sheet = 5, header = TRUE)
GO_BP_groupA_day1vs3 = read.xls("/Users/maryhoekstra/Desktop/GO_BP_EarlyVsLate.xlsx", sheet = 6, header = TRUE)
GO_BP_groupB_day1vs3 = read.xls("/Users/maryhoekstra/Desktop/GO_BP_EarlyVsLate.xlsx", sheet = 7, header = TRUE)

createVenn <- function(groupA,groupB,titleString,groupAstring,groupBstring) {
  groupATerms <- groupA$Term
  groupBTerms <- groupB$Term
  commonTerms <- intersect(groupATerms,groupBTerms)
  
  # find common pathways and those unique to each group
  commonTerms <- intersect(groupATerms,groupBTerms)
  groupATerms_unique <- setdiff(groupATerms,groupBTerms)
  groupBTerms_unique <- setdiff(groupBTerms,groupATerms)
  
  # take top 6 pathways
  topCommonTerms <- head(commonTerms)
  topGroupATerms_unique <- head(groupATerms_unique)
  topGroupBTerms_unique <- head(groupBTerms_unique)
  
  # strip GO tag prefix
  strippedCommon <- unlist(lapply(topCommonTerms,function(string) sub(".*~", "", string)))
  strippedTermsA <- unlist(lapply(topGroupATerms_unique,function(string) sub(".*~", "", string)))
  strippedTermsB <- unlist(lapply(topGroupBTerms_unique,function(string) sub(".*~", "", string)))
  
  # create venn diagram object
  v <- venn.diagram(
    x = list(groupATerms,groupBTerms),
    category.names = c(groupAstring,groupBstring),
    filename = NULL,
    fill = c('yellow', 'purple'),
    main=titleString,
    scaled = FALSE,
    fontface="bold",
    cat.fontface = "bold",
    main.fontface = "bold",
    main.cex = 2,
    cat.cex = 2,
    cex = 0.8
  )
  
  # over-write labels (5 to 7 chosen by manual check of labels)
  v[[5]]$label  <- paste(strippedTermsA, collapse="\n\n") 
  v[[6]]$label <- paste(strippedTermsB, collapse="\n\n") 
  v[[7]]$label <- paste(strippedCommon, collapse="\n\n") 
  
  # plot venn 
  grid.newpage()
  grid.draw(v)
  
}

createVenn(GO_BP_groupA_down,GO_BP_groupB_down,"Upregulated pathways in late samples","A","B")
createVenn(GO_BP_groupA_up,GO_BP_groupB_up,"Downregulated pathways in late samples","A","B")

createVenn(GO_BP_groupB,GO_BP_groupB_incDay7,"Differentially expressed pathways","W/O Day 7", "W/ Day 7")
createVenn(GO_BP_groupB,GO_BP_groupB_incDay7,"Differentially expressed pathways","W/O Day 7", "W/ Day 7")

createVenn(GO_BP_groupA_day1vs7,GO_BP_groupB_day1vs7,"Differentially expressed pathways between baseline and day 7","A","B")
createVenn(GO_BP_groupA_day1vs3,GO_BP_groupB_day1vs3,"Differentially expressed pathways between baseline and day 3","B","A") # swap labels since A is empty
