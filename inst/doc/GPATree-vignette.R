## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  # tidy = TRUE,
  prompt = TRUE,
  # split=TRUE,
  # tidy = TRUE,
  comment = ""
)

## ----setup, include=TRUE, message=FALSE, warning=FALSE------------------------
library(GPATree)

## ----eval=TRUE, include=FALSE, echo = FALSE-----------------------------------
options(digits = 4)

## ----include=TRUE, message=FALSE, warning=FALSE, split=TRUE-------------------
data(GPATreeExampleData)
dim(GPATreeExampleData$gwasPval)
head(GPATreeExampleData$gwasPval)
dim(GPATreeExampleData$annMat)
head(GPATreeExampleData$annMat)

## ----GPATree, eval=TRUE, include=TRUE, message=FALSE, warning=FALSE, message=FALSE----
fit.GPATree <- GPATree(gwasPval = GPATreeExampleData$gwasPval,
                       annMat = GPATreeExampleData$annMat,
                       initAlpha = 0.1,
                       cpTry = 0.005)

## ---- eval=TRUE, include=FALSE, echo = FALSE----------------------------------
options(digits = 4)

## ---- eval=TRUE, include=TRUE, message=FALSE, warning=FALSE-------------------
fit.GPATree

## ---- eval=FALSE, include=TRUE, message=FALSE, warning=FALSE, tidy=TRUE, comment=""----
#  ShinyGPATree(fit.GPATree)

## ---- eval=FALSE, include=TRUE, message=FALSE, warning=FALSE------------------
#  fit.GPATree.pruned <- prune(fit.GPATree, 0.001)

## ---- eval=FALSE, include=TRUE, message=FALSE, warning=FALSE------------------
#  plot(fit.GPATree)

## ---- eval=TRUE, include=TRUE, message=FALSE, warning=FALSE-------------------
leaf(fit.GPATree)

## ---- eval=TRUE, include=TRUE, message=FALSE, warning=FALSE-------------------
assoc.SNP.GPATree <- assoc(fit.GPATree, 
                           FDR = 0.05, 
                           fdrControl="global")
head(assoc.SNP.GPATree)
table(assoc.SNP.GPATree$P1)
table(assoc.SNP.GPATree$leaf)
table(assoc.SNP.GPATree$P1, assoc.SNP.GPATree$leaf)

