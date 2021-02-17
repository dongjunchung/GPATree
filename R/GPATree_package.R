#' GPATree: A package to implement the GPA-Tree method
#'
#' This package provides functions for fitting GPA-Tree, a statistical approach for integrative analysis of genome wide association studies (GWAS) data and functional annotation information within a unified framework. GPA-Tree simultaneously identifies disease risk-associated SNPs and combinations of functional annotations that potentially explain the mechanisms through which risk-associated SNPs are related with phenotypes.
#' @details
#' \itemize{
#' \item Package: GPATree
#' \item Type: Package
#' \item Version: 0.0.0.9000
#' \item Date: 2021-02-16
#' \item License: GPL(>=3)
#' \item LazyLoad: yes
#' }
#' This package contains a main class, GPATree, which represents GPATree model fit. This package contains five main methods for the GPATree framework, GPATree, plot, assoc, prune, leaf. GPATree method fits the GPATree model and assoc method implements association mapping. leaf provided information regarding functional annotations that are enriched for the leaves in the GPATree model results. plot allows plotting the GPATree model result, and prune allows further pruning the GPATree.
#' This package also contains a methods for the ShinyGPATree visualization, association mapping and functional annotation tree selection toolkit. ShinyGPATree opens the ShinyGPATree interface, which takes the results generated from GPATree method as input.
#' @seealso GPATree, assoc, leaf, prune, plot, ShinyGPATree
#' @author  Aastha Khatiwada \email{asthakhatiwada@gmail.com}
#' @examples
#' \dontrun{
#' library(GPATree)
#'
#' # load GPATree example data
#' data(GPATreeExampleData)
#'
#' # fitting the GPATree model
#' fit <- GPATree(GPATreeExampleData$gwasPval, GPATreeExampleData$annMat)
#'
#' # get functional annotation information
#' leaf(fit)
#'
#' # association mapping
#' assoc.gpatree <- assoc(fit, FDR = 0.01, fdrControl = 'global')
#'
#' # pruning the GPATree model fit
#' pruned.fit <- prune(fit, cp = 0.005)
#'
#' # plotting the GPATree model results
#' plot(fit)
#' plot(pruned.fit)
#'
#' # run the ShinyGPATree app using output from the GPATree method
#' ShinyGPATree(fit)
#' }
#' @name GPATree-package
#' @aliases GPATree-package,GPATree-method
#' @docType package
"_PACKAGE"
utils::globalVariables("quantile_reg_model")
