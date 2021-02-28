#' Fit GPA-Tree model
#'
#' This function will implement the GPA-Tree approach for integrative analysis of GWAS and functional annotation data.
#'
#' @author  Aastha Khatiwada
#'
#' @param gwasPval A matrix of M X 1 dimension, where M is the number of SNPs. The matrix includes the GWAS association p-values for the phenotype. P-values must be between 0 and 1.
#' @param annMat A matrix of binary annotations, where rows and columns correspond to SNPs and annotations, respectively.
#' @param initAlpha Initial value for alpha estimate. Default is 0.1.
#' @param cpTry Complexity parameter (cp) value to be used. cpTry can be between 0 and 1 or NULL. Default is 0.001. When cpTry is NULL, GPATree will select the optimal cp to be used.
#' @details The GPATree() function fits the GPATree model. It requires to provide GWAS p-value to gwasPval and binary annotation data to annMat.
#' It is assumed that number of rows of matrix in gwasPval and annMat are equal and correspond to the same SNP.
#'
#' The assoc() function implements association mapping.
#'
#' The plot() function takes in an object of class GPATree and will plot the functional annotation tree from the GPATree model.
#'
#' The leaf() function takes in an object of class GPATree and will provide information regarding the functional annotations that are enriched (1) or not enriched (0) for SNPs in any leaf of the GPATree model plot.
#'
#' The prune() function takes in an object of class GPATree and a cp parameter and will prune the GPATree model result. This function can be useful when the tree obtained from GPATree model is huge.
#'
#' The ShinyGPATree app provides visualization of the GPA-Tree model, identifies risk-associated SNPs, and characterizes the combinations of functional annotations that can describe the risk-associated SNPs. The app can also be utilized to improve the visualization of the GPA-Tree model fit to collate or separate layers of the model (add or remove leaves).
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
#' @return Contructs a GPATree class object
#' @importFrom methods new
#' @export

GPATree <- function(gwasPval, annMat, initAlpha = 0.1, cpTry = 0.001){

  if ( initAlpha<= 0 | initAlpha >= 1 ) {
    stop( "Inappropriate value for 'initAlpha' argument. It should be between zero and one." )
  }

  if ( !is.matrix(gwasPval) ) {
    gwasPval <- as.matrix(gwasPval, colClasses=c("numeric"))
  }

  if ( !is.null(annMat) ) {
    if ( !is.matrix(annMat) ) {
      annMat <- as.matrix(annMat)
    }
  }

  if ( !is.null(annMat) ) {
    if ( nrow(gwasPval) != nrow(annMat) ) {
      stop( "Number of SNPs are different between p-value and annotation matrices. They should coincide. Please check your p-value and annotation matrices.")
    }
  }

  if (ncol(gwasPval) > 1){
    stop("Summary statistics for more than one GWAS. gwasPval should have only one column.")
  }

  if ( any( gwasPval < 0 | gwasPval > 1 ) ) {
    stop( "Some p-values are smaller than zero or larger than one. p-value should be ranged between zero and one. Please check your p-value matrix." )
  }

  if ( !is.null(annMat) ) {
    if ( any( annMat != 0 & annMat != 1 ) ) {
      stop( "Some elements in annotation matrix has values other than zero or one. All the elements in annotation matrix should be either zero (not annotated) or one (annotated). Please check your annotation matrix." )
    }
  }

  if (!is.null(cpTry)) {
    if ( cpTry< 0 | cpTry > 1 ) {
      stop("Inappropriate value for complexity parameter (cp). cp should be between 0 and 1")
    }
  }



  nGWAS <- ncol(gwasPval)

  if ( !is.null(annMat) ) {
    nAnn <- ncol(annMat)
  }

  message( "Info: Number of GWAS data: ", nGWAS )
  if ( !is.null(annMat) ) {
    message( "Info: Number of annotation data: ", nAnn )
  }


  if (!is.null(cpTry)) {
    message( "Info: complexity parameter (cp) to be used: ", cpTry )
  }
  if (is.null(cpTry)) {
    cpTry <- 0.001
    message( "Info: complexity parameter (cp) will be determined by GPA-Tree.")
  }

  if (nGWAS == 1){


    GPATREE_S1 <- GPATreeStage1(gwasPval = gwasPval,
                                annMat = annMat,
                                initAlpha = initAlpha)

    GPATREE_S2 <- GPATreeStage2(gwasPval = gwasPval,
                                annMat = annMat,
                                alphaStage1 = GPATREE_S1$alpha_stage1,
                                initPi = GPATREE_S1$pi_stage1,
                                cpTry = cpTry)

    return_list <- list()
    return_list$numIterConvergence <- c('stage1' = GPATREE_S1$num_iter_stage1, 'stage2' =  GPATREE_S2$num_iter_stage2)
    return_list$alpha <- GPATREE_S1$alpha_stage1
    names(return_list$alpha) <- 'alpha'
    return_list$fit <- GPATREE_S2$pruned_tree
    return_list$fitSelectVar <- GPATREE_S2$var_in_tree
    return_list$Z <- as.matrix(GPATREE_S2$zi_stage2)
    return_list$Zmarg <- as.matrix(GPATREE_S2$Zmarg)
    colnames(return_list$Zmarg) <- colnames(gwasPval)[1]
    rownames(return_list$Zmarg) <- rownames(gwasPval)
    return_list$pi <- matrix(GPATREE_S2$pi_stage2, ncol = 1)
    colnames(return_list$pi) <- colnames(gwasPval)[1]
    rownames(return_list$pi) <- rownames(gwasPval)
    return_list$licVec <- GPATREE_S2$lic_list_stage2
    return_list$lcVec <- GPATREE_S2$lc
    # return_list$annMat <- annMat
    # return_list$gwasPval <- gwasPval
  }

  new( Class = "GPATree",
       fit = return_list,
       gwasPval = gwasPval,
       annMat = annMat )

}
