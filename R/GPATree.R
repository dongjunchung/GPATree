#' Implement the GPA-Tree approach
#'
#' This function will implement the GPA-Tree approach.
#'
#' @author  Aastha Khatiwada
#'
#' @param gwasPval A matrix of M X 1 dimension, where M is the number of SNPs. The matrix includes the GWAS association p-values for the phenotype. P-values must be between 0 and 1.
#' @param annMat A matrix of binary annotations, where row and column correspond to SNPs and annotations, respectively.
#' @param initAlpha Initial value for alpha estimate. Default is 0.1.
#' @param cpTry Complexity parameter (cp) value to be used. cpTry can be between 0 and 1 or NULL. Default is 0.001. When cpTry is NULL, GPATree will select the optimal cp to be used.
#'
#' @return This function returns a \code{List} including:
#' \itemize{
#' \item numIterConvergence: number of iterations taken for GPA-Tree approach to converge in Stage 1 and 2.
#' \item alpha: alpha estimated by the GPA-Tree approach.
#' \item fit: GPA-Tree model fit.
#' \item fitSelectVar: annotations included in the tree.
#' \item Z: posterior probability of being null and non-null for all SNPs.
#' \item zMarg: posterior probability of being a non-null SNP.
#' \item pi: updated prior based on GPA-Tree model.
#' \item licVec: vector of incomplete log-likelihood for all iteration of GPA-Tree approach.
#' \item lcVec: list of complete log likelihood for all iteration of GPA-Tree approach.
#' \item annMat: annMat provided by user.
#' \item gwasPval: gwasPval provided by user.
#' }
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
