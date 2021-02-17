#' Association mapping
#'
#' This function will implement association mapping for the GPA-Tree model.
#'
#' @author  Aastha Khatiwada
#'
#' @param object An object of class GPATree.
#' @param FDR FDR level. Value has to be between 0 and 1.
#' @param fdrControl Method to control FDR. Possible values are "global" (global FDR control) and "local" (local FDR control).
#'
#' @return Returns a binary matrix indicating the association between SNPs and phenotype, where its rows match those of input pvalue matrix for function GPATREE.
#' @export
#' @name assoc
#' @aliases assoc,GPATree-method
setMethod(
  f="assoc",
  signature="GPATree",
  definition=function( object, FDR=0.05, fdrControl="global") {
    
    if ( fdrControl != "global" & fdrControl != "local" ) {
      stop( "Invalid value for 'fdrControl' argument! It should be either 'global' or 'local'." )
    }
    
    if ( FDR < 0 | FDR > 1 ) {
      stop( "Invalid value for 'FDR' argument! FDR should be between zero and one." )
    }
    
    # association mapping
    fdrmat <- 1 - object@fit$Zmarg
    amat <- matrix( 0, nrow(fdrmat), ncol(fdrmat))
    colnames(amat) <- paste('P', 1:ncol(fdrmat), sep = '')
    rownames(amat) <- rownames(object@fit$Zmarg)
    
    if ( fdrControl == "local" ) {
      message( "Info: Association mapping based on the local FDR control at level ", FDR, "." )
      amat[ fdrmat <= FDR ] <- 1
    } else if ( fdrControl == "global" ) {
      message( "Info: Association mapping based on the global FDR control at level ", FDR, "." )
      # direct posterior probability approach for FDR control (Newton et al.)
      for ( j in 1:ncol(amat) ) {
        pp <- fdrmat[,j]
        pp.ordered <- sort(pp)
        pp.cum <- cumsum( pp.ordered ) / c(1:length(pp))
        cutoff <- max( pp.ordered[ pp.cum <= FDR ], 0 )
        amat[ pp <= cutoff, j ] <- 1
      }
    }
    
    # leaf mapping ####
    if (ncol(amat) >= 1){
      amat <- as.data.frame(amat)
      leaf_frame <- object@fit$fit$frame[object@fit$fit$frame$var == '<leaf>', c('var', 'n')]
      leaf_frame$var <- paste('LEAF', 1:nrow(leaf_frame))
      tree_rules <- decTree(object)
      leaf_frame$rule <- rep(NA, nrow(leaf_frame))
      for (i in 1:length(tree_rules$CART_PIs)) { leaf_frame$rule[i] <- tree_rules$CART_PIs[[i]] }
      whichnodes <- rownames(object@fit$fit$frame)[object@fit$fit$where]
      amat$leaf <- leaf_frame$var[match(whichnodes, rownames(leaf_frame))]
    }
    
    return(amat)
  
    }

)

