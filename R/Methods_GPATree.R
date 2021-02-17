#' @import graphics
#' @import methods
#' @importFrom rpart.plot rpart.plot
# generic methods for "GPATree" class
setMethod(
  f="show",
  signature="GPATree",
  definition=function( object ) {

    options(digits = 4)
    # summary of GGPA fit

    # constants

    M <- nrow(object@gwasPval)
    nAnn <- ncol(object@annMat)
    nGWAS<- ncol(object@gwasPval)
    l <- leaf(object)
    # estimates

    est_alpha = object@fit$alpha
    cat( "Summary: GPATree model results (class: GPATree)\n" )
    cat( "--------------------------------------------------\n" )
    cat( "Data summary:\n" )
    cat( "\tNumber of GWAS data: ", nGWAS , "\n", sep="" )
    cat( "\tNumber of Annotations: ", nAnn , "\n", sep="" )
    cat( "\tNumber of SNPs: ", M , "\n", sep="" )
    if(nAnn!=0) {
    cat( "\tAlpha estimate: ", est_alpha , "\n", sep="" )
    cat( "Functional annotation tree description:\n" )
    if (nrow(l) <= 15){
      print( l )
    } else {
      cat("\tNumber of leaves (terminal nodes): ", nrow(l), "\n", sep="")
    }
    } else {}
    cat( "--------------------------------------------------\n" )
  }
)

