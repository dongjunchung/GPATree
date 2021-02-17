#' Prune GPA-Tree model fit
#'
#' This function will prune the GPA-Tree model fit using the given cp value.
#'
#' @author  Aastha Khatiwada
#'
#' @param object An object of class GPATree.
#' @param cp The cp parameter to be used for pruning. cp must be between 0 and 1.
#'
#' @return GPA-Tree model output.
#' @examples
#' \dontrun{
#' library(GPATree)
#'
#' # load GPATree example data
#' data(GPATreeExampleData)
#'
#' #fitting the GPATree model
#' fit <- GPATree(GPATreeExampleData$gwasPval, GPATreeExampleData$annMat)
#'
#' # pruning the GPATree model fit
#' pruned.fit <- prune(fit, cp = 0.005)
#' }
#' @importFrom rpart prune
#' @export
#' @name prune
#' @aliases prune,GPATree-method
setMethod(
  f="prune",
  signature="GPATree",
  definition=function(object, cp=0.001) {

    if ( cp < 0 | cp >1 ) {
      stop( "Inappropriate value for 'cp' argument. It should be between zero and one." )
    }


  pruned_tree <- object

  if (ncol(object@fit$Zmarg) == 1) {
    pruned_tree@fit$fit <- rpart::prune(object@fit$fit, cp = cp)
  }

  pruned_tree@fit$fitSelectVar <- pruned_tree@fit$fit$frame$var
  pruned_tree@fit$fitSelectVar <- pruned_tree@fit$fitSelectVar[ pruned_tree@fit$fitSelectVar != '<leaf>']
  pruned_tree@fit$fitSelectVar <- ifelse(length(pruned_tree@fit$fitSelectVar) > 1,
                                     paste(as.character(sort(unique(pruned_tree@fit$fitSelectVar))), collapse = ', '),
                                     ifelse(length(pruned_tree@fit$fitSelectVar) == 1,
                                            as.character(unique(pruned_tree@fit$fitSelectVar)),
                                            'none'))
  message( "Info: GPA-Tree model pruned based on the complexity parameter (cp) value of ", round(cp, 6), "." )

  return(pruned_tree)

  }

)
