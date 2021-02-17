#' Functional annotation tree.
#'
#' This function will provide the annotation combinations relevant to risk-associated SNPs.
#'
#' @author  Aastha Khatiwada
#'
#' @param object An object of class GPATree.
#'
#' @return Returns a matrix where each row corresponds to a leaf from the GPA-Tree model fit and contains information regarding the local FDR for SNPs in the leaf, and also information regarding annotations that are enriched (1) or not enriched (0) for the leaf.
#' @examples
#' \dontrun{
#' library(GPATree)
#'
#' # load GPATree example data
#' data(GPATreeExampleData)
#'
#' #fitting the GPATree model
#' fit <- GPATree(GPATreeExampleData$gwasPval, GPATreeExampleData$annMat)
#' leaf(fit)
#' }
#' @export
#' @name leaf
#' @aliases leaf,GPATree-method
setMethod(
  f="leaf",
  signature="GPATree",
  definition=function( object ) {

    options(digits = 4)
    CARTmod <- object@fit$fit

    if ( nrow(CARTmod$frame) == 1){
      if (ncol(CARTmod$frame) >= 8 ){
        mean_localFDR <- colMeans(1-object@fit$Zmarg, na.rm = TRUE)
        select_annMat <- as.data.frame(matrix(c(mean_localFDR, "No annotations selected"), nrow = 1, ncol = length(mean_localFDR) + 1))
        rownames(select_annMat) <- 'LEAF 1'
        colnames(select_annMat) <- c('local FDR', "Note")
        select_annMat$`local FDR` <- as.numeric(select_annMat$`local FDR`)
        }
    }

    if ( nrow(CARTmod$frame) > 1){

      term_nodes <- as.numeric( rownames(CARTmod$frame)[c( which(CARTmod$frame[, 1] == "<leaf>") )  ] )
      pth <- rpart::path.rpart(tree = CARTmod, nodes = term_nodes, pretty = 0, print.it = F) #list of the PIs for CART
      anns <- setdiff(unlist(pth), 'root')
      parts <- strsplit(anns, "=")
      select_ann <- rep(NA, length(parts))

      for (i in 1:length(parts)) { select_ann[i] <- parts[[i]][1] }

      select_ann <- unique(select_ann)
      select_annMat <- matrix(NA, nrow = length(pth), ncol = length(select_ann)+ncol(object@fit$Zmarg))

      colnames(select_annMat) <- c('local FDR', select_ann)

      rownames(select_annMat) <- paste('LEAF ', 1:length(pth), sep = '')
      cartmod_frame_ind <- sort(unique(CARTmod$where))

      if (ncol(CARTmod$frame) >= 8 ){

        for (i in 1:length(pth))  {
          select_annMat[i, 1:ncol(object@fit$Zmarg)] <- round(colMeans(1-as.matrix(object@fit$Zmarg[which(CARTmod$where == cartmod_frame_ind[i]), ]), na.rm = TRUE), 4)
          curr_PI <- setdiff(pth[[i]], "root") #PI in CART format

          for (j in 1:length(curr_PI)) {#NOTE 1st element always =n root
            parts <- strsplit(curr_PI[j], "=")
            select_annMat[i, which(colnames(select_annMat) == parts[[1]][1])] <- as.numeric(parts[[1]][2])
          }
        }
      }

    }

    select_annMat <- as.data.frame(select_annMat)
    select_annMat$`local FDR` <- round(select_annMat$`local FDR`, 4)
    select_annMat[is.na(select_annMat)] <- '-'

    return(select_annMat)

  }
)

