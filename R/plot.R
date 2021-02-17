#' Plot the functional annotation tree
#'
#' This function will plot the functional annotation tree for the GPA-Tree model fit.
#'
#' @author  Aastha Khatiwada
#'
#' @param x An object of class GPATree.
#' @param y missing (not required).
#' @param ... ...
#'
#' @return Returns a plot for the functional annotation tree from the GPA-Tree model fit.
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
#' # plotting the GPATree model fit
#' plot(fit)
#' }
#' @name plot
#' @aliases plot,GPATree,missing-method
#' @export
setMethod(
  f="plot",
  signature=c(x = "GPATree", y = "missing"),
  definition=function( x, y, ... ) {

    out <- x@fit$fit
    leaves_index <- which(out$frame$var=='<leaf>')
    out$frame$leaf <- rep(NA, nrow(out$frame))
    out$frame$leaf[leaves_index] <- paste('LEAF ', 1:length(leaves_index), sep = '')
    localFDR <- vector("list", ncol(x@fit$Zmarg))

    for (j in 1:length(localFDR)) {

      localFDR[[j]] <- rep(NA, nrow(out$frame))
      for (i in 1:length(leaves_index)) {
        SNPindex <- which(out$where == leaves_index[i])
        localFDR[[j]][leaves_index[i]] <- round(mean(1 - x@fit$Zmarg[SNPindex, j], na.rm = TRUE), 4)
      }
      out$frame$v1 <- localFDR[[j]]
      colnames(out$frame)[ncol(out$frame)] <- paste('local FDR')
    }

    node.fun1 <- function(x, labs, digits, varlen) {

      p1 <- paste(out$frame$leaf, "\n \n N = ", out$frame$n, sep = '')

      if (ncol(out$frame) == 10) {
        p1 <- paste(out$frame$leaf,
                    "\n \n N = ",
                    out$frame$n,
                    "\n \n local FDR = ",
                    formatC(round(out$frame$`local FDR`, 3), 3, format = 'f'),
                    sep = '')
      }
      p1
    }
    rpart.plot::rpart.plot(out,
                           node.fun = node.fun1,
                           type = 5,
                           extra = 1,
                           font = 2,
                           split.font = 2,
                           col = 'black',
                           under.col=10,
                           compress = TRUE,
                           ycompress = TRUE,
                           roundint=FALSE)
  }

)

