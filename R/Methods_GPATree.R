# generic methods for "GPATree" class

# GPATree model fit summary

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
    cat( "Functiona annotation tree description:\n" )
    print( leaf(object) )
    } else {}
    cat( "--------------------------------------------------\n" )
  }
)

# GPATree plot
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
      colnames(out$frame)[ncol(out$frame)] <- paste('local FDR P', j, sep = '')
    }
    
    node.fun1 <- function(x, labs, digits, varlen) {
      
      p1 <- paste(out$frame$leaf, "\n \n N = ", out$frame$n, sep = '')
      
      if (ncol(out$frame) > 10) {
        for (i in 1:ncol(x@fit$Zmarg)) {
          p1 <- paste(p1,"\n \n local FDR P", i, " = ", round(out$frame[, 10+i], 3), sep = '')
        }
      }
      
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

