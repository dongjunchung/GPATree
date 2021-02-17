#' GPATree selected decision tree
#'
#' This function will extract the combinations of functional annotations selected by the decision tree in the stage 2 of GPATree method that meets the provided threshold in minPredictedProb.
#'
#' @author  Aastha Khatiwada
#' @param object An object of class GPATree.
#' 
#' @return A list containing variables in combinations selected by the decision tree (CART_PIs), the combination (CART_PIs_comb) and the predicted proportions for the selected PIs(assoc_pred) meet the provided threshold in minPredictedProb.
#' @import rpart
#' @importFrom dplyr mutate
#' @importFrom dplyr %>%
#' @name decTree
#' @aliases decTree,GPATree-method
setMethod(
  f="decTree",
  signature="GPATree",
  definition=function( object ) {
    
  
  CARTmod <- object@fit$fit
  minPredictedProb <- 0
  
  if ( nrow(CARTmod$frame) == 1){
    if (ncol(CARTmod$frame) >= 8 ){
      yvals <- CARTmod$frame$yval
      ans <- list(CART_PIs = 'none',
                  CART_PIs_comb = 'none',
                  assoc_pred = CARTmod$frame$yval)
    }
  }
  
  
  if ( nrow(CARTmod$frame) > 1){
    
    term_nodes <- as.numeric( rownames(CARTmod$frame)[c( which(CARTmod$frame[, 1] == "<leaf>") )  ] )
    pth <- rpart::path.rpart(tree = CARTmod, nodes = term_nodes, pretty = 0, print.it = F) #list of the PIs for CART
    
    if (ncol(CARTmod$frame) >= 8 ){
      yvals <- (CARTmod$frame$yval)[c(which(CARTmod$frame[, 1] == "<leaf>"))]
      y_ids <- 1:length(yvals)
      yvals <- yvals[c(y_ids)]
      num_CARTpis <- length(yvals) #number of PIs identified by CART
      CART_PIs <- vector("list", num_CARTpis)
      
      for (i in 1:num_CARTpis)  {
        curr_PI <- setdiff(pth[[y_ids[i]]], "root") #PI in CART format
        elems <- length(curr_PI) #number elements in current PI
        CPM_nm <- c()
        
        for (j in 1:elems) {#NOTE 1st element always =n root
          parts <- strsplit(curr_PI[j], "=")
          part_which <- parts[[1]][2]
          
          if (part_which == '0') {#Note these variables are compliments
            el <- strsplit(curr_PI[j], "=")[[1]][1]
            elem <- paste("!", el, sep = "")
            CPM_nm <- append(CPM_nm, elem)
          }
          
          if (part_which == '1') {
            elem <- strsplit(curr_PI[j], "=")[[1]][1]
            CPM_nm <- append(CPM_nm, elem)
          }
        }
        
        CART_PIs[i] <- paste('(',  paste(CPM_nm, collapse = ' & ', sep = ''), ')', sep = '')
        
      }
      
      paste(unlist(CART_PIs), collapse = ' OR ')
      
      
      # combining only selected combinations
      yvals_select <- which(yvals > minPredictedProb)
      
      if (length(yvals_select) > 0){
        comb_select <- paste(unlist(CART_PIs[yvals_select]), collapse = ' OR ')
      } else {
        comb_select <-  'none'
      }
      ans <- list(CART_PIs = CART_PIs[yvals > minPredictedProb],
                  CART_PIs_comb = comb_select,
                  assoc_pred = yvals[yvals > minPredictedProb])
    }
  }
  
  return(ans)
  
  }

  )

