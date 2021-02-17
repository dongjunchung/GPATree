#' Implement the Stage 2 of the GPA-Tree approach.
#'
#' This function will implement the Stage 2 of the GPA-Tree approach.
#'
#' @author  Aastha Khatiwada
#'
#' @param gwasPval A matrix of M X 1 dimension, where M is the number of SNPs. The first column is labeled 'SNPid' and contains the SNPid. The second column contains the GWAS association p-values and is called 'P1'. The values in P1 must be between 0 and 1.
#' @param annMat A matrix of binary annotations, where row and column correspond to SNPs and annotations, respectively.
#' @param alphaStage1 alpha estimated in stage 1 of the GPA-Tree approach.
#' @param initPi pi estimated in stage 1 of the GPA-Tree approach.
#' @param cpTry Complexity parameter (cp) value to be used to build annotation decision tree. cpTry can be between 0 and 1 or NULL. Default is 0.001. When cpTry is NULL, GPATree will select the optimal cp to be used.
#'
#' @return This function returns a \code{List} including:
#' \itemize{
#' \item num_iter_stage2: number of iterations taken for Stage 2 of GPA-Tree approach to converge.
#' \item lic_list_stage2: list of incomplete log-likelihood for all iterations in Stage 2 of GPA-Tree approach.
#' \item lc_list_stage2: list of complete log-likelihood for all iterations in Stage 2 of GPA-Tree approach.
#' \item zi_stage2: probability of being a null  and non-null SNP.
#' \item Zmarg: marginal posterior probability of being a non-null SNP for the phenotype.
#' \item pi_stage2: updated prior.
#' \item var_in_tree: annotations included in the tree.
#' \item pruned_tree: GPA-Tree model fit.
#' }
#' @import utils
#' @import readr
#' @importFrom rpart rpart
#' @importFrom rpart prune
#' @importFrom rpart.plot rpart.plot
#' @import quantreg
#' @importFrom rpart.utils rpart.subrules.table
#' @importFrom gtools permutations


GPATreeStage2 <- function(gwasPval, annMat, alphaStage1, initPi, cpTry){
  
  message( "Info: Processing Stage 2 of GPA-Tree...")
  
  yi <- as.vector(gwasPval[, 1])
  names(yi) <- rownames(gwasPval)
  M <- length(yi)
  nvar <- ncol(annMat)

  # percent of annotated SNPs
  est_prop_ones_all_ann <- colMeans(annMat)
  est_percent_ones <- mean(est_prop_ones_all_ann)

  # percent of overlap between annotated SNPs
  est_prop_overlap_all <- apply(gtools::permutations(nvar, 2, 1:nvar), 1, function(cm) {prop.table( table( annMat[,cm[1]], annMat[,cm[2]] ), 1 )[4]})
  names(est_prop_overlap_all) <- apply(gtools::permutations(nvar, 2, 1:nvar), 1, function(cm) {paste("A", cm, sep="", collapse = ',')})
  est_percent_overlap <- mean(est_prop_overlap_all[est_prop_overlap_all>0])

  alpha <- alphaStage1
  pi <- initPi
  pi[ pi < 1e-10 ] <- 1e-10
  pi[ pi > 1 - 1e-10 ] <- 1 - 1e-10
  
  annMat <- as.data.frame(annMat)
  for (i in 1:ncol(annMat)) {  annMat[, i] <- as.factor(annMat[, i])   }

  # initial values for EM ####
  verbose <- TRUE
  iter_max <- 20000
  iter_cur <- 2

  # determining pruning parameters ####
  if (!is.null(cpTry)) {
    cp_qreg_start <- cp_prune <- cpTry
  }
  if (is.null(cpTry)) {
    # load the quantile regression model stored in the data folder
    data("quantile_reg_model", envir=environment())
    newdat_cp <- data.frame(cbind(alpha,
                                  percent.ones = est_percent_ones,
                                  percent.overlap = est_percent_overlap))
    cp_qreg_start <- cp_prune <- 10^(quantreg::predict.rq(quantile_reg_model, newdata = newdat_cp))
  }

  # cp try list if lc doesn't increase ####
  cp_log10_prune_list <- seq(-1, -5, by =-0.01) #suggested by Dr. Chung
  cp_log10_scale_list <- 10^cp_log10_prune_list

  # storing
  liclist <- rep( NA, iter_max )
  lclist <- rep(NA, iter_max)

  # initial stores
  liclist[1] <- -10000000
  lclist[1] <- -10000000

  # EM start ####
  
  options(digits = 20)

  while( iter_cur < iter_max) {


    # E Step
    zi <- (alpha* (yi^(alpha-1)) * pi ) / ( (alpha* (yi^(alpha-1)) * pi) + (1-pi) ) # posterior probability of having signal or being non-null for each SNP
    zi[ zi < 0.01 ] <- 0.01
    zi[ zi > 0.99 ] <- 0.99

    # M Step using regression tree 
    dat <- as.data.frame(cbind(annMat, 'zi'=zi)) # data to use in regression tree
    reg_full_tree <- rpart::rpart(zi ~ . , 
                                  method = 'anova', 
                                  data=dat, 
                                  control = rpart::rpart.control(cp=0))
    reg_tree <- rpart::prune(reg_full_tree, cp = cp_prune) # prune the tree

    # update pi
    pi <- predict(reg_tree, newdata = as.data.frame(annMat))
    pi[ pi < 1e-10 ] <- 1e-10
    pi[ pi > 1 - 1e-10 ] <- 1 - 1e-10


    # Incomplete log-Likelihood
    lic <- sum(log ( (pi * alpha* (yi^(alpha-1)) ) + (1-pi) ))


    # Complete data log-likelihood
    lc <- sum( (zi * ( log(pi) + log(alpha) + ((alpha - 1) * log(yi)) )) + ((1 - zi) * log(1-pi)) )

    # if necessary, repeat M Step to make sure lic[cur.iter] > lic[cur.iter-1]
    if (lic < liclist[iter_cur-1]) {

      # print('inside the loop')
      lc_list_from_try_cp <- rep(NA, length(cp_log10_scale_list))

      for (t in 1:length(cp_log10_scale_list)) {

        # try the different cp in the list of cp to try
        reg_tree_try <- rpart::prune(reg_full_tree, cp = cp_log10_scale_list[t] )

        # update pi
        pi_try <- predict(reg_tree_try, newdata = as.data.frame(annMat))

        # Complete data log-likelihood
        lc_list_from_try_cp[t] <- sum( (zi * ( log(pi_try) + log(alpha) + ((alpha - 1) * log(yi)) )) + ((1 - zi) * log(1-pi_try)) )
      }

      cp_prune_try <- cp_log10_scale_list[which(lc_list_from_try_cp == max(lc_list_from_try_cp))[1]]

      # using the cp that maximizes the lc, get new estimates for the parameters
      reg_tree_try <- rpart::prune(reg_full_tree, cp = cp_prune_try)

      # updates
      # pi
      pi_try <- predict(reg_tree_try, newdata = as.data.frame(annMat))
      # Incomplete log-Likelihood
      lic_try <- sum(log ( (alpha* (yi^(alpha-1)) * pi_try) + (1-pi_try) ))
      # Complete data log-likelihood
      lc_try <- sum( (zi * ( log(pi_try) + log(alpha) + ((alpha - 1) * log(yi)) )) + ((1 - zi) * log(1-pi_try)) )

      if (lic_try < liclist[iter_cur-1]){
        iter_cur <- iter_cur-1
        # print('EM converged (cannot find cp that increases incomplete data log-likelihood)')
        break } else {
          pi <- pi_try
          cp_prune <- cp_prune_try
          lic <- lic_try
          lc <- lc_try
          reg_tree <- reg_tree_try
        }
    }
    # store values
    liclist[iter_cur] <- lic
    lclist[iter_cur] <- lc

    # stopping rules
    if ( abs(liclist[iter_cur] - liclist[iter_cur-1]) <= 1e-3 ) {
      break
    }

    iter_cur <- iter_cur + 1

  }
  # one extra step to make sure final zi and pi do not have restrictions
  iter_cur <- iter_cur + 1
  pi <- predict(reg_tree, newdata = as.data.frame(annMat))
  zi <- (alpha* (yi^(alpha-1)) * pi ) / ( (alpha* (yi^(alpha-1)) * pi) + (1-pi) ) # posterior probability of having signal or being non-null for each SNP
  lic <- sum(log ( (pi * alpha* (yi^(alpha-1)) ) + (1-pi) ))
  liclist[iter_cur] <- lic
  lc <- sum( (zi * ( log(pi) + log(alpha) + ((alpha - 1) * log(yi)) )) + ((1 - zi) * log(1-pi)) )
  lclist[iter_cur] <- lc
  
  
  if ( nrow(reg_tree$frame) > 1){
    vars_in_tree <- as.character(unique(rpart.utils::rpart.subrules.table(reg_tree)$Variable))
    vars_in_tree <- paste(vars_in_tree[order(vars_in_tree)], sep = '', collapse = ", ")
  } else{
    vars_in_tree <- 'none'
  }
  
  # list of variables selected by new method
  if (length(reg_tree$variable.importance) >0){ # which annotations are selected?
    var_selected_importance <- paste(names(reg_tree$variable.importance), sep ='', collapse = ', ')
    var_selected_ordered <-  paste(names(reg_tree$variable.importance)[order(names(reg_tree$variable.importance))], sep ='', collapse = ', ')
  } else {
    var_selected_importance <- var_selected_ordered <- 'none'
  }
  message('Info: Incomplete log-likelihood converged in Stage 2 of GPA-Tree.')
  message( "Info: Stage 2 of GPA-Tree completed.")
  
  EM_stage2_out <- list()
  EM_stage2_out$num_iter_stage2 = iter_cur
  EM_stage2_out$lic_list_stage2 = liclist[1:iter_cur]
  EM_stage2_out$lc_list_stage2 = lclist[1:iter_cur]
  EM_stage2_out$zi_stage2 = cbind('0' = 1-zi, '1' =  zi)
  rownames(EM_stage2_out$zi_stage2) <- names(yi)
  EM_stage2_out$Zmarg = as.matrix(zi)
  rownames(EM_stage2_out$Zmarg) <- names(yi)
  EM_stage2_out$pi_stage2 = pi
  names(EM_stage2_out$pi_stage2) <- names(yi)
  EM_stage2_out$var_in_tree = vars_in_tree
  EM_stage2_out$pruned_tree = reg_tree

  return(EM_stage2_out)

}


