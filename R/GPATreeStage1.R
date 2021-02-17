#' Implement the Stage 1 of the GPA-Tree Method
#'
#' This function will implement the stage 1 of the GPA-Tree method.
#'
#' @author  Aastha Khatiwada
#'
#' @param gwasPval A matrix of M X 1 dimension where M is the number of SNPs. The matrix contains GWAS association p-values. Values must be between 0 and 1.
#' @param annMat A matrix of binary annotations, where row and column correspond to SNPs and annotations, respectively.
#' @param initAlpha Initial value for alpha estimate. Default is 0.1.
#'
#' @import stats
#' @return This function returns a \code{List} including:
#' \itemize{
#' \item num_iter_stage1: number of iterations taken for Stage 1 of GPA-Tree to converge.
#' \item pi_stage1: predicted posterior probability of being a non-null SNP in Stage 1 of GPA-Tree Method.
#' \item alpha_stage1: estimated alpha of GPA-Tree Method.
#' \item beta_stage1: beta parameters from the linear model fitted at convergence of Stage 1 of GPA-Tree Method.
#' }

GPATreeStage1 <- function(gwasPval, annMat, initAlpha=0.1){


  message( "Info: Processing Stage 1 of GPA-Tree...")
  yi <- as.vector(gwasPval[, 1])
  names(yi) <- rownames(gwasPval)
  M <- length(yi)
  nAnn <- ncol(annMat)

  # initial values for EM ####
  verbose <- TRUE
  iter_max <- 20000
  iter_cur <- 2
  alpha <- initAlpha
  pi <- rep(NA, M)
  pi[ yi <= 1e-4 ] <- 0.9
  pi[ yi > 1e-4 ] <- 0.1
  names(pi) <- names(yi)
  
  # storing
  alphalist <- rep( NA, iter_max )
  liclist <- rep( NA, iter_max )
  lclist <- rep(NA, iter_max)
  betalist <- list()

  # initial stores
  alphalist[1] <- alpha
  liclist[1] <- -10000000
  lclist[1] <- -10000000
  betalist[1] <- NA

  # EM start ####

  options(digits = 20)

  while( iter_cur < iter_max) {


    # E Step
    zi <- (alpha* (yi^(alpha-1)) * pi ) / ( (alpha* (yi^(alpha-1)) * pi) + (1-pi) ) # posterior probability of having signal or being non-null for each SNP
    zi[ zi < 0.01 ] <- 0.01
    zi[ zi > 0.99 ] <- 0.99

    # M Step using regression tree from CART
    dat <- as.data.frame(cbind(annMat, 'zi'=zi)) # data to use in regression
    reg_fit <- lm(zi ~ . , data=dat)

    # update pi
    pi <- c(predict(reg_fit, as.data.frame(annMat)))
    pi[ pi < 0.1 ] <- 0.1 # change this to 0.1/0.25
    pi[ pi > 0.9 ] <- 0.9 # change this to 0.9/0.75

    # update alpha
    alpha <- max(min(-(sum(zi)/sum(zi*log(yi))), 0.999), 0.001) # setting 0.001 < alpha < 0.999

    # Incomplete log-Likelihood
    lic <- sum(log ( (pi * alpha* (yi^(alpha-1)) ) + (1-pi) ))

    # Complete data log-likelihood
    lc <- sum( (zi * ( log(pi) + log(alpha) + ((alpha - 1) * log(yi)) )) + ((1 - zi) * log(1-pi)) )

    # store values
    alphalist[iter_cur] <- alpha
    liclist[iter_cur] <- lic
    lclist[iter_cur] <- lc

    # stopping rules
    if ( abs(liclist[iter_cur] - liclist[iter_cur-1]) <= 1e-4 &
         abs(alphalist[iter_cur] - alphalist[iter_cur-1]) <= 1e-4) {
      message('Info: Incomplete log-likelihood and alpha converged in Stage 1 of GPA-Tree.')
      break
    }

    iter_cur <- iter_cur + 1

  }

  # store betas from linear regression fit
  beta_list <- rep(NA, (nAnn+1))
  beta_list[1] <- summary(reg_fit)$coefficient[1,1]
  for (i in 1:(nAnn)) {
    beta_list[i+1] <- summary(reg_fit)$coefficient[(i+1),1]
  }
  names(beta_list)[1:(nAnn+1)] <- c('beta_int',paste('beta_A', 1:nAnn, sep = ''))

  message( "Info: Stage 1 of GPA-Tree completed.")
  
  return_out <- list()
  return_out$num_iter_stage1 = iter_cur
  return_out$pi_stage1 = pi
  return_out$alpha_stage1 = alpha
  return_out$beta_stage1 = beta_list
  
  return(return_out)

}

