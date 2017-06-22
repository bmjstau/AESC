computeG_means_independent <- function(mean1, mean2, sd1, sd2, n1, n2){
## Computes Cohen's d and Hedge's g for independent samples
##
## Formulas taken from Introduction to Meta-Analysis by M. Borenstein,
## L. V. Hedges, J. P. T. Higgins, H. R. Rothstein (2009)
  
  #sanitize inputs
  stopifnot(is.numeric(mean1), is.numeric(mean2),
            is.numeric(sd1), is.numeric(sd2),
            is.numeric(n1), is.numeric(n2))

  #compute cohen's d plus its SD and SE, plus 95% CIs
  sd_pooled <- sqrt( ((n1-1)*sd1^2+(n2-1)*sd2^2) / (n1+n2-2) )
  d <- (mean1-mean2)/sd_pooled
  sd_d <- (n1+n2)/(n1*n2) + d^2/(2*(n1+n2))
  se_d <- sqrt(sd_d)
  d_ci <- d + c(-qnorm(.975)*se_d, qnorm(.975)*se_d)
  
  #compute Hedge's g plus its SD and SE, plus 95% CIs
  df_sd <- n1+n2-2
  J <- 1-(3/(4*df_sd-1))
  g <- J*d
  sd_g <- J^2*sd_d
  se_g <- sqrt(sd_g)
  g_ci <- g + c(-qnorm(.975)*se_g, qnorm(.975)*se_g)
  

  #print the results
  cat("Effect Sizes (independent):", "\n",
      "Cohen's d: ", d, "[", d_ci, "]", "\n",
      "sd =", sd_d, ", se =", se_d, "\n\n",
      "Hedge's g: ", g, "[", g_ci, "]", "\n",
      "sd =", sd_g, ", se =", se_g, "\n\n"
      )

}