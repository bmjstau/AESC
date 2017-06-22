computeG_means_paired <- function(mean1, mean2, n_pairs, r_pairs, sd_diff, print_d = FALSE){
## Computes Cohen's d and Hedge's g for paired samples / pre-post-designs.
##
## Formulas taken from Introduction to Meta-Analysis by M. Borenstein,
## L. V. Hedges, J. P. T. Higgins, H. R. Rothstein (2009)
  
  #sanitize inputs
  stopifnot(is.numeric(mean1), is.numeric(mean2),
            is.numeric(r_pairs), is.numeric(sd_diff),
            is.numeric(n_pairs), is.logical(print_d))

  #compute cohen's d plus its SD and SE, plus 95% CIs
  sd_pooled <- sd_diff/sqrt(2*(1-r_pairs))
  d <- (mean1-mean2)/sd_pooled
  sd_d <- (1/n_pairs+d^2/(2*n_pairs))*2*(1-r_pairs)
  se_d <- sqrt(sd_d)
  d_ci <- d + c(-qnorm(.975)*se_d, qnorm(.975)*se_d)
  
  #compute Hedge's g plus its SD and SE, plus 95% CIs
  J <- 1-(3/(4*(n_pairs-1)-1))
  g <- J*d
  sd_g <- J^2*sd_d
  se_g <- sqrt(sd_g)
  g_ci <- g + c(-qnorm(.975)*se_g, qnorm(.975)*se_g)
  

  #print the results
  if(print_d){
    cat("Effect Sizes (paired):", "\n",
        "Cohen's d: ", d, "[", d_ci, "]", "\n",
        "sd =", sd_d, ", se =", se_d, "\n\n",
        "Hedge's g: ", g, "[", g_ci, "]", "\n",
        "sd =", sd_g, ", se =", se_g, "\n\n"
    )
  } else {
    cat("Effect Sizes (paired):", "\n",
        "Hedge's g: ", g, "[", g_ci, "]", "\n",
        "sd =", sd_g, ", se =", se_g, "\n\n"
    )
  }
  

}