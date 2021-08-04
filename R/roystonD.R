#' Calculate Royston D and R2D
#'
#' @param pred predicted absolute risk at time t
#' @param data data 
#' @param time follow-up time
#' @param status status (including competing risks)
#' @param primary_event event of interest
#'
#' @return
#'
#' @author Daniele Giardiello
#'
#' @examples
#' 

library(cmprsk)
royston_R2D_cmprsk <- function(pred,
                                data, 
                                time,
                                status,
                                primary_event) {
  temp <- data
  temp$rank <- rank(pred,
                    na.last = TRUE,
                    ties.method = "first") 
  
  Kappa <- (8 / pi) ^0.5   #kappa scaling coefficient 
  rankit <- (temp$rank - 0.375) / (nrow(temp) + 0.25)
  normal_rankit <- qnorm(rankit)
  scaled_rankit <- normal_rankit / Kappa
  fg <- crr(ftime = time,
            fstatus = status,
            cov1 = scaled_rankit,
            failcode = primary_event,
            cencode = 0)
  
  alpha <- 0.05 # for confidence interval
  res <- c(
    "Royston D" = unname(fg$coef),
    "Lower .95" = fg$coef - qnorm(1 - alpha /2) * sqrt(fg$var),
    "Upper .95" = fg$coef + qnorm(1 - alpha /2) * sqrt(fg$var),
    "R2D" =  unname((fg$coef^2 / Kappa^2)/((pi^2 / 6)+(fg$coef^2 / Kappa^2)))  
  )
  
  return(res)
}