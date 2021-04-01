# 8th March 2021
# Goal: function to calculate observed and expected (OE) ratio of competing
# risk models 
# Author: Daniele Giardiello

# Input:
# fit: FG model or CSC model a riskRegression::FGR() or riskRegression:CSC() object
# newdata: validation data
# cause: code for event of interest. if the event of interest is 1 then cause=1. More details are in the help of riskRegression::predictRisk();
# thorizon: the fixed time horizon for prediction to calculate the OE ratio at the fixed time horizon t;
# obs_cif: the cumulative incidence at time horizon;
# std.error: standard error at time horizon to calculate the confidence intervals 
#            using Summary(survfit()...)$std.err. 


OE_function<-function(fit,newdata,cause,thorizon,obs_cif,
                      std.error) {
  F_t<-predictRisk(fit,newdata=newdata,cause=cause,times=thorizon)
  E_t<-mean(F_t)
  OE<-obs_cif/E_t
  OE_lower<-exp(log(OE-qnorm(0.975)*((std.error)/obs_cif)))
  OE_upper<-exp(log(OE+qnorm(0.975)*((std.error)/obs_cif)))
  res<-c(OE,OE_lower,OE_upper)
  res<-round(res,2)
  names(res)<-c('Estimate','Lower.95','Upper.95')
  return(res)
}
