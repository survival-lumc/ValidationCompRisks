# 8th March 2021
# Goal: function to calculate observed and expected ratio of competing
# risk models 
# Author: Daniele Giardiello


OE_function<-function(fit,newdata,cause,thorizon,obs_cif,obs_events) {
  F_t<-predictRisk(fit,newdata=newdata,cause=cause,times=thorizon)
  E_t<-mean(F_t)
  OE<-obs_cif/E_t
  OE_lower<-OE*exp(-qnorm(0.975)*sqrt(1/obs_events))
  OE_upper<-OE*exp(qnorm(0.975)*sqrt(1/obs_events))
  res<-c(OE,OE_lower,OE_upper)
  res<-round(res,2)
  names(res)<-c('Estimate','Lower.95','Upper.95')
  return(res)
}