#' Calculate calibration measures (ICI, E50, E90, Emax)
#'
#' @param data Validation data
#' @param thorizon time horizon for prediction
#' @param fit Model fit
#' @param Tstop observed time
#' @param status variable status
#' @param cause number indicating the event of interest. See ?riskRegression::predictRisk()
#'
#' @return
#' 
#' @author Daniele Giardiello
#' 
#' @examples
#' 
cal_measures<-function(data,thorizon,fit,Tstop,status,cause) {
  valid.df <- data.frame(data)
  pred <- riskRegression::predictRisk(fit,cause=cause,newdata=data,times=thorizon)
  cll.pred <- log(-log(1-pred))
  
  # Create grid along which to create calibration curves.
  range.pred <- quantile(pred,probs=c(0.01,0.99))
  pred.grid<- seq(from=range.pred[1],to=range.pred[2],length=100)
  cll.pred.grid <- log(-log(1-pred.grid))
  
  remove(range.pred)
  
  # Create RCS components using 3 knots.
  rcs3.mat <- Hmisc::rcspline.eval(cll.pred,nk=3,inclx=T)
  knots3 <- attr(rcs3.mat,"knots")
  valid.df$rcs3.x1 <- rcs3.mat[,1]
  valid.df$rcs3.x2 <- rcs3.mat[,2]
  
  # Create RCS terms for grid along which predictions will be obtained.
  # Use the CLL transformation of the grid. Use the knots defined above.
  grid3.rcs <- Hmisc::rcspline.eval(cll.pred.grid,knots=knots3,inclx=T)
  rcs3.x1 <- grid3.rcs[,1]
  rcs3.x2 <- grid3.rcs[,2]
  rcs3.df <- data.frame(cbind(rcs3.x1,rcs3.x2))
  
  remove(rcs3.mat,knots3,grid3.rcs,rcs3.x1,rcs3.x2)
  
  # Fit F-G model to predicted probabilities using RCS.
  # Obtain estimated probabilities from fitted model at the grid points.
  # This is for the calibration curve.
  
  # [NOTE] Expand validation data using mstate::crprep()
  valid.df.w<-mstate::crprep(Tstop=Tstop,status=status,trans=cause,
                     keep=c('rcs3.x1','rcs3.x2'),data=valid.df)
  
  model.calibrate.fg <- rms::cph(
    survival::Surv(Tstart,Tstop,status==1)~rcs3.x1+rcs3.x2,
    weights=weight.cens,
    x=T,y=T,
    surv=T,
    data=valid.df.w
  )
  
  #obs.grid.fg <- 1-survest(model.calibrate.fg,newdata=rcs3.df,time=5)$surv
  ## [NOTE] save the confidence levels of the grid useful for the calibration plot confidence bands
  #obs.grid.fg.upper<-1-survest(model.calibrate.fg,newdata=rcs3.df,time=5)$lower
  #obs.grid.fg.lower<-1-survest(model.calibrate.fg,newdata=rcs3.df,time=yr*365)$upper
  
  obs.fg <- 1 - rms::survest(model.calibrate.fg,newdata=valid.df,time=5)$surv
  ici.fg <- mean(abs(pred - obs.fg))
  E50.fg <- median(abs(pred - obs.fg))
  E90.fg <- quantile(abs(pred - obs.fg),probs=0.90)
  Emax.fg <-max(abs(pred - obs.fg))
  
  res_calmeas<-c(ici.fg,E50.fg,E90.fg,Emax.fg)
  names(res_calmeas)<-c('ICI','E50','E90','Emax')
  return(res_calmeas)
}