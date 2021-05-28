# To-do:
# - Edit source data directly so no re-leveling needed


# Load libraries and data -------------------------------------------------


library(survival)
library(pec)
library(riskRegression)
library(splines)

rdata <- readRDS("Data/rdata.rds")
vdata <- readRDS("Data/vdata.rds")

rdata$hr_status <- relevel(rdata$hr_status, ref = "ER and/or PR +")
rdata$status_num <- as.numeric(rdata$status) - 1
vdata$hr_status <- relevel(vdata$hr_status, ref = "ER and/or PR +")
vdata$status_num <- as.numeric(vdata$status) - 1


# Fit cause-specific models -----------------------------------------------


fit_csh <- CSC(
  formula = Hist(time, status_num) ~ age + size + ncat + hr_status, 
  data = rdata
)

# External validation at 5 years
horizon <- 4.99

score_vdata <- Score(
  list("csh_validation" = fit_csh),
  formula = Hist(time, status_num) ~ 1,
  cens.model = "km", 
  data = vdata, 
  conf.int = TRUE, 
  times = horizon,
  metrics = c("auc", "brier"),
  summary = c("ipa"), 
  cause = 1,
  plots = "calibration"
)

# Calculate predicted probabilities 
pred <- predictRisk(
  fit_csh, 
  cause = 1, 
  newdata = vdata, 
  times = horizon
)


# Calibration (O/E) -------------------------------------------------------


# First calculate Aalen-Johansen estimates (as 'observed')
obj <- summary(survfit(Surv(time, status) ~ 1, data = vdata), times = horizon)
aj <- list("obs" = obj$pstate[, 2], "se" = obj$std.err[, 2])

OE <- aj$obs / mean(pred)
OE_summary <- c(
  "OE" = OE,
  "lower" = exp(log(OE - qnorm(0.975) * aj$se / aj$obs)),
  "upper" = exp(log(OE + qnorm(0.975) * aj$se / aj$obs))
)
OE_summary


# Calibration plot (pseudo-obs approach) ----------------------------------


calplot_pseudo <- plotCalibration(
  score_vdata,
  brier.in.legend = FALSE,
  auc.in.legend = FALSE, 
  cens.method = "pseudo",
  bandwidth = 0.05, # or NULL to leave to prodlim::neighborhood(pred)$bandwidth; look at paper for justification
  cex = 1, 
  round = FALSE, # Important, keeps all unique risk estimates (here = 665)
  xlim = c(0, 0.6), 
  ylim = c(0, 0.6), 
  rug = TRUE
)

# We can extract predicted and observed, this will depend on smoothing
dat_pseudo <- calplot_pseudo$plotFrames$csh_validation
absdiff_pseudo <- abs(dat_pseudo$Pred - dat_pseudo$Obs)

# Collect all numerical summary measures
numsum_pseudo <- c(
  "ICI" = mean(absdiff_pseudo),
  setNames(quantile(absdiff_pseudo, c(0.5, 0.9)), c("E50", "E90")),
  "Emax" = max(absdiff_pseudo),
  "squared_bias" = mean((dat_pseudo$Pred - dat_pseudo$Obs)^2)
)
numsum_pseudo


# Calibration plot (subdistribution approach) -----------------------------


# Add predicted risk and complementary log-log of it to dataset
vdata$pred <- pred
vdata$cll_pred <- log(-log(1 - pred))

# 5 knots seems to give equivalent graph to pseudo method with bw = 0.05
rcs_vdata <- ns(vdata$cll_pred, df = 6)
class(rcs_vdata) <- "matrix"
colnames(rcs_vdata) <- paste0("basisf_", colnames(rcs_vdata))
vdata_bis <- cbind.data.frame(vdata, rcs_vdata)

# With FGR
form_fgr <- reformulate(
  termlabels = colnames(rcs_vdata),
  response = "Hist(time, status_num)"
)

# Regress subdistribution of event of interest on cloglog of predicted risks
calib_fgr <- FGR(
  form_fgr,
  cause = 1,
  data = vdata_bis
)

# Imo confidence intervals can only be obtained by resampling (discuss later)
dat_fgr <- cbind.data.frame(
  "obs" = predict(calib_fgr, times = 5, newdata = vdata_bis),
  "pred" = vdata$pred
)

dat_fgr <- dat_fgr[order(dat_fgr$pred), ]
plot(
  x = dat_fgr$pred, 
  y = dat_fgr$obs, 
  type = "l",
  xlim = c(0, 0.6), 
  ylim = c(0, 0.6),
  xlab = "Predicted risk",
  ylab = "Estimated actual risk"
)
abline(a = 0, b = 1, lty = "dashed", col = "red")

# Numerical measures
absdiff_fgr <- abs(dat_fgr$pred - dat_fgr$obs)

numsum_fgr <- c(
  "ICI" = mean(absdiff_fgr),
  setNames(quantile(absdiff_fgr, c(0.5, 0.9)), c("E50", "E90")),
  "Emax" = max(absdiff_fgr),
  "squared_bias" = mean((dat_fgr$pred - dat_fgr$obs)^2)
)
numsum_fgr

# (Plot this and pseudo together)
plot(
  x = dat_fgr$pred, 
  y = dat_fgr$obs, 
  type = "l", 
  xlim = c(0, 0.6), 
  ylim = c(0, 0.6),
  col = "blue",
  lwd = 2,
  xlab = "Predicted risk",
  ylab = "Estimated actual risk"
)
lines(x = dat_pseudo$Pred, y = dat_pseudo$Obs, col = "lightblue", lwd = 2)
abline(a = 0, b = 1, lty = "dashed", col = "red")
legend(
  x = 0, 
  y = 0.6, 
  legend = c("Subdistribution", "Pseudo-observations"),
  col = c("blue", "lightblue"),
  lty = rep(1, 2),
  lwd = rep(2, 2),
  bty = "n"
)


# Discrimination ----------------------------------------------------------


# AUC described in paper (same as AUC_2 from timeROC)
score_vdata$AUC$score

# C-index
pec::cindex(fit_csh, cause = 1, eval.times = horizon, data = vdata)$AppCindex
# Optional bootstrap for cindex here


# Prediction error --------------------------------------------------------


score_vdata$Brier$score # brier + IPA
# Add optional bootstrap for IPA here


# Clinical utility --------------------------------------------------------


# (Source Daniele's function here)
