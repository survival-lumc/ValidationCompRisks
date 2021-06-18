# Load libraries and data -------------------------------------------------

library(survival)
library(pec)
library(riskRegression)
library(splines)
library(sandwich)
library(lmtest)

rdata <- readRDS("Data/rdata.rds")
vdata <- readRDS("Data/vdata.rds")


# Fit cause-specific models -----------------------------------------------


fit_csh <- CSC(
  formula = Hist(time, status_num) ~ age + size + ncat + hr_status, 
  data = rdata
)

# External validation at 5 years
horizon <- 5

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


# Calibration plot (pseudo-obs approach) ----------------------------------


calplot_pseudo <- plotCalibration(
  score_vdata,
  brier.in.legend = FALSE,
  auc.in.legend = FALSE, 
  cens.method = "pseudo",
  bandwidth = 0.05, # or NULL to leave default prodlim::neighborhood(pred)$bandwidth; look at paper for justification
  cex = 1, 
  round = FALSE, # Important, keeps all unique risk estimates rather then rounding (here = 665)
  xlim = c(0, 0.6), 
  ylim = c(0, 0.6), 
  rug = TRUE
)

# We can extract predicted and observed, this will depend on smoothing
dat_pseudo <- calplot_pseudo$plotFrames$csh_validation

# Make sure to use all predicted risks (not just unique ones)
diff_pseudo <- pred - dat_pseudo$Obs[match(pred, dat_pseudo$Pred)]

# Collect all numerical summary measures
numsum_pseudo <- c(
  "ICI" = mean(abs(diff_pseudo)),
  setNames(quantile(abs(diff_pseudo), c(0.5, 0.9)), c("E50", "E90")),
  "Emax" = max(abs(diff_pseudo)),
  "squared_bias" = mean(diff_pseudo^2)
)
numsum_pseudo


# Calibration plot (subdistribution approach) -----------------------------


# Add predicted risk and complementary log-log of it to dataset
vdata$pred <- pred
vdata$cll_pred <- log(-log(1 - pred))

# 5 knots seems to give equivalent graph to pseudo method with bw = 0.05
n_internal_knots <- 5 # Pick between 3 (more smoothing, less flexible) -5 (less smoothing, more flexible), 
rcs_vdata <- ns(vdata$cll_pred, df = n_internal_knots + 1)
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
diff_fgr <- dat_fgr$pred - dat_fgr$obs

numsum_fgr <- c(
  "ICI" = mean(abs(diff_fgr)),
  setNames(quantile(abs(diff_fgr), c(0.5, 0.9)), c("E50", "E90")),
  "Emax" = max(abs(diff_fgr)),
  "squared_bias" = mean(diff_fgr^2)
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


# Calibration intercept/slope ---------------------------------------------


# Use pseudo values calculated by score (can also use pseudo::pseudoci)
pseudos <- data.frame(score_vdata$Calibration$plotframe)

# Notes:
# - this is dataframe with ACTUAL pseudovalues, not the smoothed ones
# - ID is not the id in data, it is just a number assigned to each row of 
# the original validation data sorted by time and event indicator
pseudos
pseudos$pseudovalue # the pseudo values
pseudos$cll_pred <- log(-log(1 - pseudos$risk)) # add the cloglog risk ests

# Fit model for calibration intercept
fit_cal_int <- glm(
  pseudovalue ~ offset(cll_pred), 
  data = pseudos,
  family = quasi(link = "cloglog"), 
  start = 0
)

# Fit model for calibration slope
# or: update(fit_cal_int, formula. = . ~ . + cll_pred, start = c(0, 0))
fit_cal_slope <- glm(
  pseudovalue ~ offset(cll_pred) + cll_pred, 
  data = pseudos,
  family = quasi(link = "cloglog"), 
  start = c(0, 0)
)

# Perform joint test 
betas <- coef(fit_cal_slope)
vcov_mat <- vcovCL(fit_cal_slope, cluster = ~ ID)
wald <- drop(betas %*% solve(vcov_mat) %*% betas)
pchisq(wald, df = 2, lower.tail = FALSE)

# Slope test + value and confidence interval
test_slope <- coeftest(fit_cal_slope, vcov. = vcovCL, cluster = ~ ID) 
c("slope" = 1 + betas[["cll_pred"]], 1 + confint(test_slope)["cll_pred", ])

# Test for intercept
test_int <- coeftest(fit_cal_int, vcov. = vcovCL, cluster = ~ ID)
c("intercept" = coef(test_int), drop(confint(test_int)))


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
