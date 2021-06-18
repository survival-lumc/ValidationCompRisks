# Load libraries and data -------------------------------------------------

library(survival)
library(pec)
library(riskRegression)
library(splines)
library(sandwich)
library(lmtest)

rdata <- readRDS("Data/rdata.rds")
vdata <- readRDS("Data/vdata.rds")

# Set seed (for bootstrapping)
set.seed(2021)


# Fit cause-specific models -----------------------------------------------


fit_csh <- CSC(
  formula = Hist(time, status_num) ~ age + size + ncat + hr_status, 
  data = rdata
)

# External validation at 5 years
horizon <- 5
primary_event <- 1 # set 2 if cause 2 is of interest 

score_vdata <- Score(
  list("csh_validation" = fit_csh),
  formula = Hist(time, status_num) ~ 1,
  cens.model = "km", 
  data = vdata, 
  conf.int = TRUE, 
  times = horizon,
  metrics = c("auc", "brier"),
  summary = c("ipa"), 
  cause = primary_event,
  plots = "calibration"
)

# Calculate predicted probabilities 
pred <- predictRisk(
  fit_csh, 
  cause = primary_event, 
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
  cause = primary_event,
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
aj <- list("obs" = obj$pstate[, primary_event + 1], "se" = obj$std.err[, primary_event + 1])

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
head(pseudos)
head(pseudos$pseudovalue) # the pseudo values
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
test_slope
c("slope" = 1 + betas[["cll_pred"]], 1 + confint(test_slope)["cll_pred", ])

# Test for intercept
test_int <- coeftest(fit_cal_int, vcov. = vcovCL, cluster = ~ ID)
test_int
c("intercept" = coef(test_int), drop(confint(test_int)))


# Discrimination ----------------------------------------------------------


# AUC described in paper (same as AUC_2 from timeROC)
score_vdata$AUC$score

# C-index
cindex_csh <- pec::cindex(
  object = fit_csh, 
  formula = Hist(time, status_num) ~ 1,
  cause = primary_event, 
  eval.times = horizon, 
  data = vdata
)$AppCindex$CauseSpecificCox
# Optional bootstrap for cindex here (or at end?)


# Prediction error --------------------------------------------------------


score_vdata$Brier$score # brier + IPA
# Add optional bootstrap for IPA here


# Clinical utility --------------------------------------------------------


# Minimal version (better to use Daniele's function):

# 1. Set grid of thresholds
thresholds <- seq(0, 0.6, by = 0.01)

# 2. Calculate Aelen johansen for all patients exceeding threshold (i.e. treat-all)
survfit_all <- summary(
  survfit(Surv(time, status) ~ 1, data = vdata), 
  times = horizon
)
f_all <- survfit_all$pstate[primary_event + 1]

# 3. Calculate Net Benefit across all thresholds
list_nb <- lapply(thresholds, function(ps) {
  
  # Treat all
  NB_all <- f_all - (1 - f_all) * (ps / (1 - ps))
  
  # Based on threshold
  p_exceed <- mean(vdata$pred > ps)
  survfit_among_exceed <- try(
    summary(
      survfit(Surv(time, status) ~ 1, data = vdata[vdata$pred > ps, ]), 
      times = horizon
    ), silent = TRUE
  )
  
  # If a) no more observations above threshold, or b) among subset exceeding..
  # ..no indiv has event time >= horizon, then NB = 0
  if (class(survfit_among_exceed) == "try-error") {
    NB <- 0
  } else {
    f_given_exceed <- survfit_among_exceed$pstate[primary_event + 1]
    TP <- f_given_exceed * p_exceed
    FP <- (1 - f_given_exceed) * p_exceed
    NB <- TP - FP * (ps / (1 - ps))
  }
  
  # Return together
  df_res <- data.frame("threshold" = ps, "NB" = NB, "treat_all" = NB_all)
  return(df_res)
})

# Combine into data frame
df_nb <- do.call(rbind.data.frame, list_nb)
head(df_nb)

# Make basic decision curve plot
par(
  xaxs = "i", 
  yaxs = "i", 
  las = 1, 
  mar = c(6.1, 5.8, 4.1, 2.1), 
  mgp = c(4.25, 1, 0)
)
plot(
  df_nb$threshold, 
  df_nb$NB,
  type = "l", 
  lwd = 2,
  ylim = c(-0.1, 0.1),
  xlim = c(0, 0.5), 
  xlab = "",
  ylab = "Net Benefit",
  bty = "n", 
  xaxt = "n"
)
lines(df_nb$threshold, df_nb$treat_all, type = "l", col = "darkgray", lwd = 2)
abline(h = 0, lty = 2, lwd = 2)
legend(
  "topright", 
  c("Treat all", "Treat none", "Prediction model"),
  lwd = c(2, 2, 2), 
  lty = c(1, 2, 1), 
  col = c("darkgray", "black", "black"), 
  bty = "n"
)
axis(side = 1, at = c(0, 0.1, 0.2, 0.3, 0.4, 0.5))
axis(
  side = 1,
  pos = -0.145, 
  at = c(0.1, 0.2, 0.3, 0.4, 0.5),
  labels = c("1:9", "1:4", "3:7", "2:3", "1:1")
)
mtext("Threshold probability", 1, line = 2)
mtext("Harm to benefit ratio", 1, line = 5)
title("Validation data")

# Restore old graphical parameters
dev.off() 

# -- End minimal script


# Extra: bootstrap confidence intervals -----------------------------------


# Repeat whole modelling process:
# - Resample training data, evaluate each model on test set
# (i.e. not resample test set with same model)
B <- 250

cindexes <- vapply(seq_len(B), function(b) {
  
  # Fit model on bootstrapped development data
  obj <- fit_csh_boot <- CSC(
    formula = Hist(time, status_num) ~ age + size + ncat + hr_status, 
    data = rdata[sample(nrow(rdata), replace = TRUE), ]
  )
  
  # Get cindex on validation data
  pec::cindex(
    object = obj, 
    formula = Hist(time, status_num) ~ 1,
    cause = 1, 
    eval.times = horizon, 
    data = vdata,
    verbose = FALSE
  )$AppCindex$CauseSpecificCox
}, FUN.VALUE = numeric(1))

hist(cindexes)
cindex_csh
quantile(cindexes, probs = c(0.025, 0.975))


# Do we do it for the calibration curves..? Probs not
