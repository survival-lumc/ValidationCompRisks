# Load libraries and data -------------------------------------------------


# General packages
library(survival)
library(riskRegression)

# Load datasets
rdata <- readRDS("Data/rdata.rds")
vdata <- readRDS("Data/vdata.rds")


# Fitting and sharing necessary parts for validation ----------------------


# Fit model
fit_csh <- CSC(
  formula = Hist(time, status_num) ~ age + size + ncat + hr_status, 
  data = rdata
)

# The bare minimum to share: coefficients, baseline hazards and model formulas
# (No data sharing is needed)
model_info <- list(
  "coefficients" = coef(fit_csh),
  "baseline_hazards" = lapply(fit_csh$models, function(mod) basehaz(mod, centered = FALSE)),
  "model_formulas" = lapply(fit_csh$models, function(mod) mod[["formula"]])
)

# Can be saved to then be shared as:
#saveRDS(model_info, file = "model_info.rds")


# Calculating predicted risks  --------------------------------------------


# Pick primary event and prediction horizon
primary_event <- 1
horizon <- 5

# Get index of causes
causes <- seq_along(model_info$coefficients)

# Calculate linear predictors for all causes in new dataset (vdata)
linpreds <- lapply(causes, function(cause) {
  mod_matrix <- model.matrix(model_info$model_formulas[[cause]], data = vdata)
  unname(drop(mod_matrix[, -1] %*% model_info$coefficients[[cause]]))
})

# Compute absolute risks for each individual
preds <- vapply(seq_len(nrow(vdata)), function(id) {
  
  # Calculate individual-specific cause-specific hazards
  hazards <- lapply(causes, function(cause) {
    cumhaz <- model_info$baseline_hazards[[cause]][["hazard"]] * exp(linpreds[[cause]][[id]])
    diff(c(0, cumhaz))
  })
  
  # Calculate event-free survival
  surv <- cumprod(1 - Reduce("+", hazards))
                  
  # Calculate cumulative incidence
  cuminc <- cumsum(hazards[[primary_event]] * c(1, surv[-length(surv)]))
  time_points <- model_info$baseline_hazards[[primary_event]][["time"]]
  cuminc_horizon <- cuminc[length(time_points[time_points <= horizon])]
  return(cuminc_horizon)
}, FUN.VALUE = numeric(1L))

# Check results against riskRegression::predictRisk()
preds_riskReg <- drop(predictRisk(fit_csh, vdata, times = horizon, cause = primary_event))
all.equal(preds, preds_riskReg)


# Validation with riskRegression ------------------------------------------


# Supply predicted probabilities, and proceed as with the minimal script..
score_vdata <- Score(
  list("csh_validation" = preds),
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

#.. We need to do all of this since predictRisk() needs coxph models to be run
# with x = TRUE, meaning data has to be passed on 