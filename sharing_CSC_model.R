# Load libraries and data -------------------------------------------------


# General packages
library(survival)
library(riskRegression)

# Load datasets
rdata <- readRDS("Data/rdata.rds")
vdata <- readRDS("Data/vdata.rds")


# Fitting and sharing necessary parts for validation ----------------------


# Fit model
fit_csh <- riskRegression::CSC(
  formula = Hist(time, status_num) ~ age + size + ncat + hr_status, 
  data = rdata
)

# The bare minimum to share: coefficients, baseline hazards and model formulas
# (No data sharing is needed)
model_info <- list(
  "coefficients" = stats::coef(fit_csh),
  "baseline_hazards" = lapply(fit_csh$models, function(mod) survival::basehaz(mod, centered = FALSE)),
  "model_formulas" = lapply(fit_csh$models, function(mod) mod[["formula"]])
)

# Can be saved to then be shared as:
#saveRDS(model_info, file = "model_info.rds")


# Calculating predicted risks  --------------------------------------------


# Function to calculate predicted risks
predictRisk_shared_CSC <- function(model_info,
                                   newdata,
                                   horizon,
                                   primary_event) {
  
  # -- Basic checks
  
  # Check model_info components
  info_names <- c("coefficients", "baseline_hazards", "model_formulas")
  if (any(!(names(model_info) %in% info_names))) {
    stop(paste0("Names of model_info components should be: ", paste(info_names, collapse = ", ")))
  }
  n_causes <- unique(vapply(model_info, length, FUN.VALUE = integer(1L)))
  if (length(n_causes) > 1) {
    stop("The elements of model_info should all be of length equal to number of competing risks!")
  } 
  
  # Check predictor variables in newdata
  predictor_vars <- Reduce(
    f = "intersect",
    x = lapply(model_info$model_formulas, function(form) all.vars(stats::update(form, 1 ~ .)))
  )
  if (any(!(predictor_vars %in% colnames(newdata)))) {
    which_missing <- predictor_vars[!(predictor_vars %in% colnames(newdata))]
    stop(paste0("newdata does not contain the following predictors: ", paste(which_missing, collapse = ", ")))
  }

  # -- Absolute risk prediction
  
  causes_ind <- seq_len(n_causes)
  
  # Calculate linear predictors for all causes in new dataset
  linpreds <- lapply(causes_ind, function(cause) {
    mod_matrix <- stats::model.matrix(model_info$model_formulas[[cause]], data = newdata)
    unname(drop(mod_matrix[, -1] %*% model_info$coefficients[[cause]]))
  })
  
  # Compute absolute risks for each individual
  preds <- vapply(seq_len(nrow(vdata)), function(id) {
    
    # Calculate individual-specific cause-specific hazards
    time_points <- model_info$baseline_hazards[[primary_event]][["time"]]
    hazards <- vapply(causes_ind, function(cause) {
      cumhaz <- model_info$baseline_hazards[[cause]][["hazard"]] * exp(linpreds[[cause]][[id]])
      diff(c(0, cumhaz))
    }, FUN.VALUE = numeric(length(time_points)))
    
    # Calculate event-free survival
    surv <- cumprod(1 - rowSums(hazards))
    
    # Calculate cumulative incidence
    cuminc <- cumsum(hazards[, primary_event] * c(1, surv[-length(surv)]))
    cuminc_horizon <- cuminc[length(time_points[time_points <= horizon])]
    return(cuminc_horizon)
  }, FUN.VALUE = numeric(1L))
  
  return(preds)
}

# Set prediction horizon and event of interest
horiz <- 5
cause_interest <- 1

# Get absolute risk
preds <- predictRisk_shared_CSC(model_info, newdata = vdata, horizon = horiz, primary_event = cause_interest)

# Check results against riskRegression::predictRisk()
preds_riskReg <- drop(predictRisk(fit_csh, vdata, times = horiz, cause = cause_interest))
all.equal(preds, preds_riskReg)


# Validation with riskRegression ------------------------------------------


# Supply predicted probabilities, and proceed as with the minimal script..
score_vdata <- Score(
  list("csh_validation" = preds),
  formula = Hist(time, status_num) ~ 1, 
  cens.model = "km", 
  data = vdata, 
  conf.int = TRUE, 
  times = horiz,
  metrics = c("auc", "brier"),
  summary = c("ipa"), 
  cause = cause_interest,
  plots = "calibration"
)

#.. We need to do all of this since predictRisk() needs coxph models to be run
# with x = TRUE, meaning data has to be passed on 