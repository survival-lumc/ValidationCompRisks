# Load libraries and data -------------------------------------------------


# General packages
pkgs <- c("survival", "mstate", "rms")
vapply(pkgs, function(pkg) {
  if (!require(pkg, character.only = TRUE)) install.packages(pkg)
  require(pkg, character.only = TRUE, quietly = TRUE)
}, FUN.VALUE = logical(length = 1L))

# Install latest development version of riskRegression
if (!require("devtools", character.only = TRUE)) install.packages("devtools")
if (!require("riskRegression", character.only = TRUE)) devtools::install_github("tagteam/riskRegression")
require("riskRegression", character.only = TRUE)

# Load datasets
rdata <- readRDS("Data/rdata.rds")


# Expand data to prepare for fitting the model using mstate::crprep()
rdata.w <- crprep(
  Tstop = "time",
  status = "status_num",
  trans = c(1, 2),
  id = "id",
  keep = c("age", "size", "ncat", "hr_status"),
  data = rdata
)

# Save extended data with weights for recurrence (failcode=1)
primary_event <- 1 # breast cancer recurrence is the primary event
rdata.w1 <- rdata.w[rdata.w$failcode == primary_event, ]


# Cumulative incidence function ------------
mfit_rdata <- survfit(
  Surv(Tstart, Tstop, status == 1) ~ 1,
  data = rdata.w1, weights = weight.cens
)

par(xaxs = "i", yaxs = "i", las = 1)
plot(mfit_rdata,
     col = 1, lwd = 2,
     xlab = "Years since BC diagnosis",
     ylab = "Cumulative incidence", bty = "n",
     ylim = c(0, 0.25), xlim = c(0, 5), fun = "event", conf.int = TRUE
)
title("Development data")

# Cumulative incidences -------------------
smfit_rdata <- summary(mfit_rdata, times = c(1, 2, 3, 4, 5))

CIF <- cbind.data.frame(
  "time" = smfit_rdata$time,
  "CI" = 1 - smfit_rdata$surv,
  "2.5 %" = 1 - smfit_rdata$upper,
  "97.5 %" = 1 - smfit_rdata$lower
)
CIF


### Check non-linearity of continuous predictors ----------------

# Defining knots of the restricted cubic splines ------------------
# Extract knots position of the restricted cubic spline based on the
# original distribution of the data

# Age at breast cancer diagnosis
rcs3_age <- rcspline.eval(rdata$age, nk = 3)
attr(rcs3_age, "dim") <- NULL
attr(rcs3_age, "knots") <- NULL
pos_knots_age <- attributes(rcspline.eval(rdata$age, nk = 3))$knots
rdata$age3 <- rcs3_age

# Size of breast cancer
rcs3_size <- rcspline.eval(rdata$size, nk = 3)
attr(rcs3_size, "dim") <- NULL
attr(rcs3_size, "knots") <- NULL
pos_knots_size <- attributes(rcspline.eval(rdata$size, nk = 3))$knots
rdata$size3 <- rcs3_size

# FG model
dd <- datadist(rdata)
options(datadist = "dd")
fit_fg_rcs <- cph(Surv(Tstart, Tstop, status == 1) ~ 
                    rcs(age, pos_knots_age) + rcs(size, pos_knots_size) +
                    ncat + hr_status, 
                  weights = weight.cens,
                  x = T, 
                  y = T, 
                  surv = T, 
                  data = rdata.w1
)
P_fg_age_rcs <- Predict(fit_fg_rcs, "age")
P_fg_size_rcs <- Predict(fit_fg_rcs, "size")

plot(P_fg_age_rcs)
plot(P_fg_size_rcs)
options(datadist = NULL)

# Fit the Fine and Gray model assuming linear relation of
# the continuous predictors
dd <- datadist(rdata, adjto.cat = "first")
options(datadist = "dd")
fit_fg <- cph(Surv(Tstart, Tstop, status == 1) ~ age + size +
                ncat + hr_status,
              weights = weight.cens,
              x = T, 
              y = T, 
              surv = T, 
              data = rdata.w1
)
options(datadist = NULL)


res_AIC <- c("Without splines" = AIC(fit_fg), 
             "With splines" = AIC(fit_fg_rcs))
res_AIC

# The AIC and the graphical check suggested 
# a potential linear relation 
# between the continuous predictors (age and size) 
# and the event of interest (breast cancer recurrence).

### Checking proportional subdistribution hazards assumption ---------
zp_fg <- cox.zph(fit_fg, transform = "identity")

par(las = 1, xaxs = "i", yaxs = "i")
# c(bottom, left, top, right)
oldpar <- par(mfrow = c(2, 2), mar = c(2, 2, 3, 2))

sub_title <- c("Age", "Size", "Lymph node status", "HR status")
for (i in 1:4) {
  plot(zp_fg[i], 
       resid = F, 
       bty = "n", 
       xlim = c(0, 5))
  abline(0, 0, lty = 3)
  title(sub_title[i])
}
mtext("Fine and Gray", 
      side = 3, 
      line = -1, 
      outer = TRUE, 
      font = 2)
par(oldpar)

zp_fg$table


# The statistical tests showed 
# a potential violation of the proportional hazards assumption 
# for size of breast cancer.
# For simplicity we ignore this violation in the remainder.


### Model development - fit the risk prediction models ---------------------

fit_fg <- coxph(Surv(Tstart, Tstop, status == 1) ~
                age + size + ncat + hr_status, 
                  weights = weight.cens,
                  x = T, 
                  y = T, 
                  data = rdata.w1
)
# NOTE: rms:cph() can also be used as before



