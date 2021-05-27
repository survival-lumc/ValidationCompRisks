External validation of the performance of competing risks prediction
models: a guide through modern methods - Cause specific hazard models
================

-   [Steps](#steps)
    -   [Installing and loading packages and import
        data](#installing-and-loading-packages-and-import-data)
    -   [Descriptive statistics](#descriptive-statistics)
-   [Goal 1 - develop a competing risks prediction
    model](#goal-1---develop-a-competing-risks-prediction-model)
    -   [1.1 Cumulative incidence
        curves](#11-cumulative-incidence-curves)
    -   [1.2 Check non-linearity of continuous
        predictors](#12-check-non-linearity-of-continuous-predictors)
    -   [1.3 Checking proportional hazards
        assumption](#13-checking-proportional-hazards-assumption)
    -   [1.4 Examine the risk fit of the
        models](#14-examine-the-risk-fit-of-the-models)
-   [Goal 2 - Assessing performance of a competing risks prediction
    model](#goal-2---assessing-performance-of-a-competing-risks-prediction-model)
    -   [2.1 Overall prediction error](#21-overall-prediction-error)
    -   [2.2 Discrimination](#22-discrimination)
    -   [2.3 Calibration](#23-calibration)
        -   [2.3.1 Numerical summaries of
            calibration](#231-numerical-summaries-of-calibration)
        -   [2.3.2 Calibration plot](#232-calibration-plot)
-   [Goal 3 - Clinical utility](#goal-3---clinical-utility)
-   [Additional references](#additional-references)
-   [Reproducibility ticket](#reproducibility-ticket)

## Steps

The steps taken in this file are:  
1. To develop a competing risks prediction model cause specific hazards
approach;  
2. To assess the performance of the model in terms of calibration,
discrimination and overall prediction error;  
3. To assess the potential clinical utility the model using decision
curve analysis;

### Installing and loading packages and import data

The following libraries are needed to achieve the outlined goals, the
code chunk below will a) check whether you already have them installed,
b) install them for you if not already present, and c) load the packages
into the session.

``` r
# Use pacman to check whether packages are installed, if not load
if (!require("pacman")) install.packages("pacman")
library(pacman)

pacman::p_load(
  rio,
  survival,
  rms,
  mstate,
  sqldf,
  pec,
  riskRegression,
  survAUC,
  survivalROC,
  timeROC,
  plotrix,
  splines,
  knitr,
  table1,
  kableExtra,
  gtsummary,
  boot,
  tidyverse,
  rsample,
  gridExtra,
  webshot
)

rdata <- readRDS(here::here("Data/rdata.rds"))
vdata <- readRDS(here::here("Data/vdata.rds"))

rdata$hr_status <- relevel(rdata$hr_status, ref = "ER and/or PR +")
vdata$hr_status <- relevel(vdata$hr_status, ref = "ER and/or PR +")
```

We loaded the development data (rdata) and the validation data (vdata).
More details about development and validation data are provided in the
manuscript.

### Descriptive statistics

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
Characteristic
</th>
<th style="text-align:left;">
Development data, N = 1,000
</th>
<th style="text-align:left;">
Validation data, N = 1,000
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Age (years)
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
</tr>
<tr>
<td style="text-align:left;padding-left: 2em;" indentlevel="1">
Mean (SD)
</td>
<td style="text-align:left;">
75 (7)
</td>
<td style="text-align:left;">
77 (6)
</td>
</tr>
<tr>
<td style="text-align:left;padding-left: 2em;" indentlevel="1">
Median (Range)
</td>
<td style="text-align:left;">
74 (65, 95)
</td>
<td style="text-align:left;">
76 (70, 96)
</td>
</tr>
<tr>
<td style="text-align:left;">
Size (cm)
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
</tr>
<tr>
<td style="text-align:left;padding-left: 2em;" indentlevel="1">
Mean (SD)
</td>
<td style="text-align:left;">
2.29 (1.31)
</td>
<td style="text-align:left;">
2.13 (1.32)
</td>
</tr>
<tr>
<td style="text-align:left;padding-left: 2em;" indentlevel="1">
Median (Range)
</td>
<td style="text-align:left;">
2.00 (0.10, 8.50)
</td>
<td style="text-align:left;">
1.80 (0.09, 11.00)
</td>
</tr>
<tr>
<td style="text-align:left;">
Nodal status
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
</tr>
<tr>
<td style="text-align:left;padding-left: 2em;" indentlevel="1">
negative
</td>
<td style="text-align:left;">
642 (64%)
</td>
<td style="text-align:left;">
688 (69%)
</td>
</tr>
<tr>
<td style="text-align:left;padding-left: 2em;" indentlevel="1">
positive
</td>
<td style="text-align:left;">
358 (36%)
</td>
<td style="text-align:left;">
312 (31%)
</td>
</tr>
<tr>
<td style="text-align:left;">
Hormon receptor status
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
</tr>
<tr>
<td style="text-align:left;padding-left: 2em;" indentlevel="1">
ER and/or PR +
</td>
<td style="text-align:left;">
822 (82%)
</td>
<td style="text-align:left;">
857 (86%)
</td>
</tr>
<tr>
<td style="text-align:left;padding-left: 2em;" indentlevel="1">
ER-/PR-
</td>
<td style="text-align:left;">
178 (18%)
</td>
<td style="text-align:left;">
143 (14%)
</td>
</tr>
</tbody>
</table>

## Goal 1 - develop a competing risks prediction model

### 1.1 Cumulative incidence curves

First, we draw the cumulative incidence curves of breast cancer
recurrence.

``` r
# Expand datasets
# Create indicator variables for the outcome
rdata$status_num <- as.numeric(rdata$status) - 1
rdata$status1[rdata$status_num == 1] <- 1
rdata$status1[rdata$status_num != 1] <- 0
rdata$status2[rdata$status_num == 2] <- 2
rdata$status2[rdata$status_num != 2] <- 0
# Create indicator variables for the outcome
vdata$status_num <- as.numeric(vdata$status) - 1
vdata$status1[vdata$status_num == 1] <- 1
vdata$status1[vdata$status_num != 1] <- 0
vdata$status2[vdata$status_num == 2] <- 2
vdata$status2[vdata$status_num != 2] <- 0

# Expand data to prepare for fitting the model 
rdata.w <- crprep(
  Tstop = "time",
  status = "status_num",
  trans = c(1, 2),
  id = "id",
  keep = c("age", "size", "ncat", "hr_status"),
  data = rdata
)
# Save extended data with weights for recurrence (failcode=1)
# and non recurrence mortality (failcode=2)
rdata.w1 <- rdata.w %>% filter(failcode == 1)
rdata.w2 <- rdata.w %>% filter(failcode == 2)
vdata.w <- crprep(
  Tstop = "time",
  status = "status_num",
  trans = c(1, 2),
  id = "id",
  keep = c("age", "size", "ncat", "hr_status"),
  data = vdata
)
vdata.w1 <- vdata.w %>% filter(failcode == 1)
vdata.w2 <- vdata.w %>% filter(failcode == 2)
# Development set
mfit3 <- survfit(
  Surv(Tstart, Tstop, status == 1) ~ 1,
  data = rdata.w1, weights = weight.cens
)
mfit4 <- survfit(
  Surv(Tstart, Tstop, status == 1) ~ 1,
  data = vdata.w1, weights = weight.cens
)
par(xaxs = "i", yaxs = "i", las = 1)
oldpar <- par(mfrow = c(1, 2), mar = c(5, 5, 1, 1))
plot(mfit3,
     col = 1, lwd = 2,
     xlab = "Years since BC diagnosis",
     ylab = "Cumulative incidence", bty = "n",
     ylim = c(0, 0.25), xlim = c(0, 5), fun = "event", conf.int = TRUE
)
title("Development data")
plot(mfit4,
     col = 1, lwd = 2,
     xlab = "Years since BC diagnosis",
     ylab = "Cumulative incidence", bty = "n",
     ylim = c(0, 0.25), xlim = c(0, 5), fun = "event", conf.int = TRUE
)
title("Validation data")
```

<img src="imgs/Prediction_CSC/cuminc-1.png" width="672" style="display: block; margin: auto;" />

``` r
par(oldpar)
# Cumulative incidences
smfit3 <- summary(mfit3, times = c(1, 2, 3, 4, 5))
smfit4 <- summary(mfit4, times = c(1, 2, 3, 4, 5))
```

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="empty-cells: hide;border-bottom:hidden;" colspan="1">
</th>
<th style="border-bottom:hidden;padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="3">

<div style="border-bottom: 1px solid #ddd; padding-bottom: 5px; ">

Development data

</div>

</th>
<th style="border-bottom:hidden;padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="3">

<div style="border-bottom: 1px solid #ddd; padding-bottom: 5px; ">

Validation data

</div>

</th>
</tr>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:right;">
Estimate
</th>
<th style="text-align:right;">
Lower .95
</th>
<th style="text-align:right;">
Upper .95
</th>
<th style="text-align:right;">
Estimate
</th>
<th style="text-align:right;">
Lower .95
</th>
<th style="text-align:right;">
Upper .95
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
1-year
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.03
</td>
</tr>
<tr>
<td style="text-align:left;">
2-year
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.06
</td>
</tr>
<tr>
<td style="text-align:left;">
3-year
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.09
</td>
</tr>
<tr>
<td style="text-align:left;">
4-year
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.11
</td>
</tr>
<tr>
<td style="text-align:left;">
5-year
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.12
</td>
</tr>
</tbody>
</table>

The 5-year cumulative incidence of breast cancer recurrence was 14% (95%
CI: 11-16%), and 10% (95%CI: 8-12%)

### 1.2 Check non-linearity of continuous predictors

Here we investigate the potential non-linear relation between continuous
predictors (i.e. age and size) and the outcomes. We apply three-knot
restricted cubic splines using `rms::rcs()` function (details are given
in e.g. Frank Harrell’s book ‘Regression Model Strategies (second
edition)’, page 27.

``` r
# Models without splines
fit_csh <- CSC(Hist(time, status_num) ~ age + size +
  ncat + hr_status, data = rdata, fitter = "cph")
fit_csc1 <- fit_csh$models$`Cause 1`
fit_csc2 <- fit_csh$models$`Cause 2`
# Models with splines
dd <- datadist(rdata)
options(datadist = "dd")
# Recurrence
fit_csc1_rcs <- cph(Surv(time, status_num == 1) ~ rcs(age, 3) + rcs(size, 3) +
  ncat + hr_status, x = T, y = T, surv = T, data = rdata)
# print(fit_csc1_rcs)
# print(summary(fit_csc1_rcs))
# print(anova(fit_csc1_rcs))
P_csc1_age_rcs <- Predict(fit_csc1_rcs, "age")
P_csc1_size_rcs <- Predict(fit_csc1_rcs, "size")
options(datadist = NULL)
# Non-recurrence mortality
dd <- datadist(rdata)
options(datadist = "dd")
fit_csc2_rcs <- cph(Surv(time, status_num == 2) ~ rcs(age, 3) + rcs(size, 3) +
  ncat + hr_status, x = T, y = T, surv = T, data = rdata)
# print(fit_csc2_rcs)
# print(summary(fit_csc2_rcs))
# print(anova(fit_csc2_rcs))
P_csc2_age_rcs <- Predict(fit_csc2_rcs, "age")
P_csc2_size_rcs <- Predict(fit_csc2_rcs, "size")
options(datadist = NULL)
oldpar <- par(mfrow = c(2, 2), mar = c(5, 5, 1, 1))
par(xaxs = "i", yaxs = "i", las = 1)
plot(P_csc1_age_rcs$age, P_csc1_age_rcs$yhat,
  type = "l", lwd = 2, col = "blue", bty = "n",
  xlab = "Age at breast cancer diagnosis", ylab = "log Relative Hazard", ylim = c(-2, 2),
  xlim = c(65, 95)
)
polygon(c(P_csc1_age_rcs$age, rev(P_csc1_age_rcs$age)),
  c(P_csc1_age_rcs$lower, rev(P_csc1_age_rcs$upper)),
  col = "grey75",
  border = FALSE
)
par(new = TRUE)
plot(P_csc1_age_rcs$age, P_csc1_age_rcs$yhat,
  type = "l", lwd = 2, col = "blue", bty = "n",
  xlab = "Age at breast cancer diagnosis", ylab = "log Relative Hazard",
  ylim = c(-2, 2), xlim = c(65, 95)
)
title("Recurrence")
# CSC 1- size
par(xaxs = "i", yaxs = "i", las = 1)
plot(P_csc1_size_rcs$size, P_csc1_size_rcs$yhat,
  type = "l", lwd = 2, col = "blue", bty = "n",
  xlab = "Size of breast cancer", ylab = "log Relative Hazard", ylim = c(-2, 2),
  xlim = c(0, 7)
)
polygon(c(P_csc1_size_rcs$size, rev(P_csc1_size_rcs$size)),
  c(P_csc1_size_rcs$lower, rev(P_csc1_size_rcs$upper)),
  col = "grey75",
  border = FALSE
)
par(new = TRUE)
plot(P_csc1_size_rcs$size, P_csc1_size_rcs$yhat,
  type = "l", lwd = 2, col = "blue", bty = "n",
  xlab = "Size of breast cancer", ylab = "log Relative Hazard",
  ylim = c(-2, 2), xlim = c(0, 7)
)
title("Recurrence")
par(xaxs = "i", yaxs = "i", las = 1)
options(datadist = NULL)
# CSC 2- age
plot(P_csc2_age_rcs$age, P_csc2_age_rcs$yhat,
  type = "l", lwd = 2, col = "blue", bty = "n",
  xlab = "Age at breast cancer diagnosis", ylab = "log Relative Hazard", ylim = c(-2, 2),
  xlim = c(65, 95)
)
polygon(c(P_csc2_age_rcs$age, rev(P_csc2_age_rcs$age)),
  c(P_csc2_age_rcs$lower, rev(P_csc2_age_rcs$upper)),
  col = "grey75",
  border = FALSE
)
par(new = TRUE)
plot(P_csc2_age_rcs$age, P_csc2_age_rcs$yhat,
  type = "l", lwd = 2, col = "blue", bty = "n",
  xlab = "Age at breast cancer diagnosis", ylab = "log Relative Hazard",
  ylim = c(-2, 2), xlim = c(65, 95)
)
title("Non recurrence mortality")
# CSC 2 - size
par(xaxs = "i", yaxs = "i", las = 1)
plot(P_csc2_size_rcs$size, P_csc2_size_rcs$yhat,
  type = "l", lwd = 2, col = "blue", bty = "n",
  xlab = "Size of breast cancer", ylab = "log Relative Hazard", ylim = c(-2, 2),
  xlim = c(0, 7)
)
polygon(c(P_csc2_size_rcs$size, rev(P_csc2_size_rcs$size)),
  c(P_csc2_size_rcs$lower, rev(P_csc2_size_rcs$upper)),
  col = "grey75",
  border = FALSE
)
par(new = TRUE)
plot(P_csc2_size_rcs$size, P_csc2_size_rcs$yhat,
  type = "l", lwd = 2, col = "blue", bty = "n",
  xlab = "Size of breast cancer", ylab = "log Relative Hazard",
  ylim = c(-2, 2), xlim = c(0, 7)
)
title("Non recurrence mortality")
```

<img src="imgs/Prediction_CSC/ff-1.png" width="672" style="display: block; margin: auto;" />

``` r
options(datadist = NULL)
par(oldpar)
```

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:right;">
AIC without splines
</th>
<th style="text-align:right;">
AIC with splines
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Recurrence specific hazard
</td>
<td style="text-align:right;">
1779.375
</td>
<td style="text-align:right;">
1780.699
</td>
</tr>
<tr>
<td style="text-align:left;">
Non recurrence mortality
</td>
<td style="text-align:right;">
2572.079
</td>
<td style="text-align:right;">
2574.830
</td>
</tr>
</tbody>
</table>

Both the graphical comparison and the AIC comparison suggested no
relevant departure from linear relations between the continuous
predictors (age and size) and the cause-specific hazards (recurrence and
non-recurrence mortality).

### 1.3 Checking proportional hazards assumption

We now examine the fits further by checking the proportionality of the
cause-specific hazards of the models.

``` r
zp_csc1 <- cox.zph(fit_csc1, transform = "identity")
par(las = 1, xaxs = "i", yaxs = "i")
# c(bottom, left, top, right)
oldpar <- par(mfrow = c(2, 2), mar = c(5, 6.1, 3.1, 1))
sub_title <- c("Age", "Size", "Lymph node status", "HR status")
for (i in 1:4) {
  plot(zp_csc1[i], resid = F, bty = "n", xlim = c(0, 5))
  abline(0, 0, lty = 3)
  title(sub_title[i])
}
mtext("Recurrence", side = 3, line = -1, outer = TRUE, font = 2)
```

<img src="imgs/Prediction_CSC/ph_csc1-1.png" width="672" style="display: block; margin: auto;" />

``` r
par(oldpar)
kable(round(zp_csc1$table, 3)) %>%
  kable_styling("striped", position = "center")
```

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:right;">
chisq
</th>
<th style="text-align:right;">
df
</th>
<th style="text-align:right;">
p
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
age
</td>
<td style="text-align:right;">
0.458
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.498
</td>
</tr>
<tr>
<td style="text-align:left;">
size
</td>
<td style="text-align:right;">
2.309
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.129
</td>
</tr>
<tr>
<td style="text-align:left;">
ncat
</td>
<td style="text-align:right;">
0.016
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.900
</td>
</tr>
<tr>
<td style="text-align:left;">
hr\_status
</td>
<td style="text-align:right;">
0.078
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.780
</td>
</tr>
<tr>
<td style="text-align:left;">
GLOBAL
</td>
<td style="text-align:right;">
2.713
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.607
</td>
</tr>
</tbody>
</table>

``` r
zp_csc2 <- cox.zph(fit_csc2, transform = "identity")
par(las = 1, xaxs = "i", yaxs = "i")
# c(bottom, left, top, right)
oldpar <- par(mfrow = c(2, 2), mar = c(5, 6.1, 3.1, 1))
sub_title <- c("Age", "Size", "Lymph node status", "HR status")
for (i in 1:4) {
  plot(zp_csc2[i], resid = F, bty = "n", xlim = c(0, 5))
  abline(0, 0, lty = 3)
  title(sub_title[i])
}
mtext("Non-recurrence mortality", side = 3, line = -1, outer = TRUE, font = 2)
```

<img src="imgs/Prediction_CSC/ph_csc2-1.png" width="672" style="display: block; margin: auto;" />

``` r
par(oldpar)
kable(round(zp_csc2$table, 3)) %>%
  kable_styling("striped", position = "center")
```

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:right;">
chisq
</th>
<th style="text-align:right;">
df
</th>
<th style="text-align:right;">
p
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
age
</td>
<td style="text-align:right;">
2.727
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.099
</td>
</tr>
<tr>
<td style="text-align:left;">
size
</td>
<td style="text-align:right;">
1.716
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.190
</td>
</tr>
<tr>
<td style="text-align:left;">
ncat
</td>
<td style="text-align:right;">
5.227
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.022
</td>
</tr>
<tr>
<td style="text-align:left;">
hr\_status
</td>
<td style="text-align:right;">
0.020
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.889
</td>
</tr>
<tr>
<td style="text-align:left;">
GLOBAL
</td>
<td style="text-align:right;">
9.374
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.052
</td>
</tr>
</tbody>
</table>

The statistical tests showed a potential violation of the proportional
hazards assumption for nodal status in the model for non-recurrence
mortality. For simplicity we ignore this violation in the remainder.

### 1.4 Examine the risk fit of the models

We show the results of the Cox proportional cause-specific hazards
regression models

-   Cox proportional hazard model for recurrence

``` r
dd <- datadist(rdata)
options(datadist = "dd")
options(prType = "html")
fit_csc1_cph <- cph(Surv(time, status_num == 1) ~ age + size +
  ncat + hr_status,
x = T, y = T, surv = T, data = rdata
)
print(fit_csc1_cph)
```

 <strong>Cox Proportional Hazards Model</strong>
 
 <pre>
 cph(formula = Surv(time, status_num == 1) ~ age + size + ncat + 
     hr_status, data = rdata, x = T, y = T, surv = T)
 </pre>
 
 <table class='gmisc_table' style='border-collapse: collapse; margin-top: 1em; margin-bottom: 1em;' >
<thead>
<tr>
<th style='font-weight: 900; border-bottom: 1px solid grey; border-top: 2px solid grey; border-left: 1px solid black; border-right: 1px solid black; text-align: center;'></th>
<th style='font-weight: 900; border-bottom: 1px solid grey; border-top: 2px solid grey; border-right: 1px solid black; text-align: center;'>Model Tests</th>
<th style='font-weight: 900; border-bottom: 1px solid grey; border-top: 2px solid grey; border-right: 1px solid black; text-align: center;'>Discrimination<br>Indexes</th>
</tr>
</thead>
<tbody>
<tr>
<td style='min-width: 9em; border-left: 1px solid black; border-right: 1px solid black; text-align: center;'>Obs 1000</td>
<td style='min-width: 9em; border-right: 1px solid black; text-align: center;'>LR χ<sup>2</sup> 46.33</td>
<td style='min-width: 9em; border-right: 1px solid black; text-align: center;'><i>R</i><sup>2</sup> 0.054</td>
</tr>
<tr>
<td style='min-width: 9em; border-left: 1px solid black; border-right: 1px solid black; text-align: center;'>Events 135</td>
<td style='min-width: 9em; border-right: 1px solid black; text-align: center;'>d.f. 4</td>
<td style='min-width: 9em; border-right: 1px solid black; text-align: center;'><i>D</i><sub>xy</sub> 0.326</td>
</tr>
<tr>
<td style='min-width: 9em; border-left: 1px solid black; border-right: 1px solid black; text-align: center;'>Center 2.0201</td>
<td style='min-width: 9em; border-right: 1px solid black; text-align: center;'>Pr(>χ<sup>2</sup>) 0.0000</td>
<td style='min-width: 9em; border-right: 1px solid black; text-align: center;'><i>g</i> 0.627</td>
</tr>
<tr>
<td style='min-width: 9em; border-left: 1px solid black; border-right: 1px solid black; text-align: center;'></td>
<td style='min-width: 9em; border-right: 1px solid black; text-align: center;'>Score χ<sup>2</sup> 53.04</td>
<td style='min-width: 9em; border-right: 1px solid black; text-align: center;'><i>g</i><sub>r</sub> 1.871</td>
</tr>
<tr>
<td style='min-width: 9em; border-bottom: 2px solid grey; border-left: 1px solid black; border-right: 1px solid black; text-align: center;'></td>
<td style='min-width: 9em; border-bottom: 2px solid grey; border-right: 1px solid black; text-align: center;'>Pr(>χ<sup>2</sup>) 0.0000</td>
<td style='min-width: 9em; border-bottom: 2px solid grey; border-right: 1px solid black; text-align: center;'></td>
</tr>
</tbody>
</table>

 
 <table class='gmisc_table' style='border-collapse: collapse; margin-top: 1em; margin-bottom: 1em;' >
<thead>
<tr><th style='border-bottom: 1px solid grey; font-weight: 900; border-top: 2px solid grey; min-width: 7em; text-align: center;'></th>
<th style='font-weight: 900; border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: right;'>β</th>
<th style='font-weight: 900; border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: right;'>S.E.</th>
<th style='font-weight: 900; border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: right;'>Wald <i>Z</i></th>
<th style='font-weight: 900; border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: right;'>Pr(>|<i>Z</i>|)</th>
</tr>
</thead>
<tbody>
<tr>
<td style='min-width: 7em; text-align: left;'>age</td>
<td style='min-width: 7em; text-align: right;'> 0.0153</td>
<td style='min-width: 7em; text-align: right;'> 0.0125</td>
<td style='min-width: 7em; text-align: right;'>1.22</td>
<td style='min-width: 7em; text-align: right;'>0.2237</td>
</tr>
<tr>
<td style='min-width: 7em; text-align: left;'>size</td>
<td style='min-width: 7em; text-align: right;'> 0.2510</td>
<td style='min-width: 7em; text-align: right;'> 0.0564</td>
<td style='min-width: 7em; text-align: right;'>4.45</td>
<td style='min-width: 7em; text-align: right;'><0.0001</td>
</tr>
<tr>
<td style='min-width: 7em; text-align: left;'>ncat=positive</td>
<td style='min-width: 7em; text-align: right;'> 0.5090</td>
<td style='min-width: 7em; text-align: right;'> 0.1769</td>
<td style='min-width: 7em; text-align: right;'>2.88</td>
<td style='min-width: 7em; text-align: right;'>0.0040</td>
</tr>
<tr>
<td style='min-width: 7em; border-bottom: 2px solid grey; text-align: left;'>hr_status=ER-/PR-</td>
<td style='min-width: 7em; border-bottom: 2px solid grey; text-align: right;'> 0.6433</td>
<td style='min-width: 7em; border-bottom: 2px solid grey; text-align: right;'> 0.1926</td>
<td style='min-width: 7em; border-bottom: 2px solid grey; text-align: right;'>3.34</td>
<td style='min-width: 7em; border-bottom: 2px solid grey; text-align: right;'>0.0008</td>
</tr>
</tbody>
</table>

``` r
# print(summary(fit_csc1_cph))
options(datadist = NULL)
```

-   Cox proportional non recurrence mortality-specific hazard model

``` r
dd <- datadist(rdata)
options(datadist = "dd")
options(prType = "html")
fit_csc2_cph <- cph(Surv(time, status_num == 2) ~ age + size +
  ncat + hr_status,
x = T, y = T, surv = T, data = rdata
)
print(fit_csc2_cph)
```

 <strong>Cox Proportional Hazards Model</strong>
 
 <pre>
 cph(formula = Surv(time, status_num == 2) ~ age + size + ncat + 
     hr_status, data = rdata, x = T, y = T, surv = T)
 </pre>
 
 <table class='gmisc_table' style='border-collapse: collapse; margin-top: 1em; margin-bottom: 1em;' >
<thead>
<tr>
<th style='font-weight: 900; border-bottom: 1px solid grey; border-top: 2px solid grey; border-left: 1px solid black; border-right: 1px solid black; text-align: center;'></th>
<th style='font-weight: 900; border-bottom: 1px solid grey; border-top: 2px solid grey; border-right: 1px solid black; text-align: center;'>Model Tests</th>
<th style='font-weight: 900; border-bottom: 1px solid grey; border-top: 2px solid grey; border-right: 1px solid black; text-align: center;'>Discrimination<br>Indexes</th>
</tr>
</thead>
<tbody>
<tr>
<td style='min-width: 9em; border-left: 1px solid black; border-right: 1px solid black; text-align: center;'>Obs 1000</td>
<td style='min-width: 9em; border-right: 1px solid black; text-align: center;'>LR χ<sup>2</sup> 170.94</td>
<td style='min-width: 9em; border-right: 1px solid black; text-align: center;'><i>R</i><sup>2</sup> 0.168</td>
</tr>
<tr>
<td style='min-width: 9em; border-left: 1px solid black; border-right: 1px solid black; text-align: center;'>Events 204</td>
<td style='min-width: 9em; border-right: 1px solid black; text-align: center;'>d.f. 4</td>
<td style='min-width: 9em; border-right: 1px solid black; text-align: center;'><i>D</i><sub>xy</sub> 0.491</td>
</tr>
<tr>
<td style='min-width: 9em; border-left: 1px solid black; border-right: 1px solid black; text-align: center;'>Center 9.06</td>
<td style='min-width: 9em; border-right: 1px solid black; text-align: center;'>Pr(>χ<sup>2</sup>) 0.0000</td>
<td style='min-width: 9em; border-right: 1px solid black; text-align: center;'><i>g</i> 1.047</td>
</tr>
<tr>
<td style='min-width: 9em; border-left: 1px solid black; border-right: 1px solid black; text-align: center;'></td>
<td style='min-width: 9em; border-right: 1px solid black; text-align: center;'>Score χ<sup>2</sup> 186.60</td>
<td style='min-width: 9em; border-right: 1px solid black; text-align: center;'><i>g</i><sub>r</sub> 2.849</td>
</tr>
<tr>
<td style='min-width: 9em; border-bottom: 2px solid grey; border-left: 1px solid black; border-right: 1px solid black; text-align: center;'></td>
<td style='min-width: 9em; border-bottom: 2px solid grey; border-right: 1px solid black; text-align: center;'>Pr(>χ<sup>2</sup>) 0.0000</td>
<td style='min-width: 9em; border-bottom: 2px solid grey; border-right: 1px solid black; text-align: center;'></td>
</tr>
</tbody>
</table>

 
 <table class='gmisc_table' style='border-collapse: collapse; margin-top: 1em; margin-bottom: 1em;' >
<thead>
<tr><th style='border-bottom: 1px solid grey; font-weight: 900; border-top: 2px solid grey; min-width: 7em; text-align: center;'></th>
<th style='font-weight: 900; border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: right;'>β</th>
<th style='font-weight: 900; border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: right;'>S.E.</th>
<th style='font-weight: 900; border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: right;'>Wald <i>Z</i></th>
<th style='font-weight: 900; border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: right;'>Pr(>|<i>Z</i>|)</th>
</tr>
</thead>
<tbody>
<tr>
<td style='min-width: 7em; text-align: left;'>age</td>
<td style='min-width: 7em; text-align: right;'> 0.1118</td>
<td style='min-width: 7em; text-align: right;'> 0.0100</td>
<td style='min-width: 7em; text-align: right;'>11.15</td>
<td style='min-width: 7em; text-align: right;'><0.0001</td>
</tr>
<tr>
<td style='min-width: 7em; text-align: left;'>size</td>
<td style='min-width: 7em; text-align: right;'> 0.2361</td>
<td style='min-width: 7em; text-align: right;'> 0.0481</td>
<td style='min-width: 7em; text-align: right;'> 4.91</td>
<td style='min-width: 7em; text-align: right;'><0.0001</td>
</tr>
<tr>
<td style='min-width: 7em; text-align: left;'>ncat=positive</td>
<td style='min-width: 7em; text-align: right;'> 0.1860</td>
<td style='min-width: 7em; text-align: right;'> 0.1442</td>
<td style='min-width: 7em; text-align: right;'> 1.29</td>
<td style='min-width: 7em; text-align: right;'>0.1971</td>
</tr>
<tr>
<td style='min-width: 7em; border-bottom: 2px solid grey; text-align: left;'>hr_status=ER-/PR-</td>
<td style='min-width: 7em; border-bottom: 2px solid grey; text-align: right;'> 0.2396</td>
<td style='min-width: 7em; border-bottom: 2px solid grey; text-align: right;'> 0.1765</td>
<td style='min-width: 7em; border-bottom: 2px solid grey; text-align: right;'> 1.36</td>
<td style='min-width: 7em; border-bottom: 2px solid grey; text-align: right;'>0.1745</td>
</tr>
</tbody>
</table>

``` r
# print(summary(fit_csc2_cph))
options(datadist = NULL)
```

The coefficients of the models indicated that larger tumor size,
positive nodal status and negative hormone receptor status status were
associated with higher risk to develop a breast cancer recurrence, while
older patients and larger tumors are associated with higher risk of non
recurrence mortality.

## Goal 2 - Assessing performance of a competing risks prediction model

Here we evaluate the performance of the risk prediction models in terms
of calibration, discrimination and overall prediction error.

### 2.1 Overall prediction error

We calculate the Brier Score, and the scaled Brier scale and the
corresponding confidence intervals..

Some confidence intervals are calculated using the bootstrap percentile
method.

``` r
# Bootstrapping data
set.seed(20201214)
rboot <- bootstraps(rdata, times = 10)
vboot <- bootstraps(vdata, times = 10)
# NOTE: B=10 to speed up the procedure, is typically set to 100 or 1000
```

``` r
# riskRegression::Score() to calculate Brier and scaled Brier (in this function called "ipa")
# Development set - apparent validation
score_rdata1 <- Score(list("CSH development" = fit_csh),
  formula = Hist(time, status_num) ~ 1,
  data = rdata, conf.int = TRUE, times = 4.99,
  cens.model = "km", metrics = "brier",
  summary = "ipa", cause = 1
)
# Validation set - external validation
score_vdata1 <-
  Score(list("CSH validation" = fit_csh),
    formula = Hist(time, status_num) ~ 1,
    data = vdata, conf.int = TRUE, times = 4.99,
    cens.model = "km", metrics = "brier",
    summary = "ipa", cause = 1
  )
# Development set - internal validation (bootstrap)
# mstate::crprep() for every bootstrap sample
crprep_boot <- function(split) {
  crprep(
    Tstop = "time", status = "status_num",
    trans = 1, data = analysis(split),
    keep = c(
      "status_num", "age", "size",
      "ncat", "hr_status"
    )
  )
}
# riskRegression::Score() to calculate Brier and scaled Brier (here called IPA) for every bootstrap sample
score_boot_1 <- function(split) {
  Score(list("CSH" = fit_csh),
    formula = Hist(time, status_num) ~ 1,
    data = analysis(split), conf.int = FALSE, times = 4.99,
    cens.model = "km", metrics = "brier", cause = 1,
    summary = "ipa"
  )$Brier[[1]]$IPA[2]
}
# Development data
rboot <- rboot %>% mutate(
  cr.prep = map(splits, crprep_boot),
  IPA1 = map_dbl(splits, score_boot_1)
)
# Validation data
vboot <- vboot %>% mutate(
  cr.prep = map(splits, crprep_boot),
  IPA1 = map_dbl(splits, score_boot_1),
)
```

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="empty-cells: hide;border-bottom:hidden;" colspan="1">
</th>
<th style="border-bottom:hidden;padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="3">

<div style="border-bottom: 1px solid #ddd; padding-bottom: 5px; ">

Development data

</div>

</th>
<th style="border-bottom:hidden;padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="3">

<div style="border-bottom: 1px solid #ddd; padding-bottom: 5px; ">

Validation data

</div>

</th>
</tr>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:right;">
Estimate
</th>
<th style="text-align:right;">
Lower .95
</th>
<th style="text-align:right;">
Upper .95
</th>
<th style="text-align:right;">
Estimate
</th>
<th style="text-align:right;">
Lower .95
</th>
<th style="text-align:right;">
Upper .95
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Brier
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.13
</td>
</tr>
<tr>
<td style="text-align:left;">
scaled Brier
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.08
</td>
</tr>
</tbody>
</table>

Note: unexpectedly, the point estimate for the Brier score is lower
(thus better) and for the scaled Brier score is higher (thus better) in
the validation data compared to the development data.

### 2.2 Discrimination

We here calculate

-   The 5-year C-index. More details are in the main manuscript and its
    references;
-   The 5-year time-dependent AUC. More details are in the manuscript
    and in its references;

We used the time horizon up to 4.99 and not 5 years since controls are
considered patients at risk after the time horizon.

``` r
# C-index
# Development set (Apparent validation)
C_rdata1_cph1 <- unlist(pec::cindex(fit_csh,
                                    cause = 1,
                                    eval.times = 4.99
)$AppCindex)
# Validdation set
C_vdata1_cph1 <- unlist(pec::cindex(fit_csh,
                                    data = vdata,
                                    cause = 1, eval.times = 4.99
)$AppCindex)

# 5-year time dependent AUC
# Development set (Apparent validation)
Uno_rdata1_CSC <-
  timeROC(
    T = rdata$time, delta = rdata$status1,
    marker = predictRisk(fit_csh, newdata = rdata, cause = 1, times = 5),
    cause = 1, weighting = "marginal", times = 4.99,
    iid = TRUE
  )
# Validdation set
Uno_vdata1_CSC <-
  timeROC(
    T = vdata$time, delta = vdata$status1,
    marker = predictRisk(fit_csh, newdata = vdata, cause = 1, times = 5),
    cause = 1, weighting = "marginal", times = 4.99,
    iid = TRUE
  )

# NOTE: if you have many observations (n > 2000), standard error computation may be really long.
# In that case, you may consider using bootstrapping to calculate confidence intervals.
# NOTE: AUC_1: controls = subjects free of any event 
# NOTE: AUC_2: controls = subjects does not experience the primary event, this is what we use here 

# Bootstraping Wolbers' C-index to calculate the bootstrap percentile confidence intervals
C_boot1_cph1 <- function(split) {
  unlist(pec::cindex(fit_csh,
                     data = analysis(split),
                     cause = 1, eval.times = 4.99
  )$AppCindex)
}
C_boot1_cph2 <- function(split) {
  unlist(pec::cindex(fit_csh,
                     data = analysis(split),
                     cause = 2, eval.times = 4.99
  )$AppCindex)
}
# Run time-dependent AUC in the bootstrapped development and validation data
# to calculate the non-parametric CI through percentile bootstrap
rboot <- rboot %>% mutate(
  C1 = map_dbl(splits, C_boot1_cph1),
  C2 = map_dbl(splits, C_boot1_cph2)
)
vboot <- vboot %>% mutate(
  C1 = map_dbl(splits, C_boot1_cph1),
  C2 = map_dbl(splits, C_boot1_cph2)
)
```

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="empty-cells: hide;border-bottom:hidden;" colspan="1">
</th>
<th style="border-bottom:hidden;padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="3">

<div style="border-bottom: 1px solid #ddd; padding-bottom: 5px; ">

Development data

</div>

</th>
<th style="border-bottom:hidden;padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="3">

<div style="border-bottom: 1px solid #ddd; padding-bottom: 5px; ">

Validation data

</div>

</th>
</tr>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:right;">
Estimate
</th>
<th style="text-align:right;">
Lower .95
</th>
<th style="text-align:right;">
Upper .95
</th>
<th style="text-align:right;">
Estimate
</th>
<th style="text-align:right;">
Lower .95
</th>
<th style="text-align:right;">
Upper .95
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Wolbers C
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
0.73
</td>
</tr>
<tr>
<td style="text-align:left;">
Uno AUC
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
0.63
</td>
<td style="text-align:right;">
0.73
</td>
<td style="text-align:right;">
0.74
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
0.79
</td>
</tr>
</tbody>
</table>

The time-dependent AUC at 5 years was 0.74 in the validation set.

### 2.3 Calibration

We assess calibration by:

-   The calibration plot as a graphical representation of calibration;

-   The observed vs expected ratio (O/E ratio);

-   The squared bias, i.e., the average squared difference between
    actual risks and risk predictions;

-   The integrated Calibration Index (ICI), i.e., the average absolute
    difference between actual risks and risk predictions;

-   E50, E90 and Emax denote the median, 90th percentile and the maximum
    of the absolute differences between actual risks and risk
    predictions;

#### 2.3.1 Numerical summaries of calibration

We calculate the O/E ratio, squared bias, ICI, E50, E90 and Emax at 5
years in the development and validation data.

``` r
# Load the function to calculate the OE ratio
source(here::here("R/OE_function.R"))
# O = estimated cumulative incidence at 5 years
# E = mean of the predicted cumulative incidence at 5 years
Po_t <- summary(
  survfit(Surv(Tstart, Tstop, status == 1) ~ 1,
          data = vdata.w1, weights = weight.cens
  ),
  times = 5
)
obs_vdata <- 1 - Po_t$surv
obs_stderror <- Po_t$std.err
# Observed/Expected ratio
OE_vdata <- OE_function(
  fit = fit_csh, newdata = vdata, cause = 1,
  thorizon = 5, obs_cif = obs_vdata,
  std.error = obs_stderror
)
res_OE <- matrix(OE_vdata,
                 ncol = 3, nrow = 1, byrow = T,
                 dimnames = list(
                   c("O/E ratio"),
                   c("Estimate", "Lower.95", "Upper.95")
                 )
)
kable(res_OE) %>%
  kable_styling("striped", position = "center")
```

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:right;">
Estimate
</th>
<th style="text-align:right;">
Lower.95
</th>
<th style="text-align:right;">
Upper.95
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
O/E ratio
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
0.99
</td>
</tr>
</tbody>
</table>

The competing risks prediction model slightly overestimates the absolute
risk to develop breast cancer recurrence in the validation data.

``` r
# Calibration measures: squared bias, ICI, E50, E90, Emax
source(here::here("R/cal_measures.R"))
calmeas_vdata <- cal_measures(vdata, 5, fit_csh,
                              Tstop = "time", status = "status_num", cause = 1
)
# Squared bias
avg_sqbias_CSC <- mean((predictRisk(fit_csh, newdata = vdata, cause = 1, times = 5)
                        - obs_vdata)**2)
res_calmeas <- matrix(c(avg_sqbias_CSC, calmeas_vdata),
                      ncol = 1, nrow = 5, byrow = T,
                      dimnames = list(
                        c("Average squared bias", "ICI", "E50", "E90", "Emax"),
                        c("Estimate")
                      )
)
res_calmeas <- round(res_calmeas, 2)
kable(res_calmeas) %>%
  kable_styling("striped", position = "center")
```

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:right;">
Estimate
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Average squared bias
</td>
<td style="text-align:right;">
0.01
</td>
</tr>
<tr>
<td style="text-align:left;">
ICI
</td>
<td style="text-align:right;">
0.02
</td>
</tr>
<tr>
<td style="text-align:left;">
E50
</td>
<td style="text-align:right;">
0.03
</td>
</tr>
<tr>
<td style="text-align:left;">
E90
</td>
<td style="text-align:right;">
0.04
</td>
</tr>
<tr>
<td style="text-align:left;">
Emax
</td>
<td style="text-align:right;">
0.04
</td>
</tr>
</tbody>
</table>

#### 2.3.2 Calibration plot

Calibration plot for the validation data is calculated using
pseudo-values.

Calibration plots reports:

-   on the *x-axis* the estimated risk by the prediction model by a
    fixed time point (e.g. at 5 years);
-   on the *y-axis* the estimated actual risk by a fixed time point
    (e.g. at 5 years);
-   The 45-degree line indicates perfect calibration. Points below the
    45-degree line indicate that the model overestimates the estimated
    actual risk. If points are above the 45-degree line, the model
    underestimates the estimated actual risk.

``` r
x <- Score(list(model1 = fit_csh), Hist(time, status_num) ~ 1,
           data = vdata,
           cause = 1, times = 5, plots = "cal"
)
oldpar <- par(
  mar = c(5.1, 5.8, 4.1, 2.1), mgp = c(4.25, 1, 0),
  xaxs = "i", yaxs = "i", las = 1
)
plotCalibration(x,
                brier.in.legend = FALSE,
                auc.in.legend = FALSE, cens.method = "pseudo",
                cex = 1, xlim = c(0, 0.5), ylim = c(0, 0.5), rug=TRUE
)
title("Cause-specific hazards models")
```

<img src="imgs/Prediction_CSC/cal_rcs-1.png" width="672" style="display: block; margin: auto;" />

``` r
par(oldpar)
```

Calibration plot suggests that the prediction model seems to
overestimate the actual risk, especially at the lower and higher values
of the estimated risk.

## Goal 3 - Clinical utility

Clinical utility can be measured by the net benefit and plotted in a
decision curve. Details about net benefit, decision curve calculation
and interpretation are provided in the manuscript (see also the
appendix) and its references.

``` r
# Run the stdca function to calculate the net benefit and the elements needed to develop decision curve analysis
source(here::here("R/stdca.R"))
# Development data
# calculation estimated risk
rdata$pred5 <- predictRisk(fit_csh, newdata = rdata, times = 5)
rdata <- as.data.frame(rdata)
dca_rdata_1 <- stdca(
  data = rdata, outcome = "status_num", ttoutcome = "time",
  timepoint = 5, predictors = "pred5", xstop = 0.35,
  ymin = -0.01, graph = FALSE, cmprsk = TRUE
)
# Decision curves plot
oldpar <- par(xaxs = "i", yaxs = "i", las = 1, mar = c(6.1, 5.8, 4.1, 2.1), mgp = c(4.25, 1, 0))
plot(dca_rdata_1$net.benefit$threshold,
     dca_rdata_1$net.benefit$pred5,
     type = "l", lwd = 2, lty = 1,
     xlab = "", ylab = "Net Benefit",
     xlim = c(0, 0.5), ylim = c(-0.10, 0.10), bty = "n", xaxt = "n"
)
legend("topright", c("Treat all", "Treat none", "Prediction model"),
       lwd = c(2, 2, 2), lty = c(1, 2, 1), col = c("darkgray", "black", "black"), bty = "n"
)
lines(dca_rdata_1$net.benefit$threshold, dca_rdata_1$net.benefit$none,
      type = "l", lwd = 2, lty = 4
)
lines(dca_rdata_1$net.benefit$threshold, dca_rdata_1$net.benefit$all,
      type = "l", lwd = 2, col = "darkgray"
)
axis(1, at = c(0, 0.1, 0.2, 0.3, 0.4, 0.5))
axis(1,
     pos = -0.145, at = c(0.1, 0.2, 0.3, 0.4, 0.5),
     labels = c("1:9", "1:4", "3:7", "2:3", "1:1")
)
mtext("Threshold probability", 1, line = 2)
mtext("Harm to benefit ratio", 1, line = 5)
title("Development data")
```

<img src="imgs/Prediction_CSC/dca-1.png" width="672" style="display: block; margin: auto;" />

``` r
par(oldpar)
# Validation data
# Predicted probability calculation
vdata$pred5 <- predictRisk(fit_csh, newdata = vdata, times = 5)
vdata <- as.data.frame(vdata)
# Run decision curve analysis
# Development data
# Model without PGR
dca_vdata_1 <- stdca(
  data = vdata, outcome = "status_num", ttoutcome = "time",
  timepoint = 5, predictors = "pred5", xstop = 0.45,
  ymin = -0.01, graph = FALSE, cmprsk = TRUE
)
# Decision curves plot
oldpar <- par(xaxs = "i", yaxs = "i", las = 1, mar = c(6.1, 5.8, 4.1, 2.1), mgp = c(4.25, 1, 0))
plot(dca_vdata_1$net.benefit$threshold,
     dca_vdata_1$net.benefit$pred5,
     type = "l", lwd = 2, lty = 1,
     xlab = "", ylab = "Net Benefit",
     xlim = c(0, 0.5), ylim = c(-0.10, 0.10), bty = "n", xaxt = "n"
)
lines(dca_vdata_1$net.benefit$threshold,
      dca_vdata_1$net.benefit$none,
      type = "l", lwd = 2, lty = 4
)
lines(dca_vdata_1$net.benefit$threshold,
      dca_vdata_1$net.benefit$all,
      type = "l", lwd = 2, col = "darkgray"
)
legend("topright", c("Treat all", "Treat none", "Prediction model"),
       lwd = c(2, 2, 2), lty = c(1, 2, 1), col = c("darkgray", "black", "black"), bty = "n"
)
axis(1, at = c(0, 0.1, 0.2, 0.3, 0.4, 0.5))
axis(1,
     pos = -0.145, at = c(0.1, 0.2, 0.3, 0.4, 0.5),
     labels = c("1:9", "1:4", "3:7", "2:3", "1:1")
)
mtext("Threshold probability", 1, line = 2)
mtext("Harm to benefit ratio", 1, line = 5)
title("Validation data")
```

<img src="imgs/Prediction_CSC/dca-2.png" width="672" style="display: block; margin: auto;" />

``` r
par(oldpar)
```

If we choose a threshold of 20%, the model had a net benefit of 0.011 in
the development data. This means that the model would identify 11
patients per 1000 who will have beast cancer recurrence within 5 years
since diagnosis where adjuvant chemotherapy is really needed. In the
validation data, the model had a net benefit of 0.014 choosing a
threshold of 20%.

## Additional references

-   Calibration  
    <https://onlinelibrary.wiley.com/doi/abs/10.1002/sim.6152>  
    <https://onlinelibrary.wiley.com/doi/full/10.1002/sim.8570>  

-   Discrimination  
    <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4059461>  
    <https://onlinelibrary.wiley.com/doi/10.1002/sim.5958>  

-   Overall prediction error
    <https://onlinelibrary.wiley.com/doi/abs/10.1002/bimj.201000073>  
    <https://diagnprognres.biomedcentral.com/articles/10.1186/s41512-018-0029-2>  
    R Vignette:
    <https://cran.r-project.org/web/packages/riskRegression/vignettes/IPA.html>  

-   Clinical utility (decision curves)  
    R/SAS/STATA code and references:
    <https://www.mskcc.org/departments/epidemiology-biostatistics/biostatistics/decision-curve-analysis>  
    More guidelines about net benefit assessment and interpretation  
    <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6261531/>  
    <https://diagnprognres.biomedcentral.com/articles/10.1186/s41512-019-0064-7>  

## Reproducibility ticket

``` r
sessionInfo()
```

    ## R version 4.0.4 (2021-02-15)
    ## Platform: x86_64-w64-mingw32/x64 (64-bit)
    ## Running under: Windows 10 x64 (build 18363)
    ## 
    ## Matrix products: default
    ## 
    ## locale:
    ## [1] LC_COLLATE=Dutch_Netherlands.1252  LC_CTYPE=Dutch_Netherlands.1252   
    ## [3] LC_MONETARY=Dutch_Netherlands.1252 LC_NUMERIC=C                      
    ## [5] LC_TIME=Dutch_Netherlands.1252    
    ## 
    ## attached base packages:
    ## [1] splines   stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] webshot_0.5.2             gridExtra_2.3            
    ##  [3] rsample_0.1.0             forcats_0.5.1            
    ##  [5] stringr_1.4.0             dplyr_1.0.6              
    ##  [7] purrr_0.3.4               readr_1.4.0              
    ##  [9] tidyr_1.1.3               tibble_3.1.1             
    ## [11] tidyverse_1.3.1           boot_1.3-28              
    ## [13] gtsummary_1.4.0           kableExtra_1.3.4         
    ## [15] table1_1.4                knitr_1.33               
    ## [17] plotrix_3.8-1             timeROC_0.4              
    ## [19] survivalROC_1.0.3         survAUC_1.0-5            
    ## [21] riskRegression_2020.12.08 pec_2020.11.17           
    ## [23] prodlim_2019.11.13        sqldf_0.4-11             
    ## [25] RSQLite_2.2.7             gsubfn_0.7               
    ## [27] proto_1.0.0               mstate_0.3.1             
    ## [29] rms_6.2-0                 SparseM_1.81             
    ## [31] Hmisc_4.5-0               ggplot2_3.3.3            
    ## [33] Formula_1.2-4             lattice_0.20-44          
    ## [35] survival_3.2-11           rio_0.5.26               
    ## [37] pacman_0.5.1             
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] readxl_1.3.1        backports_1.2.1     systemfonts_1.0.2  
    ##   [4] listenv_0.8.0       TH.data_1.0-10      digest_0.6.27      
    ##   [7] foreach_1.5.1       htmltools_0.5.1.1   fansi_0.4.2        
    ##  [10] magrittr_2.0.1      checkmate_2.0.0     memoise_2.0.0      
    ##  [13] cluster_2.1.2       openxlsx_4.2.3      globals_0.14.0     
    ##  [16] modelr_0.1.8        mets_1.2.8.1        matrixStats_0.58.0 
    ##  [19] sandwich_3.0-0      svglite_2.0.0       jpeg_0.1-8.1       
    ##  [22] colorspace_2.0-1    blob_1.2.1          rvest_1.0.0        
    ##  [25] haven_2.4.1         xfun_0.22           tcltk_4.0.4        
    ##  [28] crayon_1.4.1        jsonlite_1.7.2      zoo_1.8-9          
    ##  [31] iterators_1.0.13    glue_1.4.2          gtable_0.3.0       
    ##  [34] MatrixModels_0.5-0  scales_1.1.1        mvtnorm_1.1-1      
    ##  [37] DBI_1.1.1           Rcpp_1.0.6          viridisLite_0.4.0  
    ##  [40] cmprsk_2.2-10       htmlTable_2.1.0     foreign_0.8-81     
    ##  [43] bit_4.0.4           lava_1.6.9          htmlwidgets_1.5.3  
    ##  [46] httr_1.4.2          RColorBrewer_1.1-2  ellipsis_0.3.2     
    ##  [49] pkgconfig_2.0.3     nnet_7.3-16         dbplyr_2.1.1       
    ##  [52] here_1.0.1          utf8_1.2.1          tidyselect_1.1.1   
    ##  [55] rlang_0.4.11        munsell_0.5.0       cellranger_1.1.0   
    ##  [58] tools_4.0.4         cachem_1.0.4        cli_2.5.0          
    ##  [61] generics_0.1.0      broom_0.7.6         evaluate_0.14      
    ##  [64] fastmap_1.1.0       yaml_2.2.1          bit64_4.0.5        
    ##  [67] fs_1.5.0            timereg_1.9.8       zip_2.1.1          
    ##  [70] future_1.21.0       nlme_3.1-152        quantreg_5.85      
    ##  [73] xml2_1.3.2          compiler_4.0.4      rstudioapi_0.13    
    ##  [76] curl_4.3.1          png_0.1-7           gt_0.3.0           
    ##  [79] reprex_2.0.0        broom.helpers_1.3.0 stringi_1.6.1      
    ##  [82] highr_0.9           Matrix_1.3-3        vctrs_0.3.8        
    ##  [85] pillar_1.6.0        lifecycle_1.0.0     furrr_0.2.2        
    ##  [88] data.table_1.14.0   conquer_1.0.2       R6_2.5.0           
    ##  [91] latticeExtra_0.6-29 KernSmooth_2.23-20  parallelly_1.25.0  
    ##  [94] codetools_0.2-18    polspline_1.1.19    MASS_7.3-54        
    ##  [97] assertthat_0.2.1    chron_2.3-56        rprojroot_2.0.2    
    ## [100] withr_2.4.2         multcomp_1.4-17     parallel_4.0.4     
    ## [103] hms_1.0.0           grid_4.0.4          rpart_4.1-15       
    ## [106] rmarkdown_2.8       numDeriv_2016.8-1.1 lubridate_1.7.10   
    ## [109] base64enc_0.1-3
