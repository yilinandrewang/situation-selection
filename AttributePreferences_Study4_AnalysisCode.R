# ***************** Summarized Attribute Preferences Have *****************
# ****************** Unique Antecedents and Consequences ******************
# Aline da Silva Frost, Y. Andre Wang, Paul W. Eastwick, and Alison Ledgerwood


# Study 4 Analysis Code
# Last update: 2022/06/19



# Prepare packages and functions ------------------------------------------

# Name the packages needed
packages <- c("lm.beta", "lavaan", "semTools", "apaTables", "multicon", "dplyr",
              "data.table")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == F)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = T))

# Create function that converts log odds ratios to correlation coefficients
# Formula from Borenstein, Hedges, Higgins, and Rothstein (2009)
lor2r <- function(x) {
  x*(sqrt(3) / pi)/
    sqrt((x*(sqrt(3)/pi))^2 + 4)
}



# Load data ---------------------------------------------------------------

# Load main data file
ss <- read.csv("AttributePreferences_Study4_data.csv")

# Calculate correlations between summarized preference (SP) and 
# functional preference (FP) for use below

# SP-FP correlation for intelligence
cor_sp <- cor(ss$spcomp_intel, ss$spcomp_conf)

# SP-FP correlation for confidence
cor_fp <- cor(ss$fpcomp_intel, ss$fpcomp_conf)

# Load stimuli data files for calculating reliabilities of FPs
mal <- read.csv("AttributePreferences_Study4_TargetRatingMale.csv")
fem <- read.csv("AttributePreferences_Study4_TargetRatingFemale.csv")
rom <- read.csv("AttributePreferences_Study4_RomanticRatings.csv")



# Correlations among observed variables -----------------------------------

# Reverse-code situation selection with experience (dv2_ci)
# so that the correlations are easier to interpret
ss$dv2_ci_r <- 1 - ss$dv2_ci

# Subset variables needed
key_var <- ss[, c(
  "spcomp_intel", # SP for intelligence
  "spcomp_conf",  # SP for confidence
  "fpcomp_intel", # FP for intelligence
  "fpcomp_conf",  # FP for confidence
  "dv1_intel",    # situation selection (SS) at a distance (intelligence)
  "dv1_conf",     # SS at a distance (confidence)
  "dv2_ci_r")]    # SS with experience (higher value = confidence)

# Create correlation tables
obs_cor_main <- cor(key_var)
round(obs_cor_main, 2)
apaTables::apa.cor.table(key_var)



# Reliability -------------------------------------------------------------


# *- SP -------------------------------------------------------------------

# Because we conducted structural equation modeling (SEM), reliabilities of
# SP are estimated from the measurement models using confirmatory factor
# analysis (CFA). 


# *--- Intelligence -------------------------------------------------------

# Specify CFA
cfa_sp_intel <- "
sp_intel_lat =~ sp_intelligent + sp_smart + sp_sharp
"

# Fit model
fit_sp_intel <- sem(model = cfa_sp_intel, data = ss, std.lv = T)

# Reliability
reliability(fit_sp_intel) # alpha = 0.89



# *--- Confidence ---------------------------------------------------------

# Specify CFA
cfa_sp_conf <- "
sp_conf_lat =~ a*sp_confident + a*sp_assured
"

# Fit model
fit_sp_conf <- sem(model = cfa_sp_conf, data = ss, std.lv = T)

# Reliability
reliability(fit_sp_conf) # alpha = 0.79



# *- FP -------------------------------------------------------------------

# Reliabilities of FP are estimated using the method discussed
# in Sherman & Wood (2014). We calculate them separately for male
# and female targets.

# Load romantic desire ratings for male targets
ratings_m <- rom[, 110:209]
set_ratings_m <- transpose(ratings_m)
set_ratings_m <- set_ratings_m[, colSums(
  is.na(set_ratings_m)) != nrow(set_ratings_m)]

# Load romantic desire ratings for female targets
ratings_f <- rom[, 10:109]
set_ratings_f <- transpose(ratings_f)
set_ratings_f <- set_ratings_f[, colSums(
  is.na(set_ratings_f)) != nrow(set_ratings_f)]



# *--- Intelligence ----------------------------------------------------

# Trait ratings of intelligence for male targets (norming data)
int_m <- pull(mal, 3)

# Vectorize intelligence ratings for male targets
vec_int_m <- (int_m - mean(int_m))/sd(int_m)

# Calculate alphaCP for intelligence ratings for male targets
alpha_mal_int <- vector.alpha(vec_int_m, set_ratings_m, type = "cor", 
                              CI = 0.95, CItype = "xci", minval = -1)
alpha_mal_int


# Trait ratings of intelligence for female targets(norming data)
int_f <- pull(fem, 3)

# Vectorize intelligence ratings for female targets
vec_int_f <- (int_f - mean(int_f))/sd(int_f)

# Calculate alphaCP for intelligence ratings for female targets
alpha_fem_int <- vector.alpha(vec_int_f, set_ratings_f, type = "cor", 
                              CI = 0.95, CItype = "xci", minval = -1)
alpha_fem_int



# *--- Confidence ---------------------------------------------------------

# Trait ratings of confidence for male targets (norming data)
con_m <- pull(mal, 2)

# Vectorize confidence ratings for male targets
vec_con_m <- (con_m - mean(con_m))/sd(con_m)

# Calculate alphaCP for confidence ratings for male targets
alpha_mal_con <- vector.alpha(vec_con_m, set_ratings_m, type = "cor", 
                              CI = 0.95, CItype = "xci", minval = -1)
alpha_mal_con


# Trait ratings of confidence for female targets (norming data)
con_f <- pull(fem, 2)

# Vectorize confidence ratings for female targets
vec_con_f <- (con_f - mean(con_f))/sd(con_f)

# Calculate alphaCP for confidence ratings for female targets
alpha_fem_con <- vector.alpha(vec_con_f, set_ratings_f, type = "cor", 
                              CI = 0.95, CItype = "xci", minval = -1)
alpha_fem_con



# Hypotheses 2 and 3: preregistered analyses, primary DVs -----------------


# *- SEM ------------------------------------------------------------------


# *--- SP predicting SS at a distance: intelligence -----------------------

# Specify model
mod_dv1_sp_intel <- "
sp_intel_lat =~ sp_intelligent + sp_smart + sp_sharp
sp_conf_lat =~ sp_confident + sp_assured

dv1_intel ~ sp_intel_lat + sp_conf_lat
"

# Fit model
fit_dv1_sp_intel <- sem(model = mod_dv1_sp_intel, data = ss, std.lv = T)

# Estimated model
summary(fit_dv1_sp_intel, fit.measures = T, standardized = T)

# Fit indices
round(fitMeasures(fit_dv1_sp_intel)[c("chisq", "df", "cfi", "tli", "rmsea")], 2)
round(fitMeasures(fit_dv1_sp_intel)["pvalue"], 3)
# Chi-squared(7) = 10.34, p = .170, CFI = 1.00, TLI = 1.00, RMSEA = 0.03

# Get parameter estimates
par_dv1_sp_intel <- parameterestimates(fit_dv1_sp_intel, standardized = T)
std_dv1_sp_intel <- standardizedsolution(fit_dv1_sp_intel)

# Estimate coefficients
round(par_dv1_sp_intel[6, c("est", "se", "std.all")], 2) # b, SE, beta
round(par_dv1_sp_intel[6, "pvalue"], 3) # p-value
round(par_dv1_sp_intel[6, "std.all"]*
        sqrt(1 - par_dv1_sp_intel[16, "std.all"]^2), 2) # r
round(std_dv1_sp_intel[6, "ci.lower"]*
        sqrt(1 - std_dv1_sp_intel[16, "est.std"]^2), 2) # r lower bound
round(std_dv1_sp_intel[6, "ci.upper"]*
        sqrt(1 - std_dv1_sp_intel[16, "est.std"]^2), 2) # r upper bound
# b = 0.90, SE = 0.13, p < .001, beta = 0.39, r = .32, 95% CI [.24, .41]



# *--- SP predicting SS at a distance: confidence -------------------------

# Specify model
mod_dv1_sp_conf <- "
sp_intel_lat =~ sp_intelligent + sp_smart + sp_sharp
sp_conf_lat =~ sp_confident + sp_assured

dv1_conf ~ sp_intel_lat + sp_conf_lat
"

# Fit model
fit_dv1_sp_conf <- sem(model = mod_dv1_sp_conf, data = ss, std.lv = T)

# Estimated model
summary(fit_dv1_sp_conf, fit.measures = T, standardized = T)

# Fit indices
round(fitMeasures(fit_dv1_sp_conf)[c("chisq", "df", "cfi", "tli", "rmsea")], 2)
round(fitMeasures(fit_dv1_sp_conf)["pvalue"], 3)
# Chi-squared(7) = 8.96, p = .256, CFI = 1.00, TLI = 1.00, RMSEA = 0.02

# Get parameter estimates
par_dv1_sp_conf <- parameterestimates(fit_dv1_sp_conf, standardized = T)
std_dv1_sp_conf <- standardizedsolution(fit_dv1_sp_conf)

# Estimate coefficients
round(par_dv1_sp_conf[7, c("est", "se", "std.all")], 2) # b, SE, beta
round(par_dv1_sp_conf[7, "pvalue"], 3) # p-value
round(par_dv1_sp_conf[7, "std.all"]*
        sqrt(1 - par_dv1_sp_conf[16, "std.all"]^2), 2) # r
round(std_dv1_sp_conf[7, "ci.lower"]*
        sqrt(1 - std_dv1_sp_conf[16, "est.std"]^2), 2) # r lower bound
round(std_dv1_sp_conf[7, "ci.upper"]*
        sqrt(1 - std_dv1_sp_conf[16, "est.std"]^2), 2) # r upper bound
# b = 0.86, SE = 0.14, p < .001, beta = 0.38, r = .31, 95% CI [.22, .40]



# *--- FP predicting SS at a distance: intelligence -----------------------

# Specify model
mod_dv1_fp_intel <- "
fp_intel_lat =~ fp_intel_p1 + fp_intel_p2 + fp_intel_p3 + fp_intel_p4
fp_conf_lat =~ fp_conf_p1 + fp_conf_p2 + fp_conf_p3 + fp_conf_p4

dv1_intel ~ fp_intel_lat + fp_conf_lat

fp_intel_p1 ~~ fp_conf_p1
fp_intel_p2 ~~ fp_conf_p2
fp_intel_p3 ~~ fp_conf_p3
fp_intel_p4 ~~ fp_conf_p4
"

# Fit model
fit_dv1_fp_intel <- sem(model = mod_dv1_fp_intel, data = ss, std.lv = T)

# Estimated model
summary(fit_dv1_fp_intel, fit.measures = T, standardized = T)

# Fit indices
round(fitMeasures(fit_dv1_fp_intel)[c("chisq", "df", "cfi", "tli", "rmsea")], 2)
round(fitMeasures(fit_dv1_fp_intel)["pvalue"], 3)
# Chi-squared(21) = 265.30, p < .001, CFI = 0.90, TLI = 0.84, RMSEA = 0.14

# Get parameter estimates
par_dv1_fp_intel <- parameterestimates(fit_dv1_fp_intel, standardized = T)
std_dv1_fp_intel <- standardizedsolution(fit_dv1_fp_intel)

# Estimate coefficients
round(par_dv1_fp_intel[9, c("est", "se", "std.all")], 2) # b, SE, beta
round(par_dv1_fp_intel[9, "pvalue"], 3) # p-value
round(par_dv1_fp_intel[9, "std.all"]*
        sqrt(1 - par_dv1_fp_intel[26, "std.all"]^2), 2) # Calculate r
round(std_dv1_fp_intel[9, "ci.lower"]*
        sqrt(1 - std_dv1_fp_intel[26, "est.std"]^2), 2) # r lower bound
round(std_dv1_fp_intel[9, "ci.upper"]*
        sqrt(1 - std_dv1_fp_intel[26, "est.std"]^2), 2) # r upper bound
# b = 0.40, SE = 0.18, p = .026, beta = 0.17, r = .11, 95% CI [.01, .21]



# *--- FP predicting SS at a distance: confidence -------------------------

# Specify model
mod_dv1_fp_conf <- "
fp_intel_lat =~ fp_intel_p1 + fp_intel_p2 + fp_intel_p3 + fp_intel_p4
fp_conf_lat =~ fp_conf_p1 + fp_conf_p2 + fp_conf_p3 + fp_conf_p4

dv1_conf ~ fp_intel_lat + fp_conf_lat

fp_intel_p1 ~~ fp_conf_p1
fp_intel_p2 ~~ fp_conf_p2
fp_intel_p3 ~~ fp_conf_p3
fp_intel_p4 ~~ fp_conf_p4
"

# Fit model
fit_dv1_fp_conf <- sem(model = mod_dv1_fp_conf, data = ss, std.lv = T)

# Estimated model
summary(fit_dv1_fp_conf, fit.measures = T, standardized = T)

# Fit indices
round(fitMeasures(fit_dv1_fp_conf)[c("chisq", "df", "cfi", "tli", "rmsea")], 2)
round(fitMeasures(fit_dv1_fp_conf)["pvalue"], 3)
# Chi-squared(21) = 263.68, p < .001, CFI = 0.90, TLI = 0.84, RMSEA = 0.14

# Get parameter estimates
par_dv1_fp_conf <- parameterestimates(fit_dv1_fp_conf, standardized = T)
std_dv1_fp_conf <- standardizedsolution(fit_dv1_fp_conf)

# Estimate coefficient
round(par_dv1_fp_conf[10, c("est", "se", "std.all")], 2) # b, SE, beta
round(par_dv1_fp_conf[10, "pvalue"], 3) # p-value
round(par_dv1_fp_conf[10, "std.all"]*
        sqrt(1 - par_dv1_fp_conf[26, "std.all"]^2), 2) # Calculate r
round(std_dv1_fp_conf[10, "ci.lower"]*
        sqrt(1 - std_dv1_fp_conf[26, "est.std"]^2), 2) # r lower bound
round(std_dv1_fp_conf[10, "ci.upper"]*
        sqrt(1 - std_dv1_fp_conf[26, "est.std"]^2), 2) # r upper bound
# b = 0.24, SE = 0.18, p = .174, beta = 0.10, r = .07, 95% CI [-.03, .17]



# *--- SP predicting SS with experience: intelligence & confidence --------

# Specify model
mod_dv2_sp <- "
sp_intel_lat =~ sp_intelligent + sp_smart + sp_sharp
sp_conf_lat =~ sp_confident + sp_assured

dv2_ci ~ sp_intel_lat + sp_conf_lat
"

# Fit model
fit_dv2_sp <- sem(model = mod_dv2_sp, ordered = "dv2_ci", data = ss, std.lv = T)

# Estimated model
summary(fit_dv2_sp, fit.measures = T, standardized = T)

# Fit indices
round(fitMeasures(fit_dv2_sp)[c("chisq", "df", "cfi", "tli", "rmsea")], 2)
round(fitMeasures(fit_dv2_sp)["pvalue"], 3)
# Chi-squared(7) = 1.23, p = .990, CFI = 1.00, TLI = 1.01, RMSEA = 0.00

# Get parameter estimates
par_dv2_sp <- parameterestimates(fit_dv2_sp)

# Estimate coefficients: Intelligence
round(par_dv2_sp[6, c("est", "se")]*1.7, 2) # b, SE (logit = probit*1.7)
par_dv2_sp[6, "pvalue"] # p-value
round(exp(par_dv2_sp[6, "est"]*1.7), 2) # odds ratio
round(lor2r(par_dv2_sp[6, "est"]*1.7), 2) # r
round(lor2r(par_dv2_sp[6, "ci.lower"]*1.7), 2) # r lower
round(lor2r(par_dv2_sp[6, "ci.upper"]*1.7), 2) # r upper
# b = 0.44, SE = 0.13, p < .001, OR = 1.55, r = .12, 95% CI [.05, .19]

# Estimate coefficients: Confidence
# Note that because choice of "confident" website was coded as 0,
# the signs of all coefficients need to be flipped
round(par_dv2_sp[7, "est"]*1.7*(-1), 2) # b (logit = probit*1.7)
round(par_dv2_sp[7, "se"]*1.7, 2) # SE (logit = probit*1.7)
round(par_dv2_sp[7, "pvalue"], 3) # p-value
round(exp(par_dv2_sp[7, "est"]*1.7*(-1)), 2) # odds ratio
round(lor2r(par_dv2_sp[7, "est"]*1.7*(-1)), 2) # r
round(lor2r(par_dv2_sp[7, "ci.upper"]*1.7*(-1)), 2) # r lower (b/c sign flip)
round(lor2r(par_dv2_sp[7, "ci.lower"]*1.7*(-1)), 2) # r upper (b/c sign flip)
# b = 0.54, SE = 0.13, p < .001, OR = 1.71, r = .15, 95% CI [.08, .21]



# *--- FP predicting SS with experience: intelligence & confidence --------

# Specify model
mod_dv2_fp <- "
fp_intel_lat =~ fp_intel_p1 + fp_intel_p2 + fp_intel_p3 + fp_intel_p4
fp_conf_lat =~ fp_conf_p1 + fp_conf_p2 + fp_conf_p3 + fp_conf_p4

dv2_ci ~ fp_intel_lat + fp_conf_lat

fp_intel_p1 ~~ fp_conf_p1
fp_intel_p2 ~~ fp_conf_p2
fp_intel_p3 ~~ fp_conf_p3
fp_intel_p4 ~~ fp_conf_p4
"

# Fit model
fit_dv2_fp <- sem(model = mod_dv2_fp, ordered = "dv2_ci", data = ss, std.lv = T)

# Estimated model
summary(fit_dv2_fp, fit.measures = T, standardized = T)

# Fit indices
round(fitMeasures(fit_dv2_fp)[c("chisq", "df", "cfi", "tli", "rmsea")], 2)
round(fitMeasures(fit_dv2_fp)["pvalue"], 3)
# Chi-squared(21) = 125.86, p < .001, CFI = 0.96, TLI = 0.94, RMSEA = 0.09

# Get parameter estimates
par_dv2_fp <- parameterestimates(fit_dv2_fp)

# Estimate coefficient: Intelligence
round(par_dv2_fp[9, c("est", "se")]*1.7, 2) # b, SE (logit = probit*1.7)
round(par_dv2_fp[9, "pvalue"], 3) # p-value
round(exp(par_dv2_fp[9, "est"]*1.7), 2) # odds ratio
round(lor2r(par_dv2_fp[9, "est"]*1.7), 2) # r
round(lor2r(par_dv2_fp[9, "ci.lower"]*1.7), 2) # r lower
round(lor2r(par_dv2_fp[9, "ci.upper"]*1.7), 2) # r upper
# b = 1.06, SE = 0.13, p < .001, OR = 2.87, r = .28, 95% CI [.22, .34]

# Estimate coefficient: Confidence
# Note that because choice of "confident" website was coded as 0,
# the signs of all coefficients need to be flipped
round(par_dv2_fp[10, "est"]*1.7*(-1), 2) # b (logit = probit*1.7)
round(par_dv2_fp[10, "se"]*1.7, 2) # SE (logit = probit*1.7)
round(par_dv2_fp[10, "pvalue"], 3) # p-value
round(exp(par_dv2_fp[10, "est"]*1.7*(-1)), 2) # odds ratio
round(lor2r(par_dv2_fp[10, "est"]*1.7*(-1)), 2) # r
round(lor2r(par_dv2_fp[10, "ci.upper"]*1.7*(-1)), 2) # r lower (b/c sign flip)
round(lor2r(par_dv2_fp[10, "ci.lower"]*1.7*(-1)), 2) # r upper (b/c sign flip)
# b = 1.44, SE = 0.13, p < .001, OR = 4.21, r = .37, 95% CI [.31, .42]



# *- Bivariate Regression -------------------------------------------------


# *--- SP predicting SS at a distance: intelligence -----------------------

# Fit model
br_dv1_sp_intel <- lm.beta(lm(dv1_intel ~ spcomp_intel, data = ss))

# Estimate coefficients
round(summary(br_dv1_sp_intel)$
        coefficients["spcomp_intel", 1:3], 2) # b, beta ( = r), SE
round(summary(br_dv1_sp_intel)$coefficients["spcomp_intel", 5], 3) # p-value
round(confint(br_dv1_sp_intel), 2)
# b = 0.73, beta = 0.33, SE = 0.09, p < .001, 95% CI [.15, .50]



# *--- SP predicting SS at a distance: confidence -------------------------

# Fit model
br_dv1_sp_conf <- lm.beta(lm(dv1_conf ~ spcomp_conf, data = ss))

# Estimate coefficients
round(summary(br_dv1_sp_conf)$
        coefficients["spcomp_conf", 1:3], 2) # b, beta ( = r), SE
round(summary(br_dv1_sp_conf)$coefficients["spcomp_conf", 5], 3) # p-value
round(confint(br_dv1_sp_conf), 2)
# b = 0.59, beta = 0.29, SE = 0.08, p < .001, 95% CI [.13, .45]



# *--- FP predicting SS at a distance: intelligence -----------------------

# Fit model
br_dv1_fp_intel <- lm.beta(lm(dv1_intel ~ fpcomp_intel, data = ss))

# Estimate coefficients
round(summary(br_dv1_fp_intel)$
        coefficients["fpcomp_intel", 1:3], 2) # b, beta ( = r), SE
round(summary(br_dv1_fp_intel)$coefficients["fpcomp_intel", 5], 3) # p-value
cor.test(ss$dv1_intel, ss$fpcomp_intel) # 95% CI of correlation
# Note: confint(br_dv1_fp_intel) gives inprobable values
# b = 1.06, beta = 0.09, SE = 0.49, p = .031, 95% CI [.008, .17]



# *--- FP predicting SS at a distance: confidence -----------------------

# Fit model
br_dv1_fp_conf <- lm.beta(lm(dv1_conf ~ fpcomp_conf, data = ss))

# Estimate coefficients
round(summary(br_dv1_fp_conf)$
        coefficients["fpcomp_conf", 1:3], 2) # b, beta ( = r), SE
round(summary(br_dv1_fp_conf)$coefficients["fpcomp_conf", 5], 3) # p-value
cor.test(ss$dv1_conf, ss$fpcomp_conf) # 95% CI of correlation
# b = 0.39, beta = 0.04, SE = 0.43, p = .366, 95% CI [-0.04, .12]



# *--- SP predicting SS with experience: intelligence ---------------------

# Fit model
br_dv2_sp_intel <- glm(dv2_ci ~ spcomp_intel, data = ss, family = "binomial")

# Estimate coefficients
round(summary(br_dv2_sp_intel)$coefficients["spcomp_intel", 1:2], 2) # b, SE
round(summary(br_dv2_sp_intel)$coefficients["spcomp_intel", 4], 3) # p-value
round(exp(coef(br_dv2_sp_intel)["spcomp_intel"]), 2) # odds ratio
round(lor2r(coef(br_dv2_sp_intel)["spcomp_intel"]), 2) # r
round(lor2r(confint(br_dv2_sp_intel)["spcomp_intel", ]), 2) # 95% CI
# b = 0.13, SE = 0.10, p = .192, OR = 1.13, r = .03, 95% CI [-.02, .09]



# *--- SP predicting SS with experience: confidence -----------------------

# Fit model
br_dv2_sp_conf <- glm(dv2_ci ~ spcomp_conf, data = ss, family = "binomial")

# Estimate coefficients
# Note that because choice of "confident" website was coded as 0,
# the signs of all coefficients need to be flipped
round(summary(br_dv2_sp_conf)$coefficients["spcomp_conf", 1]*(-1), 2) # b
round(summary(br_dv2_sp_conf)$coefficients["spcomp_conf", 2], 2) # SE
round(summary(br_dv2_sp_conf)$coefficients["spcomp_conf", 4], 3) # p-value
round(exp(coef(br_dv2_sp_conf)["spcomp_conf"]*(-1)), 2) # odds ratio
round(lor2r(coef(br_dv2_sp_conf)["spcomp_conf"]*(-1)), 2) # r
round(lor2r(confint(br_dv2_sp_conf)["spcomp_conf", 2])*(-1), 2) # r lower bound
round(lor2r(confint(br_dv2_sp_conf)["spcomp_conf", 1])*(-1), 2) # r upper bound
# b = 0.23, SE = 0.08, p = .006, OR = 1.25, r = .06, 95% CI [.02, .11]



# *--- FP predicting SS with experience: intelligence ---------------------

# Fit model
br_dv2_fp_intel <- glm(dv2_ci ~ fpcomp_intel, data = ss, family = "binomial")

# Estimate coefficients
round(summary(br_dv2_fp_intel)$coefficients["fpcomp_intel", 1:2], 2) # b, SE
round(summary(br_dv2_fp_intel)$coefficients["fpcomp_intel", 4], 3) # p-value
round(exp(coef(br_dv2_fp_intel)["fpcomp_intel"]), 2) # odds ratio
round(lor2r(coef(br_dv2_fp_intel)["fpcomp_intel"]), 2) # r
round(lor2r(confint(br_dv2_fp_intel)), 2)["fpcomp_intel", ] # 95% CI
# b = 0.46, SE = 0.48, p = .339, OR = 1.59, r = .13, 95% CI [-.13, .36]



# *--- FP predicting SS with experience: confidence -----------------------

# Fit model
br_dv2_fp_conf <- glm(dv2_ci ~ fpcomp_conf, data = ss, family = "binomial")

# Estimate coefficients
# Note that because choice of "confident" website was coded as 0,
# the signs of all coefficients need to be flipped
round(summary(br_dv2_fp_conf)$coefficients["fpcomp_conf", 1]*(-1), 2) # b
round(summary(br_dv2_fp_conf)$coefficients["fpcomp_conf", 2], 2) # SE
round(summary(br_dv2_fp_conf)$coefficients["fpcomp_conf", 4], 3) # p-value
round(exp(coef(br_dv2_fp_conf)["fpcomp_conf"]*(-1)), 2) # odds ratio
round(lor2r(coef(br_dv2_fp_conf)["fpcomp_conf"]*(-1)), 2) # r
round(lor2r(confint(br_dv2_fp_conf)["fpcomp_conf", 2])*(-1), 2) # r lower bound
round(lor2r(confint(br_dv2_fp_conf)["fpcomp_conf", 1])*(-1), 2) # r upper bound
# b = 2.95, SE = 0.46, p < .001, OR = 19.08, r = .63, 95% CI [.50, .73]



# *- Multiple Regression --------------------------------------------------


# *--- SP predicting SS at a distance: intelligence -----------------------

# Fit model
mr_dv1_sp_intel <- lm.beta(lm(dv1_intel ~ spcomp_intel + spcomp_conf, 
                              data = ss))

# Estimate coefficients
round(summary(mr_dv1_sp_intel)$
        coefficients["spcomp_intel", 1:3], 2) # b, beta, SE
round(summary(mr_dv1_sp_intel)$coefficients["spcomp_intel", 5], 3) # p-value
round(summary(mr_dv1_sp_intel)$coefficients["spcomp_intel", 2]*
        sqrt(1 - cor_sp^2), 2) # r
round(confint(mr_dv1_sp_intel)*sqrt(1 - cor_sp^2), 2)["spcomp_intel", ]
# b = 0.75, SE = 0.10, beta = 0.34, p < .001, r = .30, 95% CI [.12, .47]



# *--- SP predicting SS at a distance: confidence ------------------------

# Fit model
mr_dv1_sp_conf <- lm.beta(lm(dv1_conf ~ spcomp_conf + spcomp_intel, 
                             data = ss))

# Estimate coefficients
round(summary(mr_dv1_sp_conf)$
        coefficients["spcomp_conf", 1:3], 2) # b, beta, SE
round(summary(mr_dv1_sp_conf)$coefficients["spcomp_conf", 5], 3) # p-value
round(summary(mr_dv1_sp_conf)$coefficients["spcomp_conf", 2]*
        sqrt(1 - cor_sp^2), 2) # r
round(confint(mr_dv1_sp_conf)*sqrt(1 - cor_sp^2), 2)["spcomp_conf", ]
# b = 0.62, SE = 0.09, beta = 0.31, p < .001, r = .27, 95% CI [.11, .43]



# *--- FP predicting SS at a distance: intelligence -----------------------

# Fit model
mr_dv1_fp_intel <- lm.beta(lm(dv1_intel ~ fpcomp_intel + fpcomp_conf, 
                              data = ss))

# Estimate coefficients
round(summary(mr_dv1_fp_intel)$
        coefficients["fpcomp_intel", 1:3], 2) # b, beta, SE
round(summary(mr_dv1_fp_intel)$coefficients["fpcomp_intel", 5], 3) # p-value
round(summary(mr_dv1_fp_intel)$coefficients["fpcomp_intel", 2]*
        sqrt(1- cor_fp^2), 2) # r

# Manually standardize each variable to get accurate CIs
# Note that this model yields identical betas and p-values as the model above
ss$dv1_intel_std <- scale(ss$dv1_intel)
ss$fpcomp_intel_std <- scale(ss$fpcomp_intel)
ss$fpcomp_conf_std <- scale(ss$fpcomp_conf)
mr_dv1_fp_intel_std <- lm.beta(lm(
  dv1_intel_std ~ fpcomp_intel_std + fpcomp_conf_std, data = ss))
confint(mr_dv1_fp_intel_std)["fpcomp_intel_std", ]*sqrt(1 - cor_fp^2)
# b = 1.36, SE = 0.66, beta = 0.12, p = .039, r = .09, 95% CI [.004, .17]



# *--- FP predicting SS at a distance: confidence -------------------------

# Fit model
mr_dv1_fp_conf <- lm.beta(lm(dv1_conf ~ fpcomp_conf + fpcomp_intel, 
                             data = ss))

# Estimate coefficients
round(summary(mr_dv1_fp_conf)$
        coefficients["fpcomp_conf", 1:3], 2) # b, beta, SE
round(summary(mr_dv1_fp_conf)$coefficients["fpcomp_conf", 5], 3) # p-value
round(summary(mr_dv1_fp_conf)$coefficients["fpcomp_conf", 2]*
        sqrt(1- cor_fp^2), 2) # r

# Manually standardize each variable to get accurate CIs
# Note that this model yields identical betas and p-values as the model above
ss$dv1_conf_std <- scale(ss$dv1_conf)
mr_dv1_fp_conf_std <- lm.beta(lm(
  dv1_conf_std ~ fpcomp_intel_std + fpcomp_conf_std, data = ss))
round(confint(mr_dv1_fp_conf_std)["fpcomp_conf_std", ]*sqrt(1 - cor_fp^2), 2)
# b = 0.74, SE = 0.57, beta = 0.07, p = .193, r = .06, 95% CI [-.03, .14]



# *--- SP predicting SS with experience: intelligence & confidence --------

# Fit model
mr_dv2_sp <- glm(dv2_ci ~ spcomp_intel + spcomp_conf, 
                 data = ss, family = "binomial")

# Estimate coefficients: Intelligence
round(summary(mr_dv2_sp)$coefficients["spcomp_intel", 1:2], 2) # b, SE
round(summary(mr_dv2_sp)$coefficients["spcomp_intel", 4], 3) # p-value
round(exp(coef(mr_dv2_sp)["spcomp_intel"]), 2) # odds ratio
round(lor2r(coef(mr_dv2_sp)["spcomp_intel"]), 2) # r
round(lor2r(confint(mr_dv2_sp)), 2)["spcomp_intel", ]
# b = 0.32, SE = 0.11, p = .004, OR = 1.38, r = .09, 95% CI [.03, .15]

# Estimate coefficients: Confidence
# Note that because choice of "confident" website was coded as 0,
# the signs of all coefficients need to be flipped
round(summary(mr_dv2_sp)$coefficients["spcomp_conf", 1]*(-1), 2) # b
round(summary(mr_dv2_sp)$coefficients["spcomp_conf", 2], 2) # SE
round(summary(mr_dv2_sp)$coefficients["spcomp_conf", 4], 3) # p-value
round(exp(coef(mr_dv2_sp)["spcomp_conf"]*(-1)), 2) # odds ratio
round(lor2r(coef(mr_dv2_sp)["spcomp_conf"]*(-1)), 2) # r 
round(lor2r(confint(mr_dv2_sp)["spcomp_conf", 2]*(-1)), 2) # r lower bound
round(lor2r(confint(mr_dv2_sp)["spcomp_conf", 1]*(-1)), 2) # r upper bound
# b = 0.36, SE = 0.10, p < 0.001, OR = 1.44, r = .10, 95% CI [.05, .15]



# *--- FP predicting SS with experience: intelligence & confidence --------

# Fit model
mr_dv2_fp <- glm(dv2_ci ~ fpcomp_intel + fpcomp_conf, 
                 data = ss, family = "binomial")

# Estimate coefficients: Intelligence
round(summary(mr_dv2_fp)$coefficients["fpcomp_intel", 1:2], 2) # b, SE
round(summary(mr_dv2_fp)$coefficients["fpcomp_intel", 4], 3) # p-value
round(exp(coef(mr_dv2_fp)["fpcomp_intel"]), 2) # odds ratio
round(lor2r(coef(mr_dv2_fp)["fpcomp_intel"]), 2) # r
round(lor2r(confint(mr_dv2_fp)), 2)
# b = 5.48, SE = 0.78, p < 0.001, OR = 239.86, r = .83, 95% CI [.74, .89]

# Estimate coefficients: Confidence
# Note that because choice of "confident" website was coded as 0,
# the signs of all coefficients need to be flipped
round(summary(mr_dv2_fp)$coefficients["fpcomp_conf", 1]*(-1), 2) # b
round(summary(mr_dv2_fp)$coefficients["fpcomp_conf", 2], 2) # SE
round(summary(mr_dv2_fp)$coefficients["fpcomp_conf", 4], 3) # p-value
round(exp(coef(mr_dv2_fp)["fpcomp_conf"]*(-1)), 2) # odds ratio
round(lor2r(coef(mr_dv2_fp)["fpcomp_conf"]*(-1)), 2) # r
round(lor2r(confint(mr_dv2_fp)["fpcomp_conf", 2]*(-1)), 2) # r lower bound
round(lor2r(confint(mr_dv2_fp)["fpcomp_conf", 1]*(-1)), 2) # r upper bound
# b = 6.41, SE = 0.72, p < 0.001, OR = 605.52, r = .87, 95% CI [.81, .91]



# Hypotheses 2 and 3: exploring a full model of double dissociation -------


# *- Fit model ------------------------------------------------------------

mod_SEM_all <- "
sp_intel_lat =~ sp_intelligent + sp_smart + sp_sharp
sp_conf_lat =~ a*sp_confident + a*sp_assured

fp_intel_lat =~ fp_intel_p1 + fp_intel_p2 + fp_intel_p3 + fp_intel_p4
fp_conf_lat =~ fp_conf_p1 + fp_conf_p2 + fp_conf_p3 + fp_conf_p4

fp_intel_p1 ~~ fp_conf_p1
fp_intel_p2 ~~ fp_conf_p2
fp_intel_p3 ~~ fp_conf_p3
fp_intel_p4 ~~ fp_conf_p4

dv1_intel ~ sp_intel_lat + fp_intel_lat + sp_conf_lat + fp_conf_lat
dv1_conf ~ sp_intel_lat + fp_intel_lat + sp_conf_lat + fp_conf_lat
dv2_ci ~ sp_intel_lat + fp_intel_lat + sp_conf_lat + fp_conf_lat
"

# Fit model
fit_SEM_all <- sem(model = mod_SEM_all, ordered = "dv2_ci", 
                   data = ss, std.lv = T)

# Estimated model
summary(fit_SEM_all, fit.measures = T, standardized = T, n = 4)

# Fit indices
round(fitMeasures(fit_SEM_all)[c("chisq", "df", "cfi", "tli", "rmsea")], 2)
round(fitMeasures(fit_SEM_all)["pvalue"], 3)
# Chi-squared(83) = 178.16, p < .001, CFI = 0.98, TLI = 0.97, RMSEA = 0.05



# *- Parameter estimates ------------------------------------------------

par_SEM_all <- parameterestimates(fit_SEM_all, standardized = T)
std_SEM_all <- standardizedsolution(fit_SEM_all)

# Get logit coefficients (note: only read cofficients from experience DV)
par_SEM_all_logit <- par_SEM_all[26:29, ] # probit
par_SEM_all_logit$est <- round(par_SEM_all_logit$est * 1.7, 2) # convert to logit
par_SEM_all_logit$se <- round(par_SEM_all_logit$se * 1.7, 2) # convert to logit
par_SEM_all_logit

# SP predicting SS at a distance: intelligence
par_SEM_all[18, 1:3] # Parameter 18
round(par_SEM_all[18, c("std.all")], 2)
round(par_SEM_all[18, "pvalue"], 3)
round(std_SEM_all[18, "ci.lower"], 2); round(std_SEM_all[18, "ci.upper"], 2)
# beta = 0.37, p < .001, 95% CI [0.28, 0.46]

# SP predicting SS at a distance: confidence
par_SEM_all[24, 1:3] # Parameter 24
round(par_SEM_all[24, c("std.all")], 2)
round(par_SEM_all[24, "pvalue"], 3)
round(std_SEM_all[24, "ci.lower"], 2); round(std_SEM_all[24, "ci.upper"], 2)
# beta = 0.39, p < .001, 95% CI [0.29, 0.49]

# FP predicting SS at a distance: intelligence
par_SEM_all[19, 1:3] # Parameter 19
round(par_SEM_all[19, c("std.all")], 2)
round(par_SEM_all[19, "pvalue"], 3)
round(std_SEM_all[19, "ci.lower"], 2); round(std_SEM_all[19, "ci.upper"], 2)
# beta = 0.06, p = .323, 95% CI [-0.06, 0.19]

# FP predicting SS at a distance: confidence
par_SEM_all[25, 1:3] # Parameter 25
round(par_SEM_all[25, c("std.all")], 2)
round(par_SEM_all[25, "pvalue"], 3)
round(std_SEM_all[25, "ci.lower"], 2); round(std_SEM_all[25, "ci.upper"], 2)
# beta = 0.02, p = .769, 95% CI [-0.11, 0.14]

# SP predicting SS with experience: intelligence
par_SEM_all[26, 1:3] # Parameter 26
par_SEM_all[26, "std.all"] # probit (est = std.lv = std.all)
par_SEM_all[26, "std.all"]*1.7 # convert to logit
round(exp(par_SEM_all[26, "std.all"]*1.7), 2) # OR
round(par_SEM_all[26, "pvalue"], 3) # p-value
round(exp(par_SEM_all[26, "ci.lower"]*1.7), 2) # OR lower bound
round(exp(par_SEM_all[26, "ci.upper"]*1.7), 2) # OR upper bound
# OR = 1.45, p = .006, 95% CI = [1.11, 1.90]

# SP predicting SS with experience: confidence (signs need to be flipped)
par_SEM_all[28, 1:3] # Parameter 28
par_SEM_all[28, "std.all"]*(-1) # probit coeffient (est = std.lv = std.all)
par_SEM_all[28, "std.all"]*1.7*(-1) # convert to logit
round(exp(par_SEM_all[28, "std.all"]*1.7*(-1)), 2) # OR
round(par_SEM_all[28, "pvalue"], 3) # p-value
round(exp(par_SEM_all[28, "ci.upper"]*1.7*(-1)), 2) # OR upper bound
round(exp(par_SEM_all[28, "ci.lower"]*1.7*(-1)), 2) # OR lower bound
# OR = 1.42, p = .011, 95% CI = [1.08, 1.85]

# FP predicting SS with experience: intelligence
par_SEM_all[27, 1:3] # Parameter 27
par_SEM_all[27, "std.all"] # probit coeffient (est = std.lv = std.all)
par_SEM_all[27, "std.all"]*1.7 # convert to logit
round(exp(par_SEM_all[27, "std.all"]*1.7), 2) # OR
round(par_SEM_all[27, "pvalue"], 3) # p-value
round(exp(par_SEM_all[27, "ci.lower"]*1.7), 2) # OR lower bound
round(exp(par_SEM_all[27, "ci.upper"]*1.7), 2) # OR upper bound
# OR = 2.51, p < .001, 95% CI = [1.91, 3.32]

# FP predicting SS with experience: confidence
par_SEM_all[29, 1:3] # Parameter 29
par_SEM_all[29, "std.all"]*(-1) # probit coeffient (est = std.lv = std.all)
par_SEM_all[29, "std.all"]*1.7*(-1) # convert to logit
round(exp(par_SEM_all[29, "std.all"]*1.7*(-1)), 2) # OR
round(par_SEM_all[29, "pvalue"], 3) # p-value
round(exp(par_SEM_all[29, "ci.upper"]*1.7*(-1)), 2) # OR upper bound
round(exp(par_SEM_all[29, "ci.lower"]*1.7*(-1)), 2) # OR lower bound
# OR = 3.94, p < .001, 95% CI = [3.05, 5.08]



# *- Testing significance of difference -----------------------------------

# Specify measurement model to be used below
mod_SEM_measure <- "
sp_intel_lat =~ sp_intelligent + sp_smart + sp_sharp
sp_conf_lat =~ a*sp_confident + a*sp_assured

fp_intel_lat =~ fp_intel_p1 + fp_intel_p2 + fp_intel_p3 + fp_intel_p4
fp_conf_lat =~ fp_conf_p1 + fp_conf_p2 + fp_conf_p3 + fp_conf_p4

fp_intel_p1 ~~ fp_conf_p1
fp_intel_p2 ~~ fp_conf_p2
fp_intel_p3 ~~ fp_conf_p3
fp_intel_p4 ~~ fp_conf_p4
"



# *--- SP-FP intelligence on distance -------------------------------------

# Specify model with equality constraint
mod_SPFPintel_dv1 <- "
dv1_intel ~ b*sp_intel_lat + b*fp_intel_lat + sp_conf_lat + fp_conf_lat
dv1_conf ~ sp_intel_lat + fp_intel_lat + sp_conf_lat + fp_conf_lat
dv2_ci ~ sp_intel_lat + fp_intel_lat + sp_conf_lat + fp_conf_lat
"
mod_SPFPintel_dv1_all <- paste(mod_SEM_measure, mod_SPFPintel_dv1)

# Fit constrained model
fit_SPFPintel_dv1_all <- sem(model = mod_SPFPintel_dv1_all, ordered = "dv2_ci", 
                             data = ss, std.lv = T)

# Estimate constrained model
summary(fit_SPFPintel_dv1_all, fit.measures = T, standardized = T)

# Compare constrained model to unconstrained model
anova(fit_SEM_all, fit_SPFPintel_dv1_all) # chi-squared(1) = 10.26, p = .001



# *--- SP-FP confidence on distance ---------------------------------------

# Specify model with equality constraint
mod_SPFPconf_dv1 <- "
dv1_intel ~ sp_intel_lat + fp_intel_lat + sp_conf_lat + fp_conf_lat
dv1_conf ~ sp_intel_lat + fp_intel_lat + b*sp_conf_lat + b*fp_conf_lat
dv2_ci ~ sp_intel_lat + fp_intel_lat + sp_conf_lat + fp_conf_lat
"
mod_SPFPconf_dv1_all <- paste(mod_SEM_measure, mod_SPFPconf_dv1)

# Fit constrained model
fit_SPFPconf_dv1_all <- sem(model = mod_SPFPconf_dv1_all, ordered = "dv2_ci", 
                            data = ss, std.lv = T)

# Estimate constrained model
summary(fit_SPFPconf_dv1_all, fit.measures = T, standardized = T)

# Compare constrained model to unconstrained model
anova(fit_SEM_all, fit_SPFPconf_dv1_all) # chi-squared(1) = 17.77, p < .001



# *--- SP-FP intelligence on experience -----------------------------------

# Specify model with equality constraint
mod_SPFPintel_dv2 <- "
dv1_intel ~ sp_intel_lat + fp_intel_lat + sp_conf_lat + fp_conf_lat
dv1_conf ~ sp_intel_lat + fp_intel_lat + sp_conf_lat + fp_conf_lat
dv2_ci ~ b*sp_intel_lat + b*fp_intel_lat + sp_conf_lat + fp_conf_lat
"
mod_SPFPintel_dv2_all <- paste(mod_SEM_measure, mod_SPFPintel_dv2)

# Fit constrained model
fit_SPFPintel_dv2_all <- sem(model = mod_SPFPintel_dv2_all, ordered = "dv2_ci", 
                             data = ss, std.lv = T)

# Estimate constrained model
summary(fit_SPFPintel_dv2_all, fit.measures = T, standardized = T)

# Compare constrained model to unconstrained model
anova(fit_SEM_all, fit_SPFPintel_dv2_all) # chi-squared(1) = 5.28, p = .022



# *--- SP-FP confidence on experience -----------------------------------

# Specify model with equality constraint
mod_SPFPconf_dv2 <- "
dv1_intel ~ sp_intel_lat + fp_intel_lat + sp_conf_lat + fp_conf_lat
dv1_conf ~ sp_intel_lat + fp_intel_lat + sp_conf_lat + fp_conf_lat
dv2_ci ~ sp_intel_lat + fp_intel_lat + b*sp_conf_lat + b*fp_conf_lat
"

mod_SPFPconf_dv2_all <- paste(mod_SEM_measure, mod_SPFPconf_dv2)

# Fit constrained model
fit_SPFPconf_dv2_all <- sem(model = mod_SPFPconf_dv2_all, ordered = "dv2_ci", 
                            data = ss, std.lv = T)

# Estimate constrained model
summary(fit_SPFPconf_dv2_all, fit.measures = T, standardized = T)

# Compare constrained model to unconstrained model
anova(fit_SEM_all, fit_SPFPconf_dv2_all) # chi-squared(1) = 20.46, p < .001



# Hypotheses 2 and 3: preregistered analyses, secondary DVs ---------------


# *- SEM ------------------------------------------------------------------


# *--- SP predicting SS at a distance (tradeoff): intelligence & confidence ----

# Specify model
mod_dv3_sp <- "
sp_intel_lat =~ sp_intelligent + sp_smart + sp_sharp
sp_conf_lat =~ sp_confident + sp_assured

dv3_ic ~ sp_intel_lat + sp_conf_lat
"

# Fit model
fit_dv3_sp <- sem(model = mod_dv3_sp, data = ss, std.lv = T)

# Estimated model
summary(fit_dv3_sp, fit.measures = T, standardized = T)

# Fit indices
round(fitMeasures(fit_dv3_sp)[c("chisq", "df", "cfi", "tli", "rmsea")], 2)
round(fitMeasures(fit_dv3_sp)["pvalue"], 3)
# Chi-squared(7) = 7.81, p = .350, CFI = 1.00, TLI = 1.00, RMSEA = 0.01

# Get parameter estimates
par_dv3_sp <- parameterestimates(fit_dv3_sp, standardized = T)
std_dv3_sp <- standardizedsolution(fit_dv3_sp)

# Estimate coefficients: Intelligence
# Note that lower value of DV = greater preference of "intelligent" website
# so the signs of all coefficients need to be flipped
round(par_dv3_sp[6, c("est", "std.all")]*(-1), 2) # b, beta
round(par_dv3_sp[6, "se"], 2) # SE
round(par_dv3_sp[6, "pvalue"], 3) # p-value
round(par_dv3_sp[6, "std.all"]*(-1)*
        sqrt(1 - (par_dv3_sp[16, "std.all"]*(-1))^2), 2) # r
round(std_dv3_sp[6, "ci.upper"]*
        sqrt(1 - std_dv3_sp[16, "est.std"]^2)*(-1), 2) # r lower bound
round(std_dv3_sp[6, "ci.lower"]*
        sqrt(1 - std_dv3_sp[16, "est.std"]^2)*(-1), 2) # r upper bound
# b = 1.27, SE = 0.12, p < .001, beta = 0.58, r = .48, 95% CI [.40, .57]

# Estimate coefficients: Confidence
round(par_dv3_sp[7, c("est", "se", "std.all")], 2) # b, SE, beta
round(par_dv3_sp[7, "pvalue"], 3) # p-value
round(par_dv3_sp[7, "std.all"]*
        sqrt(1 - par_dv3_sp[16, "std.all"]^2), 2) # r
round(std_dv3_sp[7, "ci.lower"]*
        sqrt(1 - std_dv3_sp[16, "est.std"]^2), 2) # r lower bound
round(std_dv3_sp[7, "ci.upper"]*
        sqrt(1 - std_dv3_sp[16, "est.std"]^2), 2) # r upper bound
# b = 0.98, SE = 0.12, p < .001, beta = 0.45, r = .38, 95% CI [.29, .46]



# *--- FP predicting SS at a distance (tradeoff): intelligence & confidence ----

# Specify model
mod_dv3_fp <- "
fp_intel_lat =~ fp_intel_p1 + fp_intel_p2 + fp_intel_p3 + fp_intel_p4
fp_conf_lat =~ fp_conf_p1 + fp_conf_p2 + fp_conf_p3 + fp_conf_p4

dv3_ic ~ fp_intel_lat + fp_conf_lat

fp_intel_p1 ~~ fp_conf_p1
fp_intel_p2 ~~ fp_conf_p2
fp_intel_p3 ~~ fp_conf_p3
fp_intel_p4 ~~ fp_conf_p4
"

# Fit model
fit_dv3_fp <- sem(model = mod_dv3_fp, data = ss, std.lv = T)

# Estimated model
summary(fit_dv3_fp, fit.measures = T, standardized = T)

# Fit indices
round(fitMeasures(fit_dv3_fp)[c("chisq", "df", "cfi", "tli", "rmsea")], 2)
round(fitMeasures(fit_dv3_fp)["pvalue"], 3)
# Chi-squared(21) = 264.36, p < .001, CFI = 0.91, TLI = 0.84, RMSEA = 0.14

# Get parameter estimates
par_dv3_fp <- parameterestimates(fit_dv3_fp, standardized = T)

# Estimate coefficients: Intelligence
# Note that lower value of DV = greater preference of "intelligent" website
# so the signs of all coefficients need to be flipped
round(par_dv3_fp[9, c("est", "std.all")]*(-1), 2) # b, beta
round(par_dv3_fp[9, "se"], 2) # SE
round(par_dv3_fp[9, "pvalue"], 3) # p-value
round(par_dv3_fp[9, "std.all"]*(-1)*
        sqrt(1 - (par_dv3_fp[26, "std.all"]*(-1))^2), 2) # r
# b = 0.69, SE = 0.17, p < .001, beta = 0.32, r = .20

# Estimate coefficient: Confidence
round(par_dv3_fp[10, c("est", "se", "std.all")], 2) # b, SE, beta
round(par_dv3_fp[10, "pvalue"], 3) # p-value
round(par_dv3_fp[10, "std.all"]*
        sqrt(1 - par_dv3_fp[26, "std.all"]^2), 2) # r
# b = 0.44, SE = 0.17, p = .009, beta = 0.20, r = .13



# *--- SP predicting SS at a distance (choice): intelligence & confidence ----

# Specify model
mod_dv4_sp <- "
sp_intel_lat =~ sp_intelligent + sp_smart + sp_sharp
sp_conf_lat =~ sp_confident + sp_assured

dv4_ic ~ sp_intel_lat + sp_conf_lat
"

# Fit model
fit_dv4_sp <- sem(model = mod_dv4_sp, ordered = "dv4_ic", data = ss, std.lv = T)

# Estimated model
summary(fit_dv4_sp, fit.measures = T, standardized = T)

# Fit indices
round(fitMeasures(fit_dv4_sp)[c("chisq", "df", "cfi", "tli", "rmsea")], 2)
round(fitMeasures(fit_dv4_sp)["pvalue"], 3)
# Chi-squared(7) = 2.55, p = .923, CFI = 1.00, TLI = 1.01, RMSEA = 0.00

# Get parameter estimates
par_dv4_sp <- parameterestimates(fit_dv4_sp)

# Estimate coefficients: Intelligence
# Note that because choice of "intelligent" website was coded as 0,
# the signs of all coefficients need to be flipped
round(par_dv4_sp[6, "est"]*1.7*(-1), 2) # b (probit to logit)
round(par_dv4_sp[6, "se"]*1.7, 2) # SE (probit to logit)
round(par_dv4_sp[6, "pvalue"], 3) # p-value
round(exp(par_dv4_sp[6, "est"]*1.7*(-1)), 2) # odds ratio
round(lor2r(par_dv4_sp[6, "est"]*1.7*(-1)), 2) # r
round(lor2r(par_dv4_sp[6, "ci.upper"]*1.7*(-1)), 2) # r lower (b/c sign flip)
round(lor2r(par_dv4_sp[6, "ci.lower"]*1.7*(-1)), 2) # r upper (b/c sign flip)
# b = 1.04, SE = 0.10, p < .001, OR = 2.84, r = .28, 95% CI [.23, .32]

# Estimate coefficients: Confidence
round(par_dv4_sp[7, c("est", "se")]*1.7, 2) # b, SE (probit to logit)
round(par_dv4_sp[7, "pvalue"], 3) # p-value
round(exp(par_dv4_sp[7, "est"]*1.7), 2) # odds ratio
round(lor2r(par_dv4_sp[7, "est"]*1.7), 2) # r
round(lor2r(par_dv4_sp[7, "ci.lower"]*1.7), 2) # r lower
round(lor2r(par_dv4_sp[7, "ci.upper"]*1.7), 2) # r upper
# b = 1.00, SE = 0.13, p < .001, OR = 2.71, r = .26, 95% CI [.20, .33]



# *--- FP predicting SS at a distance (choice): intelligence & confidence ----

# Specify model
mod_dv4_fp <- "
fp_intel_lat =~ fp_intel_p1 + fp_intel_p2 + fp_intel_p3 + fp_intel_p4
fp_conf_lat =~ fp_conf_p1 + fp_conf_p2 + fp_conf_p3 + fp_conf_p4

dv4_ic ~ fp_intel_lat + fp_conf_lat

fp_intel_p1 ~~ fp_conf_p1
fp_intel_p2 ~~ fp_conf_p2
fp_intel_p3 ~~ fp_conf_p3
fp_intel_p4 ~~ fp_conf_p4
"

# Fit model
fit_dv4_fp <- sem(model = mod_dv4_fp, ordered = "dv4_ic", data = ss, std.lv = T)

# Estimated model
summary(fit_dv4_fp, fit.measures = T, standardized = T)

# Fit indices
round(fitMeasures(fit_dv4_fp)[c("chisq", "df", "cfi", "tli", "rmsea")], 2)
round(fitMeasures(fit_dv4_fp)["pvalue"], 3)
# Chi-squared(21) = 53.89, p < .001, CFI = 0.99, TLI = 0.98, RMSEA = 0.05

# Get parameter estimates
par_dv4_fp <- parameterestimates(fit_dv4_fp)

# Estimate coefficients: Intelligence
# Note that because choice of "intelligent" website was coded as 0,
# the signs of all coefficients need to be flipped
round(par_dv4_fp[9, "est"]*1.7*(-1), 2) # b (probit to logit)
round(par_dv4_fp[9, "se"]*1.7, 2) # SE (probit to logit)
round(par_dv4_fp[9, "pvalue"], 3) # p-value
round(exp(par_dv4_fp[9, "est"]*1.7*(-1)), 2) # odds ratio
round(lor2r(par_dv4_fp[9, "est"]*1.7*(-1)), 2) # r
round(lor2r(par_dv4_fp[9, "ci.upper"]*1.7*(-1)), 2) # r lower (b/c sign flip)
round(lor2r(par_dv4_fp[9, "ci.lower"]*1.7*(-1)), 2) # r upper (b/c sign flip)
# b = 0.44, SE = 0.16, p = .005, OR = 1.56, r = .12, 95% CI [.04, .20]

# Estimate coefficients: Confidence
round(par_dv4_fp[10, c("est", "se")]*1.7, 2) # b, SE (probit to logit)
round(par_dv4_fp[10, "pvalue"], 3) # p-value
round(exp(par_dv4_fp[10, "est"]*1.7), 2) # odds ratio
round(lor2r(par_dv4_fp[10, "est"]*1.7), 2) # r
round(lor2r(par_dv4_fp[10, "ci.lower"]*1.7), 2) # r lower
round(lor2r(par_dv4_fp[10, "ci.upper"]*1.7), 2) # r upper
# b = 0.45, SE = 0.16, p = .004, OR = 1.57, r = .12, 95% CI [.04, .21]



# *- Bivariate Regression -------------------------------------------------


# *--- SP predicting SS at a distance (tradeoff): intelligence ------------

# Fit model
br_dv3_sp_intel <- lm.beta(lm(dv3_ic ~ spcomp_intel, data = ss))

# Estimate coefficients
# Note that lower value of DV = greater preference of "intelligent" website
# so the signs of all coefficients need to be flipped
round(summary(br_dv3_sp_intel)$
        coefficients["spcomp_intel", 1:2]*(-1), 2) # b, beta ( = r)
round(summary(br_dv3_sp_intel)$coefficients["spcomp_intel", 3], 2) # SE
round(summary(br_dv3_sp_intel)$coefficients["spcomp_intel", 5], 3) # p-value
# b = 0.66, beta = r = 0.31, SE = 0.09, p < .001



# *--- SP predicting SS at a distance (tradeoff): confidence --------------

# Fit model
br_dv3_sp_conf <- lm.beta(lm(dv3_ic ~ spcomp_conf, data = ss))

# Estimate coefficients
round(summary(br_dv3_sp_conf)$
        coefficients["spcomp_conf", 1:3], 2) # b, beta ( = r), SE
round(summary(br_dv3_sp_conf)$coefficients["spcomp_conf", 5], 3) # p-value
# b = 0.21, beta = r = 0.11, SE = 0.08, p = .010



# *--- FP predicting SS at a distance (tradeoff): intelligence ------------

# Fit model
br_dv3_fp_intel <- lm.beta(lm(dv3_ic ~ fpcomp_intel, data = ss))

# Estimate coefficients
# Note that lower value of DV = greater preference of "intelligent" website
# so the signs of all coefficients need to be flipped
round(summary(br_dv3_fp_intel)$
        coefficients["fpcomp_intel", 1:2]*(-1), 2) # b, beta ( = r)
round(summary(br_dv3_fp_intel)$coefficients["fpcomp_intel", 3], 2) # SE
summary(br_dv3_fp_intel)$coefficients["fpcomp_intel", 5] # p-value
# b = 1.62, beta = r = 0.15, SE = 0.46, p < .001



# *--- FP predicting SS at a distance (tradeoff): confidence --------------

# Fit model
br_dv3_fp_conf <- lm.beta(lm(dv3_ic ~ fpcomp_conf, data = ss))

# Estimate coefficients
round(summary(br_dv3_fp_conf)$
        coefficients["fpcomp_conf", 1:3], 2) # b, beta ( = r), SE
round(summary(br_dv3_fp_conf)$coefficients["fpcomp_conf", 5], 3) # p-value
# b = -0.37, beta = r = -0.04, SE = 0.41, p = .358



# *--- SP predicting SS at a distance (choice): intelligence --------------

# Fit model
br_dv4_sp_intel <- glm(dv4_ic ~ spcomp_intel, data = ss, family = "binomial")

# Estimate coefficients
# Note that because choice of "intelligent" website was coded as 0,
# the signs of all coefficients need to be flipped
round(summary(br_dv4_sp_intel)$coefficients["spcomp_intel", 1]*(-1), 2) # b
round(summary(br_dv4_sp_intel)$coefficients["spcomp_intel", 2], 2) # SE
round(summary(br_dv4_sp_intel)$coefficients["spcomp_intel", 4], 3) # p-value
round(exp(coef(br_dv4_sp_intel)["spcomp_intel"]*(-1)), 2) # odds ratio
round(lor2r(coef(br_dv4_sp_intel)["spcomp_intel"]*(-1)), 2) # r
# b = 0.44, SE = 0.09, p < .001, OR = 1.55, r = .12



# *--- SP predicting SS at a distance (choice): confidence ----------------

# Fit model
br_dv4_sp_conf <- glm(dv4_ic ~ spcomp_conf, data = ss, family = "binomial")

# Estimate coefficients
round(summary(br_dv4_sp_conf)$coefficients["spcomp_conf", 1:2], 2) # b, SE
summary(br_dv4_sp_conf)$coefficients["spcomp_conf", 4] # p-value
# Note that this p-value rounds up to .001 but is actually < .001
round(exp(coef(br_dv4_sp_conf)["spcomp_conf"]), 2) # odds ratio
round(lor2r(coef(br_dv4_sp_conf)["spcomp_conf"]), 2) # r
# b = 0.33, SE = 0.10, p < .001, OR = 1.39, r = .09



# *--- FP predicting SS at a distance (choice): intelligence --------------

# Fit model
br_dv4_fp_intel <- glm(dv4_ic ~ fpcomp_intel, data = ss, family = "binomial")

# Estimate coefficients
# Note that because choice of "intelligent" website was coded as 0,
# the signs of all coefficients need to be flipped
round(summary(br_dv4_fp_intel)$coefficients["fpcomp_intel", 1]*(-1), 2) # b
round(summary(br_dv4_fp_intel)$coefficients["fpcomp_intel", 2], 2) # SE
round(summary(br_dv4_fp_intel)$coefficients["fpcomp_intel", 4], 3) # p-value
round(exp(coef(br_dv4_fp_intel)["fpcomp_intel"]*(-1)), 2) # odds ratio
round(lor2r(coef(br_dv4_fp_intel)["fpcomp_intel"]*(-1)), 2) # r
# b = 0.65, SE = 0.49, p = .190, OR = 1.91, r = .18



# *--- FP predicting SS at a distance (choice): confidence ----------------

# Fit model
br_dv4_fp_conf <- glm(dv4_ic ~ fpcomp_conf, data = ss, family = "binomial")

# Estimate coefficients
round(summary(br_dv4_fp_conf)$coefficients["fpcomp_conf", 1:2], 2) # b, SE
round(summary(br_dv4_fp_conf)$coefficients["fpcomp_conf", 4], 3) # p-value
round(exp(coef(br_dv4_fp_conf)["fpcomp_conf"]), 2) # odds ratio
round(lor2r(coef(br_dv4_fp_conf)["fpcomp_conf"]), 2) # r
# b = 0.47, SE = 0.43, p = .277, OR = 1.60, r = .13



# *- Multiple Regression --------------------------------------------------


# *--- SP predicting SS at a distance (tradeoff): intelligence & confidence ----

# Fit model
mr_dv3_sp <- lm.beta(lm(dv3_ic ~ spcomp_intel + spcomp_conf, data = ss))

# Estimate coefficient: Intelligence
# Note that sign of coefficient needs to be flipped
# because DV with lower value = greater preference of "intelligent" website
round(summary(mr_dv3_sp)$coefficients["spcomp_intel", 1:2]*(-1), 2) # b, beta
round(summary(mr_dv3_sp)$coefficients["spcomp_intel", 3], 2) # SE
round(summary(mr_dv3_sp)$coefficients["spcomp_intel", 5], 3) # p-value
round(coef(mr_dv3_sp)["spcomp_intel"]*(-1)*sqrt(1 - cor_sp^2), 2) # r
# b = 0.99, beta = 0.47, SE = 0.09, p < .001, r = .41

# Estimate coefficient: Confidence
round(summary(mr_dv3_sp)$coefficients["spcomp_conf", 1:3], 2) # b, beta, SE
round(summary(mr_dv3_sp)$coefficients["spcomp_conf", 5], 3) # p-value
round(coef(mr_dv3_sp)["spcomp_conf"]*sqrt(1- cor_sp^2), 2) # r
# b = 0.63, beta = 0.33, SE = 0.08, p < .001, r = .29



# *--- FP predicting SS at a distance (tradeoff): intelligence & confidence ----

# Fit model
mr_dv3_fp <- lm.beta(lm(dv3_ic ~ fpcomp_intel + fpcomp_conf, data = ss))

# Estimate coefficient: Intelligence
# Note that sign of coefficient needs to be flipped
# because DV with lower value = greater preference of "intelligent" website
round(summary(mr_dv3_fp)$coefficients["fpcomp_intel", 1:2]*(-1), 2) # b, beta
round(summary(mr_dv3_fp)$coefficients["fpcomp_intel", 3], 2) # SE
round(summary(mr_dv3_fp)$coefficients["fpcomp_intel", 5], 3) # p-value
round(coef(mr_dv3_fp)["fpcomp_intel"]*(-1)*sqrt(1- cor_fp^2), 2) # r
# b = 2.36, beta = 0.21, SE = 0.61, p < .001, r = .16

# Estimate coefficient: Confidence
round(summary(mr_dv3_fp)$coefficients["fpcomp_conf", 1:3], 2) # b, beta, SE
round(summary(mr_dv3_fp)$coefficients["fpcomp_conf", 5], 3) # p-value
round(coef(mr_dv3_fp)["fpcomp_conf"]*sqrt(1- cor_fp^2), 2) # r
# b = 0.98, beta = 0.10, SE = 0.53, p = .067, r = .08



# *--- SP predicting SS at a distance (choice): intelligence & confidence ----

# Fit model
mr_dv4_sp <- glm(dv4_ic ~ spcomp_intel + spcomp_conf, 
                 data = ss, family = "binomial")

# Estimate coefficient: Intelligence
# Note that because choice of "intelligent" website was coded as 0,
# the signs of all coefficients need to be flipped
round(summary(mr_dv4_sp)$coefficients["spcomp_intel", 1]*(-1), 2) # b
round(summary(mr_dv4_sp)$coefficients["spcomp_intel", 2], 2) # SE
round(summary(mr_dv4_sp)$coefficients["spcomp_intel", 4], 3) # p-value
round(exp(coef(mr_dv4_sp)["spcomp_intel"]*(-1)), 2) # odds ratio
round(lor2r(coef(mr_dv4_sp)["spcomp_intel"]*(-1)), 2) # r
# b = 0.94, SE = 0.13, p < .001, OR = 2.55, r = .25

# Estimate coefficient: Confidence
round(summary(mr_dv4_sp)$coefficients["spcomp_conf", 1:2], 2) # b, SE
round(summary(mr_dv4_sp)$coefficients["spcomp_conf", 4], 3) # p-value
round(exp(coef(mr_dv4_sp)["spcomp_conf"]), 2) # odds ratio
round(lor2r(coef(mr_dv4_sp)["spcomp_conf"]), 2) # r
# b = 0.84, SE = 0.13, p < .001, OR = 2.33, r = .23



# *--- FP predicting SS at a distance (choice): intelligence & confidence ----

# Fit model
mr_dv4_fp <- glm(dv4_ic ~ fpcomp_intel + fpcomp_conf, 
                 data = ss, family = "binomial")

# Estimate coefficients: Intelligence
# Note that because choice of "intelligent" website was coded as 0,
# the signs of all coefficients need to be flipped
round(summary(mr_dv4_fp)$coefficients["fpcomp_intel", 1]*(-1), 2) # b
round(summary(mr_dv4_fp)$coefficients["fpcomp_intel", 2], 2) # SE
round(summary(mr_dv4_fp)$coefficients["fpcomp_intel", 4], 3) # p-value
round(exp(coef(mr_dv4_fp)["fpcomp_intel"]*(-1)), 2) # odds ratio
round(lor2r(coef(mr_dv4_fp)["fpcomp_intel"]*(-1)), 2) # r
# b = 1.81, SE = 0.67, p = .007, OR = 6.11, r = .45

# Estimate coefficients: Confidence
round(summary(mr_dv4_fp)$coefficients["fpcomp_conf", 1:2], 2) # b, SE
round(summary(mr_dv4_fp)$coefficients["fpcomp_conf", 4], 3) # p-value
round(exp(coef(mr_dv4_fp)["fpcomp_conf"]), 2) # odds ratio
round(lor2r(coef(mr_dv4_fp)["fpcomp_conf"]), 2) # r
# b = 1.52, SE = 0.59, p = .010, OR = 4.55, r = .39

