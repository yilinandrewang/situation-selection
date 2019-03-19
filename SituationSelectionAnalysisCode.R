# ******* IDEAS ABOUT LIKING PREDICT SITUATION SELECTION AT A DISTANCE *******
# *******     Y. ANDRE WANG, PAUL W. EASTWICK, & ALISON LEDGERWOOD     *******
# ****************************************************************************

# ANALYSIS CODE
# WRITTEN BY Y. ANDRE WANG
# LAST UPDATE: 03/17/2019


# ********************************* INTRO ************************************

# This code file will produce all results reported in the results section of
# Wang, Eastwick, and Ledgerwood (under review), as well as all results
# reported in Table S2 in the supplemental materials.
# To obtain code for additional analyses, please contact the first author
# Email: ylawang@ucdavis.edu

# This code file is organized into two sections:
# Section I reports all tests with summarized preferences (SP) as predictors.
# Section II reports all tests with functional preferences (FP) as predictors.

# Each code section includes tests conducted on the four DVs.
# The four DVs are reported in the same order they appear in the manuscript:
# DV 1: Situation selection at a distance
# DV 2: Situation selection with experience
# DV 3: Situation selection at a distance (tradeoff)
# DV 4: Situation selection at a distance (choice)

# With each DV, three analytic approaches were conducted and reported:
# Analytic approach 1: Structural Equation Modeling (SEM)
# Analytic approach 2: Bivariate Regressions (BR)
# Analytic approach 3: Multiple Regressions (MR)


############################## DATA PREPARATION ###########################

# Load R packages used for data analysis
# Use the function install.packages() to install packages if not yet in R
library(lm.beta); library(lavaan); library(semPlot)

# Load data file
ss <- read.csv("SituationSelectionData.csv")

# Create function that converts log odds ratios to correlation coefficients
# Formula from Borenstein, Hedges, Higgins, and Rothstein (2009)
lor2r <- function(x) {
  x*(sqrt(3) / pi)/
    sqrt((x*(sqrt(3)/pi))^2 + 4)
}

# Calculate raw correlations between SPs and between FPs for use below
cor_sp <- cor(ss$spcomp_intel, ss$spcomp_conf)
cor_fp <- cor(ss$fpcomp_intel, ss$fpcomp_conf)

# Label indicators in SEM for graphical purposes
lab_sp <- c("intelligent", "smart", "sharp", "confident", "assured")
lab_fp <- c("intel1", "intel2", "intel3", "intel4", 
            "conf1", "conf2", "conf3", "conf4")

########################## SECTION I: SP PREDICTIONS ######################
################# DV 1: Situation selection at a distance #################

# *-- SEM ####

# *---- SP: Intelligence ####

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
# Chi-squared(7) = 10.34, p = 0.170, CFI = 1.00, TLI = 1.00, RMSEA = 0.03

# Regression coefficient
par_dv1_sp_intel <- parameterestimates(fit_dv1_sp_intel, standardized = T)
round(par_dv1_sp_intel[6, c("est", "std.all", "se")], 2) # b, beta, SE
round(par_dv1_sp_intel[6, "pvalue"], 3) # p-value
round(par_dv1_sp_intel[6, "std.all"]*
        sqrt(1 - par_dv1_sp_intel[16, "std.all"]^2), 2) # r
# b = 0.90, SE = 0.13, p < .001, beta = 0.39, r = 0.32

# Visualize model
semPaths(fit_dv1_sp_intel, "std", fade = F, esize = T, edge.color = 'black', 
         edge.label.cex = 1, sizeMan = 10, sizeLat = 10, rotation = 2,
         nodeLabels = c(lab_sp, "DV1_intel", "SP_intel", "SP_conf"),
         label.norm = "OOOOO")


# *---- SP: Confidence ####

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
# Chi-squared(7) = 8.96, p = 0.256, CFI = 1.00, TLI = 1.00, RMSEA = 0.02

# Regression coefficient
par_dv1_sp_conf <- parameterestimates(fit_dv1_sp_conf, standardized = T)
round(par_dv1_sp_conf[7, c("est", "std.all", "se")], 2) # b, beta, SE
round(par_dv1_sp_conf[7, "pvalue"], 3) # p-value
round(par_dv1_sp_conf[7, "std.all"]*
        sqrt(1 - par_dv1_sp_conf[16, "std.all"]^2), 2) # r
# b = 0.86, SE = 0.14, p < .001, beta = 0.38, r = 0.31

# Visualize model
semPaths(fit_dv1_sp_conf, "std", fade = F, esize = T, edge.color = 'black', 
         edge.label.cex = 1, sizeMan = 10, sizeLat = 10, rotation = 2,
         nodeLabels = c(lab_sp, "DV1_conf", "SP_intel", "SP_conf"),
         label.norm = "OOOOO")



# *-- Bivariate Regression ####

# *---- SP: Intelligence ####
br_dv1_sp_intel <- lm.beta(lm(dv1_intel ~ spcomp_intel, data = ss))
round(summary(br_dv1_sp_intel)$
        coefficients["spcomp_intel", 1:3], 2) # b, beta ( = r), SE
round(summary(br_dv1_sp_intel)$coefficients["spcomp_intel", 5], 3) # p-value
# b = 0.73, beta = 0.33, SE = 0.09, p < .001


# *---- SP: Confidence ####
br_dv1_sp_conf <- lm.beta(lm(dv1_conf ~ spcomp_conf, data = ss))
round(summary(br_dv1_sp_conf)$
        coefficients["spcomp_conf", 1:3], 2) # b, beta ( = r), SE
round(summary(br_dv1_sp_conf)$coefficients["spcomp_conf", 5], 3) # p-value
# b = 0.59, beta = 0.29, SE = 0.08, p < .001



# *-- Multiple Regression ####

# *---- SP: Intelligence ####
mr_dv1_sp_intel <- lm.beta(lm(dv1_intel ~ spcomp_intel + spcomp_conf, 
                              data = ss))
round(summary(mr_dv1_sp_intel)$
        coefficients["spcomp_intel", 1:3], 2) # b, beta, SE
round(summary(mr_dv1_sp_intel)$coefficients["spcomp_intel", 5], 3) # p-value
round(summary(mr_dv1_sp_intel)$coefficients["spcomp_intel", 2]*
        sqrt(1- cor_sp^2), 2) # r
# b = 0.75, beta = 0.34, SE = 0.10, p < .001, r = .30


# *---- SP: Confidence ####
mr_dv1_sp_conf <- lm.beta(lm(dv1_conf ~ spcomp_conf + spcomp_intel, 
                              data = ss))
round(summary(mr_dv1_sp_conf)$
        coefficients["spcomp_conf", 1:3], 2) # b, beta, SE
round(summary(mr_dv1_sp_conf)$coefficients["spcomp_conf", 5], 3) # p-value
round(summary(mr_dv1_sp_conf)$coefficients["spcomp_conf", 2]*
        sqrt(1- cor_sp^2), 2) # r
# b = 0.62, beta = 0.31, SE = 0.09, p < .001, r = .27




################# DV 2: Situation selection with experience #################

# *-- SEM ####

# *---- SP: Intelligence & Confidence ####

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
# Chi-squared(7) = 1.23, p = 0.990, CFI = 1.00, TLI = 1.01, RMSEA = 0.00

# Regression coefficient: Intelligence
# Transform probit coefficients into logit coefficients
par_dv2_sp <- parameterestimates(fit_dv2_sp)
round(par_dv2_sp[6, c("est", "se")]*1.7, 2) # b, SE (probit to logit)
par_dv2_sp[6, "pvalue"] # p-value
round(exp(par_dv2_sp[6, "est"]*1.7), 2) # odds ratio
round(lor2r(par_dv2_sp[6, "est"]*1.7), 2) # r
# b = 0.44, SE = 0.13, p < .001, OR = 1.55, r = 0.12

# Regression coefficient: Confidence
# Note that sign of coefficient needs to be flipped
# because choice of "confident" website was coded as 0
round(par_dv2_sp[7, "est"]*1.7*(-1), 2) # b (probit to logit)
round(par_dv2_sp[7, "se"]*1.7, 2) # SE (probit to logit)
round(par_dv2_sp[7, "pvalue"], 3) # p-value
round(exp(par_dv2_sp[7, "est"]*1.7*(-1)), 2) # odds ratio
round(lor2r(par_dv2_sp[7, "est"]*1.7*(-1)), 2) # r
# b = 0.54, SE = 0.13, p < .001, OR = 1.71, r = 0.15

# Visualize model (note that the regression coefficients are in probit)
semPaths(fit_dv2_sp, "std", fade = F, esize = T, edge.color = 'black', 
         edge.label.cex = 1, sizeMan = 10, sizeLat = 10, rotation = 2,
         nodeLabels = c(lab_sp, "DV2", "SP_intel", "SP_conf"),
         label.norm = "OOOOO", intercepts = F, thresholds = F)



# *-- Bivariate Regression ####

# *---- SP: Intelligence ####
br_dv2_sp_intel <- glm(dv2_ci ~ spcomp_intel, data = ss, family = "binomial")
round(summary(br_dv2_sp_intel)$coefficients["spcomp_intel", 1:2], 2) # b, SE
round(summary(br_dv2_sp_intel)$coefficients["spcomp_intel", 4], 3) # p-value
round(exp(coef(br_dv2_sp_intel)["spcomp_intel"]), 2) # odds ratio
round(lor2r(coef(br_dv2_sp_intel)["spcomp_intel"]), 2) # r
# b = 0.13, SE = 0.10, p = .192, OR = 1.13, r = 0.03


# *---- SP: Confidence ####
# Note that sign of coefficient needs to be flipped
# because choice of "confident" website was coded as 0
br_dv2_sp_conf <- glm(dv2_ci ~ spcomp_conf, data = ss, family = "binomial")
round(summary(br_dv2_sp_conf)$coefficients["spcomp_conf", 1]*(-1), 2) # b
round(summary(br_dv2_sp_conf)$coefficients["spcomp_conf", 2], 2) # SE
round(summary(br_dv2_sp_conf)$coefficients["spcomp_conf", 4], 3) # p-value
round(exp(coef(br_dv2_sp_conf)["spcomp_conf"]*(-1)), 2) # odds ratio
round(lor2r(coef(br_dv2_sp_conf)["spcomp_conf"]*(-1)), 2) # r
# b = 0.23, SE = 0.08, p = 0.006, OR = 1.25, r = 0.06



# *-- Multiple Regression ####

# *---- SP: Intelligence & Confidence ####
mr_dv2_sp <- glm(dv2_ci ~ spcomp_intel + spcomp_conf, 
                 data = ss, family = "binomial")

# Regression coefficient: Intelligence
round(summary(mr_dv2_sp)$coefficients["spcomp_intel", 1:2], 2) # b, SE
round(summary(mr_dv2_sp)$coefficients["spcomp_intel", 4], 3) # p-value
round(exp(coef(mr_dv2_sp)["spcomp_intel"]), 2) # odds ratio
round(lor2r(coef(mr_dv2_sp)["spcomp_intel"]), 2) # r
# b = 0.32, SE = 0.11, p = 0.004, OR = 1.38, r = 0.09

# Regression coefficient: Confidence
# Note that sign of coefficient needs to be flipped
# because choice of "confident" website was coded as 0
round(summary(mr_dv2_sp)$coefficients["spcomp_conf", 1]*(-1), 2) # b
round(summary(mr_dv2_sp)$coefficients["spcomp_conf", 2], 2) # SE
round(summary(mr_dv2_sp)$coefficients["spcomp_conf", 4], 3) # p-value
round(exp(coef(mr_dv2_sp)["spcomp_conf"]*(-1)), 2) # odds ratio
round(lor2r(coef(mr_dv2_sp)["spcomp_conf"]*(-1)), 2) # r
# b = 0.36, SE = 0.10, p < 0.001, OR = 1.44, r = 0.10




################# DV 3: Situation selection at a distance (tradeoff) ##########

# *-- SEM ####

# *---- SP: Intelligence & Confidence ####

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
# Chi-squared(7) = 7.81, p = 0.350, CFI = 1.00, TLI = 1.00, RMSEA = 0.01

# Regression coefficient: Intelligence
# Note that sign of coefficient needs to be flipped
# because DV with lower value = greater preference of "intelligent" website
par_dv3_sp <- parameterestimates(fit_dv3_sp, standardized = T)
round(par_dv3_sp[6, c("est", "std.all")]*(-1), 2) # b, beta
round(par_dv3_sp[6, "se"], 2) # SE
round(par_dv3_sp[6, "pvalue"], 3) # p-value
round(par_dv3_sp[6, "std.all"]*(-1)*
        sqrt(1 - (par_dv3_sp[16, "std.all"]*(-1))^2), 2) # r
# b = 1.27, SE = 0.12, p < .001, beta = 0.58, r = 0.48

# Regression coefficient: Confidence
round(par_dv3_sp[7, c("est", "std.all", "se")], 2) # b, beta, SE
round(par_dv3_sp[7, "pvalue"], 3) # p-value
round(par_dv3_sp[7, "std.all"]*
        sqrt(1 - par_dv3_sp[16, "std.all"]^2), 2) # r
# b = 0.98, SE = 0.12, p < .001, beta = 0.45, r = 0.38

# Visualize model
semPaths(fit_dv3_sp, "std", fade = F, esize = T, edge.color = 'black', 
         edge.label.cex = 1, sizeMan = 10, sizeLat = 10, rotation = 2,
         nodeLabels = c(lab_sp, "DV3", "SP_intel", "SP_conf"),
         label.norm = "OOOOO")



# *-- Bivariate Regression ####

# *---- SP: Intelligence ####
# Note that sign of coefficient needs to be flipped
# because DV with lower value = greater preference of "intelligent" website
br_dv3_sp_intel <- lm.beta(lm(dv3_ic ~ spcomp_intel, data = ss))
round(summary(br_dv3_sp_intel)$
        coefficients["spcomp_intel", 1:2]*(-1), 2) # b, beta ( = r)
round(summary(br_dv3_sp_intel)$coefficients["spcomp_intel", 3], 2) # SE
round(summary(br_dv3_sp_intel)$coefficients["spcomp_intel", 5], 3) # p-value
# b = 0.66, beta = 0.31, SE = 0.09, p < .001


# *---- SP: Confidence ####
br_dv3_sp_conf <- lm.beta(lm(dv3_ic ~ spcomp_conf, data = ss))
round(summary(br_dv3_sp_conf)$
        coefficients["spcomp_conf", 1:3], 2) # b, beta ( = r), SE
round(summary(br_dv3_sp_conf)$coefficients["spcomp_conf", 5], 3) # p-value
# b = 0.21, beta = 0.11, SE = 0.08, p = .010



# *-- Multiple Regression ####

# *---- SP: Intelligence & Confidence ####
mr_dv3_sp <- lm.beta(lm(dv3_ic ~ spcomp_intel + spcomp_conf, data = ss))

# Regression coefficient: Intelligence
# Note that sign of coefficient needs to be flipped
# because DV with lower value = greater preference of "intelligent" website
round(summary(mr_dv3_sp)$coefficients["spcomp_intel", 1:2]*(-1), 2) # b, beta
round(summary(mr_dv3_sp)$coefficients["spcomp_intel", 3], 2) # SE
round(summary(mr_dv3_sp)$coefficients["spcomp_intel", 5], 3) # p-value
round(coef(mr_dv3_sp)["spcomp_intel"]*(-1)*sqrt(1- cor_sp^2), 2) # r
# b = 0.99, beta = 0.47, SE = 0.09, p < 0.001, r = .41

# Regression coefficient: Confidence
round(summary(mr_dv3_sp)$coefficients["spcomp_conf", 1:3], 2) # b, beta, SE
round(summary(mr_dv3_sp)$coefficients["spcomp_conf", 5], 3) # p-value
round(coef(mr_dv3_sp)["spcomp_conf"]*sqrt(1- cor_sp^2), 2) # r
# b = 0.63, beta = 0.33, SE = 0.08, p < 0.001, r = .29




################# DV 4: Situation selection at a distance (choice) ############

# *-- SEM ####

# *---- SP: Intelligence & Confidence ####

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
# Chi-squared(7) = 2.55, p = 0.923, CFI = 1.00, TLI = 1.01, RMSEA = 0.00

# Regression coefficient: Intelligence
# Note that sign of coefficient needs to be flipped
# because choice of "intelligent" website was coded as 0
par_dv4_sp <- parameterestimates(fit_dv4_sp)
round(par_dv4_sp[6, "est"]*1.7*(-1), 2) # b (probit to logit)
round(par_dv4_sp[6, "se"]*1.7, 2) # SE (probit to logit)
round(par_dv4_sp[6, "pvalue"], 3) # p-value
round(exp(par_dv4_sp[6, "est"]*1.7*(-1)), 2) # odds ratio
round(lor2r(par_dv4_sp[6, "est"]*1.7*(-1)), 2) # r
# b = 1.04, SE = 0.10, p < .001, OR = 2.84, r = 0.28

# Regression coefficient: Confidence
round(par_dv4_sp[7, c("est", "se")]*1.7, 2) # b, SE (probit to logit)
round(par_dv4_sp[7, "pvalue"], 3) # p-value
round(exp(par_dv4_sp[7, "est"]*1.7), 2) # odds ratio
round(lor2r(par_dv4_sp[7, "est"]*1.7), 2) # r
# b = 1.00, SE = 0.13, p < .001, OR = 2.71, r = 0.26

# Visualize model (note that the regression coefficients are in probit)
semPaths(fit_dv4_sp, "std", fade = F, esize = T, edge.color = 'black', 
         edge.label.cex = 1, sizeMan = 10, sizeLat = 10, rotation = 2,
         nodeLabels = c(lab_sp, "DV4", "SP_intel", "SP_conf"),
         label.norm = "OOOOO", intercepts = F, thresholds = F)



# *-- Bivariate Regression ####

# *---- SP: Intelligence ####
# Note that sign of coefficient needs to be flipped
# because choice of "intelligent" website was coded as 0
br_dv4_sp_intel <- glm(dv4_ic ~ spcomp_intel, data = ss, family = "binomial")
round(summary(br_dv4_sp_intel)$coefficients["spcomp_intel", 1]*(-1), 2) # b
round(summary(br_dv4_sp_intel)$coefficients["spcomp_intel", 2], 2) # SE
round(summary(br_dv4_sp_intel)$coefficients["spcomp_intel", 4], 3) # p-value
round(exp(coef(br_dv4_sp_intel)["spcomp_intel"]*(-1)), 2) # odds ratio
round(lor2r(coef(br_dv4_sp_intel)["spcomp_intel"]*(-1)), 2) # r
# b = 0.44, SE = 0.09, p < .001, OR = 1.55, r = 0.12


# *---- SP: Confidence ####
br_dv4_sp_conf <- glm(dv4_ic ~ spcomp_conf, data = ss, family = "binomial")
round(summary(br_dv4_sp_conf)$coefficients["spcomp_conf", 1:2], 2) # b, SE
summary(br_dv4_sp_conf)$coefficients["spcomp_conf", 4] # p-value
# Note that this p-value rounds up to .001 but is actually < .001
round(exp(coef(br_dv4_sp_conf)["spcomp_conf"]), 2) # odds ratio
round(lor2r(coef(br_dv4_sp_conf)["spcomp_conf"]), 2) # r
# b = 0.33, SE = 0.10, p < 0.001, OR = 1.39, r = 0.09



# *-- Multiple Regression ####

# *---- SP: Intelligence & Confidence ####
mr_dv4_sp <- glm(dv4_ic ~ spcomp_intel + spcomp_conf, 
                 data = ss, family = "binomial")

# Regression coefficient: Intelligence
# Note that sign of coefficient needs to be flipped
# because choice of "intelligent" website was coded as 0
round(summary(mr_dv4_sp)$coefficients["spcomp_intel", 1]*(-1), 2) # b
round(summary(mr_dv4_sp)$coefficients["spcomp_intel", 2], 2) # SE
round(summary(mr_dv4_sp)$coefficients["spcomp_intel", 4], 3) # p-value
round(exp(coef(mr_dv4_sp)["spcomp_intel"]*(-1)), 2) # odds ratio
round(lor2r(coef(mr_dv4_sp)["spcomp_intel"]*(-1)), 2) # r
# b = 0.94, SE = 0.13, p < 0.001, OR = 2.55, r = 0.25

# Regression coefficient: Confidence
round(summary(mr_dv4_sp)$coefficients["spcomp_conf", 1:2], 2) # b, SE
round(summary(mr_dv4_sp)$coefficients["spcomp_conf", 4], 3) # p-value
round(exp(coef(mr_dv4_sp)["spcomp_conf"]), 2) # odds ratio
round(lor2r(coef(mr_dv4_sp)["spcomp_conf"]), 2) # r
# b = 0.84, SE = 0.13, p < 0.001, OR = 2.33, r = 0.23




########################## SECTION II: FP PREDICTIONS #####################
################# DV 1: Situation selection at a distance #################

# *-- SEM ####

# *---- FP: Intelligence ####

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
# Chi-squared(21) = 265.30, p < 0.001, CFI = 0.90, TLI = 0.84, RMSEA = 0.14

# Regression coefficient
par_dv1_fp_intel <- parameterestimates(fit_dv1_fp_intel, standardized = T)
round(par_dv1_fp_intel[9, c("est", "se", "std.all")], 2)
round(par_dv1_fp_intel[9, "pvalue"], 3)
round(par_dv1_fp_intel[9, "std.all"]*
        sqrt(1 - par_dv1_fp_intel[26, "std.all"]^2), 2) # Calculate r
# b = 0.40, SE = 0.18, p = 0.026, beta = 0.17, r = 0.11

# Visualize model
semPaths(fit_dv1_fp_intel, "std", fade = F, esize = T, edge.color = 'black', 
         edge.label.cex = 1, sizeMan = 7, sizeLat = 10, rotation = 2,
         nodeLabels = c(lab_fp, "DV1_intel", "FP_intel", "FP_conf"), 
         label.norm = "OOOO", curvature = 3)


# *---- FP: Confidence ####

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
# Chi-squared(21) = 263.68, p < 0.001, CFI = 0.90, TLI = 0.84, RMSEA = 0.14

# Regression coefficient
par_dv1_fp_conf <- parameterestimates(fit_dv1_fp_conf, standardized = T)
round(par_dv1_fp_conf[10, c("est", "se", "std.all")], 2)
round(par_dv1_fp_conf[10, "pvalue"], 3)
round(par_dv1_fp_conf[10, "std.all"]*
        sqrt(1 - par_dv1_fp_conf[26, "std.all"]^2), 2) # Calculate r
# b = 0.24, SE = 0.18, p = 0.174, beta = 0.10, r = 0.07

# Visualize model
semPaths(fit_dv1_fp_conf, "std", fade = F, esize = T, edge.color = 'black', 
         edge.label.cex = 1, sizeMan = 7, sizeLat = 10, rotation = 2,
         nodeLabels = c(lab_fp, "DV1_conf", "FP_intel", "FP_conf"), 
         label.norm = "OOOO", curvature = 3)



# *-- Bivariate Regression ####

# *---- FP: Intelligence ####
br_dv1_fp_intel <- lm.beta(lm(dv1_intel ~ fpcomp_intel, data = ss))
round(summary(br_dv1_fp_intel)$
        coefficients["fpcomp_intel", 1:3], 2) # b, beta ( = r), SE
round(summary(br_dv1_fp_intel)$coefficients["fpcomp_intel", 5], 3) # p-value
# b = 1.06, beta = 0.09, SE = 0.49, p = 0.031


# *---- FP: Confidence ####
br_dv1_fp_conf <- lm.beta(lm(dv1_conf ~ fpcomp_conf, data = ss))
round(summary(br_dv1_fp_conf)$
        coefficients["fpcomp_conf", 1:3], 2) # b, beta ( = r), SE
round(summary(br_dv1_fp_conf)$coefficients["fpcomp_conf", 5], 3) # p-value
# b = 0.39, beta = 0.04, SE = 0.43, p = 0.366



# *-- Multiple Regression ####

# *---- FP: Intelligence ####
mr_dv1_fp_intel <- lm.beta(lm(dv1_intel ~ fpcomp_intel + fpcomp_conf, 
                              data = ss))
round(summary(mr_dv1_fp_intel)$
        coefficients["fpcomp_intel", 1:3], 2) # b, beta, SE
round(summary(mr_dv1_fp_intel)$coefficients["fpcomp_intel", 5], 3) # p-value
round(summary(mr_dv1_fp_intel)$coefficients["fpcomp_intel", 2]*
        sqrt(1- cor_fp^2), 2) # r
# b = 1.36, beta = 0.12, SE = 0.66, p = 0.039, r = .09


# *---- FP: Confidence ####
mr_dv1_fp_conf <- lm.beta(lm(dv1_conf ~ fpcomp_conf + fpcomp_intel, 
                             data = ss))
round(summary(mr_dv1_fp_conf)$
        coefficients["fpcomp_conf", 1:3], 2) # b, beta, SE
round(summary(mr_dv1_fp_conf)$coefficients["fpcomp_conf", 5], 3) # p-value
round(summary(mr_dv1_fp_conf)$coefficients["fpcomp_conf", 2]*
        sqrt(1- cor_fp^2), 2) # r
# b = 0.74, beta = 0.07, SE = 0.57, p = 0.193, r = .06




################# DV 2: Situation selection with experience #################

# *-- SEM ####

# *---- FP: Intelligence & Confidence ####

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
# Chi-squared(21) = 125.86, p < 0.001, CFI = 0.96, TLI = 0.94, RMSEA = 0.09

# Regression coefficient: Intelligence
par_dv2_fp <- parameterestimates(fit_dv2_fp)
round(par_dv2_fp[9, c("est", "se")]*1.7, 2) # b, SE (probit to logit)
round(par_dv2_fp[9, "pvalue"], 3) # p-value
round(exp(par_dv2_fp[9, "est"]*1.7), 2) # odds ratio
round(lor2r(par_dv2_fp[9, "est"]*1.7), 2) # r
# b = 1.06, SE = 0.13, p < .001, OR = 2.87, r = 0.28

# Regression coefficient: Confidence
# Note that sign of coefficient needs to be flipped
# because choice of "confident" website was coded as 0
round(par_dv2_fp[10, "est"]*1.7*(-1), 2) # b (probit to logit)
round(par_dv2_fp[10, "se"]*1.7, 2) # SE (probit to logit)
round(par_dv2_fp[10, "pvalue"], 3) # p-value
round(exp(par_dv2_fp[10, "est"]*1.7*(-1)), 2) # odds ratio
round(lor2r(par_dv2_fp[10, "est"]*1.7*(-1)), 2) # r
# b = 1.44, SE = 0.13, p < .001, OR = 4.21, r = 0.37

# Visualize model (note that the regression coefficients are in probit)
semPaths(fit_dv2_fp, "std", fade = F, esize = T, edge.color = 'black',
         edge.label.cex = 1, sizeMan = 7, sizeLat = 10, rotation = 2,
         nodeLabels = c(lab_fp, "DV2", "FP_intel", "FP_conf"), 
         label.norm = "OOOO", curvature = 3, intercepts = F, thresholds = F)



# *-- Bivariate Regression ####

# *---- FP: Intelligence ####
br_dv2_fp_intel <- glm(dv2_ci ~ fpcomp_intel, data = ss, family = "binomial")
round(summary(br_dv2_fp_intel)$coefficients["fpcomp_intel", 1:2], 2) # b, SE
round(summary(br_dv2_fp_intel)$coefficients["fpcomp_intel", 4], 3) # p-value
round(exp(coef(br_dv2_fp_intel)["fpcomp_intel"]), 2) # odds ratio
round(lor2r(coef(br_dv2_fp_intel)["fpcomp_intel"]), 2) # r
# b = 0.46, SE = 0.48, p = .339, OR = 1.59, r = 0.13


# *---- FP: Confidence ####
# Note that sign of coefficient needs to be flipped
# because choice of "confident" website was coded as 0
br_dv2_fp_conf <- glm(dv2_ci ~ fpcomp_conf, data = ss, family = "binomial")
round(summary(br_dv2_fp_conf)$coefficients["fpcomp_conf", 1]*(-1), 2) # b
round(summary(br_dv2_fp_conf)$coefficients["fpcomp_conf", 2], 2) # SE
round(summary(br_dv2_fp_conf)$coefficients["fpcomp_conf", 4], 3) # p-value
round(exp(coef(br_dv2_fp_conf)["fpcomp_conf"]*(-1)), 2) # odds ratio
round(lor2r(coef(br_dv2_fp_conf)["fpcomp_conf"]*(-1)), 2) # r
# b = 2.95, SE = 0.46, p < 0.001, OR = 19.08, r = 0.63



# *-- Multiple Regression ####

# *---- FP: Intelligence & Confidence ####
mr_dv2_fp <- glm(dv2_ci ~ fpcomp_intel + fpcomp_conf, 
                 data = ss, family = "binomial")

# Regression coefficient: Intelligence
round(summary(mr_dv2_fp)$coefficients["fpcomp_intel", 1:2], 2) # b, SE
round(summary(mr_dv2_fp)$coefficients["fpcomp_intel", 4], 3) # p-value
round(exp(coef(mr_dv2_fp)["fpcomp_intel"]), 2) # odds ratio
round(lor2r(coef(mr_dv2_fp)["fpcomp_intel"]), 2) # r
# b = 5.48, SE = 0.78, p < 0.001, OR = 239.86, r = 0.83

# Regression coefficient: Confidence
# Note that sign of coefficient needs to be flipped
# because choice of "confident" website was coded as 0
round(summary(mr_dv2_fp)$coefficients["fpcomp_conf", 1]*(-1), 2) # b
round(summary(mr_dv2_fp)$coefficients["fpcomp_conf", 2], 2) # SE
round(summary(mr_dv2_fp)$coefficients["fpcomp_conf", 4], 3) # p-value
round(exp(coef(mr_dv2_fp)["fpcomp_conf"]*(-1)), 2) # odds ratio
round(lor2r(coef(mr_dv2_fp)["fpcomp_conf"]*(-1)), 2) # r
# b = 6.41, SE = 0.72, p < 0.001, OR = 605.52, r = 0.87




################# DV 3: Situation selection at a distance (tradeoff) ##########

# *-- SEM ####

# *---- FP: Intelligence & Confidence ####

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
# Chi-squared(21) = 264.36, p < 0.001, CFI = 0.91, TLI = 0.84, RMSEA = 0.14

# Regression coefficient: Intelligence
# Note that sign of coefficient needs to be flipped
# because DV with lower value = greater preference of "intelligent" website
par_dv3_fp <- parameterestimates(fit_dv3_fp, standardized = T)
round(par_dv3_fp[9, c("est", "std.all")]*(-1), 2) # b, beta
round(par_dv3_fp[9, "se"], 2) # SE
round(par_dv3_fp[9, "pvalue"], 3) # p-value
round(par_dv3_fp[9, "std.all"]*(-1)*
        sqrt(1 - (par_dv3_fp[26, "std.all"]*(-1))^2), 2) # r
# b = 0.69, SE = 0.17, p < .001, beta = 0.32, r = 0.20

# Regression coefficient: Confidence
round(par_dv3_fp[10, c("est", "std.all", "se")], 2) # b, beta, SE
round(par_dv3_fp[10, "pvalue"], 3) # p-value
round(par_dv3_fp[10, "std.all"]*
        sqrt(1 - par_dv3_fp[26, "std.all"]^2), 2) # r
# b = 0.44, SE = 0.17, p = 0.009, beta = 0.20, r = 0.13

# Visualize model
semPaths(fit_dv3_fp, "std", fade = F, esize = T, edge.color = 'black', 
         edge.label.cex = 1, sizeMan = 7, sizeLat = 10, rotation = 2,
         nodeLabels = c(lab_fp, "DV3", "FP_intel", "FP_conf"), 
         label.norm = "OOOO", curvature = 3)



# *-- Bivariate Regression ####

# *---- FP: Intelligence ####
# Note that sign of coefficient needs to be flipped
# because DV with lower value = greater preference of "intelligent" website
br_dv3_fp_intel <- lm.beta(lm(dv3_ic ~ fpcomp_intel, data = ss))
round(summary(br_dv3_fp_intel)$
        coefficients["fpcomp_intel", 1:2]*(-1), 2) # b, beta ( = r)
round(summary(br_dv3_fp_intel)$coefficients["fpcomp_intel", 3], 2) # SE
summary(br_dv3_fp_intel)$coefficients["fpcomp_intel", 5] # p-value
# b = 1.62, beta = 0.15, SE = 0.46, p < 0.001


# *---- FP: Confidence ####
br_dv3_fp_conf <- lm.beta(lm(dv3_ic ~ fpcomp_conf, data = ss))
round(summary(br_dv3_fp_conf)$
        coefficients["fpcomp_conf", 1:3], 2) # b, beta ( = r), SE
round(summary(br_dv3_fp_conf)$coefficients["fpcomp_conf", 5], 3) # p-value
# b = -0.37, beta = -0.04, SE = 0.41, p = 0.358



# *-- Multiple Regression ####

# *---- FP: Intelligence & Confidence ####
mr_dv3_fp <- lm.beta(lm(dv3_ic ~ fpcomp_intel + fpcomp_conf, data = ss))

# Regression coefficient: Intelligence
# Note that sign of coefficient needs to be flipped
# because DV with lower value = greater preference of "intelligent" website
round(summary(mr_dv3_fp)$coefficients["fpcomp_intel", 1:2]*(-1), 2) # b, beta
round(summary(mr_dv3_fp)$coefficients["fpcomp_intel", 3], 2) # SE
round(summary(mr_dv3_fp)$coefficients["fpcomp_intel", 5], 3) # p-value
round(coef(mr_dv3_fp)["fpcomp_intel"]*(-1)*sqrt(1- cor_fp^2), 2) # r
# b = 2.36, beta = 0.21, SE = 0.61, p < 0.001, r = 0.16

# Regression coefficient: Confidence
round(summary(mr_dv3_fp)$coefficients["fpcomp_conf", 1:3], 2) # b, beta, SE
round(summary(mr_dv3_fp)$coefficients["fpcomp_conf", 5], 3) # p-value
round(coef(mr_dv3_fp)["fpcomp_conf"]*sqrt(1- cor_fp^2), 2) # r
# b = 0.98, beta = 0.10, SE = 0.53, p = 0.067, r = 0.08




################# DV 4: Situation selection at a distance (choice) ############

# *-- SEM ####

# *---- FP: Intelligence & Confidence ####

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
# Chi-squared(21) = 53.89, p < 0.001, CFI = 0.99, TLI = 0.98, RMSEA = 0.05

# Regression coefficient: Intelligence
# Note that sign of coefficient needs to be flipped
# because choice of "intelligent" website was coded as 0
par_dv4_fp <- parameterestimates(fit_dv4_fp)
round(par_dv4_fp[9, "est"]*1.7*(-1), 2) # b (probit to logit)
round(par_dv4_fp[9, "se"]*1.7, 2) # SE (probit to logit)
round(par_dv4_fp[9, "pvalue"], 3) # p-value
round(exp(par_dv4_fp[9, "est"]*1.7*(-1)), 2) # odds ratio
round(lor2r(par_dv4_fp[9, "est"]*1.7*(-1)), 2) # r
# b = 0.44, SE = 0.16, p = 0.005, OR = 1.56, r = 0.12

# Regression coefficient: Confidence
round(par_dv4_fp[10, c("est", "se")]*1.7, 2) # b, SE (probit to logit)
round(par_dv4_fp[10, "pvalue"], 3) # p-value
round(exp(par_dv4_fp[10, "est"]*1.7), 2) # odds ratio
round(lor2r(par_dv4_fp[10, "est"]*1.7), 2) # r
# b = 0.45, SE = 0.16, p = 0.004, OR = 1.57, r = 0.12

# Visualize model (note that the regression coefficients are in probit)
semPaths(fit_dv4_fp, "std", fade = F, esize = T, edge.color = 'black', 
         edge.label.cex = 1, sizeMan = 7, sizeLat = 10, rotation = 2,
         nodeLabels = c(lab_fp, "DV4", "FP_intel", "FP_conf"),
         label.norm = "OOOO", curvature = 3, intercepts = F, thresholds = F)



# *-- Bivariate Regression ####

# *---- FP: Intelligence ####
# Note that sign of coefficient needs to be flipped
# because choice of "intelligent" website was coded as 0
br_dv4_fp_intel <- glm(dv4_ic ~ fpcomp_intel, data = ss, family = "binomial")
round(summary(br_dv4_fp_intel)$coefficients["fpcomp_intel", 1]*(-1), 2) # b
round(summary(br_dv4_fp_intel)$coefficients["fpcomp_intel", 2], 2) # SE
round(summary(br_dv4_fp_intel)$coefficients["fpcomp_intel", 4], 3) # p-value
round(exp(coef(br_dv4_fp_intel)["fpcomp_intel"]*(-1)), 2) # odds ratio
round(lor2r(coef(br_dv4_fp_intel)["fpcomp_intel"]*(-1)), 2) # r
# b = 0.65, SE = 0.49, p = 0.190, OR = 1.91, r = 0.18


# *---- FP: Confidence ####
br_dv4_fp_conf <- glm(dv4_ic ~ fpcomp_conf, data = ss, family = "binomial")
round(summary(br_dv4_fp_conf)$coefficients["fpcomp_conf", 1:2], 2) # b, SE
round(summary(br_dv4_fp_conf)$coefficients["fpcomp_conf", 4], 3) # p-value
round(exp(coef(br_dv4_fp_conf)["fpcomp_conf"]), 2) # odds ratio
round(lor2r(coef(br_dv4_fp_conf)["fpcomp_conf"]), 2) # r
# b = 0.47, SE = 0.43, p = 0.277, OR = 1.60, r = 0.13



# *-- Multiple Regression ####

# *---- FP: Intelligence & Confidence ####
mr_dv4_fp <- glm(dv4_ic ~ fpcomp_intel + fpcomp_conf, 
                 data = ss, family = "binomial")

# Regression coefficient: Intelligence
# Note that sign of coefficient needs to be flipped
# because choice of "intelligent" website was coded as 0
round(summary(mr_dv4_fp)$coefficients["fpcomp_intel", 1]*(-1), 2) # b
round(summary(mr_dv4_fp)$coefficients["fpcomp_intel", 2], 2) # SE
round(summary(mr_dv4_fp)$coefficients["fpcomp_intel", 4], 3) # p-value
round(exp(coef(mr_dv4_fp)["fpcomp_intel"]*(-1)), 2) # odds ratio
round(lor2r(coef(mr_dv4_fp)["fpcomp_intel"]*(-1)), 2) # r
# b = 1.81, SE = 0.67, p = 0.007, OR = 6.11, r = 0.45

# Regression coefficient: Confidence
round(summary(mr_dv4_fp)$coefficients["fpcomp_conf", 1:2], 2) # b, SE
round(summary(mr_dv4_fp)$coefficients["fpcomp_conf", 4], 3) # p-value
round(exp(coef(mr_dv4_fp)["fpcomp_conf"]), 2) # odds ratio
round(lor2r(coef(mr_dv4_fp)["fpcomp_conf"]), 2) # r
# b = 1.52, SE = 0.59, p = 0.010, OR = 4.55, r = 0.39
