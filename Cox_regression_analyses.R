#--COX REGRESSION ANALYSES------------------------------------------------------
# Author: Pascal B. Meyre
# Date: 01/28/24
# Location: Reinach, Baselland, Switzerland

# Use cleaned dataset of SWISS.BEAT.biomarker.cleaned.csv
#-------------------------------------------------------------------------------

# Read the cleaned dataset
library(data.table)
dat <- fread("/Users/pascalmeyre/Desktop/Research/1_Projects_Analysis/18_Biomarkers_MACE_bleeding/analysis/datasets/SWISS.BEAT.biomarker.cleaned.imputed.csv")

# Calculate eGFR using Cockcroft-Gault formula
# dat$crep2 <- ifelse(dat$crep2 == 1, NA, dat$crep2)

# Define a function for Cockcroft-Gault equation
cockcroft_gault <- function(age, weight, serum_creatinine) {
  crcl <- (140 - age) * weight / (72 * serum_creatinine)
  return(crcl)
}

dat$creatinine_mg_dl <- dat$crep2/88.4

# Calculate CrCl using Cockcroft-Gault equation
dat$crcl_cg <- cockcroft_gault(age = dat$age.bl, weight = dat$gewicht, serum_creatinine = dat$creatinine_mg_dl)

# Print the eGFR values
print(dat$crcl_cg)

# List of biomarkers
biomarkers <- c("ang2", "ddi2h", "cysc", "alat", "gdf.15", "crphs", "igfbp7",
                "il6", "probnpii", "opn", "tnt.hs", "crcl_cg")

# Log-transform biomarkers
for (col in biomarkers) {
  dat[[paste0(col, "_log")]] <- log1p(dat[[col]])
}

# Standardize per 1-SD increase
biomarkers_log <- c("ang2_log", "ddi2h_log", "cysc_log", "alat_log", 
                    "gdf.15_log", "crphs_log", "igfbp7_log", "il6_log", "probnpii_log", 
                    "opn_log", "tnt.hs_log", "crcl_cg_log")

# Loop through each log-transformed biomarker
for (biomarker in biomarkers_log) {
  biomarker_std <- paste0(biomarker, "_std")
  
  # Compute mean and standard deviation
  biomarker_mean <- mean(dat[[biomarker]], na.rm = TRUE)
  biomarker_sd <- sd(dat[[biomarker]], na.rm = TRUE)
  
  # Standardize the log-transformed biomarker per 1-SD increase
  dat[[biomarker_std]] <- (dat[[biomarker]] - biomarker_mean) / biomarker_sd
}


################################################################################
# Composite outcome
################################################################################

#-------------------------------------------------------------------------------
# Check whether there are missing "time" variables 
num_missing_timeto1.death <- sum(is.na(dat$timeto1.death))
print(num_missing_timeto1.death)

num_missing_timeto1.isch.stroke <- sum(is.na(dat$timeto1.isch.stroke))
print(num_missing_timeto1.isch.stroke)

num_missing_timeto1.sys.embolism <- sum(is.na(dat$timeto1.sys.embolism))
print(num_missing_timeto1.sys.embolism)

num_missing_timeto1.mi <- sum(is.na(dat$timeto1.mi))
print(num_missing_timeto1.mi)
#-------------------------------------------------------------------------------

################### COX REGRESSION MODELS ######################################

# Create a single time and outcome variable for the composite outcome
# Calculate the shortest time variable for composite outcome
dat$timeto1.composite <- pmin(dat$timeto1.death, dat$timeto1.isch.stroke, dat$timeto1.sys.embolism, dat$timeto1.mi)

# Create a composite variable indicating whether an event occurred (1) or not (0)
dat$composite <- ifelse((dat$timeto1.composite == dat$timeto1.death) & dat$death.cardiac, 1,
                        ifelse((dat$timeto1.composite == dat$timeto1.isch.stroke) & dat$ischemic.stroke, 1, 
                               ifelse((dat$timeto1.composite == dat$timeto1.mi) & dat$mi, 1, 
                                      ifelse((dat$timeto1.composite == dat$timeto1.sys.embolism) & dat$sys.embolism, 1, 0))))

# Fit the age + sex adjusted Cox model
library(survival)

# Create empty lists to store results from age+sex and multivariable adjusted models
results.age_summary <- list()
results.multi_summary <- list()
results.age <- list()
results.multi <- list()

# Age+sex adjusted: Iterate over each biomarker
for (biomarker in biomarkers) {
  
  # Fit the Cox model
  cox_model <- coxph(Surv(timeto1.composite, composite) ~
                       get(paste0(biomarker, "_log_std")) + age.bl + pat.sex, data = dat)
  
  # Store results in the list
  results.age_summary[[biomarker]] <- summary(cox_model)
  
  # Store results in the list for volcano plots
  results.age[[biomarker]] <- coef(summary(cox_model))[, "exp(coef)"]
}

for (biomarker in biomarkers) {

  # Fit the Cox model
  cox_model <- coxph(Surv(timeto1.composite, composite) ~
                       get(paste0(biomarker, "_log_std")) + age.bl + pat.sex + bmi +
                       current.smoker + rr.sys.liegend + prev.diabetes +
                       prev.stroke.tia + prev.heart.failure + prev.niereninsuff +
                       coronary.heart.disease, data = dat)
  
  # Store results in the list
  results.multi_summary[[biomarker]] <- summary(cox_model)
  
  # Extract HR and p-value for _log_std
  coef_summary <- coef(summary(cox_model))
  hr <- coef_summary["get(paste0(biomarker, \"_log_std\"))", "exp(coef)"]
  p_value <- coef_summary["get(paste0(biomarker, \"_log_std\"))", "Pr(>|z|)"]
  
  # Store HR and p-value in the list for volcano plots
  results.multi[[biomarker]] <- c(HR = hr, P_Value = p_value)
}

# Combined model
cox_model <- coxph(Surv(timeto1.composite, composite) ~ ang2_log_std + ddi2h_log_std + 
                    cysc_log_std + alat_log_std + gdf.15_log_std + crphs_log_std + 
                     igfbp7_log_std + il6_log_std + probnpii_log_std + opn_log_std + tnt.hs_log_std + 
                     crcl_cg_log_std + age.bl + pat.sex + bmi + current.smoker + rr.sys.liegend + 
                     prev.diabetes + prev.stroke.tia + prev.heart.failure + crcl_cg_log_std +
                     prev.niereninsuff + coronary.heart.disease, data = dat)
summary(cox_model)

################### BACKWARD BIOMARKER SELECTION ###############################

# Perform a backward and forward selection of biomarkers
# Create a new dataset excluding patients with missing values in specified variables
# Exclude patients with missing variables
library(dplyr)

data <- dat %>% 
  select(ends_with("_log_std"), age.bl, pat.sex, bmi, current.smoker, rr.sys.liegend,
         prev.diabetes, prev.stroke.tia, prev.heart.failure, prev.niereninsuff,
         coronary.heart.disease, timeto1.composite, composite)

data <- na.omit(data)

# Fit the Cox proportional hazards model with the filtered dataset
full_model <- coxph(Surv(timeto1.composite, composite) ~ ang2_log_std + ddi2h_log_std + 
                      cysc_log_std + alat_log_std + gdf.15_log_std + crphs_log_std + 
                      igfbp7_log_std + il6_log_std + probnpii_log_std + opn_log_std + tnt.hs_log_std + 
                      crcl_cg_log_std + age.bl + pat.sex + bmi + current.smoker + rr.sys.liegend + 
                      prev.diabetes + prev.stroke.tia + prev.heart.failure +
                      prev.niereninsuff + coronary.heart.disease, data = data)

summary(full_model)

library(car)
vif(full_model)

# Define the range of models for stepwise selection (only consider changes in biomarkers)
scope <- list(
  lower = ~ age.bl + pat.sex + bmi + 
                current.smoker + rr.sys.liegend + prev.diabetes + prev.stroke.tia + 
                prev.heart.failure + prev.niereninsuff + coronary.heart.disease,
  upper = ~ ang2_log_std + ddi2h_log_std + 
                cysc_log_std + alat_log_std + gdf.15_log_std + crphs_log_std + 
                igfbp7_log_std + il6_log_std + probnpii_log_std + opn_log_std + tnt.hs_log_std + 
                crcl_cg_log_std + age.bl + pat.sex + bmi + 
                current.smoker + rr.sys.liegend + prev.diabetes + prev.stroke.tia + 
                prev.heart.failure + prev.niereninsuff + coronary.heart.disease
  )

# Perform stepwise selection for model refinement
reduced_model <- step(full_model, direction = "backward", scope = scope)

# Perform stepwise selection using stepAIC from MASS package
library(MASS) 
reduced_model <- stepAIC(full_model, scope = scope, direction = "backward")

# View the final selected model
summary(reduced_model)

# Exclude CRP
cox_model <- coxph(Surv(timeto1.composite, composite) ~ ang2_log_std + ddi2h_log_std + 
                     alat_log_std + gdf.15_log_std + il6_log_std + 
                     probnpii_log_std + tnt.hs_log_std + age.bl + pat.sex + bmi + 
                     current.smoker + rr.sys.liegend + prev.diabetes + 
                     prev.stroke.tia + prev.heart.failure + prev.niereninsuff + 
                     coronary.heart.disease, data = data)

summary(cox_model)

################ AUC FOR BASE MODEL AND BASE MODEL + BIOMARKERS ################

# Load necessary libraries
library(survival)
library(timeROC)
library(pROC)

# Fit the Cox proportional hazards model on the complete cases
cox_model_base <- coxph(Surv(timeto1.composite, composite) ~ age.bl + pat.sex + bmi + 
                     current.smoker + rr.sys.liegend + 
                     prev.diabetes + prev.stroke.tia + prev.heart.failure + 
                     prev.niereninsuff + coronary.heart.disease, data = data)

cox_model_biomarkers <- coxph(Surv(timeto1.composite, composite) ~ ang2_log_std + ddi2h_log_std + 
                                alat_log_std + gdf.15_log_std + il6_log_std + 
                                probnpii_log_std + tnt.hs_log_std + age.bl + pat.sex + bmi + 
                                current.smoker + rr.sys.liegend + prev.diabetes + 
                                prev.stroke.tia + prev.heart.failure + prev.niereninsuff + 
                                coronary.heart.disease, data = data)

# Obtain the risk scores (linear predictors) for each model
risk_scores_base <- predict(cox_model_base, type = "lp", newdata = data)
risk_scores_biomarkers <- predict(cox_model_biomarkers, type = "lp", newdata = data)

# Define time points at which to calculate the AUC
time_points <- c(1, 2, 3) * 365.25  # 1, 2, and 3 years in days

# Calculate time-dependent AUC for the base model
roc_results_base <- timeROC(T = data$timeto1.composite,
                            delta = data$composite,
                            marker = risk_scores_base,
                            cause = 1,
                            weighting = "marginal",
                            times = time_points,
                            ROC = TRUE)

# Calculate time-dependent AUC for the biomarker model
roc_results_biomarkers <- timeROC(T = data$timeto1.composite,
                                  delta = data$composite,
                                  marker = risk_scores_biomarkers,
                                  cause = 1,
                                  weighting = "marginal",
                                  times = time_points,
                                  ROC = TRUE)

# Print AUC values for both models
print(roc_results_base$AUC)
print(roc_results_biomarkers$AUC)

# Extract predicted probabilities at the specific time point (1 year, 365.25 days)
pred_probs_base <- predict(cox_model_base, type = "risk", newdata = data)
pred_probs_biomarkers <- predict(cox_model_biomarkers, type = "risk", newdata = data)

# Use pROC to create ROC curves at the 1-year mark
roc_base <- roc(data$composite, pred_probs_base, plot = FALSE, ci = TRUE, quiet = TRUE)
roc_biomarkers <- roc(data$composite, pred_probs_biomarkers, plot = FALSE, ci = TRUE, quiet = TRUE)

# Perform DeLong's test to compare the AUCs
delong_test <- roc.test(roc_base, roc_biomarkers)
print(delong_test)

# Plot ROC Curves at different time points
plot(roc_results, time = time_points[1], col = "blue", title = "ROC Curve at 1 Year")
plot(roc_results, time = time_points[2], col = "red", title = "ROC Curve at 2 Years", add = TRUE)
plot(roc_results, time = time_points[3], col = "green", title = "ROC Curve at 3 Years", add = TRUE)
legend("bottomright", legend = c("1 Year", "2 Years", "3 Years"), col = c("blue", "red", "green"), lwd = 2)

##############################  VOLCANO PLOTS ##################################

# Create an empty dataframe to store results
volcano_data <- data.frame(Biomarker = character(),
                           HR = numeric(),
                           P_Value = numeric(),
                           stringsAsFactors = FALSE)

# Iterate through each biomarker
for (biomarker in biomarkers) {
  
  # Fit the Cox model
  cox_model <- coxph(Surv(timeto1.composite, composite) ~
                       get(paste0(biomarker, "_log_std")) + age.bl + pat.sex + bmi +
                       current.smoker + rr.sys.liegend + prev.diabetes +
                       prev.stroke.tia + prev.heart.failure + prev.niereninsuff +
                       coronary.heart.disease, data = dat)
  
  # Extract HR and p-value for _log_std
  coef_summary <- coef(summary(cox_model))
  hr <- coef_summary["get(paste0(biomarker, \"_log_std\"))", "exp(coef)"]
  p_value <- coef_summary["get(paste0(biomarker, \"_log_std\"))", "Pr(>|z|)"]
  
  # Add results to the dataframe
  volcano_data <- rbind(volcano_data, data.frame(biomarker = biomarker, HR = hr, P_Value = p_value))
}

# Volcano plot from the combined model after stepwise biomarker selection
volcano_data <- data.frame(Biomarker = character(),
                           HR = numeric(),
                           P_Value = numeric(),
                           stringsAsFactors = FALSE)

biomarker_names <- c("ang2", "ddi2h", "alat", "gdf.15", 
                     "il6", "probnpii", "tnt.hs")

cox_model <- coxph(Surv(timeto1.composite, composite) ~ ang2_log_std + ddi2h_log_std + 
                     alat_log_std + gdf.15_log_std + il6_log_std + 
                     probnpii_log_std + tnt.hs_log_std + age.bl + pat.sex + bmi + 
                     current.smoker + rr.sys.liegend + prev.diabetes + 
                     prev.stroke.tia + prev.heart.failure + prev.niereninsuff + 
                     coronary.heart.disease, data = dat)

# Obtain summary of coefficients
coef_summary <- summary(cox_model)

hr <- coef_summary$coefficients[, "exp(coef)"][paste0(biomarker_names, "_log_std")]
p_value <- coef_summary$coefficients[, "Pr(>|z|)"][paste0(biomarker_names, "_log_std")]

# Combine results into a data frame
volcano_data <- rbind(volcano_data, data.frame(biomarker = biomarker_names, HR = hr, P_Value = p_value))

# Display the dataframe
print(volcano_data)

new_names <- c("ANG-2", "D-dimer", "ALAT", "GDF-15", "IL-6", "NT-proBNP", "hsTropT")

# Replace biomarker names in coef_data
volcano_data$biomarker <- ifelse(volcano_data$biomarker == "ang2", new_names[1],
                                 ifelse(volcano_data$biomarker == "ddi2h", new_names[2],
                                        ifelse(volcano_data$biomarker == "alat", new_names[3],
                                               ifelse(volcano_data$biomarker == "gdf.15", new_names[4],
                                                      ifelse(volcano_data$biomarker == "il6", new_names[5],
                                                             ifelse(volcano_data$biomarker == "probnpii", new_names[6],
                                                                    ifelse(volcano_data$biomarker == "tnt.hs", new_names[7],
                                                                           volcano_data$biomarker)))))))

library(ggplot2)
library(ggrepel)

ggplot(volcano_data, aes(x = HR, y = -log10(P_Value), label = biomarker)) +
  geom_point() +
  labs(x = "Standardized hazard ratio", y = "-Log10 p-value") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "red") +
  theme_bw() + 
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 9, color = "black"),
        axis.text.x = element_text(hjust = 0.5)) +
  coord_cartesian(xlim = c(0.5, 2.0)) +
  scale_x_continuous(trans = "log10", breaks = c(0.5, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2.0), labels = c(0.5, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2.0)) +
  geom_text_repel(size = 2.75)

########################### VARIABLE IMPORTANCE PLOTS ##########################

library(dplyr)

features <- dat %>% 
  select(ends_with("_log_std"), age.bl, pat.sex, bmi, current.smoker, rr.sys.liegend,
         prev.diabetes, prev.stroke.tia, prev.heart.failure, prev.niereninsuff,
         coronary.heart.disease, timeto1.composite, composite)

features <- na.omit(features)

cox_model <- coxph(Surv(timeto1.composite, composite) ~ ang2_log_std + ddi2h_log_std + 
                     alat_log_std + gdf.15_log_std + il6_log_std + 
                     probnpii_log_std + tnt.hs_log_std + age.bl + pat.sex + bmi + 
                     current.smoker + rr.sys.liegend + prev.diabetes + 
                     prev.stroke.tia + prev.heart.failure + prev.niereninsuff + 
                     coronary.heart.disease, data = features)
summary(cox_model)

# Extract coefficient estimates and standard errors
coef_est <- coef(cox_model)
coef_se <- sqrt(diag(vcov(cox_model)))

# Calculate partial chi-squared statistic for each predictor
partial_chisq <- (coef_est / coef_se)^2

# Extract the total degrees of freedom for the model
num_params <- length(coef(cox_model))
df_total <- num_params

# Adjust partial chi-squared statistic by subtracting total degrees of freedom
partial_chisq_adj <- partial_chisq - df_total

# Display the results
result <- data.frame(Variable = names(coef_est), Partial_ChiSq = partial_chisq)
print(result)

# Order the result by partial chi-squared statistic in descending order
result <- result[order(-result$Partial_ChiSq), ]

library(ggplot2)

# Define the desired label changes
label_changes <- c("ang2_log_std" = "ANG-2",
                   "ddi2h_log_std" = "D-dimer",
                   "alat_log_std" = "ALAT",
                   "gdf.15_log_std" = "GDF-15",
                   "il6_log_std" = "IL-6",
                   "probnpii_log_std" = "NT-proBNP",
                   "tnt.hs_log_std" = "hsTropT",
                   "age.bl" = "Age",
                   "pat.sex" = "Sex",
                   "bmi" = "BMI",
                   "current.smoker" = "Smoker",
                   "rr.sys.liegend" = "SBP",
                   "prev.diabetes" = "Diabetes",
                   "prev.stroke.tia" = "Prior stroke/TIA",
                   "prev.heart.failure" = "Heart failure",
                   "prev.niereninsuff" = "CKD",
                   "coronary.heart.disease" = "CAD")

# Create the dot plot
p <- ggplot(result, aes(x = Partial_ChiSq, y = reorder(Variable, Partial_ChiSq))) +
  geom_point(shape = 1, size = 2.5, color = "black") +  # Change shape to open circle and set color to black
  labs(x = "Partial χ2-df", size = 2, y = "", color = "black") +  # Fixed x-axis label and removed y-axis title
  scale_y_discrete(labels = label_changes) +  # Apply label changes to y-axis
  theme_bw() +  # Minimal theme
  theme(axis.text.x = element_text(size = 10),  # Adjust x-axis text size
        axis.text.y = element_text(size = 10),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(color = "black"),  # Set x-axis ticks color to black
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.major.y = element_line(linetype = "dotted", color = "grey")  # Specify major horizontal grid lines as dotted and grey
  )

print(p)

################################# FOREST PLOTS #################################

# Combined model
library(randomForest)
library(survminer)

cox_model <- coxph(Surv(timeto1.composite, composite) ~ ang2_log_std + ddi2h_log_std + 
                     alat_log_std + gdf.15_log_std + il6_log_std + 
                     probnpii_log_std + tnt.hs_log_std + age.bl + pat.sex + bmi + 
                     current.smoker + rr.sys.liegend + prev.diabetes + 
                     prev.stroke.tia + prev.heart.failure + prev.niereninsuff + 
                     coronary.heart.disease, data = dat)
summary(cox_model)

biomarkers <- c("ang2_log_std", "ddi2h_log_std", "alat_log_std", "gdf.15_log_std", "il6_log_std", "probnpii_log_std", "tnt.hs_log_std")

# Create a dataframe to store hazard ratios and 95% CIs
coef_data <- data.frame(
  Biomarker = biomarkers,
  HR = exp(coef(cox_model)[biomarkers]),             # Extracting hazard ratios
  Lower_CI = exp(confint(cox_model)[biomarkers, 1]),  # Extracting lower CIs
  Upper_CI = exp(confint(cox_model)[biomarkers, 2]),   # Extracting upper CIs
  p_value = summary(cox_model)$coefficients[biomarkers, "Pr(>|z|)"]  # Extracting p-values
)

# View the dataframe
print(coef_data)

new_names <- c("ANG-2", "D-dimer", "ALAT", "GDF-15", "IL-6", "NT-proBNP", "hsTropT")

# Replace biomarker names in coef_data
coef_data$Biomarker <- ifelse(coef_data$Biomarker == "ang2_log_std", new_names[1],
                              ifelse(coef_data$Biomarker == "ddi2h_log_std", new_names[2],
                                     ifelse(coef_data$Biomarker == "alat_log_std", new_names[3],
                                            ifelse(coef_data$Biomarker == "gdf.15_log_std", new_names[4],
                                                   ifelse(coef_data$Biomarker == "il6_log_std", new_names[5],
                                                          ifelse(coef_data$Biomarker == "probnpii_log_std", new_names[6],
                                                                 ifelse(coef_data$Biomarker == "tnt.hs_log_std", new_names[7],
                                                                        coef_data$Biomarker)))))))

coef_data$p.signif <- ifelse(coef_data$p_value < 0.05, "*", "")

library(ggplot2)

# Create coefficient plot using ggplot2
coefficient_plot <- ggplot(coef_data, aes(x = Biomarker, y = HR, color = Biomarker)) +
  geom_point(size = 2.5) + 
  geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 0.3, linetype = "solid") + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") + 
  labs(x = "Biomarker", y = "Hazard ratio (95% CI) per 1 SD increase in biomarker levels", size=11, color = "black") + 
  theme_bw() +
  theme(legend.position="none",
        axis.title = element_text(face="bold"),
        axis.text.x = element_text(size=9, color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(size=9, color = "black")) +
  scale_y_continuous(limits = c(0.8, 1.6), breaks = seq(0.8, 1.6, by = 0.1)) +
  geom_text(aes(label = paste0(sprintf("%.2f", HR))), 
            hjust = -0.3, size = 3)
print(coefficient_plot)


################################################################################
# Heart failure hospitalization
################################################################################

#-------------------------------------------------------------------------------
# Check whether there are missing "time" variables 
num_missing_timeto1.hf <- sum(is.na(dat$timeto1.hf))
print(num_missing_timeto1.hf)
#-------------------------------------------------------------------------------

############################ COX REGRESSION MODELS #############################

# Fit the age + sex adjusted Cox model
library(survival)

# Create empty lists to store results from age+sex and multivariable adjusted models
results.age_summary <- list()
results.multi_summary <- list()
results.age <- list()
results.multi <- list()

# Age+sex adjusted: Iterate over each biomarker
for (biomarker in biomarkers) {
  
  # Fit the Cox model
  cox_model <- coxph(Surv(timeto1.hf, heart.failure) ~
                       get(paste0(biomarker, "_log_std")) + age.bl + pat.sex, data = dat)
  
  # Store results in the list
  results.age_summary[[biomarker]] <- summary(cox_model)
  
  # Store results in the list for volcano plots
  results.age[[biomarker]] <- coef(summary(cox_model))[, "exp(coef)"]
}

for (biomarker in biomarkers) {
  
  # Fit the Cox model
  cox_model <- coxph(Surv(timeto1.hf, heart.failure) ~
                       get(paste0(biomarker, "_log_std")) + age.bl + pat.sex + bmi +
                       current.smoker + rr.sys.liegend + prev.diabetes +
                       prev.stroke.tia + prev.heart.failure + prev.niereninsuff +
                       coronary.heart.disease, data = dat)
  
  # Store results in the list
  results.multi_summary[[biomarker]] <- summary(cox_model)
  
  # Extract HR and p-value for _log_std
  coef_summary <- coef(summary(cox_model))
  hr <- coef_summary["get(paste0(biomarker, \"_log_std\"))", "exp(coef)"]
  p_value <- coef_summary["get(paste0(biomarker, \"_log_std\"))", "Pr(>|z|)"]
  
  # Store HR and p-value in the list for volcano plots
  results.multi[[biomarker]] <- c(HR = hr, P_Value = p_value)
}

# Combined model 1 (only significant biomarkers)
cox_model <- coxph(Surv(timeto1.hf, heart.failure) ~ ang2_log_std + ddi2h_log_std + 
                     crcl_cg_log_std + cysc_log_std + gdf.15_log_std + crphs_log_std + 
                     igfbp7_log_std + il6_log_std + probnpii_log_std + opn_log_std + tnt.hs_log_std + 
                     age.bl + pat.sex + bmi + current.smoker + rr.sys.liegend + 
                     prev.diabetes + prev.stroke.tia + prev.heart.failure + 
                     prev.niereninsuff + coronary.heart.disease, data = dat)
summary(cox_model)

# Combined model 2 (all biomarkers)
cox_model <- coxph(Surv(timeto1.hf, heart.failure) ~ ang2_log_std + ddi2h_log_std + 
                     crcl_cg_log_std + cysc_log_std + alat_log_std + gdf.15_log_std + crphs_log_std + 
                     igfbp7_log_std + il6_log_std + probnpii_log_std + opn_log_std + tnt.hs_log_std + 
                     age.bl + pat.sex + bmi + current.smoker + rr.sys.liegend + 
                     prev.diabetes + prev.stroke.tia + prev.heart.failure + 
                     prev.niereninsuff + coronary.heart.disease, data = dat)
summary(cox_model)

################### BACKWARD BIOMARKER SELECTION ###############################

# Perform a backward and forward selection of biomarkers
# Create a new dataset excluding patients with missing values in specified variables
# Exclude patients with missing variables
library(dplyr)

data <- dat %>% 
  select(ends_with("_log_std"), age.bl, pat.sex, bmi, current.smoker, rr.sys.liegend,
         prev.diabetes, prev.stroke.tia, prev.heart.failure, prev.niereninsuff,
         coronary.heart.disease, timeto1.hf, heart.failure)

data <- na.omit(data)

# Fit the Cox proportional hazards model with the filtered dataset
full_model <- coxph(Surv(timeto1.hf, heart.failure) ~ ang2_log_std + ddi2h_log_std + 
                      cysc_log_std + alat_log_std + gdf.15_log_std + crphs_log_std + 
                      igfbp7_log_std + il6_log_std + probnpii_log_std + opn_log_std + tnt.hs_log_std + 
                      crcl_cg_log_std + age.bl + pat.sex + bmi + current.smoker + rr.sys.liegend + 
                      prev.diabetes + prev.stroke.tia + prev.heart.failure + crcl_cg_log_std +
                      prev.niereninsuff + coronary.heart.disease, 
                    data = data)

summary(full_model)

library(car)
vif(full_model)

# Define the range of models for stepwise selection (only consider changes in biomarkers)
scope <- list(lower = ~ age.bl + pat.sex + bmi + 
                current.smoker + rr.sys.liegend + prev.diabetes + prev.stroke.tia + 
                prev.heart.failure + prev.niereninsuff + coronary.heart.disease,
              upper = ~ ang2_log_std + ddi2h_log_std + 
                cysc_log_std + alat_log_std + gdf.15_log_std + crphs_log_std + 
                igfbp7_log_std + il6_log_std + probnpii_log_std + opn_log_std + tnt.hs_log_std + 
                crcl_cg_log_std + age.bl + pat.sex + bmi + 
                current.smoker + rr.sys.liegend + prev.diabetes + prev.stroke.tia + 
                prev.heart.failure + prev.niereninsuff + coronary.heart.disease)

# Perform stepwise selection for model refinement
reduced_model <- step(full_model, direction = "backward", scope = scope)

# Perform stepwise selection using stepAIC from MASS package
library(MASS) 
reduced_model <- stepAIC(full_model, scope = scope, direction = "backward")

# View the final selected model
summary(reduced_model)

cox_model <- coxph(Surv(timeto1.hf, heart.failure) ~ alat_log_std + gdf.15_log_std + 
                     igfbp7_log_std + probnpii_log_std + tnt.hs_log_std + crcl_cg_log_std +
                     age.bl + pat.sex + bmi + current.smoker + rr.sys.liegend + 
                     prev.diabetes + prev.stroke.tia + prev.heart.failure + 
                     prev.niereninsuff + coronary.heart.disease, data = data)

summary(cox_model)

################ AUC FOR BASE MODEL AND BASE MODEL + BIOMARKERS ################

# Load necessary libraries
library(survival)
library(timeROC)
library(pROC)

# Fit the Cox proportional hazards model on the complete cases
cox_model_base <- coxph(Surv(timeto1.hf, heart.failure) ~ age.bl + pat.sex + bmi + 
                          current.smoker + rr.sys.liegend + 
                          prev.diabetes + prev.stroke.tia + prev.heart.failure + 
                          prev.niereninsuff + coronary.heart.disease, data = data)

cox_model_biomarkers <- coxph(Surv(timeto1.hf, heart.failure) ~ alat_log_std + gdf.15_log_std + 
                                igfbp7_log_std + probnpii_log_std + tnt.hs_log_std + crcl_cg_log_std +
                                age.bl + pat.sex + bmi + current.smoker + rr.sys.liegend + 
                                prev.diabetes + prev.stroke.tia + prev.heart.failure + 
                                prev.niereninsuff + coronary.heart.disease, data = data)

# Obtain the risk scores (linear predictors) for each model
risk_scores_base <- predict(cox_model_base, type = "lp", newdata = data)
risk_scores_biomarkers <- predict(cox_model_biomarkers, type = "lp", newdata = data)

# Define time points at which to calculate the AUC
time_points <- c(1, 2, 3) * 365.25  # 1, 2, and 3 years in days

# Calculate time-dependent AUC for the base model
roc_results_base <- timeROC(T = data$timeto1.hf,
                            delta = data$heart.failure,
                            marker = risk_scores_base,
                            cause = 1,
                            weighting = "marginal",
                            times = time_points,
                            ROC = TRUE)

# Calculate time-dependent AUC for the biomarker model
roc_results_biomarkers <- timeROC(T = data$timeto1.hf,
                                  delta = data$heart.failure,
                                  marker = risk_scores_biomarkers,
                                  cause = 1,
                                  weighting = "marginal",
                                  times = time_points,
                                  ROC = TRUE)

# Print AUC values for both models
print(roc_results_base$AUC)
print(roc_results_biomarkers$AUC)

# Extract predicted probabilities at the specific time point (1 year, 365.25 days)
pred_probs_base <- predict(cox_model_base, type = "risk", newdata = data)
pred_probs_biomarkers <- predict(cox_model_biomarkers, type = "risk", newdata = data)

# Use pROC to create ROC curves at the 1-year mark
roc_base <- roc(data$heart.failure, pred_probs_base, plot = FALSE, ci = TRUE, quiet = TRUE)
roc_biomarkers <- roc(data$heart.failure, pred_probs_biomarkers, plot = FALSE, ci = TRUE, quiet = TRUE)

# Perform DeLong's test to compare the AUCs
delong_test <- roc.test(roc_base, roc_biomarkers)
print(delong_test)

# Plot ROC Curves at different time points
plot(roc_results_biomarkers, time = time_points[1], col = "blue", title = "ROC Curve at 1 Year")
plot(roc_results_biomarkers, time = time_points[2], col = "red", title = "ROC Curve at 2 Years", add = TRUE)
plot(roc_results_biomarkers, time = time_points[3], col = "green", title = "ROC Curve at 3 Years", add = TRUE)
legend("bottomright", legend = c("1 Year", "2 Years", "3 Years"), col = c("blue", "red", "green"), lwd = 2)

############################### VOLCANO PLOTS ##################################

# Create an empty dataframe to store results
volcano_data <- data.frame(Biomarker = character(),
                           HR = numeric(),
                           P_Value = numeric(),
                           stringsAsFactors = FALSE)

# Iterate through each biomarker
for (biomarker in biomarkers) {
  
  # Fit the Cox model
  cox_model <- coxph(Surv(timeto1.hf, heart.failure) ~
                       get(paste0(biomarker, "_log_std")) + age.bl + pat.sex + bmi +
                       current.smoker + rr.sys.liegend + prev.diabetes +
                       prev.stroke.tia + prev.heart.failure + prev.niereninsuff +
                       coronary.heart.disease, data = dat)
  
  # Extract HR and p-value for _log_std
  coef_summary <- coef(summary(cox_model))
  hr <- coef_summary["get(paste0(biomarker, \"_log_std\"))", "exp(coef)"]
  p_value <- coef_summary["get(paste0(biomarker, \"_log_std\"))", "Pr(>|z|)"]
  
  # Add results to the dataframe
  volcano_data <- rbind(volcano_data, data.frame(biomarker = biomarker, HR = hr, P_Value = p_value))
}

# Volcano plot from the combined model after stepwise biomarker selection
volcano_data <- data.frame(Biomarker = character(),
                           HR = numeric(),
                           P_Value = numeric(),
                           stringsAsFactors = FALSE)

biomarker_names <- c("alat", "gdf.15", "igfbp7", "probnpii", "tnt.hs")

cox_model <- coxph(Surv(timeto1.hf, heart.failure) ~ alat_log_std + gdf.15_log_std + 
                     igfbp7_log_std + probnpii_log_std + tnt.hs_log_std + crcl_cg_log_std +
                     age.bl + pat.sex + bmi + current.smoker + rr.sys.liegend + 
                     prev.diabetes + prev.stroke.tia + prev.heart.failure + 
                     prev.niereninsuff + coronary.heart.disease, data = dat)

# Obtain summary of coefficients
coef_summary <- summary(cox_model)

hr <- coef_summary$coefficients[, "exp(coef)"][paste0(biomarker_names, "_log_std")]
p_value <- coef_summary$coefficients[, "Pr(>|z|)"][paste0(biomarker_names, "_log_std")]

# Combine results into a data frame
volcano_data <- rbind(volcano_data, data.frame(biomarker = biomarker_names, HR = hr, P_Value = p_value))

# Display the dataframe
print(volcano_data)

new_names <- c("ALAT", "GDF-15", "IGFBP-7", "NT-proBNP", "hsTropT")

# Replace biomarker names in coef_data
volcano_data$biomarker <- ifelse(volcano_data$biomarker == "alat", new_names[1],
                                 ifelse(volcano_data$biomarker == "gdf.15", new_names[2],
                                        ifelse(volcano_data$biomarker == "igfbp7", new_names[3],
                                               ifelse(volcano_data$biomarker == "probnpii", new_names[4],
                                                      ifelse(volcano_data$biomarker == "tnt.hs", new_names[5],
                                                             volcano_data$biomarker)))))

library(ggplot2)
library(ggrepel)

ggplot(volcano_data, aes(x = HR, y = -log10(P_Value), label = biomarker)) +
  geom_point() +
  labs(x = "Standardized hazard ratio", y = "-Log10 P-value") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "red") +
  theme_bw() + 
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 9, color = "black"),
        axis.text.x = element_text(hjust = 0.5)) +
  coord_cartesian(xlim = c(0.5, 2.0)) +
  scale_x_continuous(trans = "log10", breaks = c(0.5, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2.0), labels = c(0.5, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2.0)) +
  geom_text_repel(size = 2.75)

########################## VARIABLE IMPORTANCE PLOTS ###########################

library(dplyr)

features <- dat %>% 
  select(ends_with("_log_std"), age.bl, pat.sex, bmi, current.smoker, rr.sys.liegend,
         prev.diabetes, prev.stroke.tia, prev.heart.failure, prev.niereninsuff,
         coronary.heart.disease, timeto1.hf, heart.failure)

features <- na.omit(features)

cox_model <- coxph(Surv(timeto1.hf, heart.failure) ~ alat_log_std + gdf.15_log_std + 
                     igfbp7_log_std + probnpii_log_std + tnt.hs_log_std + crcl_cg_log_std +
                     age.bl + pat.sex + bmi + current.smoker + rr.sys.liegend + 
                     prev.diabetes + prev.stroke.tia + prev.heart.failure + 
                     prev.niereninsuff + coronary.heart.disease, data = features)
summary(cox_model)

# Extract coefficient estimates and standard errors
coef_est <- coef(cox_model)
coef_se <- sqrt(diag(vcov(cox_model)))

# Calculate partial chi-squared statistic for each predictor
partial_chisq <- (coef_est / coef_se)^2

# Extract the total degrees of freedom for the model
num_params <- length(coef(cox_model))
df_total <- num_params

# Adjust partial chi-squared statistic by subtracting total degrees of freedom
partial_chisq_adj <- partial_chisq - df_total

# Display the results
result <- data.frame(Variable = names(coef_est), Partial_ChiSq = partial_chisq)
print(result)

# Order the result by partial chi-squared statistic in descending order
result <- result[order(-result$Partial_ChiSq), ]

library(ggplot2)

# Define the desired label changes
label_changes <- c("alat_log_std" = "ALAT",
                   "gdf.15_log_std" = "GDF-15",
                   "igfbp7_log_std" = "IGFBP-7",
                   "probnpii_log_std" = "NT-proBNP",
                   "crcl_cg_log_std" = "eGFR",
                   "tnt.hs_log_std" = "hsTropT",
                   "age.bl" = "Age",
                   "pat.sex" = "Sex",
                   "bmi" = "BMI",
                   "current.smoker" = "Smoker",
                   "rr.sys.liegend" = "SBP",
                   "prev.diabetes" = "Diabetes",
                   "prev.stroke.tia" = "Prior stroke/TIA",
                   "prev.heart.failure" = "Heart failure",
                   "prev.niereninsuff" = "CKD",
                   "coronary.heart.disease" = "CAD")

# Create the dot plot
p <- ggplot(result, aes(x = Partial_ChiSq, y = reorder(Variable, Partial_ChiSq))) +
  geom_point(shape = 1, size = 2.5, color = "black") +  # Change shape to open circle and set color to black
  labs(x = "Partial χ2-df", size = 2, y = "", color = "black") +  # Fixed x-axis label and removed y-axis title
  scale_y_discrete(labels = label_changes) +  # Apply label changes to y-axis
  theme_bw() +  # Minimal theme
  theme(axis.text.x = element_text(size = 10),  # Adjust x-axis text size
        axis.text.y = element_text(size = 10),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(color = "black"),  # Set x-axis ticks color to black
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.major.y = element_line(linetype = "dotted", color = "grey")  # Specify major horizontal grid lines as dotted and grey
  )

print(p)

################################ FOREST PLOTS ##################################

# Combined model
library(randomForest)
library(survminer)

cox_model <- coxph(Surv(timeto1.hf, heart.failure) ~ alat_log_std + gdf.15_log_std + 
                     igfbp7_log_std + probnpii_log_std + tnt.hs_log_std + crcl_cg_log_std +
                     age.bl + pat.sex + bmi + current.smoker + rr.sys.liegend + 
                     prev.diabetes + prev.stroke.tia + prev.heart.failure + 
                     prev.niereninsuff + coronary.heart.disease, data = dat)
summary(cox_model)

biomarkers <- c("alat_log_std", "gdf.15_log_std", "igfbp7_log_std", "probnpii_log_std", "tnt.hs_log_std")

# Create a dataframe to store hazard ratios and 95% CIs
coef_data <- data.frame(
  Biomarker = biomarkers,
  HR = exp(coef(cox_model)[biomarkers]),             # Extracting hazard ratios
  Lower_CI = exp(confint(cox_model)[biomarkers, 1]),  # Extracting lower CIs
  Upper_CI = exp(confint(cox_model)[biomarkers, 2]),   # Extracting upper CIs
  p_value = summary(cox_model)$coefficients[biomarkers, "Pr(>|z|)"]  # Extracting p-values
)

# View the dataframe
print(coef_data)

new_names <- c("ALAT", "GDF-15", "IGFBP-7", "NT-proBNP", "hsTropT")

# Replace biomarker names in coef_data
coef_data$Biomarker <- ifelse(coef_data$Biomarker == "alat_log_std", new_names[1],
                              ifelse(coef_data$Biomarker == "gdf.15_log_std", new_names[2],
                                     ifelse(coef_data$Biomarker == "igfbp7_log_std", new_names[3],
                                            ifelse(coef_data$Biomarker == "probnpii_log_std", new_names[4],
                                                   ifelse(coef_data$Biomarker == "tnt.hs_log_std", new_names[5],
                                                          coef_data$Biomarker)))))

coef_data$p.signif <- ifelse(coef_data$p_value < 0.05, "*", "")

library(ggplot2)

# Create coefficient plot using ggplot2
coefficient_plot <- ggplot(coef_data, aes(x = Biomarker, y = HR, color = Biomarker)) +
  geom_point(size = 2.5) + 
  geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 0.3, linetype = "solid") + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") + 
  labs(x = "Biomarker", y = "Hazard ratio (95% CI) per 1 SD increase in biomarker levels", size=11, color = "black") + 
  theme_bw() +
  theme(legend.position="none",
        axis.title = element_text(face="bold"),
        axis.text.x = element_text(size=9, color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(size=9, color = "black")) +
  scale_y_continuous(limits = c(0.8, 1.8), breaks = seq(0.8, 1.8, by = 0.1)) +
  geom_text(aes(label = paste0(sprintf("%.2f", HR))), 
            hjust = -0.3, size = 3)
print(coefficient_plot)


################################################################################
# Major bleeding
################################################################################

#-------------------------------------------------------------------------------
# Check whether there are missing "time" variables 
num_missing_timeto1.major.bleed <- sum(is.na(dat$timeto1.major.bleed))
print(num_missing_timeto1.major.bleed)
#-------------------------------------------------------------------------------

############################## COX REGRESSION MODELS ###########################

# Fit the age + sex adjusted Cox model
library(survival)

# Create empty lists to store results from age+sex and multivariable adjusted models
results.age_summary <- list()
results.multi_summary <- list()
results.age <- list()
results.multi <- list()

# Age+sex adjusted: Iterate over each biomarker
for (biomarker in biomarkers) {
  
  # Fit the Cox model
  cox_model <- coxph(Surv(timeto1.major.bleed, major.bleed) ~
                       get(paste0(biomarker, "_log_std")) + age.bl + pat.sex, data = dat)
  
  # Store results in the list
  results.age_summary[[biomarker]] <- summary(cox_model)
  
  # Store results in the list for volcano plots
  results.age[[biomarker]] <- coef(summary(cox_model))[, "exp(coef)"]
}

# Multivariable: Iterate over each biomarker
for (biomarker in biomarkers) {
  
  # Fit the Cox model
  cox_model <- coxph(Surv(timeto1.major.bleed, major.bleed) ~
                       get(paste0(biomarker, "_log_std")) + age.bl + pat.sex + bmi +
                       current.smoker + rr.sys.liegend + prev.diabetes +
                       prev.stroke.tia + prev.heart.failure + prev.niereninsuff +
                       coronary.heart.disease, data = dat)
  
  # Store results in the list
  results.multi_summary[[biomarker]] <- summary(cox_model)
  
  # Store results in the list for volcano plots
  results.multi[[biomarker]] <- coef(summary(cox_model))[, "exp(coef)"]
}

# Combined model 1 (only significant biomarkers)
cox_model <- coxph(Surv(timeto1.major.bleed, major.bleed) ~ ang2_log_std + 
                     cysc_log_std + gdf.15_log_std + igfbp7_log_std + 
                     il6_log_std + probnpii_log_std + opn_log_std + tnt.hs_log_std + 
                     age.bl + pat.sex + bmi + current.smoker + rr.sys.liegend + 
                     prev.diabetes + prev.stroke.tia + prev.heart.failure + 
                     prev.niereninsuff + coronary.heart.disease, data = dat)
summary(cox_model)

# Combined model 2 (all biomarkers)
cox_model <- coxph(Surv(timeto1.major.bleed, major.bleed) ~ ang2_log_std + ddi2h_log_std + 
                     crcl_cg_log_std + cysc_log_std + alat_log_std + gdf.15_log_std + crphs_log_std + 
                     igfbp7_log_std + il6_log_std + probnpii_log_std + opn_log_std + tnt.hs_log_std + 
                     age.bl + pat.sex + bmi + current.smoker + rr.sys.liegend + 
                     prev.diabetes + prev.stroke.tia + prev.heart.failure + 
                     prev.niereninsuff + coronary.heart.disease, data = dat)
summary(cox_model)

######################### BACKWARD BIOMARKER SELECTION #########################

# Perform a backward and forward selection of biomarkers
# Create a new dataset excluding patients with missing values in specified variables
# Exclude patients with missing variables
library(dplyr)

data <- dat %>% 
  select(ends_with("_log_std"), age.bl, pat.sex, bmi, current.smoker, rr.sys.liegend,
         prev.diabetes, prev.stroke.tia, prev.heart.failure, prev.niereninsuff,
         coronary.heart.disease, timeto1.major.bleed, major.bleed)

data <- na.omit(data)

# Fit the Cox proportional hazards model with the filtered dataset
full_model <- coxph(Surv(timeto1.major.bleed, major.bleed) ~ ang2_log_std + ddi2h_log_std + 
                      crcl_cg_log_std + cysc_log_std + alat_log_std + gdf.15_log_std + crphs_log_std + 
                      igfbp7_log_std + il6_log_std + probnpii_log_std + opn_log_std + tnt.hs_log_std + 
                      age.bl + pat.sex + bmi + current.smoker + rr.sys.liegend + 
                      prev.diabetes + prev.stroke.tia + prev.heart.failure + 
                      prev.niereninsuff + coronary.heart.disease, data = data)

summary(full_model)

# Define the range of models for stepwise selection (only consider changes in biomarkers)
scope <- list(lower = ~ age.bl + pat.sex + bmi + 
                current.smoker + rr.sys.liegend + prev.diabetes + prev.stroke.tia + 
                prev.heart.failure + prev.niereninsuff + coronary.heart.disease,
              upper = ~ ang2_log_std + ddi2h_log_std + 
                crcl_cg_log_std + cysc_log_std + alat_log_std + gdf.15_log_std + 
                crphs_log_std + igfbp7_log_std + il6_log_std + probnpii_log_std + 
                opn_log_std + tnt.hs_log_std + age.bl + pat.sex + bmi + 
                current.smoker + rr.sys.liegend + prev.diabetes + prev.stroke.tia + 
                prev.heart.failure + prev.niereninsuff + coronary.heart.disease)

# Perform stepwise selection for model refinement
reduced_model <- step(full_model, direction = "backward", scope = scope)

# Perform stepwise selection using stepAIC from MASS package
library(MASS) 
reduced_model <- stepAIC(full_model, scope = scope, direction = "backward")

# View the final selected model
summary(reduced_model)

# Use the biomarkers included in the step_model_filtered for the combined model
cox_model <- coxph(Surv(timeto1.major.bleed, major.bleed) ~ crcl_cg_log_std + gdf.15_log_std + 
                     igfbp7_log_std + il6_log_std + tnt.hs_log_std + 
                     age.bl + pat.sex + bmi + current.smoker + rr.sys.liegend + 
                     prev.diabetes + prev.stroke.tia + prev.heart.failure + 
                     prev.niereninsuff + coronary.heart.disease, data = dat)
summary(cox_model)

################ AUC FOR BASE MODEL AND BASE MODEL + BIOMARKERS ################

# Load necessary libraries
library(survival)
library(timeROC)
library(pROC)

# Fit the Cox proportional hazards model on the complete cases
cox_model_base <- coxph(Surv(timeto1.major.bleed, major.bleed) ~ age.bl + pat.sex + bmi + 
                          current.smoker + rr.sys.liegend + 
                          prev.diabetes + prev.stroke.tia + prev.heart.failure + 
                          prev.niereninsuff + coronary.heart.disease, data = data)

cox_model_biomarkers <- coxph(Surv(timeto1.major.bleed, major.bleed) ~ crcl_cg_log_std + gdf.15_log_std + 
                                igfbp7_log_std + il6_log_std + tnt.hs_log_std + 
                                age.bl + pat.sex + bmi + current.smoker + rr.sys.liegend + 
                                prev.diabetes + prev.stroke.tia + prev.heart.failure + 
                                prev.niereninsuff + coronary.heart.disease, data = data)

# Obtain the risk scores (linear predictors) for each model
risk_scores_base <- predict(cox_model_base, type = "lp", newdata = data)
risk_scores_biomarkers <- predict(cox_model_biomarkers, type = "lp", newdata = data)

# Define time points at which to calculate the AUC
time_points <- c(1, 2, 3) * 365.25  # 1, 2, and 3 years in days

# Calculate time-dependent AUC for the base model
roc_results_base <- timeROC(T = data$timeto1.major.bleed,
                            delta = data$major.bleed,
                            marker = risk_scores_base,
                            cause = 1,
                            weighting = "marginal",
                            times = time_points,
                            ROC = TRUE)

# Calculate time-dependent AUC for the biomarker model
roc_results_biomarkers <- timeROC(T = data$timeto1.major.bleed,
                                  delta = data$major.bleed,
                                  marker = risk_scores_biomarkers,
                                  cause = 1,
                                  weighting = "marginal",
                                  times = time_points,
                                  ROC = TRUE)

# Print AUC values for both models
print(roc_results_base$AUC)
print(roc_results_biomarkers$AUC)

# Extract predicted probabilities at the specific time point (1 year, 365.25 days)
pred_probs_base <- predict(cox_model_base, type = "risk", newdata = data)
pred_probs_biomarkers <- predict(cox_model_biomarkers, type = "risk", newdata = data)

# Use pROC to create ROC curves at the 1-year mark
roc_base <- roc(data$major.bleed, pred_probs_base, plot = FALSE, ci = TRUE, quiet = TRUE)
roc_biomarkers <- roc(data$major.bleed, pred_probs_biomarkers, plot = FALSE, ci = TRUE, quiet = TRUE)

# Perform DeLong's test to compare the AUCs
delong_test <- roc.test(roc_base, roc_biomarkers)
print(delong_test)

# Plot ROC Curves at different time points
plot(roc_biomarkers, time = time_points[1], col = "blue", title = "ROC Curve at 1 Year")
plot(roc_biomarkers, time = time_points[2], col = "red", title = "ROC Curve at 2 Years", add = TRUE)
plot(roc_biomarkers, time = time_points[3], col = "green", title = "ROC Curve at 3 Years", add = TRUE)
legend("bottomright", legend = c("1 Year", "2 Years", "3 Years"), col = c("blue", "red", "green"), lwd = 2)

############################### VOLCANO PLOTS ##################################

# Create an empty dataframe to store results
volcano_data <- data.frame(Biomarker = character(),
                           HR = numeric(),
                           P_Value = numeric(),
                           stringsAsFactors = FALSE)

# Iterate through each biomarker
for (biomarker in biomarkers) {
  
  # Fit the Cox model
  cox_model <- coxph(Surv(timeto1.major.bleed, major.bleed) ~
                       get(paste0(biomarker, "_log_std")) + age.bl + pat.sex + bmi +
                       current.smoker + rr.sys.liegend + prev.diabetes +
                       prev.stroke.tia + prev.heart.failure + prev.niereninsuff +
                       coronary.heart.disease, data = dat)
  
  # Extract HR and p-value for _log_std
  coef_summary <- coef(summary(cox_model))
  hr <- coef_summary["get(paste0(biomarker, \"_log_std\"))", "exp(coef)"]
  p_value <- coef_summary["get(paste0(biomarker, \"_log_std\"))", "Pr(>|z|)"]
  
  # Add results to the dataframe
  volcano_data <- rbind(volcano_data, data.frame(biomarker = biomarker, HR = hr, P_Value = p_value))
}

# Volcano plot from the combined model after stepwise biomarker selection
volcano_data <- data.frame(Biomarker = character(),
                           HR = numeric(),
                           P_Value = numeric(),
                           stringsAsFactors = FALSE)

biomarker_names <- c("gdf.15", "igfbp7", "il6", "tnt.hs")

cox_model <- coxph(Surv(timeto1.major.bleed, major.bleed) ~ crcl_cg_log_std + gdf.15_log_std + 
                     igfbp7_log_std + il6_log_std + tnt.hs_log_std + 
                     age.bl + pat.sex + bmi + current.smoker + rr.sys.liegend + 
                     prev.diabetes + prev.stroke.tia + prev.heart.failure + 
                     prev.niereninsuff + coronary.heart.disease, data = dat)

# Obtain summary of coefficients
coef_summary <- summary(cox_model)

hr <- coef_summary$coefficients[, "exp(coef)"][paste0(biomarker_names, "_log_std")]
p_value <- coef_summary$coefficients[, "Pr(>|z|)"][paste0(biomarker_names, "_log_std")]

# Combine results into a data frame
volcano_data <- rbind(volcano_data, data.frame(biomarker = biomarker_names, HR = hr, P_Value = p_value))

# Display the dataframe
print(volcano_data)

new_names <- c("GDF-15", "IGFBP-7", "IL-6", "hsTropT")

# Replace biomarker names in coef_data
volcano_data$biomarker <- ifelse(volcano_data$biomarker == "gdf.15", new_names[1],
                                 ifelse(volcano_data$biomarker == "igfbp7", new_names[2],
                                        ifelse(volcano_data$biomarker == "il6", new_names[3],
                                               ifelse(volcano_data$biomarker == "tnt.hs", new_names[4],
                                                      volcano_data$biomarker))))

library(ggplot2)
library(ggrepel)

ggplot(volcano_data, aes(x = HR, y = -log10(P_Value), label = biomarker)) +
  geom_point() +
  labs(x = "Standardized hazard ratio", y = "-Log10 P-value") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "red") +
  theme_bw() + 
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 9, color = "black"),
        axis.text.x = element_text(hjust = 0.5)) +
  coord_cartesian(xlim = c(0.5, 2.0)) +
  scale_x_continuous(trans = "log10", breaks = c(0.5, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2.0), labels = c(0.5, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2.0)) +
  geom_text_repel(size = 2.75)

########################## VARIABLE IMPORTANCE PLOTS ###########################

library(dplyr)

features <- dat %>% 
  select(ends_with("_log_std"), age.bl, pat.sex, bmi, current.smoker, rr.sys.liegend,
         prev.diabetes, prev.stroke.tia, prev.heart.failure, prev.niereninsuff,
         coronary.heart.disease, timeto1.major.bleed, major.bleed)

features <- na.omit(features)

cox_model <- coxph(Surv(timeto1.major.bleed, major.bleed) ~ crcl_cg_log_std + gdf.15_log_std + 
                     igfbp7_log_std + il6_log_std + tnt.hs_log_std + 
                     age.bl + pat.sex + bmi + current.smoker + rr.sys.liegend + 
                     prev.diabetes + prev.stroke.tia + prev.heart.failure + 
                     prev.niereninsuff + coronary.heart.disease, data = features)
summary(cox_model)

# Extract coefficient estimates and standard errors
coef_est <- coef(cox_model)
coef_se <- sqrt(diag(vcov(cox_model)))

# Calculate partial chi-squared statistic for each predictor
partial_chisq <- (coef_est / coef_se)^2

# Extract the total degrees of freedom for the model
num_params <- length(coef(cox_model))
df_total <- num_params

# Adjust partial chi-squared statistic by subtracting total degrees of freedom
partial_chisq_adj <- partial_chisq - df_total

# Display the results
result <- data.frame(Variable = names(coef_est), Partial_ChiSq = partial_chisq)
print(result)

# Order the result by partial chi-squared statistic in descending order
result <- result[order(-result$Partial_ChiSq), ]

library(ggplot2)

# Define the desired label changes
label_changes <- c("crcl_cg_log_std" = "eGFR",
                   "gdf.15_log_std" = "GDF-15",
                   "igfbp7_log_std" = "IGFBP-7",
                   "il6_log_std" = "IL-6",
                   "tnt.hs_log_std" = "hsTropT",
                   "age.bl" = "Age",
                   "pat.sex" = "Sex",
                   "bmi" = "BMI",
                   "current.smoker" = "Smoker",
                   "rr.sys.liegend" = "SBP",
                   "prev.diabetes" = "Diabetes",
                   "prev.stroke.tia" = "Prior stroke/TIA",
                   "prev.heart.failure" = "Heart failure",
                   "prev.niereninsuff" = "CKD",
                   "coronary.heart.disease" = "CAD")

# Create the dot plot
p <- ggplot(result, aes(x = Partial_ChiSq, y = reorder(Variable, Partial_ChiSq))) +
  geom_point(shape = 1, size = 2.5, color = "black") +  # Change shape to open circle and set color to black
  labs(x = "Partial χ2-df", size = 2, y = "", color = "black") +  # Fixed x-axis label and removed y-axis title
  scale_y_discrete(labels = label_changes) +  # Apply label changes to y-axis
  theme_bw() +  # Minimal theme
  theme(axis.text.x = element_text(size = 10),  # Adjust x-axis text size
        axis.text.y = element_text(size = 10),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(color = "black"),  # Set x-axis ticks color to black
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.major.y = element_line(linetype = "dotted", color = "grey")  # Specify major horizontal grid lines as dotted and grey
  )

print(p)

################################ FOREST PLOTS ##################################

# Combined model
library(randomForest)
library(survminer)

cox_model <- coxph(Surv(timeto1.major.bleed, major.bleed) ~ crcl_cg_log_std + gdf.15_log_std + 
                     igfbp7_log_std + il6_log_std + tnt.hs_log_std + 
                     age.bl + pat.sex + bmi + current.smoker + rr.sys.liegend + 
                     prev.diabetes + prev.stroke.tia + prev.heart.failure + 
                     prev.niereninsuff + coronary.heart.disease, data = dat)
summary(cox_model)

biomarkers <- c("gdf.15_log_std", "igfbp7_log_std", 
                "il6_log_std", "tnt.hs_log_std")

# Create a dataframe to store hazard ratios and 95% CIs
coef_data <- data.frame(
  Biomarker = biomarkers,
  HR = exp(coef(cox_model)[biomarkers]),             # Extracting hazard ratios
  Lower_CI = exp(confint(cox_model)[biomarkers, 1]),  # Extracting lower CIs
  Upper_CI = exp(confint(cox_model)[biomarkers, 2]),   # Extracting upper CIs
  p_value = summary(cox_model)$coefficients[biomarkers, "Pr(>|z|)"]  # Extracting p-values
)

# View the dataframe
print(coef_data)

new_names <- c("GDF-15", "IGFBP-7", "IL-6", "hsTropT")

# Replace biomarker names in coef_data
coef_data$Biomarker <- ifelse(coef_data$Biomarker == "gdf.15_log_std", new_names[1],
                              ifelse(coef_data$Biomarker == "igfbp7_log_std", new_names[2],
                                     ifelse(coef_data$Biomarker == "il6_log_std", new_names[3],
                                            ifelse(coef_data$Biomarker == "tnt.hs_log_std", new_names[4],
                                                   coef_data$Biomarker))))

coef_data$p.signif <- ifelse(coef_data$p_value < 0.05, "*", "")

library(ggplot2)

# Create coefficient plot using ggplot2
coefficient_plot <- ggplot(coef_data, aes(x = Biomarker, y = HR, color = Biomarker)) +
  geom_point(size = 2.5) + 
  geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 0.3, linetype = "solid") + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") + 
  labs(x = "Biomarker", y = "Hazard ratio (95% CI) per 1 SD increase in biomarker levels", size=11, color = "black") + 
  theme_bw() +
  theme(legend.position="none",
        axis.title = element_text(face="bold"),
        axis.text.x = element_text(size=9, color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(size=9, color = "black")) +
  scale_y_continuous(limits = c(0.8, 1.6), breaks = seq(0.8, 1.6, by = 0.1)) +
  geom_text(aes(label = paste0(sprintf("%.2f", HR))), 
            hjust = -0.3, size = 3)
print(coefficient_plot)


################################################################################
# Ischemic stroke
################################################################################

#-------------------------------------------------------------------------------
# Check whether there are missing "time" variables 
num_missing_timeto1.isch.stroke <- sum(is.na(dat$timeto1.isch.stroke))
print(num_missing_timeto1.isch.stroke)
#-------------------------------------------------------------------------------

############################ COX REGRESSION MODELS #############################

# Fit the age + sex adjusted Cox model
library(survival)

# Create empty lists to store results from age+sex and multivariable adjusted models
results.age_summary <- list()
results.multi_summary <- list()
results.age <- list()
results.multi <- list()

# Age+sex adjusted: Iterate over each biomarker
for (biomarker in biomarkers) {
  
  # Fit the Cox model
  cox_model <- coxph(Surv(timeto1.isch.stroke, ischemic.stroke) ~
                       get(paste0(biomarker, "_log_std")) + age.bl + pat.sex, data = dat)
  
  # Store results in the list
  results.age_summary[[biomarker]] <- summary(cox_model)
  
  # Store results in the list for volcano plots
  results.age[[biomarker]] <- coef(summary(cox_model))[, "exp(coef)"]
}

for (biomarker in biomarkers) {
  
  # Fit the Cox model
  cox_model <- coxph(Surv(timeto1.isch.stroke, ischemic.stroke) ~
                       get(paste0(biomarker, "_log_std")) + age.bl + pat.sex + bmi +
                       current.smoker + rr.sys.liegend + prev.diabetes +
                       prev.stroke.tia + prev.heart.failure + prev.niereninsuff +
                       coronary.heart.disease, data = dat)
  
  # Store results in the list
  results.multi_summary[[biomarker]] <- summary(cox_model)
  
  # Extract HR and p-value for _log_std
  coef_summary <- coef(summary(cox_model))
  hr <- coef_summary["get(paste0(biomarker, \"_log_std\"))", "exp(coef)"]
  p_value <- coef_summary["get(paste0(biomarker, \"_log_std\"))", "Pr(>|z|)"]
  
  # Store HR and p-value in the list for volcano plots
  results.multi[[biomarker]] <- c(HR = hr, P_Value = p_value)
}

# Combined model 2 (all biomarkers)
cox_model <- coxph(Surv(timeto1.isch.stroke, ischemic.stroke) ~ ang2_log_std + ddi2h_log_std + 
                     crcl_cg_log_std + cysc_log_std + alat_log_std + gdf.15_log_std + crphs_log_std + 
                     igfbp7_log_std + il6_log_std + probnpii_log_std + opn_log_std + tnt.hs_log_std + 
                     age.bl + pat.sex + bmi + current.smoker + rr.sys.liegend + 
                     prev.diabetes + prev.stroke.tia + prev.heart.failure + 
                     prev.niereninsuff + coronary.heart.disease, data = dat)
summary(cox_model)

################### BACKWARD BIOMARKER SELECTION ###############################

# Perform a backward and forward selection of biomarkers
# Create a new dataset excluding patients with missing values in specified variables
# Exclude patients with missing variables
library(dplyr)

data <- dat %>% 
  select(ends_with("_log_std"), age.bl, pat.sex, bmi, current.smoker, rr.sys.liegend,
         prev.diabetes, prev.stroke.tia, prev.heart.failure, prev.niereninsuff,
         coronary.heart.disease, timeto1.isch.stroke, ischemic.stroke)

data <- na.omit(data)

# Fit the Cox proportional hazards model with the filtered dataset
full_model <- coxph(Surv(timeto1.isch.stroke, ischemic.stroke) ~ ang2_log_std + ddi2h_log_std + 
                      crcl_cg_log_std + cysc_log_std + alat_log_std + gdf.15_log_std + 
                      crphs_log_std + igfbp7_log_std + il6_log_std + probnpii_log_std + 
                      opn_log_std + tnt.hs_log_std + age.bl + pat.sex + bmi + 
                      current.smoker + rr.sys.liegend + prev.diabetes + prev.stroke.tia + 
                      prev.heart.failure + prev.niereninsuff + coronary.heart.disease, 
                    data = data)

summary(full_model)

# Define the range of models for stepwise selection (only consider changes in biomarkers)
scope <- list(lower = ~ age.bl + pat.sex + bmi + 
                current.smoker + rr.sys.liegend + prev.diabetes + prev.stroke.tia + 
                prev.heart.failure + prev.niereninsuff + coronary.heart.disease,
              upper = ~ ang2_log_std + ddi2h_log_std + 
                crcl_cg_log_std + cysc_log_std + alat_log_std + gdf.15_log_std + 
                crphs_log_std + igfbp7_log_std + il6_log_std + probnpii_log_std + 
                opn_log_std + tnt.hs_log_std + age.bl + pat.sex + bmi + 
                current.smoker + rr.sys.liegend + prev.diabetes + prev.stroke.tia + 
                prev.heart.failure + prev.niereninsuff + coronary.heart.disease)

# Perform stepwise selection for model refinement
reduced_model <- step(full_model, direction = "backward", scope = scope)

# Perform stepwise selection using stepAIC from MASS package
library(MASS) 
reduced_model <- stepAIC(full_model, scope = scope, direction = "backward")

# View the final selected model
summary(reduced_model)

# Use the biomarkers included in the step_model_filtered for the combined model
cox_model <- coxph(Surv(timeto1.isch.stroke, ischemic.stroke) ~ alat_log_std + 
                     probnpii_log_std + opn_log_std +
                     age.bl + pat.sex + bmi + current.smoker + rr.sys.liegend + 
                     prev.diabetes + prev.stroke.tia + prev.heart.failure + 
                     prev.niereninsuff + coronary.heart.disease, data = dat)
summary(cox_model)

########## BACKWARD BIOMARKER SELECTION USING CHADSVASC AS BASE MODEL ##########

# Perform a backward and forward selection of biomarkers
# Create a new dataset excluding patients with missing values in specified variables
# Exclude patients with missing variables
library(dplyr)

data <- dat %>% 
  select(ends_with("_log"), timeto1.isch.stroke, ischemic.stroke, chads)

data <- na.omit(data)

# Fit the Cox proportional hazards model with the filtered dataset
full_model <- coxph(Surv(timeto1.isch.stroke, ischemic.stroke) ~ ang2_log + ddi2h_log + 
                      crcl_cg_log + cysc_log + alat_log + gdf.15_log + 
                      crphs_log + igfbp7_log + il6_log + probnpii_log + 
                      opn_log + tnt.hs_log + chads, 
                    data = data)

summary(full_model)

# Define the range of models for stepwise selection (only consider changes in biomarkers)
scope <- list(lower = ~ chads,
              upper = ~ ang2_log + ddi2h_log + 
                crcl_cg_log + cysc_log + alat_log + gdf.15_log + 
                crphs_log + igfbp7_log + il6_log + probnpii_log + 
                opn_log + tnt.hs_log + chads)

# Perform stepwise selection for model refinement
reduced_model <- step(full_model, direction = "backward", scope = scope)

# Perform stepwise selection using stepAIC from MASS package
library(MASS) 
reduced_model <- stepAIC(full_model, scope = scope, direction = "backward")

# View the final selected model
summary(reduced_model)

# Use the biomarkers included in the step_model_filtered for the combined model
cox_model <- coxph(Surv(timeto1.isch.stroke, ischemic.stroke) ~ probnpii_log + chads, data = dat)
summary(cox_model)

cox_zph <- cox.zph(cox_model)
print(cox_zph)
plot(cox_zph)

################ AUC FOR BASE MODEL AND BASE MODEL + BIOMARKERS ################

# Load necessary libraries
library(survival)
library(timeROC)
library(pROC)

# Fit the Cox proportional hazards model on the complete cases
cox_model_base <- coxph(Surv(timeto1.isch.stroke, ischemic.stroke) ~ age.bl + pat.sex + bmi + 
                          current.smoker + rr.sys.liegend + 
                          prev.diabetes + prev.stroke.tia + prev.heart.failure + 
                          prev.niereninsuff + coronary.heart.disease, data = data)

cox_model_biomarkers <- coxph(Surv(timeto1.isch.stroke, ischemic.stroke) ~ alat_log_std + 
                                probnpii_log_std + opn_log_std +
                                age.bl + pat.sex + bmi + current.smoker + rr.sys.liegend + 
                                prev.diabetes + prev.stroke.tia + prev.heart.failure + 
                                prev.niereninsuff + coronary.heart.disease, data = data)

# Obtain the risk scores (linear predictors) for each model
risk_scores_base <- predict(cox_model_base, type = "lp", newdata = data)
risk_scores_biomarkers <- predict(cox_model_biomarkers, type = "lp", newdata = data)

# Define time points at which to calculate the AUC
time_points <- c(1, 2, 3) * 365.25  # 1, 2, and 3 years in days

# Calculate time-dependent AUC for the base model
roc_results_base <- timeROC(T = data$timeto1.isch.stroke,
                            delta = data$ischemic.stroke,
                            marker = risk_scores_base,
                            cause = 1,
                            weighting = "marginal",
                            times = time_points,
                            ROC = TRUE)

# Calculate time-dependent AUC for the biomarker model
roc_results_biomarkers <- timeROC(T = data$timeto1.isch.stroke,
                                  delta = data$ischemic.stroke,
                                  marker = risk_scores_biomarkers,
                                  cause = 1,
                                  weighting = "marginal",
                                  times = time_points,
                                  ROC = TRUE)

# Print AUC values for both models
print(roc_results_base$AUC)
print(roc_results_biomarkers$AUC)

# Extract predicted probabilities at the specific time point (1 year, 365.25 days)
pred_probs_base <- predict(cox_model_base, type = "risk", newdata = data)
pred_probs_biomarkers <- predict(cox_model_biomarkers, type = "risk", newdata = data)

# Use pROC to create ROC curves at the 1-year mark
roc_base <- roc(data$ischemic.stroke, pred_probs_base, plot = FALSE, ci = TRUE, quiet = TRUE)
roc_biomarkers <- roc(data$ischemic.stroke, pred_probs_biomarkers, plot = FALSE, ci = TRUE, quiet = TRUE)

# Perform DeLong's test to compare the AUCs
delong_test <- roc.test(roc_base, roc_biomarkers)
print(delong_test)

# Plot ROC Curves at different time points
plot(roc_results_biomarkers, time = time_points[1], col = "blue", title = "ROC Curve at 1 Year")
plot(roc_results_biomarkers, time = time_points[2], col = "red", title = "ROC Curve at 2 Years", add = TRUE)
plot(roc_results_biomarkers, time = time_points[3], col = "green", title = "ROC Curve at 3 Years", add = TRUE)
legend("bottomright", legend = c("1 Year", "2 Years", "3 Years"), col = c("blue", "red", "green"), lwd = 2)

####################### AUC of CHADS-VASC + biomarkers #########################
library(survival)
library(timeROC)
library(pROC)

# Fit the Cox proportional hazards model on the complete cases
cox_model_base <- coxph(Surv(timeto1.isch.stroke, ischemic.stroke) ~ chads, data = dat)
cox_model_biomarkers <- coxph(Surv(timeto1.isch.stroke, ischemic.stroke) ~ probnpii_log_std + chads, data = dat)

# Obtain the risk scores (linear predictors) for each model
risk_scores_base <- predict(cox_model_base, type = "lp", newdata = dat)
risk_scores_biomarkers <- predict(cox_model_biomarkers, type = "lp", newdata = dat)

time_points <- c(1, 2, 3, 4, 5, 6) * 365.25  # 1 to 6 years in days

# Calculate time-dependent AUC for the base model
roc_results_base <- timeROC(T = dat$timeto1.isch.stroke,
                            delta = dat$ischemic.stroke,
                            marker = risk_scores_base,
                            cause = 1,
                            weighting = "marginal",
                            times = time_points,
                            ROC = TRUE)

# Calculate time-dependent AUC for the biomarker model
roc_results_biomarkers <- timeROC(T = dat$timeto1.isch.stroke,
                                  delta = dat$ischemic.stroke,
                                  marker = risk_scores_biomarkers,
                                  cause = 1,
                                  weighting = "marginal",
                                  times = time_points,
                                  ROC = TRUE)

# Print AUC values for both models
print(roc_results_base$AUC)
print(roc_results_biomarkers$AUC)

# Extract predicted probabilities using the entire follow-up duration
pred_probs_base <- predict(cox_model_base, type = "risk", newdata = dat)
pred_probs_biomarkers <- predict(cox_model_biomarkers, type = "risk", newdata = dat)

# Use pROC to create ROC curves (without specific time point argument)
roc_base <- roc(dat$ischemic.stroke, pred_probs_base, plot = FALSE, ci = TRUE, quiet = TRUE)
roc_biomarkers <- roc(dat$ischemic.stroke, pred_probs_biomarkers, plot = FALSE, ci = TRUE, quiet = TRUE)

# Perform DeLong's test to compare the AUCs
delong_test <- roc.test(roc_base, roc_biomarkers)
print(delong_test)

# Extract ROC curve data for 1-year mark (365.25 days)
roc_data_base <- data.frame(
  FPR = 1 - roc_results_base$FP[, 1],
  TPR = roc_results_base$TP[, 1],
  model = "CHA2DS2-VASc Score"
)

roc_data_biomarkers <- data.frame(
  FPR = 1 - roc_results_biomarkers$FP[, 1],
  TPR = roc_results_biomarkers$TP[, 1],
  model = "CHA2DS2-VASc Score + biomarkers"
)

plot(roc_results_base, time = time_points[1], col = "blue", title = "CHA2DS2-VASc Score")
plot(roc_results_biomarkers, time = time_points[1], col = "red", title = "CHA2DS2-VASc Score + biomarkers", add = TRUE)

# Add a legend
legend("topright", 
       legend = c("CHA2DS2-VASc Score", "CHA2DS2-VASc Score + biomarkers"), 
       col = c("blue", "red"), 
       lwd = 2,
       inset = c(-0.4, -0.3),  # adjust the inset to move the legend further inside
       xpd = TRUE)


















############################### VOLCANO PLOTS ##################################

# Create an empty dataframe to store results
volcano_data <- data.frame(Biomarker = character(),
                           HR = numeric(),
                           P_Value = numeric(),
                           stringsAsFactors = FALSE)

# Iterate through each biomarker
for (biomarker in biomarkers) {
  
  # Fit the Cox model
  cox_model <- coxph(Surv(timeto1.isch.stroke, ischemic.stroke) ~
                       get(paste0(biomarker, "_log_std")) + age.bl + pat.sex + bmi +
                       current.smoker + rr.sys.liegend + prev.diabetes +
                       prev.stroke.tia + prev.heart.failure + prev.niereninsuff +
                       coronary.heart.disease, data = dat)
  
  # Extract HR and p-value for _log_std
  coef_summary <- coef(summary(cox_model))
  hr <- coef_summary["get(paste0(biomarker, \"_log_std\"))", "exp(coef)"]
  p_value <- coef_summary["get(paste0(biomarker, \"_log_std\"))", "Pr(>|z|)"]
  
  # Add results to the dataframe
  volcano_data <- rbind(volcano_data, data.frame(biomarker = biomarker, HR = hr, P_Value = p_value))
}

# Volcano plot from the combined model after stepwise biomarker selection
volcano_data <- data.frame(Biomarker = character(),
                           HR = numeric(),
                           P_Value = numeric(),
                           stringsAsFactors = FALSE)

biomarker_names <- c("alat", "probnpii", "opn")

cox_model <- coxph(Surv(timeto1.isch.stroke, ischemic.stroke) ~ alat_log_std + 
                     probnpii_log_std + opn_log_std + age.bl + pat.sex + bmi + 
                     current.smoker + rr.sys.liegend + 
                     prev.diabetes + prev.stroke.tia + prev.heart.failure + 
                     prev.niereninsuff + coronary.heart.disease, data = dat)

# Obtain summary of coefficients
coef_summary <- summary(cox_model)

hr <- coef_summary$coefficients[, "exp(coef)"][paste0(biomarker_names, "_log_std")]
p_value <- coef_summary$coefficients[, "Pr(>|z|)"][paste0(biomarker_names, "_log_std")]

# Combine results into a data frame
volcano_data <- rbind(volcano_data, data.frame(biomarker = biomarker_names, HR = hr, P_Value = p_value))

# Display the dataframe
print(volcano_data)

new_names <- c("ALAT", "NT-proBNP", "OPN")

# Replace biomarker names in coef_data
volcano_data$biomarker <- ifelse(volcano_data$biomarker == "alat", new_names[1],
                                 ifelse(volcano_data$biomarker == "probnpii", new_names[2],
                                        ifelse(volcano_data$biomarker == "opn", new_names[3],
                                               volcano_data$biomarker)))

library(ggplot2)
library(ggrepel)

ggplot(volcano_data, aes(x = HR, y = -log10(P_Value), label = biomarker)) +
  geom_point() +
  labs(x = "Standardized hazard ratio", y = "-Log10 P-value") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "red") +
  theme_bw() + 
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 9, color = "black"),
        axis.text.x = element_text(hjust = 0.5)) +
  coord_cartesian(xlim = c(0.5, 2.0)) +
  scale_x_continuous(trans = "log10", breaks = c(0.5, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2.0), labels = c(0.5, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2.0)) +
  geom_text_repel(size = 2.75)

########################## VARIABLE IMPORTANCE PLOTS ###########################

library(dplyr)

features <- dat %>% 
  select(ends_with("_log_std"), age.bl, pat.sex, bmi, current.smoker, rr.sys.liegend,
         prev.diabetes, prev.stroke.tia, prev.heart.failure, prev.niereninsuff,
         coronary.heart.disease, timeto1.isch.stroke, ischemic.stroke)

features <- na.omit(features)

cox_model <- coxph(Surv(timeto1.isch.stroke, ischemic.stroke) ~ alat_log_std + 
                     probnpii_log_std + opn_log_std + age.bl + pat.sex + bmi + 
                     current.smoker + rr.sys.liegend + 
                     prev.diabetes + prev.stroke.tia + prev.heart.failure + 
                     prev.niereninsuff + coronary.heart.disease, data = features)
summary(cox_model)

# Extract coefficient estimates and standard errors
coef_est <- coef(cox_model)
coef_se <- sqrt(diag(vcov(cox_model)))

# Calculate partial chi-squared statistic for each predictor
partial_chisq <- (coef_est / coef_se)^2

# Extract the total degrees of freedom for the model
num_params <- length(coef(cox_model))
df_total <- num_params

# Adjust partial chi-squared statistic by subtracting total degrees of freedom
partial_chisq_adj <- partial_chisq - df_total

# Display the results
result <- data.frame(Variable = names(coef_est), Partial_ChiSq = partial_chisq)
print(result)

# Order the result by partial chi-squared statistic in descending order
result <- result[order(-result$Partial_ChiSq), ]

library(ggplot2)

# Define the desired label changes
label_changes <- c("alat_log_std" = "ALAT",
                   "probnpii_log_std" = "NT-proBNP",
                   "opn_log_std" = "OPN",
                   "age.bl" = "Age",
                   "pat.sex" = "Sex",
                   "bmi" = "BMI",
                   "current.smoker" = "Smoker",
                   "rr.sys.liegend" = "SBP",
                   "prev.diabetes" = "Diabetes",
                   "prev.stroke.tia" = "Prior stroke/TIA",
                   "prev.heart.failure" = "Heart failure",
                   "prev.niereninsuff" = "CKD",
                   "coronary.heart.disease" = "CAD")

# Create the dot plot
p <- ggplot(result, aes(x = Partial_ChiSq, y = reorder(Variable, Partial_ChiSq))) +
  geom_point(shape = 1, size = 2.5, color = "black") +  # Change shape to open circle and set color to black
  labs(x = "Partial χ2-df", size = 2, y = "", color = "black") +  # Fixed x-axis label and removed y-axis title
  scale_y_discrete(labels = label_changes) +  # Apply label changes to y-axis
  theme_bw() +  # Minimal theme
  theme(axis.text.x = element_text(size = 10),  # Adjust x-axis text size
        axis.text.y = element_text(size = 10),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(color = "black"),  # Set x-axis ticks color to black
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.major.y = element_line(linetype = "dotted", color = "grey")  # Specify major horizontal grid lines as dotted and grey
  )

print(p)

################################ FOREST PLOTS ##################################

# Combined model
library(randomForest)
library(survminer)

cox_model <- coxph(Surv(timeto1.isch.stroke, ischemic.stroke) ~ alat_log_std + 
                     probnpii_log_std + opn_log_std + age.bl + pat.sex + bmi + 
                     current.smoker + rr.sys.liegend + 
                     prev.diabetes + prev.stroke.tia + prev.heart.failure + 
                     prev.niereninsuff + coronary.heart.disease, data = dat)
summary(cox_model)

biomarkers <- c("alat_log_std", "probnpii_log_std", "opn_log_std")

# Create a dataframe to store hazard ratios and 95% CIs
coef_data <- data.frame(
  Biomarker = biomarkers,
  HR = exp(coef(cox_model)[biomarkers]),             # Extracting hazard ratios
  Lower_CI = exp(confint(cox_model)[biomarkers, 1]),  # Extracting lower CIs
  Upper_CI = exp(confint(cox_model)[biomarkers, 2]),   # Extracting upper CIs
  p_value = summary(cox_model)$coefficients[biomarkers, "Pr(>|z|)"]  # Extracting p-values
)

# View the dataframe
print(coef_data)

new_names <- c("ALAT", "NT-proBNP", "OPN")

# Replace biomarker names in coef_data
coef_data$Biomarker <- ifelse(coef_data$Biomarker == "alat_log_std", new_names[1],
                              ifelse(coef_data$Biomarker == "probnpii_log_std", new_names[2],
                                     ifelse(coef_data$Biomarker == "opn_log_std", new_names[3], 
                                            coef_data$Biomarker)))

coef_data$p.signif <- ifelse(coef_data$p_value < 0.05, "*", "")

library(ggplot2)

# Create coefficient plot using ggplot2
coefficient_plot <- ggplot(coef_data, aes(x = Biomarker, y = HR, color = Biomarker)) +
  geom_point(size = 2.5) + 
  geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 0.3, linetype = "solid") + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") + 
  labs(x = "Biomarker", y = "Hazard ratio (95% CI) per 1 SD increase in biomarker levels", size=11, color = "black") + 
  theme_bw() +
  theme(legend.position="none",
        axis.title = element_text(face="bold"),
        axis.text.x = element_text(size=9, color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(size=9, color = "black")) +
  scale_y_continuous(limits = c(0.6, 2.0), breaks = seq(0.6, 2.0, by = 0.1)) +
  geom_text(aes(label = paste0(sprintf("%.2f", HR))), 
            hjust = -0.3, size = 3)
print(coefficient_plot)


################################################################################
# All strokes
################################################################################

#-------------------------------------------------------------------------------
# Check whether there are missing "time" variables 
num_missing_timeto1.stroke <- sum(is.na(dat$timeto1.stroke))
print(num_missing_timeto1.stroke)
#-------------------------------------------------------------------------------

############################ COX REGRESSION MODELS #############################

# Fit the age + sex adjusted Cox model
library(survival)

# Create empty lists to store results from age+sex and multivariable adjusted models
results.age_summary <- list()
results.multi_summary <- list()
results.age <- list()
results.multi <- list()

# Age+sex adjusted: Iterate over each biomarker
for (biomarker in biomarkers) {
  
  # Fit the Cox model
  cox_model <- coxph(Surv(timeto1.stroke, stroke) ~
                       get(paste0(biomarker, "_log_std")) + age.bl + pat.sex, data = dat)
  
  # Store results in the list
  results.age_summary[[biomarker]] <- summary(cox_model)
  
  # Store results in the list for volcano plots
  results.age[[biomarker]] <- coef(summary(cox_model))[, "exp(coef)"]
}

# Multivariable: Iterate over each biomarker
for (biomarker in biomarkers) {
  
  # Fit the Cox model
  cox_model <- coxph(Surv(timeto1.stroke, stroke) ~
                       get(paste0(biomarker, "_log_std")) + age.bl + pat.sex + bmi +
                       current.smoker + rr.sys.liegend + prev.diabetes +
                       prev.stroke.tia + prev.heart.failure + prev.niereninsuff +
                       coronary.heart.disease, data = dat)
  
  # Store results in the list
  results.multi_summary[[biomarker]] <- summary(cox_model)
  
  # Store results in the list for volcano plots
  results.multi[[biomarker]] <- coef(summary(cox_model))[, "exp(coef)"]
}

# Combined model 1 (only significant biomarkers)
cox_model <- coxph(Surv(timeto1.stroke, stroke) ~ ang2_log_std + il6_log_std + 
                     probnpii_log_std + tnt.hs_log_std + 
                     age.bl + pat.sex + bmi + current.smoker + rr.sys.liegend + 
                     prev.diabetes + prev.stroke.tia + prev.heart.failure + 
                     prev.niereninsuff + coronary.heart.disease, data = dat)
summary(cox_model)

# Combined model 2 (all biomarkers)
cox_model <- coxph(Surv(timeto1.stroke, stroke) ~ ang2_log_std + ddi2h_log_std + 
                     crcl_cg_log_std + cysc_log_std + alat_log_std + gdf.15_log_std + crphs_log_std + 
                     igfbp7_log_std + il6_log_std + probnpii_log_std + opn_log_std + tnt.hs_log_std + 
                     age.bl + pat.sex + bmi + current.smoker + rr.sys.liegend + 
                     prev.diabetes + prev.stroke.tia + prev.heart.failure + 
                     prev.niereninsuff + coronary.heart.disease, data = dat)
summary(cox_model)

################### BACKWARD BIOMARKER SELECTION ###############################

# Perform a backward and forward selection of biomarkers
# Create a new dataset excluding patients with missing values in specified variables
# Exclude patients with missing variables
library(dplyr)

data <- dat %>% 
  select(ends_with("_log_std"), age.bl, pat.sex, bmi, current.smoker, rr.sys.liegend,
         prev.diabetes, prev.stroke.tia, prev.heart.failure, prev.niereninsuff,
         coronary.heart.disease, timeto1.stroke, stroke, chads)

data <- na.omit(data)

library(dplyr)

data <- dat %>% 
  select(ends_with("_log"), age.bl, pat.sex, bmi, current.smoker, rr.sys.liegend,
         prev.diabetes, prev.stroke.tia, prev.heart.failure, prev.niereninsuff,
         coronary.heart.disease, timeto1.stroke, stroke, chads)

data <- na.omit(data)

# Fit the Cox proportional hazards model with the filtered dataset
full_model <- coxph(Surv(timeto1.stroke, stroke) ~ ang2_log_std + ddi2h_log_std + 
                      crcl_cg_log_std + cysc_log_std + alat_log_std + gdf.15_log_std + 
                      crphs_log_std + igfbp7_log_std + il6_log_std + probnpii_log_std + 
                      opn_log_std + tnt.hs_log_std + age.bl + pat.sex + bmi + 
                      current.smoker + rr.sys.liegend + prev.diabetes + prev.stroke.tia + 
                      prev.heart.failure + prev.niereninsuff + coronary.heart.disease, 
                    data = data)

summary(full_model)

# Define the range of models for stepwise selection (only consider changes in biomarkers)
scope <- list(lower = ~ age.bl + pat.sex + bmi + 
                current.smoker + rr.sys.liegend + prev.diabetes + prev.stroke.tia + 
                prev.heart.failure + prev.niereninsuff + coronary.heart.disease,
              upper = ~ ang2_log_std + ddi2h_log_std + 
                crcl_cg_log_std + cysc_log_std + alat_log_std + gdf.15_log_std + 
                crphs_log_std + igfbp7_log_std + il6_log_std + probnpii_log_std + 
                opn_log_std + tnt.hs_log_std + age.bl + pat.sex + bmi + 
                current.smoker + rr.sys.liegend + prev.diabetes + prev.stroke.tia + 
                prev.heart.failure + prev.niereninsuff + coronary.heart.disease)

# Perform stepwise selection for model refinement
reduced_model <- step(full_model, direction = "backward", scope = scope)

# Perform stepwise selection using stepAIC from MASS package
library(MASS) 
reduced_model <- stepAIC(full_model, scope = scope, direction = "backward")

# View the final selected model
summary(reduced_model)

# Use the biomarkers included in the step_model_filtered for the combined model
cox_model <- coxph(Surv(timeto1.stroke, stroke) ~ alat_log_std + 
                     il6_log_std + probnpii_log_std + opn_log_std + age.bl + pat.sex + 
                     bmi + current.smoker + rr.sys.liegend + prev.diabetes + prev.stroke.tia + 
                     prev.heart.failure + prev.niereninsuff + coronary.heart.disease, data = dat)
summary(cox_model)

cox_zph <- cox.zph(cox_model)
print(cox_zph)
plot(cox_zph)

########## BACKWARD BIOMARKER SELECTION USING CHADSVASC AS BASE MODEL ##########

# Perform a backward and forward selection of biomarkers
# Create a new dataset excluding patients with missing values in specified variables
# Exclude patients with missing variables
library(dplyr)

data <- dat %>% 
  select(ends_with("_log"), timeto1.stroke, stroke, chads)

data <- na.omit(data)

# Fit the Cox proportional hazards model with the filtered dataset
full_model <- coxph(Surv(timeto1.stroke, stroke) ~ ang2_log + ddi2h_log + 
                      crcl_cg_log + cysc_log + alat_log + gdf.15_log + 
                      crphs_log + igfbp7_log + il6_log + probnpii_log + 
                      opn_log + tnt.hs_log + chads, 
                    data = data)

summary(full_model)

# Define the range of models for stepwise selection (only consider changes in biomarkers)
scope <- list(lower = ~ chads,
              upper = ~ ang2_log + ddi2h_log + 
                crcl_cg_log + cysc_log + alat_log + gdf.15_log + 
                crphs_log + igfbp7_log + il6_log + probnpii_log + 
                opn_log + tnt.hs_log + chads)

# Perform stepwise selection for model refinement
reduced_model <- step(full_model, direction = "backward", scope = scope)

# Perform stepwise selection using stepAIC from MASS package
library(MASS) 
reduced_model <- stepAIC(full_model, scope = scope, direction = "backward")

# View the final selected model
summary(reduced_model)

# Use the biomarkers included in the step_model_filtered for the combined model
cox_model <- coxph(Surv(timeto1.stroke, stroke) ~ il6_log + probnpii_log + 
                     tnt.hs_log + chads, data = dat)
summary(cox_model)

cox_zph <- cox.zph(cox_model)
print(cox_zph)
plot(cox_zph)

################### CALCULATE AUC FOR BASE MODEL AND BASE MODEL + BIOMARKERS ###

# Load necessary libraries
library(survival)
library(timeROC)
library(pROC)

# Fit the Cox proportional hazards model on the complete cases
cox_model_base <- coxph(Surv(timeto1.stroke, stroke) ~ age.bl + pat.sex + bmi + 
                          current.smoker + rr.sys.liegend + 
                          prev.diabetes + prev.stroke.tia + prev.heart.failure + 
                          prev.niereninsuff + coronary.heart.disease, data = data)

cox_model_biomarkers <- coxph(Surv(timeto1.stroke, stroke) ~ alat_log_std + 
                                il6_log_std + probnpii_log_std + opn_log_std + age.bl + pat.sex + 
                                bmi + current.smoker + rr.sys.liegend + prev.diabetes + prev.stroke.tia + 
                                prev.heart.failure + prev.niereninsuff + coronary.heart.disease, data = data)

# Obtain the risk scores (linear predictors) for each model
risk_scores_base <- predict(cox_model_base, type = "lp", newdata = data)
risk_scores_biomarkers <- predict(cox_model_biomarkers, type = "lp", newdata = data)

# Define time points at which to calculate the AUC
time_points <- c(1, 2, 3) * 365.25  # 1, 2, and 3 years in days

# Calculate time-dependent AUC for the base model
roc_results_base <- timeROC(T = data$timeto1.stroke,
                            delta = data$stroke,
                            marker = risk_scores_base,
                            cause = 1,
                            weighting = "marginal",
                            times = time_points,
                            ROC = TRUE)

# Calculate time-dependent AUC for the biomarker model
roc_results_biomarkers <- timeROC(T = data$timeto1.stroke,
                                  delta = data$stroke,
                                  marker = risk_scores_biomarkers,
                                  cause = 1,
                                  weighting = "marginal",
                                  times = time_points,
                                  ROC = TRUE)

# Print AUC values for both models
print(roc_results_base$AUC)
print(roc_results_biomarkers$AUC)

# Extract predicted probabilities at the specific time point (1 year, 365.25 days)
pred_probs_base <- predict(cox_model_base, type = "risk", newdata = data)
pred_probs_biomarkers <- predict(cox_model_biomarkers, type = "risk", newdata = data)

# Use pROC to create ROC curves at the 1-year mark
roc_base <- roc(data$stroke, pred_probs_base, plot = FALSE, ci = TRUE, quiet = TRUE)
roc_biomarkers <- roc(data$stroke, pred_probs_biomarkers, plot = FALSE, ci = TRUE, quiet = TRUE)

# Perform DeLong's test to compare the AUCs
delong_test <- roc.test(roc_base, roc_biomarkers)
print(delong_test)

####################### AUC of CHADS-VASC + biomarkers #########################
library(survival)
library(timeROC)
library(pROC)

# Fit the Cox proportional hazards model on the complete cases
cox_model_base <- coxph(Surv(timeto1.stroke, stroke) ~ chads, data = dat)
cox_model_biomarkers <- coxph(Surv(timeto1.stroke, stroke) ~ il6_log + 
                                probnpii_log + tnt.hs_log + chads, data = dat)

# Obtain the risk scores (linear predictors) for each model
risk_scores_base <- predict(cox_model_base, type = "lp", newdata = dat)
risk_scores_biomarkers <- predict(cox_model_biomarkers, type = "lp", newdata = dat)

# Define time points at which to calculate the AUC
time_points <- c(1, 2, 3) * 365.25  # 1, 2, and 3 years in days

# Calculate time-dependent AUC for the base model
roc_results_base <- timeROC(T = dat$timeto1.stroke,
                            delta = dat$stroke,
                            marker = risk_scores_base,
                            cause = 1,
                            weighting = "marginal",
                            times = time_points,
                            ROC = TRUE)

# Calculate time-dependent AUC for the biomarker model
roc_results_biomarkers <- timeROC(T = dat$timeto1.stroke,
                                  delta = dat$stroke,
                                  marker = risk_scores_biomarkers,
                                  cause = 1,
                                  weighting = "marginal",
                                  times = time_points,
                                  ROC = TRUE)

# Print AUC values for both models
print(roc_results_base$AUC)
print(roc_results_biomarkers$AUC)

# Extract predicted probabilities at the specific time point (1 year, 365.25 days)
pred_probs_base <- predict(cox_model_base, type = "risk", newdata = dat)
pred_probs_biomarkers <- predict(cox_model_biomarkers, type = "risk", newdata = dat)

# Use pROC to create ROC curves at the 1-year mark
roc_base <- roc(dat$stroke, pred_probs_base, plot = FALSE, ci = TRUE, quiet = TRUE)
roc_biomarkers <- roc(dat$stroke, pred_probs_biomarkers, plot = FALSE, ci = TRUE, quiet = TRUE)

# Perform DeLong's test to compare the AUCs
delong_test <- roc.test(roc_base, roc_biomarkers)
print(delong_test)

# Extract ROC curve data for 1-year mark (365.25 days)
roc_data_base <- data.frame(
  FPR = 1 - roc_results_base$FP[, 1],
  TPR = roc_results_base$TP[, 1],
  model = "CHA2DS2-VASc Score"
)

roc_data_biomarkers <- data.frame(
  FPR = 1 - roc_results_biomarkers$FP[, 1],
  TPR = roc_results_biomarkers$TP[, 1],
  model = "CHA2DS2-VASc Score + IL6, NT-proBNP, TNTHS"
)

plot(roc_results_base, time = time_points[1], col = "blue", title = "CHA2DS2-VASc Score")
plot(roc_results_biomarkers, time = time_points[1], col = "red", title = "CHA2DS2-VASc Score + biomarkers", add = TRUE)

# Add a legend
legend("topright", 
       legend = c("CHA2DS2-VASc Score", "CHA2DS2-VASc Score + biomarkers"), 
       col = c("blue", "red"), 
       lwd = 2,
       inset = c(-0.4, -0.3),  # adjust the inset to move the legend further inside
       xpd = TRUE)







# Load necessary libraries
library(survival)
library(timeROC)
library(pROC)

# Fit the Cox proportional hazards model on the complete cases
cox_model_base <- coxph(Surv(timeto1.stroke, stroke) ~ chads, data = data)

cox_model_biomarkers <- coxph(Surv(timeto1.stroke, stroke) ~ chads + probnpii_log + il6_log, data = data)

# Obtain the risk scores (linear predictors) for each model
risk_scores_base <- predict(cox_model_base, type = "lp", newdata = data)
risk_scores_biomarkers <- predict(cox_model_biomarkers, type = "lp", newdata = data)

# Define time points at which to calculate the AUC
time_points <- 2 * 365.25

# Calculate time-dependent AUC for the base model
roc_results_base <- timeROC(T = data$timeto1.stroke,
                            delta = data$stroke,
                            marker = risk_scores_base,
                            cause = 1,
                            weighting = "marginal",
                            times = time_points,
                            ROC = TRUE)

# Calculate time-dependent AUC for the biomarker model
roc_results_biomarkers <- timeROC(T = data$timeto1.stroke,
                                  delta = data$stroke,
                                  marker = risk_scores_biomarkers,
                                  cause = 1,
                                  weighting = "marginal",
                                  times = time_points,
                                  ROC = TRUE)

# Print AUC values for both models
print(roc_results_base$AUC)
print(roc_results_biomarkers$AUC)

# Extract predicted probabilities at the specific time point (1 year, 365.25 days)
pred_probs_base <- predict(cox_model_base, type = "risk", newdata = data)
pred_probs_biomarkers <- predict(cox_model_biomarkers, type = "risk", newdata = data)

# Use pROC to create ROC curves at the 1-year mark
roc_base <- roc(data$stroke, pred_probs_base, plot = FALSE, ci = TRUE, quiet = TRUE)
roc_biomarkers <- roc(data$stroke, pred_probs_biomarkers, plot = FALSE, ci = TRUE, quiet = TRUE)

# Perform DeLong's test to compare the AUCs
delong_test <- roc.test(roc_base, roc_biomarkers)
print(delong_test)





















############################### VOLCANO PLOTS ##################################

# Create an empty dataframe to store results
volcano_data <- data.frame(Biomarker = character(),
                           HR = numeric(),
                           P_Value = numeric(),
                           stringsAsFactors = FALSE)

# Iterate through each biomarker
for (biomarker in biomarkers) {
  
  # Fit the Cox model
  cox_model <- coxph(Surv(timeto1.stroke, stroke) ~
                       get(paste0(biomarker, "_log_std")) + age.bl + pat.sex + bmi +
                       current.smoker + rr.sys.liegend + prev.diabetes +
                       prev.stroke.tia + prev.heart.failure + prev.niereninsuff +
                       coronary.heart.disease, data = dat)
  
  # Extract HR and p-value for _log_std
  coef_summary <- coef(summary(cox_model))
  hr <- coef_summary["get(paste0(biomarker, \"_log_std\"))", "exp(coef)"]
  p_value <- coef_summary["get(paste0(biomarker, \"_log_std\"))", "Pr(>|z|)"]
  
  # Add results to the dataframe
  volcano_data <- rbind(volcano_data, data.frame(biomarker = biomarker, HR = hr, P_Value = p_value))
}

# Volcano plot from the combined model after stepwise biomarker selection
volcano_data <- data.frame(Biomarker = character(),
                           HR = numeric(),
                           P_Value = numeric(),
                           stringsAsFactors = FALSE)

biomarker_names <- c("alat", "il6", "probnpii", "opn")

cox_model <- coxph(Surv(timeto1.stroke, stroke) ~ alat_log_std + 
                     il6_log_std + probnpii_log_std + opn_log_std + age.bl + pat.sex + 
                     bmi + current.smoker + rr.sys.liegend + prev.diabetes + prev.stroke.tia + 
                     prev.heart.failure + prev.niereninsuff + coronary.heart.disease, data = dat)

# Obtain summary of coefficients
coef_summary <- summary(cox_model)

hr <- coef_summary$coefficients[, "exp(coef)"][paste0(biomarker_names, "_log_std")]
p_value <- coef_summary$coefficients[, "Pr(>|z|)"][paste0(biomarker_names, "_log_std")]

# Combine results into a data frame
volcano_data <- rbind(volcano_data, data.frame(biomarker = biomarker_names, HR = hr, P_Value = p_value))

# Display the dataframe
print(volcano_data)

new_names <- c("ALAT", "IL-6", "NT-proBNP", "OPN")

# Replace biomarker names in coef_data
volcano_data$biomarker <- ifelse(volcano_data$biomarker == "alat", new_names[1],
                                 ifelse(volcano_data$biomarker == "il6", new_names[2],
                                        ifelse(volcano_data$biomarker == "probnpii", new_names[3], 
                                               ifelse(volcano_data$biomarker == "opn", new_names[4],
                                                      volcano_data$biomarker))))

library(ggplot2)
library(ggrepel)

ggplot(volcano_data, aes(x = HR, y = -log10(P_Value), label = biomarker)) +
  geom_point() +
  labs(x = "Standardized hazard ratio", y = "-Log10 P-value") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "red") +
  theme_bw() + 
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 9, color = "black"),
        axis.text.x = element_text(hjust = 0.5)) +
  coord_cartesian(xlim = c(0.5, 2.0)) +
  scale_x_continuous(trans = "log10", breaks = c(0.5, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2.0), labels = c(0.5, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2.0)) +
  geom_text_repel(size = 2.75)

########################## VARIABLE IMPORTANCE PLOTS ###########################

library(dplyr)

features <- dat %>% 
  select(ends_with("_log_std"), age.bl, pat.sex, bmi, current.smoker, rr.sys.liegend,
         prev.diabetes, prev.stroke.tia, prev.heart.failure, prev.niereninsuff,
         coronary.heart.disease, timeto1.stroke, stroke)

features <- na.omit(features)

cox_model <- coxph(Surv(timeto1.stroke, stroke) ~ alat_log_std + 
                     il6_log_std + probnpii_log_std + opn_log_std + age.bl + pat.sex + 
                     bmi + current.smoker + rr.sys.liegend + prev.diabetes + prev.stroke.tia + 
                     prev.heart.failure + prev.niereninsuff + coronary.heart.disease, data = features)
summary(cox_model)

# Extract coefficient estimates and standard errors
coef_est <- coef(cox_model)
coef_se <- sqrt(diag(vcov(cox_model)))

# Calculate partial chi-squared statistic for each predictor
partial_chisq <- (coef_est / coef_se)^2

# Extract the total degrees of freedom for the model
num_params <- length(coef(cox_model))
df_total <- num_params

# Adjust partial chi-squared statistic by subtracting total degrees of freedom
partial_chisq_adj <- partial_chisq - df_total

# Display the results
result <- data.frame(Variable = names(coef_est), Partial_ChiSq = partial_chisq)
print(result)

# Order the result by partial chi-squared statistic in descending order
result <- result[order(-result$Partial_ChiSq), ]

library(ggplot2)

# Define the desired label changes
label_changes <- c("opn_log_std" = "OPN",
                   "il6_log_std" = "IL-6",
                   "probnpii_log_std" = "NT-proBNP",
                   "alat_log_std" = "ALAT",
                   "age.bl" = "Age",
                   "pat.sex" = "Sex",
                   "bmi" = "BMI",
                   "current.smoker" = "Smoker",
                   "rr.sys.liegend" = "SBP",
                   "prev.diabetes" = "Diabetes",
                   "prev.stroke.tia" = "Prior stroke/TIA",
                   "prev.heart.failure" = "Heart failure",
                   "prev.niereninsuff" = "CKD",
                   "coronary.heart.disease" = "CAD")

# Create the dot plot
p <- ggplot(result, aes(x = Partial_ChiSq, y = reorder(Variable, Partial_ChiSq))) +
  geom_point(shape = 1, size = 2.5, color = "black") +  # Change shape to open circle and set color to black
  labs(x = "Partial χ2-df", size = 2, y = "", color = "black") +  # Fixed x-axis label and removed y-axis title
  scale_y_discrete(labels = label_changes) +  # Apply label changes to y-axis
  theme_bw() +  # Minimal theme
  theme(axis.text.x = element_text(size = 10),  # Adjust x-axis text size
        axis.text.y = element_text(size = 10),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(color = "black"),  # Set x-axis ticks color to black
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.major.y = element_line(linetype = "dotted", color = "grey")  # Specify major horizontal grid lines as dotted and grey
  )

print(p)

################################ FOREST PLOTS ##################################

# Combined model
library(randomForest)
library(survminer)

cox_model <- coxph(Surv(timeto1.stroke, stroke) ~ alat_log_std + 
                     il6_log_std + probnpii_log_std + opn_log_std + age.bl + pat.sex + 
                     bmi + current.smoker + rr.sys.liegend + prev.diabetes + prev.stroke.tia + 
                     prev.heart.failure + prev.niereninsuff + coronary.heart.disease, data = dat)
summary(cox_model)

biomarkers <- c("alat_log_std", "il6_log_std", "probnpii_log_std", "opn_log_std")

# Create a dataframe to store hazard ratios and 95% CIs
coef_data <- data.frame(
  Biomarker = biomarkers,
  HR = exp(coef(cox_model)[biomarkers]),             # Extracting hazard ratios
  Lower_CI = exp(confint(cox_model)[biomarkers, 1]),  # Extracting lower CIs
  Upper_CI = exp(confint(cox_model)[biomarkers, 2]),   # Extracting upper CIs
  p_value = summary(cox_model)$coefficients[biomarkers, "Pr(>|z|)"]  # Extracting p-values
)

# View the dataframe
print(coef_data)

new_names <- c("ALAT", "IL-6", "NT-proBNP", "OPN")

# Replace biomarker names in coef_data
coef_data$Biomarker <- ifelse(coef_data$Biomarker == "alat_log_std", new_names[1],
                              ifelse(coef_data$Biomarker == "il6_log_std", new_names[2],
                                     ifelse(coef_data$Biomarker == "probnpii_log_std", new_names[3],
                                            ifelse(coef_data$Biomarker == "opn_log_std", new_names[4], 
                                                   coef_data$Biomarker))))

coef_data$p.signif <- ifelse(coef_data$p_value < 0.05, "*", "")

library(ggplot2)

# Create coefficient plot using ggplot2
coefficient_plot <- ggplot(coef_data, aes(x = Biomarker, y = HR, color = Biomarker)) +
  geom_point(size = 2.5) + 
  geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 0.3, linetype = "solid") + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") + 
  labs(x = "Biomarker", y = "Hazard ratio (95% CI) per 1 SD increase in biomarker levels", size=11, color = "black") + 
  theme_bw() +
  theme(legend.position="none",
        axis.title = element_text(face="bold"),
        axis.text.x = element_text(size=9, color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(size=9, color = "black")) +
  scale_y_continuous(limits = c(0.7, 1.9), breaks = seq(0.7, 1.9, by = 0.1)) +
  geom_text(aes(label = paste0(sprintf("%.2f", HR))), 
            hjust = -0.3, size = 3)
print(coefficient_plot)

################################################################################
# Myocardial infarction
################################################################################

#-------------------------------------------------------------------------------
# Check whether there are missing "time" variables 
num_missing_timeto1.mi <- sum(is.na(dat$timeto1.mi))
print(num_missing_timeto1.mi)
#-------------------------------------------------------------------------------

################### COX REGRESSION MODELS ######################################

# Fit the age + sex adjusted Cox model
library(survival)

# Create empty lists to store results from age+sex and multivariable adjusted models
results.age_summary <- list()
results.multi_summary <- list()
results.age <- list()
results.multi <- list()

# Age+sex adjusted: Iterate over each biomarker
for (biomarker in biomarkers) {
  
  # Fit the Cox model
  cox_model <- coxph(Surv(timeto1.mi, mi) ~
                       get(paste0(biomarker, "_log_std")) + age.bl + pat.sex, data = dat)
  
  # Store results in the list
  results.age_summary[[biomarker]] <- summary(cox_model)
  
  # Store results in the list for volcano plots
  results.age[[biomarker]] <- coef(summary(cox_model))[, "exp(coef)"]
}

# Multivariable: Iterate over each biomarker
for (biomarker in biomarkers) {
  
  # Fit the Cox model
  cox_model <- coxph(Surv(timeto1.mi, mi) ~
                       get(paste0(biomarker, "_log_std")) + age.bl + pat.sex + bmi +
                       current.smoker + rr.sys.liegend + prev.diabetes +
                       prev.stroke.tia + prev.heart.failure + prev.niereninsuff +
                       coronary.heart.disease, data = dat)
  
  # Store results in the list
  results.multi_summary[[biomarker]] <- summary(cox_model)
  
  # Store results in the list for volcano plots
  results.multi[[biomarker]] <- coef(summary(cox_model))[, "exp(coef)"]
}

# Combined model 1 (only significant biomarkers)
cox_model <- coxph(Surv(timeto1.mi, mi) ~ gdf.15_log_std + il6_log_std + 
                     probnpii_log_std + tnt.hs_log_std + age.bl + pat.sex + 
                     bmi + current.smoker + rr.sys.liegend + 
                     prev.diabetes + prev.stroke.tia + prev.heart.failure + 
                     prev.niereninsuff + coronary.heart.disease, data = dat)
summary(cox_model)

# Combined model 2 (all biomarkers)
cox_model <- coxph(Surv(timeto1.mi, mi) ~ ang2_log_std + ddi2h_log_std + 
                     crcl_cg_log_std + cysc_log_std + alat_log_std + gdf.15_log_std + crphs_log_std + 
                     igfbp7_log_std + il6_log_std + probnpii_log_std + opn_log_std + tnt.hs_log_std + 
                     age.bl + pat.sex + bmi + current.smoker + rr.sys.liegend + 
                     prev.diabetes + prev.stroke.tia + prev.heart.failure + 
                     prev.niereninsuff + coronary.heart.disease, data = dat)
summary(cox_model)

################### BACKWARD BIOMARKER SELECTION ###############################

# Perform a backward and forward selection of biomarkers
# Create a new dataset excluding patients with missing values in specified variables
# Exclude patients with missing variables
library(dplyr)

data <- dat %>% 
  select(ends_with("_log_std"), age.bl, pat.sex, bmi, current.smoker, rr.sys.liegend,
         prev.diabetes, prev.stroke.tia, prev.heart.failure, prev.niereninsuff,
         coronary.heart.disease, timeto1.mi, mi)

data <- na.omit(data)

# Fit the Cox proportional hazards model with the filtered dataset
full_model <- coxph(Surv(timeto1.mi, mi) ~ ang2_log_std + ddi2h_log_std + 
                      crcl_cg_log_std + cysc_log_std + alat_log_std + gdf.15_log_std + 
                      crphs_log_std + igfbp7_log_std + il6_log_std + probnpii_log_std + 
                      opn_log_std + tnt.hs_log_std + age.bl + pat.sex + bmi + 
                      current.smoker + rr.sys.liegend + prev.diabetes + prev.stroke.tia + 
                      prev.heart.failure + prev.niereninsuff + coronary.heart.disease, 
                    data = data)

summary(full_model)

# Define the range of models for stepwise selection (only consider changes in biomarkers)
scope <- list(lower = ~ age.bl + pat.sex + bmi + 
                current.smoker + rr.sys.liegend + prev.diabetes + prev.stroke.tia + 
                prev.heart.failure + prev.niereninsuff + coronary.heart.disease,
              upper = ~ ang2_log_std + ddi2h_log_std + 
                crcl_cg_log_std + cysc_log_std + alat_log_std + gdf.15_log_std + 
                crphs_log_std + igfbp7_log_std + il6_log_std + probnpii_log_std + 
                opn_log_std + tnt.hs_log_std + age.bl + pat.sex + bmi + 
                current.smoker + rr.sys.liegend + prev.diabetes + prev.stroke.tia + 
                prev.heart.failure + prev.niereninsuff + coronary.heart.disease)

# Perform stepwise selection for model refinement
reduced_model <- step(full_model, direction = "backward", scope = scope)

# Perform stepwise selection using stepAIC from MASS package
library(MASS) 
reduced_model <- stepAIC(full_model, scope = scope, direction = "backward")

# View the final selected model
summary(reduced_model)

# Use the biomarkers included in the step_model_filtered for the combined model
cox_model <- coxph(Surv(timeto1.mi, mi) ~ il6_log_std + tnt.hs_log_std + 
                     age.bl + pat.sex + bmi + current.smoker + rr.sys.liegend + 
                     prev.diabetes + prev.stroke.tia + prev.heart.failure + 
                     prev.niereninsuff + coronary.heart.disease, data = dat)
summary(cox_model)

################ AUC FOR BASE MODEL AND BASE MODEL + BIOMARKERS ################

# Load necessary libraries
library(survival)
library(timeROC)
library(pROC)

# Fit the Cox proportional hazards model on the complete cases
cox_model_base <- coxph(Surv(timeto1.mi, mi) ~ age.bl + pat.sex + bmi + 
                          current.smoker + rr.sys.liegend + 
                          prev.diabetes + prev.stroke.tia + prev.heart.failure + 
                          prev.niereninsuff + coronary.heart.disease, data = data)

cox_model_biomarkers <- coxph(Surv(timeto1.mi, mi) ~ il6_log_std + tnt.hs_log_std + 
                                age.bl + pat.sex + bmi + current.smoker + rr.sys.liegend + 
                                prev.diabetes + prev.stroke.tia + prev.heart.failure + 
                                prev.niereninsuff + coronary.heart.disease, data = data)

# Obtain the risk scores (linear predictors) for each model
risk_scores_base <- predict(cox_model_base, type = "lp", newdata = data)
risk_scores_biomarkers <- predict(cox_model_biomarkers, type = "lp", newdata = data)

# Define time points at which to calculate the AUC
time_points <- c(1, 2, 3) * 365.25  # 1, 2, and 3 years in days

# Calculate time-dependent AUC for the base model
roc_results_base <- timeROC(T = data$timeto1.mi,
                            delta = data$mi,
                            marker = risk_scores_base,
                            cause = 1,
                            weighting = "marginal",
                            times = time_points,
                            ROC = TRUE)

# Calculate time-dependent AUC for the biomarker model
roc_results_biomarkers <- timeROC(T = data$timeto1.mi,
                                  delta = data$mi,
                                  marker = risk_scores_biomarkers,
                                  cause = 1,
                                  weighting = "marginal",
                                  times = time_points,
                                  ROC = TRUE)

# Print AUC values for both models
print(roc_results_base$AUC)
print(roc_results_biomarkers$AUC)

# Extract predicted probabilities at the specific time point (1 year, 365.25 days)
pred_probs_base <- predict(cox_model_base, type = "risk", newdata = data)
pred_probs_biomarkers <- predict(cox_model_biomarkers, type = "risk", newdata = data)

# Use pROC to create ROC curves at the 1-year mark
roc_base <- roc(data$mi, pred_probs_base, plot = FALSE, ci = TRUE, quiet = TRUE)
roc_biomarkers <- roc(data$mi, pred_probs_biomarkers, plot = FALSE, ci = TRUE, quiet = TRUE)

# Perform DeLong's test to compare the AUCs
delong_test <- roc.test(roc_base, roc_biomarkers)
print(delong_test)

# Plot ROC Curves at different time points
plot(roc_results, time = time_points[1], col = "blue", title = "ROC Curve at 1 Year")
plot(roc_results, time = time_points[2], col = "red", title = "ROC Curve at 2 Years", add = TRUE)
plot(roc_results, time = time_points[3], col = "green", title = "ROC Curve at 3 Years", add = TRUE)
legend("bottomright", legend = c("1 Year", "2 Years", "3 Years"), col = c("blue", "red", "green"), lwd = 2)

############################### VOLCANO PLOTS ##################################

# Create an empty dataframe to store results
volcano_data <- data.frame(Biomarker = character(),
                           HR = numeric(),
                           P_Value = numeric(),
                           stringsAsFactors = FALSE)

# Iterate through each biomarker
for (biomarker in biomarkers) {
  
  # Fit the Cox model
  cox_model <- coxph(Surv(timeto1.mi, mi) ~
                       get(paste0(biomarker, "_log_std")) + age.bl + pat.sex + bmi +
                       current.smoker + rr.sys.liegend + prev.diabetes +
                       prev.stroke.tia + prev.heart.failure + prev.niereninsuff +
                       coronary.heart.disease, data = dat)
  
  # Extract HR and p-value for _log_std
  coef_summary <- coef(summary(cox_model))
  hr <- coef_summary["get(paste0(biomarker, \"_log_std\"))", "exp(coef)"]
  p_value <- coef_summary["get(paste0(biomarker, \"_log_std\"))", "Pr(>|z|)"]
  
  # Add results to the dataframe
  volcano_data <- rbind(volcano_data, data.frame(biomarker = biomarker, HR = hr, P_Value = p_value))
}

# Volcano plot from the combined model after stepwise biomarker selection
volcano_data <- data.frame(Biomarker = character(),
                           HR = numeric(),
                           P_Value = numeric(),
                           stringsAsFactors = FALSE)

biomarker_names <- c("il6", "tnt.hs")

cox_model <- coxph(Surv(timeto1.mi, mi) ~ il6_log_std + tnt.hs_log_std + 
                     age.bl + pat.sex + bmi + current.smoker + rr.sys.liegend + 
                     prev.diabetes + prev.stroke.tia + prev.heart.failure + 
                     prev.niereninsuff + coronary.heart.disease, data = dat)

# Obtain summary of coefficients
coef_summary <- summary(cox_model)

hr <- coef_summary$coefficients[, "exp(coef)"][paste0(biomarker_names, "_log_std")]
p_value <- coef_summary$coefficients[, "Pr(>|z|)"][paste0(biomarker_names, "_log_std")]

# Combine results into a data frame
volcano_data <- rbind(volcano_data, data.frame(biomarker = biomarker_names, HR = hr, P_Value = p_value))

# Display the dataframe
print(volcano_data)

new_names <- c("IL-6", "hsTropT")

# Replace biomarker names in coef_data
volcano_data$biomarker <- ifelse(volcano_data$biomarker == "il6", new_names[1],
                                 ifelse(volcano_data$biomarker == "tnt.hs", new_names[2],
                                        volcano_data$biomarker))

library(ggplot2)
library(ggrepel)

ggplot(volcano_data, aes(x = HR, y = -log10(P_Value), label = biomarker)) +
  geom_point() +
  labs(x = "Standardized hazard ratio", y = "-Log10 P-value") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "red") +
  theme_bw() + 
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 9, color = "black"),
        axis.text.x = element_text(hjust = 0.5)) +
  coord_cartesian(xlim = c(0.5, 2.0)) +
  scale_x_continuous(trans = "log10", breaks = c(0.5, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2.0), labels = c(0.5, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2.0)) +
  geom_text_repel(size = 2.75)

########################## VARIABLE IMPORTANCE PLOTS ###########################

library(dplyr)

features <- dat %>% 
  select(ends_with("_log_std"), age.bl, pat.sex, bmi, current.smoker, rr.sys.liegend,
         prev.diabetes, prev.stroke.tia, prev.heart.failure, prev.niereninsuff,
         coronary.heart.disease, timeto1.mi, mi)

features <- na.omit(features)

cox_model <- coxph(Surv(timeto1.mi, mi) ~ il6_log_std + tnt.hs_log_std + 
                     age.bl + pat.sex + bmi + current.smoker + rr.sys.liegend + 
                     prev.diabetes + prev.stroke.tia + prev.heart.failure + 
                     prev.niereninsuff + coronary.heart.disease, data = data)
summary(cox_model)

# Extract coefficient estimates and standard errors
coef_est <- coef(cox_model)
coef_se <- sqrt(diag(vcov(cox_model)))

# Calculate partial chi-squared statistic for each predictor
partial_chisq <- (coef_est / coef_se)^2

# Extract the total degrees of freedom for the model
num_params <- length(coef(cox_model))
df_total <- num_params

# Adjust partial chi-squared statistic by subtracting total degrees of freedom
partial_chisq_adj <- partial_chisq - df_total

# Display the results
result <- data.frame(Variable = names(coef_est), Partial_ChiSq = partial_chisq)
print(result)

# Order the result by partial chi-squared statistic in descending order
result <- result[order(-result$Partial_ChiSq), ]

library(ggplot2)

# Define the desired label changes
label_changes <- c("il6_log_std" = "IL-6",
                   "tnt.hs_log_std" = "hsTropT",
                   "age.bl" = "Age",
                   "pat.sex" = "Sex",
                   "bmi" = "BMI",
                   "current.smoker" = "Smoker",
                   "rr.sys.liegend" = "SBP",
                   "prev.diabetes" = "Diabetes",
                   "prev.stroke.tia" = "Prior stroke/TIA",
                   "prev.heart.failure" = "Heart failure",
                   "prev.niereninsuff" = "CKD",
                   "coronary.heart.disease" = "CAD")

# Create the dot plot
p <- ggplot(result, aes(x = Partial_ChiSq, y = reorder(Variable, Partial_ChiSq))) +
  geom_point(shape = 1, size = 2.5, color = "black") +  # Change shape to open circle and set color to black
  labs(x = "Partial χ2-df", size = 2, y = "", color = "black") +  # Fixed x-axis label and removed y-axis title
  scale_y_discrete(labels = label_changes) +  # Apply label changes to y-axis
  theme_bw() +  # Minimal theme
  theme(axis.text.x = element_text(size = 10),  # Adjust x-axis text size
        axis.text.y = element_text(size = 10),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(color = "black"),  # Set x-axis ticks color to black
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.major.y = element_line(linetype = "dotted", color = "grey")  # Specify major horizontal grid lines as dotted and grey
  )

print(p)

################################ FOREST PLOTS ##################################

# Combined model
library(randomForest)
library(survminer)

cox_model <- coxph(Surv(timeto1.mi, mi) ~ il6_log_std + tnt.hs_log_std + 
                     age.bl + pat.sex + bmi + current.smoker + rr.sys.liegend + 
                     prev.diabetes + prev.stroke.tia + prev.heart.failure + 
                     prev.niereninsuff + coronary.heart.disease, data = dat)
summary(cox_model)

biomarkers <- c("il6_log_std", "tnt.hs_log_std")

# Create a dataframe to store hazard ratios and 95% CIs
coef_data <- data.frame(
  Biomarker = biomarkers,
  HR = exp(coef(cox_model)[biomarkers]),             # Extracting hazard ratios
  Lower_CI = exp(confint(cox_model)[biomarkers, 1]),  # Extracting lower CIs
  Upper_CI = exp(confint(cox_model)[biomarkers, 2]),   # Extracting upper CIs
  p_value = summary(cox_model)$coefficients[biomarkers, "Pr(>|z|)"]  # Extracting p-values
)

# View the dataframe
print(coef_data)

new_names <- c("IL-6", "hsTropT")

# Replace biomarker names in coef_data
coef_data$Biomarker <- ifelse(coef_data$Biomarker == "il6_log_std", new_names[1],
                              ifelse(coef_data$Biomarker == "tnt.hs_log_std", new_names[2],
                                     coef_data$Biomarker))

coef_data$p.signif <- ifelse(coef_data$p_value < 0.05, "*", "")

library(ggplot2)

# Create coefficient plot using ggplot2
coefficient_plot <- ggplot(coef_data, aes(x = Biomarker, y = HR, color = Biomarker)) +
  geom_point(size = 2.5) + 
  geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 0.3, linetype = "solid") + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") + 
  labs(x = "Biomarker", y = "Hazard ratio (95% CI) per 1 SD increase in biomarker levels", size=11, color = "black") + 
  theme_bw() +
  theme(legend.position="none",
        axis.title = element_text(face="bold"),
        axis.text.x = element_text(size=9, color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(size=9, color = "black")) +
  scale_y_continuous(limits = c(0.7, 1.7), breaks = seq(0.7, 1.7, by = 0.1)) +
  geom_text(aes(label = paste0(sprintf("%.2f", HR))), 
            hjust = -0.3, size = 3)
print(coefficient_plot)


################################################################################
# Cardiovascular death
################################################################################

#-------------------------------------------------------------------------------
# Check whether there are missing "time" variables 
num_missing_timeto1.death <- sum(is.na(dat$timeto1.death))
print(num_missing_timeto1.death)
#-------------------------------------------------------------------------------

############################ COX REGRESSION MODELS #############################

# Fit the age + sex adjusted Cox model
library(survival)

# Create empty lists to store results from age+sex and multivariable adjusted models
results.age_summary <- list()
results.multi_summary <- list()
results.age <- list()
results.multi <- list()

# Age+sex adjusted: Iterate over each biomarker
for (biomarker in biomarkers) {
  
  # Fit the Cox model
  cox_model <- coxph(Surv(timeto1.death, death.cardiac) ~
                       get(paste0(biomarker, "_log_std")) + age.bl + pat.sex, data = dat)
  
  # Store results in the list
  results.age_summary[[biomarker]] <- summary(cox_model)
  
  # Store results in the list for volcano plots
  results.age[[biomarker]] <- coef(summary(cox_model))[, "exp(coef)"]
}

# Multivariable: Iterate over each biomarker
for (biomarker in biomarkers) {
  
  # Fit the Cox model
  cox_model <- coxph(Surv(timeto1.death, death.cardiac) ~
                       get(paste0(biomarker, "_log_std")) + age.bl + pat.sex + bmi +
                       current.smoker + rr.sys.liegend + prev.diabetes +
                       prev.stroke.tia + prev.heart.failure + prev.niereninsuff +
                       coronary.heart.disease, data = dat)
  
  # Store results in the list
  results.multi_summary[[biomarker]] <- summary(cox_model)
  
  # Store results in the list for volcano plots
  results.multi[[biomarker]] <- coef(summary(cox_model))[, "exp(coef)"]
}

# Combined model 2 (all biomarkers)
cox_model <- coxph(Surv(timeto1.death, death.cardiac) ~ ang2_log_std + ddi2h_log_std + 
                     crep2_log_std + cysc_log_std + alat_log_std + gdf.15_log_std + crphs_log_std + 
                     igfbp7_log_std + il6_log_std + probnpii_log_std + opn_log_std + tnt.hs_log_std + 
                     age.bl + pat.sex + bmi + current.smoker + rr.sys.liegend + 
                     prev.diabetes + prev.stroke.tia + prev.heart.failure + 
                     prev.niereninsuff + coronary.heart.disease, data = dat)
summary(cox_model)

################### BACKWARD BIOMARKER SELECTION ###############################

# Perform a backward and forward selection of biomarkers
# Create a new dataset excluding patients with missing values in specified variables
# Exclude patients with missing variables
library(dplyr)

data <- dat %>% 
  select(ends_with("_log_std"), age.bl, pat.sex, bmi, current.smoker, rr.sys.liegend,
         prev.diabetes, prev.stroke.tia, prev.heart.failure, prev.niereninsuff,
         coronary.heart.disease, timeto1.death, death.cardiac)

data <- na.omit(data)

# Fit the Cox proportional hazards model with the filtered dataset
full_model <- coxph(Surv(timeto1.death, death.cardiac) ~ ang2_log_std + ddi2h_log_std + 
                      crcl_cg_log_std + cysc_log_std + alat_log_std + gdf.15_log_std + 
                      crphs_log_std + igfbp7_log_std + il6_log_std + probnpii_log_std + 
                      opn_log_std + tnt.hs_log_std + age.bl + pat.sex + bmi + 
                      current.smoker + rr.sys.liegend + prev.diabetes + prev.stroke.tia + 
                      prev.heart.failure + prev.niereninsuff + coronary.heart.disease, 
                    data = data)

summary(full_model)

# Define the range of models for stepwise selection (only consider changes in biomarkers)
scope <- list(lower = ~ age.bl + pat.sex + bmi + 
                current.smoker + rr.sys.liegend + prev.diabetes + prev.stroke.tia + 
                prev.heart.failure + prev.niereninsuff + coronary.heart.disease,
              upper = ~ ang2_log_std + ddi2h_log_std + 
                crcl_cg_log_std + cysc_log_std + alat_log_std + gdf.15_log_std + 
                crphs_log_std + igfbp7_log_std + il6_log_std + probnpii_log_std + 
                opn_log_std + tnt.hs_log_std + age.bl + pat.sex + bmi + 
                current.smoker + rr.sys.liegend + prev.diabetes + prev.stroke.tia + 
                prev.heart.failure + prev.niereninsuff + coronary.heart.disease)

# Perform stepwise selection for model refinement
reduced_model <- step(full_model, direction = "backward", scope = scope)

# Perform stepwise selection using stepAIC from MASS package
library(MASS) 
reduced_model <- stepAIC(full_model, scope = scope, direction = "backward")

# View the final selected model
summary(reduced_model)

# Use the biomarkers included in the step_model_filtered for the combined model
# Exclude CRP
cox_model <- coxph(Surv(timeto1.death, death.cardiac) ~ ddi2h_log_std + 
                     alat_log_std + gdf.15_log_std + 
                     il6_log_std + probnpii_log_std + tnt.hs_log_std + 
                     age.bl + pat.sex + bmi + current.smoker + rr.sys.liegend + 
                     prev.diabetes + prev.stroke.tia + prev.heart.failure + 
                     prev.niereninsuff + coronary.heart.disease, data = dat)
summary(cox_model)

################ AUC FOR BASE MODEL AND BASE MODEL + BIOMARKERS ################

# Load necessary libraries
library(survival)
library(timeROC)
library(pROC)

# Fit the Cox proportional hazards model on the complete cases
cox_model_base <- coxph(Surv(timeto1.death, death.cardiac) ~ age.bl + pat.sex + bmi + 
                          current.smoker + rr.sys.liegend + 
                          prev.diabetes + prev.stroke.tia + prev.heart.failure + 
                          prev.niereninsuff + coronary.heart.disease, data = data)

cox_model_biomarkers <- coxph(Surv(timeto1.death, death.cardiac) ~ ddi2h_log_std + 
                                alat_log_std + gdf.15_log_std + 
                                il6_log_std + probnpii_log_std + tnt.hs_log_std + 
                                age.bl + pat.sex + bmi + current.smoker + rr.sys.liegend + 
                                prev.diabetes + prev.stroke.tia + prev.heart.failure + 
                                prev.niereninsuff + coronary.heart.disease, data = data)

# Obtain the risk scores (linear predictors) for each model
risk_scores_base <- predict(cox_model_base, type = "lp", newdata = data)
risk_scores_biomarkers <- predict(cox_model_biomarkers, type = "lp", newdata = data)

# Define time points at which to calculate the AUC
time_points <- c(1, 2, 3) * 365.25  # 1, 2, and 3 years in days

# Calculate time-dependent AUC for the base model
roc_results_base <- timeROC(T = data$timeto1.death,
                            delta = data$death.cardiac,
                            marker = risk_scores_base,
                            cause = 1,
                            weighting = "marginal",
                            times = time_points,
                            ROC = TRUE)

# Calculate time-dependent AUC for the biomarker model
roc_results_biomarkers <- timeROC(T = data$timeto1.death,
                                  delta = data$death.cardiac,
                                  marker = risk_scores_biomarkers,
                                  cause = 1,
                                  weighting = "marginal",
                                  times = time_points,
                                  ROC = TRUE)

# Print AUC values for both models
print(roc_results_base$AUC)
print(roc_results_biomarkers$AUC)

# Extract predicted probabilities at the specific time point (1 year, 365.25 days)
pred_probs_base <- predict(cox_model_base, type = "risk", newdata = data)
pred_probs_biomarkers <- predict(cox_model_biomarkers, type = "risk", newdata = data)

# Use pROC to create ROC curves at the 1-year mark
roc_base <- roc(data$death.cardiac, pred_probs_base, plot = FALSE, ci = TRUE, quiet = TRUE)
roc_biomarkers <- roc(data$death.cardiac, pred_probs_biomarkers, plot = FALSE, ci = TRUE, quiet = TRUE)

# Perform DeLong's test to compare the AUCs
delong_test <- roc.test(roc_base, roc_biomarkers)
print(delong_test)

# Plot ROC Curves at different time points
plot(roc_results, time = time_points[1], col = "blue", title = "ROC Curve at 1 Year")
plot(roc_results, time = time_points[2], col = "red", title = "ROC Curve at 2 Years", add = TRUE)
plot(roc_results, time = time_points[3], col = "green", title = "ROC Curve at 3 Years", add = TRUE)
legend("bottomright", legend = c("1 Year", "2 Years", "3 Years"), col = c("blue", "red", "green"), lwd = 2)

############################### VOLCANO PLOTS ##################################

# Create an empty dataframe to store results
volcano_data <- data.frame(Biomarker = character(),
                           HR = numeric(),
                           P_Value = numeric(),
                           stringsAsFactors = FALSE)

# Iterate through each biomarker
for (biomarker in biomarkers) {
  
  # Fit the Cox model
  cox_model <- coxph(Surv(timeto1.death, death.cardiac) ~
                       get(paste0(biomarker, "_log_std")) + age.bl + pat.sex + bmi +
                       current.smoker + rr.sys.liegend + prev.diabetes +
                       prev.stroke.tia + prev.heart.failure + prev.niereninsuff +
                       coronary.heart.disease, data = dat)
  
  # Extract HR and p-value for _log_std
  coef_summary <- coef(summary(cox_model))
  hr <- coef_summary["get(paste0(biomarker, \"_log_std\"))", "exp(coef)"]
  p_value <- coef_summary["get(paste0(biomarker, \"_log_std\"))", "Pr(>|z|)"]
  
  # Add results to the dataframe
  volcano_data <- rbind(volcano_data, data.frame(biomarker = biomarker, HR = hr, P_Value = p_value))
}

# Volcano plot from the combined model after stepwise biomarker selection
volcano_data <- data.frame(Biomarker = character(),
                           HR = numeric(),
                           P_Value = numeric(),
                           stringsAsFactors = FALSE)

biomarker_names <- c("ddi2h", "alat", "gdf.15",
                     "il6", "probnpii", "tnt.hs")

cox_model <- coxph(Surv(timeto1.death, death.cardiac) ~ ddi2h_log_std + 
                     alat_log_std + gdf.15_log_std + il6_log_std + probnpii_log_std + 
                     tnt.hs_log_std + age.bl + pat.sex + bmi + current.smoker + 
                     rr.sys.liegend + prev.diabetes + prev.stroke.tia + prev.heart.failure + 
                     prev.niereninsuff + coronary.heart.disease, data = dat)

# Obtain summary of coefficients
coef_summary <- summary(cox_model)

hr <- coef_summary$coefficients[, "exp(coef)"][paste0(biomarker_names, "_log_std")]
p_value <- coef_summary$coefficients[, "Pr(>|z|)"][paste0(biomarker_names, "_log_std")]

# Combine results into a data frame
volcano_data <- rbind(volcano_data, data.frame(biomarker = biomarker_names, HR = hr, P_Value = p_value))

# Display the dataframe
print(volcano_data)

new_names <- c("D-dimer", "ALAT", "GDF-15", "IL-6", "NT-proBNP", "hsTropT")

# Replace biomarker names in coef_data
volcano_data$biomarker <- ifelse(volcano_data$biomarker == "ddi2h", new_names[1],
                                 ifelse(volcano_data$biomarker == "alat", new_names[2],
                                        ifelse(volcano_data$biomarker == "gdf.15", new_names[3],
                                               ifelse(volcano_data$biomarker == "il6", new_names[4],
                                                      ifelse(volcano_data$biomarker == "probnpii", new_names[5],
                                                             ifelse(volcano_data$biomarker == "tnt.hs", new_names[6],
                                                                    volcano_data$biomarker))))))

library(ggplot2)
library(ggrepel)

ggplot(volcano_data, aes(x = HR, y = -log10(P_Value), label = biomarker)) +
  geom_point() +
  labs(x = "Standardized hazard ratio", y = "-Log10 P-value") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "red") +
  theme_bw() + 
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 9, color = "black"),
        axis.text.x = element_text(hjust = 0.5)) +
  coord_cartesian(xlim = c(0.5, 2.0)) +
  scale_x_continuous(trans = "log10", breaks = c(0.5, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2.0), labels = c(0.5, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2.0)) +
  geom_text_repel(size = 2.75)

########################## VARIABLE IMPORTANCE PLOTS ###########################

library(dplyr)

features <- dat %>% 
  select(ends_with("_log_std"), age.bl, pat.sex, bmi, current.smoker, rr.sys.liegend,
         prev.diabetes, prev.stroke.tia, prev.heart.failure, prev.niereninsuff,
         coronary.heart.disease, timeto1.death, death.cardiac)

features <- na.omit(features)

cox_model <- coxph(Surv(timeto1.death, death.cardiac) ~ ddi2h_log_std + 
                     alat_log_std + gdf.15_log_std + il6_log_std + probnpii_log_std + 
                     tnt.hs_log_std + age.bl + pat.sex + bmi + current.smoker + 
                     rr.sys.liegend + prev.diabetes + prev.stroke.tia + prev.heart.failure + 
                     prev.niereninsuff + coronary.heart.disease, data = data)
summary(cox_model)

# Extract coefficient estimates and standard errors
coef_est <- coef(cox_model)
coef_se <- sqrt(diag(vcov(cox_model)))

# Calculate partial chi-squared statistic for each predictor
partial_chisq <- (coef_est / coef_se)^2

# Extract the total degrees of freedom for the model
num_params <- length(coef(cox_model))
df_total <- num_params

# Adjust partial chi-squared statistic by subtracting total degrees of freedom
partial_chisq_adj <- partial_chisq - df_total

# Display the results
result <- data.frame(Variable = names(coef_est), Partial_ChiSq = partial_chisq)
print(result)

# Order the result by partial chi-squared statistic in descending order
result <- result[order(-result$Partial_ChiSq), ]

library(ggplot2)

# Define the desired label changes
label_changes <- c("ddi2h_log_std" = "D-dimer",
                   "alat_log_std" = "ALAT",
                   "gdf.15_log_std" = "GDF-15",
                   "il6_log_std" = "IL-6",
                   "probnpii_log_std" = "NT-proBNP",
                   "tnt.hs_log_std" = "hsTropT",
                   "age.bl" = "Age",
                   "pat.sex" = "Sex",
                   "bmi" = "BMI",
                   "current.smoker" = "Smoker",
                   "rr.sys.liegend" = "SBP",
                   "prev.diabetes" = "Diabetes",
                   "prev.stroke.tia" = "Prior stroke/TIA",
                   "prev.heart.failure" = "Heart failure",
                   "prev.niereninsuff" = "CKD",
                   "coronary.heart.disease" = "CAD")

# Create the dot plot
p <- ggplot(result, aes(x = Partial_ChiSq, y = reorder(Variable, Partial_ChiSq))) +
  geom_point(shape = 1, size = 2.5, color = "black") +  # Change shape to open circle and set color to black
  labs(x = "Partial χ2-df", size = 2, y = "", color = "black") +  # Fixed x-axis label and removed y-axis title
  scale_y_discrete(labels = label_changes) +  # Apply label changes to y-axis
  theme_bw() +  # Minimal theme
  theme(axis.text.x = element_text(size = 10),  # Adjust x-axis text size
        axis.text.y = element_text(size = 10),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(color = "black"),  # Set x-axis ticks color to black
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.major.y = element_line(linetype = "dotted", color = "grey")  # Specify major horizontal grid lines as dotted and grey
  )

print(p)

################################ FOREST PLOTS ##################################

# Combined model
library(randomForest)
library(survminer)

cox_model <- coxph(Surv(timeto1.death, death.cardiac) ~ ddi2h_log_std + 
                     alat_log_std + gdf.15_log_std + il6_log_std + probnpii_log_std + 
                     tnt.hs_log_std + age.bl + pat.sex + bmi + current.smoker + 
                     rr.sys.liegend + prev.diabetes + prev.stroke.tia + prev.heart.failure + 
                     prev.niereninsuff + coronary.heart.disease, data = dat)
summary(cox_model)

biomarkers <- c("ddi2h_log_std", "alat_log_std", "gdf.15_log_std", 
                "il6_log_std", "probnpii_log_std", "tnt.hs_log_std")

# Create a dataframe to store hazard ratios and 95% CIs
coef_data <- data.frame(
  Biomarker = biomarkers,
  HR = exp(coef(cox_model)[biomarkers]),             # Extracting hazard ratios
  Lower_CI = exp(confint(cox_model)[biomarkers, 1]),  # Extracting lower CIs
  Upper_CI = exp(confint(cox_model)[biomarkers, 2]),   # Extracting upper CIs
  p_value = summary(cox_model)$coefficients[biomarkers, "Pr(>|z|)"]  # Extracting p-values
)

# View the dataframe
print(coef_data)

new_names <- c("D-dimer", "ALAT", "GDF-15", "IL-6", "NT-proBNP", "hsTropT")

# Replace biomarker names in coef_data
coef_data$Biomarker <- ifelse(coef_data$Biomarker == "ddi2h_log_std", new_names[1],
                              ifelse(coef_data$Biomarker == "alat_log_std", new_names[2],
                                     ifelse(coef_data$Biomarker == "gdf.15_log_std", new_names[3],
                                            ifelse(coef_data$Biomarker == "il6_log_std", new_names[4],
                                                   ifelse(coef_data$Biomarker == "probnpii_log_std", new_names[5],
                                                          ifelse(coef_data$Biomarker == "tnt.hs_log_std", new_names[6],
                                                                 coef_data$Biomarker))))))
coef_data$p.signif <- ifelse(coef_data$p_value < 0.05, "*", "")

library(ggplot2)

# Create coefficient plot using ggplot2
coefficient_plot <- ggplot(coef_data, aes(x = Biomarker, y = HR, color = Biomarker)) +
  geom_point(size = 2.5) + 
  geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 0.3, linetype = "solid") + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") + 
  labs(x = "Biomarker", y = "Hazard ratio (95% CI) per 1 SD increase in biomarker levels", size=11, color = "black") + 
  theme_bw() +
  theme(legend.position="none",
        axis.title = element_text(face="bold"),
        axis.text.x = element_text(size=9, color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(size=9, color = "black")) +
  scale_y_continuous(limits = c(0.7, 1.7), breaks = seq(0.7, 1.7, by = 0.1)) +
  geom_text(aes(label = paste0(sprintf("%.2f", HR))), 
            hjust = -0.3, size = 3)
print(coefficient_plot)


################################################################################
# All-cause death
################################################################################

#-------------------------------------------------------------------------------
# Check whether there are missing "time" variables 
num_missing_timeto1.death <- sum(is.na(dat$timeto1.death))
print(num_missing_timeto1.death)
#-------------------------------------------------------------------------------

########################### COX REGRESSION MODELS ##############################

# Fit the age + sex adjusted Cox model
library(survival)

timeto1.death_years <- dat$timeto1.death / 365.25
median_follow_up_years <- median(timeto1.death_years, na.rm = TRUE)
iqr_25_75 <- quantile(timeto1.death_years, probs = c(0.25, 0.75), na.rm = TRUE)

# Create empty lists to store results from age+sex and multivariable adjusted models
results.age_summary <- list()
results.multi_summary <- list()
results.age <- list()
results.multi <- list()

# Age+sex adjusted: Iterate over each biomarker
for (biomarker in biomarkers) {
  
  # Fit the Cox model
  cox_model <- coxph(Surv(timeto1.death, death.any) ~
                       get(paste0(biomarker, "_log_std")) + age.bl + pat.sex, data = dat)
  
  # Store results in the list
  results.age_summary[[biomarker]] <- summary(cox_model)
  
  # Store results in the list for volcano plots
  results.age[[biomarker]] <- coef(summary(cox_model))[, "exp(coef)"]
}

# Multivariable: Iterate over each biomarker
for (biomarker in biomarkers) {
  
  # Fit the Cox model
  cox_model <- coxph(Surv(timeto1.death, death.any) ~
                       get(paste0(biomarker, "_log_std")) + age.bl + pat.sex + bmi +
                       current.smoker + rr.sys.liegend + prev.diabetes +
                       prev.stroke.tia + prev.heart.failure + prev.niereninsuff +
                       coronary.heart.disease, data = dat)
  
  # Store results in the list
  results.multi_summary[[biomarker]] <- summary(cox_model)
  
  # Store results in the list for volcano plots
  results.multi[[biomarker]] <- coef(summary(cox_model))[, "exp(coef)"]
}

# Combined model 2 (all biomarkers)
cox_model <- coxph(Surv(timeto1.death, death.any) ~ ang2_log_std + ddi2h_log_std + 
                     crcl_cg_log_std + cysc_log_std + alat_log_std + gdf.15_log_std + crphs_log_std + 
                     igfbp7_log_std + il6_log_std + probnpii_log_std + opn_log_std + tnt.hs_log_std + 
                     age.bl + pat.sex + bmi + current.smoker + rr.sys.liegend + 
                     prev.diabetes + prev.stroke.tia + prev.heart.failure + 
                     prev.niereninsuff + coronary.heart.disease, data = dat)
summary(cox_model)

################### BACKWARD BIOMARKER SELECTION ###############################

# Perform a backward and forward selection of biomarkers
# Create a new dataset excluding patients with missing values in specified variables
# Exclude patients with missing variables
library(dplyr)

data <- dat %>% 
  select(ends_with("_log_std"), age.bl, pat.sex, bmi, current.smoker, rr.sys.liegend,
         prev.diabetes, prev.stroke.tia, prev.heart.failure, prev.niereninsuff,
         coronary.heart.disease, timeto1.death, death.any)

data <- na.omit(data)

# Fit the Cox proportional hazards model with the filtered dataset
full_model <- coxph(Surv(timeto1.death, death.any) ~ ang2_log_std + ddi2h_log_std + 
                      crcl_cg_log_std + cysc_log_std + alat_log_std + gdf.15_log_std + 
                      crphs_log_std + igfbp7_log_std + il6_log_std + probnpii_log_std + 
                      opn_log_std + tnt.hs_log_std + age.bl + pat.sex + bmi + 
                      current.smoker + rr.sys.liegend + prev.diabetes + prev.stroke.tia + 
                      prev.heart.failure + prev.niereninsuff + coronary.heart.disease, 
                    data = data)

summary(full_model)

# Define the range of models for stepwise selection (only consider changes in biomarkers)
scope <- list(lower = ~ age.bl + pat.sex + bmi + 
                current.smoker + rr.sys.liegend + prev.diabetes + prev.stroke.tia + 
                prev.heart.failure + prev.niereninsuff + coronary.heart.disease,
              upper = ~ ang2_log_std + ddi2h_log_std + 
                crcl_cg_log_std + cysc_log_std + alat_log_std + gdf.15_log_std + 
                crphs_log_std + igfbp7_log_std + il6_log_std + probnpii_log_std + 
                opn_log_std + tnt.hs_log_std + age.bl + pat.sex + bmi + 
                current.smoker + rr.sys.liegend + prev.diabetes + prev.stroke.tia + 
                prev.heart.failure + prev.niereninsuff + coronary.heart.disease)

# Perform stepwise selection for model refinement
reduced_model <- step(full_model, direction = "backward", scope = scope)

# Perform stepwise selection using stepAIC from MASS package
library(MASS) 
reduced_model <- stepAIC(full_model, scope = scope, direction = "backward")

# View the final selected model
summary(reduced_model)

# Use the biomarkers included in the step_model_filtered for the combined model
# Exclude CRP
cox_model <- coxph(Surv(timeto1.death, death.any) ~ ddi2h_log_std + 
                     alat_log_std + gdf.15_log_std + 
                     igfbp7_log_std + il6_log_std + probnpii_log_std + opn_log_std +
                     tnt.hs_log_std + age.bl + pat.sex + bmi + current.smoker + rr.sys.liegend + 
                     prev.diabetes + prev.stroke.tia + prev.heart.failure + 
                     prev.niereninsuff + coronary.heart.disease, data = dat)
summary(cox_model)

################ AUC FOR BASE MODEL AND BASE MODEL + BIOMARKERS ################

# Load necessary libraries
library(survival)
library(timeROC)
library(pROC)

# Fit the Cox proportional hazards model on the complete cases
cox_model_base <- coxph(Surv(timeto1.death, death.any) ~ age.bl + pat.sex + bmi + 
                          current.smoker + rr.sys.liegend + 
                          prev.diabetes + prev.stroke.tia + prev.heart.failure + 
                          prev.niereninsuff + coronary.heart.disease, data = data)

cox_model_biomarkers <- coxph(Surv(timeto1.death, death.any) ~ ddi2h_log_std + 
                                alat_log_std + gdf.15_log_std + 
                                igfbp7_log_std + il6_log_std + probnpii_log_std + opn_log_std +
                                tnt.hs_log_std + age.bl + pat.sex + bmi + current.smoker + rr.sys.liegend + 
                                prev.diabetes + prev.stroke.tia + prev.heart.failure + 
                                prev.niereninsuff + coronary.heart.disease, data = data)

# Obtain the risk scores (linear predictors) for each model
risk_scores_base <- predict(cox_model_base, type = "lp", newdata = data)
risk_scores_biomarkers <- predict(cox_model_biomarkers, type = "lp", newdata = data)

# Define time points at which to calculate the AUC
time_points <- c(1, 2, 3) * 365.25  # 1, 2, and 3 years in days

# Calculate time-dependent AUC for the base model
roc_results_base <- timeROC(T = data$timeto1.death,
                            delta = data$death.any,
                            marker = risk_scores_base,
                            cause = 1,
                            weighting = "marginal",
                            times = time_points,
                            ROC = TRUE)

# Calculate time-dependent AUC for the biomarker model
roc_results_biomarkers <- timeROC(T = data$timeto1.death,
                                  delta = data$death.any,
                                  marker = risk_scores_biomarkers,
                                  cause = 1,
                                  weighting = "marginal",
                                  times = time_points,
                                  ROC = TRUE)

# Print AUC values for both models
print(roc_results_base$AUC)
print(roc_results_biomarkers$AUC)

# Extract predicted probabilities at the specific time point (1 year, 365.25 days)
pred_probs_base <- predict(cox_model_base, type = "risk", newdata = data)
pred_probs_biomarkers <- predict(cox_model_biomarkers, type = "risk", newdata = data)

# Use pROC to create ROC curves at the 1-year mark
roc_base <- roc(data$death.any, pred_probs_base, plot = FALSE, ci = TRUE, quiet = TRUE)
roc_biomarkers <- roc(data$death.any, pred_probs_biomarkers, plot = FALSE, ci = TRUE, quiet = TRUE)

# Perform DeLong's test to compare the AUCs
delong_test <- roc.test(roc_base, roc_biomarkers)
print(delong_test)

# Plot ROC Curves at different time points
plot(roc_results, time = time_points[1], col = "blue", title = "ROC Curve at 1 Year")
plot(roc_results, time = time_points[2], col = "red", title = "ROC Curve at 2 Years", add = TRUE)
plot(roc_results, time = time_points[3], col = "green", title = "ROC Curve at 3 Years", add = TRUE)
legend("bottomright", legend = c("1 Year", "2 Years", "3 Years"), col = c("blue", "red", "green"), lwd = 2)

############################### VOLCANO PLOTS ##################################

# Create an empty dataframe to store results
volcano_data <- data.frame(Biomarker = character(),
                           HR = numeric(),
                           P_Value = numeric(),
                           stringsAsFactors = FALSE)

# Iterate through each biomarker
for (biomarker in biomarkers) {
  
  # Fit the Cox model
  cox_model <- coxph(Surv(timeto1.death, death.any) ~
                       get(paste0(biomarker, "_log_std")) + age.bl + pat.sex + bmi +
                       current.smoker + rr.sys.liegend + prev.diabetes +
                       prev.stroke.tia + prev.heart.failure + prev.niereninsuff +
                       coronary.heart.disease, data = dat)
  
  # Extract HR and p-value for _log_std
  coef_summary <- coef(summary(cox_model))
  hr <- coef_summary["get(paste0(biomarker, \"_log_std\"))", "exp(coef)"]
  p_value <- coef_summary["get(paste0(biomarker, \"_log_std\"))", "Pr(>|z|)"]
  
  # Add results to the dataframe
  volcano_data <- rbind(volcano_data, data.frame(biomarker = biomarker, HR = hr, P_Value = p_value))
}

# Volcano plot from the combined model after stepwise biomarker selection
volcano_data <- data.frame(Biomarker = character(),
                           HR = numeric(),
                           P_Value = numeric(),
                           stringsAsFactors = FALSE)

biomarker_names <- c("ddi2h", "alat", "gdf.15", 
                     "igfbp7", "il6", "probnpii", "opn", "tnt.hs")

cox_model <- coxph(Surv(timeto1.death, death.any) ~ ddi2h_log_std + 
                     alat_log_std + gdf.15_log_std + 
                     igfbp7_log_std + il6_log_std + probnpii_log_std + opn_log_std +
                     tnt.hs_log_std + age.bl + pat.sex + bmi + current.smoker + rr.sys.liegend + 
                     prev.diabetes + prev.stroke.tia + prev.heart.failure + 
                     prev.niereninsuff + coronary.heart.disease, data = dat)

# Obtain summary of coefficients
coef_summary <- summary(cox_model)

hr <- coef_summary$coefficients[, "exp(coef)"][paste0(biomarker_names, "_log_std")]
p_value <- coef_summary$coefficients[, "Pr(>|z|)"][paste0(biomarker_names, "_log_std")]

# Combine results into a data frame
volcano_data <- rbind(volcano_data, data.frame(biomarker = biomarker_names, HR = hr, P_Value = p_value))

# Display the dataframe
print(volcano_data)

new_names <- c("D-dimer", "ALAT", "GDF-15", "IGFBP-7", "IL-6", "NT-proBNP", "OPN", "hsTropT")

# Replace biomarker names in coef_data
volcano_data$biomarker <- ifelse(volcano_data$biomarker == "ddi2h", new_names[1],
                                 ifelse(volcano_data$biomarker == "alat", new_names[2],
                                        ifelse(volcano_data$biomarker == "gdf.15", new_names[3],
                                               ifelse(volcano_data$biomarker == "igfbp7", new_names[4],
                                                      ifelse(volcano_data$biomarker == "il6", new_names[5],
                                                             ifelse(volcano_data$biomarker == "probnpii", new_names[6],
                                                                    ifelse(volcano_data$biomarker == "opn", new_names[7],
                                                                           ifelse(volcano_data$biomarker == "tnt.hs", new_names[8],
                                                                                  volcano_data$biomarker))))))))

library(ggplot2)
library(ggrepel)

ggplot(volcano_data, aes(x = HR, y = -log10(P_Value), label = biomarker)) +
  geom_point() +
  labs(x = "Standardized hazard ratio", y = "-Log10 P-value") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "red") +
  theme_bw() + 
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 9, color = "black"),
        axis.text.x = element_text(hjust = 0.5)) +
  coord_cartesian(xlim = c(0.5, 2.0)) +
  scale_x_continuous(trans = "log10", breaks = c(0.5, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2.0), labels = c(0.5, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2.0)) +
  geom_text_repel(size = 2.75)

########################## VARIABLE IMPORTANCE PLOTS ###########################

library(dplyr)

features <- dat %>% 
  select(ends_with("_log_std"), age.bl, pat.sex, bmi, current.smoker, rr.sys.liegend,
         prev.diabetes, prev.stroke.tia, prev.heart.failure, prev.niereninsuff,
         coronary.heart.disease, timeto1.death, death.any)

features <- na.omit(features)

cox_model <- coxph(Surv(timeto1.death, death.any) ~ ddi2h_log_std + 
                     alat_log_std + gdf.15_log_std + 
                     igfbp7_log_std + il6_log_std + probnpii_log_std + opn_log_std +
                     tnt.hs_log_std + age.bl + pat.sex + bmi + current.smoker + rr.sys.liegend + 
                     prev.diabetes + prev.stroke.tia + prev.heart.failure + 
                     prev.niereninsuff + coronary.heart.disease, data = features)
summary(cox_model)

# Extract coefficient estimates and standard errors
coef_est <- coef(cox_model)
coef_se <- sqrt(diag(vcov(cox_model)))

# Calculate partial chi-squared statistic for each predictor
partial_chisq <- (coef_est / coef_se)^2

# Extract the total degrees of freedom for the model
num_params <- length(coef(cox_model))
df_total <- num_params

# Adjust partial chi-squared statistic by subtracting total degrees of freedom
partial_chisq_adj <- partial_chisq - df_total

# Display the results
result <- data.frame(Variable = names(coef_est), Partial_ChiSq = partial_chisq)
print(result)

# Order the result by partial chi-squared statistic in descending order
result <- result[order(-result$Partial_ChiSq), ]

library(ggplot2)

# Define the desired label changes
label_changes <- c("ddi2h_log_std" = "D-dimer",
                   "alat_log_std" = "ALAT",
                   "gdf.15_log_std" = "GDF-15",
                   "igfbp7_log_std" = "IGFBP-7",
                   "il6_log_std" = "IL-6",
                   "probnpii_log_std" = "NT-proBNP",
                   "opn_log_std" = "OPN",
                   "tnt.hs_log_std" = "hsTropT",
                   "age.bl" = "Age",
                   "pat.sex" = "Sex",
                   "bmi" = "BMI",
                   "current.smoker" = "Smoker",
                   "rr.sys.liegend" = "SBP",
                   "prev.diabetes" = "Diabetes",
                   "prev.stroke.tia" = "Prior stroke/TIA",
                   "prev.heart.failure" = "Heart failure",
                   "prev.niereninsuff" = "CKD",
                   "coronary.heart.disease" = "CAD")

# Create the dot plot
p <- ggplot(result, aes(x = Partial_ChiSq, y = reorder(Variable, Partial_ChiSq))) +
  geom_point(shape = 1, size = 2.5, color = "black") +  # Change shape to open circle and set color to black
  labs(x = "Partial χ2-df", size = 2, y = "", color = "black") +  # Fixed x-axis label and removed y-axis title
  scale_y_discrete(labels = label_changes) +  # Apply label changes to y-axis
  theme_bw() +  # Minimal theme
  theme(axis.text.x = element_text(size = 10),  # Adjust x-axis text size
        axis.text.y = element_text(size = 10),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(color = "black"),  # Set x-axis ticks color to black
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.major.y = element_line(linetype = "dotted", color = "grey")  # Specify major horizontal grid lines as dotted and grey
  )

print(p)

################################ FOREST PLOTS ##################################

# Combined model
library(randomForest)
library(survminer)

cox_model <- coxph(Surv(timeto1.death, death.any) ~ ddi2h_log_std + 
                     alat_log_std + gdf.15_log_std + 
                     igfbp7_log_std + il6_log_std + probnpii_log_std + opn_log_std +
                     tnt.hs_log_std + age.bl + pat.sex + bmi + current.smoker + rr.sys.liegend + 
                     prev.diabetes + prev.stroke.tia + prev.heart.failure + 
                     prev.niereninsuff + coronary.heart.disease, data = dat)
summary(cox_model)

biomarkers <- c("ddi2h_log_std", "alat_log_std", "gdf.15_log_std", "igfbp7_log_std", 
                "il6_log_std", "probnpii_log_std", "opn_log_std", "tnt.hs_log_std")

# Create a dataframe to store hazard ratios and 95% CIs
coef_data <- data.frame(
  Biomarker = biomarkers,
  HR = exp(coef(cox_model)[biomarkers]),             # Extracting hazard ratios
  Lower_CI = exp(confint(cox_model)[biomarkers, 1]),  # Extracting lower CIs
  Upper_CI = exp(confint(cox_model)[biomarkers, 2]),   # Extracting upper CIs
  p_value = summary(cox_model)$coefficients[biomarkers, "Pr(>|z|)"]  # Extracting p-values
)

# View the dataframe
print(coef_data)

new_names <- c("D-dimer", "ALAT", "GDF-15", "IGFBP-7", "IL-6", "NT-proBNP", "OPN", "hsTropT")

# Replace biomarker names in coef_data
coef_data$Biomarker <- ifelse(coef_data$Biomarker == "ddi2h_log_std", new_names[1],
                              ifelse(coef_data$Biomarker == "alat_log_std", new_names[2],
                                     ifelse(coef_data$Biomarker == "gdf.15_log_std", new_names[3],
                                            ifelse(coef_data$Biomarker == "igfbp7_log_std", new_names[4],
                                                   ifelse(coef_data$Biomarker == "il6_log_std", new_names[5],
                                                          ifelse(coef_data$Biomarker == "probnpii_log_std", new_names[6],
                                                                 ifelse(coef_data$Biomarker == "opn_log_std", new_names[7],
                                                                        ifelse(coef_data$Biomarker == "tnt.hs_log_std", new_names[8],
                                                                               coef_data$Biomarker))))))))

coef_data$p.signif <- ifelse(coef_data$p_value < 0.05, "*", "")

library(ggplot2)

# Create coefficient plot using ggplot2
coefficient_plot <- ggplot(coef_data, aes(x = Biomarker, y = HR, color = Biomarker)) +
  geom_point(size = 2.5) + 
  geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 0.3, linetype = "solid") + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") + 
  labs(x = "Biomarker", y = "Hazard ratio (95% CI) per 1 SD increase in biomarker levels", size=11, color = "black") + 
  theme_bw() +
  theme(legend.position="none",
        axis.title = element_text(face="bold"),
        axis.text.x = element_text(size=9, color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(size=9, color = "black")) +
  scale_y_continuous(limits = c(0.7, 1.7), breaks = seq(0.7, 1.7, by = 0.1)) +
  geom_text(aes(label = paste0(sprintf("%.2f", HR))), 
            hjust = -0.3, size = 3)
print(coefficient_plot)


################################################################################
# Composite major and clinically relevant non-major bleeding
################################################################################

#-------------------------------------------------------------------------------
# Check whether there are missing "time" variables 
num_missing_timeto1.anybleed <- sum(is.na(dat$timeto1.anybleed))
print(num_missing_timeto1.anybleed)
#-------------------------------------------------------------------------------

################### COX REGRESSION MODELS ######################################

# Fit the age + sex adjusted Cox model
library(survival)

# Create empty lists to store results from age+sex and multivariable adjusted models
results.age_summary <- list()
results.multi_summary <- list()
results.age <- list()
results.multi <- list()

# Age+sex adjusted: Iterate over each biomarker
for (biomarker in biomarkers) {
  
  # Fit the Cox model
  cox_model <- coxph(Surv(timeto1.anybleed, bleed.any) ~
                       get(paste0(biomarker, "_log_std")) + age.bl + pat.sex, data = dat)
  
  # Store results in the list
  results.age_summary[[biomarker]] <- summary(cox_model)
  
  # Store results in the list for volcano plots
  results.age[[biomarker]] <- coef(summary(cox_model))[, "exp(coef)"]
}

# Multivariable: Iterate over each biomarker
for (biomarker in biomarkers) {
  
  # Fit the Cox model
  cox_model <- coxph(Surv(timeto1.anybleed, bleed.any) ~
                       get(paste0(biomarker, "_log_std")) + age.bl + pat.sex + bmi +
                       current.smoker + rr.sys.liegend + prev.diabetes +
                       prev.stroke.tia + prev.heart.failure + prev.niereninsuff +
                       coronary.heart.disease, data = dat)
  
  # Store results in the list
  results.multi_summary[[biomarker]] <- summary(cox_model)
  
  # Store results in the list for volcano plots
  results.multi[[biomarker]] <- coef(summary(cox_model))[, "exp(coef)"]
}

# Combined model 1 (significant biomarkers)
cox_model <- coxph(Surv(timeto1.anybleed, bleed.any) ~ ang2_log_std + ddi2h_log_std + 
                     crcl_cg_log_std + gdf.15_log_std + 
                     igfbp7_log_std + il6_log_std + probnpii_log_std + opn_log_std + tnt.hs_log_std + 
                     age.bl + pat.sex + bmi + current.smoker + rr.sys.liegend + 
                     prev.diabetes + prev.stroke.tia + prev.heart.failure + 
                     prev.niereninsuff + coronary.heart.disease, data = dat)
summary(cox_model)

# Combined model 2 (all biomarkers)
cox_model <- coxph(Surv(timeto1.anybleed, bleed.any) ~ ang2_log_std + ddi2h_log_std + 
                     crcl_cg_log_std + cysc_log_std + alat_log_std + gdf.15_log_std + crphs_log_std + 
                     igfbp7_log_std + il6_log_std + probnpii_log_std + opn_log_std + tnt.hs_log_std + 
                     age.bl + pat.sex + bmi + current.smoker + rr.sys.liegend + 
                     prev.diabetes + prev.stroke.tia + prev.heart.failure + 
                     prev.niereninsuff + coronary.heart.disease, data = dat)
summary(cox_model)

################### BACKWARD BIOMARKER SELECTION ###############################

# Create a new dataset excluding patients with missing values in specified variables
# Exclude patients with missing variables
library(dplyr)

data <- dat %>% 
  select(ends_with("_log_std"), age.bl, pat.sex, bmi, current.smoker, rr.sys.liegend,
         prev.diabetes, prev.stroke.tia, prev.heart.failure, prev.niereninsuff,
         coronary.heart.disease, timeto1.anybleed, bleed.any)

data <- na.omit(data)

# Fit the Cox proportional hazards model with the filtered dataset
full_model <- coxph(Surv(timeto1.anybleed, bleed.any) ~ ang2_log_std + ddi2h_log_std + 
                      cysc_log_std + alat_log_std + gdf.15_log_std + crphs_log_std + 
                      igfbp7_log_std + il6_log_std + probnpii_log_std + opn_log_std + tnt.hs_log_std + 
                      crcl_cg_log_std + age.bl + pat.sex + bmi + current.smoker + rr.sys.liegend + 
                      prev.diabetes + prev.stroke.tia + prev.heart.failure + crcl_cg_log_std +
                      prev.niereninsuff + coronary.heart.disease, 
                    data = data)

summary(full_model)
vif(full_model)

# Define the range of models for stepwise selection (only consider changes in biomarkers)
scope <- list(lower = ~ age.bl + pat.sex + bmi + 
                current.smoker + rr.sys.liegend + prev.diabetes + prev.stroke.tia + 
                prev.heart.failure + prev.niereninsuff + coronary.heart.disease,
              upper = ~ ang2_log_std + ddi2h_log_std + 
                cysc_log_std + alat_log_std + gdf.15_log_std + crphs_log_std + 
                igfbp7_log_std + il6_log_std + probnpii_log_std + opn_log_std + tnt.hs_log_std + 
                crcl_cg_log_std + age.bl + pat.sex + bmi + 
                current.smoker + rr.sys.liegend + prev.diabetes + prev.stroke.tia + 
                prev.heart.failure + prev.niereninsuff + coronary.heart.disease)

# Perform stepwise selection for model refinement
reduced_model <- step(full_model, direction = "backward", scope = scope)

# Perform stepwise selection using stepAIC from MASS package
library(MASS) 
reduced_model <- stepAIC(full_model, scope = scope, direction = "backward")

# View the final selected model
summary(reduced_model)

# Exclude CRP
cox_model <- coxph(Surv(timeto1.anybleed, bleed.any) ~ gdf.15_log_std + 
                     igfbp7_log_std + il6_log_std + probnpii_log_std + 
                     age.bl + pat.sex + bmi + current.smoker + 
                     rr.sys.liegend + prev.diabetes + prev.stroke.tia + prev.heart.failure + 
                     prev.niereninsuff + coronary.heart.disease, data = data)

summary(cox_model)

################ AUC FOR BASE MODEL AND BASE MODEL + BIOMARKERS ################

library(survival)
library(timeROC)
library(pROC)

# Fit the Cox proportional hazards model on the complete cases
cox_model_base <- coxph(Surv(timeto1.anybleed, bleed.any) ~ age.bl + pat.sex + bmi + 
                          current.smoker + rr.sys.liegend + 
                          prev.diabetes + prev.stroke.tia + prev.heart.failure + 
                          prev.niereninsuff + coronary.heart.disease, data = data)

cox_model_biomarkers <- coxph(Surv(timeto1.anybleed, bleed.any) ~ gdf.15_log_std + 
                                igfbp7_log_std + il6_log_std + probnpii_log_std + 
                                age.bl + pat.sex + bmi + current.smoker + 
                                rr.sys.liegend + prev.diabetes + prev.stroke.tia + prev.heart.failure + 
                                prev.niereninsuff + coronary.heart.disease, data = data)

# Obtain the risk scores (linear predictors) for each model
risk_scores_base <- predict(cox_model_base, type = "lp", newdata = data)
risk_scores_biomarkers <- predict(cox_model_biomarkers, type = "lp", newdata = data)

# Define time points at which to calculate the AUC
time_points <- c(1, 2, 3) * 365.25  # 1, 2, and 3 years in days

# Calculate time-dependent AUC for the base model
roc_results_base <- timeROC(T = data$timeto1.anybleed,
                            delta = data$bleed.any,
                            marker = risk_scores_base,
                            cause = 1,
                            weighting = "marginal",
                            times = time_points,
                            ROC = TRUE)

# Calculate time-dependent AUC for the biomarker model
roc_results_biomarkers <- timeROC(T = data$timeto1.anybleed,
                                  delta = data$bleed.any,
                                  marker = risk_scores_biomarkers,
                                  cause = 1,
                                  weighting = "marginal",
                                  times = time_points,
                                  ROC = TRUE)

# Print AUC values for both models
print(roc_results_base$AUC)
print(roc_results_biomarkers$AUC)

# Extract predicted probabilities at the specific time point (1 year, 365.25 days)
pred_probs_base <- predict(cox_model_base, type = "risk", newdata = data)
pred_probs_biomarkers <- predict(cox_model_biomarkers, type = "risk", newdata = data)

# Use pROC to create ROC curves at the 1-year mark
roc_base <- roc(data$bleed.any, pred_probs_base, plot = FALSE, ci = TRUE, quiet = TRUE)
roc_biomarkers <- roc(data$bleed.any, pred_probs_biomarkers, plot = FALSE, ci = TRUE, quiet = TRUE)

# Perform DeLong's test to compare the AUCs
delong_test <- roc.test(roc_base, roc_biomarkers)
print(delong_test)

# Plot ROC Curves at different time points
plot(roc_results, time = time_points[1], col = "blue", title = "ROC Curve at 1 Year")
plot(roc_results, time = time_points[2], col = "red", title = "ROC Curve at 2 Years", add = TRUE)
plot(roc_results, time = time_points[3], col = "green", title = "ROC Curve at 3 Years", add = TRUE)
legend("bottomright", legend = c("1 Year", "2 Years", "3 Years"), col = c("blue", "red", "green"), lwd = 2)


################################ VOLCANO PLOTS #################################

# Create an empty dataframe to store results
volcano_data <- data.frame(Biomarker = character(),
                           HR = numeric(),
                           P_Value = numeric(),
                           stringsAsFactors = FALSE)

# Iterate through each biomarker
for (biomarker in biomarkers) {
  
  # Fit the Cox model
  cox_model <- coxph(Surv(timeto1.anybleed, bleed.any) ~
                       get(paste0(biomarker, "_log_std")) + age.bl + pat.sex + bmi +
                       current.smoker + rr.sys.liegend + prev.diabetes +
                       prev.stroke.tia + prev.heart.failure + prev.niereninsuff +
                       coronary.heart.disease, data = dat)
  
  # Extract HR and p-value for _log_std
  coef_summary <- coef(summary(cox_model))
  hr <- coef_summary["get(paste0(biomarker, \"_log_std\"))", "exp(coef)"]
  p_value <- coef_summary["get(paste0(biomarker, \"_log_std\"))", "Pr(>|z|)"]
  
  # Add results to the dataframe
  volcano_data <- rbind(volcano_data, data.frame(biomarker = biomarker, HR = hr, P_Value = p_value))
}

# Volcano plot from the combined model after stepwise biomarker selection
volcano_data <- data.frame(Biomarker = character(),
                           HR = numeric(),
                           P_Value = numeric(),
                           stringsAsFactors = FALSE)

biomarker_names <- c("gdf.15", "il6", 
                     "igfbp7", "probnpii")

cox_model <- coxph(Surv(timeto1.anybleed, bleed.any) ~ gdf.15_log_std + 
                     igfbp7_log_std + il6_log_std + probnpii_log_std + 
                     age.bl + pat.sex + bmi + current.smoker + 
                     rr.sys.liegend + prev.diabetes + prev.stroke.tia + prev.heart.failure + 
                     prev.niereninsuff + coronary.heart.disease, data = dat)

# Obtain summary of coefficients
coef_summary <- summary(cox_model)

hr <- coef_summary$coefficients[, "exp(coef)"][paste0(biomarker_names, "_log_std")]
p_value <- coef_summary$coefficients[, "Pr(>|z|)"][paste0(biomarker_names, "_log_std")]

# Combine results into a data frame
volcano_data <- rbind(volcano_data, data.frame(biomarker = biomarker_names, HR = hr, P_Value = p_value))

# Display the dataframe
print(volcano_data)

new_names <- c("GDF-15", "IL-6", "IGFBP-7", "NT-proBNP")

# Replace biomarker names in coef_data
volcano_data$biomarker <- ifelse(volcano_data$biomarker == "gdf.15", new_names[1],
                                 ifelse(volcano_data$biomarker == "il6", new_names[2],
                                        ifelse(volcano_data$biomarker == "igfbp7", new_names[3],
                                               ifelse(volcano_data$biomarker == "probnpii", new_names[4],
                                                      volcano_data$biomarker))))

library(ggplot2)
library(ggrepel)

ggplot(volcano_data, aes(x = HR, y = -log10(P_Value), label = biomarker)) +
  geom_point() +
  labs(x = "Standardized hazard ratio", y = "-Log10 P-value") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "red") +
  theme_bw() + 
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 9, color = "black"),
        axis.text.x = element_text(hjust = 0.5)) +
  coord_cartesian(xlim = c(0.5, 2.0)) +
  scale_x_continuous(trans = "log10", breaks = c(0.5, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2.0), labels = c(0.5, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2.0)) +
  geom_text_repel(size = 2.75)

########################## VARIABLE IMPORTANCE PLOTS ###########################

library(dplyr)

features <- dat %>% 
  select(ends_with("_log_std"), age.bl, pat.sex, bmi, current.smoker, rr.sys.liegend,
         prev.diabetes, prev.stroke.tia, prev.heart.failure, prev.niereninsuff,
         coronary.heart.disease, timeto1.anybleed, bleed.any)

features <- na.omit(features)

cox_model <- coxph(Surv(timeto1.anybleed, bleed.any) ~ gdf.15_log_std + 
                     igfbp7_log_std + il6_log_std + probnpii_log_std + 
                     age.bl + pat.sex + bmi + current.smoker + 
                     rr.sys.liegend + prev.diabetes + prev.stroke.tia + prev.heart.failure + 
                     prev.niereninsuff + coronary.heart.disease, data = data)
summary(cox_model)

# Extract coefficient estimates and standard errors
coef_est <- coef(cox_model)
coef_se <- sqrt(diag(vcov(cox_model)))

# Calculate partial chi-squared statistic for each predictor
partial_chisq <- (coef_est / coef_se)^2

# Extract the total degrees of freedom for the model
num_params <- length(coef(cox_model))
df_total <- num_params

# Adjust partial chi-squared statistic by subtracting total degrees of freedom
partial_chisq_adj <- partial_chisq - df_total

# Display the results
result <- data.frame(Variable = names(coef_est), Partial_ChiSq = partial_chisq)
print(result)

# Order the result by partial chi-squared statistic in descending order
result <- result[order(-result$Partial_ChiSq), ]

library(ggplot2)

# Define the desired label changes
label_changes <- c("gdf.15_log_std" = "GDF-15",
                   "igfbp7_log_std" = "IGFBP-7",
                   "il6_log_std" = "IL-6",
                   "probnpii_log_std" = "NT-proBNP",
                   "age.bl" = "Age",
                   "pat.sex" = "Sex",
                   "bmi" = "BMI",
                   "current.smoker" = "Smoker",
                   "rr.sys.liegend" = "SBP",
                   "prev.diabetes" = "Diabetes",
                   "prev.stroke.tia" = "Prior stroke/TIA",
                   "prev.heart.failure" = "Heart failure",
                   "prev.niereninsuff" = "CKD",
                   "coronary.heart.disease" = "CAD")

# Create the dot plot
p <- ggplot(result, aes(x = Partial_ChiSq, y = reorder(Variable, Partial_ChiSq))) +
  geom_point(shape = 1, size = 2.5, color = "black") +  # Change shape to open circle and set color to black
  labs(x = "Partial χ2-df", size = 2, y = "", color = "black") +  # Fixed x-axis label and removed y-axis title
  scale_y_discrete(labels = label_changes) +  # Apply label changes to y-axis
  theme_bw() +  # Minimal theme
  theme(axis.text.x = element_text(size = 10),  # Adjust x-axis text size
        axis.text.y = element_text(size = 10),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(color = "black"),  # Set x-axis ticks color to black
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.major.y = element_line(linetype = "dotted", color = "grey")  # Specify major horizontal grid lines as dotted and grey
  )

print(p)

############################## FOREST PLOTS ####################################

# Combined model
library(randomForest)
library(survminer)

cox_model <- coxph(Surv(timeto1.anybleed, bleed.any) ~ gdf.15_log_std + 
                     igfbp7_log_std + il6_log_std + probnpii_log_std + 
                     age.bl + pat.sex + bmi + current.smoker + 
                     rr.sys.liegend + prev.diabetes + prev.stroke.tia + prev.heart.failure + 
                     prev.niereninsuff + coronary.heart.disease, data = data)
summary(cox_model)

biomarkers <- c("gdf.15_log_std", "igfbp7_log_std", "il6_log_std", "probnpii_log_std")

# Create a dataframe to store hazard ratios and 95% CIs
coef_data <- data.frame(
  Biomarker = biomarkers,
  HR = exp(coef(cox_model)[biomarkers]),             # Extracting hazard ratios
  Lower_CI = exp(confint(cox_model)[biomarkers, 1]),  # Extracting lower CIs
  Upper_CI = exp(confint(cox_model)[biomarkers, 2]),   # Extracting upper CIs
  p_value = summary(cox_model)$coefficients[biomarkers, "Pr(>|z|)"]  # Extracting p-values
)

# View the dataframe
print(coef_data)

new_names <- c("GDF-15", "IGFBP-7", "IL-6", "NT-proBNP")

# Replace biomarker names in coef_data
coef_data$Biomarker <- ifelse(coef_data$Biomarker == "gdf.15_log_std", new_names[1],
                              ifelse(coef_data$Biomarker == "igfbp7_log_std", new_names[2],
                                     ifelse(coef_data$Biomarker == "il6_log_std", new_names[3],
                                            ifelse(coef_data$Biomarker == "probnpii_log_std", new_names[4],
                                                   coef_data$Biomarker))))

coef_data$p.signif <- ifelse(coef_data$p_value < 0.05, "*", "")

library(ggplot2)

# Create coefficient plot using ggplot2
coefficient_plot <- ggplot(coef_data, aes(x = Biomarker, y = HR, color = Biomarker)) +
  geom_point(size = 2.5) + 
  geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 0.3, linetype = "solid") + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") + 
  labs(x = "Biomarker", y = "Hazard ratio (95% CI) per 1 SD increase in biomarker levels", size=11, color = "black") + 
  theme_bw() +
  theme(legend.position="none",
        axis.title = element_text(face="bold"),
        axis.text.x = element_text(size=9, color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(size=9, color = "black")) +
  scale_y_continuous(limits = c(0.8, 1.5), breaks = seq(0.8, 1.5, by = 0.1)) +
  geom_text(aes(label = paste0(sprintf("%.2f", HR))), 
            hjust = -0.3, size = 3)
print(coefficient_plot)


################################################################################
# Clinically relevant non-major bleeding
################################################################################

#-------------------------------------------------------------------------------
# Check whether there are missing "time" variables 
num_missing_timeto1.minor.bleed <- sum(is.na(dat$timeto1.minor.bleed))
print(num_missing_timeto1.minor.bleed)
#-------------------------------------------------------------------------------

########################### COX REGRESSION MODELS ##############################

# Fit the age + sex adjusted Cox model
library(survival)

# Create empty lists to store results from age+sex and multivariable adjusted models
results.age_summary <- list()
results.multi_summary <- list()
results.age <- list()
results.multi <- list()

# Age+sex adjusted: Iterate over each biomarker
for (biomarker in biomarkers) {
  
  # Fit the Cox model
  cox_model <- coxph(Surv(timeto1.minor.bleed, minor.bleed) ~
                       get(paste0(biomarker, "_log_std")) + age.bl + pat.sex, data = dat)
  
  # Store results in the list
  results.age_summary[[biomarker]] <- summary(cox_model)
  
  # Store results in the list for volcano plots
  results.age[[biomarker]] <- coef(summary(cox_model))[, "exp(coef)"]
}

# Multivariable: Iterate over each biomarker
for (biomarker in biomarkers) {
  
  # Fit the Cox model
  cox_model <- coxph(Surv(timeto1.minor.bleed, minor.bleed) ~
                       get(paste0(biomarker, "_log_std")) + age.bl + pat.sex + bmi +
                       current.smoker + rr.sys.liegend + prev.diabetes +
                       prev.stroke.tia + prev.heart.failure + prev.niereninsuff +
                       coronary.heart.disease, data = dat)
  
  # Store results in the list
  results.multi_summary[[biomarker]] <- summary(cox_model)
  
  # Store results in the list for volcano plots
  results.multi[[biomarker]] <- coef(summary(cox_model))[, "exp(coef)"]
}

# Combined model 1 (significant biomarkers)
cox_model <- coxph(Surv(timeto1.minor.bleed, minor.bleed) ~ ang2_log_std + ddi2h_log_std + 
                    cysc_log_std + gdf.15_log_std + 
                     igfbp7_log_std + il6_log_std + probnpii_log_std + opn_log_std + tnt.hs_log_std + 
                     age.bl + pat.sex + bmi + current.smoker + rr.sys.liegend + 
                     prev.diabetes + prev.stroke.tia + prev.heart.failure + 
                     prev.niereninsuff + coronary.heart.disease, data = dat)
summary(cox_model)

# Combined model 2 (all biomarkers)
cox_model <- coxph(Surv(timeto1.minor.bleed, minor.bleed) ~ ang2_log_std + ddi2h_log_std + 
                     crcl_cg_log_std + cysc_log_std + alat_log_std + gdf.15_log_std + crphs_log_std + 
                     igfbp7_log_std + il6_log_std + probnpii_log_std + opn_log_std + tnt.hs_log_std + 
                     age.bl + pat.sex + bmi + current.smoker + rr.sys.liegend + 
                     prev.diabetes + prev.stroke.tia + prev.heart.failure + 
                     prev.niereninsuff + coronary.heart.disease, data = dat)
summary(cox_model)

######################## BACKWARD BIOMARKER SELECTION ##########################

# Create a new dataset excluding patients with missing values in specified variables
# Exclude patients with missing variables
library(dplyr)

data <- dat %>% 
  select(ends_with("_log_std"), age.bl, pat.sex, bmi, current.smoker, rr.sys.liegend,
         prev.diabetes, prev.stroke.tia, prev.heart.failure, prev.niereninsuff,
         coronary.heart.disease, timeto1.minor.bleed, minor.bleed)

data <- na.omit(data)

# Fit the Cox proportional hazards model with the filtered dataset
full_model <- coxph(Surv(timeto1.minor.bleed, minor.bleed) ~ ang2_log_std + ddi2h_log_std + 
                      crcl_cg_log_std + cysc_log_std + alat_log_std + gdf.15_log_std + 
                      crphs_log_std + igfbp7_log_std + il6_log_std + probnpii_log_std + 
                      opn_log_std + tnt.hs_log_std + age.bl + pat.sex + bmi + 
                      current.smoker + rr.sys.liegend + prev.diabetes + prev.stroke.tia + 
                      prev.heart.failure + prev.niereninsuff + coronary.heart.disease, 
                    data = data)

summary(full_model)

# Define the range of models for stepwise selection (only consider changes in biomarkers)
scope <- list(lower = ~ age.bl + pat.sex + bmi + 
                current.smoker + rr.sys.liegend + prev.diabetes + prev.stroke.tia + 
                prev.heart.failure + prev.niereninsuff + coronary.heart.disease,
              upper = ~ ang2_log_std + ddi2h_log_std + 
                crcl_cg_log_std + cysc_log_std + alat_log_std + gdf.15_log_std + 
                crphs_log_std + igfbp7_log_std + il6_log_std + probnpii_log_std + 
                opn_log_std + tnt.hs_log_std + age.bl + pat.sex + bmi + 
                current.smoker + rr.sys.liegend + prev.diabetes + prev.stroke.tia + 
                prev.heart.failure + prev.niereninsuff + coronary.heart.disease)

# Perform stepwise selection for model refinement
reduced_model <- step(full_model, direction = "backward", scope = scope)

# Perform stepwise selection using stepAIC from MASS package
library(MASS) 
reduced_model <- stepAIC(full_model, scope = scope, direction = "backward")

# View the final selected model
summary(reduced_model)

# Use the biomarkers included in the step_model_filtered for the combined model
cox_model <- coxph(Surv(timeto1.minor.bleed, minor.bleed) ~ gdf.15_log_std + 
                     il6_log_std + probnpii_log_std + 
                     age.bl + pat.sex + bmi + current.smoker + rr.sys.liegend + 
                     prev.diabetes + prev.stroke.tia + prev.heart.failure + 
                     prev.niereninsuff + coronary.heart.disease, data = dat)
summary(cox_model)


################ AUC FOR BASE MODEL AND BASE MODEL + BIOMARKERS ################

library(survival)
library(timeROC)
library(pROC)

# Fit the Cox proportional hazards model on the complete cases
cox_model_base <- coxph(Surv(timeto1.minor.bleed, minor.bleed) ~ age.bl + pat.sex + bmi + 
                          current.smoker + rr.sys.liegend + 
                          prev.diabetes + prev.stroke.tia + prev.heart.failure + 
                          prev.niereninsuff + coronary.heart.disease, data = data)

cox_model_biomarkers <- coxph(Surv(timeto1.minor.bleed, minor.bleed) ~ gdf.15_log_std + 
                                il6_log_std + probnpii_log_std + 
                                age.bl + pat.sex + bmi + current.smoker + rr.sys.liegend + 
                                prev.diabetes + prev.stroke.tia + prev.heart.failure + 
                                prev.niereninsuff + coronary.heart.disease, data = data)

# Obtain the risk scores (linear predictors) for each model
risk_scores_base <- predict(cox_model_base, type = "lp", newdata = data)
risk_scores_biomarkers <- predict(cox_model_biomarkers, type = "lp", newdata = data)

# Define time points at which to calculate the AUC
time_points <- c(1, 2, 3) * 365.25  # 1, 2, and 3 years in days

# Calculate time-dependent AUC for the base model
roc_results_base <- timeROC(T = data$timeto1.minor.bleed,
                            delta = data$minor.bleed,
                            marker = risk_scores_base,
                            cause = 1,
                            weighting = "marginal",
                            times = time_points,
                            ROC = TRUE)

# Calculate time-dependent AUC for the biomarker model
roc_results_biomarkers <- timeROC(T = data$timeto1.minor.bleed,
                                  delta = data$minor.bleed,
                                  marker = risk_scores_biomarkers,
                                  cause = 1,
                                  weighting = "marginal",
                                  times = time_points,
                                  ROC = TRUE)

# Print AUC values for both models
print(roc_results_base$AUC)
print(roc_results_biomarkers$AUC)

# Extract predicted probabilities at the specific time point (1 year, 365.25 days)
pred_probs_base <- predict(cox_model_base, type = "risk", newdata = data)
pred_probs_biomarkers <- predict(cox_model_biomarkers, type = "risk", newdata = data)

# Use pROC to create ROC curves at the 1-year mark
roc_base <- roc(data$minor.bleed, pred_probs_base, plot = FALSE, ci = TRUE, quiet = TRUE)
roc_biomarkers <- roc(data$minor.bleed, pred_probs_biomarkers, plot = FALSE, ci = TRUE, quiet = TRUE)

# Perform DeLong's test to compare the AUCs
delong_test <- roc.test(roc_base, roc_biomarkers)
print(delong_test)

# Plot ROC Curves at different time points
plot(roc_results, time = time_points[1], col = "blue", title = "ROC Curve at 1 Year")
plot(roc_results, time = time_points[2], col = "red", title = "ROC Curve at 2 Years", add = TRUE)
plot(roc_results, time = time_points[3], col = "green", title = "ROC Curve at 3 Years", add = TRUE)
legend("bottomright", legend = c("1 Year", "2 Years", "3 Years"), col = c("blue", "red", "green"), lwd = 2)

################################ VOLCANO PLOTS #################################

# Create an empty dataframe to store results
volcano_data <- data.frame(Biomarker = character(),
                           HR = numeric(),
                           P_Value = numeric(),
                           stringsAsFactors = FALSE)

# Iterate through each biomarker
for (biomarker in biomarkers) {
  
  # Fit the Cox model
  cox_model <- coxph(Surv(timeto1.minor.bleed, minor.bleed) ~
                       get(paste0(biomarker, "_log_std")) + age.bl + pat.sex + bmi +
                       current.smoker + rr.sys.liegend + prev.diabetes +
                       prev.stroke.tia + prev.heart.failure + prev.niereninsuff +
                       coronary.heart.disease, data = dat)
  
  # Extract HR and p-value for _log_std
  coef_summary <- coef(summary(cox_model))
  hr <- coef_summary["get(paste0(biomarker, \"_log_std\"))", "exp(coef)"]
  p_value <- coef_summary["get(paste0(biomarker, \"_log_std\"))", "Pr(>|z|)"]
  
  # Add results to the dataframe
  volcano_data <- rbind(volcano_data, data.frame(biomarker = biomarker, HR = hr, P_Value = p_value))
}

# Volcano plot from the combined model after stepwise biomarker selection
volcano_data <- data.frame(Biomarker = character(),
                           HR = numeric(),
                           P_Value = numeric(),
                           stringsAsFactors = FALSE)

biomarker_names <- c("gdf.15", "il6", "probnpii")

cox_model <- coxph(Surv(timeto1.minor.bleed, minor.bleed) ~ gdf.15_log_std + 
                     il6_log_std + probnpii_log_std + 
                     age.bl + pat.sex + bmi + current.smoker + rr.sys.liegend + 
                     prev.diabetes + prev.stroke.tia + prev.heart.failure + 
                     prev.niereninsuff + coronary.heart.disease, data = dat)

# Obtain summary of coefficients
coef_summary <- summary(cox_model)

hr <- coef_summary$coefficients[, "exp(coef)"][paste0(biomarker_names, "_log_std")]
p_value <- coef_summary$coefficients[, "Pr(>|z|)"][paste0(biomarker_names, "_log_std")]

# Combine results into a data frame
volcano_data <- rbind(volcano_data, data.frame(biomarker = biomarker_names, HR = hr, P_Value = p_value))

# Display the dataframe
print(volcano_data)

new_names <- c("GDF-15", "IL-6", "NT-proBNP")

# Replace biomarker names in coef_data
volcano_data$biomarker <- ifelse(volcano_data$biomarker == "gdf.15", new_names[1],
                                 ifelse(volcano_data$biomarker == "il6", new_names[2],
                                        ifelse(volcano_data$biomarker == "probnpii", new_names[3],
                                               volcano_data$biomarker)))

library(ggplot2)
library(ggrepel)

ggplot(volcano_data, aes(x = HR, y = -log10(P_Value), label = biomarker)) +
  geom_point() +
  labs(x = "Standardized hazard ratio", y = "-Log10 P-value") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "red") +
  theme_bw() + 
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 9, color = "black"),
        axis.text.x = element_text(hjust = 0.5)) +
  coord_cartesian(xlim = c(0.5, 2.0)) +
  scale_x_continuous(trans = "log10", breaks = c(0.5, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2.0), labels = c(0.5, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2.0)) +
  geom_text_repel(size = 2.75)

########################## VARIABLE IMPORTANCE PLOTS ###########################

library(dplyr)

features <- dat %>% 
  select(ends_with("_log_std"), age.bl, pat.sex, bmi, current.smoker, rr.sys.liegend,
         prev.diabetes, prev.stroke.tia, prev.heart.failure, prev.niereninsuff,
         coronary.heart.disease, timeto1.minor.bleed, minor.bleed)

features <- na.omit(features)

cox_model <- coxph(Surv(timeto1.minor.bleed, minor.bleed) ~ gdf.15_log_std + 
                     il6_log_std + probnpii_log_std + 
                     age.bl + pat.sex + bmi + current.smoker + rr.sys.liegend + 
                     prev.diabetes + prev.stroke.tia + prev.heart.failure + 
                     prev.niereninsuff + coronary.heart.disease, data = features)
summary(cox_model)

# Extract coefficient estimates and standard errors
coef_est <- coef(cox_model)
coef_se <- sqrt(diag(vcov(cox_model)))

# Calculate partial chi-squared statistic for each predictor
partial_chisq <- (coef_est / coef_se)^2

# Extract the total degrees of freedom for the model
num_params <- length(coef(cox_model))
df_total <- num_params

# Adjust partial chi-squared statistic by subtracting total degrees of freedom
partial_chisq_adj <- partial_chisq - df_total

# Display the results
result <- data.frame(Variable = names(coef_est), Partial_ChiSq = partial_chisq)
print(result)

# Order the result by partial chi-squared statistic in descending order
result <- result[order(-result$Partial_ChiSq), ]

library(ggplot2)

# Define the desired label changes
label_changes <- c("gdf.15_log_std" = "GDF-15",
                   "il6_log_std" = "IL-6",
                   "probnpii_log_std" = "NT-proBNP",
                   "age.bl" = "Age",
                   "pat.sex" = "Sex",
                   "bmi" = "BMI",
                   "current.smoker" = "Smoker",
                   "rr.sys.liegend" = "SBP",
                   "prev.diabetes" = "Diabetes",
                   "prev.stroke.tia" = "Prior stroke/TIA",
                   "prev.heart.failure" = "Heart failure",
                   "prev.niereninsuff" = "CKD",
                   "coronary.heart.disease" = "CAD")

# Create the dot plot
p <- ggplot(result, aes(x = Partial_ChiSq, y = reorder(Variable, Partial_ChiSq))) +
  geom_point(shape = 1, size = 2.5, color = "black") +  # Change shape to open circle and set color to black
  labs(x = "Partial χ2-df", size = 2, y = "", color = "black") +  # Fixed x-axis label and removed y-axis title
  scale_y_discrete(labels = label_changes) +  # Apply label changes to y-axis
  theme_bw() +  # Minimal theme
  theme(axis.text.x = element_text(size = 10),  # Adjust x-axis text size
        axis.text.y = element_text(size = 10),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(color = "black"),  # Set x-axis ticks color to black
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.major.y = element_line(linetype = "dotted", color = "grey")  # Specify major horizontal grid lines as dotted and grey
  )

print(p)

############################## FOREST PLOTS ####################################

# Combined model
library(randomForest)
library(survminer)

cox_model <- coxph(Surv(timeto1.minor.bleed, minor.bleed) ~ gdf.15_log_std + 
                     il6_log_std + probnpii_log_std + 
                     age.bl + pat.sex + bmi + current.smoker + rr.sys.liegend + 
                     prev.diabetes + prev.stroke.tia + prev.heart.failure + 
                     prev.niereninsuff + coronary.heart.disease, data = dat)
summary(cox_model)

biomarkers <- c("gdf.15_log_std", "il6_log_std", "probnpii_log_std")

# Create a dataframe to store hazard ratios and 95% CIs
coef_data <- data.frame(
  Biomarker = biomarkers,
  HR = exp(coef(cox_model)[biomarkers]),             # Extracting hazard ratios
  Lower_CI = exp(confint(cox_model)[biomarkers, 1]),  # Extracting lower CIs
  Upper_CI = exp(confint(cox_model)[biomarkers, 2]),   # Extracting upper CIs
  p_value = summary(cox_model)$coefficients[biomarkers, "Pr(>|z|)"]  # Extracting p-values
)

# View the dataframe
print(coef_data)

new_names <- c("GDF-15", "IL-6", "NT-proBNP")

# Replace biomarker names in coef_data
coef_data$Biomarker <- ifelse(coef_data$Biomarker == "gdf.15_log_std", new_names[1],
                              ifelse(coef_data$Biomarker == "il6_log_std", new_names[2],
                                     ifelse(coef_data$Biomarker == "probnpii_log_std", new_names[3],
                                            coef_data$Biomarker)))

coef_data$p.signif <- ifelse(coef_data$p_value < 0.05, "*", "")

library(ggplot2)

# Create coefficient plot using ggplot2
coefficient_plot <- ggplot(coef_data, aes(x = Biomarker, y = HR, color = Biomarker)) +
  geom_point(size = 2.5) + 
  geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 0.3, linetype = "solid") + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") + 
  labs(x = "Biomarker", y = "Hazard ratio (95% CI) per 1 SD increase in biomarker levels", size=11, color = "black") + 
  theme_bw() +
  theme(legend.position="none",
        axis.title = element_text(face="bold"),
        axis.text.x = element_text(size=9, color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(size=9, color = "black")) +
  scale_y_continuous(limits = c(0.8, 1.5), breaks = seq(0.8, 1.5, by = 0.1)) +
  geom_text(aes(label = paste0(sprintf("%.2f", HR))), 
            hjust = -0.3, size = 3)
print(coefficient_plot)

################################# END ##########################################
