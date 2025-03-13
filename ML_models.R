#--ML MODELS--------------------------------------------------------------------
# Author: Pascal B. Meyre
# Date: 05/18/24, latest update: 03/13/25
# Location: Reinach, Baselland, Switzerland

# Use cleaned and imputed dataset of SWISS.BEAT.biomarker.cleaned.imputed.csv
#-------------------------------------------------------------------------------

# Read the cleaned dataset
library(data.table)
dat <- fread("/Users/pascalmeyre/Desktop/Research/1_Projects_Analysis/18_Biomarkers_MACE_bleeding/analysis/datasets/SWISS.BEAT.biomarker.cleaned.imputed.csv")

dat$crep2 <- ifelse(dat$crep2 == 1, NA, dat$crep2)
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

# Relabel variables
# Center
# Define the mapping of center names to numeric values
center_map <- c(
  "Universitätsspital Basel" = 1,
  "Inselspital, Universitätsspital Bern" = 2,
  "Triemli Spital Zürich" = 3,
  "Kantonsspital Baden" = 4,
  "Kantonsspital St. Gallen" = 5,
  "Cardiocentro Lugano" = 6,
  "EOC Lugano" = 7,
  "Luzerner Kantonsspital" = 8,
  "Hôpital Cantonal Fribourg" = 9,
  "HUG" = 10,
  "Bürgerspital Solothurn" = 11,
  "CHUV" = 12,
  "Kantonsspital Luzern" = 8,
  "Regionalspital Lugano" = 7,
  "EOC Bellinzona" = 13,
  "Regionalspital Bellinzona" = 13,
  "HUG Genève" = 11,
  "Hirslanden Klinik St. Anna, Luzern" = 14,
  "CHUV Lausanne" = 12,
  "Universitätsspital Zürich" = 15
)

# Relabel the 'center' variable using the mapping
dat$center <- center_map[dat$center]

# Herkunft
# Define the mapping of 'herkunft' values to numeric values
herkunft_map <- c(
  "Mitteleuropa" = 1,
  "Südeuropa" = 2,
  "Nordeuropa" = 3,
  "Osteuropa" = 4,
  "andere" = 5,
  "Mittel-/Südamerika" = 6,
  "Nordamerika/weiss" = 7
)

# Relabel the 'herkunft' variable using the mapping
dat$herkunft <- herkunft_map[dat$herkunft]

# AF type
dat$vhf.typ.aktuell.bl <- ifelse(dat$vhf.typ.aktuell.bl == "paroxysmal", 1,
                                 ifelse(dat$vhf.typ.aktuell.bl == "persistierend (>7 Tage, EKV)", 2,
                                        ifelse(dat$vhf.typ.aktuell.bl == "permanent", 3, NA)))

# Alcohol beer
# Define the mapping of alcohol consumption categories to numeric values
alkohol_bier_map <- c(
  "6+ pro Tag" = 1,
  "4-5 pro Tag" = 2,
  "2-3 pro Tag" = 3,
  "1 pro Tag" = 4,
  "5-6 pro Woche" = 5,
  "2-4 pro Woche" = 6,
  "1 pro Woche" = 7,
  "1-3 pro Monat" = 8,
  "nie oder weniger als 1 Monat" = 9
)

# Relabel the 'alkohol.bier' variable using the mapping
dat$alkohol.bier <- alkohol_bier_map[dat$alkohol.bier]

# Alcohol red wine
# Define the mapping of alcohol consumption categories to numeric values for 'rotwein'
alkohol_rotwein_map <- c(
  "6+ pro Tag" = 1,
  "4-5 pro Tag" = 2,
  "2-3 pro Tag" = 3,
  "1 pro Tag" = 4,
  "5-6 pro Woche" = 5,
  "2-4 pro Woche" = 6,
  "1 pro Woche" = 7,
  "1-3 pro Monat" = 8,
  "nie oder weniger als 1 Monat" = 9
)

# Relabel the 'alkohol.rotwein' variable using the mapping
dat$alkohol.rotwein <- alkohol_rotwein_map[dat$alkohol.rotwein]

# Alcohol white wine
# Define the mapping of alcohol consumption categories to numeric values for 'weisswein'
alkohol_weisswein_map <- c(
  "6+ pro Tag" = 1,
  "4-5 pro Tag" = 2,
  "2-3 pro Tag" = 3,
  "1 pro Tag" = 4,
  "5-6 pro Woche" = 5,
  "2-4 pro Woche" = 6,
  "1 pro Woche" = 7,
  "1-3 pro Monat" = 8,
  "nie oder weniger als 1 Monat" = 9
)

# Relabel the 'alkohol.weisswein' variable using the mapping
dat$alkohol.weisswein <- alkohol_weisswein_map[dat$alkohol.weisswein]

# Alcohol liquor
# Define the mapping of alcohol consumption categories to numeric values for 'schnaps'
alkohol_schnaps_map <- c(
  "6+ pro Tag" = 1,
  "2-3 pro Tag" = 2,
  "1 pro Tag" = 3,
  "5-6 pro Woche" = 4,
  "2-4 pro Woche" = 5,
  "1 pro Woche" = 6,
  "1-3 pro Monat" = 7,
  "nie oder weniger als 1 Monat" = 8
)

# Relabel the 'alkohol.schnaps' variable using the mapping
dat$alkohol.schnaps <- alkohol_schnaps_map[dat$alkohol.schnaps]

# AF duration
# Define the mapping of 'vhf.episoden.dauer' categories to numeric values
vhf_episoden_dauer_map <- c(
  "dauerhaft" = 1,
  "Tage" = 2,
  "Stunden" = 3,
  "Minuten" = 4,
  "keine mehr" = 5
)

# Relabel the 'vhf.episoden.dauer' variable using the mapping
dat$vhf.episoden.dauer <- vhf_episoden_dauer_map[dat$vhf.episoden.dauer]

# AF nr. of episodes
# Define the mapping of 'vhf.episoden.anz' categories to numeric values
vhf_episoden_anz_map <- c(
  "dauerhaft" = 1,
  ">=1x pro Woche" = 2,
  "=1x pro Woche" = 3,
  "<1x pro Woche aber >1x pro Monat" = 4,
  "<1x pro Monat" = 5,
  "keine mehr" = 6
)

# Relabel the 'vhf.episoden.anz' variable using the mapping
dat$vhf.episoden.anz <- vhf_episoden_anz_map[dat$vhf.episoden.anz]

# AFlutter
dat$vorhofflattern <- ifelse(dat$vorhofflattern == "Yes", 1,
                             ifelse(dat$vorhofflattern == "No", 0, NA))

# Aspirin
dat$med.aspirin <- ifelse(dat$med.aspirin == "Yes", 1,
                          ifelse(dat$med.aspirin == "No", 0, NA))

# TCA
dat$med.tca.yn <- ifelse(dat$med.tca.yn == "Yes", 1,
                         ifelse(dat$med.tca.yn == "No", 0, NA))

# Antiplatelets
dat$med.antiplatelet.yn <- ifelse(dat$med.antiplatelet.yn == "Yes", 1,
                                  ifelse(dat$med.antiplatelet.yn == "No", 0, NA))

# Ilicit drugs
dat$drogen <- ifelse(dat$drogen == "Yes", 1,
                     ifelse(dat$drogen == "No", 0, NA))

# Cancer
dat$krk.malignom <- ifelse(dat$krk.malignom == "Yes", 1,
                           ifelse(dat$krk.malignom == "No", 0, NA))

# PTE
dat$krk.tvt <- ifelse(dat$krk.tvt == "Yes", 1,
                      ifelse(dat$krk.tvt == "No", 0, NA))

# Paternal AF history
dat$fam.vhf.vater <- ifelse(dat$fam.vhf.vater == "Ja", 1,
                            ifelse(dat$fam.vhf.vater == "Nein", 0,
                                   ifelse(dat$fam.vhf.vater == "unbekannt", 0, NA)))

# Brother AF history
dat$fam.vhf.bruder <- ifelse(dat$fam.vhf.bruder == "Ja", 1,
                             ifelse(dat$fam.vhf.bruder == "Nein", 0,
                                    ifelse(dat$fam.vhf.bruder == "unbekannt", 0,
                                           ifelse(dat$fam.vhf.bruder == "kein Bruder", 2, NA))))

# Maternal AF history
dat$fam.vhf.mutter <- ifelse(dat$fam.vhf.mutter == "Ja", 1,
                             ifelse(dat$fam.vhf.mutter == "Nein", 0,
                                    ifelse(dat$fam.vhf.mutter == "unbekannt", 0, NA)))

# Sister AF history
dat$fam.vhf.schwester <- ifelse(dat$fam.vhf.schwester == "Ja", 1,
                                ifelse(dat$fam.vhf.schwester == "Nein", 0,
                                       ifelse(dat$fam.vhf.schwester == "unbekannt", 0,
                                              ifelse(dat$fam.vhf.schwester == "kein Bruder", 2, NA))))

# Family history hypertension
dat$fam.hypertonie <- ifelse(dat$fam.hypertonie == "Ja", 1,
                             ifelse(dat$fam.hypertonie == "Nein", 0,
                                    ifelse(dat$fam.hypertonie == "unbekannt", 0, NA)))

# Family history diabetes
dat$fam.diabetes <- ifelse(dat$fam.diabetes == "Ja", 1,
                           ifelse(dat$fam.diabetes == "Nein", 0,
                                  ifelse(dat$fam.diabetes == "unbekannt", 0, NA)))

# Family history obesity
dat$fam.uebergewicht <- ifelse(dat$fam.uebergewicht == "Ja", 1,
                               ifelse(dat$fam.uebergewicht == "Nein", 0,
                                      ifelse(dat$fam.uebergewicht == "unbekannt", 0, NA)))

# Family history CAD
dat$fam.khk <- ifelse(dat$fam.khk == "Ja", 1,
                      ifelse(dat$fam.khk == "Nein", 0,
                             ifelse(dat$fam.khk == "unbekannt", 0, NA)))

# Rhythm on ECG
# Define the mapping of 'ecg.rhythm.algo' categories to numeric values
ecg_rhythm_map <- c(
  "AF" = 1,
  "Fibrillation" = 1,
  "Flutter" = 2,
  "Sinus" = 3,
  "Other" = 4,
  "unclear" = 4
)

# Relabel the 'ecg.rhythm.algo' variable using the mapping
dat$ecg.rhythm.algo <- ecg_rhythm_map[dat$ecg.rhythm.algo]

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

# Create a single time and outcome variable for the composite outcome
# Calculate the shortest time variable for composite outcome
dat$timeto1.composite <- pmin(dat$timeto1.death, dat$timeto1.isch.stroke, dat$timeto1.sys.embolism, dat$timeto1.mi)

# Create a composite variable indicating whether an event occurred (1) or not (0)
dat$composite <- ifelse((dat$timeto1.composite == dat$timeto1.death) & dat$death.cardiac, 1,
                        ifelse((dat$timeto1.composite == dat$timeto1.isch.stroke) & dat$ischemic.stroke, 1, 
                               ifelse((dat$timeto1.composite == dat$timeto1.mi) & dat$mi, 1, 
                                      ifelse((dat$timeto1.composite == dat$timeto1.sys.embolism) & dat$sys.embolism, 1, 0))))

# Split data into features (biomarkers) and target variable
library(caret)
library(randomForest)
library(glmnet)
library(dplyr)
library(pROC)
library(xgboost)

# Exclude patients with missing variables
data <- dat %>% 
  select(ends_with("_log"), age.bl, pat.sex, bmi, current.smoker, rr.sys.liegend,
         prev.diabetes, prev.stroke.tia, prev.heart.failure, prev.niereninsuff,
         coronary.heart.disease, center, herkunft, vhf.typ.aktuell.bl, alkohol.bier, 
         alkohol.rotwein, alkohol.weisswein, alkohol.schnaps, vhf.episoden.dauer,
         vhf.episoden.anz, vorhofflattern, med.aspirin, med.tca.yn, med.antiplatelet.yn,
         drogen, prev.akb, prev.hyperthyreose, prev.hypothyreose, krk.malignom, krk.tvt,
         fam.vhf.vater, fam.vhf.bruder, fam.vhf.mutter, fam.vhf.schwester, fam.hypertonie,
         fam.diabetes, fam.uebergewicht, fam.khk, groesse, gewicht, heart.rate,
         rr.dia.liegend, ecg.rhythm.algo, prev.schlaf.apnoe, prev.hypertonie, prev.pavk,
         prev.mi, prev.sys.embolism, prev.major.bleed, composite)

data <- na.omit(data)

# Make feature and target outcome data frames"
features <- data %>% 
  select(age.bl, pat.sex, bmi, current.smoker, rr.sys.liegend,
         prev.diabetes, prev.stroke.tia, prev.heart.failure, prev.niereninsuff,
         coronary.heart.disease, center, herkunft, vhf.typ.aktuell.bl, alkohol.bier, 
         alkohol.rotwein, alkohol.weisswein, alkohol.schnaps, vhf.episoden.dauer,
         vhf.episoden.anz, vorhofflattern, med.aspirin, med.tca.yn, med.antiplatelet.yn,
         drogen, prev.akb, prev.hyperthyreose, prev.hypothyreose, krk.malignom, krk.tvt,
         fam.vhf.vater, fam.vhf.bruder, fam.vhf.mutter, fam.vhf.schwester, fam.hypertonie,
         fam.diabetes, fam.uebergewicht, fam.khk, groesse, gewicht, heart.rate,
         rr.dia.liegend, ecg.rhythm.algo, prev.schlaf.apnoe, prev.hypertonie, prev.pavk,
         prev.mi, prev.sys.embolism, prev.major.bleed, composite)

features$composite <- as.factor(features$composite)

#### Base model ################################################################
#### Random forest ####---------------------------------------------------------
# Train random forest model
set.seed(123)

# Train a Random Forest model
rf_model <- randomForest(composite ~ ., data = features, importance = TRUE)

# Get variable importance
importance <- importance(rf_model)
var_importance <- data.frame(Variables = rownames(importance), Importance = importance[, 'MeanDecreaseGini'])

# Get top biomarkers
top_biomarkers <- var_importance %>% arrange(desc(Importance)) %>% head(12)
print(top_biomarkers)

# Get predicted probabilities
pred_probs <- predict(rf_model, type = "prob")[, "1"]

# Get true outcomes
true_outcomes <- ifelse(features$composite == 1, TRUE, FALSE)

# Calculate AUC ROC
roc_curve <- roc(true_outcomes, pred_probs)
auc_base <- roc_curve$auc
print(paste("AUC ROC for Random Forest model:", auc_base))

# Calculate 95% CI for the AUC
ci_auc <- ci.auc(roc_curve)
print(paste("95% CI for AUC:", ci_auc[1], "-", ci_auc[3]))

#### LASSO ####-----------------------------------------------------------------
# Prepare data for LASSO
x <- model.matrix(composite ~ ., features)[, -1]
y <- features$composite

# Train LASSO model
lasso_model <- cv.glmnet(x, y, family = "binomial", alpha = 1)

# Get coefficients from the best model
lasso_coef <- coef(lasso_model, s = "lambda.min")
selected_biomarkers <- rownames(lasso_coef)[lasso_coef[,1] != 0]
print(selected_biomarkers)

# Predict probabilities using the selected features
pred_probs <- predict(lasso_model, newx = x, s = "lambda.min", type = "response")

# Calculate AUC ROC
roc_curve <- roc(y, as.vector(pred_probs))
auc_base <- roc_curve$auc
print(paste("AUC ROC for LASSO model:", auc_base))

# Calculate 95% CI for the AUC
ci_auc <- ci.auc(roc_curve)
print(paste("95% CI for AUC:", ci_auc[1], "-", ci_auc[3]))

#### XGBoost ####---------------------------------------------------------------
new_colnames <- c("Age", "Sex", "BMI", "Smoker", "Systolic blood pressure", "Diabetes", 
                  "Prior stroke/TIA", "Heart failure", "Renal failure", "CAD", "Study center", 
                  "Ethnicity", "AF type", "Beer drinker", "Red wine drinker", "White wine drinker", 
                  "Liquor drinker", "AF episode duration", "Nr. of AF episodes", "Atrial flutter", 
                  "Aspirin", "TCA", "Antiplatelet therapy", "Illicit drugs", "Prior CABG", 
                  "Hyperthyroidism", "Hypothyroidism", "Cancer", "VTE", "Paternal AF", 
                  "Brother AF", "Maternal AF", "Sister AF", "Familial hypertension", "Familial diabetes", 
                  "Familial obesity", "Familial CAD", "Height", "Weight", "Heart rate", 
                  "Diastolic blood pressure", "Rhythm at baseline", "Sleep apnea", "Hypertension", 
                  "PAD", "Prior MI", "Systemic embolism", "Prior major bleeding")

# Assign the new column names to the features data frame
colnames(features) <- new_colnames

features$composite <- as.numeric(features$composite == 1)

# Train XGBoost model
set.seed(123)
xgb_model <- xgboost(data = as.matrix(features[, -c("composite")]), 
                     label = features$composite, 
                     objective = "binary:logistic", 
                     nrounds = 10)

# Get variable importance
var_importance_xgb <- xgb.importance(feature_names = colnames(features[, -c("composite")]), model = xgb_model)

# Get top biomarkers
top_biomarkers_xgb <- var_importance_xgb$Feature[1:12]
print(top_biomarkers_xgb)

# Get predicted probabilities
pred_probs_xgb <- as.numeric(predict(xgb_model, as.matrix(features[, -c("composite")])))

# Get true outcomes
true_outcomes <- ifelse(features$composite == 1, TRUE, FALSE)

# Calculate AUC ROC for XGBoost
roc_curve_xgb <- roc(true_outcomes, pred_probs_xgb)
auc_xgb_base <- roc_curve_xgb$auc
print(paste("AUC ROC for XGBoost model:", auc_xgb_base))

# Calculate 95% CI for the AUC
ci_auc_xgb <- ci.auc(roc_curve_xgb)
print(paste("95% CI for AUC:", ci_auc_xgb[1], "-", ci_auc_xgb[3]))

# Plot variable importance
plot(var_importance_xgb, col = "black", pch = 20, cex = 1.25, top = 20)
xgb.plot.importance(var_importance_xgb, top_n = 20, measure = "Gain", rel_to_first = TRUE, xlab = "Relative Importance")

ggplot(var_importance_xgb, aes(x = reorder(Feature, Gain), y = Gain)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Variable Importance from XGBoost Model",
       x = "Features",
       y = "Gain") +
  theme_minimal()

ggplot(var_importance_xgb, aes(x = Gain, y = reorder(Feature, Gain))) +
  geom_point(shape = 1, size = 2.5, color = "black") +
  labs(x = "Importance", size = 2, y = "", color = "black") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(color = "black"),  # Set x-axis ticks color to black
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.major.y = element_line(linetype = "dotted", color = "grey"))

#### Model with biomarkers #####################################################
#### Random forest ####---------------------------------------------------------
library(caret)
library(randomForest)
library(glmnet)
library(dplyr)
library(pROC)
library(xgboost)

# Exclude patients with missing variables
data <- dat %>% 
  select(ends_with("_log"), age.bl, pat.sex, bmi, current.smoker, rr.sys.liegend,
         prev.diabetes, prev.stroke.tia, prev.heart.failure, prev.niereninsuff,
         coronary.heart.disease, center, herkunft, vhf.typ.aktuell.bl, alkohol.bier, 
         alkohol.rotwein, alkohol.weisswein, alkohol.schnaps, vhf.episoden.dauer,
         vhf.episoden.anz, vorhofflattern, med.aspirin, med.tca.yn, med.antiplatelet.yn,
         drogen, prev.akb, prev.hyperthyreose, prev.hypothyreose, krk.malignom, krk.tvt,
         fam.vhf.vater, fam.vhf.bruder, fam.vhf.mutter, fam.vhf.schwester, fam.hypertonie,
         fam.diabetes, fam.uebergewicht, fam.khk, groesse, gewicht, heart.rate,
         rr.dia.liegend, ecg.rhythm.algo, prev.schlaf.apnoe, prev.hypertonie, prev.pavk,
         prev.mi, prev.sys.embolism, prev.major.bleed, composite)

data <- na.omit(data)

# Make feature and target outcome data frames"
features <- data %>% 
  select(age.bl, pat.sex, bmi, current.smoker, rr.sys.liegend,
         prev.diabetes, prev.stroke.tia, prev.heart.failure, prev.niereninsuff,
         coronary.heart.disease, center, herkunft, vhf.typ.aktuell.bl, alkohol.bier, 
         alkohol.rotwein, alkohol.weisswein, alkohol.schnaps, vhf.episoden.dauer,
         vhf.episoden.anz, vorhofflattern, med.aspirin, med.tca.yn, med.antiplatelet.yn,
         drogen, prev.akb, prev.hyperthyreose, prev.hypothyreose, krk.malignom, krk.tvt,
         fam.vhf.vater, fam.vhf.bruder, fam.vhf.mutter, fam.vhf.schwester, fam.hypertonie,
         fam.diabetes, fam.uebergewicht, fam.khk, groesse, gewicht, heart.rate,
         rr.dia.liegend, ecg.rhythm.algo, prev.schlaf.apnoe, prev.hypertonie, prev.pavk,
         prev.mi, prev.sys.embolism, prev.major.bleed, ends_with("_log"), composite)

features$composite <- as.factor(features$composite)

new_colnames <- c("Age", "Sex", "BMI", "Smoker", "SBP", "Diabetes", "Prior.stroke.TIA", 
                  "Heart.failure", "Renal.failure", "CAD", "Study.center", "Ethnicity", 
                  "AF.type", "Beer.drinker", "Red.wine.drinker", "White.wine.drinker", 
                  "Liquor.drinker", "AF.episode.duration", "Nr.of.AF.episodes", "Atrial.flutter", 
                  "Aspirin", "TCA", "Antiplatelet.therapy", "Illicit.drugs", "Prior.CABG", 
                  "Hyperthyroidism", "Hypothyroidism", "Cancer", "VTE", "Paternal.AF", 
                  "Brother.AF", "Maternal.AF", "Sister.AF", "Familial.hypertension", "Familial.diabetes", 
                  "Familial.obesity", "Familial.CAD", "Height", "Weight", "Heart.rate", "DBP", 
                  "Rhythm.at.baseline", "Sleep.apnea", "Hypertension", "PAD", "Prior.MI", 
                  "Systemic.embolism", "Prior.major.bleeding", "ANG.2", "d.dimer", "cystatin.C", 
                  "ALAT", "GDF.15", "hs.CRP", "IGFBP.7", "IL.6", "NT.proBNP", "OPN", "hsTropT", "eGFR")

# Assign the new column names to the features data frame
colnames(features) <- new_colnames

# Train random forest model
set.seed(123)

# Train a Random Forest model
rf_model <- randomForest(composite ~ ., data = features, importance = TRUE)

# Get variable importance
importance <- importance(rf_model)
var_importance <- data.frame(Variables = rownames(importance), Importance = importance[, 'MeanDecreaseGini'])

# Get top biomarkers
top_biomarkers <- var_importance %>% arrange(desc(Importance)) %>% head(12)
print(top_biomarkers)

# Get predicted probabilities
pred_probs <- predict(rf_model, type = "prob")[, "1"]

# Get true outcomes
true_outcomes <- ifelse(features$composite == 1, TRUE, FALSE)

# Calculate AUC ROC
roc_curve <- roc(true_outcomes, pred_probs)
auc_biomarkers <- roc_curve$auc
print(paste("AUC ROC for Random Forest model:", auc_biomarkers))

# Calculate 95% CI for the AUC
ci_auc <- ci.auc(roc_curve)
print(paste("95% CI for AUC:", ci_auc[1], "-", ci_auc[3]))

# Perform DeLong's test to compare the AUCs
delong_test <- roc.test(auc_base, auc_biomarkers)
print(delong_test)

# Plot variable importance of the random forest model
ggplot(var_importance, aes(x = Importance, y = reorder(Variables, Importance))) +
  geom_point(shape = 1, size = 2.5, color = "black") +
  labs(x = "Importance", y = "Variables") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(color = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.major.y = element_line(linetype = "dotted", color = "grey"))

#### LASSO ####
# Prepare data for LASSO
x <- model.matrix(composite ~ ., features)[, -1]
y <- features$composite

# Train LASSO model
lasso_model <- cv.glmnet(x, y, family = "binomial", alpha = 1)

# Get coefficients from the best model
lasso_coef <- coef(lasso_model, s = "lambda.min")
selected_biomarkers <- rownames(lasso_coef)[lasso_coef[,1] != 0]
print(selected_biomarkers)

# Get coefficients from the best model
lasso_coef <- coef(lasso_model, s = "lambda.min")
lasso_coef <- as.matrix(lasso_coef)

# Create a data frame of coefficients (excluding the intercept)
var_importance <- data.frame(
  Variables = rownames(lasso_coef)[-1],  # Exclude the intercept
  Importance = abs(lasso_coef[-1, 1])  # Take the absolute value of the coefficients
)

# Predict probabilities using the selected features
pred_probs <- predict(lasso_model, newx = x, s = "lambda.min", type = "response")

# Calculate AUC ROC
roc_curve <- roc(true_outcomes, pred_probs)
auc_biomarkers <- roc_curve$auc
print(paste("AUC ROC for Random Forest model:", auc_biomarkers))

# Calculate 95% CI for the AUC
ci_auc <- ci.auc(roc_curve)
print(paste("95% CI for AUC:", ci_auc[1], "-", ci_auc[3]))

# Perform DeLong's test to compare the AUCs
delong_test <- roc.test(auc_base, auc_biomarkers)
print(delong_test)

# Plot variable importance of the LASSO model
ggplot(var_importance, aes(x = Importance, y = reorder(Variables, Importance))) +
  geom_point(shape = 1, size = 2.5, color = "black") +
  labs(x = "Importance", y = "Variables") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(color = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.major.y = element_line(linetype = "dotted", color = "grey"))

#### XGBoost ####---------------------------------------------------------------
new_colnames <- c("Age", "Sex", "BMI", "Smoker", "SBP", "Diabetes", "Prior.stroke.TIA", 
                  "Heart.failure", "Renal.failure", "CAD", "Study.center", "Ethnicity", 
                  "AF.type", "Beer.drinker", "Red.wine.drinker", "White.wine.drinker", 
                  "Liquor.drinker", "AF.episode.duration", "Nr.of.AF.episodes", "Atrial.flutter", 
                  "Aspirin", "TCA", "Antiplatelet.therapy", "Illicit.drugs", "Prior.CABG", 
                  "Hyperthyroidism", "Hypothyroidism", "Cancer", "VTE", "Paternal.AF", 
                  "Brother.AF", "Maternal.AF", "Sister.AF", "Familial.hypertension", "Familial.diabetes", 
                  "Familial.obesity", "Familial.CAD", "Height", "Weight", "Heart.rate", "DBP", 
                  "Rhythm.at.baseline", "Sleep.apnea", "Hypertension", "PAD", "Prior.MI", 
                  "Systemic.embolism", "Prior.major.bleeding", "ANG.2", "d.dimer", "cystatin.C", 
                  "ALAT", "GDF.15", "hs.CRP", "IGFBP.7", "IL.6", "NT.proBNP", "OPN", "hsTropT", "eGFR")

# Assign the new column names to the features data frame
colnames(features) <- new_colnames

features$composite <- as.numeric(features$composite == 1)

# Train XGBoost model
set.seed(123)
xgb_model <- xgboost(data = as.matrix(features[, -c("composite")]), 
                     label = features$composite, 
                     objective = "binary:logistic", 
                     nrounds = 10)

# Get variable importance
var_importance_xgb <- xgb.importance(feature_names = colnames(features[, -c("composite")]), model = xgb_model)

# Get top biomarkers
top_biomarkers_xgb <- var_importance_xgb$Feature[1:12]
print(top_biomarkers_xgb)

# Get predicted probabilities
pred_probs_xgb <- as.numeric(predict(xgb_model, as.matrix(features[, -c("composite")])))

# Get true outcomes
true_outcomes <- ifelse(features$composite == 1, TRUE, FALSE)

# Calculate AUC ROC for XGBoost
roc_curve_xgb <- roc(true_outcomes, pred_probs_xgb)
auc_xgb_biomarkers <- roc_curve_xgb$auc
print(paste("AUC ROC for XGBoost model:", auc_xgb_biomarkers))

# Calculate 95% CI for the AUC
ci_auc_xgb <- ci.auc(roc_curve_xgb)
print(paste("95% CI for AUC:", ci_auc_xgb[1], "-", ci_auc_xgb[3]))

# Perform DeLong's test to compare the AUCs
delong_test <- roc.test(auc_xgb_base, auc_xgb_biomarkers)
print(delong_test)

# Plot variable importance
plot(var_importance_xgb, col = "black", pch = 20, cex = 1.25, top = 20)

ggplot(var_importance_xgb, aes(x = Gain, y = reorder(Feature, Gain))) +
  geom_point(shape = 1, size = 2.5, color = "black") +
  labs(x = "Importance", size = 2, y = "Variables", color = "black") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(color = "black"),  # Set x-axis ticks color to black
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.major.y = element_line(linetype = "dotted", color = "grey"))

################################################################################
# Heart failure hospitalization
################################################################################

#-------------------------------------------------------------------------------
# Check whether there are missing "time" variables 
num_missing_timeto1.hf <- sum(is.na(dat$timeto1.hf))
print(num_missing_timeto1.hf)
#-------------------------------------------------------------------------------

# Split data into features (biomarkers) and target variable
library(caret)
library(randomForest)
library(glmnet)
library(dplyr)
library(pROC)
library(xgboost)

# Exclude patients with missing variables
data <- dat %>% 
  select(ends_with("_log"), age.bl, pat.sex, bmi, current.smoker, rr.sys.liegend,
         prev.diabetes, prev.stroke.tia, prev.heart.failure, prev.niereninsuff,
         coronary.heart.disease, center, herkunft, vhf.typ.aktuell.bl, alkohol.bier, 
         alkohol.rotwein, alkohol.weisswein, alkohol.schnaps, vhf.episoden.dauer,
         vhf.episoden.anz, vorhofflattern, med.aspirin, med.tca.yn, med.antiplatelet.yn,
         drogen, prev.akb, prev.hyperthyreose, prev.hypothyreose, krk.malignom, krk.tvt,
         fam.vhf.vater, fam.vhf.bruder, fam.vhf.mutter, fam.vhf.schwester, fam.hypertonie,
         fam.diabetes, fam.uebergewicht, fam.khk, groesse, gewicht, heart.rate,
         rr.dia.liegend, ecg.rhythm.algo, prev.schlaf.apnoe, prev.hypertonie, prev.pavk,
         prev.mi, prev.sys.embolism, prev.major.bleed, heart.failure)

data <- na.omit(data)

# Make feature and target outcome data frames"
features <- data %>% 
  select(age.bl, pat.sex, bmi, current.smoker, rr.sys.liegend,
         prev.diabetes, prev.stroke.tia, prev.heart.failure, prev.niereninsuff,
         coronary.heart.disease, center, herkunft, vhf.typ.aktuell.bl, alkohol.bier, 
         alkohol.rotwein, alkohol.weisswein, alkohol.schnaps, vhf.episoden.dauer,
         vhf.episoden.anz, vorhofflattern, med.aspirin, med.tca.yn, med.antiplatelet.yn,
         drogen, prev.akb, prev.hyperthyreose, prev.hypothyreose, krk.malignom, krk.tvt,
         fam.vhf.vater, fam.vhf.bruder, fam.vhf.mutter, fam.vhf.schwester, fam.hypertonie,
         fam.diabetes, fam.uebergewicht, fam.khk, groesse, gewicht, heart.rate,
         rr.dia.liegend, ecg.rhythm.algo, prev.schlaf.apnoe, prev.hypertonie, prev.pavk,
         prev.mi, prev.sys.embolism, prev.major.bleed, heart.failure)

features$heart.failure <- as.factor(features$heart.failure)

#### Base model ################################################################
#### Random forest ####---------------------------------------------------------
# Train random forest model
set.seed(123)

# Train a Random Forest model
rf_model <- randomForest(heart.failure ~ ., data = features, importance = TRUE)

# Get variable importance
importance <- importance(rf_model)
var_importance <- data.frame(Variables = rownames(importance), Importance = importance[, 'MeanDecreaseGini'])

# Get top biomarkers
top_biomarkers <- var_importance %>% arrange(desc(Importance)) %>% head(12)
print(top_biomarkers)

# Get predicted probabilities
pred_probs <- predict(rf_model, type = "prob")[, "1"]

# Get true outcomes
true_outcomes <- ifelse(features$heart.failure == 1, TRUE, FALSE)

# Calculate AUC ROC
roc_curve <- roc(true_outcomes, pred_probs)
auc_base <- roc_curve$auc
print(paste("AUC ROC for Random Forest model:", auc_base))

# Calculate 95% CI for the AUC
ci_auc <- ci.auc(roc_curve)
print(paste("95% CI for AUC:", ci_auc[1], "-", ci_auc[3]))

#### LASSO ####-----------------------------------------------------------------
# Prepare data for LASSO
x <- model.matrix(heart.failure ~ ., features)[, -1]
y <- features$heart.failure

# Train LASSO model
lasso_model <- cv.glmnet(x, y, family = "binomial", alpha = 1)

# Get coefficients from the best model
lasso_coef <- coef(lasso_model, s = "lambda.min")
selected_biomarkers <- rownames(lasso_coef)[lasso_coef[,1] != 0]
print(selected_biomarkers)

# Predict probabilities using the selected features
pred_probs <- predict(lasso_model, newx = x, s = "lambda.min", type = "response")

# Calculate AUC ROC
roc_curve <- roc(y, as.vector(pred_probs))
auc_base <- roc_curve$auc
print(paste("AUC ROC for Random Forest model:", auc_base))

# Calculate 95% CI for the AUC
ci_auc <- ci.auc(roc_curve)
print(paste("95% CI for AUC:", ci_auc[1], "-", ci_auc[3]))

#### XGBoost ####---------------------------------------------------------------
new_colnames <- c("Age", "Sex", "BMI", "Smoker", "Systolic blood pressure", "Diabetes", 
                  "Prior stroke/TIA", "Heart failure", "Renal failure", "CAD", "Study center", 
                  "Ethnicity", "AF type", "Beer drinker", "Red wine drinker", "White wine drinker", 
                  "Liquor drinker", "AF episode duration", "Nr. of AF episodes", "Atrial flutter", 
                  "Aspirin", "TCA", "Antiplatelet therapy", "Illicit drugs", "Prior CABG", 
                  "Hyperthyroidism", "Hypothyroidism", "Cancer", "VTE", "Paternal AF", "Brother AF", 
                  "Maternal AF", "Sister AF", "Familial hypertension", "Familial diabetes", 
                  "Familial obesity", "Familial CAD", "Height", "Weight", "Heart rate", 
                  "Diastolic blood pressure", "Rhythm at baseline", "Sleep apnea", "Hypertension", 
                  "PAD", "Prior MI", "Systemic embolism", "Prior major bleeding")

# Assign the new column names to the features data frame
colnames(features) <- new_colnames

features$heart.failure <- as.numeric(features$heart.failure == 1)

# Train XGBoost model
set.seed(123)
xgb_model <- xgboost(data = as.matrix(features[, -c("heart.failure")]), 
                     label = features$heart.failure, 
                     objective = "binary:logistic", 
                     nrounds = 10)

# Get variable importance
var_importance_xgb <- xgb.importance(feature_names = colnames(features[, -c("heart.failure")]), model = xgb_model)

# Get top biomarkers
top_biomarkers_xgb <- var_importance_xgb$Feature[1:12]
print(top_biomarkers_xgb)

# Get predicted probabilities
pred_probs_xgb <- as.numeric(predict(xgb_model, as.matrix(features[, -c("heart.failure")])))

# Get true outcomes
true_outcomes <- ifelse(features$heart.failure == 1, TRUE, FALSE)

# Calculate AUC ROC for XGBoost
roc_curve_xgb <- roc(true_outcomes, pred_probs_xgb)
auc_xgb_base <- roc_curve_xgb$auc
print(paste("AUC ROC for XGBoost model:", auc_xgb_base))

# Calculate 95% CI for the AUC
ci_auc_xgb <- ci.auc(roc_curve_xgb)
print(paste("95% CI for AUC:", ci_auc_xgb[1], "-", ci_auc_xgb[3]))

# Plot variable importance
plot(var_importance_xgb, col = "black", pch = 20, cex = 1.25, top = 20)
xgb.plot.importance(var_importance_xgb, top_n = 20, measure = "Gain", rel_to_first = TRUE, xlab = "Relative Importance")

ggplot(var_importance_xgb, aes(x = reorder(Feature, Gain), y = Gain)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Variable Importance from XGBoost Model",
       x = "Features",
       y = "Gain") +
  theme_minimal()

ggplot(var_importance_xgb, aes(x = Gain, y = reorder(Feature, Gain))) +
  geom_point(shape = 1, size = 2.5, color = "black") +
  labs(x = "Importance", size = 2, y = "", color = "black") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(color = "black"),  # Set x-axis ticks color to black
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.major.y = element_line(linetype = "dotted", color = "grey"))

#### Model with biomarkers #####################################################
#### Random forest ####---------------------------------------------------------
library(caret)
library(randomForest)
library(glmnet)
library(dplyr)
library(pROC)
library(xgboost)

# Exclude patients with missing variables
data <- dat %>% 
  select(ends_with("_log"), age.bl, pat.sex, bmi, current.smoker, rr.sys.liegend,
         prev.diabetes, prev.stroke.tia, prev.heart.failure, prev.niereninsuff,
         coronary.heart.disease, center, herkunft, vhf.typ.aktuell.bl, alkohol.bier, 
         alkohol.rotwein, alkohol.weisswein, alkohol.schnaps, vhf.episoden.dauer,
         vhf.episoden.anz, vorhofflattern, med.aspirin, med.tca.yn, med.antiplatelet.yn,
         drogen, prev.akb, prev.hyperthyreose, prev.hypothyreose, krk.malignom, krk.tvt,
         fam.vhf.vater, fam.vhf.bruder, fam.vhf.mutter, fam.vhf.schwester, fam.hypertonie,
         fam.diabetes, fam.uebergewicht, fam.khk, groesse, gewicht, heart.rate,
         rr.dia.liegend, ecg.rhythm.algo, prev.schlaf.apnoe, prev.hypertonie, prev.pavk,
         prev.mi, prev.sys.embolism, prev.major.bleed, heart.failure)

data <- na.omit(data)

# Make feature and target outcome data frames"
features <- data %>% 
  select(age.bl, pat.sex, bmi, current.smoker, rr.sys.liegend,
         prev.diabetes, prev.stroke.tia, prev.heart.failure, prev.niereninsuff,
         coronary.heart.disease, center, herkunft, vhf.typ.aktuell.bl, alkohol.bier, 
         alkohol.rotwein, alkohol.weisswein, alkohol.schnaps, vhf.episoden.dauer,
         vhf.episoden.anz, vorhofflattern, med.aspirin, med.tca.yn, med.antiplatelet.yn,
         drogen, prev.akb, prev.hyperthyreose, prev.hypothyreose, krk.malignom, krk.tvt,
         fam.vhf.vater, fam.vhf.bruder, fam.vhf.mutter, fam.vhf.schwester, fam.hypertonie,
         fam.diabetes, fam.uebergewicht, fam.khk, groesse, gewicht, heart.rate,
         rr.dia.liegend, ecg.rhythm.algo, prev.schlaf.apnoe, prev.hypertonie, prev.pavk,
         prev.mi, prev.sys.embolism, prev.major.bleed, ends_with("_log"), heart.failure)

features$heart.failure <- as.factor(features$heart.failure)

# Train random forest model
set.seed(123)

# Train a Random Forest model
rf_model <- randomForest(heart.failure ~ ., data = features, importance = TRUE)

# Get variable importance
importance <- importance(rf_model)
var_importance <- data.frame(Variables = rownames(importance), Importance = importance[, 'MeanDecreaseGini'])

# Get top biomarkers
top_biomarkers <- var_importance %>% arrange(desc(Importance)) %>% head(12)
print(top_biomarkers)

# Get predicted probabilities
pred_probs <- predict(rf_model, type = "prob")[, "1"]

# Get true outcomes
true_outcomes <- ifelse(features$heart.failure == 1, TRUE, FALSE)

# Calculate AUC ROC
roc_curve <- roc(true_outcomes, pred_probs)
auc_biomarkers <- roc_curve$auc
print(paste("AUC ROC for Random Forest model:", auc_biomarkers))

# Calculate 95% CI for the AUC
ci_auc <- ci.auc(roc_curve)
print(paste("95% CI for AUC:", ci_auc[1], "-", ci_auc[3]))

# Perform DeLong's test to compare the AUCs
delong_test <- roc.test(auc_base, auc_biomarkers)
print(delong_test)

#### LASSO ####
# Prepare data for LASSO
x <- model.matrix(heart.failure ~ ., features)[, -1]
y <- features$heart.failure

# Train LASSO model
lasso_model <- cv.glmnet(x, y, family = "binomial", alpha = 1)

# Get coefficients from the best model
lasso_coef <- coef(lasso_model, s = "lambda.min")
selected_biomarkers <- rownames(lasso_coef)[lasso_coef[,1] != 0]
print(selected_biomarkers)

# Predict probabilities using the selected features
pred_probs <- predict(lasso_model, newx = x, s = "lambda.min", type = "response")

# Calculate AUC ROC
roc_curve <- roc(true_outcomes, pred_probs)
auc_biomarkers <- roc_curve$auc
print(paste("AUC ROC for Random Forest model:", auc_biomarkers))

# Calculate 95% CI for the AUC
ci_auc <- ci.auc(roc_curve)
print(paste("95% CI for AUC:", ci_auc[1], "-", ci_auc[3]))

# Perform DeLong's test to compare the AUCs
delong_test <- roc.test(auc_base, auc_biomarkers)
print(delong_test)

#### XGBoost ####---------------------------------------------------------------
new_colnames <- c("Age", "Sex", "BMI", "Smoker", "Systolic blood pressure", "Diabetes", 
                  "Prior stroke/TIA", "Heart failure", "Renal failure", "CAD", "Study center", 
                  "Ethnicity", "AF type", "Beer drinker", "Red wine drinker", "White wine drinker", 
                  "Liquor drinker", "AF episode duration", "Nr. of AF episodes", "Atrial flutter", 
                  "Aspirin", "TCA", "Antiplatelet therapy", "Illicit drugs", "Prior CABG", 
                  "Hyperthyroidism", "Hypothyroidism", "Cancer", "VTE", "Paternal AF", "Brother AF", 
                  "Maternal AF", "Sister AF", "Familial hypertension", "Familial diabetes", 
                  "Familial obesity", "Familial CAD", "Height", "Weight", "Heart rate", 
                  "Diastolic blood pressure", "Rhythm at baseline", "Sleep apnea", "Hypertension", 
                  "PAD", "Prior MI", "Systemic embolism", "Prior major bleeding", "ANG2", "DDI2H", 
                  "CYSC", "ALAT", "GDF15", "CRPHS", "IGFBP7", "IL6", "PROBNPII", "OPN", "TNTHS", "eGFR")

# Assign the new column names to the features data frame
colnames(features) <- new_colnames

features$heart.failure <- as.numeric(features$heart.failure == 1)

# Train XGBoost model
set.seed(123)
xgb_model <- xgboost(data = as.matrix(features[, -c("heart.failure")]), 
                     label = features$heart.failure, 
                     objective = "binary:logistic", 
                     nrounds = 10)

# Get variable importance
var_importance_xgb <- xgb.importance(feature_names = colnames(features[, -c("heart.failure")]), model = xgb_model)

# Get top biomarkers
top_biomarkers_xgb <- var_importance_xgb$Feature[1:12]
print(top_biomarkers_xgb)

# Get predicted probabilities
pred_probs_xgb <- as.numeric(predict(xgb_model, as.matrix(features[, -c("heart.failure")])))

# Get true outcomes
true_outcomes <- ifelse(features$heart.failure == 1, TRUE, FALSE)

# Calculate AUC ROC for XGBoost
roc_curve_xgb <- roc(true_outcomes, pred_probs_xgb)
auc_xgb_biomarkers <- roc_curve_xgb$auc
print(paste("AUC ROC for XGBoost model:", auc_xgb_biomarkers))

# Calculate 95% CI for the AUC
ci_auc_xgb <- ci.auc(roc_curve_xgb)
print(paste("95% CI for AUC:", ci_auc_xgb[1], "-", ci_auc_xgb[3]))

# Perform DeLong's test to compare the AUCs
delong_test <- roc.test(auc_xgb_base, auc_xgb_biomarkers)
print(delong_test)

# Plot variable importance
plot(var_importance_xgb, col = "black", pch = 20, cex = 1.25, top = 20)

ggplot(var_importance_xgb, aes(x = Gain, y = reorder(Feature, Gain))) +
  geom_point(shape = 1, size = 2.5, color = "black") +
  labs(x = "Importance", size = 2, y = "", color = "black") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(color = "black"),  # Set x-axis ticks color to black
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.major.y = element_line(linetype = "dotted", color = "grey"))


################################################################################
# Major bleeding
################################################################################

#-------------------------------------------------------------------------------
# Check whether there are missing "time" variables 
num_missing_timeto1.major.bleed <- sum(is.na(dat$timeto1.major.bleed))
print(num_missing_timeto1.major.bleed)
#-------------------------------------------------------------------------------

# Split data into features (biomarkers) and target variable
library(caret)
library(randomForest)
library(glmnet)
library(dplyr)
library(pROC)
library(xgboost)

# Exclude patients with missing variables
data <- dat %>% 
  select(ends_with("_log"), age.bl, pat.sex, bmi, current.smoker, rr.sys.liegend,
         prev.diabetes, prev.stroke.tia, prev.heart.failure, prev.niereninsuff,
         coronary.heart.disease, center, herkunft, vhf.typ.aktuell.bl, alkohol.bier, 
         alkohol.rotwein, alkohol.weisswein, alkohol.schnaps, vhf.episoden.dauer,
         vhf.episoden.anz, vorhofflattern, med.aspirin, med.tca.yn, med.antiplatelet.yn,
         drogen, prev.akb, prev.hyperthyreose, prev.hypothyreose, krk.malignom, krk.tvt,
         fam.vhf.vater, fam.vhf.bruder, fam.vhf.mutter, fam.vhf.schwester, fam.hypertonie,
         fam.diabetes, fam.uebergewicht, fam.khk, groesse, gewicht, heart.rate,
         rr.dia.liegend, ecg.rhythm.algo, prev.schlaf.apnoe, prev.hypertonie, prev.pavk,
         prev.mi, prev.sys.embolism, prev.major.bleed, major.bleed)

data <- na.omit(data)

# Make feature and target outcome data frames"
features <- data %>% 
  select(age.bl, pat.sex, bmi, current.smoker, rr.sys.liegend,
         prev.diabetes, prev.stroke.tia, prev.heart.failure, prev.niereninsuff,
         coronary.heart.disease, center, herkunft, vhf.typ.aktuell.bl, alkohol.bier, 
         alkohol.rotwein, alkohol.weisswein, alkohol.schnaps, vhf.episoden.dauer,
         vhf.episoden.anz, vorhofflattern, med.aspirin, med.tca.yn, med.antiplatelet.yn,
         drogen, prev.akb, prev.hyperthyreose, prev.hypothyreose, krk.malignom, krk.tvt,
         fam.vhf.vater, fam.vhf.bruder, fam.vhf.mutter, fam.vhf.schwester, fam.hypertonie,
         fam.diabetes, fam.uebergewicht, fam.khk, groesse, gewicht, heart.rate,
         rr.dia.liegend, ecg.rhythm.algo, prev.schlaf.apnoe, prev.hypertonie, prev.pavk,
         prev.mi, prev.sys.embolism, prev.major.bleed, major.bleed)

features$major.bleed <- as.factor(features$major.bleed)

#### Base model ################################################################
#### Random forest ####---------------------------------------------------------
# Train random forest model
set.seed(123)

# Train a Random Forest model
rf_model <- randomForest(major.bleed ~ ., data = features, importance = TRUE)

# Get variable importance
importance <- importance(rf_model)
var_importance <- data.frame(Variables = rownames(importance), Importance = importance[, 'MeanDecreaseGini'])

# Get top biomarkers
top_biomarkers <- var_importance %>% arrange(desc(Importance)) %>% head(12)
print(top_biomarkers)

# Get predicted probabilities
pred_probs <- predict(rf_model, type = "prob")[, "1"]

# Get true outcomes
true_outcomes <- ifelse(features$major.bleed == 1, TRUE, FALSE)

# Calculate AUC ROC
roc_curve <- roc(true_outcomes, pred_probs)
auc_base <- roc_curve$auc
print(paste("AUC ROC for Random Forest model:", auc_base))

# Calculate 95% CI for the AUC
ci_auc <- ci.auc(roc_curve)
print(paste("95% CI for AUC:", ci_auc[1], "-", ci_auc[3]))

#### LASSO ####-----------------------------------------------------------------
# Prepare data for LASSO
x <- model.matrix(major.bleed ~ ., features)[, -1]
y <- features$major.bleed

# Train LASSO model
lasso_model <- cv.glmnet(x, y, family = "binomial", alpha = 1)

# Get coefficients from the best model
lasso_coef <- coef(lasso_model, s = "lambda.min")
selected_biomarkers <- rownames(lasso_coef)[lasso_coef[,1] != 0]
print(selected_biomarkers)

# Predict probabilities using the selected features
pred_probs <- predict(lasso_model, newx = x, s = "lambda.min", type = "response")

# Calculate AUC ROC
roc_curve <- roc(y, as.vector(pred_probs))
auc_base <- roc_curve$auc
print(paste("AUC ROC for Random Forest model:", auc_base))

# Calculate 95% CI for the AUC
ci_auc <- ci.auc(roc_curve)
print(paste("95% CI for AUC:", ci_auc[1], "-", ci_auc[3]))

#### XGBoost ####---------------------------------------------------------------
new_colnames <- c("Age", "Sex", "BMI", "Smoker", "Systolic blood pressure", "Diabetes", 
                  "Prior stroke/TIA", "Heart failure", "Renal failure", "CAD", "Study center", 
                  "Ethnicity", "AF type", "Beer drinker", "Red wine drinker", "White wine drinker", 
                  "Liquor drinker", "AF episode duration", "Nr. of AF episodes", "Atrial flutter", 
                  "Aspirin", "TCA", "Antiplatelet therapy", "Illicit drugs", "Prior CABG", 
                  "Hyperthyroidism", "Hypothyroidism", "Cancer", "VTE", "Paternal AF", "Brother AF", 
                  "Maternal AF", "Sister AF", "Familial hypertension", "Familial diabetes", 
                  "Familial obesity", "Familial CAD", "Height", "Weight", "Heart rate", 
                  "Diastolic blood pressure", "Rhythm at baseline", "Sleep apnea", "Hypertension", 
                  "PAD", "Prior MI", "Systemic embolism", "Prior major bleeding")

# Assign the new column names to the features data frame
colnames(features) <- new_colnames

features$major.bleed <- as.numeric(features$major.bleed == 1)

# Train XGBoost model
set.seed(123)
xgb_model <- xgboost(data = as.matrix(features[, -c("major.bleed")]), 
                     label = features$major.bleed, 
                     objective = "binary:logistic", 
                     nrounds = 10)

# Get variable importance
var_importance_xgb <- xgb.importance(feature_names = colnames(features[, -c("major.bleed")]), model = xgb_model)

# Get top biomarkers
top_biomarkers_xgb <- var_importance_xgb$Feature[1:12]
print(top_biomarkers_xgb)

# Get predicted probabilities
pred_probs_xgb <- as.numeric(predict(xgb_model, as.matrix(features[, -c("major.bleed")])))

# Get true outcomes
true_outcomes <- ifelse(features$major.bleed == 1, TRUE, FALSE)

# Calculate AUC ROC for XGBoost
roc_curve_xgb <- roc(true_outcomes, pred_probs_xgb)
auc_xgb_base <- roc_curve_xgb$auc
print(paste("AUC ROC for XGBoost model:", auc_xgb_base))

# Calculate 95% CI for the AUC
ci_auc_xgb <- ci.auc(roc_curve_xgb)
print(paste("95% CI for AUC:", ci_auc_xgb[1], "-", ci_auc_xgb[3]))

# Plot variable importance
plot(var_importance_xgb, col = "black", pch = 20, cex = 1.25, top = 20)
xgb.plot.importance(var_importance_xgb, top_n = 20, measure = "Gain", rel_to_first = TRUE, xlab = "Relative Importance")

ggplot(var_importance_xgb, aes(x = reorder(Feature, Gain), y = Gain)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Variable Importance from XGBoost Model",
       x = "Features",
       y = "Gain") +
  theme_minimal()

ggplot(var_importance_xgb, aes(x = Gain, y = reorder(Feature, Gain))) +
  geom_point(shape = 1, size = 2.5, color = "black") +
  labs(x = "Importance", size = 2, y = "", color = "black") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(color = "black"),  # Set x-axis ticks color to black
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.major.y = element_line(linetype = "dotted", color = "grey"))

#### Model with biomarkers #####################################################
#### Random forest ####---------------------------------------------------------
library(caret)
library(randomForest)
library(glmnet)
library(dplyr)
library(pROC)
library(xgboost)

# Exclude patients with missing variables
data <- dat %>% 
  select(ends_with("_log"), age.bl, pat.sex, bmi, current.smoker, rr.sys.liegend,
         prev.diabetes, prev.stroke.tia, prev.heart.failure, prev.niereninsuff,
         coronary.heart.disease, center, herkunft, vhf.typ.aktuell.bl, alkohol.bier, 
         alkohol.rotwein, alkohol.weisswein, alkohol.schnaps, vhf.episoden.dauer,
         vhf.episoden.anz, vorhofflattern, med.aspirin, med.tca.yn, med.antiplatelet.yn,
         drogen, prev.akb, prev.hyperthyreose, prev.hypothyreose, krk.malignom, krk.tvt,
         fam.vhf.vater, fam.vhf.bruder, fam.vhf.mutter, fam.vhf.schwester, fam.hypertonie,
         fam.diabetes, fam.uebergewicht, fam.khk, groesse, gewicht, heart.rate,
         rr.dia.liegend, ecg.rhythm.algo, prev.schlaf.apnoe, prev.hypertonie, prev.pavk,
         prev.mi, prev.sys.embolism, prev.major.bleed, major.bleed)

data <- na.omit(data)

# Make feature and target outcome data frames"
features <- data %>% 
  select(age.bl, pat.sex, bmi, current.smoker, rr.sys.liegend,
         prev.diabetes, prev.stroke.tia, prev.heart.failure, prev.niereninsuff,
         coronary.heart.disease, center, herkunft, vhf.typ.aktuell.bl, alkohol.bier, 
         alkohol.rotwein, alkohol.weisswein, alkohol.schnaps, vhf.episoden.dauer,
         vhf.episoden.anz, vorhofflattern, med.aspirin, med.tca.yn, med.antiplatelet.yn,
         drogen, prev.akb, prev.hyperthyreose, prev.hypothyreose, krk.malignom, krk.tvt,
         fam.vhf.vater, fam.vhf.bruder, fam.vhf.mutter, fam.vhf.schwester, fam.hypertonie,
         fam.diabetes, fam.uebergewicht, fam.khk, groesse, gewicht, heart.rate,
         rr.dia.liegend, ecg.rhythm.algo, prev.schlaf.apnoe, prev.hypertonie, prev.pavk,
         prev.mi, prev.sys.embolism, prev.major.bleed, ends_with("_log"), major.bleed)

features$major.bleed <- as.factor(features$major.bleed)

# Train random forest model
set.seed(123)

# Train a Random Forest model
rf_model <- randomForest(major.bleed ~ ., data = features, importance = TRUE)

# Get variable importance
importance <- importance(rf_model)
var_importance <- data.frame(Variables = rownames(importance), Importance = importance[, 'MeanDecreaseGini'])

# Get top biomarkers
top_biomarkers <- var_importance %>% arrange(desc(Importance)) %>% head(12)
print(top_biomarkers)

# Get predicted probabilities
pred_probs <- predict(rf_model, type = "prob")[, "1"]

# Get true outcomes
true_outcomes <- ifelse(features$major.bleed == 1, TRUE, FALSE)

# Calculate AUC ROC
roc_curve <- roc(true_outcomes, pred_probs)
auc_biomarkers <- roc_curve$auc
print(paste("AUC ROC for Random Forest model:", auc_biomarkers))

# Calculate 95% CI for the AUC
ci_auc <- ci.auc(roc_curve)
print(paste("95% CI for AUC:", ci_auc[1], "-", ci_auc[3]))

# Perform DeLong's test to compare the AUCs
delong_test <- roc.test(auc_base, auc_biomarkers)
print(delong_test)

#### LASSO ####
# Prepare data for LASSO
x <- model.matrix(major.bleed ~ ., features)[, -1]
y <- features$major.bleed

# Train LASSO model
lasso_model <- cv.glmnet(x, y, family = "binomial", alpha = 1)

# Get coefficients from the best model
lasso_coef <- coef(lasso_model, s = "lambda.min")
selected_biomarkers <- rownames(lasso_coef)[lasso_coef[,1] != 0]
print(selected_biomarkers)

# Predict probabilities using the selected features
pred_probs <- predict(lasso_model, newx = x, s = "lambda.min", type = "response")

# Calculate AUC ROC
roc_curve <- roc(true_outcomes, pred_probs)
auc_biomarkers <- roc_curve$auc
print(paste("AUC ROC for Random Forest model:", auc_biomarkers))

# Calculate 95% CI for the AUC
ci_auc <- ci.auc(roc_curve)
print(paste("95% CI for AUC:", ci_auc[1], "-", ci_auc[3]))

# Perform DeLong's test to compare the AUCs
delong_test <- roc.test(auc_base, auc_biomarkers)
print(delong_test)

#### XGBoost ####---------------------------------------------------------------
new_colnames <- c("Age", "Sex", "BMI", "Smoker", "Systolic blood pressure", "Diabetes", 
                  "Prior stroke/TIA", "Heart failure", "Renal failure", "CAD", "Study center", 
                  "Ethnicity", "AF type", "Beer drinker", "Red wine drinker", "White wine drinker", 
                  "Liquor drinker", "AF episode duration", "Nr. of AF episodes", "Atrial flutter", 
                  "Aspirin", "TCA", "Antiplatelet therapy", "Illicit drugs", "Prior CABG", 
                  "Hyperthyroidism", "Hypothyroidism", "Cancer", "VTE", "Paternal AF", "Brother AF", 
                  "Maternal AF", "Sister AF", "Familial hypertension", "Familial diabetes", 
                  "Familial obesity", "Familial CAD", "Height", "Weight", "Heart rate", 
                  "Diastolic blood pressure", "Rhythm at baseline", "Sleep apnea", "Hypertension", 
                  "PAD", "Prior MI", "Systemic embolism", "Prior major bleeding", 
                  "ANG2", "DDI2H", "CYSC", "ALAT", "GDF15", "CRPHS", "IGFBP7", "IL6", 
                  "PROBNPII", "OPN", "TNTHS", "eGFR")

# Assign the new column names to the features data frame
colnames(features) <- new_colnames

features$major.bleed <- as.numeric(features$major.bleed == 1)

# Train XGBoost model
set.seed(123)
xgb_model <- xgboost(data = as.matrix(features[, -c("major.bleed")]), 
                     label = features$major.bleed, 
                     objective = "binary:logistic", 
                     nrounds = 10)

# Get variable importance
var_importance_xgb <- xgb.importance(feature_names = colnames(features[, -c("major.bleed")]), model = xgb_model)

# Get top biomarkers
top_biomarkers_xgb <- var_importance_xgb$Feature[1:12]
print(top_biomarkers_xgb)

# Get predicted probabilities
pred_probs_xgb <- as.numeric(predict(xgb_model, as.matrix(features[, -c("major.bleed")])))

# Get true outcomes
true_outcomes <- ifelse(features$major.bleed == 1, TRUE, FALSE)

# Calculate AUC ROC for XGBoost
roc_curve_xgb <- roc(true_outcomes, pred_probs_xgb)
auc_xgb_biomarkers <- roc_curve_xgb$auc
print(paste("AUC ROC for XGBoost model:", auc_xgb_biomarkers))

# Calculate 95% CI for the AUC
ci_auc_xgb <- ci.auc(roc_curve_xgb)
print(paste("95% CI for AUC:", ci_auc_xgb[1], "-", ci_auc_xgb[3]))

# Perform DeLong's test to compare the AUCs
delong_test <- roc.test(auc_xgb_base, auc_xgb_biomarkers)
print(delong_test)

# Plot variable importance
plot(var_importance_xgb, col = "black", pch = 20, cex = 1.25, top = 20)

ggplot(var_importance_xgb, aes(x = Gain, y = reorder(Feature, Gain))) +
  geom_point(shape = 1, size = 2.5, color = "black") +
  labs(x = "Importance", size = 2, y = "", color = "black") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(color = "black"),  # Set x-axis ticks color to black
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.major.y = element_line(linetype = "dotted", color = "grey"))


################################################################################
# Ischemic stroke    
################################################################################

#-------------------------------------------------------------------------------
# Check whether there are missing "time" variables 
num_missing_timeto1.isch.stroke <- sum(is.na(dat$timeto1.isch.stroke))
print(num_missing_timeto1.isch.stroke)
#-------------------------------------------------------------------------------

# Split data into features (biomarkers) and target variable
library(caret)
library(randomForest)
library(glmnet)
library(dplyr)
library(pROC)
library(xgboost)

# Exclude patients with missing variables
data <- dat %>% 
  select(ends_with("_log"), age.bl, pat.sex, bmi, current.smoker, rr.sys.liegend,
         prev.diabetes, prev.stroke.tia, prev.heart.failure, prev.niereninsuff,
         coronary.heart.disease, center, herkunft, vhf.typ.aktuell.bl, alkohol.bier, 
         alkohol.rotwein, alkohol.weisswein, alkohol.schnaps, vhf.episoden.dauer,
         vhf.episoden.anz, vorhofflattern, med.aspirin, med.tca.yn, med.antiplatelet.yn,
         drogen, prev.akb, prev.hyperthyreose, prev.hypothyreose, krk.malignom, krk.tvt,
         fam.vhf.vater, fam.vhf.bruder, fam.vhf.mutter, fam.vhf.schwester, fam.hypertonie,
         fam.diabetes, fam.uebergewicht, fam.khk, groesse, gewicht, heart.rate,
         rr.dia.liegend, ecg.rhythm.algo, prev.schlaf.apnoe, prev.hypertonie, prev.pavk,
         prev.mi, prev.sys.embolism, prev.major.bleed, ischemic.stroke)

data <- na.omit(data)

# Make feature and target outcome data frames"
features <- data %>% 
  select(age.bl, pat.sex, bmi, current.smoker, rr.sys.liegend,
         prev.diabetes, prev.stroke.tia, prev.heart.failure, prev.niereninsuff,
         coronary.heart.disease, center, herkunft, vhf.typ.aktuell.bl, alkohol.bier, 
         alkohol.rotwein, alkohol.weisswein, alkohol.schnaps, vhf.episoden.dauer,
         vhf.episoden.anz, vorhofflattern, med.aspirin, med.tca.yn, med.antiplatelet.yn,
         drogen, prev.akb, prev.hyperthyreose, prev.hypothyreose, krk.malignom, krk.tvt,
         fam.vhf.vater, fam.vhf.bruder, fam.vhf.mutter, fam.vhf.schwester, fam.hypertonie,
         fam.diabetes, fam.uebergewicht, fam.khk, groesse, gewicht, heart.rate,
         rr.dia.liegend, ecg.rhythm.algo, prev.schlaf.apnoe, prev.hypertonie, prev.pavk,
         prev.mi, prev.sys.embolism, prev.major.bleed, ischemic.stroke)

features$ischemic.stroke <- as.factor(features$ischemic.stroke)

#### Base model ################################################################
#### Random forest ####---------------------------------------------------------
# Train random forest model
set.seed(123)

# Train a Random Forest model
rf_model <- randomForest(ischemic.stroke ~ ., data = features, importance = TRUE)

# Get variable importance
importance <- importance(rf_model)
var_importance <- data.frame(Variables = rownames(importance), Importance = importance[, 'MeanDecreaseGini'])

# Get top biomarkers
top_biomarkers <- var_importance %>% arrange(desc(Importance)) %>% head(12)
print(top_biomarkers)

# Get predicted probabilities
pred_probs <- predict(rf_model, type = "prob")[, "1"]

# Get true outcomes
true_outcomes <- ifelse(features$ischemic.stroke == 1, TRUE, FALSE)

# Calculate AUC ROC
roc_curve <- roc(true_outcomes, pred_probs)
auc_base <- roc_curve$auc
print(paste("AUC ROC for Random Forest model:", auc_base))

# Calculate 95% CI for the AUC
ci_auc <- ci.auc(roc_curve)
print(paste("95% CI for AUC:", ci_auc[1], "-", ci_auc[3]))

#### LASSO ####-----------------------------------------------------------------
# Prepare data for LASSO
x <- model.matrix(ischemic.stroke ~ ., features)[, -1]
y <- features$ischemic.stroke

# Train LASSO model
lasso_model <- cv.glmnet(x, y, family = "binomial", alpha = 1)

# Get coefficients from the best model
lasso_coef <- coef(lasso_model, s = "lambda.min")
selected_biomarkers <- rownames(lasso_coef)[lasso_coef[,1] != 0]
print(selected_biomarkers)

# Predict probabilities using the selected features
pred_probs <- predict(lasso_model, newx = x, s = "lambda.min", type = "response")

# Calculate AUC ROC
roc_curve <- roc(y, as.vector(pred_probs))
auc_base <- roc_curve$auc
print(paste("AUC ROC for Random Forest model:", auc_base))

# Calculate 95% CI for the AUC
ci_auc <- ci.auc(roc_curve)
print(paste("95% CI for AUC:", ci_auc[1], "-", ci_auc[3]))

#### XGBoost ####---------------------------------------------------------------
new_colnames <- c("Age", "Sex", "BMI", "Smoker", "Systolic blood pressure", "Diabetes", 
                  "Prior stroke/TIA", "Heart failure", "Renal failure", "CAD", "Study center", 
                  "Ethnicity", "AF type", "Beer drinker", "Red wine drinker", "White wine drinker", 
                  "Liquor drinker", "AF episode duration", "Nr. of AF episodes", "Atrial flutter", 
                  "Aspirin", "TCA", "Antiplatelet therapy", "Illicit drugs", "Prior CABG", 
                  "Hyperthyroidism", "Hypothyroidism", "Cancer", "VTE", "Paternal AF", "Brother AF", 
                  "Maternal AF", "Sister AF", "Familial hypertension", "Familial diabetes", 
                  "Familial obesity", "Familial CAD", "Height", "Weight", "Heart rate", 
                  "Diastolic blood pressure", "Rhythm at baseline", "Sleep apnea", "Hypertension", 
                  "PAD", "Prior MI", "Systemic embolism", "Prior major bleeding")

# Assign the new column names to the dataframe
colnames(features) <- new_colnames

features$ischemic.stroke <- as.numeric(features$ischemic.stroke == 1)

# Train XGBoost model
set.seed(123)
xgb_model <- xgboost(data = as.matrix(features[, -c("ischemic.stroke")]), 
                     label = features$ischemic.stroke, 
                     objective = "binary:logistic", 
                     nrounds = 10)

# Get variable importance
var_importance_xgb <- xgb.importance(feature_names = colnames(features[, -c("ischemic.stroke")]), model = xgb_model)

# Get top biomarkers
top_biomarkers_xgb <- var_importance_xgb$Feature[1:12]
print(top_biomarkers_xgb)

# Get predicted probabilities
pred_probs_xgb <- as.numeric(predict(xgb_model, as.matrix(features[, -c("ischemic.stroke")])))

# Get true outcomes
true_outcomes <- ifelse(features$ischemic.stroke == 1, TRUE, FALSE)

# Calculate AUC ROC for XGBoost
roc_curve_xgb <- roc(true_outcomes, pred_probs_xgb)
auc_xgb_base <- roc_curve_xgb$auc
print(paste("AUC ROC for XGBoost model:", auc_xgb_base))

# Calculate 95% CI for the AUC
ci_auc_xgb <- ci.auc(roc_curve_xgb)
print(paste("95% CI for AUC:", ci_auc_xgb[1], "-", ci_auc_xgb[3]))

# Plot variable importance
plot(var_importance_xgb, col = "black", pch = 20, cex = 1.25, top = 20)
xgb.plot.importance(var_importance_xgb, top_n = 20, measure = "Gain", rel_to_first = TRUE, xlab = "Relative Importance")

ggplot(var_importance_xgb, aes(x = reorder(Feature, Gain), y = Gain)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Variable Importance from XGBoost Model",
       x = "Features",
       y = "Gain") +
  theme_minimal()

ggplot(var_importance_xgb, aes(x = Gain, y = reorder(Feature, Gain))) +
  geom_point(shape = 1, size = 2.5, color = "black") +
  labs(x = "Importance", size = 2, y = "", color = "black") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(color = "black"),  # Set x-axis ticks color to black
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.major.y = element_line(linetype = "dotted", color = "grey"))

#### Model with biomarkers #####################################################
#### Random forest ####---------------------------------------------------------
library(caret)
library(randomForest)
library(glmnet)
library(dplyr)
library(pROC)
library(xgboost)

# Exclude patients with missing variables
data <- dat %>% 
  select(ends_with("_log"), age.bl, pat.sex, bmi, current.smoker, rr.sys.liegend,
         prev.diabetes, prev.stroke.tia, prev.heart.failure, prev.niereninsuff,
         coronary.heart.disease, center, herkunft, vhf.typ.aktuell.bl, alkohol.bier, 
         alkohol.rotwein, alkohol.weisswein, alkohol.schnaps, vhf.episoden.dauer,
         vhf.episoden.anz, vorhofflattern, med.aspirin, med.tca.yn, med.antiplatelet.yn,
         drogen, prev.akb, prev.hyperthyreose, prev.hypothyreose, krk.malignom, krk.tvt,
         fam.vhf.vater, fam.vhf.bruder, fam.vhf.mutter, fam.vhf.schwester, fam.hypertonie,
         fam.diabetes, fam.uebergewicht, fam.khk, groesse, gewicht, heart.rate,
         rr.dia.liegend, ecg.rhythm.algo, prev.schlaf.apnoe, prev.hypertonie, prev.pavk,
         prev.mi, prev.sys.embolism, prev.major.bleed, ischemic.stroke)

data <- na.omit(data)

# Make feature and target outcome data frames"
features <- data %>% 
  select(age.bl, pat.sex, bmi, current.smoker, rr.sys.liegend,
         prev.diabetes, prev.stroke.tia, prev.heart.failure, prev.niereninsuff,
         coronary.heart.disease, center, herkunft, vhf.typ.aktuell.bl, alkohol.bier, 
         alkohol.rotwein, alkohol.weisswein, alkohol.schnaps, vhf.episoden.dauer,
         vhf.episoden.anz, vorhofflattern, med.aspirin, med.tca.yn, med.antiplatelet.yn,
         drogen, prev.akb, prev.hyperthyreose, prev.hypothyreose, krk.malignom, krk.tvt,
         fam.vhf.vater, fam.vhf.bruder, fam.vhf.mutter, fam.vhf.schwester, fam.hypertonie,
         fam.diabetes, fam.uebergewicht, fam.khk, groesse, gewicht, heart.rate,
         rr.dia.liegend, ecg.rhythm.algo, prev.schlaf.apnoe, prev.hypertonie, prev.pavk,
         prev.mi, prev.sys.embolism, prev.major.bleed, ends_with("_log"), ischemic.stroke)

features$ischemic.stroke <- as.factor(features$ischemic.stroke)

# Train random forest model
set.seed(123)

# Train a Random Forest model
rf_model <- randomForest(ischemic.stroke ~ ., data = features, importance = TRUE)

# Get variable importance
importance <- importance(rf_model)
var_importance <- data.frame(Variables = rownames(importance), Importance = importance[, 'MeanDecreaseGini'])

# Get top biomarkers
top_biomarkers <- var_importance %>% arrange(desc(Importance)) %>% head(12)
print(top_biomarkers)

# Get predicted probabilities
pred_probs <- predict(rf_model, type = "prob")[, "1"]

# Get true outcomes
true_outcomes <- ifelse(features$ischemic.stroke == 1, TRUE, FALSE)

# Calculate AUC ROC
roc_curve <- roc(true_outcomes, pred_probs)
auc_biomarkers <- roc_curve$auc
print(paste("AUC ROC for Random Forest model:", auc_biomarkers))

# Calculate 95% CI for the AUC
ci_auc <- ci.auc(roc_curve)
print(paste("95% CI for AUC:", ci_auc[1], "-", ci_auc[3]))

# Perform DeLong's test to compare the AUCs
delong_test <- roc.test(auc_base, auc_biomarkers)
print(delong_test)

#### LASSO ####
# Prepare data for LASSO
x <- model.matrix(ischemic.stroke ~ ., features)[, -1]
y <- features$ischemic.stroke

# Train LASSO model
lasso_model <- cv.glmnet(x, y, family = "binomial", alpha = 1)

# Get coefficients from the best model
lasso_coef <- coef(lasso_model, s = "lambda.min")
selected_biomarkers <- rownames(lasso_coef)[lasso_coef[,1] != 0]
print(selected_biomarkers)

# Predict probabilities using the selected features
pred_probs <- predict(lasso_model, newx = x, s = "lambda.min", type = "response")

# Calculate AUC ROC
roc_curve <- roc(true_outcomes, pred_probs)
auc_biomarkers <- roc_curve$auc
print(paste("AUC ROC for Random Forest model:", auc_biomarkers))

# Calculate 95% CI for the AUC
ci_auc <- ci.auc(roc_curve)
print(paste("95% CI for AUC:", ci_auc[1], "-", ci_auc[3]))

# Perform DeLong's test to compare the AUCs
delong_test <- roc.test(auc_base, auc_biomarkers)
print(delong_test)

#### XGBoost ####---------------------------------------------------------------
new_colnames <- c("Age", "Sex", "BMI", "Smoker", "Systolic blood pressure", "Diabetes", 
                  "Prior stroke/TIA", "Heart failure", "Renal failure", "CAD", "Study center", 
                  "Ethnicity", "AF type", "Beer drinker", "Red wine drinker", "White wine drinker", 
                  "Liquor drinker", "AF episode duration", "Nr. of AF episodes", "Atrial flutter", 
                  "Aspirin", "TCA", "Antiplatelet therapy", "Illicit drugs", "Prior CABG", 
                  "Hyperthyroidism", "Hypothyroidism", "Cancer", "VTE", "Paternal AF", "Brother AF", 
                  "Maternal AF", "Sister AF", "Familial hypertension", "Familial diabetes", 
                  "Familial obesity", "Familial CAD", "Height", "Weight", "Heart rate", 
                  "Diastolic blood pressure", "Rhythm at baseline", "Sleep apnea", "Hypertension", 
                  "PAD", "Prior MI", "Systemic embolism", "Prior major bleeding", "ANG2", "DDI2H", 
                  "CYSC", "ALAT", "GDF15", "CRPHS", "IGFBP7", "IL6", "PROBNPII", "OPN", "TNTHS", "eGFR")

# Assign the new column names to the dataframe
colnames(features) <- new_colnames

features$ischemic.stroke <- as.numeric(features$ischemic.stroke == 1)

# Train XGBoost model
set.seed(123)
xgb_model <- xgboost(data = as.matrix(features[, -c("ischemic.stroke")]), 
                     label = features$ischemic.stroke, 
                     objective = "binary:logistic", 
                     nrounds = 10)

# Get variable importance
var_importance_xgb <- xgb.importance(feature_names = colnames(features[, -c("ischemic.stroke")]), model = xgb_model)

# Get top biomarkers
top_biomarkers_xgb <- var_importance_xgb$Feature[1:12]
print(top_biomarkers_xgb)

# Get predicted probabilities
pred_probs_xgb <- as.numeric(predict(xgb_model, as.matrix(features[, -c("ischemic.stroke")])))

# Get true outcomes
true_outcomes <- ifelse(features$ischemic.stroke == 1, TRUE, FALSE)

# Calculate AUC ROC for XGBoost
roc_curve_xgb <- roc(true_outcomes, pred_probs_xgb)
auc_xgb_biomarkers <- roc_curve_xgb$auc
print(paste("AUC ROC for XGBoost model:", auc_xgb_biomarkers))

# Calculate 95% CI for the AUC
ci_auc_xgb <- ci.auc(roc_curve_xgb)
print(paste("95% CI for AUC:", ci_auc_xgb[1], "-", ci_auc_xgb[3]))

# Perform DeLong's test to compare the AUCs
delong_test <- roc.test(auc_xgb_base, auc_xgb_biomarkers)
print(delong_test)

# Plot variable importance
plot(var_importance_xgb, col = "black", pch = 20, cex = 1.25, top = 20)

ggplot(var_importance_xgb, aes(x = Gain, y = reorder(Feature, Gain))) +
  geom_point(shape = 1, size = 2.5, color = "black") +
  labs(x = "Importance", size = 2, y = "", color = "black") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(color = "black"),  # Set x-axis ticks color to black
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.major.y = element_line(linetype = "dotted", color = "grey"))


################################################################################
# All strokes
################################################################################

#-------------------------------------------------------------------------------
# Check whether there are missing "time" variables 
num_missing_timeto1.stroke <- sum(is.na(dat$timeto1.stroke))
print(num_missing_timeto1.stroke)
#-------------------------------------------------------------------------------

# Split data into features (biomarkers) and target variable
library(caret)
library(randomForest)
library(glmnet)
library(dplyr)
library(pROC)
library(xgboost)

# Exclude patients with missing variables
data <- dat %>% 
  select(ends_with("_log"), age.bl, pat.sex, bmi, current.smoker, rr.sys.liegend,
         prev.diabetes, prev.stroke.tia, prev.heart.failure, prev.niereninsuff,
         coronary.heart.disease, center, herkunft, vhf.typ.aktuell.bl, alkohol.bier, 
         alkohol.rotwein, alkohol.weisswein, alkohol.schnaps, vhf.episoden.dauer,
         vhf.episoden.anz, vorhofflattern, med.aspirin, med.tca.yn, med.antiplatelet.yn,
         drogen, prev.akb, prev.hyperthyreose, prev.hypothyreose, krk.malignom, krk.tvt,
         fam.vhf.vater, fam.vhf.bruder, fam.vhf.mutter, fam.vhf.schwester, fam.hypertonie,
         fam.diabetes, fam.uebergewicht, fam.khk, groesse, gewicht, heart.rate,
         rr.dia.liegend, ecg.rhythm.algo, prev.schlaf.apnoe, prev.hypertonie, prev.pavk,
         prev.mi, prev.sys.embolism, prev.major.bleed, stroke)

data <- na.omit(data)

# Make feature and target outcome data frames"
features <- data %>% 
  select(age.bl, pat.sex, bmi, current.smoker, rr.sys.liegend,
         prev.diabetes, prev.stroke.tia, prev.heart.failure, prev.niereninsuff,
         coronary.heart.disease, center, herkunft, vhf.typ.aktuell.bl, alkohol.bier, 
         alkohol.rotwein, alkohol.weisswein, alkohol.schnaps, vhf.episoden.dauer,
         vhf.episoden.anz, vorhofflattern, med.aspirin, med.tca.yn, med.antiplatelet.yn,
         drogen, prev.akb, prev.hyperthyreose, prev.hypothyreose, krk.malignom, krk.tvt,
         fam.vhf.vater, fam.vhf.bruder, fam.vhf.mutter, fam.vhf.schwester, fam.hypertonie,
         fam.diabetes, fam.uebergewicht, fam.khk, groesse, gewicht, heart.rate,
         rr.dia.liegend, ecg.rhythm.algo, prev.schlaf.apnoe, prev.hypertonie, prev.pavk,
         prev.mi, prev.sys.embolism, prev.major.bleed, stroke)

features$stroke <- as.factor(features$stroke)

#### Base model ################################################################
#### Random forest ####---------------------------------------------------------
# Train random forest model
set.seed(123)

# Train a Random Forest model
rf_model <- randomForest(stroke ~ ., data = features, importance = TRUE)

# Get variable importance
importance <- importance(rf_model)
var_importance <- data.frame(Variables = rownames(importance), Importance = importance[, 'MeanDecreaseGini'])

# Get top biomarkers
top_biomarkers <- var_importance %>% arrange(desc(Importance)) %>% head(12)
print(top_biomarkers)

# Get predicted probabilities
pred_probs <- predict(rf_model, type = "prob")[, "1"]

# Get true outcomes
true_outcomes <- ifelse(features$stroke == 1, TRUE, FALSE)

# Calculate AUC ROC
roc_curve <- roc(true_outcomes, pred_probs)
auc_base <- roc_curve$auc
print(paste("AUC ROC for Random Forest model:", auc_base))

# Calculate 95% CI for the AUC
ci_auc <- ci.auc(roc_curve)
print(paste("95% CI for AUC:", ci_auc[1], "-", ci_auc[3]))

#### LASSO ####-----------------------------------------------------------------
# Prepare data for LASSO
x <- model.matrix(stroke ~ ., features)[, -1]
y <- features$stroke

# Train LASSO model
lasso_model <- cv.glmnet(x, y, family = "binomial", alpha = 1)

# Get coefficients from the best model
lasso_coef <- coef(lasso_model, s = "lambda.min")
selected_biomarkers <- rownames(lasso_coef)[lasso_coef[,1] != 0]
print(selected_biomarkers)

# Predict probabilities using the selected features
pred_probs <- predict(lasso_model, newx = x, s = "lambda.min", type = "response")

# Calculate AUC ROC
roc_curve <- roc(y, as.vector(pred_probs))
auc_base <- roc_curve$auc
print(paste("AUC ROC for Random Forest model:", auc_base))

# Calculate 95% CI for the AUC
ci_auc <- ci.auc(roc_curve)
print(paste("95% CI for AUC:", ci_auc[1], "-", ci_auc[3]))

#### XGBoost ####---------------------------------------------------------------
new_colnames <- c(
  "Age", "Sex", "BMI", "Smoker", "Systolic blood pressure", "Diabetes", 
  "Prior stroke/TIA", "Heart failure", "Renal failure", "CAD", "Study center", 
  "Ethnicity", "AF type", "Beer drinker", "Red wine drinker", "White wine drinker", 
  "Liquor drinker", "AF episode duration", "Nr. of AF episodes", "Atrial flutter", 
  "Aspirin", "TCA", "Antiplatelet therapy", "Illicit drugs", "Prior CABG", 
  "Hyperthyroidism", "Hypothyroidism", "Cancer", "VTE", "Paternal AF", "Brother AF", 
  "Maternal AF", "Sister AF", "Familial hypertension", "Familial diabetes", 
  "Familial obesity", "Familial CAD", "Height", "Weight", "Heart rate", 
  "Diastolic blood pressure", "Rhythm at baseline", "Sleep apnea", "Hypertension", 
  "PAD", "Prior MI", "Systemic embolism", "Prior major bleeding"
)

# Assign the new column names to the dataframe
colnames(features) <- new_colnames

features$stroke <- as.numeric(features$stroke == 1)

# Train XGBoost model
set.seed(123)
xgb_model <- xgboost(data = as.matrix(features[, -c("stroke")]), 
                     label = features$stroke, 
                     objective = "binary:logistic", 
                     nrounds = 10)

# Get variable importance
var_importance_xgb <- xgb.importance(feature_names = colnames(features[, -c("stroke")]), model = xgb_model)

# Get top biomarkers
top_biomarkers_xgb <- var_importance_xgb$Feature[1:12]
print(top_biomarkers_xgb)

# Get predicted probabilities
pred_probs_xgb <- as.numeric(predict(xgb_model, as.matrix(features[, -c("stroke")])))

# Get true outcomes
true_outcomes <- ifelse(features$stroke == 1, TRUE, FALSE)

# Calculate AUC ROC for XGBoost
roc_curve_xgb <- roc(true_outcomes, pred_probs_xgb)
auc_xgb_base <- roc_curve_xgb$auc
print(paste("AUC ROC for XGBoost model:", auc_xgb_base))

# Calculate 95% CI for the AUC
ci_auc_xgb <- ci.auc(roc_curve_xgb)
print(paste("95% CI for AUC:", ci_auc_xgb[1], "-", ci_auc_xgb[3]))

# Plot variable importance
plot(var_importance_xgb, col = "black", pch = 20, cex = 1.25, top = 20)
xgb.plot.importance(var_importance_xgb, top_n = 20, measure = "Gain", rel_to_first = TRUE, xlab = "Relative Importance")

ggplot(var_importance_xgb, aes(x = reorder(Feature, Gain), y = Gain)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Variable Importance from XGBoost Model",
       x = "Features",
       y = "Gain") +
  theme_minimal()

ggplot(var_importance_xgb, aes(x = Gain, y = reorder(Feature, Gain))) +
  geom_point(shape = 1, size = 2.5, color = "black") +
  labs(x = "Importance", size = 2, y = "", color = "black") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(color = "black"),  # Set x-axis ticks color to black
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.major.y = element_line(linetype = "dotted", color = "grey"))

#### Model with biomarkers #####################################################
#### Random forest ####---------------------------------------------------------
library(caret)
library(randomForest)
library(glmnet)
library(dplyr)
library(pROC)
library(xgboost)

# Exclude patients with missing variables
data <- dat %>% 
  select(ends_with("_log"), age.bl, pat.sex, bmi, current.smoker, rr.sys.liegend,
         prev.diabetes, prev.stroke.tia, prev.heart.failure, prev.niereninsuff,
         coronary.heart.disease, center, herkunft, vhf.typ.aktuell.bl, alkohol.bier, 
         alkohol.rotwein, alkohol.weisswein, alkohol.schnaps, vhf.episoden.dauer,
         vhf.episoden.anz, vorhofflattern, med.aspirin, med.tca.yn, med.antiplatelet.yn,
         drogen, prev.akb, prev.hyperthyreose, prev.hypothyreose, krk.malignom, krk.tvt,
         fam.vhf.vater, fam.vhf.bruder, fam.vhf.mutter, fam.vhf.schwester, fam.hypertonie,
         fam.diabetes, fam.uebergewicht, fam.khk, groesse, gewicht, heart.rate,
         rr.dia.liegend, ecg.rhythm.algo, prev.schlaf.apnoe, prev.hypertonie, prev.pavk,
         prev.mi, prev.sys.embolism, prev.major.bleed, stroke)

data <- na.omit(data)

# Make feature and target outcome data frames"
features <- data %>% 
  select(age.bl, pat.sex, bmi, current.smoker, rr.sys.liegend,
         prev.diabetes, prev.stroke.tia, prev.heart.failure, prev.niereninsuff,
         coronary.heart.disease, center, herkunft, vhf.typ.aktuell.bl, alkohol.bier, 
         alkohol.rotwein, alkohol.weisswein, alkohol.schnaps, vhf.episoden.dauer,
         vhf.episoden.anz, vorhofflattern, med.aspirin, med.tca.yn, med.antiplatelet.yn,
         drogen, prev.akb, prev.hyperthyreose, prev.hypothyreose, krk.malignom, krk.tvt,
         fam.vhf.vater, fam.vhf.bruder, fam.vhf.mutter, fam.vhf.schwester, fam.hypertonie,
         fam.diabetes, fam.uebergewicht, fam.khk, groesse, gewicht, heart.rate,
         rr.dia.liegend, ecg.rhythm.algo, prev.schlaf.apnoe, prev.hypertonie, prev.pavk,
         prev.mi, prev.sys.embolism, prev.major.bleed, ends_with("_log"), stroke)

features$stroke <- as.factor(features$stroke)

# Train random forest model
set.seed(123)

# Train a Random Forest model
rf_model <- randomForest(stroke ~ ., data = features, importance = TRUE)

# Get variable importance
importance <- importance(rf_model)
var_importance <- data.frame(Variables = rownames(importance), Importance = importance[, 'MeanDecreaseGini'])

# Get top biomarkers
top_biomarkers <- var_importance %>% arrange(desc(Importance)) %>% head(12)
print(top_biomarkers)

# Get predicted probabilities
pred_probs <- predict(rf_model, type = "prob")[, "1"]

# Get true outcomes
true_outcomes <- ifelse(features$stroke == 1, TRUE, FALSE)

# Calculate AUC ROC
roc_curve <- roc(true_outcomes, pred_probs)
auc_biomarkers <- roc_curve$auc
print(paste("AUC ROC for Random Forest model:", auc_biomarkers))

# Calculate 95% CI for the AUC
ci_auc <- ci.auc(roc_curve)
print(paste("95% CI for AUC:", ci_auc[1], "-", ci_auc[3]))

# Perform DeLong's test to compare the AUCs
delong_test <- roc.test(auc_base, auc_biomarkers)
print(delong_test)

#### LASSO ####
# Prepare data for LASSO
x <- model.matrix(stroke ~ ., features)[, -1]
y <- features$stroke

# Train LASSO model
lasso_model <- cv.glmnet(x, y, family = "binomial", alpha = 1)

# Get coefficients from the best model
lasso_coef <- coef(lasso_model, s = "lambda.min")
selected_biomarkers <- rownames(lasso_coef)[lasso_coef[,1] != 0]
print(selected_biomarkers)

# Predict probabilities using the selected features
pred_probs <- predict(lasso_model, newx = x, s = "lambda.min", type = "response")

# Calculate AUC ROC
roc_curve <- roc(true_outcomes, pred_probs)
auc_biomarkers <- roc_curve$auc
print(paste("AUC ROC for Random Forest model:", auc_biomarkers))

# Calculate 95% CI for the AUC
ci_auc <- ci.auc(roc_curve)
print(paste("95% CI for AUC:", ci_auc[1], "-", ci_auc[3]))

# Perform DeLong's test to compare the AUCs
delong_test <- roc.test(auc_base, auc_biomarkers)
print(delong_test)

#### XGBoost ####---------------------------------------------------------------
new_colnames <- c(
  "Age", "Sex", "BMI", "Smoker", "Systolic blood pressure", "Diabetes", 
  "Prior stroke/TIA", "Heart failure", "Renal failure", "CAD", "Study center", 
  "Ethnicity", "AF type", "Beer drinker", "Red wine drinker", "White wine drinker", 
  "Liquor drinker", "AF episode duration", "Nr. of AF episodes", "Atrial flutter", 
  "Aspirin", "TCA", "Antiplatelet therapy", "Illicit drugs", "Prior CABG", 
  "Hyperthyroidism", "Hypothyroidism", "Cancer", "VTE", "Paternal AF", "Brother AF", 
  "Maternal AF", "Sister AF", "Familial hypertension", "Familial diabetes", 
  "Familial obesity", "Familial CAD", "Height", "Weight", "Heart rate", 
  "Diastolic blood pressure", "Rhythm at baseline", "Sleep apnea", "Hypertension", 
  "PAD", "Prior MI", "Systemic embolism", "Prior major bleeding", 
  "ANG2", "DDI2H", "CYSC", "ALAT", "GDF15", "CRPHS", "IGFBP7", "IL6", 
  "PROBNPII", "OPN", "TNTHS", "eGFR"
)

# Assign the new column names to the dataframe
colnames(features) <- new_colnames

features$stroke <- as.numeric(features$stroke == 1)

# Train XGBoost model
set.seed(123)
xgb_model <- xgboost(data = as.matrix(features[, -c("stroke")]), 
                     label = features$stroke, 
                     objective = "binary:logistic", 
                     nrounds = 10)

# Get variable importance
var_importance_xgb <- xgb.importance(feature_names = colnames(features[, -c("stroke")]), model = xgb_model)

# Get top biomarkers
top_biomarkers_xgb <- var_importance_xgb$Feature[1:12]
print(top_biomarkers_xgb)

# Get predicted probabilities
pred_probs_xgb <- as.numeric(predict(xgb_model, as.matrix(features[, -c("stroke")])))

# Get true outcomes
true_outcomes <- ifelse(features$stroke == 1, TRUE, FALSE)

# Calculate AUC ROC for XGBoost
roc_curve_xgb <- roc(true_outcomes, pred_probs_xgb)
auc_xgb_biomarkers <- roc_curve_xgb$auc
print(paste("AUC ROC for XGBoost model:", auc_xgb_biomarkers))

# Calculate 95% CI for the AUC
ci_auc_xgb <- ci.auc(roc_curve_xgb)
print(paste("95% CI for AUC:", ci_auc_xgb[1], "-", ci_auc_xgb[3]))

# Perform DeLong's test to compare the AUCs
delong_test <- roc.test(auc_xgb_base, auc_xgb_biomarkers)
print(delong_test)

# Plot variable importance
plot(var_importance_xgb, col = "black", pch = 20, cex = 1.25, top = 20)

ggplot(var_importance_xgb, aes(x = Gain, y = reorder(Feature, Gain))) +
  geom_point(shape = 1, size = 2.5, color = "black") +
  labs(x = "Importance", size = 2, y = "", color = "black") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(color = "black"),  # Set x-axis ticks color to black
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.major.y = element_line(linetype = "dotted", color = "grey"))


################################################################################
# Myocardial infarction
################################################################################

#-------------------------------------------------------------------------------
# Check whether there are missing "time" variables 
num_missing_timeto1.mi <- sum(is.na(dat$timeto1.mi))
print(num_missing_timeto1.mi)
#-------------------------------------------------------------------------------

# Split data into features (biomarkers) and target variable
library(caret)
library(randomForest)
library(glmnet)
library(dplyr)
library(pROC)
library(xgboost)

# Exclude patients with missing variables
data <- dat %>% 
  select(ends_with("_log"), age.bl, pat.sex, bmi, current.smoker, rr.sys.liegend,
         prev.diabetes, prev.stroke.tia, prev.heart.failure, prev.niereninsuff,
         coronary.heart.disease, center, herkunft, vhf.typ.aktuell.bl, alkohol.bier, 
         alkohol.rotwein, alkohol.weisswein, alkohol.schnaps, vhf.episoden.dauer,
         vhf.episoden.anz, vorhofflattern, med.aspirin, med.tca.yn, med.antiplatelet.yn,
         drogen, prev.akb, prev.hyperthyreose, prev.hypothyreose, krk.malignom, krk.tvt,
         fam.vhf.vater, fam.vhf.bruder, fam.vhf.mutter, fam.vhf.schwester, fam.hypertonie,
         fam.diabetes, fam.uebergewicht, fam.khk, groesse, gewicht, heart.rate,
         rr.dia.liegend, ecg.rhythm.algo, prev.schlaf.apnoe, prev.hypertonie, prev.pavk,
         prev.mi, prev.sys.embolism, prev.major.bleed, mi)

data <- na.omit(data)

# Make feature and target outcome data frames"
features <- data %>% 
  select(age.bl, pat.sex, bmi, current.smoker, rr.sys.liegend,
         prev.diabetes, prev.stroke.tia, prev.heart.failure, prev.niereninsuff,
         coronary.heart.disease, center, herkunft, vhf.typ.aktuell.bl, alkohol.bier, 
         alkohol.rotwein, alkohol.weisswein, alkohol.schnaps, vhf.episoden.dauer,
         vhf.episoden.anz, vorhofflattern, med.aspirin, med.tca.yn, med.antiplatelet.yn,
         drogen, prev.akb, prev.hyperthyreose, prev.hypothyreose, krk.malignom, krk.tvt,
         fam.vhf.vater, fam.vhf.bruder, fam.vhf.mutter, fam.vhf.schwester, fam.hypertonie,
         fam.diabetes, fam.uebergewicht, fam.khk, groesse, gewicht, heart.rate,
         rr.dia.liegend, ecg.rhythm.algo, prev.schlaf.apnoe, prev.hypertonie, prev.pavk,
         prev.mi, prev.sys.embolism, prev.major.bleed, mi)

features$mi <- as.factor(features$mi)

#### Base model ################################################################
#### Random forest ####---------------------------------------------------------
# Train random forest model
set.seed(123)

# Train a Random Forest model
rf_model <- randomForest(mi ~ ., data = features, importance = TRUE)

# Get variable importance
importance <- importance(rf_model)
var_importance <- data.frame(Variables = rownames(importance), Importance = importance[, 'MeanDecreaseGini'])

# Get top biomarkers
top_biomarkers <- var_importance %>% arrange(desc(Importance)) %>% head(12)
print(top_biomarkers)

# Get predicted probabilities
pred_probs <- predict(rf_model, type = "prob")[, "1"]

# Get true outcomes
true_outcomes <- ifelse(features$mi == 1, TRUE, FALSE)

# Calculate AUC ROC
roc_curve <- roc(true_outcomes, pred_probs)
auc_base <- roc_curve$auc
print(paste("AUC ROC for Random Forest model:", auc_base))

# Calculate 95% CI for the AUC
ci_auc <- ci.auc(roc_curve)
print(paste("95% CI for AUC:", ci_auc[1], "-", ci_auc[3]))

#### LASSO ####-----------------------------------------------------------------
# Prepare data for LASSO
x <- model.matrix(mi ~ ., features)[, -1]
y <- features$mi

# Train LASSO model
lasso_model <- cv.glmnet(x, y, family = "binomial", alpha = 1)

# Get coefficients from the best model
lasso_coef <- coef(lasso_model, s = "lambda.min")
selected_biomarkers <- rownames(lasso_coef)[lasso_coef[,1] != 0]
print(selected_biomarkers)

# Predict probabilities using the selected features
pred_probs <- predict(lasso_model, newx = x, s = "lambda.min", type = "response")

# Calculate AUC ROC
roc_curve <- roc(y, as.vector(pred_probs))
auc_base <- roc_curve$auc
print(paste("AUC ROC for Random Forest model:", auc_base))

# Calculate 95% CI for the AUC
ci_auc <- ci.auc(roc_curve)
print(paste("95% CI for AUC:", ci_auc[1], "-", ci_auc[3]))

#### XGBoost ####---------------------------------------------------------------
new_colnames <- c(
  "Age", "Sex", "BMI", "Smoker", "Systolic blood pressure", "Diabetes", 
  "Prior stroke/TIA", "Heart failure", "Renal failure", "CAD", "Study center", 
  "Ethnicity", "AF type", "Beer drinker", "Red wine drinker", "White wine drinker", 
  "Liquor drinker", "AF episode duration", "Nr. of AF episodes", "Atrial flutter", 
  "Aspirin", "TCA", "Antiplatelet therapy", "Illicit drugs", "Prior CABG", 
  "Hyperthyroidism", "Hypothyroidism", "Cancer", "VTE", "Paternal AF", "Brother AF", 
  "Maternal AF", "Sister AF", "Familial hypertension", "Familial diabetes", 
  "Familial obesity", "Familial CAD", "Height", "Weight", "Heart rate", 
  "Diastolic blood pressure", "Rhythm at baseline", "Sleep apnea", "Hypertension", 
  "PAD", "Prior MI", "Systemic embolism", "Prior major bleeding"
)

# Assign the new column names to the dataframe
colnames(features) <- new_colnames

features$mi <- as.numeric(features$mi == 1)

# Train XGBoost model
set.seed(123)
xgb_model <- xgboost(data = as.matrix(features[, -c("mi")]), 
                     label = features$mi, 
                     objective = "binary:logistic", 
                     nrounds = 10)

# Get variable importance
var_importance_xgb <- xgb.importance(feature_names = colnames(features[, -c("mi")]), model = xgb_model)

# Get top biomarkers
top_biomarkers_xgb <- var_importance_xgb$Feature[1:12]
print(top_biomarkers_xgb)

# Get predicted probabilities
pred_probs_xgb <- as.numeric(predict(xgb_model, as.matrix(features[, -c("mi")])))

# Get true outcomes
true_outcomes <- ifelse(features$mi == 1, TRUE, FALSE)

# Calculate AUC ROC for XGBoost
roc_curve_xgb <- roc(true_outcomes, pred_probs_xgb)
auc_xgb_base <- roc_curve_xgb$auc
print(paste("AUC ROC for XGBoost model:", auc_xgb_base))

# Calculate 95% CI for the AUC
ci_auc_xgb <- ci.auc(roc_curve_xgb)
print(paste("95% CI for AUC:", ci_auc_xgb[1], "-", ci_auc_xgb[3]))

# Plot variable importance
plot(var_importance_xgb, col = "black", pch = 20, cex = 1.25, top = 20)
xgb.plot.importance(var_importance_xgb, top_n = 20, measure = "Gain", rel_to_first = TRUE, xlab = "Relative Importance")

ggplot(var_importance_xgb, aes(x = reorder(Feature, Gain), y = Gain)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Variable Importance from XGBoost Model",
       x = "Features",
       y = "Gain") +
  theme_minimal()

ggplot(var_importance_xgb, aes(x = Gain, y = reorder(Feature, Gain))) +
  geom_point(shape = 1, size = 2.5, color = "black") +
  labs(x = "Importance", size = 2, y = "", color = "black") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(color = "black"),  # Set x-axis ticks color to black
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.major.y = element_line(linetype = "dotted", color = "grey"))

#### Model with biomarkers #####################################################
#### Random forest ####---------------------------------------------------------
library(caret)
library(randomForest)
library(glmnet)
library(dplyr)
library(pROC)
library(xgboost)

# Exclude patients with missing variables
data <- dat %>% 
  select(ends_with("_log"), age.bl, pat.sex, bmi, current.smoker, rr.sys.liegend,
         prev.diabetes, prev.stroke.tia, prev.heart.failure, prev.niereninsuff,
         coronary.heart.disease, center, herkunft, vhf.typ.aktuell.bl, alkohol.bier, 
         alkohol.rotwein, alkohol.weisswein, alkohol.schnaps, vhf.episoden.dauer,
         vhf.episoden.anz, vorhofflattern, med.aspirin, med.tca.yn, med.antiplatelet.yn,
         drogen, prev.akb, prev.hyperthyreose, prev.hypothyreose, krk.malignom, krk.tvt,
         fam.vhf.vater, fam.vhf.bruder, fam.vhf.mutter, fam.vhf.schwester, fam.hypertonie,
         fam.diabetes, fam.uebergewicht, fam.khk, groesse, gewicht, heart.rate,
         rr.dia.liegend, ecg.rhythm.algo, prev.schlaf.apnoe, prev.hypertonie, prev.pavk,
         prev.mi, prev.sys.embolism, prev.major.bleed, mi)

data <- na.omit(data)

# Make feature and target outcome data frames"
features <- data %>% 
  select(age.bl, pat.sex, bmi, current.smoker, rr.sys.liegend,
         prev.diabetes, prev.stroke.tia, prev.heart.failure, prev.niereninsuff,
         coronary.heart.disease, center, herkunft, vhf.typ.aktuell.bl, alkohol.bier, 
         alkohol.rotwein, alkohol.weisswein, alkohol.schnaps, vhf.episoden.dauer,
         vhf.episoden.anz, vorhofflattern, med.aspirin, med.tca.yn, med.antiplatelet.yn,
         drogen, prev.akb, prev.hyperthyreose, prev.hypothyreose, krk.malignom, krk.tvt,
         fam.vhf.vater, fam.vhf.bruder, fam.vhf.mutter, fam.vhf.schwester, fam.hypertonie,
         fam.diabetes, fam.uebergewicht, fam.khk, groesse, gewicht, heart.rate,
         rr.dia.liegend, ecg.rhythm.algo, prev.schlaf.apnoe, prev.hypertonie, prev.pavk,
         prev.mi, prev.sys.embolism, prev.major.bleed, ends_with("_log"), mi)

features$mi <- as.factor(features$mi)

# Train random forest model
set.seed(123)

# Train a Random Forest model
rf_model <- randomForest(mi ~ ., data = features, importance = TRUE)

# Get variable importance
importance <- importance(rf_model)
var_importance <- data.frame(Variables = rownames(importance), Importance = importance[, 'MeanDecreaseGini'])

# Get top biomarkers
top_biomarkers <- var_importance %>% arrange(desc(Importance)) %>% head(12)
print(top_biomarkers)

# Get predicted probabilities
pred_probs <- predict(rf_model, type = "prob")[, "1"]

# Get true outcomes
true_outcomes <- ifelse(features$mi == 1, TRUE, FALSE)

# Calculate AUC ROC
roc_curve <- roc(true_outcomes, pred_probs)
auc_biomarkers <- roc_curve$auc
print(paste("AUC ROC for Random Forest model:", auc_biomarkers))

# Calculate 95% CI for the AUC
ci_auc <- ci.auc(roc_curve)
print(paste("95% CI for AUC:", ci_auc[1], "-", ci_auc[3]))

# Perform DeLong's test to compare the AUCs
delong_test <- roc.test(auc_base, auc_biomarkers)
print(delong_test)

#### LASSO ####
# Prepare data for LASSO
x <- model.matrix(mi ~ ., features)[, -1]
y <- features$mi

# Train LASSO model
lasso_model <- cv.glmnet(x, y, family = "binomial", alpha = 1)

# Get coefficients from the best model
lasso_coef <- coef(lasso_model, s = "lambda.min")
selected_biomarkers <- rownames(lasso_coef)[lasso_coef[,1] != 0]
print(selected_biomarkers)

# Predict probabilities using the selected features
pred_probs <- predict(lasso_model, newx = x, s = "lambda.min", type = "response")

# Calculate AUC ROC
roc_curve <- roc(true_outcomes, pred_probs)
auc_biomarkers <- roc_curve$auc
print(paste("AUC ROC for Random Forest model:", auc_biomarkers))

# Calculate 95% CI for the AUC
ci_auc <- ci.auc(roc_curve)
print(paste("95% CI for AUC:", ci_auc[1], "-", ci_auc[3]))

# Perform DeLong's test to compare the AUCs
delong_test <- roc.test(auc_base, auc_biomarkers)
print(delong_test)

#### XGBoost ####---------------------------------------------------------------
new_colnames <- c(
  "Age", "Sex", "BMI", "Smoker", "Systolic blood pressure", "Diabetes", 
  "Prior stroke/TIA", "Heart failure", "Renal failure", "CAD", "Study center", 
  "Ethnicity", "AF type", "Beer drinker", "Red wine drinker", "White wine drinker", 
  "Liquor drinker", "AF episode duration", "Nr. of AF episodes", "Atrial flutter", 
  "Aspirin", "TCA", "Antiplatelet therapy", "Illicit drugs", "Prior CABG", 
  "Hyperthyroidism", "Hypothyroidism", "Cancer", "VTE", "Paternal AF", "Brother AF", 
  "Maternal AF", "Sister AF", "Familial hypertension", "Familial diabetes", 
  "Familial obesity", "Familial CAD", "Height", "Weight", "Heart rate", 
  "Diastolic blood pressure", "Rhythm at baseline", "Sleep apnea", "Hypertension", 
  "PAD", "Prior MI", "Systemic embolism", "Prior major bleeding", 
  "ANG2", "DDI2H", "CYSC", "ALAT", "GDF15", "CRPHS", "IGFBP7", "IL6", 
  "PROBNPII", "OPN", "TNTHS", "eGFR"
)

# Assign the new column names to the dataframe
colnames(features) <- new_colnames

features$mi <- as.numeric(features$mi == 1)

# Train XGBoost model
set.seed(123)
xgb_model <- xgboost(data = as.matrix(features[, -c("mi")]), 
                     label = features$mi, 
                     objective = "binary:logistic", 
                     nrounds = 10)

# Get variable importance
var_importance_xgb <- xgb.importance(feature_names = colnames(features[, -c("mi")]), model = xgb_model)

# Get top biomarkers
top_biomarkers_xgb <- var_importance_xgb$Feature[1:12]
print(top_biomarkers_xgb)

# Get predicted probabilities
pred_probs_xgb <- as.numeric(predict(xgb_model, as.matrix(features[, -c("mi")])))

# Get true outcomes
true_outcomes <- ifelse(features$mi == 1, TRUE, FALSE)

# Calculate AUC ROC for XGBoost
roc_curve_xgb <- roc(true_outcomes, pred_probs_xgb)
auc_xgb_biomarkers <- roc_curve_xgb$auc
print(paste("AUC ROC for XGBoost model:", auc_xgb_biomarkers))

# Calculate 95% CI for the AUC
ci_auc_xgb <- ci.auc(roc_curve_xgb)
print(paste("95% CI for AUC:", ci_auc_xgb[1], "-", ci_auc_xgb[3]))

# Perform DeLong's test to compare the AUCs
delong_test <- roc.test(auc_xgb_base, auc_xgb_biomarkers)
print(delong_test)

# Plot variable importance
plot(var_importance_xgb, col = "black", pch = 20, cex = 1.25, top = 20)

ggplot(var_importance_xgb, aes(x = Gain, y = reorder(Feature, Gain))) +
  geom_point(shape = 1, size = 2.5, color = "black") +
  labs(x = "Importance", size = 2, y = "", color = "black") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(color = "black"),  # Set x-axis ticks color to black
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.major.y = element_line(linetype = "dotted", color = "grey"))


################################################################################
# Cardiovascular death
################################################################################

#-------------------------------------------------------------------------------
# Check whether there are missing "time" variables 
num_missing_timeto1.death <- sum(is.na(dat$timeto1.death))
print(num_missing_timeto1.death)
#-------------------------------------------------------------------------------

# Split data into features (biomarkers) and target variable
library(caret)
library(randomForest)
library(glmnet)
library(dplyr)
library(pROC)
library(xgboost)

# Exclude patients with missing variables
data <- dat %>% 
  select(ends_with("_log"), age.bl, pat.sex, bmi, current.smoker, rr.sys.liegend,
         prev.diabetes, prev.stroke.tia, prev.heart.failure, prev.niereninsuff,
         coronary.heart.disease, center, herkunft, vhf.typ.aktuell.bl, alkohol.bier, 
         alkohol.rotwein, alkohol.weisswein, alkohol.schnaps, vhf.episoden.dauer,
         vhf.episoden.anz, vorhofflattern, med.aspirin, med.tca.yn, med.antiplatelet.yn,
         drogen, prev.akb, prev.hyperthyreose, prev.hypothyreose, krk.malignom, krk.tvt,
         fam.vhf.vater, fam.vhf.bruder, fam.vhf.mutter, fam.vhf.schwester, fam.hypertonie,
         fam.diabetes, fam.uebergewicht, fam.khk, groesse, gewicht, heart.rate,
         rr.dia.liegend, ecg.rhythm.algo, prev.schlaf.apnoe, prev.hypertonie, prev.pavk,
         prev.mi, prev.sys.embolism, prev.major.bleed, death.cardiac)

data <- na.omit(data)

# Make feature and target outcome data frames"
features <- data %>% 
  select(age.bl, pat.sex, bmi, current.smoker, rr.sys.liegend,
         prev.diabetes, prev.stroke.tia, prev.heart.failure, prev.niereninsuff,
         coronary.heart.disease, center, herkunft, vhf.typ.aktuell.bl, alkohol.bier, 
         alkohol.rotwein, alkohol.weisswein, alkohol.schnaps, vhf.episoden.dauer,
         vhf.episoden.anz, vorhofflattern, med.aspirin, med.tca.yn, med.antiplatelet.yn,
         drogen, prev.akb, prev.hyperthyreose, prev.hypothyreose, krk.malignom, krk.tvt,
         fam.vhf.vater, fam.vhf.bruder, fam.vhf.mutter, fam.vhf.schwester, fam.hypertonie,
         fam.diabetes, fam.uebergewicht, fam.khk, groesse, gewicht, heart.rate,
         rr.dia.liegend, ecg.rhythm.algo, prev.schlaf.apnoe, prev.hypertonie, prev.pavk,
         prev.mi, prev.sys.embolism, prev.major.bleed, death.cardiac)

features$death.cardiac <- as.factor(features$death.cardiac)

#### Base model ################################################################
#### Random forest ####---------------------------------------------------------
# Train random forest model
set.seed(123)

# Train a Random Forest model
rf_model <- randomForest(death.cardiac ~ ., data = features, importance = TRUE)

# Get variable importance
importance <- importance(rf_model)
var_importance <- data.frame(Variables = rownames(importance), Importance = importance[, 'MeanDecreaseGini'])

# Get top biomarkers
top_biomarkers <- var_importance %>% arrange(desc(Importance)) %>% head(12)
print(top_biomarkers)

# Get predicted probabilities
pred_probs <- predict(rf_model, type = "prob")[, "1"]

# Get true outcomes
true_outcomes <- ifelse(features$death.cardiac == 1, TRUE, FALSE)

# Calculate AUC ROC
roc_curve <- roc(true_outcomes, pred_probs)
auc_base <- roc_curve$auc
print(paste("AUC ROC for Random Forest model:", auc_base))

# Calculate 95% CI for the AUC
ci_auc <- ci.auc(roc_curve)
print(paste("95% CI for AUC:", ci_auc[1], "-", ci_auc[3]))

#### LASSO ####-----------------------------------------------------------------
# Prepare data for LASSO
x <- model.matrix(death.cardiac ~ ., features)[, -1]
y <- features$death.cardiac

# Train LASSO model
lasso_model <- cv.glmnet(x, y, family = "binomial", alpha = 1)

# Get coefficients from the best model
lasso_coef <- coef(lasso_model, s = "lambda.min")
selected_biomarkers <- rownames(lasso_coef)[lasso_coef[,1] != 0]
print(selected_biomarkers)

# Predict probabilities using the selected features
pred_probs <- predict(lasso_model, newx = x, s = "lambda.min", type = "response")

# Calculate AUC ROC
roc_curve <- roc(y, as.vector(pred_probs))
auc_base <- roc_curve$auc
print(paste("AUC ROC for Random Forest model:", auc_base))

# Calculate 95% CI for the AUC
ci_auc <- ci.auc(roc_curve)
print(paste("95% CI for AUC:", ci_auc[1], "-", ci_auc[3]))

#### XGBoost ####---------------------------------------------------------------
new_colnames <- c(
  "Age", "Sex", "BMI", "Smoker", "Systolic blood pressure", "Diabetes", 
  "Prior stroke/TIA", "Heart failure", "Renal failure", "CAD", "Study center", 
  "Ethnicity", "AF type", "Beer drinker", "Red wine drinker", "White wine drinker", 
  "Liquor drinker", "AF episode duration", "Nr. of AF episodes", "Atrial flutter", 
  "Aspirin", "TCA", "Antiplatelet therapy", "Illicit drugs", "Prior CABG", 
  "Hyperthyroidism", "Hypothyroidism", "Cancer", "VTE", "Paternal AF", "Brother AF", 
  "Maternal AF", "Sister AF", "Familial hypertension", "Familial diabetes", 
  "Familial obesity", "Familial CAD", "Height", "Weight", "Heart rate", 
  "Diastolic blood pressure", "Rhythm at baseline", "Sleep apnea", "Hypertension", 
  "PAD", "Prior MI", "Systemic embolism", "Prior major bleeding"
)

# Assign the new column names to the dataframe
colnames(features) <- new_colnames

features$death.cardiac <- as.numeric(features$death.cardiac == 1)

# Train XGBoost model
set.seed(123)
xgb_model <- xgboost(data = as.matrix(features[, -c("death.cardiac")]), 
                     label = features$death.cardiac, 
                     objective = "binary:logistic", 
                     nrounds = 10)

# Get variable importance
var_importance_xgb <- xgb.importance(feature_names = colnames(features[, -c("death.cardiac")]), model = xgb_model)

# Get top biomarkers
top_biomarkers_xgb <- var_importance_xgb$Feature[1:12]
print(top_biomarkers_xgb)

# Get predicted probabilities
pred_probs_xgb <- as.numeric(predict(xgb_model, as.matrix(features[, -c("death.cardiac")])))

# Get true outcomes
true_outcomes <- ifelse(features$death.cardiac == 1, TRUE, FALSE)

# Calculate AUC ROC for XGBoost
roc_curve_xgb <- roc(true_outcomes, pred_probs_xgb)
auc_xgb_base <- roc_curve_xgb$auc
print(paste("AUC ROC for XGBoost model:", auc_xgb_base))

# Calculate 95% CI for the AUC
ci_auc_xgb <- ci.auc(roc_curve_xgb)
print(paste("95% CI for AUC:", ci_auc_xgb[1], "-", ci_auc_xgb[3]))

# Plot variable importance
plot(var_importance_xgb, col = "black", pch = 20, cex = 1.25, top = 20)
xgb.plot.importance(var_importance_xgb, top_n = 20, measure = "Gain", rel_to_first = TRUE, xlab = "Relative Importance")

ggplot(var_importance_xgb, aes(x = reorder(Feature, Gain), y = Gain)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Variable Importance from XGBoost Model",
       x = "Features",
       y = "Gain") +
  theme_minimal()

ggplot(var_importance_xgb, aes(x = Gain, y = reorder(Feature, Gain))) +
  geom_point(shape = 1, size = 2.5, color = "black") +
  labs(x = "Importance", size = 2, y = "", color = "black") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(color = "black"),  # Set x-axis ticks color to black
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.major.y = element_line(linetype = "dotted", color = "grey"))

#### Model with biomarkers #####################################################
#### Random forest ####---------------------------------------------------------
library(caret)
library(randomForest)
library(glmnet)
library(dplyr)
library(pROC)
library(xgboost)

# Exclude patients with missing variables
data <- dat %>% 
  select(ends_with("_log"), age.bl, pat.sex, bmi, current.smoker, rr.sys.liegend,
         prev.diabetes, prev.stroke.tia, prev.heart.failure, prev.niereninsuff,
         coronary.heart.disease, center, herkunft, vhf.typ.aktuell.bl, alkohol.bier, 
         alkohol.rotwein, alkohol.weisswein, alkohol.schnaps, vhf.episoden.dauer,
         vhf.episoden.anz, vorhofflattern, med.aspirin, med.tca.yn, med.antiplatelet.yn,
         drogen, prev.akb, prev.hyperthyreose, prev.hypothyreose, krk.malignom, krk.tvt,
         fam.vhf.vater, fam.vhf.bruder, fam.vhf.mutter, fam.vhf.schwester, fam.hypertonie,
         fam.diabetes, fam.uebergewicht, fam.khk, groesse, gewicht, heart.rate,
         rr.dia.liegend, ecg.rhythm.algo, prev.schlaf.apnoe, prev.hypertonie, prev.pavk,
         prev.mi, prev.sys.embolism, prev.major.bleed, death.cardiac)

data <- na.omit(data)

# Make feature and target outcome data frames"
features <- data %>% 
  select(age.bl, pat.sex, bmi, current.smoker, rr.sys.liegend,
         prev.diabetes, prev.stroke.tia, prev.heart.failure, prev.niereninsuff,
         coronary.heart.disease, center, herkunft, vhf.typ.aktuell.bl, alkohol.bier, 
         alkohol.rotwein, alkohol.weisswein, alkohol.schnaps, vhf.episoden.dauer,
         vhf.episoden.anz, vorhofflattern, med.aspirin, med.tca.yn, med.antiplatelet.yn,
         drogen, prev.akb, prev.hyperthyreose, prev.hypothyreose, krk.malignom, krk.tvt,
         fam.vhf.vater, fam.vhf.bruder, fam.vhf.mutter, fam.vhf.schwester, fam.hypertonie,
         fam.diabetes, fam.uebergewicht, fam.khk, groesse, gewicht, heart.rate,
         rr.dia.liegend, ecg.rhythm.algo, prev.schlaf.apnoe, prev.hypertonie, prev.pavk,
         prev.mi, prev.sys.embolism, prev.major.bleed, ends_with("_log"), death.cardiac)

features$death.cardiac <- as.factor(features$death.cardiac)

# Train random forest model
set.seed(123)

# Train a Random Forest model
rf_model <- randomForest(death.cardiac ~ ., data = features, importance = TRUE)

# Get variable importance
importance <- importance(rf_model)
var_importance <- data.frame(Variables = rownames(importance), Importance = importance[, 'MeanDecreaseGini'])

# Get top biomarkers
top_biomarkers <- var_importance %>% arrange(desc(Importance)) %>% head(12)
print(top_biomarkers)

# Get predicted probabilities
pred_probs <- predict(rf_model, type = "prob")[, "1"]

# Get true outcomes
true_outcomes <- ifelse(features$death.cardiac == 1, TRUE, FALSE)

# Calculate AUC ROC
roc_curve <- roc(true_outcomes, pred_probs)
auc_biomarkers <- roc_curve$auc
print(paste("AUC ROC for Random Forest model:", auc_biomarkers))

# Calculate 95% CI for the AUC
ci_auc <- ci.auc(roc_curve)
print(paste("95% CI for AUC:", ci_auc[1], "-", ci_auc[3]))

# Perform DeLong's test to compare the AUCs
delong_test <- roc.test(auc_base, auc_biomarkers)
print(delong_test)

#### LASSO ####
# Prepare data for LASSO
x <- model.matrix(death.cardiac ~ ., features)[, -1]
y <- features$death.cardiac

# Train LASSO model
lasso_model <- cv.glmnet(x, y, family = "binomial", alpha = 1)

# Get coefficients from the best model
lasso_coef <- coef(lasso_model, s = "lambda.min")
selected_biomarkers <- rownames(lasso_coef)[lasso_coef[,1] != 0]
print(selected_biomarkers)

# Predict probabilities using the selected features
pred_probs <- predict(lasso_model, newx = x, s = "lambda.min", type = "response")

# Calculate AUC ROC
roc_curve <- roc(true_outcomes, pred_probs)
auc_biomarkers <- roc_curve$auc
print(paste("AUC ROC for Random Forest model:", auc_biomarkers))

# Calculate 95% CI for the AUC
ci_auc <- ci.auc(roc_curve)
print(paste("95% CI for AUC:", ci_auc[1], "-", ci_auc[3]))

# Perform DeLong's test to compare the AUCs
delong_test <- roc.test(auc_base, auc_biomarkers)
print(delong_test)

#### XGBoost ####---------------------------------------------------------------
new_colnames <- c("Age", "Sex", "BMI", "Smoker", "Systolic blood pressure", "Diabetes", 
                  "Prior stroke/TIA", "Heart failure", "Renal failure", "CAD", "Study center", 
                  "Ethnicity", "AF type", "Beer drinker", "Red wine drinker", "White wine drinker", 
                  "Liquor drinker", "AF episode duration", "Nr. of AF episodes", "Atrial flutter", 
                  "Aspirin", "TCA", "Antiplatelet therapy", "Illicit drugs", "Prior CABG", 
                  "Hyperthyroidism", "Hypothyroidism", "Cancer", "VTE", "Paternal AF", "Brother AF", 
                  "Maternal AF", "Sister AF", "Familial hypertension", "Familial diabetes", 
                  "Familial obesity", "Familial CAD", "Height", "Weight", "Heart rate", 
                  "Diastolic blood pressure", "Rhythm at baseline", "Sleep apnea", "Hypertension", 
                  "PAD", "Prior MI", "Systemic embolism", "Prior major bleeding", "ANG2", "DDI2H", 
                  "CYSC", "ALAT", "GDF15", "CRPHS", "IGFBP7", "IL6", "PROBNPII", "OPN", "TNTHS", "eGFR")

# Assign the new column names to the dataframe
colnames(features) <- new_colnames

features$death.cardiac <- as.numeric(features$death.cardiac == 1)

# Train XGBoost model
set.seed(123)
xgb_model <- xgboost(data = as.matrix(features[, -c("death.cardiac")]), 
                     label = features$death.cardiac, 
                     objective = "binary:logistic", 
                     nrounds = 10)

# Get variable importance
var_importance_xgb <- xgb.importance(feature_names = colnames(features[, -c("death.cardiac")]), model = xgb_model)

# Get top biomarkers
top_biomarkers_xgb <- var_importance_xgb$Feature[1:12]
print(top_biomarkers_xgb)

# Get predicted probabilities
pred_probs_xgb <- as.numeric(predict(xgb_model, as.matrix(features[, -c("death.cardiac")])))

# Get true outcomes
true_outcomes <- ifelse(features$death.cardiac == 1, TRUE, FALSE)

# Calculate AUC ROC for XGBoost
roc_curve_xgb <- roc(true_outcomes, pred_probs_xgb)
auc_xgb_biomarkers <- roc_curve_xgb$auc
print(paste("AUC ROC for XGBoost model:", auc_xgb_biomarkers))

# Calculate 95% CI for the AUC
ci_auc_xgb <- ci.auc(roc_curve_xgb)
print(paste("95% CI for AUC:", ci_auc_xgb[1], "-", ci_auc_xgb[3]))

# Perform DeLong's test to compare the AUCs
delong_test <- roc.test(auc_xgb_base, auc_xgb_biomarkers)
print(delong_test)

# Plot variable importance
plot(var_importance_xgb, col = "black", pch = 20, cex = 1.25, top = 20)

ggplot(var_importance_xgb, aes(x = Gain, y = reorder(Feature, Gain))) +
  geom_point(shape = 1, size = 2.5, color = "black") +
  labs(x = "Importance", size = 2, y = "", color = "black") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(color = "black"),  # Set x-axis ticks color to black
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.major.y = element_line(linetype = "dotted", color = "grey"))


################################################################################
# All-cause death
################################################################################

#-------------------------------------------------------------------------------
# Check whether there are missing "time" variables 
num_missing_timeto1.death <- sum(is.na(dat$timeto1.death))
print(num_missing_timeto1.death)
#-------------------------------------------------------------------------------

# Split data into features (biomarkers) and target variable
library(caret)
library(randomForest)
library(glmnet)
library(dplyr)
library(pROC)
library(xgboost)

# Exclude patients with missing variables
data <- dat %>% 
  select(ends_with("_log"), age.bl, pat.sex, bmi, current.smoker, rr.sys.liegend,
         prev.diabetes, prev.stroke.tia, prev.heart.failure, prev.niereninsuff,
         coronary.heart.disease, center, herkunft, vhf.typ.aktuell.bl, alkohol.bier, 
         alkohol.rotwein, alkohol.weisswein, alkohol.schnaps, vhf.episoden.dauer,
         vhf.episoden.anz, vorhofflattern, med.aspirin, med.tca.yn, med.antiplatelet.yn,
         drogen, prev.akb, prev.hyperthyreose, prev.hypothyreose, krk.malignom, krk.tvt,
         fam.vhf.vater, fam.vhf.bruder, fam.vhf.mutter, fam.vhf.schwester, fam.hypertonie,
         fam.diabetes, fam.uebergewicht, fam.khk, groesse, gewicht, heart.rate,
         rr.dia.liegend, ecg.rhythm.algo, prev.schlaf.apnoe, prev.hypertonie, prev.pavk,
         prev.mi, prev.sys.embolism, prev.major.bleed, death.any)

data <- na.omit(data)

# Make feature and target outcome data frames"
features <- data %>% 
  select(age.bl, pat.sex, bmi, current.smoker, rr.sys.liegend,
         prev.diabetes, prev.stroke.tia, prev.heart.failure, prev.niereninsuff,
         coronary.heart.disease, center, herkunft, vhf.typ.aktuell.bl, alkohol.bier, 
         alkohol.rotwein, alkohol.weisswein, alkohol.schnaps, vhf.episoden.dauer,
         vhf.episoden.anz, vorhofflattern, med.aspirin, med.tca.yn, med.antiplatelet.yn,
         drogen, prev.akb, prev.hyperthyreose, prev.hypothyreose, krk.malignom, krk.tvt,
         fam.vhf.vater, fam.vhf.bruder, fam.vhf.mutter, fam.vhf.schwester, fam.hypertonie,
         fam.diabetes, fam.uebergewicht, fam.khk, groesse, gewicht, heart.rate,
         rr.dia.liegend, ecg.rhythm.algo, prev.schlaf.apnoe, prev.hypertonie, prev.pavk,
         prev.mi, prev.sys.embolism, prev.major.bleed, death.any)

features$death.any <- as.factor(features$death.any)

#### Base model ################################################################
#### Random forest ####---------------------------------------------------------
# Train random forest model
set.seed(123)

# Train a Random Forest model
rf_model <- randomForest(death.any ~ ., data = features, importance = TRUE)

# Get variable importance
importance <- importance(rf_model)
var_importance <- data.frame(Variables = rownames(importance), Importance = importance[, 'MeanDecreaseGini'])

# Get top biomarkers
top_biomarkers <- var_importance %>% arrange(desc(Importance)) %>% head(12)
print(top_biomarkers)

# Get predicted probabilities
pred_probs <- predict(rf_model, type = "prob")[, "1"]

# Get true outcomes
true_outcomes <- ifelse(features$death.any == 1, TRUE, FALSE)

# Calculate AUC ROC
roc_curve <- roc(true_outcomes, pred_probs)
auc_base <- roc_curve$auc
print(paste("AUC ROC for Random Forest model:", auc_base))

# Calculate 95% CI for the AUC
ci_auc <- ci.auc(roc_curve)
print(paste("95% CI for AUC:", ci_auc[1], "-", ci_auc[3]))

#### LASSO ####-----------------------------------------------------------------
# Prepare data for LASSO
x <- model.matrix(death.any ~ ., features)[, -1]
y <- features$death.any

# Train LASSO model
lasso_model <- cv.glmnet(x, y, family = "binomial", alpha = 1)

# Get coefficients from the best model
lasso_coef <- coef(lasso_model, s = "lambda.min")
selected_biomarkers <- rownames(lasso_coef)[lasso_coef[,1] != 0]
print(selected_biomarkers)

# Predict probabilities using the selected features
pred_probs <- predict(lasso_model, newx = x, s = "lambda.min", type = "response")

# Calculate AUC ROC
roc_curve <- roc(y, as.vector(pred_probs))
auc_base <- roc_curve$auc
print(paste("AUC ROC for Random Forest model:", auc_base))

# Calculate 95% CI for the AUC
ci_auc <- ci.auc(roc_curve)
print(paste("95% CI for AUC:", ci_auc[1], "-", ci_auc[3]))

#### XGBoost ####---------------------------------------------------------------
new_colnames <- c("Age", "Sex", "BMI", "Smoker", "Systolic blood pressure", "Diabetes", 
                  "Prior stroke/TIA", "Heart failure", "Renal failure", "CAD", "Study center", 
                  "Ethnicity", "AF type", "Beer drinker", "Red wine drinker", "White wine drinker", 
                  "Liquor drinker", "AF episode duration", "Nr. of AF episodes", "Atrial flutter", 
                  "Aspirin", "TCA", "Antiplatelet therapy", "Illicit drugs", "Prior CABG", 
                  "Hyperthyroidism", "Hypothyroidism", "Cancer", "VTE", "Paternal AF", "Brother AF", 
                  "Maternal AF", "Sister AF", "Familial hypertension", "Familial diabetes", 
                  "Familial obesity", "Familial CAD", "Height", "Weight", "Heart rate", 
                  "Diastolic blood pressure", "Rhythm at baseline", "Sleep apnea", "Hypertension", 
                  "PAD", "Prior MI", "Systemic embolism", "Prior major bleeding")

# Assign the new column names to the dataframe
colnames(features) <- new_colnames

features$death.any <- as.numeric(features$death.any == 1)

# Train XGBoost model
set.seed(123)
xgb_model <- xgboost(data = as.matrix(features[, -c("death.any")]), 
                     label = features$death.any, 
                     objective = "binary:logistic", 
                     nrounds = 10)

# Get variable importance
var_importance_xgb <- xgb.importance(feature_names = colnames(features[, -c("death.any")]), model = xgb_model)

# Get top biomarkers
top_biomarkers_xgb <- var_importance_xgb$Feature[1:12]
print(top_biomarkers_xgb)

# Get predicted probabilities
pred_probs_xgb <- as.numeric(predict(xgb_model, as.matrix(features[, -c("death.any")])))

# Get true outcomes
true_outcomes <- ifelse(features$death.any == 1, TRUE, FALSE)

# Calculate AUC ROC for XGBoost
roc_curve_xgb <- roc(true_outcomes, pred_probs_xgb)
auc_xgb_base <- roc_curve_xgb$auc
print(paste("AUC ROC for XGBoost model:", auc_xgb_base))

# Calculate 95% CI for the AUC
ci_auc_xgb <- ci.auc(roc_curve_xgb)
print(paste("95% CI for AUC:", ci_auc_xgb[1], "-", ci_auc_xgb[3]))

# Plot variable importance
plot(var_importance_xgb, col = "black", pch = 20, cex = 1.25, top = 20)
xgb.plot.importance(var_importance_xgb, top_n = 20, measure = "Gain", rel_to_first = TRUE, xlab = "Relative Importance")

ggplot(var_importance_xgb, aes(x = reorder(Feature, Gain), y = Gain)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Variable Importance from XGBoost Model",
       x = "Features",
       y = "Gain") +
  theme_minimal()

ggplot(var_importance_xgb, aes(x = Gain, y = reorder(Feature, Gain))) +
  geom_point(shape = 1, size = 2.5, color = "black") +
  labs(x = "Importance", size = 2, y = "", color = "black") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(color = "black"),  # Set x-axis ticks color to black
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.major.y = element_line(linetype = "dotted", color = "grey"))

#### Model with biomarkers #####################################################
#### Random forest ####---------------------------------------------------------
library(caret)
library(randomForest)
library(glmnet)
library(dplyr)
library(pROC)
library(xgboost)

# Exclude patients with missing variables
data <- dat %>% 
  select(ends_with("_log"), age.bl, pat.sex, bmi, current.smoker, rr.sys.liegend,
         prev.diabetes, prev.stroke.tia, prev.heart.failure, prev.niereninsuff,
         coronary.heart.disease, center, herkunft, vhf.typ.aktuell.bl, alkohol.bier, 
         alkohol.rotwein, alkohol.weisswein, alkohol.schnaps, vhf.episoden.dauer,
         vhf.episoden.anz, vorhofflattern, med.aspirin, med.tca.yn, med.antiplatelet.yn,
         drogen, prev.akb, prev.hyperthyreose, prev.hypothyreose, krk.malignom, krk.tvt,
         fam.vhf.vater, fam.vhf.bruder, fam.vhf.mutter, fam.vhf.schwester, fam.hypertonie,
         fam.diabetes, fam.uebergewicht, fam.khk, groesse, gewicht, heart.rate,
         rr.dia.liegend, ecg.rhythm.algo, prev.schlaf.apnoe, prev.hypertonie, prev.pavk,
         prev.mi, prev.sys.embolism, prev.major.bleed, death.any)

data <- na.omit(data)

# Make feature and target outcome data frames"
features <- data %>% 
  select(age.bl, pat.sex, bmi, current.smoker, rr.sys.liegend,
         prev.diabetes, prev.stroke.tia, prev.heart.failure, prev.niereninsuff,
         coronary.heart.disease, center, herkunft, vhf.typ.aktuell.bl, alkohol.bier, 
         alkohol.rotwein, alkohol.weisswein, alkohol.schnaps, vhf.episoden.dauer,
         vhf.episoden.anz, vorhofflattern, med.aspirin, med.tca.yn, med.antiplatelet.yn,
         drogen, prev.akb, prev.hyperthyreose, prev.hypothyreose, krk.malignom, krk.tvt,
         fam.vhf.vater, fam.vhf.bruder, fam.vhf.mutter, fam.vhf.schwester, fam.hypertonie,
         fam.diabetes, fam.uebergewicht, fam.khk, groesse, gewicht, heart.rate,
         rr.dia.liegend, ecg.rhythm.algo, prev.schlaf.apnoe, prev.hypertonie, prev.pavk,
         prev.mi, prev.sys.embolism, prev.major.bleed, ends_with("_log"), death.any)

features$death.any <- as.factor(features$death.any)

# Train random forest model
set.seed(123)

# Train a Random Forest model
rf_model <- randomForest(death.any ~ ., data = features, importance = TRUE)

# Get variable importance
importance <- importance(rf_model)
var_importance <- data.frame(Variables = rownames(importance), Importance = importance[, 'MeanDecreaseGini'])

# Get top biomarkers
top_biomarkers <- var_importance %>% arrange(desc(Importance)) %>% head(12)
print(top_biomarkers)

# Get predicted probabilities
pred_probs <- predict(rf_model, type = "prob")[, "1"]

# Get true outcomes
true_outcomes <- ifelse(features$death.any == 1, TRUE, FALSE)

# Calculate AUC ROC
roc_curve <- roc(true_outcomes, pred_probs)
auc_biomarkers <- roc_curve$auc
print(paste("AUC ROC for Random Forest model:", auc_biomarkers))

# Calculate 95% CI for the AUC
ci_auc <- ci.auc(roc_curve)
print(paste("95% CI for AUC:", ci_auc[1], "-", ci_auc[3]))

# Perform DeLong's test to compare the AUCs
delong_test <- roc.test(auc_base, auc_biomarkers)
print(delong_test)

#### LASSO ####
# Prepare data for LASSO
x <- model.matrix(death.any ~ ., features)[, -1]
y <- features$death.any

# Train LASSO model
lasso_model <- cv.glmnet(x, y, family = "binomial", alpha = 1)

# Get coefficients from the best model
lasso_coef <- coef(lasso_model, s = "lambda.min")
selected_biomarkers <- rownames(lasso_coef)[lasso_coef[,1] != 0]
print(selected_biomarkers)

# Predict probabilities using the selected features
pred_probs <- predict(lasso_model, newx = x, s = "lambda.min", type = "response")

# Calculate AUC ROC
roc_curve <- roc(true_outcomes, pred_probs)
auc_biomarkers <- roc_curve$auc
print(paste("AUC ROC for Random Forest model:", auc_biomarkers))

# Calculate 95% CI for the AUC
ci_auc <- ci.auc(roc_curve)
print(paste("95% CI for AUC:", ci_auc[1], "-", ci_auc[3]))

# Perform DeLong's test to compare the AUCs
delong_test <- roc.test(auc_base, auc_biomarkers)
print(delong_test)

#### XGBoost ####---------------------------------------------------------------
new_colnames <- c("Age", "Sex", "BMI", "Smoker", "Systolic blood pressure", "Diabetes", 
                  "Prior stroke/TIA", "Heart failure", "Renal failure", "CAD", "Study center", 
                  "Ethnicity", "AF type", "Beer drinker", "Red wine drinker", "White wine drinker", 
                  "Liquor drinker", "AF episode duration", "Nr. of AF episodes", "Atrial flutter", 
                  "Aspirin", "TCA", "Antiplatelet therapy", "Illicit drugs", "Prior CABG", 
                  "Hyperthyroidism", "Hypothyroidism", "Cancer", "VTE", "Paternal AF", "Brother AF", 
                  "Maternal AF", "Sister AF", "Familial hypertension", "Familial diabetes", 
                  "Familial obesity", "Familial CAD", "Height", "Weight", "Heart rate", 
                  "Diastolic blood pressure", "Rhythm at baseline", "Sleep apnea", "Hypertension", 
                  "PAD", "Prior MI", "Systemic embolism", "Prior major bleeding", "ANG2", "DDI2H", 
                  "CYSC", "ALAT", "GDF15", "CRPHS", "IGFBP7", "IL6", "PROBNPII", "OPN", "TNTHS", "eGFR")

# Assign the new column names to the dataframe
colnames(features) <- new_colnames

features$death.any <- as.numeric(features$death.any == 1)

# Train XGBoost model
set.seed(123)
xgb_model <- xgboost(data = as.matrix(features[, -c("death.any")]), 
                     label = features$death.any, 
                     objective = "binary:logistic", 
                     nrounds = 10)

# Get variable importance
var_importance_xgb <- xgb.importance(feature_names = colnames(features[, -c("death.any")]), model = xgb_model)

# Get top biomarkers
top_biomarkers_xgb <- var_importance_xgb$Feature[1:12]
print(top_biomarkers_xgb)

# Get predicted probabilities
pred_probs_xgb <- as.numeric(predict(xgb_model, as.matrix(features[, -c("death.any")])))

# Get true outcomes
true_outcomes <- ifelse(features$death.any == 1, TRUE, FALSE)

# Calculate AUC ROC for XGBoost
roc_curve_xgb <- roc(true_outcomes, pred_probs_xgb)
auc_xgb_biomarkers <- roc_curve_xgb$auc
print(paste("AUC ROC for XGBoost model:", auc_xgb_biomarkers))

# Calculate 95% CI for the AUC
ci_auc_xgb <- ci.auc(roc_curve_xgb)
print(paste("95% CI for AUC:", ci_auc_xgb[1], "-", ci_auc_xgb[3]))

# Perform DeLong's test to compare the AUCs
delong_test <- roc.test(auc_xgb_base, auc_xgb_biomarkers)
print(delong_test)

# Plot variable importance
plot(var_importance_xgb, col = "black", pch = 20, cex = 1.25, top = 20)

ggplot(var_importance_xgb, aes(x = Gain, y = reorder(Feature, Gain))) +
  geom_point(shape = 1, size = 2.5, color = "black") +
  labs(x = "Importance", size = 2, y = "", color = "black") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(color = "black"),  # Set x-axis ticks color to black
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.major.y = element_line(linetype = "dotted", color = "grey"))


################################################################################
# Composite major and clinically relevant non-major bleeding
################################################################################

#-------------------------------------------------------------------------------
# Check whether there are missing "time" variables 
num_missing_timeto1.anybleed <- sum(is.na(dat$timeto1.anybleed))
print(num_missing_timeto1.anybleed)
#-------------------------------------------------------------------------------

# Split data into features (biomarkers) and target variable
library(caret)
library(randomForest)
library(glmnet)
library(dplyr)
library(pROC)
library(xgboost)

# Exclude patients with missing variables
data <- dat %>% 
  select(ends_with("_log"), age.bl, pat.sex, bmi, current.smoker, rr.sys.liegend,
         prev.diabetes, prev.stroke.tia, prev.heart.failure, prev.niereninsuff,
         coronary.heart.disease, center, herkunft, vhf.typ.aktuell.bl, alkohol.bier, 
         alkohol.rotwein, alkohol.weisswein, alkohol.schnaps, vhf.episoden.dauer,
         vhf.episoden.anz, vorhofflattern, med.aspirin, med.tca.yn, med.antiplatelet.yn,
         drogen, prev.akb, prev.hyperthyreose, prev.hypothyreose, krk.malignom, krk.tvt,
         fam.vhf.vater, fam.vhf.bruder, fam.vhf.mutter, fam.vhf.schwester, fam.hypertonie,
         fam.diabetes, fam.uebergewicht, fam.khk, groesse, gewicht, heart.rate,
         rr.dia.liegend, ecg.rhythm.algo, prev.schlaf.apnoe, prev.hypertonie, prev.pavk,
         prev.mi, prev.sys.embolism, prev.major.bleed, bleed.any)

data <- na.omit(data)

# Make feature and target outcome data frames"
features <- data %>% 
  select(age.bl, pat.sex, bmi, current.smoker, rr.sys.liegend,
         prev.diabetes, prev.stroke.tia, prev.heart.failure, prev.niereninsuff,
         coronary.heart.disease, center, herkunft, vhf.typ.aktuell.bl, alkohol.bier, 
         alkohol.rotwein, alkohol.weisswein, alkohol.schnaps, vhf.episoden.dauer,
         vhf.episoden.anz, vorhofflattern, med.aspirin, med.tca.yn, med.antiplatelet.yn,
         drogen, prev.akb, prev.hyperthyreose, prev.hypothyreose, krk.malignom, krk.tvt,
         fam.vhf.vater, fam.vhf.bruder, fam.vhf.mutter, fam.vhf.schwester, fam.hypertonie,
         fam.diabetes, fam.uebergewicht, fam.khk, groesse, gewicht, heart.rate,
         rr.dia.liegend, ecg.rhythm.algo, prev.schlaf.apnoe, prev.hypertonie, prev.pavk,
         prev.mi, prev.sys.embolism, prev.major.bleed, bleed.any)

features$bleed.any <- as.factor(features$bleed.any)

#### Base model ################################################################
#### Random forest ####---------------------------------------------------------
# Train random forest model
set.seed(123)

# Train a Random Forest model
rf_model <- randomForest(bleed.any ~ ., data = features, importance = TRUE)

# Get variable importance
importance <- importance(rf_model)
var_importance <- data.frame(Variables = rownames(importance), Importance = importance[, 'MeanDecreaseGini'])

# Get top biomarkers
top_biomarkers <- var_importance %>% arrange(desc(Importance)) %>% head(12)
print(top_biomarkers)

# Get predicted probabilities
pred_probs <- predict(rf_model, type = "prob")[, "1"]

# Get true outcomes
true_outcomes <- ifelse(features$bleed.any == 1, TRUE, FALSE)

# Calculate AUC ROC
roc_curve <- roc(true_outcomes, pred_probs)
auc_base <- roc_curve$auc
print(paste("AUC ROC for Random Forest model:", auc_base))

# Calculate 95% CI for the AUC
ci_auc <- ci.auc(roc_curve)
print(paste("95% CI for AUC:", ci_auc[1], "-", ci_auc[3]))

#### LASSO ####-----------------------------------------------------------------
# Prepare data for LASSO
x <- model.matrix(bleed.any ~ ., features)[, -1]
y <- features$bleed.any

# Train LASSO model
lasso_model <- cv.glmnet(x, y, family = "binomial", alpha = 1)

# Get coefficients from the best model
lasso_coef <- coef(lasso_model, s = "lambda.min")
selected_biomarkers <- rownames(lasso_coef)[lasso_coef[,1] != 0]
print(selected_biomarkers)

# Predict probabilities using the selected features
pred_probs <- predict(lasso_model, newx = x, s = "lambda.min", type = "response")

# Calculate AUC ROC
roc_curve <- roc(y, as.vector(pred_probs))
auc_base <- roc_curve$auc
print(paste("AUC ROC for Random Forest model:", auc_base))

# Calculate 95% CI for the AUC
ci_auc <- ci.auc(roc_curve)
print(paste("95% CI for AUC:", ci_auc[1], "-", ci_auc[3]))

#### XGBoost ####---------------------------------------------------------------
new_colnames <- c("Age", "Sex", "BMI", "Smoker", "Systolic blood pressure", "Diabetes", 
                  "Prior stroke/TIA", "Heart failure", "Renal failure", "CAD", "Study center", 
                  "Ethnicity", "AF type", "Beer drinker", "Red wine drinker", "White wine drinker", 
                  "Liquor drinker", "AF episode duration", "Nr. of AF episodes", "Atrial flutter", 
                  "Aspirin", "TCA", "Antiplatelet therapy", "Illicit drugs", "Prior CABG", 
                  "Hyperthyroidism", "Hypothyroidism", "Cancer", "VTE", "Paternal AF", "Brother AF", 
                  "Maternal AF", "Sister AF", "Familial hypertension", "Familial diabetes", 
                  "Familial obesity", "Familial CAD", "Height", "Weight", "Heart rate", 
                  "Diastolic blood pressure", "Rhythm at baseline", "Sleep apnea", "Hypertension", 
                  "PAD", "Prior MI", "Systemic embolism", "Prior major bleeding")

# Assign the new column names to the dataframe
colnames(features) <- new_colnames

features$bleed.any <- as.numeric(features$bleed.any == 1)

# Train XGBoost model
set.seed(123)
xgb_model <- xgboost(data = as.matrix(features[, -c("bleed.any")]), 
                     label = features$bleed.any, 
                     objective = "binary:logistic", 
                     nrounds = 10)

# Get variable importance
var_importance_xgb <- xgb.importance(feature_names = colnames(features[, -c("bleed.any")]), model = xgb_model)

# Get top biomarkers
top_biomarkers_xgb <- var_importance_xgb$Feature[1:12]
print(top_biomarkers_xgb)

# Get predicted probabilities
pred_probs_xgb <- as.numeric(predict(xgb_model, as.matrix(features[, -c("bleed.any")])))

# Get true outcomes
true_outcomes <- ifelse(features$bleed.any == 1, TRUE, FALSE)

# Calculate AUC ROC for XGBoost
roc_curve_xgb <- roc(true_outcomes, pred_probs_xgb)
auc_xgb_base <- roc_curve_xgb$auc
print(paste("AUC ROC for XGBoost model:", auc_xgb_base))

# Calculate 95% CI for the AUC
ci_auc_xgb <- ci.auc(roc_curve_xgb)
print(paste("95% CI for AUC:", ci_auc_xgb[1], "-", ci_auc_xgb[3]))

# Plot variable importance
plot(var_importance_xgb, col = "black", pch = 20, cex = 1.25, top = 20)
xgb.plot.importance(var_importance_xgb, top_n = 20, measure = "Gain", rel_to_first = TRUE, xlab = "Relative Importance")

ggplot(var_importance_xgb, aes(x = reorder(Feature, Gain), y = Gain)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Variable Importance from XGBoost Model",
       x = "Features",
       y = "Gain") +
  theme_minimal()

ggplot(var_importance_xgb, aes(x = Gain, y = reorder(Feature, Gain))) +
  geom_point(shape = 1, size = 2.5, color = "black") +
  labs(x = "Importance", size = 2, y = "", color = "black") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(color = "black"),  # Set x-axis ticks color to black
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.major.y = element_line(linetype = "dotted", color = "grey"))

#### Model with biomarkers #####################################################
#### Random forest ####---------------------------------------------------------
library(caret)
library(randomForest)
library(glmnet)
library(dplyr)
library(pROC)
library(xgboost)

# Exclude patients with missing variables
data <- dat %>% 
  select(ends_with("_log"), age.bl, pat.sex, bmi, current.smoker, rr.sys.liegend,
         prev.diabetes, prev.stroke.tia, prev.heart.failure, prev.niereninsuff,
         coronary.heart.disease, center, herkunft, vhf.typ.aktuell.bl, alkohol.bier, 
         alkohol.rotwein, alkohol.weisswein, alkohol.schnaps, vhf.episoden.dauer,
         vhf.episoden.anz, vorhofflattern, med.aspirin, med.tca.yn, med.antiplatelet.yn,
         drogen, prev.akb, prev.hyperthyreose, prev.hypothyreose, krk.malignom, krk.tvt,
         fam.vhf.vater, fam.vhf.bruder, fam.vhf.mutter, fam.vhf.schwester, fam.hypertonie,
         fam.diabetes, fam.uebergewicht, fam.khk, groesse, gewicht, heart.rate,
         rr.dia.liegend, ecg.rhythm.algo, prev.schlaf.apnoe, prev.hypertonie, prev.pavk,
         prev.mi, prev.sys.embolism, prev.major.bleed, bleed.any)

data <- na.omit(data)

# Make feature and target outcome data frames"
features <- data %>% 
  select(age.bl, pat.sex, bmi, current.smoker, rr.sys.liegend,
         prev.diabetes, prev.stroke.tia, prev.heart.failure, prev.niereninsuff,
         coronary.heart.disease, center, herkunft, vhf.typ.aktuell.bl, alkohol.bier, 
         alkohol.rotwein, alkohol.weisswein, alkohol.schnaps, vhf.episoden.dauer,
         vhf.episoden.anz, vorhofflattern, med.aspirin, med.tca.yn, med.antiplatelet.yn,
         drogen, prev.akb, prev.hyperthyreose, prev.hypothyreose, krk.malignom, krk.tvt,
         fam.vhf.vater, fam.vhf.bruder, fam.vhf.mutter, fam.vhf.schwester, fam.hypertonie,
         fam.diabetes, fam.uebergewicht, fam.khk, groesse, gewicht, heart.rate,
         rr.dia.liegend, ecg.rhythm.algo, prev.schlaf.apnoe, prev.hypertonie, prev.pavk,
         prev.mi, prev.sys.embolism, prev.major.bleed, ends_with("_log"), bleed.any)

features$bleed.any <- as.factor(features$bleed.any)

# Train random forest model
set.seed(123)

# Train a Random Forest model
rf_model <- randomForest(bleed.any ~ ., data = features, importance = TRUE)

# Get variable importance
importance <- importance(rf_model)
var_importance <- data.frame(Variables = rownames(importance), Importance = importance[, 'MeanDecreaseGini'])

# Get top biomarkers
top_biomarkers <- var_importance %>% arrange(desc(Importance)) %>% head(12)
print(top_biomarkers)

# Get predicted probabilities
pred_probs <- predict(rf_model, type = "prob")[, "1"]

# Get true outcomes
true_outcomes <- ifelse(features$bleed.any == 1, TRUE, FALSE)

# Calculate AUC ROC
roc_curve <- roc(true_outcomes, pred_probs)
auc_biomarkers <- roc_curve$auc
print(paste("AUC ROC for Random Forest model:", auc_biomarkers))

# Calculate 95% CI for the AUC
ci_auc <- ci.auc(roc_curve)
print(paste("95% CI for AUC:", ci_auc[1], "-", ci_auc[3]))

# Perform DeLong's test to compare the AUCs
delong_test <- roc.test(auc_base, auc_biomarkers)
print(delong_test)

#### LASSO ####
# Prepare data for LASSO
x <- model.matrix(bleed.any ~ ., features)[, -1]
y <- features$bleed.any

# Train LASSO model
lasso_model <- cv.glmnet(x, y, family = "binomial", alpha = 1)

# Get coefficients from the best model
lasso_coef <- coef(lasso_model, s = "lambda.min")
selected_biomarkers <- rownames(lasso_coef)[lasso_coef[,1] != 0]
print(selected_biomarkers)

# Predict probabilities using the selected features
pred_probs <- predict(lasso_model, newx = x, s = "lambda.min", type = "response")

# Calculate AUC ROC
roc_curve <- roc(true_outcomes, pred_probs)
auc_biomarkers <- roc_curve$auc
print(paste("AUC ROC for Random Forest model:", auc_biomarkers))

# Calculate 95% CI for the AUC
ci_auc <- ci.auc(roc_curve)
print(paste("95% CI for AUC:", ci_auc[1], "-", ci_auc[3]))

# Perform DeLong's test to compare the AUCs
delong_test <- roc.test(auc_base, auc_biomarkers)
print(delong_test)

#### XGBoost ####---------------------------------------------------------------
new_colnames <- c("Age", "Sex", "BMI", "Smoker", "Systolic blood pressure", "Diabetes", 
                  "Prior stroke/TIA", "Heart failure", "Renal failure", "CAD", "Study center", 
                  "Ethnicity", "AF type", "Beer drinker", "Red wine drinker", "White wine drinker", 
                  "Liquor drinker", "AF episode duration", "Nr. of AF episodes", "Atrial flutter", 
                  "Aspirin", "TCA", "Antiplatelet therapy", "Illicit drugs", "Prior CABG", 
                  "Hyperthyroidism", "Hypothyroidism", "Cancer", "VTE", "Paternal AF", "Brother AF", 
                  "Maternal AF", "Sister AF", "Familial hypertension", "Familial diabetes", 
                  "Familial obesity", "Familial CAD", "Height", "Weight", "Heart rate", 
                  "Diastolic blood pressure", "Rhythm at baseline", "Sleep apnea", "Hypertension", 
                  "PAD", "Prior MI", "Systemic embolism", "Prior major bleeding", 
                  "ANG2", "DDI2H", "CYSC", "ALAT", "GDF15", "CRPHS", "IGFBP7", "IL6", 
                  "PROBNPII", "OPN", "TNTHS", "eGFR")

# Assign the new column names to the features data frame
colnames(features) <- new_colnames

features$bleed.any <- as.numeric(features$bleed.any == 1)

# Train XGBoost model
set.seed(123)
xgb_model <- xgboost(data = as.matrix(features[, -c("bleed.any")]), 
                     label = features$bleed.any, 
                     objective = "binary:logistic", 
                     nrounds = 10)

# Get variable importance
var_importance_xgb <- xgb.importance(feature_names = colnames(features[, -c("bleed.any")]), model = xgb_model)

# Get top biomarkers
top_biomarkers_xgb <- var_importance_xgb$Feature[1:12]
print(top_biomarkers_xgb)

# Get predicted probabilities
pred_probs_xgb <- as.numeric(predict(xgb_model, as.matrix(features[, -c("bleed.any")])))

# Get true outcomes
true_outcomes <- ifelse(features$bleed.any == 1, TRUE, FALSE)

# Calculate AUC ROC for XGBoost
roc_curve_xgb <- roc(true_outcomes, pred_probs_xgb)
auc_xgb_biomarkers <- roc_curve_xgb$auc
print(paste("AUC ROC for XGBoost model:", auc_xgb_biomarkers))

# Calculate 95% CI for the AUC
ci_auc_xgb <- ci.auc(roc_curve_xgb)
print(paste("95% CI for AUC:", ci_auc_xgb[1], "-", ci_auc_xgb[3]))

# Perform DeLong's test to compare the AUCs
delong_test <- roc.test(auc_xgb_base, auc_xgb_biomarkers)
print(delong_test)

# Plot variable importance
plot(var_importance_xgb, col = "black", pch = 20, cex = 1.25, top = 20)

ggplot(var_importance_xgb, aes(x = Gain, y = reorder(Feature, Gain))) +
  geom_point(shape = 1, size = 2.5, color = "black") +
  labs(x = "Importance", size = 2, y = "", color = "black") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(color = "black"),  # Set x-axis ticks color to black
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.major.y = element_line(linetype = "dotted", color = "grey"))


################################################################################
# Clinically relevant non-major bleeding
################################################################################

#-------------------------------------------------------------------------------
# Check whether there are missing "time" variables 
num_missing_timeto1.minor.bleed <- sum(is.na(dat$timeto1.minor.bleed))
print(num_missing_timeto1.minor.bleed)
#-------------------------------------------------------------------------------

# Split data into features (biomarkers) and target variable
library(caret)
library(randomForest)
library(glmnet)
library(dplyr)
library(pROC)
library(xgboost)

# Exclude patients with missing variables
data <- dat %>% 
  select(ends_with("_log"), age.bl, pat.sex, bmi, current.smoker, rr.sys.liegend,
         prev.diabetes, prev.stroke.tia, prev.heart.failure, prev.niereninsuff,
         coronary.heart.disease, center, herkunft, vhf.typ.aktuell.bl, alkohol.bier, 
         alkohol.rotwein, alkohol.weisswein, alkohol.schnaps, vhf.episoden.dauer,
         vhf.episoden.anz, vorhofflattern, med.aspirin, med.tca.yn, med.antiplatelet.yn,
         drogen, prev.akb, prev.hyperthyreose, prev.hypothyreose, krk.malignom, krk.tvt,
         fam.vhf.vater, fam.vhf.bruder, fam.vhf.mutter, fam.vhf.schwester, fam.hypertonie,
         fam.diabetes, fam.uebergewicht, fam.khk, groesse, gewicht, heart.rate,
         rr.dia.liegend, ecg.rhythm.algo, prev.schlaf.apnoe, prev.hypertonie, prev.pavk,
         prev.mi, prev.sys.embolism, prev.major.bleed, minor.bleed)

data <- na.omit(data)

# Make feature and target outcome data frames"
features <- data %>% 
  select(age.bl, pat.sex, bmi, current.smoker, rr.sys.liegend,
         prev.diabetes, prev.stroke.tia, prev.heart.failure, prev.niereninsuff,
         coronary.heart.disease, center, herkunft, vhf.typ.aktuell.bl, alkohol.bier, 
         alkohol.rotwein, alkohol.weisswein, alkohol.schnaps, vhf.episoden.dauer,
         vhf.episoden.anz, vorhofflattern, med.aspirin, med.tca.yn, med.antiplatelet.yn,
         drogen, prev.akb, prev.hyperthyreose, prev.hypothyreose, krk.malignom, krk.tvt,
         fam.vhf.vater, fam.vhf.bruder, fam.vhf.mutter, fam.vhf.schwester, fam.hypertonie,
         fam.diabetes, fam.uebergewicht, fam.khk, groesse, gewicht, heart.rate,
         rr.dia.liegend, ecg.rhythm.algo, prev.schlaf.apnoe, prev.hypertonie, prev.pavk,
         prev.mi, prev.sys.embolism, prev.major.bleed, minor.bleed)

features$minor.bleed <- as.factor(features$minor.bleed)

#### Base model ################################################################
#### Random forest ####---------------------------------------------------------
# Train random forest model
set.seed(123)

# Train a Random Forest model
rf_model <- randomForest(minor.bleed ~ ., data = features, importance = TRUE)

# Get variable importance
importance <- importance(rf_model)
var_importance <- data.frame(Variables = rownames(importance), Importance = importance[, 'MeanDecreaseGini'])

# Get top biomarkers
top_biomarkers <- var_importance %>% arrange(desc(Importance)) %>% head(12)
print(top_biomarkers)

# Get predicted probabilities
pred_probs <- predict(rf_model, type = "prob")[, "1"]

# Get true outcomes
true_outcomes <- ifelse(features$minor.bleed == 1, TRUE, FALSE)

# Calculate AUC ROC
roc_curve <- roc(true_outcomes, pred_probs)
auc_base <- roc_curve$auc
print(paste("AUC ROC for Random Forest model:", auc_base))

# Calculate 95% CI for the AUC
ci_auc <- ci.auc(roc_curve)
print(paste("95% CI for AUC:", ci_auc[1], "-", ci_auc[3]))

#### LASSO ####-----------------------------------------------------------------
# Prepare data for LASSO
x <- model.matrix(minor.bleed ~ ., features)[, -1]
y <- features$minor.bleed

# Train LASSO model
lasso_model <- cv.glmnet(x, y, family = "binomial", alpha = 1)

# Get coefficients from the best model
lasso_coef <- coef(lasso_model, s = "lambda.min")
selected_biomarkers <- rownames(lasso_coef)[lasso_coef[,1] != 0]
print(selected_biomarkers)

# Predict probabilities using the selected features
pred_probs <- predict(lasso_model, newx = x, s = "lambda.min", type = "response")

# Calculate AUC ROC
roc_curve <- roc(y, as.vector(pred_probs))
auc_base <- roc_curve$auc
print(paste("AUC ROC for Random Forest model:", auc_base))

# Calculate 95% CI for the AUC
ci_auc <- ci.auc(roc_curve)
print(paste("95% CI for AUC:", ci_auc[1], "-", ci_auc[3]))

#### XGBoost ####---------------------------------------------------------------
new_colnames <- c("Age", "Sex", "BMI", "Smoker", "Systolic blood pressure", "Diabetes", 
                  "Prior stroke/TIA", "Heart failure", "Renal failure", "CAD", "Study center", 
                  "Ethnicity", "AF type", "Beer drinker", "Red wine drinker", "White wine drinker", 
                  "Liquor drinker", "AF episode duration", "Nr. of AF episodes", "Atrial flutter", 
                  "Aspirin", "TCA", "Antiplatelet therapy", "Illicit drugs", "Prior CABG", 
                  "Hyperthyroidism", "Hypothyroidism", "Cancer", "VTE", "Paternal AF", "Brother AF", 
                  "Maternal AF", "Sister AF", "Familial hypertension", "Familial diabetes", 
                  "Familial obesity", "Familial CAD", "Height", "Weight", "Heart rate", 
                  "Diastolic blood pressure", "Rhythm at baseline", "Sleep apnea", "Hypertension", 
                  "PAD", "Prior MI", "Systemic embolism", "Prior major bleeding")

# Assign the new column names to the dataframe
colnames(features) <- new_colnames

features$minor.bleed <- as.numeric(features$minor.bleed == 1)

# Train XGBoost model
set.seed(123)
xgb_model <- xgboost(data = as.matrix(features[, -c("minor.bleed")]), 
                     label = features$minor.bleed, 
                     objective = "binary:logistic", 
                     nrounds = 10)

# Get variable importance
var_importance_xgb <- xgb.importance(feature_names = colnames(features[, -c("minor.bleed")]), model = xgb_model)

# Get top biomarkers
top_biomarkers_xgb <- var_importance_xgb$Feature[1:12]
print(top_biomarkers_xgb)

# Get predicted probabilities
pred_probs_xgb <- as.numeric(predict(xgb_model, as.matrix(features[, -c("minor.bleed")])))

# Get true outcomes
true_outcomes <- ifelse(features$minor.bleed == 1, TRUE, FALSE)

# Calculate AUC ROC for XGBoost
roc_curve_xgb <- roc(true_outcomes, pred_probs_xgb)
auc_xgb_base <- roc_curve_xgb$auc
print(paste("AUC ROC for XGBoost model:", auc_xgb_base))

# Calculate 95% CI for the AUC
ci_auc_xgb <- ci.auc(roc_curve_xgb)
print(paste("95% CI for AUC:", ci_auc_xgb[1], "-", ci_auc_xgb[3]))

# Plot variable importance
plot(var_importance_xgb, col = "black", pch = 20, cex = 1.25, top = 20)
xgb.plot.importance(var_importance_xgb, top_n = 20, measure = "Gain", rel_to_first = TRUE, xlab = "Relative Importance")

ggplot(var_importance_xgb, aes(x = reorder(Feature, Gain), y = Gain)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Variable Importance from XGBoost Model",
       x = "Features",
       y = "Gain") +
  theme_minimal()

ggplot(var_importance_xgb, aes(x = Gain, y = reorder(Feature, Gain))) +
  geom_point(shape = 1, size = 2.5, color = "black") +
  labs(x = "Importance", size = 2, y = "", color = "black") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(color = "black"),  # Set x-axis ticks color to black
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.major.y = element_line(linetype = "dotted", color = "grey"))

#### Model with biomarkers #####################################################
#### Random forest ####---------------------------------------------------------
library(caret)
library(randomForest)
library(glmnet)
library(dplyr)
library(pROC)
library(xgboost)

# Exclude patients with missing variables
data <- dat %>% 
  select(ends_with("_log"), age.bl, pat.sex, bmi, current.smoker, rr.sys.liegend,
         prev.diabetes, prev.stroke.tia, prev.heart.failure, prev.niereninsuff,
         coronary.heart.disease, center, herkunft, vhf.typ.aktuell.bl, alkohol.bier, 
         alkohol.rotwein, alkohol.weisswein, alkohol.schnaps, vhf.episoden.dauer,
         vhf.episoden.anz, vorhofflattern, med.aspirin, med.tca.yn, med.antiplatelet.yn,
         drogen, prev.akb, prev.hyperthyreose, prev.hypothyreose, krk.malignom, krk.tvt,
         fam.vhf.vater, fam.vhf.bruder, fam.vhf.mutter, fam.vhf.schwester, fam.hypertonie,
         fam.diabetes, fam.uebergewicht, fam.khk, groesse, gewicht, heart.rate,
         rr.dia.liegend, ecg.rhythm.algo, prev.schlaf.apnoe, prev.hypertonie, prev.pavk,
         prev.mi, prev.sys.embolism, prev.major.bleed, minor.bleed)

data <- na.omit(data)

# Make feature and target outcome data frames"
features <- data %>% 
  select(age.bl, pat.sex, bmi, current.smoker, rr.sys.liegend,
         prev.diabetes, prev.stroke.tia, prev.heart.failure, prev.niereninsuff,
         coronary.heart.disease, center, herkunft, vhf.typ.aktuell.bl, alkohol.bier, 
         alkohol.rotwein, alkohol.weisswein, alkohol.schnaps, vhf.episoden.dauer,
         vhf.episoden.anz, vorhofflattern, med.aspirin, med.tca.yn, med.antiplatelet.yn,
         drogen, prev.akb, prev.hyperthyreose, prev.hypothyreose, krk.malignom, krk.tvt,
         fam.vhf.vater, fam.vhf.bruder, fam.vhf.mutter, fam.vhf.schwester, fam.hypertonie,
         fam.diabetes, fam.uebergewicht, fam.khk, groesse, gewicht, heart.rate,
         rr.dia.liegend, ecg.rhythm.algo, prev.schlaf.apnoe, prev.hypertonie, prev.pavk,
         prev.mi, prev.sys.embolism, prev.major.bleed, ends_with("_log"), minor.bleed)

features$minor.bleed <- as.factor(features$minor.bleed)

# Train random forest model
set.seed(123)

# Train a Random Forest model
rf_model <- randomForest(minor.bleed ~ ., data = features, importance = TRUE)

# Get variable importance
importance <- importance(rf_model)
var_importance <- data.frame(Variables = rownames(importance), Importance = importance[, 'MeanDecreaseGini'])

# Get top biomarkers
top_biomarkers <- var_importance %>% arrange(desc(Importance)) %>% head(12)
print(top_biomarkers)

# Get predicted probabilities
pred_probs <- predict(rf_model, type = "prob")[, "1"]

# Get true outcomes
true_outcomes <- ifelse(features$minor.bleed == 1, TRUE, FALSE)

# Calculate AUC ROC
roc_curve <- roc(true_outcomes, pred_probs)
auc_biomarkers <- roc_curve$auc
print(paste("AUC ROC for Random Forest model:", auc_biomarkers))

# Calculate 95% CI for the AUC
ci_auc <- ci.auc(roc_curve)
print(paste("95% CI for AUC:", ci_auc[1], "-", ci_auc[3]))

# Perform DeLong's test to compare the AUCs
delong_test <- roc.test(auc_base, auc_biomarkers)
print(delong_test)

#### LASSO ####
# Prepare data for LASSO
x <- model.matrix(minor.bleed ~ ., features)[, -1]
y <- features$minor.bleed

# Train LASSO model
lasso_model <- cv.glmnet(x, y, family = "binomial", alpha = 1)

# Get coefficients from the best model
lasso_coef <- coef(lasso_model, s = "lambda.min")
selected_biomarkers <- rownames(lasso_coef)[lasso_coef[,1] != 0]
print(selected_biomarkers)

# Predict probabilities using the selected features
pred_probs <- predict(lasso_model, newx = x, s = "lambda.min", type = "response")

# Calculate AUC ROC
roc_curve <- roc(true_outcomes, pred_probs)
auc_biomarkers <- roc_curve$auc
print(paste("AUC ROC for Random Forest model:", auc_biomarkers))

# Calculate 95% CI for the AUC
ci_auc <- ci.auc(roc_curve)
print(paste("95% CI for AUC:", ci_auc[1], "-", ci_auc[3]))

# Perform DeLong's test to compare the AUCs
delong_test <- roc.test(auc_base, auc_biomarkers)
print(delong_test)

#### XGBoost ####---------------------------------------------------------------
new_colnames <- c("Age", "Sex", "BMI", "Smoker", "Systolic blood pressure", "Diabetes", 
                  "Prior stroke/TIA", "Heart failure", "Renal failure", "CAD", "Study center", 
                  "Ethnicity", "AF type", "Beer drinker", "Red wine drinker", "White wine drinker", 
                  "Liquor drinker", "AF episode duration", "Nr. of AF episodes", "Atrial flutter", 
                  "Aspirin", "TCA", "Antiplatelet therapy", "Illicit drugs", "Prior CABG", 
                  "Hyperthyroidism", "Hypothyroidism", "Cancer", "VTE", "Paternal AF", "Brother AF", 
                  "Maternal AF", "Sister AF", "Familial hypertension", "Familial diabetes", 
                  "Familial obesity", "Familial CAD", "Height", "Weight", "Heart rate", 
                  "Diastolic blood pressure", "Rhythm at baseline", "Sleep apnea", "Hypertension", 
                  "PAD", "Prior MI", "Systemic embolism", "Prior major bleeding", 
                  "ANG2", "DDI2H", "CYSC", "ALAT", "GDF15", "CRPHS", "IGFBP7", "IL6", 
                  "PROBNPII", "OPN", "TNTHS", "eGFR")

# Assign the new column names to the features data frame
colnames(features) <- new_colnames

features$minor.bleed <- as.numeric(features$minor.bleed == 1)

# Train XGBoost model
set.seed(123)
xgb_model <- xgboost(data = as.matrix(features[, -c("minor.bleed")]), 
                     label = features$minor.bleed, 
                     objective = "binary:logistic", 
                     nrounds = 10)

# Get variable importance
var_importance_xgb <- xgb.importance(feature_names = colnames(features[, -c("minor.bleed")]), model = xgb_model)

# Get top biomarkers
top_biomarkers_xgb <- var_importance_xgb$Feature[1:12]
print(top_biomarkers_xgb)

# Get predicted probabilities
pred_probs_xgb <- as.numeric(predict(xgb_model, as.matrix(features[, -c("minor.bleed")])))

# Get true outcomes
true_outcomes <- ifelse(features$minor.bleed == 1, TRUE, FALSE)

# Calculate AUC ROC for XGBoost
roc_curve_xgb <- roc(true_outcomes, pred_probs_xgb)
auc_xgb_biomarkers <- roc_curve_xgb$auc
print(paste("AUC ROC for XGBoost model:", auc_xgb_biomarkers))

# Calculate 95% CI for the AUC
ci_auc_xgb <- ci.auc(roc_curve_xgb)
print(paste("95% CI for AUC:", ci_auc_xgb[1], "-", ci_auc_xgb[3]))

# Perform DeLong's test to compare the AUCs
delong_test <- roc.test(auc_xgb_base, auc_xgb_biomarkers)
print(delong_test)

# Plot variable importance
plot(var_importance_xgb, col = "black", pch = 20, cex = 1.25, top = 20)

ggplot(var_importance_xgb, aes(x = Gain, y = reorder(Feature, Gain))) +
  geom_point(shape = 1, size = 2.5, color = "black") +
  labs(x = "Importance", size = 2, y = "", color = "black") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(color = "black"),  # Set x-axis ticks color to black
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.major.y = element_line(linetype = "dotted", color = "grey"))

################################# END ##########################################
