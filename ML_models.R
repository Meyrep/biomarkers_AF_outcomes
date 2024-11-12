#--ML MODELS--------------------------------------------------------------------
# Author: Pascal B. Meyre
# Date: 05/18/24
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
dat$center <- ifelse(dat$center == "Universit\xe4tsspital Basel", 1,
                     ifelse(dat$center == "Inselspital, Universit\xe4tsspital Bern", 2,
                            ifelse(dat$center == "Triemli Spital Z\xfcrich", 3, 
                                   ifelse(dat$center == "Kantonsspital Baden", 4,
                                          ifelse(dat$center == "Kantonsspital St. Gallen", 5,
                                                 ifelse(dat$center == "Cardiocentro Lugano", 6,
                                                        ifelse(dat$center == "EOC Lugano", 7,
                                                               ifelse(dat$center == "Luzerner Kantonsspital", 8,
                                                                      ifelse(dat$center == "H\xf4pital Cantonal Fribourg", 9,
                                                                             ifelse(dat$center == "HUG", 10,
                                                                                    ifelse(dat$center == "B\xfcrgerspital Solothurn", 11,
                                                                                           ifelse(dat$center == "CHUV", 12,
                                                                                                  ifelse(dat$center == "Kantonsspital Luzern", 8,
                                                                                                         ifelse(dat$center == "Regionalspital Lugano", 7,
                                                                                                                ifelse(dat$center == "EOC Bellinzona", 13,
                                                                                                                       ifelse(dat$center == "Regionalspital Bellinzona", 13,
                                                                                                                              ifelse(dat$center == "HUG Gen\xe8ve", 11,
                                                                                                                                     ifelse(dat$center == "Hirslanden Klinik St. Anna, Luzern", 14,
                                                                                                                                            ifelse(dat$center == "CHUV Lausanne", 12,
                                                                                                                                                   ifelse(dat$center == "Universit\xe4tsspital Z\xfcrich", 15, NA))))))))))))))))))))

# Herkunft
dat$herkunft <- ifelse(dat$herkunft == "Mitteleuropa", 1,
                       ifelse(dat$herkunft == "S\xfcdeuropa", 2,
                              ifelse(dat$herkunft == "Nordeuropa", 3,
                                     ifelse(dat$herkunft == "Osteuropa", 4,
                                            ifelse(dat$herkunft == "andere", 5,
                                                   ifelse(dat$herkunft == "Mittel-/S\xfcdamerika", 6,
                                                          ifelse(dat$herkunft == "Nordamerika/weiss", 7,NA)))))))

# AF type
dat$vhf.typ.aktuell.bl <- ifelse(dat$vhf.typ.aktuell.bl == "paroxysmal", 1,
                                 ifelse(dat$vhf.typ.aktuell.bl == "persistierend (>7 Tage, EKV)", 2,
                                        ifelse(dat$vhf.typ.aktuell.bl == "permanent", 3, NA)))

# Alcohol beer
dat$alkohol.bier <- ifelse(dat$alkohol.bier == "6+ pro Tag", 1,
                           ifelse(dat$alkohol.bier == "4-5 pro Tag", 2,
                                  ifelse(dat$alkohol.bier == "2-3 pro Tag", 3,
                                         ifelse(dat$alkohol.bier == "1 pro Tag", 4,
                                                ifelse(dat$alkohol.bier == "5-6 pro Woche", 5,
                                                       ifelse(dat$alkohol.bier == "2-4 pro Woche", 6,
                                                              ifelse(dat$alkohol.bier == "1 pro Woche", 7,
                                                                     ifelse(dat$alkohol.bier == "1-3 pro Monat", 8,
                                                                            ifelse(dat$alkohol.bier == "nie oder weniger als 1 Monat", 9, NA)))))))))

# Alcohol red wine
dat$alkohol.rotwein <- ifelse(dat$alkohol.rotwein == "6+ pro Tag", 1,
                              ifelse(dat$alkohol.rotwein == "4-5 pro Tag", 2,
                                     ifelse(dat$alkohol.rotwein == "2-3 pro Tag", 3,
                                            ifelse(dat$alkohol.rotwein == "1 pro Tag", 4,
                                                   ifelse(dat$alkohol.rotwein == "5-6 pro Woche", 5,
                                                          ifelse(dat$alkohol.rotwein == "2-4 pro Woche", 6,
                                                                 ifelse(dat$alkohol.rotwein == "1 pro Woche", 7,
                                                                        ifelse(dat$alkohol.rotwein == "1-3 pro Monat", 8,
                                                                               ifelse(dat$alkohol.rotwein == "nie oder weniger als 1 Monat", 9, NA)))))))))

# Alcohol white wine
dat$alkohol.weisswein <- ifelse(dat$alkohol.weisswein == "6+ pro Tag", 1,
                                ifelse(dat$alkohol.weisswein == "4-5 pro Tag", 2,
                                       ifelse(dat$alkohol.weisswein == "2-3 pro Tag", 3,
                                              ifelse(dat$alkohol.weisswein == "1 pro Tag", 4,
                                                     ifelse(dat$alkohol.weisswein == "5-6 pro Woche", 5,
                                                            ifelse(dat$alkohol.weisswein == "2-4 pro Woche", 6,
                                                                   ifelse(dat$alkohol.weisswein == "1 pro Woche", 7,
                                                                          ifelse(dat$alkohol.weisswein == "1-3 pro Monat", 8,
                                                                                 ifelse(dat$alkohol.weisswein == "nie oder weniger als 1 Monat", 9, NA)))))))))

# Alcohol liquor
dat$alkohol.schnaps <- ifelse(dat$alkohol.schnaps == "6+ pro Tag", 1,
                              ifelse(dat$alkohol.schnaps == "2-3 pro Tag", 2,
                                     ifelse(dat$alkohol.schnaps == "1 pro Tag", 3,
                                            ifelse(dat$alkohol.schnaps == "5-6 pro Woche", 4,
                                                   ifelse(dat$alkohol.schnaps == "2-4 pro Woche", 5,
                                                          ifelse(dat$alkohol.schnaps == "1 pro Woche", 6,
                                                                 ifelse(dat$alkohol.schnaps == "1-3 pro Monat", 7,
                                                                        ifelse(dat$alkohol.schnaps == "nie oder weniger als 1 Monat", 8, NA))))))))

# AF duration
dat$vhf.episoden.dauer <- ifelse(dat$vhf.episoden.dauer == "dauerhaft", 1,
                                 ifelse(dat$vhf.episoden.dauer == "Tage", 2,
                                        ifelse(dat$vhf.episoden.dauer == "Stunden", 3,
                                               ifelse(dat$vhf.episoden.dauer == "Minuten", 4,
                                                      ifelse(dat$vhf.episoden.dauer == "keine mehr", 5, NA)))))

# AF nr. of episodes
dat$vhf.episoden.anz <- ifelse(dat$vhf.episoden.anz == "dauerhaft", 1,
                               ifelse(dat$vhf.episoden.anz == ">=1x pro Woche", 2,
                                      ifelse(dat$vhf.episoden.anz == "=1x pro Woche", 3,
                                             ifelse(dat$vhf.episoden.anz == "<1x pro Woche aber >1x pro Monat", 4,
                                                    ifelse(dat$vhf.episoden.anz == "<1x pro Monat", 5,
                                                           ifelse(dat$vhf.episoden.anz == "keine mehr", 6, NA))))))

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
dat$ecg.rhythm.algo <- ifelse(dat$ecg.rhythm.algo == "AF", 1,
                              ifelse(dat$ecg.rhythm.algo == "Fibrillation", 1,
                                     ifelse(dat$ecg.rhythm.algo == "Flutter", 2, 
                                            ifelse(dat$ecg.rhythm.algo == "Sinus", 3,
                                                   ifelse(dat$ecg.rhythm.algo == "Other", 4, 
                                                          ifelse(dat$ecg.rhythm.algo == "unclear", 4, NA))))))

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
colnames(features)[1] <- "Age"
colnames(features)[2] <- "Sex"
colnames(features)[3] <- "BMI"
colnames(features)[4] <- "Smoker"
colnames(features)[5] <- "Systolic blood pressure"
colnames(features)[6] <- "Diabetes"
colnames(features)[7] <- "Prior stroke/TIA"
colnames(features)[8] <- "Heart failure"
colnames(features)[9] <- "Renal failure"
colnames(features)[10] <- "CAD"
colnames(features)[11] <- "Study center"
colnames(features)[12] <- "Ethnicity"
colnames(features)[13] <- "AF type"
colnames(features)[14] <- "Beer drinker"
colnames(features)[15] <- "Red wine drinker"
colnames(features)[16] <- "White wine drinker"
colnames(features)[17] <- "Liquor drinker"
colnames(features)[18] <- "AF episode duration"
colnames(features)[19] <- "Nr. of AF episodes"
colnames(features)[20] <- "Atrial flutter"
colnames(features)[21] <- "Aspirin"
colnames(features)[22] <- "TCA"
colnames(features)[23] <- "Antiplatelet therapy"
colnames(features)[24] <- "Illicit drugs"
colnames(features)[25] <- "Prior CABG"
colnames(features)[26] <- "Hyperthyroidism"
colnames(features)[27] <- "Hypothyroidism"
colnames(features)[28] <- "Cancer"
colnames(features)[29] <- "VTE"
colnames(features)[30] <- "Paternal AF"
colnames(features)[31] <- "Brother AF"
colnames(features)[32] <- "Maternal AF"
colnames(features)[33] <- "Sister AF"
colnames(features)[34] <- "Familial hypertension"
colnames(features)[35] <- "Familial diabetes"
colnames(features)[36] <- "Familial obesity"
colnames(features)[37] <- "Familial CAD"
colnames(features)[38] <- "Height"
colnames(features)[39] <- "Weight"
colnames(features)[40] <- "Heart rate"
colnames(features)[41] <- "Diastolic blood pressure"
colnames(features)[42] <- "Rhythm at baseline"
colnames(features)[43] <- "Sleep apnea"
colnames(features)[44] <- "Hypertension"
colnames(features)[45] <- "PAD"
colnames(features)[46] <- "Prior MI"
colnames(features)[47] <- "Systemic embolism"
colnames(features)[48] <- "Prior major bleeding"

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

colnames(features)[1] <- "Age"
colnames(features)[2] <- "Sex"
colnames(features)[3] <- "BMI"
colnames(features)[4] <- "Smoker"
colnames(features)[5] <- "SBP"
colnames(features)[6] <- "Diabetes"
colnames(features)[7] <- "Prior.stroke.TIA"
colnames(features)[8] <- "Heart.failure"
colnames(features)[9] <- "Renal.failure"
colnames(features)[10] <- "CAD"
colnames(features)[11] <- "Study.center"
colnames(features)[12] <- "Ethnicity"
colnames(features)[13] <- "AF.type"
colnames(features)[14] <- "Beer.drinker"
colnames(features)[15] <- "Red.wine.drinker"
colnames(features)[16] <- "White.wine.drinker"
colnames(features)[17] <- "Liquor.drinker"
colnames(features)[18] <- "AF.episode.duration"
colnames(features)[19] <- "Nr.of.AF.episodes"
colnames(features)[20] <- "Atrial.flutter"
colnames(features)[21] <- "Aspirin"
colnames(features)[22] <- "TCA"
colnames(features)[23] <- "Antiplatelet.therapy"
colnames(features)[24] <- "Illicit.drugs"
colnames(features)[25] <- "Prior.CABG"
colnames(features)[26] <- "Hyperthyroidism"
colnames(features)[27] <- "Hypothyroidism"
colnames(features)[28] <- "Cancer"
colnames(features)[29] <- "VTE"
colnames(features)[30] <- "Paternal.AF"
colnames(features)[31] <- "Brother.AF"
colnames(features)[32] <- "Maternal.AF"
colnames(features)[33] <- "Sister.AF"
colnames(features)[34] <- "Familial.hypertension"
colnames(features)[35] <- "Familial.diabetes"
colnames(features)[36] <- "Familial.obesity"
colnames(features)[37] <- "Familial.CAD"
colnames(features)[38] <- "Height"
colnames(features)[39] <- "Weight"
colnames(features)[40] <- "Heart.rate"
colnames(features)[41] <- "DBP"
colnames(features)[42] <- "Rhythm.at.baseline"
colnames(features)[43] <- "Sleep.apnea"
colnames(features)[44] <- "Hypertension"
colnames(features)[45] <- "PAD"
colnames(features)[46] <- "Prior.MI"
colnames(features)[47] <- "Systemic.embolism"
colnames(features)[48] <- "Prior.major.bleeding"
colnames(features)[49] <- "ANG.2"
colnames(features)[50] <- "d.dimer"
colnames(features)[51] <- "cystatin.C"
colnames(features)[52] <- "ALAT"
colnames(features)[53] <- "GDF.15"
colnames(features)[54] <- "hs.CRP"
colnames(features)[55] <- "IGFBP.7"
colnames(features)[56] <- "IL.6"
colnames(features)[57] <- "NT.proBNP"
colnames(features)[58] <- "OPN"
colnames(features)[59] <- "hsTropT"
colnames(features)[60] <- "eGFR"

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
colnames(features)[1] <- "Age"
colnames(features)[2] <- "Sex"
colnames(features)[3] <- "BMI"
colnames(features)[4] <- "Smoker"
colnames(features)[5] <- "SBP"
colnames(features)[6] <- "Diabetes"
colnames(features)[7] <- "Prior stroke/TIA"
colnames(features)[8] <- "Heart failure"
colnames(features)[9] <- "Renal failure"
colnames(features)[10] <- "CAD"
colnames(features)[11] <- "Study center"
colnames(features)[12] <- "Ethnicity"
colnames(features)[13] <- "AF type"
colnames(features)[14] <- "Beer drinker"
colnames(features)[15] <- "Red wine drinker"
colnames(features)[16] <- "White wine drinker"
colnames(features)[17] <- "Liquor drinker"
colnames(features)[18] <- "AF episode duration"
colnames(features)[19] <- "Nr. of AF episodes"
colnames(features)[20] <- "Atrial flutter"
colnames(features)[21] <- "Aspirin"
colnames(features)[22] <- "TCA"
colnames(features)[23] <- "Antiplatelet therapy"
colnames(features)[24] <- "Illicit drugs"
colnames(features)[25] <- "Prior CABG"
colnames(features)[26] <- "Hyperthyroidism"
colnames(features)[27] <- "Hypothyroidism"
colnames(features)[28] <- "Cancer"
colnames(features)[29] <- "VTE"
colnames(features)[30] <- "Paternal AF"
colnames(features)[31] <- "Brother AF"
colnames(features)[32] <- "Maternal AF"
colnames(features)[33] <- "Sister AF"
colnames(features)[34] <- "Familial hypertension"
colnames(features)[35] <- "Familial diabetes"
colnames(features)[36] <- "Familial obesity"
colnames(features)[37] <- "Familial CAD"
colnames(features)[38] <- "Height"
colnames(features)[39] <- "Weight"
colnames(features)[40] <- "Heart rate"
colnames(features)[41] <- "Diastolic blood pressure"
colnames(features)[42] <- "Rhythm at baseline"
colnames(features)[43] <- "Sleep apnea"
colnames(features)[44] <- "Hypertension"
colnames(features)[45] <- "PAD"
colnames(features)[46] <- "Prior MI"
colnames(features)[47] <- "Systemic embolism"
colnames(features)[48] <- "Prior major bleeding"
colnames(features)[49] <- "ANG-2"
colnames(features)[50] <- "d-dimer"
colnames(features)[51] <- "cystatin C"
colnames(features)[52] <- "ALAT"
colnames(features)[53] <- "GDF-15"
colnames(features)[54] <- "hs-CRP"
colnames(features)[55] <- "IGFBP-7"
colnames(features)[56] <- "IL-6"
colnames(features)[57] <- "NT-proBNP"
colnames(features)[58] <- "OPN"
colnames(features)[59] <- "hsTropT"
colnames(features)[60] <- "eGFR"

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
colnames(features)[1] <- "Age"
colnames(features)[2] <- "Sex"
colnames(features)[3] <- "BMI"
colnames(features)[4] <- "Smoker"
colnames(features)[5] <- "Systolic blood pressure"
colnames(features)[6] <- "Diabetes"
colnames(features)[7] <- "Prior stroke/TIA"
colnames(features)[8] <- "Heart failure"
colnames(features)[9] <- "Renal failure"
colnames(features)[10] <- "CAD"
colnames(features)[11] <- "Study center"
colnames(features)[12] <- "Ethnicity"
colnames(features)[13] <- "AF type"
colnames(features)[14] <- "Beer drinker"
colnames(features)[15] <- "Red wine drinker"
colnames(features)[16] <- "White wine drinker"
colnames(features)[17] <- "Liquor drinker"
colnames(features)[18] <- "AF episode duration"
colnames(features)[19] <- "Nr. of AF episodes"
colnames(features)[20] <- "Atrial flutter"
colnames(features)[21] <- "Aspirin"
colnames(features)[22] <- "TCA"
colnames(features)[23] <- "Antiplatelet therapy"
colnames(features)[24] <- "Illicit drugs"
colnames(features)[25] <- "Prior CABG"
colnames(features)[26] <- "Hyperthyroidism"
colnames(features)[27] <- "Hypothyroidism"
colnames(features)[28] <- "Cancer"
colnames(features)[29] <- "VTE"
colnames(features)[30] <- "Paternal AF"
colnames(features)[31] <- "Brother AF"
colnames(features)[32] <- "Maternal AF"
colnames(features)[33] <- "Sister AF"
colnames(features)[34] <- "Familial hypertension"
colnames(features)[35] <- "Familial diabetes"
colnames(features)[36] <- "Familial obesity"
colnames(features)[37] <- "Familial CAD"
colnames(features)[38] <- "Height"
colnames(features)[39] <- "Weight"
colnames(features)[40] <- "Heart rate"
colnames(features)[41] <- "Diastolic blood pressure"
colnames(features)[42] <- "Rhythm at baseline"
colnames(features)[43] <- "Sleep apnea"
colnames(features)[44] <- "Hypertension"
colnames(features)[45] <- "PAD"
colnames(features)[46] <- "Prior MI"
colnames(features)[47] <- "Systemic embolism"
colnames(features)[48] <- "Prior major bleeding"

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
colnames(features)[1] <- "Age"
colnames(features)[2] <- "Sex"
colnames(features)[3] <- "BMI"
colnames(features)[4] <- "Smoker"
colnames(features)[5] <- "Systolic blood pressure"
colnames(features)[6] <- "Diabetes"
colnames(features)[7] <- "Prior stroke/TIA"
colnames(features)[8] <- "Heart failure"
colnames(features)[9] <- "Renal failure"
colnames(features)[10] <- "CAD"
colnames(features)[11] <- "Study center"
colnames(features)[12] <- "Ethnicity"
colnames(features)[13] <- "AF type"
colnames(features)[14] <- "Beer drinker"
colnames(features)[15] <- "Red wine drinker"
colnames(features)[16] <- "White wine drinker"
colnames(features)[17] <- "Liquor drinker"
colnames(features)[18] <- "AF episode duration"
colnames(features)[19] <- "Nr. of AF episodes"
colnames(features)[20] <- "Atrial flutter"
colnames(features)[21] <- "Aspirin"
colnames(features)[22] <- "TCA"
colnames(features)[23] <- "Antiplatelet therapy"
colnames(features)[24] <- "Illicit drugs"
colnames(features)[25] <- "Prior CABG"
colnames(features)[26] <- "Hyperthyroidism"
colnames(features)[27] <- "Hypothyroidism"
colnames(features)[28] <- "Cancer"
colnames(features)[29] <- "VTE"
colnames(features)[30] <- "Paternal AF"
colnames(features)[31] <- "Brother AF"
colnames(features)[32] <- "Maternal AF"
colnames(features)[33] <- "Sister AF"
colnames(features)[34] <- "Familial hypertension"
colnames(features)[35] <- "Familial diabetes"
colnames(features)[36] <- "Familial obesity"
colnames(features)[37] <- "Familial CAD"
colnames(features)[38] <- "Height"
colnames(features)[39] <- "Weight"
colnames(features)[40] <- "Heart rate"
colnames(features)[41] <- "Diastolic blood pressure"
colnames(features)[42] <- "Rhythm at baseline"
colnames(features)[43] <- "Sleep apnea"
colnames(features)[44] <- "Hypertension"
colnames(features)[45] <- "PAD"
colnames(features)[46] <- "Prior MI"
colnames(features)[47] <- "Systemic embolism"
colnames(features)[48] <- "Prior major bleeding"
colnames(features)[49] <- "ANG2"
colnames(features)[50] <- "DDI2H"
colnames(features)[51] <- "CYSC"
colnames(features)[52] <- "ALAT"
colnames(features)[53] <- "GDF15"
colnames(features)[54] <- "CRPHS"
colnames(features)[55] <- "IGFBP7"
colnames(features)[56] <- "IL6"
colnames(features)[57] <- "PROBNPII"
colnames(features)[58] <- "OPN"
colnames(features)[59] <- "TNTHS"
colnames(features)[60] <- "eGFR"

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
colnames(features)[1] <- "Age"
colnames(features)[2] <- "Sex"
colnames(features)[3] <- "BMI"
colnames(features)[4] <- "Smoker"
colnames(features)[5] <- "Systolic blood pressure"
colnames(features)[6] <- "Diabetes"
colnames(features)[7] <- "Prior stroke/TIA"
colnames(features)[8] <- "Heart failure"
colnames(features)[9] <- "Renal failure"
colnames(features)[10] <- "CAD"
colnames(features)[11] <- "Study center"
colnames(features)[12] <- "Ethnicity"
colnames(features)[13] <- "AF type"
colnames(features)[14] <- "Beer drinker"
colnames(features)[15] <- "Red wine drinker"
colnames(features)[16] <- "White wine drinker"
colnames(features)[17] <- "Liquor drinker"
colnames(features)[18] <- "AF episode duration"
colnames(features)[19] <- "Nr. of AF episodes"
colnames(features)[20] <- "Atrial flutter"
colnames(features)[21] <- "Aspirin"
colnames(features)[22] <- "TCA"
colnames(features)[23] <- "Antiplatelet therapy"
colnames(features)[24] <- "Illicit drugs"
colnames(features)[25] <- "Prior CABG"
colnames(features)[26] <- "Hyperthyroidism"
colnames(features)[27] <- "Hypothyroidism"
colnames(features)[28] <- "Cancer"
colnames(features)[29] <- "VTE"
colnames(features)[30] <- "Paternal AF"
colnames(features)[31] <- "Brother AF"
colnames(features)[32] <- "Maternal AF"
colnames(features)[33] <- "Sister AF"
colnames(features)[34] <- "Familial hypertension"
colnames(features)[35] <- "Familial diabetes"
colnames(features)[36] <- "Familial obesity"
colnames(features)[37] <- "Familial CAD"
colnames(features)[38] <- "Height"
colnames(features)[39] <- "Weight"
colnames(features)[40] <- "Heart rate"
colnames(features)[41] <- "Diastolic blood pressure"
colnames(features)[42] <- "Rhythm at baseline"
colnames(features)[43] <- "Sleep apnea"
colnames(features)[44] <- "Hypertension"
colnames(features)[45] <- "PAD"
colnames(features)[46] <- "Prior MI"
colnames(features)[47] <- "Systemic embolism"
colnames(features)[48] <- "Prior major bleeding"

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
colnames(features)[1] <- "Age"
colnames(features)[2] <- "Sex"
colnames(features)[3] <- "BMI"
colnames(features)[4] <- "Smoker"
colnames(features)[5] <- "Systolic blood pressure"
colnames(features)[6] <- "Diabetes"
colnames(features)[7] <- "Prior stroke/TIA"
colnames(features)[8] <- "Heart failure"
colnames(features)[9] <- "Renal failure"
colnames(features)[10] <- "CAD"
colnames(features)[11] <- "Study center"
colnames(features)[12] <- "Ethnicity"
colnames(features)[13] <- "AF type"
colnames(features)[14] <- "Beer drinker"
colnames(features)[15] <- "Red wine drinker"
colnames(features)[16] <- "White wine drinker"
colnames(features)[17] <- "Liquor drinker"
colnames(features)[18] <- "AF episode duration"
colnames(features)[19] <- "Nr. of AF episodes"
colnames(features)[20] <- "Atrial flutter"
colnames(features)[21] <- "Aspirin"
colnames(features)[22] <- "TCA"
colnames(features)[23] <- "Antiplatelet therapy"
colnames(features)[24] <- "Illicit drugs"
colnames(features)[25] <- "Prior CABG"
colnames(features)[26] <- "Hyperthyroidism"
colnames(features)[27] <- "Hypothyroidism"
colnames(features)[28] <- "Cancer"
colnames(features)[29] <- "VTE"
colnames(features)[30] <- "Paternal AF"
colnames(features)[31] <- "Brother AF"
colnames(features)[32] <- "Maternal AF"
colnames(features)[33] <- "Sister AF"
colnames(features)[34] <- "Familial hypertension"
colnames(features)[35] <- "Familial diabetes"
colnames(features)[36] <- "Familial obesity"
colnames(features)[37] <- "Familial CAD"
colnames(features)[38] <- "Height"
colnames(features)[39] <- "Weight"
colnames(features)[40] <- "Heart rate"
colnames(features)[41] <- "Diastolic blood pressure"
colnames(features)[42] <- "Rhythm at baseline"
colnames(features)[43] <- "Sleep apnea"
colnames(features)[44] <- "Hypertension"
colnames(features)[45] <- "PAD"
colnames(features)[46] <- "Prior MI"
colnames(features)[47] <- "Systemic embolism"
colnames(features)[48] <- "Prior major bleeding"
colnames(features)[49] <- "ANG2"
colnames(features)[50] <- "DDI2H"
colnames(features)[51] <- "CYSC"
colnames(features)[52] <- "ALAT"
colnames(features)[53] <- "GDF15"
colnames(features)[54] <- "CRPHS"
colnames(features)[55] <- "IGFBP7"
colnames(features)[56] <- "IL6"
colnames(features)[57] <- "PROBNPII"
colnames(features)[58] <- "OPN"
colnames(features)[59] <- "TNTHS"
colnames(features)[60] <- "eGFR"

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
colnames(features)[1] <- "Age"
colnames(features)[2] <- "Sex"
colnames(features)[3] <- "BMI"
colnames(features)[4] <- "Smoker"
colnames(features)[5] <- "Systolic blood pressure"
colnames(features)[6] <- "Diabetes"
colnames(features)[7] <- "Prior stroke/TIA"
colnames(features)[8] <- "Heart failure"
colnames(features)[9] <- "Renal failure"
colnames(features)[10] <- "CAD"
colnames(features)[11] <- "Study center"
colnames(features)[12] <- "Ethnicity"
colnames(features)[13] <- "AF type"
colnames(features)[14] <- "Beer drinker"
colnames(features)[15] <- "Red wine drinker"
colnames(features)[16] <- "White wine drinker"
colnames(features)[17] <- "Liquor drinker"
colnames(features)[18] <- "AF episode duration"
colnames(features)[19] <- "Nr. of AF episodes"
colnames(features)[20] <- "Atrial flutter"
colnames(features)[21] <- "Aspirin"
colnames(features)[22] <- "TCA"
colnames(features)[23] <- "Antiplatelet therapy"
colnames(features)[24] <- "Illicit drugs"
colnames(features)[25] <- "Prior CABG"
colnames(features)[26] <- "Hyperthyroidism"
colnames(features)[27] <- "Hypothyroidism"
colnames(features)[28] <- "Cancer"
colnames(features)[29] <- "VTE"
colnames(features)[30] <- "Paternal AF"
colnames(features)[31] <- "Brother AF"
colnames(features)[32] <- "Maternal AF"
colnames(features)[33] <- "Sister AF"
colnames(features)[34] <- "Familial hypertension"
colnames(features)[35] <- "Familial diabetes"
colnames(features)[36] <- "Familial obesity"
colnames(features)[37] <- "Familial CAD"
colnames(features)[38] <- "Height"
colnames(features)[39] <- "Weight"
colnames(features)[40] <- "Heart rate"
colnames(features)[41] <- "Diastolic blood pressure"
colnames(features)[42] <- "Rhythm at baseline"
colnames(features)[43] <- "Sleep apnea"
colnames(features)[44] <- "Hypertension"
colnames(features)[45] <- "PAD"
colnames(features)[46] <- "Prior MI"
colnames(features)[47] <- "Systemic embolism"
colnames(features)[48] <- "Prior major bleeding"

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
colnames(features)[1] <- "Age"
colnames(features)[2] <- "Sex"
colnames(features)[3] <- "BMI"
colnames(features)[4] <- "Smoker"
colnames(features)[5] <- "Systolic blood pressure"
colnames(features)[6] <- "Diabetes"
colnames(features)[7] <- "Prior stroke/TIA"
colnames(features)[8] <- "Heart failure"
colnames(features)[9] <- "Renal failure"
colnames(features)[10] <- "CAD"
colnames(features)[11] <- "Study center"
colnames(features)[12] <- "Ethnicity"
colnames(features)[13] <- "AF type"
colnames(features)[14] <- "Beer drinker"
colnames(features)[15] <- "Red wine drinker"
colnames(features)[16] <- "White wine drinker"
colnames(features)[17] <- "Liquor drinker"
colnames(features)[18] <- "AF episode duration"
colnames(features)[19] <- "Nr. of AF episodes"
colnames(features)[20] <- "Atrial flutter"
colnames(features)[21] <- "Aspirin"
colnames(features)[22] <- "TCA"
colnames(features)[23] <- "Antiplatelet therapy"
colnames(features)[24] <- "Illicit drugs"
colnames(features)[25] <- "Prior CABG"
colnames(features)[26] <- "Hyperthyroidism"
colnames(features)[27] <- "Hypothyroidism"
colnames(features)[28] <- "Cancer"
colnames(features)[29] <- "VTE"
colnames(features)[30] <- "Paternal AF"
colnames(features)[31] <- "Brother AF"
colnames(features)[32] <- "Maternal AF"
colnames(features)[33] <- "Sister AF"
colnames(features)[34] <- "Familial hypertension"
colnames(features)[35] <- "Familial diabetes"
colnames(features)[36] <- "Familial obesity"
colnames(features)[37] <- "Familial CAD"
colnames(features)[38] <- "Height"
colnames(features)[39] <- "Weight"
colnames(features)[40] <- "Heart rate"
colnames(features)[41] <- "Diastolic blood pressure"
colnames(features)[42] <- "Rhythm at baseline"
colnames(features)[43] <- "Sleep apnea"
colnames(features)[44] <- "Hypertension"
colnames(features)[45] <- "PAD"
colnames(features)[46] <- "Prior MI"
colnames(features)[47] <- "Systemic embolism"
colnames(features)[48] <- "Prior major bleeding"
colnames(features)[49] <- "ANG2"
colnames(features)[50] <- "DDI2H"
colnames(features)[51] <- "CYSC"
colnames(features)[52] <- "ALAT"
colnames(features)[53] <- "GDF15"
colnames(features)[54] <- "CRPHS"
colnames(features)[55] <- "IGFBP7"
colnames(features)[56] <- "IL6"
colnames(features)[57] <- "PROBNPII"
colnames(features)[58] <- "OPN"
colnames(features)[59] <- "TNTHS"
colnames(features)[60] <- "eGFR"

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
colnames(features)[1] <- "Age"
colnames(features)[2] <- "Sex"
colnames(features)[3] <- "BMI"
colnames(features)[4] <- "Smoker"
colnames(features)[5] <- "Systolic blood pressure"
colnames(features)[6] <- "Diabetes"
colnames(features)[7] <- "Prior stroke/TIA"
colnames(features)[8] <- "Heart failure"
colnames(features)[9] <- "Renal failure"
colnames(features)[10] <- "CAD"
colnames(features)[11] <- "Study center"
colnames(features)[12] <- "Ethnicity"
colnames(features)[13] <- "AF type"
colnames(features)[14] <- "Beer drinker"
colnames(features)[15] <- "Red wine drinker"
colnames(features)[16] <- "White wine drinker"
colnames(features)[17] <- "Liquor drinker"
colnames(features)[18] <- "AF episode duration"
colnames(features)[19] <- "Nr. of AF episodes"
colnames(features)[20] <- "Atrial flutter"
colnames(features)[21] <- "Aspirin"
colnames(features)[22] <- "TCA"
colnames(features)[23] <- "Antiplatelet therapy"
colnames(features)[24] <- "Illicit drugs"
colnames(features)[25] <- "Prior CABG"
colnames(features)[26] <- "Hyperthyroidism"
colnames(features)[27] <- "Hypothyroidism"
colnames(features)[28] <- "Cancer"
colnames(features)[29] <- "VTE"
colnames(features)[30] <- "Paternal AF"
colnames(features)[31] <- "Brother AF"
colnames(features)[32] <- "Maternal AF"
colnames(features)[33] <- "Sister AF"
colnames(features)[34] <- "Familial hypertension"
colnames(features)[35] <- "Familial diabetes"
colnames(features)[36] <- "Familial obesity"
colnames(features)[37] <- "Familial CAD"
colnames(features)[38] <- "Height"
colnames(features)[39] <- "Weight"
colnames(features)[40] <- "Heart rate"
colnames(features)[41] <- "Diastolic blood pressure"
colnames(features)[42] <- "Rhythm at baseline"
colnames(features)[43] <- "Sleep apnea"
colnames(features)[44] <- "Hypertension"
colnames(features)[45] <- "PAD"
colnames(features)[46] <- "Prior MI"
colnames(features)[47] <- "Systemic embolism"
colnames(features)[48] <- "Prior major bleeding"

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
colnames(features)[1] <- "Age"
colnames(features)[2] <- "Sex"
colnames(features)[3] <- "BMI"
colnames(features)[4] <- "Smoker"
colnames(features)[5] <- "Systolic blood pressure"
colnames(features)[6] <- "Diabetes"
colnames(features)[7] <- "Prior stroke/TIA"
colnames(features)[8] <- "Heart failure"
colnames(features)[9] <- "Renal failure"
colnames(features)[10] <- "CAD"
colnames(features)[11] <- "Study center"
colnames(features)[12] <- "Ethnicity"
colnames(features)[13] <- "AF type"
colnames(features)[14] <- "Beer drinker"
colnames(features)[15] <- "Red wine drinker"
colnames(features)[16] <- "White wine drinker"
colnames(features)[17] <- "Liquor drinker"
colnames(features)[18] <- "AF episode duration"
colnames(features)[19] <- "Nr. of AF episodes"
colnames(features)[20] <- "Atrial flutter"
colnames(features)[21] <- "Aspirin"
colnames(features)[22] <- "TCA"
colnames(features)[23] <- "Antiplatelet therapy"
colnames(features)[24] <- "Illicit drugs"
colnames(features)[25] <- "Prior CABG"
colnames(features)[26] <- "Hyperthyroidism"
colnames(features)[27] <- "Hypothyroidism"
colnames(features)[28] <- "Cancer"
colnames(features)[29] <- "VTE"
colnames(features)[30] <- "Paternal AF"
colnames(features)[31] <- "Brother AF"
colnames(features)[32] <- "Maternal AF"
colnames(features)[33] <- "Sister AF"
colnames(features)[34] <- "Familial hypertension"
colnames(features)[35] <- "Familial diabetes"
colnames(features)[36] <- "Familial obesity"
colnames(features)[37] <- "Familial CAD"
colnames(features)[38] <- "Height"
colnames(features)[39] <- "Weight"
colnames(features)[40] <- "Heart rate"
colnames(features)[41] <- "Diastolic blood pressure"
colnames(features)[42] <- "Rhythm at baseline"
colnames(features)[43] <- "Sleep apnea"
colnames(features)[44] <- "Hypertension"
colnames(features)[45] <- "PAD"
colnames(features)[46] <- "Prior MI"
colnames(features)[47] <- "Systemic embolism"
colnames(features)[48] <- "Prior major bleeding"
colnames(features)[49] <- "ANG2"
colnames(features)[50] <- "DDI2H"
colnames(features)[51] <- "CYSC"
colnames(features)[52] <- "ALAT"
colnames(features)[53] <- "GDF15"
colnames(features)[54] <- "CRPHS"
colnames(features)[55] <- "IGFBP7"
colnames(features)[56] <- "IL6"
colnames(features)[57] <- "PROBNPII"
colnames(features)[58] <- "OPN"
colnames(features)[59] <- "TNTHS"
colnames(features)[60] <- "eGFR"

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
colnames(features)[1] <- "Age"
colnames(features)[2] <- "Sex"
colnames(features)[3] <- "BMI"
colnames(features)[4] <- "Smoker"
colnames(features)[5] <- "Systolic blood pressure"
colnames(features)[6] <- "Diabetes"
colnames(features)[7] <- "Prior stroke/TIA"
colnames(features)[8] <- "Heart failure"
colnames(features)[9] <- "Renal failure"
colnames(features)[10] <- "CAD"
colnames(features)[11] <- "Study center"
colnames(features)[12] <- "Ethnicity"
colnames(features)[13] <- "AF type"
colnames(features)[14] <- "Beer drinker"
colnames(features)[15] <- "Red wine drinker"
colnames(features)[16] <- "White wine drinker"
colnames(features)[17] <- "Liquor drinker"
colnames(features)[18] <- "AF episode duration"
colnames(features)[19] <- "Nr. of AF episodes"
colnames(features)[20] <- "Atrial flutter"
colnames(features)[21] <- "Aspirin"
colnames(features)[22] <- "TCA"
colnames(features)[23] <- "Antiplatelet therapy"
colnames(features)[24] <- "Illicit drugs"
colnames(features)[25] <- "Prior CABG"
colnames(features)[26] <- "Hyperthyroidism"
colnames(features)[27] <- "Hypothyroidism"
colnames(features)[28] <- "Cancer"
colnames(features)[29] <- "VTE"
colnames(features)[30] <- "Paternal AF"
colnames(features)[31] <- "Brother AF"
colnames(features)[32] <- "Maternal AF"
colnames(features)[33] <- "Sister AF"
colnames(features)[34] <- "Familial hypertension"
colnames(features)[35] <- "Familial diabetes"
colnames(features)[36] <- "Familial obesity"
colnames(features)[37] <- "Familial CAD"
colnames(features)[38] <- "Height"
colnames(features)[39] <- "Weight"
colnames(features)[40] <- "Heart rate"
colnames(features)[41] <- "Diastolic blood pressure"
colnames(features)[42] <- "Rhythm at baseline"
colnames(features)[43] <- "Sleep apnea"
colnames(features)[44] <- "Hypertension"
colnames(features)[45] <- "PAD"
colnames(features)[46] <- "Prior MI"
colnames(features)[47] <- "Systemic embolism"
colnames(features)[48] <- "Prior major bleeding"

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
colnames(features)[1] <- "Age"
colnames(features)[2] <- "Sex"
colnames(features)[3] <- "BMI"
colnames(features)[4] <- "Smoker"
colnames(features)[5] <- "Systolic blood pressure"
colnames(features)[6] <- "Diabetes"
colnames(features)[7] <- "Prior stroke/TIA"
colnames(features)[8] <- "Heart failure"
colnames(features)[9] <- "Renal failure"
colnames(features)[10] <- "CAD"
colnames(features)[11] <- "Study center"
colnames(features)[12] <- "Ethnicity"
colnames(features)[13] <- "AF type"
colnames(features)[14] <- "Beer drinker"
colnames(features)[15] <- "Red wine drinker"
colnames(features)[16] <- "White wine drinker"
colnames(features)[17] <- "Liquor drinker"
colnames(features)[18] <- "AF episode duration"
colnames(features)[19] <- "Nr. of AF episodes"
colnames(features)[20] <- "Atrial flutter"
colnames(features)[21] <- "Aspirin"
colnames(features)[22] <- "TCA"
colnames(features)[23] <- "Antiplatelet therapy"
colnames(features)[24] <- "Illicit drugs"
colnames(features)[25] <- "Prior CABG"
colnames(features)[26] <- "Hyperthyroidism"
colnames(features)[27] <- "Hypothyroidism"
colnames(features)[28] <- "Cancer"
colnames(features)[29] <- "VTE"
colnames(features)[30] <- "Paternal AF"
colnames(features)[31] <- "Brother AF"
colnames(features)[32] <- "Maternal AF"
colnames(features)[33] <- "Sister AF"
colnames(features)[34] <- "Familial hypertension"
colnames(features)[35] <- "Familial diabetes"
colnames(features)[36] <- "Familial obesity"
colnames(features)[37] <- "Familial CAD"
colnames(features)[38] <- "Height"
colnames(features)[39] <- "Weight"
colnames(features)[40] <- "Heart rate"
colnames(features)[41] <- "Diastolic blood pressure"
colnames(features)[42] <- "Rhythm at baseline"
colnames(features)[43] <- "Sleep apnea"
colnames(features)[44] <- "Hypertension"
colnames(features)[45] <- "PAD"
colnames(features)[46] <- "Prior MI"
colnames(features)[47] <- "Systemic embolism"
colnames(features)[48] <- "Prior major bleeding"
colnames(features)[49] <- "ANG2"
colnames(features)[50] <- "DDI2H"
colnames(features)[51] <- "CYSC"
colnames(features)[52] <- "ALAT"
colnames(features)[53] <- "GDF15"
colnames(features)[54] <- "CRPHS"
colnames(features)[55] <- "IGFBP7"
colnames(features)[56] <- "IL6"
colnames(features)[57] <- "PROBNPII"
colnames(features)[58] <- "OPN"
colnames(features)[59] <- "TNTHS"
colnames(features)[60] <- "eGFR"

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
colnames(features)[1] <- "Age"
colnames(features)[2] <- "Sex"
colnames(features)[3] <- "BMI"
colnames(features)[4] <- "Smoker"
colnames(features)[5] <- "Systolic blood pressure"
colnames(features)[6] <- "Diabetes"
colnames(features)[7] <- "Prior stroke/TIA"
colnames(features)[8] <- "Heart failure"
colnames(features)[9] <- "Renal failure"
colnames(features)[10] <- "CAD"
colnames(features)[11] <- "Study center"
colnames(features)[12] <- "Ethnicity"
colnames(features)[13] <- "AF type"
colnames(features)[14] <- "Beer drinker"
colnames(features)[15] <- "Red wine drinker"
colnames(features)[16] <- "White wine drinker"
colnames(features)[17] <- "Liquor drinker"
colnames(features)[18] <- "AF episode duration"
colnames(features)[19] <- "Nr. of AF episodes"
colnames(features)[20] <- "Atrial flutter"
colnames(features)[21] <- "Aspirin"
colnames(features)[22] <- "TCA"
colnames(features)[23] <- "Antiplatelet therapy"
colnames(features)[24] <- "Illicit drugs"
colnames(features)[25] <- "Prior CABG"
colnames(features)[26] <- "Hyperthyroidism"
colnames(features)[27] <- "Hypothyroidism"
colnames(features)[28] <- "Cancer"
colnames(features)[29] <- "VTE"
colnames(features)[30] <- "Paternal AF"
colnames(features)[31] <- "Brother AF"
colnames(features)[32] <- "Maternal AF"
colnames(features)[33] <- "Sister AF"
colnames(features)[34] <- "Familial hypertension"
colnames(features)[35] <- "Familial diabetes"
colnames(features)[36] <- "Familial obesity"
colnames(features)[37] <- "Familial CAD"
colnames(features)[38] <- "Height"
colnames(features)[39] <- "Weight"
colnames(features)[40] <- "Heart rate"
colnames(features)[41] <- "Diastolic blood pressure"
colnames(features)[42] <- "Rhythm at baseline"
colnames(features)[43] <- "Sleep apnea"
colnames(features)[44] <- "Hypertension"
colnames(features)[45] <- "PAD"
colnames(features)[46] <- "Prior MI"
colnames(features)[47] <- "Systemic embolism"
colnames(features)[48] <- "Prior major bleeding"

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
colnames(features)[1] <- "Age"
colnames(features)[2] <- "Sex"
colnames(features)[3] <- "BMI"
colnames(features)[4] <- "Smoker"
colnames(features)[5] <- "Systolic blood pressure"
colnames(features)[6] <- "Diabetes"
colnames(features)[7] <- "Prior stroke/TIA"
colnames(features)[8] <- "Heart failure"
colnames(features)[9] <- "Renal failure"
colnames(features)[10] <- "CAD"
colnames(features)[11] <- "Study center"
colnames(features)[12] <- "Ethnicity"
colnames(features)[13] <- "AF type"
colnames(features)[14] <- "Beer drinker"
colnames(features)[15] <- "Red wine drinker"
colnames(features)[16] <- "White wine drinker"
colnames(features)[17] <- "Liquor drinker"
colnames(features)[18] <- "AF episode duration"
colnames(features)[19] <- "Nr. of AF episodes"
colnames(features)[20] <- "Atrial flutter"
colnames(features)[21] <- "Aspirin"
colnames(features)[22] <- "TCA"
colnames(features)[23] <- "Antiplatelet therapy"
colnames(features)[24] <- "Illicit drugs"
colnames(features)[25] <- "Prior CABG"
colnames(features)[26] <- "Hyperthyroidism"
colnames(features)[27] <- "Hypothyroidism"
colnames(features)[28] <- "Cancer"
colnames(features)[29] <- "VTE"
colnames(features)[30] <- "Paternal AF"
colnames(features)[31] <- "Brother AF"
colnames(features)[32] <- "Maternal AF"
colnames(features)[33] <- "Sister AF"
colnames(features)[34] <- "Familial hypertension"
colnames(features)[35] <- "Familial diabetes"
colnames(features)[36] <- "Familial obesity"
colnames(features)[37] <- "Familial CAD"
colnames(features)[38] <- "Height"
colnames(features)[39] <- "Weight"
colnames(features)[40] <- "Heart rate"
colnames(features)[41] <- "Diastolic blood pressure"
colnames(features)[42] <- "Rhythm at baseline"
colnames(features)[43] <- "Sleep apnea"
colnames(features)[44] <- "Hypertension"
colnames(features)[45] <- "PAD"
colnames(features)[46] <- "Prior MI"
colnames(features)[47] <- "Systemic embolism"
colnames(features)[48] <- "Prior major bleeding"
colnames(features)[49] <- "ANG2"
colnames(features)[50] <- "DDI2H"
colnames(features)[51] <- "CYSC"
colnames(features)[52] <- "ALAT"
colnames(features)[53] <- "GDF15"
colnames(features)[54] <- "CRPHS"
colnames(features)[55] <- "IGFBP7"
colnames(features)[56] <- "IL6"
colnames(features)[57] <- "PROBNPII"
colnames(features)[58] <- "OPN"
colnames(features)[59] <- "TNTHS"
colnames(features)[60] <- "eGFR"

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
colnames(features)[1] <- "Age"
colnames(features)[2] <- "Sex"
colnames(features)[3] <- "BMI"
colnames(features)[4] <- "Smoker"
colnames(features)[5] <- "Systolic blood pressure"
colnames(features)[6] <- "Diabetes"
colnames(features)[7] <- "Prior stroke/TIA"
colnames(features)[8] <- "Heart failure"
colnames(features)[9] <- "Renal failure"
colnames(features)[10] <- "CAD"
colnames(features)[11] <- "Study center"
colnames(features)[12] <- "Ethnicity"
colnames(features)[13] <- "AF type"
colnames(features)[14] <- "Beer drinker"
colnames(features)[15] <- "Red wine drinker"
colnames(features)[16] <- "White wine drinker"
colnames(features)[17] <- "Liquor drinker"
colnames(features)[18] <- "AF episode duration"
colnames(features)[19] <- "Nr. of AF episodes"
colnames(features)[20] <- "Atrial flutter"
colnames(features)[21] <- "Aspirin"
colnames(features)[22] <- "TCA"
colnames(features)[23] <- "Antiplatelet therapy"
colnames(features)[24] <- "Illicit drugs"
colnames(features)[25] <- "Prior CABG"
colnames(features)[26] <- "Hyperthyroidism"
colnames(features)[27] <- "Hypothyroidism"
colnames(features)[28] <- "Cancer"
colnames(features)[29] <- "VTE"
colnames(features)[30] <- "Paternal AF"
colnames(features)[31] <- "Brother AF"
colnames(features)[32] <- "Maternal AF"
colnames(features)[33] <- "Sister AF"
colnames(features)[34] <- "Familial hypertension"
colnames(features)[35] <- "Familial diabetes"
colnames(features)[36] <- "Familial obesity"
colnames(features)[37] <- "Familial CAD"
colnames(features)[38] <- "Height"
colnames(features)[39] <- "Weight"
colnames(features)[40] <- "Heart rate"
colnames(features)[41] <- "Diastolic blood pressure"
colnames(features)[42] <- "Rhythm at baseline"
colnames(features)[43] <- "Sleep apnea"
colnames(features)[44] <- "Hypertension"
colnames(features)[45] <- "PAD"
colnames(features)[46] <- "Prior MI"
colnames(features)[47] <- "Systemic embolism"
colnames(features)[48] <- "Prior major bleeding"

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
colnames(features)[1] <- "Age"
colnames(features)[2] <- "Sex"
colnames(features)[3] <- "BMI"
colnames(features)[4] <- "Smoker"
colnames(features)[5] <- "Systolic blood pressure"
colnames(features)[6] <- "Diabetes"
colnames(features)[7] <- "Prior stroke/TIA"
colnames(features)[8] <- "Heart failure"
colnames(features)[9] <- "Renal failure"
colnames(features)[10] <- "CAD"
colnames(features)[11] <- "Study center"
colnames(features)[12] <- "Ethnicity"
colnames(features)[13] <- "AF type"
colnames(features)[14] <- "Beer drinker"
colnames(features)[15] <- "Red wine drinker"
colnames(features)[16] <- "White wine drinker"
colnames(features)[17] <- "Liquor drinker"
colnames(features)[18] <- "AF episode duration"
colnames(features)[19] <- "Nr. of AF episodes"
colnames(features)[20] <- "Atrial flutter"
colnames(features)[21] <- "Aspirin"
colnames(features)[22] <- "TCA"
colnames(features)[23] <- "Antiplatelet therapy"
colnames(features)[24] <- "Illicit drugs"
colnames(features)[25] <- "Prior CABG"
colnames(features)[26] <- "Hyperthyroidism"
colnames(features)[27] <- "Hypothyroidism"
colnames(features)[28] <- "Cancer"
colnames(features)[29] <- "VTE"
colnames(features)[30] <- "Paternal AF"
colnames(features)[31] <- "Brother AF"
colnames(features)[32] <- "Maternal AF"
colnames(features)[33] <- "Sister AF"
colnames(features)[34] <- "Familial hypertension"
colnames(features)[35] <- "Familial diabetes"
colnames(features)[36] <- "Familial obesity"
colnames(features)[37] <- "Familial CAD"
colnames(features)[38] <- "Height"
colnames(features)[39] <- "Weight"
colnames(features)[40] <- "Heart rate"
colnames(features)[41] <- "Diastolic blood pressure"
colnames(features)[42] <- "Rhythm at baseline"
colnames(features)[43] <- "Sleep apnea"
colnames(features)[44] <- "Hypertension"
colnames(features)[45] <- "PAD"
colnames(features)[46] <- "Prior MI"
colnames(features)[47] <- "Systemic embolism"
colnames(features)[48] <- "Prior major bleeding"
colnames(features)[49] <- "ANG2"
colnames(features)[50] <- "DDI2H"
colnames(features)[51] <- "CYSC"
colnames(features)[52] <- "ALAT"
colnames(features)[53] <- "GDF15"
colnames(features)[54] <- "CRPHS"
colnames(features)[55] <- "IGFBP7"
colnames(features)[56] <- "IL6"
colnames(features)[57] <- "PROBNPII"
colnames(features)[58] <- "OPN"
colnames(features)[59] <- "TNTHS"
colnames(features)[60] <- "eGFR"

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
colnames(features)[1] <- "Age"
colnames(features)[2] <- "Sex"
colnames(features)[3] <- "BMI"
colnames(features)[4] <- "Smoker"
colnames(features)[5] <- "Systolic blood pressure"
colnames(features)[6] <- "Diabetes"
colnames(features)[7] <- "Prior stroke/TIA"
colnames(features)[8] <- "Heart failure"
colnames(features)[9] <- "Renal failure"
colnames(features)[10] <- "CAD"
colnames(features)[11] <- "Study center"
colnames(features)[12] <- "Ethnicity"
colnames(features)[13] <- "AF type"
colnames(features)[14] <- "Beer drinker"
colnames(features)[15] <- "Red wine drinker"
colnames(features)[16] <- "White wine drinker"
colnames(features)[17] <- "Liquor drinker"
colnames(features)[18] <- "AF episode duration"
colnames(features)[19] <- "Nr. of AF episodes"
colnames(features)[20] <- "Atrial flutter"
colnames(features)[21] <- "Aspirin"
colnames(features)[22] <- "TCA"
colnames(features)[23] <- "Antiplatelet therapy"
colnames(features)[24] <- "Illicit drugs"
colnames(features)[25] <- "Prior CABG"
colnames(features)[26] <- "Hyperthyroidism"
colnames(features)[27] <- "Hypothyroidism"
colnames(features)[28] <- "Cancer"
colnames(features)[29] <- "VTE"
colnames(features)[30] <- "Paternal AF"
colnames(features)[31] <- "Brother AF"
colnames(features)[32] <- "Maternal AF"
colnames(features)[33] <- "Sister AF"
colnames(features)[34] <- "Familial hypertension"
colnames(features)[35] <- "Familial diabetes"
colnames(features)[36] <- "Familial obesity"
colnames(features)[37] <- "Familial CAD"
colnames(features)[38] <- "Height"
colnames(features)[39] <- "Weight"
colnames(features)[40] <- "Heart rate"
colnames(features)[41] <- "Diastolic blood pressure"
colnames(features)[42] <- "Rhythm at baseline"
colnames(features)[43] <- "Sleep apnea"
colnames(features)[44] <- "Hypertension"
colnames(features)[45] <- "PAD"
colnames(features)[46] <- "Prior MI"
colnames(features)[47] <- "Systemic embolism"
colnames(features)[48] <- "Prior major bleeding"

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
colnames(features)[1] <- "Age"
colnames(features)[2] <- "Sex"
colnames(features)[3] <- "BMI"
colnames(features)[4] <- "Smoker"
colnames(features)[5] <- "Systolic blood pressure"
colnames(features)[6] <- "Diabetes"
colnames(features)[7] <- "Prior stroke/TIA"
colnames(features)[8] <- "Heart failure"
colnames(features)[9] <- "Renal failure"
colnames(features)[10] <- "CAD"
colnames(features)[11] <- "Study center"
colnames(features)[12] <- "Ethnicity"
colnames(features)[13] <- "AF type"
colnames(features)[14] <- "Beer drinker"
colnames(features)[15] <- "Red wine drinker"
colnames(features)[16] <- "White wine drinker"
colnames(features)[17] <- "Liquor drinker"
colnames(features)[18] <- "AF episode duration"
colnames(features)[19] <- "Nr. of AF episodes"
colnames(features)[20] <- "Atrial flutter"
colnames(features)[21] <- "Aspirin"
colnames(features)[22] <- "TCA"
colnames(features)[23] <- "Antiplatelet therapy"
colnames(features)[24] <- "Illicit drugs"
colnames(features)[25] <- "Prior CABG"
colnames(features)[26] <- "Hyperthyroidism"
colnames(features)[27] <- "Hypothyroidism"
colnames(features)[28] <- "Cancer"
colnames(features)[29] <- "VTE"
colnames(features)[30] <- "Paternal AF"
colnames(features)[31] <- "Brother AF"
colnames(features)[32] <- "Maternal AF"
colnames(features)[33] <- "Sister AF"
colnames(features)[34] <- "Familial hypertension"
colnames(features)[35] <- "Familial diabetes"
colnames(features)[36] <- "Familial obesity"
colnames(features)[37] <- "Familial CAD"
colnames(features)[38] <- "Height"
colnames(features)[39] <- "Weight"
colnames(features)[40] <- "Heart rate"
colnames(features)[41] <- "Diastolic blood pressure"
colnames(features)[42] <- "Rhythm at baseline"
colnames(features)[43] <- "Sleep apnea"
colnames(features)[44] <- "Hypertension"
colnames(features)[45] <- "PAD"
colnames(features)[46] <- "Prior MI"
colnames(features)[47] <- "Systemic embolism"
colnames(features)[48] <- "Prior major bleeding"
colnames(features)[49] <- "ANG2"
colnames(features)[50] <- "DDI2H"
colnames(features)[51] <- "CYSC"
colnames(features)[52] <- "ALAT"
colnames(features)[53] <- "GDF15"
colnames(features)[54] <- "CRPHS"
colnames(features)[55] <- "IGFBP7"
colnames(features)[56] <- "IL6"
colnames(features)[57] <- "PROBNPII"
colnames(features)[58] <- "OPN"
colnames(features)[59] <- "TNTHS"
colnames(features)[60] <- "eGFR"

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
colnames(features)[1] <- "Age"
colnames(features)[2] <- "Sex"
colnames(features)[3] <- "BMI"
colnames(features)[4] <- "Smoker"
colnames(features)[5] <- "Systolic blood pressure"
colnames(features)[6] <- "Diabetes"
colnames(features)[7] <- "Prior stroke/TIA"
colnames(features)[8] <- "Heart failure"
colnames(features)[9] <- "Renal failure"
colnames(features)[10] <- "CAD"
colnames(features)[11] <- "Study center"
colnames(features)[12] <- "Ethnicity"
colnames(features)[13] <- "AF type"
colnames(features)[14] <- "Beer drinker"
colnames(features)[15] <- "Red wine drinker"
colnames(features)[16] <- "White wine drinker"
colnames(features)[17] <- "Liquor drinker"
colnames(features)[18] <- "AF episode duration"
colnames(features)[19] <- "Nr. of AF episodes"
colnames(features)[20] <- "Atrial flutter"
colnames(features)[21] <- "Aspirin"
colnames(features)[22] <- "TCA"
colnames(features)[23] <- "Antiplatelet therapy"
colnames(features)[24] <- "Illicit drugs"
colnames(features)[25] <- "Prior CABG"
colnames(features)[26] <- "Hyperthyroidism"
colnames(features)[27] <- "Hypothyroidism"
colnames(features)[28] <- "Cancer"
colnames(features)[29] <- "VTE"
colnames(features)[30] <- "Paternal AF"
colnames(features)[31] <- "Brother AF"
colnames(features)[32] <- "Maternal AF"
colnames(features)[33] <- "Sister AF"
colnames(features)[34] <- "Familial hypertension"
colnames(features)[35] <- "Familial diabetes"
colnames(features)[36] <- "Familial obesity"
colnames(features)[37] <- "Familial CAD"
colnames(features)[38] <- "Height"
colnames(features)[39] <- "Weight"
colnames(features)[40] <- "Heart rate"
colnames(features)[41] <- "Diastolic blood pressure"
colnames(features)[42] <- "Rhythm at baseline"
colnames(features)[43] <- "Sleep apnea"
colnames(features)[44] <- "Hypertension"
colnames(features)[45] <- "PAD"
colnames(features)[46] <- "Prior MI"
colnames(features)[47] <- "Systemic embolism"
colnames(features)[48] <- "Prior major bleeding"

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
colnames(features)[1] <- "Age"
colnames(features)[2] <- "Sex"
colnames(features)[3] <- "BMI"
colnames(features)[4] <- "Smoker"
colnames(features)[5] <- "Systolic blood pressure"
colnames(features)[6] <- "Diabetes"
colnames(features)[7] <- "Prior stroke/TIA"
colnames(features)[8] <- "Heart failure"
colnames(features)[9] <- "Renal failure"
colnames(features)[10] <- "CAD"
colnames(features)[11] <- "Study center"
colnames(features)[12] <- "Ethnicity"
colnames(features)[13] <- "AF type"
colnames(features)[14] <- "Beer drinker"
colnames(features)[15] <- "Red wine drinker"
colnames(features)[16] <- "White wine drinker"
colnames(features)[17] <- "Liquor drinker"
colnames(features)[18] <- "AF episode duration"
colnames(features)[19] <- "Nr. of AF episodes"
colnames(features)[20] <- "Atrial flutter"
colnames(features)[21] <- "Aspirin"
colnames(features)[22] <- "TCA"
colnames(features)[23] <- "Antiplatelet therapy"
colnames(features)[24] <- "Illicit drugs"
colnames(features)[25] <- "Prior CABG"
colnames(features)[26] <- "Hyperthyroidism"
colnames(features)[27] <- "Hypothyroidism"
colnames(features)[28] <- "Cancer"
colnames(features)[29] <- "VTE"
colnames(features)[30] <- "Paternal AF"
colnames(features)[31] <- "Brother AF"
colnames(features)[32] <- "Maternal AF"
colnames(features)[33] <- "Sister AF"
colnames(features)[34] <- "Familial hypertension"
colnames(features)[35] <- "Familial diabetes"
colnames(features)[36] <- "Familial obesity"
colnames(features)[37] <- "Familial CAD"
colnames(features)[38] <- "Height"
colnames(features)[39] <- "Weight"
colnames(features)[40] <- "Heart rate"
colnames(features)[41] <- "Diastolic blood pressure"
colnames(features)[42] <- "Rhythm at baseline"
colnames(features)[43] <- "Sleep apnea"
colnames(features)[44] <- "Hypertension"
colnames(features)[45] <- "PAD"
colnames(features)[46] <- "Prior MI"
colnames(features)[47] <- "Systemic embolism"
colnames(features)[48] <- "Prior major bleeding"
colnames(features)[49] <- "ANG2"
colnames(features)[50] <- "DDI2H"
colnames(features)[51] <- "CYSC"
colnames(features)[52] <- "ALAT"
colnames(features)[53] <- "GDF15"
colnames(features)[54] <- "CRPHS"
colnames(features)[55] <- "IGFBP7"
colnames(features)[56] <- "IL6"
colnames(features)[57] <- "PROBNPII"
colnames(features)[58] <- "OPN"
colnames(features)[59] <- "TNTHS"
colnames(features)[60] <- "eGFR"

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
