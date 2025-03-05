#--MULTIPLE IMPUTATION----------------------------------------------------------
# Author: Pascal B. Meyre
# Date: 01/28/24
# Location: Reinach, Baselland, Switzerland

# Use the merged SWISS.BEAT-biomarker dataset.
#-------------------------------------------------------------------------------

################################################################################
# Perform multiple imputation with mice
################################################################################

pMiss <- function(x){sum(is.na(x))/length(x)*100}
apply(biomarker_data,2,pMiss)
apply(biomarker_data,1,pMiss)

# alat is missing 10.6% of the data points, the rest has between 0.3-3.3% missing data points

# Exclude 'pat.id' variable before imputation
biomarker_data_for_imputation <- biomarker_data[, -1]  # Exclude the first column (pat.id)

# Using mice package for looking at missing data pattern
library(mice)

par(cex=0.6)
md.pattern(biomarker_data_for_imputation, plot = TRUE, rotate.names = TRUE)

# Use mice to impute the missing data
tempData <- mice(biomarker_data_for_imputation, m = 5, maxit = 50, meth = 'pmm', seed = 500)
summary(tempData)

# Visualization of the imputed dataset
densityplot(tempData)
stripplot(tempData, pch = 20, cex = 1.2, xlab = "Imputation number")

completedData <- complete(tempData,1)

# Add 'pat.id' variable back to the imputed dataset
imputed_data <- cbind(biomarker_data[, "pat.id", drop = FALSE], completedData)

# Now 'imputed_data' contains the imputed data with 'pat.id' variable appended

################################# END ##########################################
