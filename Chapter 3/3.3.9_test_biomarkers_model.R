X3_cleanDat <- read.csv("3.cleanDat.csv")
X0_Traits <- read.csv("0.Traits.csv")

X3_cleanDat$X=gsub("\\|.*", "", X3_cleanDat$X)
X0_Traits=X0_Traits[which(X0_Traits$ClinicalGroup!="GIS"),]



X3_cleanDat=X3_cleanDat[which(X3_cleanDat$X %in% c(neuro,micro,astro)),]
rownames(X3_cleanDat)=X3_cleanDat$X
X3_cleanDat=t(X3_cleanDat)
X3_cleanDat=data.frame(X3_cleanDat)
X3_cleanDat=X3_cleanDat[X0_Traits$batch.channel,]



X0_Traits=cbind(X0_Traits,X3_cleanDat)

X0_Traits[, 13:23] <- lapply(X0_Traits[, 13:23], as.numeric)
X0_Traits[, 13:23] <- lapply(X0_Traits[, 13:23], scale)

X0_Traits$micro=X0_Traits$C3+X0_Traits$CD14
X0_Traits$neuro=X0_Traits$GAP43+X0_Traits$NPY
X0_Traits$astro=X0_Traits$AZGP1
X0_Traits$glia=X0_Traits$micro+X0_Traits$astro
X0_Traits$neuro_glia=X0_Traits$neuro-X0_Traits$glia

X0_Traits$tTau.ELISA=as.numeric(X0_Traits$Tau)
X0_Traits$AB42.ELISA=as.numeric(X0_Traits$Ab42)
X0_Traits$MoCA=as.numeric(X0_Traits$Moca)
table(X0_Traits$ClinicalGroup)

cohort2 <- read_csv("models/cohort2.csv")
X0_Traits$ClinicalGroup2=cohort2$`Biomarker-Based Group`[match(X0_Traits$batch.channel,cohort2$`TMT Batch.Channel`)]
X0_Traits$ClinicalGroupOG=X0_Traits$ClinicalGroup
table(X0_Traits$ClinicalGroup)
X0_Traits$ClinicalGroup=X0_Traits$ClinicalGroup2

saveRDS(X0_Traits,"X0_Traits_C2.rds")
#32 cognitively healthy, 31 AsymAD and 33 mild cognitive impairment/AD
X0_Traits <- readRDS("~/aim2/ADNI/models/X0_Traits_C2.rds")
predictors <- c("tTau.ELISA", "AB42.ELISA", "neuro_glia", "Age", "Sex","MoCA")
model_cv <- readRDS("~/aim2/ADNI/models/tTau.ELISA_AB42.ELISA_Age_Sex_MoCA_neuro_glia.rds")

predictors <- c("tTau.ELISA", "AB42.ELISA", "Age", "Sex","MoCA")
model_cv <- readRDS("~/aim2/ADNI/models/tTau.ELISA_AB42.ELISA_Age_Sex_MoCA.rds")

predictors <- c("tTau.ELISA", "AB42.ELISA", "neuro_glia", "Age", "Sex")
model_cv <- readRDS("~/aim2/ADNI/models/tTau.ELISA_AB42.ELISA_Age_Sex_neuro_glia.rds")


predictors <- c("tTau.ELISA", "AB42.ELISA", "Age", "Sex")
model_cv <- readRDS("~/aim2/ADNI/models/tTau.ELISA_AB42.ELISA_Age_Sex.rds")




# Ensure that the outcome is treated as a factor with 4 levels
X0_Traits <- na.omit(X0_Traits[, c( predictors,"batch.channel","ClinicalGroup")])  # Remove missing data rows 296 to 264

# Apply the trained model to the independent dataset
# Ensure the independent dataset has the same columns as the training data (without the outcome variable)
predictions <- predict(model_cv, newdata = X0_Traits[, predictors])

# View the predicted classes
print(predictions)

X0_Traits$predictions=predictions
#saveRDS(X0_Traits,"X0_Traits_PREDICTED-baselinemodelmoca.rds")

table(X0_Traits$predictions,X0_Traits$ClinicalGroup)


# Create simplified predicted classes
predicted_class <- ifelse(predictions %in% c("AD_Resilient", "AD_Vulnerable"), "AD", "Control")

# Create simplified actual classes
actual_class <- ifelse(X0_Traits$ClinicalGroup %in% c("AD"), "AD", "Control")
table(actual_class)
# Create a contingency table
# Create a data frame
data <- data.frame(
  Actual = actual_class,
  Predicted = predicted_class
)

# Convert to factors
data$Actual <- factor(data$Actual, levels = c("AD", "Control"))
data$Predicted <- factor(data$Predicted, levels = c("AD", "Control"))
table(data$Predicted,data$Actual)
# Load caret package
library(caret)

# Create confusion matrix
conf_matrix <- confusionMatrix(data = data$Predicted, reference = data$Actual, positive = "AD")
print(conf_matrix)


