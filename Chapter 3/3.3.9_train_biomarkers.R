
### old script
# Load necessary libraries
detach("package:Omix", unload=TRUE)
detach("package:caret", unload=TRUE)
library(caret)


predictors <- c("tTau.ELISA", "AB42.ELISA", "neuro_glia", "Age", "Sex","MoCA")
predictors <- c("tTau.ELISA", "AB42.ELISA", "Age", "Sex","MoCA")

predictors <- c("tTau.ELISA", "AB42.ELISA", "neuro_glia", "Age", "Sex")
predictors <- c("tTau.ELISA", "AB42.ELISA", "Age", "Sex")
X0_Traits <- readRDS("~/aim2/ADNI/models/X0_Traits.rds")
X0_Traits$outcome <- as.factor(X0_Traits$outcome)

# Ensure that the outcome is treated as a factor with 4 levels
X0_Traits <- na.omit(X0_Traits[, c("outcome", predictors)])  # Remove missing data rows 296 to 264

# Set seed for reproducibility
set.seed(123)

# Define predictors for the models

# Control object for 10-fold cross-validation with ROC as the evaluation metric
control <- trainControl(method="cv", number=10, classProbs=TRUE, summaryFunction=multiClassSummary, savePredictions=TRUE,sampling = "smote")

# Train logistic regression model using the selected predictors
model_cv <- train(outcome ~ ., data = X0_Traits[, c("outcome", predictors)], method = "rf", trControl = control)


# Aggregate predictions
predictions <- model_cv$pred
library(pROC)
# Calculate ROC curve for each outcome class
roc_control_resilient <- roc(predictions$obs == "Control_Resilient", predictions$Control_Resilient)
roc_control_vulnerable <- roc(predictions$obs == "Control_Vulnerable", predictions$Control_Vulnerable)
roc_ad_resilient <- roc(predictions$obs == "AD_Resilient", predictions$AD_Resilient)
roc_ad_vulnerable <- roc(predictions$obs == "AD_Vulnerable", predictions$AD_Vulnerable)

# Calculate AUC for each class
auc_control_resilient <- auc(roc_control_resilient)
auc_control_vulnerable <- auc(roc_control_vulnerable)
auc_ad_resilient <- auc(roc_ad_resilient)
auc_ad_vulnerable <- auc(roc_ad_vulnerable)

# Plot ROC curves for all outcome classes
plot(roc_control_resilient, col = "blue", main = "ROC Curves for Cognitive Resilience Classification")
text(x = 0.3, y = 0.2, labels = paste("Control_Resilient AUC:", round(auc_control_resilient, 2)), col = "blue")

plot(roc_control_vulnerable, add = TRUE, col = "red")
text(x = 0.3, y = 0.15, labels = paste("Control_Vulnerable AUC:", round(auc_control_vulnerable, 2)), col = "red")

plot(roc_ad_resilient, add = TRUE, col = "green")
text(x = 0.3, y = 0.10, labels = paste("AD_Resilient AUC:", round(auc_ad_resilient, 2)), col = "green")

plot(roc_ad_vulnerable, add = TRUE, col = "orange")
text(x = 0.3, y = 0.05, labels = paste("AD_Vulnerable AUC:", round(auc_ad_vulnerable, 2)), col = "orange")



saveRDS(model_cv,"tTau.ELISA_AB42.ELISA_Age_Sex.rds")
