###### Sample Size Calculation ######
event_rate <- sum(as.numeric(train_data$composite_outcome))/sum(train_data$survival_time)
mean(train_data$survival_time)
#R2 calculation
time_person <- 4.128735*160870
expectedevents <- event_rate*time_person
lnnull<- (expectedevents*(log(expectedevents/160870)))-expectedevents
max_r2a <- (1-exp((2*lnnull)/160870))
R2<- max_r2a*0.15
R2
pmsampsize(type="s", rsquared = R2, parameters = 78, rate = event_rate, timepoint = 1, meanfup = 4.128735)



###### Table 1 ######
library(dplyr)
install.packages("table1")
library(table1)
dtable <- read.csv("ready_toimpute.csv")
dtable <- dtable %>%
  mutate(survival_time = as.numeric(as.Date(composite_outcome_date)-as.Date(indexdate)))%>%
  filter(survival_time < 7000) %>%
  mutate(survival_time = round(survival_time / 365, 2))%>%
  filter(survival_time >0)



tbl <- table1(~X_age+ X_time_since_HF +as.factor(region)+gender+ ethnicity+as.factor(imd10_pat)+serum_creatinine+bmi_cat+smoking_cat+
                + e_gfr+as.factor(aki_hx)+ me_sbp + me_dbp+ me_haemogl + me_heartrate + me_mcv+ me_albumin+me_potassium+
                +me_urea+ me_sodium+as.factor(co_afib)+as.factor(co_cardiomyop)+as.factor(co_liver_disease)+as.factor(co_diab_t1)+as.factor(co_ischaemic)+as.factor(co_diab_t2)+as.factor(co_diab_unspec)+
                as.factor(co_diab_any)+as.factor(co_glomer)+as.factor(co_ischaemic)+as.factor(co_lupus)+as.factor(co_nephritis)+as.factor(co_proteinuria)+as.factor(co_pvd)+as.factor(co_valvular)+as.factor(co_ventric_hypertr)+
                as.factor(ev_nephrect_hx)+as.factor(ev_stones_hx)+as.factor(ev_nephrot_synd_hx)+as.factor(ev_ventric_tachy_hx)+as.factor(ep_hyperten)+as.factor(ep_hyperkal)+as.factor(dr_antimicrob_cat)+
                as.factor(dr_entresto_cat)+as.factor(dr_ace_cat)+as.factor(dr_arb_cat)+as.factor(dr_alpha_cat)+as.factor(dr_beta_cat)+as.factor(dr_calciumch_cat)+as.factor(dr_central_cat)+as.factor(dr_renin_cat)+as.factor(dr_loop_diur_cat)+
                as.factor(dr_pai_cat)+as.factor(dr_potass_diur_cat)+as.factor(dr_thiaz_diur_cat)+as.factor(dr_vasodil_cat)+as.factor(dr_immuno_cat)+as.factor(dr_nsaid_cat)+as.factor(dr_proton_cat)+as.factor(dr_sglt2i_cat)+survival_time| as.factor(composite_outcome),
              data= dtable)
print(tbl)



###### PH Assumption ######
library(dplyr)
library(survival)

## Prepare Data
data <- read.csv("fpsm.csv")
data <- data %>%
  filter(survival_time < 7000) %>%
  mutate(survival_time = round(survival_time / 365, 2))%>%
  filter(survival_time >0)

# Splitting the dataset into training and validation sets
set.seed(123)
train_index <- createDataPartition(data$composite_outcome, p = 0.8, list = FALSE)
train_data <- data[train_index, ]
validation_data <- data[-train_index, ]

#test Proportional hazards assumption 
# Fit the Cox proportional hazards model using the original survival object
cox_fit <- coxph(Surv(survival_time, composite_outcome == 1) ~ X_age + e_gfr + me_mcv + X_time_since_HF +
                   me_urea + me_heartrate + me_haemogl + serum_creatinine + me_sbp + me_dbp +
                   me_potassium + me_albumin + me_sodium + as.factor(imd10_pat) + as.factor(region) + numeric_bmi +
                   numeric_smoking + co_ischaemic + co_afib + numeric_gender + co_valvular + ep_hyperten +
                   numeric_ethnicity + aki_hx + co_proteinuria,
                 data = train_data)
# Calculate Schoenfeld residuals
schoenfeld_res <- cox.zph(cox_fit)
# Print the Schoenfeld test results
print(schoenfeld_res)

##martingale residuals - linearity
X <- train_data$X_age
Y <- resid(cox_fit, type = "martingale")
plot(X,Y,pch=20, col= "darkgrey",
     xlab="Age", ylab= "Martingale Residual",
     main = "Residuals vs. Predictor")
abline(h=0)
line(smooth.spline(X,Y, df=3), lty=2, lwd=2)

# Reduce figure margins
par(mfrow = c(3, 5))
par(mar = c(2, 2, 2, 2))  # Adjust the margins as needed

# List of variables in your model
variables <- c("X_age", "e_gfr", "me_mcv", "X_time_since_HF", "me_urea",
               "me_heartrate", "me_haemogl", "serum_creatinine", "me_sbp",
               "me_dbp", "me_potassium", "me_albumin", "me_sodium")
# Customize Schoenfeld residuals plots with triangles
# Customize Schoenfeld residuals plots with triangles
for (var in variables) {
  plot(schoenfeld_res, var = var, main = var, xlab = "Time", ylab = "Schoenfeld Residuals",
       col = ifelse(var %in% c("X_age", "me_albumin", "X_time_since_HF", "region"), "red", "blue"),
       pch = ifelse(var %in% c("X_age", "me_albumin", "X_time_since_HF", "region"), 3, 2),
       lty = ifelse(var %in% c("X_age", "me_albumin", "X_time_since_HF", "region"), 3, 2))
  abline(0,0, col=1, lty=3, lwd=2)
  abline(h = cox_fit$coefficients[2], col = 3, lwd=2,lty = 2)  # Add horizontal line at 0
  grid()  # Add grid lines
}
legend("bottomright",
       legend= c('Reference line for null effect',
                 'Average hazard over time',
                 'Time-Varying hazard'),
       lty = c(3,2,1), col= c(1,3,"red"),lwd=2)

###### MODEL DEVELOPMENT ######
## utility function to row bind from a list
Rbind <- function(object) do.call(rbind,object)
out <- lapply(1:15, function(i) {
  fitaic <- stpm2(Surv(survival_time,composite_outcome==1)~X_age + e_gfr + me_mcv + X_time_since_HF +
                    me_urea + me_heartrate + me_haemogl + serum_creatinine + me_sbp + me_dbp +
                    me_potassium + me_albumin + me_sodium + imd10_pat ,
                  data = train_data, df=i)
  data.frame(
    i,
    AIC=AIC(fitaic),
    BIC=BIC(fitaic),
    beta=as.numeric(coef(fitaic)[2]),
    se=coef(summary(fitaic))[2,2])
})
out %>% Rbind

# Create a plot of AIC and BIC values
results <- list()
# Iterate through degrees of freedom
for (i in 1:15) {
  fitaic <- stpm2(Surv(survival_time, composite_outcome == 1) ~
                    X_age + e_gfr + me_mcv + X_time_since_HF +
                    me_urea + me_heartrate + me_haemogl + serum_creatinine +
                    me_sbp + me_dbp + me_potassium + me_albumin + me_sodium +
                    imd10_pat  ,
                  data = train_data, df = i)
  
  results[[i]] <- data.frame(
    i,
    AIC = AIC(fitaic),
    BIC = BIC(fitaic),
    beta = as.numeric(coef(fitaic)[2]),
    se = coef(summary(fitaic))[2, 2])
}
# Combine the list of data frames into one data frame
out <- do.call(rbind, results)
ggplot(out, aes(x = i)) +
  geom_line(aes(y = AIC, color = "AIC"), size = 1.5) +
  geom_point(aes(y = AIC, color = "AIC"), size = 4, shape = 16) +
  geom_line(aes(y = BIC, color = "BIC"), size = 1.5, linetype = "dashed") +
  geom_point(aes(y = BIC, color = "BIC"), size = 4, shape = 15) +
  labs(x = "Degrees of Freedom", y = "AIC / BIC") +
  scale_x_continuous(breaks = seq(min(out$i), max(out$i), by = 1)) +
  scale_color_manual(values = c("AIC" = "#1f77b4", "BIC" = "#d62728")) +
  theme_minimal() +
  theme(legend.position = "top",
        legend.background = element_rect(fill = "white"),
        legend.key = element_rect(fill = "white"),
        legend.key.size = unit(1.2, "lines"),
        legend.spacing.x = unit(0.5, "lines"),
        legend.title.align = 0.5,
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        plot.title = element_text(size = 16, hjust = 0.5))


# Fit the Royston-Parmar survival model with time-dependent effects using tvc=list()
fit <- stpm2(Surv(survival_time, composite_outcome == 1) ~
               X_age + e_gfr + me_mcv + X_time_since_HF +me_potassium+
               me_urea + me_heartrate + me_haemogl + serum_creatinine +
               me_sbp + me_dbp + me_albumin + me_sodium +
               factor(imd10_pat),
             data = train_data,
             df = 3,
             tvc = list(X_age = 3)
)

# Display summary of the model
summary(fit)
eform(fit)


###### Baseline hazard plot ######
# Create a copy of the original dataset
zero_covariates <- train_data

# Set the values of all covariates to zero 
zero_covariates$X_age <- 0
zero_covariates$e_gfr <- 0
zero_covariates$me_mcv <- 0
zero_covariates$X_time_since_HF <- 0
zero_covariates$me_potassium  <- 0
zero_covariates$me_urea <- 0
zero_covariates$me_heartrate <- 0
zero_covariates$me_haemogl <- 0
zero_covariates$serum_creatinine <- 0
zero_covariates$me_sbp <- 0
zero_covariates$me_dbp <- 0
zero_covariates$me_albumin <- 0
zero_covariates$me_sodium <- 0
zero_covariates$imd10_pat <- 1
# Filter the dataset to include only observations with events (deaths)
events_data <- zero_covariates[zero_covariates$composite_outcome == 1, ]

# Extract the survival times for events
survival_times <- events_data$survival_time

# Predict the cumulative hazard at the survival times for events
predicted_risk <- predict(fit, newdata = events_data, times = survival_times, type = "cumhaz")

# Create a data frame to store the results
result_table <- data.frame(Time = survival_times, Hazard = predicted_risk)

# Print the result table
print(result_table)

# Load required libraries for plotting
library(ggplot2)

# Create a more visually appealing plot using ggplot2
ggplot(result_table, aes(x = Time, y = Hazard)) +
  geom_step(color = "blue") +  # Adjust line color
  labs(
    x = "Time (years)",
    y = "Cumulative Hazard",
    title = "Cumulative Hazard Plot"
  ) +
  theme_minimal() +  # Use a minimal theme
  theme(plot.title = element_text(hjust = 0.5)) +  # Center the title
  scale_x_continuous(breaks = seq(0, 17, by = 1)) +  # Set x-axis breaks
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1))  # Set y-axis range and intervals




###### Model Performance ######
# Create bins for predicted probabilities
validation_data <- validation_data %>%
  mutate(survival_6mo = ifelse(survival_time <= 0.5 & composite_outcome == 1, 1, 0))
validation_data$predicted_prob <- predict(fit, newdata = validation_data, type = "hazard",times = 0.5)

# Create 10 quantile-based groups ordered by predicted risks
validation_data$predicted_bin <- cut(validation_data$predicted_prob, quantile(validation_data$predicted_prob, probs = seq(0, 1, by = 0.1)))

# Calculate mean predicted probability and observed event rate using aggregate
calibration_summary <- aggregate(cbind(predicted_prob, survival_6mo) ~ predicted_bin, 
                                 data = validation_data, 
                                 FUN = function(x) mean(x, na.rm = TRUE))  # Add na.rm = TRUE to handle NAs

# Rename columns
colnames(calibration_summary) <- c("predicted_bin", "mean_predicted_prob", "observed_event_rate")

# Order the bins by the mean predicted probability within each bin
calibration_summary <- calibration_summary %>%
  arrange(mean_predicted_prob)

# Plot calibration curve
calibration_plot <- ggplot(calibration_summary, aes(x = mean_predicted_prob, y = observed_event_rate)) +
  geom_point() +
  geom_line(aes(group = 1)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(x = "Mean Predicted Probability", y = "Observed Event Rate") +
  theme(plot.title = element_text(hjust = 0.5))
# Set the title using ggtitle and center it
calibration_plot <- calibration_plot + ggtitle("Calibration Curve")
# Print the plot
print(calibration_plot)


# Calculate the calibration slope using individual participant data
calibration_slope <- with(validation_data, cov(predicted_prob, survival_6mo) / var(predicted_prob))
# Print the calibration slope
cat("Calibration Slope:", calibration_slope)
# Set the number of bootstrap iterations
n_iterations <- 1000  # You can adjust this number as needed

# Create an empty vector to store the bootstrap calibration slopes
bootstrap_slopes <- numeric(n_iterations)

# Perform bootstrapping
for (i in 1:n_iterations) {
  # Resample the data with replacement
  bootstrap_sample <- validation_data[sample(nrow(validation_data), replace = TRUE), ]
  
  # Calculate the calibration slope for the bootstrap sample
  bootstrap_slope <- with(bootstrap_sample, cov(predicted_prob, survival_6mo) / var(predicted_prob))
  
  # Store the bootstrap slope
  bootstrap_slopes[i] <- bootstrap_slope
}

# Calculate the confidence interval
confidence_interval <- quantile(bootstrap_slopes, c(0.025, 0.975))

# Print the confidence interval
cat("Calibration Slope:", calibration_slope, "\n")
cat("95% Confidence Interval:", confidence_interval[1], "-", confidence_interval[2])

## Generate ROC curve and calculate C-index
roc_curve <- roc(validation_data$survival_6mo, validation_data$predicted_prob)
c_index <- roc_curve$auc

# Create an empty vector to store AUC values from bootstrapping
auc_values <- numeric(1000)  # You can adjust the number of bootstrap samples (e.g., 1000)
# Number of bootstrap samples
n_bootstrap <- length(auc_values)
# Perform bootstrapping to calculate AUC values
for (i in 1:n_bootstrap) {
  # Sample with replacement from your validation data
  boot_sample <- validation_data[sample(nrow(validation_data), replace = TRUE), ]
  
  # Calculate ROC curve for the bootstrap sample
  roc_curve_boot <- roc(boot_sample$survival_6mo, boot_sample$predicted_prob)
  
  # Store the AUC value
  auc_values[i] <- roc_curve_boot$auc
}

# Calculate the confidence interval
ci_lower <- quantile(auc_values, probs = 0.025)
ci_upper <- quantile(auc_values, probs = 0.975)

# Display the C-index (AUC) and its confidence interval
cat("C-index (AUC):", c_index, "\n")
cat("95% Confidence Interval for C-index (AUC):", ci_lower, "-", ci_upper, "\n")
# Plot ROC curve
plot(roc_curve, print.auc = TRUE, auc.polygon = TRUE,
     auc.polygon.col = "skyblue", max.auc.polygon = TRUE,
     grid = TRUE, main = "ROC Curve")
cat("C-Index:", c_index, "\n")

