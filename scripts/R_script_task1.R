

# Install methylKit if not already installed
if (!requireNamespace("methylKit", quietly = TRUE)) {
  install.packages("methylKit")
}

BiocManager::install("methylKit")
library(methylKit)
install.packages("reshape2")  # Run this if you haven't installed it yet
library(tidyr)

data <- read.csv("PupilBioTest_PMP_revA.csv")


# Check for missing values
colSums(is.na(data))

# Remove rows with missing values
data <- na.omit(data)

library(dplyr)
grouped_data <- data %>% group_by(Tissue)

head(data)
str(data)

colnames(data)

median_coverage <- grouped_data %>% 
  summarise(across(starts_with("X"), ~ median(.x, na.rm = TRUE), .names = "median_{col}"))

cv_coverage <- grouped_data %>% 
  summarise(across(starts_with("X"), ~ (sd(.x, na.rm = TRUE) / mean(.x, na.rm = TRUE)) * 100, .names = "cv_{col}"))

ggplot(median_coverage, aes(x = Tissue, y = median_X.000)) +
  geom_boxplot() +
  labs(title = "Median Coverage for X.000 by Tissue", y = "Median Coverage", x = "Tissue") +
  theme_minimal()

ggplot(median_coverage, aes(x = Tissue, y = median_X.001)) +
  geom_boxplot() +
  labs(title = "Median Coverage for X.000 by Tissue", y = "Median Coverage", x = "Tissue") +
  theme_minimal()

ggplot(median_coverage, aes(x = Tissue, y = median_X.010)) +
  geom_boxplot() +
  labs(title = "Median Coverage for X.000 by Tissue", y = "Median Coverage", x = "Tissue") +
  theme_minimal()


ggplot(median_coverage, aes(x = Tissue, y = median_X.111)) +
  geom_boxplot() +
  labs(title = "Median Coverage for X.000 by Tissue", y = "Median Coverage", x = "Tissue") +
  theme_minimal()

# Reshape data for bar plot (median)
median_long <- median_coverage %>% 
  pivot_longer(cols = starts_with("median"), names_to = "Methylation_Pattern", values_to = "Median_Coverage")

ggplot(median_long, aes(x = Tissue, y = Median_Coverage, fill = Methylation_Pattern)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Median Coverage by Tissue and Methylation Pattern", y = "Median Coverage", x = "Tissue") +
  theme_minimal()




# Combine median and CV data
combined_data <- median_coverage %>%
  left_join(cv_coverage, by = "Tissue")

# Example for scatter plot of median vs. CV for X.000
ggplot(combined_data, aes(x = starts_with("median") , y = starts_with("cv") , color = Tissue)) +
  geom_point(size = 3) +
  labs(title = "Median Coverage vs. CV for All patterns", x = "Median Coverage", y = "CV (%)") +
  theme_minimal()

#heatmap 
library(reshape2)

# Reshape data for heatmap (using median coverage as an example)
heatmap_data <- median_coverage %>% 
  pivot_longer(cols = starts_with("median"), names_to = "Methylation_Pattern", values_to = "Median_Coverage")

ggplot(heatmap_data, aes(x = Methylation_Pattern, y = Tissue, fill = Median_Coverage)) +
  geom_tile() +
  labs(title = "Heatmap of Median Coverage by Tissue and Methylation Pattern", x = "Methylation Pattern", y = "Tissue") +
  theme_minimal() +
  scale_fill_gradient(low = "blue", high = "red")

########################################################################################################################
#Biomarker identification 

# Load necessary libraries
library(dplyr)

# Perform t-test for each methylation pattern column, grouping by Tissue
t_test_results <- data %>%
  group_by(Tissue) %>%
  summarise(across(starts_with("X"), 
                   ~ broom::tidy(t.test(. ~ Tissue, data = data))$p.value[1], 
                   .names = "p_{col}"),
            .groups = "drop")

# View the results
print(t_test_results)

summary(data)

data %>%
  group_by(Tissue) %>%
  summarise(across(starts_with("X"), var))

#################################################################################################################
# Calculate Total Reads as the sum of all X.* columns
data <- data %>%
  mutate(Total_Reads = rowSums(select(., starts_with("X"))))

# Calculate VRF for each PMP
data <- data %>%
  mutate(across(starts_with("X"), ~ . / Total_Reads, .names = "vrf_{col}"))

# Calculate mean VRF by Tissue
mean_vrf <- data %>%
  group_by(Tissue) %>%
  summarise(across(starts_with("vrf_"), mean, na.rm = TRUE))

library(ggplot2)

ggplot(mean_vrf, aes(x = Tissue, y = vrf_X.000)) +
  geom_bar(stat = "identity") +
  labs(title = "Mean Variant Read Fraction (VRF) by Tissue", 
       y = "Mean VRF", 
       x = "Tissue") +
  theme_minimal()



#################################################################################################################

# Select top 10 PMPs based on specificity
top_pmps <- data %>% 
  arrange(desc(specificity_score)) %>% 
  head(10)

# Calculate specificity for individual CpG sites
cpg_specificity <- data %>% 
  summarise(across(starts_with("X"), ~ mean(. == Tissue)))

# Statistical comparison
t_test_results <- t.test(top_pmps$specificity_score, cpg_specificity)


#################################################################################################################

# Example code to calculate specificity (simplified approach)
calculate_specificity_function <- function(row) {
  # Example: Calculate specificity for each row
  # Modify this function based on your actual data structure and calculation needs
  TN <- sum(row == "Tissue1" & actual == "Tissue1")  # True Negatives
  FP <- sum(row == "Tissue2" & actual == "Tissue1")  # False Positives
  specificity <- TN / (TN + FP)
  return(specificity)
}
# Assuming you have a binary outcome (presence/absence) for PMPs in each tissue
data <- data %>%
  mutate(specificity_score = calculate_specificity_function(.))

# Replace calculate_specificity_function with your actual specificity calculation logic



# Select top 10 PMPs based on specificity score
top_pmps <- data %>% 
  arrange(desc(specificity_score)) %>% 
  head(10)

# View the top 10 PMPs
print(top_pmps)

###############################################################################################################################################
#15.01.25 Agnik Sir approach

library(dplyr)

# Step 1: Calculate mean for each pattern by tissue
mean_by_tissue <- data %>%
  group_by(Tissue) %>%
  summarise(across(starts_with("X"), mean, na.rm = TRUE))

# Step 2: Calculate the proportion of non-zero entries for each pattern
non_zero_proportion <- data %>%
  summarise(across(starts_with("X"), ~ mean(. > 0)))

# Step 3: Combine means and non-zero proportions
mean_diff <- mean_by_tissue %>%
  pivot_wider(names_from = Tissue, values_from = starts_with("X")) %>%
  mutate(across(starts_with("X"), ~ abs(.[["cfDNA"]] - .[["Islet"]]), .names = "mean_diff_{col}"))

# Step 4: Weight mean difference by non-zero proportion
weighted_mean_diff <- mean_diff %>%
  mutate(across(starts_with("mean_diff_"), ~ . * non_zero_proportion[[sub("mean_diff_", "", cur_column())]], .names = "weighted_{col}"))

# Step 5: Rank by weighted mean difference
top_pmps <- weighted_mean_diff %>%
  arrange(desc(across(starts_with("weighted_")))) %>%
  slice(1:10)

print(top_pmps)


#########################################################################################33
library(dplyr)

mean_by_tissue <- data %>%
  group_by(Tissue) %>%
  summarise(across(starts_with("X"), mean, na.rm = TRUE), .groups = "drop")


cfDNA_means <- filter(mean_by_tissue, Tissue == "cfDNA")
islet_means <- filter(mean_by_tissue, Tissue == "Islet")

mean_diff <- abs(cfDNA_means[,-1] - islet_means[,-1])  # Removing the Tissue column
mean_diff_df <- data.frame(
  Methylation_Pattern = colnames(mean_diff),
  Mean_Difference = as.numeric(mean_diff)
)


# Step 4: Weight mean difference by non-zero proportion
weighted_mean_diff <- mean_diff %>%
  mutate(across(starts_with("mean_diff_"), ~ . * non_zero_proportion[[sub("mean_diff_", "", cur_column())]], .names = "weighted_{col}"))

# Step 5: Rank by weighted mean difference
top_pmps <- weighted_mean_diff %>%
  arrange(desc(across(starts_with("weighted_")))) %>%
  slice(1:10)

print(top_pmps)




top_pmps <- mean_diff_df %>%
  arrange(desc(Mean_Difference)) %>%
  slice(1:10)

print(top_pmps)


###############################################################################################################

# Calculate Total Reads if not available
data <- data %>%
  rowwise() %>%
  mutate(Total_Reads = sum(c_across(starts_with("X"))))


# Calculate VRF for each PMP
vrf_data <- data %>%
  mutate(across(starts_with("X"), ~ . / Total_Reads, .names = "vrf_{col}"))








