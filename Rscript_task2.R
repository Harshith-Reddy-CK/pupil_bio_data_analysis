
# Load necessary libraries
library(dplyr)
library(ggplot2)

# Read the data
data <- read.table("normal_quality_depth.tsv", header = FALSE, col.names = c("QUAL", "DP"))

# Calculate median quality score
median_quality <- median(data$QUAL)
print(paste("Median Quality Score:", median_quality))

# Plot the distribution of quality scores
ggplot(data, aes(x = QUAL)) +
  geom_histogram(binwidth = 1, fill = "blue", color = "black") +
  labs(title = "Distribution of Quality Scores", x = "Quality Score", y = "Frequency") +
  theme_minimal()

# Total number of reads (Assuming you have extracted this using samtools flagstat)
total_reads <- 5116483   # Replace this with the actual total read count

# Calculate RPM for each variant
rmp_data <- data %>%
  mutate(RPM = (DP / total_reads) * 1e6)

# Determine RPM threshold based on a given percentile (e.g., 95th percentile)
rpm_threshold <- quantile(rmp_data$RPM, 0.95)
print(paste("RPM Threshold (95th percentile):", rpm_threshold))

