# Load necessary libraries
library(ggplot2)
library(dplyr)

# Read the CSV
data <- read.csv("/Users/talimass/Documents/Documents - MacBook Pro/GitHub/Cayman-translocation/microCT/Porosity_python_2025.csv")
names(data) <- trimws(names(data))  # Clean column names

# Clean column names (remove trailing spaces)
names(data) <- trimws(names(data))

# Extract treatment from the 'sample' column
data$treatment <- sub(".*_", "", data$sample)

# Calculate mean and SEM for each treatment
summary_data <- data %>%
  group_by(treatment) %>%
  summarise(
    mean_porosity = mean(total_porosity_perc, na.rm = TRUE),
    sem_porosity = sd(total_porosity_perc, na.rm = TRUE) / sqrt(n())
  )

# Plot
porosity_plot <- ggplot(summary_data, aes(x = treatment, y = mean_porosity)) +
  geom_bar(stat = "identity", fill = "skyblue", width = 0.6) +
  geom_errorbar(aes(ymin = mean_porosity - sem_porosity, ymax = mean_porosity + sem_porosity),
                width = 0.2) +
  labs(title = "Average Porosity by Treatment",
       x = "Treatment",
       y = "Porosity% (mean ± SEM)") +
  theme_minimal()

# Define output path
output_path <- "~/Documents/Documents - MacBook Pro/GitHub/Cayman-translocation/Ranalysis/output"

# Save plot
ggsave(filename = file.path(output_path, "Average_Porosity_by_Treatment.png"),
       porosity_plot, width = 8, height = 6)
ggsave(filename = file.path(output_path, "Average_Porosity_by_Treatment.pdf"),
       porosity_plot, width = 8, height = 6)

# Run ANOVA
anova_result <- aov(total_porosity_perc ~ treatment, data = data)

# Define path to save the ANOVA results
anova_output_path <- file.path(output_path, "anova_results.txt")

# Save ANOVA summary to text file
capture.output(summary(anova_result), file = anova_output_path)
