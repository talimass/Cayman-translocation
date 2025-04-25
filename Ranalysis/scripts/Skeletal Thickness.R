# Load necessary libraries
library(tidyverse)

# Read the CSV file
df <- read.csv("/Users/talimass/Documents/Documents - MacBook Pro/GitHub/Cayman-translocation/microCT/Skeletal thickness 3Dum.csv")

# Reshape and process
data_long <- df %>%
  pivot_longer(-slice, names_to = "Sample", values_to = "Thickness") %>%
  mutate(
    Treatment = sapply(strsplit(Sample, "_"), function(x) x[2])
  )

# Compute mean and SEM per slice per treatment
summary_df <- data_long %>%
  group_by(slice, Treatment) %>%
  summarise(
    Mean_Thickness = mean(Thickness, na.rm = TRUE),
    SD = sd(Thickness, na.rm = TRUE),
    n = n(),
    SEM = SD / sqrt(n),
    .groups = "drop"
  )

# Plot WITHOUT error bars
plot_no_error <- ggplot(summary_df, aes(x = slice, y = Mean_Thickness, color = Treatment)) +
  geom_point(alpha = 0.6) +
  theme_minimal() +
  labs(title = "Skeletal Thickness by Slice and Treatment (No Error Bars)",
       x = "Slice",
       y = "Skeletal Thickness (3 µm)") +
  scale_color_brewer(palette = "Set2") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 9))

# Plot WITH SEM error bars
plot_with_sem <- ggplot(summary_df, aes(x = slice, y = Mean_Thickness, color = Treatment)) +
  geom_point(alpha = 0.6) +
  geom_errorbar(aes(ymin = Mean_Thickness - SEM, ymax = Mean_Thickness + SEM),
                width = 0.5, alpha = 0.5) +
  theme_minimal() +
  labs(title = "Skeletal Thickness by Slice and Treatment (With SEM)",
       x = "Slice",
       y = "Skeletal Thickness (3 µm)") +
  scale_color_brewer(palette = "Set2") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 9))

# Define output directory
output_path <- "/Users/talimass/Documents/Documents - MacBook Pro/GitHub/Cayman-translocation/Ranalysis/output"

# Save plot WITHOUT error bars
ggsave(filename = file.path(output_path, "skeletal_thickness_no_errorbars.png"),
       plot = plot_no_error, width = 8, height = 5, dpi = 300)
ggsave(filename = file.path(output_path, "skeletal_thickness_no_errorbars.pdf"),
       plot = plot_no_error, width = 8, height = 5)

# Save plot WITH SEM error bars
ggsave(filename = file.path(output_path, "skeletal_thickness_with_sem.png"),
       plot = plot_with_sem, width = 8, height = 5, dpi = 300)
ggsave(filename = file.path(output_path, "skeletal_thickness_with_sem.pdf"),
       plot = plot_with_sem, width = 8, height = 5)

