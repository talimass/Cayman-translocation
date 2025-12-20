# Load necessary libraries
library(tidyverse)

# Read the CSV file
df <- read.csv("/Users/talimass/Documents/Documents - MacBook Pro/GitHub/Cayman-translocation/microCT/Skeletal thickness 3Dum.csv")

# Reshape and process
data_long <- df %>%
  pivot_longer(-slice, names_to = "Sample", values_to = "Thickness") %>%
  mutate(
    Treatment = sapply(strsplit(Sample, "_"), function(x) x[2]),
    um = slice * 15)  # Convert slice to microns

# Compute mean and SEM per slice per treatment
summary_df <- data_long %>%
  group_by(slice, um, Treatment) %>%
  summarise(
    Mean_Thickness = mean(Thickness, na.rm = TRUE),
    SD = sd(Thickness, na.rm = TRUE),
    n = n(),
    SEM = SD / sqrt(n),
    .groups = "drop"
  )

# Plot WITHOUT error bars
plot_no_error <- ggplot(summary_df, aes(x = um, y = Mean_Thickness, color = Treatment)) +
  geom_point(alpha = 0.6) +
  theme_minimal() +
  labs(title = "Skeletal Thickness by Slice and Treatment (No Error Bars)",
       x = "Depth (µm)",
       y = "Skeletal Thickness ( µm)") +
  scale_color_brewer(palette = "Set2") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 9))
-------------------------------------------------------------------
# flip X and Y plot without error bar

  plot_no_error <- ggplot(summary_df,
                          aes(x = Mean_Thickness, y = um, color = Treatment)) +
  geom_point(alpha = 0.6) +
  
  geom_hline(yintercept = 1500,
             color = "black",
             linetype = "dashed",
             linewidth = 0.8) +
  
  scale_y_reverse() +
  theme_minimal() +
  labs(
    title = NULL,
    x = "Skeletal Thickness (µm)",
    y = "Depth (µm)"
  ) +
  scale_color_brewer(palette = "Set2") +
  theme(
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )



        
-------------------------------------------------------------------        
# Plot WITH SEM error bars
plot_with_sem <- ggplot(summary_df, aes(x = um, y = Mean_Thickness, color = Treatment)) +
  geom_point(alpha = 0.6) +
  geom_errorbar(aes(ymin = Mean_Thickness - SEM, ymax = Mean_Thickness + SEM),
                width = 0.5, alpha = 0.5) +
  theme_minimal() +
  labs(title = "Skeletal Thickness by Slice and Treatment (With SEM)",
       x = "Depth (µm)",
       y = "Skeletal Thickness ( µm)") +
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

#-----------------------------------------------------------------------------
# plot the first 100 slices, X-axis in um

# Reshape and process
data_long <- df %>%
  pivot_longer(-slice, names_to = "Sample", values_to = "Thickness") %>%
  mutate(
    Treatment = sapply(strsplit(Sample, "_"), function(x) x[2]),
    um = slice * 15  # Convert slice to microns
  )

# Filter for first 100 slices only
data_100 <- data_long %>% filter(slice <= 100)

# Compute mean and SEM per slice per treatment
summary_100 <- data_100 %>%
  group_by(slice, um, Treatment) %>%
  summarise(
    Mean_Thickness = mean(Thickness, na.rm = TRUE),
    SD = sd(Thickness, na.rm = TRUE),
    n = n(),
    SEM = SD / sqrt(n),
    .groups = "drop"
  )

# Plot WITHOUT error bars
plot_no_error_100 <- ggplot(summary_100, aes(x = um, y = Mean_Thickness, color = Treatment)) +
  geom_point(alpha = 0.6) +
  theme_minimal() +
  labs(title = "Skeletal Thickness (First 100 Slices)",
       x = "Depth (µm)",
       y = "Skeletal Thickness ( µm)") +
  scale_color_brewer(palette = "Set2") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 9))

# Plot WITH SEM error bars
plot_with_sem_100 <- ggplot(summary_100, aes(x = um, y = Mean_Thickness, color = Treatment)) +
  geom_point(alpha = 0.6) +
  geom_errorbar(aes(ymin = Mean_Thickness - SEM, ymax = Mean_Thickness + SEM),
                width = 10, alpha = 0.5) +
  theme_minimal() +
  labs(title = "Skeletal Thickness with SEM (First 100 Slices)",
       x = "Depth (µm)",
       y = "Skeletal Thickness ( µm)") +
  scale_color_brewer(palette = "Set2") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 9))

# Define output path
output_path <- "/Users/talimass/Documents/Documents - MacBook Pro/GitHub/Cayman-translocation/Ranalysis/output"

# Save plots
ggsave(file.path(output_path, "skeletal_thickness_100slices_no_errorbars.png"),
       plot = plot_no_error_100, width = 8, height = 5, dpi = 300)
ggsave(file.path(output_path, "skeletal_thickness_100slices_no_errorbars.pdf"),
       plot = plot_no_error_100, width = 8, height = 5)

ggsave(file.path(output_path, "skeletal_thickness_100slices_with_sem.png"),
       plot = plot_with_sem_100, width = 8, height = 5, dpi = 300)
ggsave(file.path(output_path, "skeletal_thickness_100slices_with_sem.pdf"),
       plot = plot_with_sem_100, width = 8, height = 5)

# Filter for DD and SS treatments and first 100 slices
data_dd_ss <- data_long %>%
  filter(Treatment %in% c("DD", "SS"), slice <= 100)

# Compute mean and SEM per slice per treatment
summary_dd_ss <- data_dd_ss %>%
  group_by(slice, um, Treatment) %>%
  summarise(
    Mean_Thickness = mean(Thickness, na.rm = TRUE),
    SD = sd(Thickness, na.rm = TRUE),
    n = n(),
    SEM = SD / sqrt(n),
    .groups = "drop"
  )

# Plot WITHOUT error bars
plot_no_error_dd_ss <- ggplot(summary_dd_ss, aes(x = um, y = Mean_Thickness, color = Treatment)) +
  geom_point(alpha = 0.7) +
  theme_minimal() +
  labs(title = "Skeletal Thickness: DD vs SS (No Error Bars)",
       x = "Depth (µm)",
       y = "Skeletal Thickness ( µm)") +
  scale_color_brewer(palette = "Set1") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 9))

# Plot WITH SEM error bars
plot_with_sem_dd_ss <- ggplot(summary_dd_ss, aes(x = um, y = Mean_Thickness, color = Treatment)) +
  geom_point(alpha = 0.7) +
  geom_errorbar(aes(ymin = Mean_Thickness - SEM, ymax = Mean_Thickness + SEM),
                width = 10, alpha = 0.5) +
  theme_minimal() +
  labs(title = "Skeletal Thickness: DD vs SS (With SEM)",
       x = "Depth (µm)",
       y = "Skeletal Thickness ( µm)") +
  scale_color_brewer(palette = "Set1") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 9))

# Define output path
output_path <- "/Users/talimass/Documents/Documents - MacBook Pro/GitHub/Cayman-translocation/Ranalysis/output"

# Save plots
ggsave(file.path(output_path, "DD_SS_100slices_no_errorbars.png"),
       plot = plot_no_error_dd_ss, width = 8, height = 5, dpi = 300)
ggsave(file.path(output_path, "DD_SS_100slices_no_errorbars.pdf"),
       plot = plot_no_error_dd_ss, width = 8, height = 5)

ggsave(file.path(output_path, "DD_SS_100slices_with_sem.png"),
       plot = plot_with_sem_dd_ss, width = 8, height = 5, dpi = 300)
ggsave(file.path(output_path, "DD_SS_100slices_with_sem.pdf"),
       plot = plot_with_sem_dd_ss, width = 8, height = 5)

# Filter for DD and DS treatments and first 100 slices
data_dd_ds <- data_long %>%
  filter(Treatment %in% c("DD", "DS"), slice <= 100)

# Compute mean and SEM per slice per treatment
summary_dd_ds <- data_dd_ds %>%
  group_by(slice, um, Treatment) %>%
  summarise(
    Mean_Thickness = mean(Thickness, na.rm = TRUE),
    SD = sd(Thickness, na.rm = TRUE),
    n = n(),
    SEM = SD / sqrt(n),
    .groups = "drop"
  )

# Plot WITHOUT error bars
plot_no_error_dd_ds <- ggplot(summary_dd_ds, aes(x = um, y = Mean_Thickness, color = Treatment)) +
  geom_point(alpha = 0.7) +
  theme_minimal() +
  labs(title = "Skeletal Thickness: DD vs DS (No Error Bars)",
       x = "Depth (µm)",
       y = "Skeletal Thickness ( µm)") +
  scale_color_brewer(palette = "Set1") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 9))

# Plot WITH SEM error bars
plot_with_sem_dd_ds <- ggplot(summary_dd_ds, aes(x = um, y = Mean_Thickness, color = Treatment)) +
  geom_point(alpha = 0.7) +
  geom_errorbar(aes(ymin = Mean_Thickness - SEM, ymax = Mean_Thickness + SEM),
                width = 10, alpha = 0.5) +
  theme_minimal() +
  labs(title = "Skeletal Thickness: DD vs DS (With SEM)",
       x = "Depth (µm)",
       y = "Skeletal Thickness ( µm)") +
  scale_color_brewer(palette = "Set1") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 9))
# Define output path
output_path <- "/Users/talimass/Documents/Documents - MacBook Pro/GitHub/Cayman-translocation/Ranalysis/output"

# Save plots
ggsave(file.path(output_path, "DD_DS_100slices_no_errorbars.png"),
       plot = plot_no_error_dd_ds, width = 8, height = 5, dpi = 300)
ggsave(file.path(output_path, "DD_DS_100slices_no_errorbars.pdf"),
       plot = plot_no_error_dd_ds, width = 8, height = 5)

ggsave(file.path(output_path, "DD_DS_100slices_with_sem.png"),
       plot = plot_with_sem_dd_ds, width = 8, height = 5, dpi = 300)
ggsave(file.path(output_path, "DD_DS_100slices_with_sem.pdf"),
       plot = plot_with_sem_dd_ds, width = 8, height = 5)

# Filter for SS and SD treatments and first 100 slices
data_ss_sd <- data_long %>%
  filter(Treatment %in% c("SS", "SD"), slice <= 100)

# Compute mean and SEM per slice per treatment
summary_ss_sd <- data_ss_sd %>%
  group_by(slice, um, Treatment) %>%
  summarise(
    Mean_Thickness = mean(Thickness, na.rm = TRUE),
    SD = sd(Thickness, na.rm = TRUE),
    n = n(),
    SEM = SD / sqrt(n),
    .groups = "drop"
  )

# Plot WITHOUT error bars
plot_no_error_ss_sd <- ggplot(summary_ss_sd, aes(x = um, y = Mean_Thickness, color = Treatment)) +
  geom_point(alpha = 0.7) +
  theme_minimal() +
  labs(title = "Skeletal Thickness: SS vs SD (No Error Bars)",
       x = "Depth (µm)",
       y = "Skeletal Thickness ( µm)") +
  scale_color_brewer(palette = "Set1") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 9))

# Plot WITH SEM error bars
plot_with_sem_ss_sd <- ggplot(summary_ss_sd, aes(x = um, y = Mean_Thickness, color = Treatment)) +
  geom_point(alpha = 0.7) +
  geom_errorbar(aes(ymin = Mean_Thickness - SEM, ymax = Mean_Thickness + SEM),
                width = 10, alpha = 0.5) +
  theme_minimal() +
  labs(title = "Skeletal Thickness: SS vs SD (With SEM)",
       x = "Depth (µm)",
       y = "Skeletal Thickness ( µm)") +
  scale_color_brewer(palette = "Set1") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 9))

# Define output path
output_path <- "/Users/talimass/Documents/Documents - MacBook Pro/GitHub/Cayman-translocation/Ranalysis/output"

# Save plots
ggsave(file.path(output_path, "SS_SD_100slices_no_errorbars.png"),
       plot = plot_no_error_ss_sd, width = 8, height = 5, dpi = 300)
ggsave(file.path(output_path, "SS_SD_100slices_no_errorbars.pdf"),
       plot = plot_no_error_ss_sd, width = 8, height = 5)

ggsave(file.path(output_path, "SS_SD_100slices_with_sem.png"),
       plot = plot_with_sem_ss_sd, width = 8, height = 5, dpi = 300)
ggsave(file.path(output_path, "SS_SD_100slices_with_sem.pdf"),
       plot = plot_with_sem_ss_sd, width = 8, height = 5)
