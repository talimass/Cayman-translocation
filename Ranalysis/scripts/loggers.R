library(readxl)
library(dplyr)
library(lubridate)
library(ggplot2)
library(scales)
library(tidyr)
library(patchwork)
library(ggbreak)
library(emmeans)
library(ggpattern)

Sys.setlocale("LC_TIME", "en_US.UTF-8")

setwd("/home/gospozha/haifa/cayman/P.astreoides_physiology/loggers/")


#### Temperature ####

# Function to process temperature data
process_temp_file <- function(filename, source_label, linetype_value, color_value) {
  read.csv(filename, skip = 3) %>%
    mutate(
      datetime = mdy_hms(.[[2]]),
      Temp = as.numeric(.[[3]]),
      date = as_date(datetime)
    ) %>%
    filter(!is.na(date), !is.na(Temp)) %>%
    group_by(date) %>%
    summarise(mean_temp = mean(Temp, na.rm = TRUE)) %>%
    mutate(Location = source_label, line_type = linetype_value, color = color_value)
}

# Replace with your actual file names and location labels
df1 <- process_temp_file("CoralCity_10m_temp.csv",  "Coral City, 10m",'solid', '#f75f55')
df2 <- process_temp_file("CoralCity_45m_temp.csv",  "Coral City, 40m", 'solid', "#00A9FF")
df3 <- process_temp_file("Marthas_10m_temp.csv", "Martha’s Finyard, 10m", 'twodash', "#f75f55")
df4 <- process_temp_file("Marthas40m_temp.csv", "Martha’s Finyard, 40m", 'twodash', "#00A9FF")

# Combine all
all_temp <- bind_rows(df1, df2, df3, df4)

# Plot
temp <- ggplot(all_temp, aes(x = date, y = mean_temp, group = Location)) +
  geom_line(aes(color = Location, linetype = Location), size = 0.3) +
  scale_color_manual(values = setNames(all_temp$color, all_temp$Location)) +
  scale_linetype_manual(values = setNames(all_temp$line_type, all_temp$Location)) +
  scale_x_date(date_labels = "%m/%y", date_breaks = "1 month") +
  labs(title = "Daily average temperature",
       x = "Date (MM/YY)",
       y = "Temperature (°C)",
       color = "Location",
       linetype = "Location") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

temp
# Save
#ggsave("temp.jpg", temp, width = 6, height = 4)

# Mean temp for July
# Extract Location and Depth from the combined column
all_data <- all_temp %>%
  separate(Location, into = c("Location", "Depth"), sep = ", ") %>%
  mutate(Depth = trimws(Depth))

# Filter only July 2023
temp_july <- all_data %>%
  filter(month(date) == 7 & year(date) == 2023)

# Compute mean PAR per Location and Depth
mean_temp <- temp_july %>%
  group_by(Location, Depth) %>%
  summarise(mean_temp = mean(mean_temp, na.rm = TRUE))

# View result
print(mean_temp)

# A tibble: 4 × 3
# Groups:   Location [2]
#Location         Depth mean_temp
#<chr>            <chr>     <dbl>
# 1 Coral City       10m        30.8
# 2 Coral City       45m        29.6
# 3 Martha’s Finyard 10m        30.2
# 4 Martha’s Finyard 40m        29.7


#### additional PAR and temp for 2022 ####

#### Temperature ####

# Function to process temperature data
process_temp_file2 <- function(filename, source_label, linetype_value, color_value) {
  read.csv(filename, skip = 3) %>%
    mutate(
      datetime = mdy_hm(.[[2]]),
      Temp = as.numeric(.[[3]]),
      date = as_date(datetime)
    ) %>%
    filter(!is.na(date), !is.na(Temp)) %>%
    group_by(date) %>%
    summarise(mean_temp = mean(Temp, na.rm = TRUE)) %>%
    mutate(Location = source_label, line_type = linetype_value, color = color_value)
}

# Replace with your actual file names and location labels
df5 <- process_temp_file2("Marthas_12m_070122-092122.csv", "Martha’s Finyard, 10m", 'twodash', "#f75f55")
df6 <- process_temp_file2("Marthas_36m_070122-092222.csv", "Martha’s Finyard, 40m", 'twodash', "#00A9FF")

# Combine all
all_temp2 <- bind_rows(df5, df6)

# Plot
temp <- ggplot(all_temp2, aes(x = date, y = mean_temp, group = Location)) +
  geom_line(aes(color = Location, linetype = Location), size = 0.3) +
  scale_color_manual(values = setNames(all_temp2$color, all_temp2$Location)) +
  scale_linetype_manual(values = setNames(all_temp2$line_type, all_temp2$Location)) +
  scale_x_date(date_labels = "%m/%y", date_breaks = "1 month") +
  labs(title = "Daily average temperature",
       x = "Date (MM/YY)",
       y = "Temperature (°C)",
       color = "Location",
       linetype = "Location") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

temp
# Save
#ggsave("temp.jpg", temp, width = 6, height = 4)

# Mean temp for September

# Extract Location and Depth from the combined column
all_temp2 <- all_temp2 %>%
  separate(Location, into = c("Location", "Depth"), sep = ", ") %>%
  mutate(Depth = trimws(Depth))

# Filter only July 2023
temp_sep <- all_temp2 %>%
  filter(month(date) == 9 & year(date) == 2022)

# Compute mean PAR per Location and Depth
mean_temp <- temp_sep%>%
  group_by(Location, Depth) %>%
  summarise(mean_temp = mean(mean_temp, na.rm = TRUE))

# View result
print(mean_temp)
# A tibble: 2 × 3
# Groups:   Location [1]
#Location         Depth mean_temp
#<chr>            <chr>     <dbl>
#  1 Martha’s Finyard 10m        30.1
#2 Martha’s Finyard 40m        29.9


#### midday average par ####

process_par_file <- function(filename, source_label, linetype_value, color_value,
                             start_h = 12, end_h = 14) {
  read_excel(filename, skip = 5) %>%
    slice(-1) %>%
    mutate(
      datetime = ymd_hms(`Colombia Time`),
      PAR = as.numeric(PAR),
      date = as_date(datetime),
      hour = hour(datetime)
    ) %>%
    filter(!is.na(date), !is.na(PAR)) %>%
    # keep only the midday window [start_h, end_h)
    filter(hour >= start_h, hour < end_h) %>%
    group_by(date) %>%
    summarise(midday_mean_PAR = mean(PAR, na.rm = TRUE), .groups = "drop") %>%
    mutate(Location = source_label, line_type = linetype_value, color = color_value)
}


# Process all four files (12:00–13:00)
df1 <- process_par_file("CC_10m_PAR.xlsx", "Coral City, 10m", 'solid', '#f75f55', 12, 14)
df2 <- process_par_file("CC_40m_PAR.xlsx", "Coral City, 40m", 'solid', "#00A9FF", 12, 14)
df3 <- process_par_file("MF_10m_PAR.xlsx", "Martha’s Finyard, 10m", 'twodash', "#f75f55", 12, 14)
df4 <- process_par_file("MF_40m_PAR.xlsx", "Martha’s Finyard, 40m", 'twodash', "#00A9FF", 12, 14)

# Process all four files
df5 <- process_par_file("PAR_12m_July22_Sept22.xlsx", "Martha’s Finyard, 10m", 'twodash', "#f75f55", 12, 13)
df6 <- process_par_file("PAR_36m_july22_Sept22.xlsx", "Martha’s Finyard, 40m", 'twodash', "#00A9FF", 12, 13)

# Combine them
all_data2 <- bind_rows( df5, df6)

all_data <- bind_rows(df1, df2, df3, df4)

PAR_midday <- ggplot(all_data, aes(x = date, y = midday_mean_PAR, color = Location, linetype = Location)) +
  geom_line(size = 0.3) +
  scale_color_manual(values = setNames(all_data$color, all_data$Location)) +
  scale_linetype_manual(values = setNames(all_data$line_type, all_data$Location)) +
  scale_x_date(date_labels = "%m/%y", date_breaks = "1 month") +
  labs(title = "Midday average PAR",
       x = "Date (MM/YY)",
       y = "PAR (µmol/(s·m²))") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

PAR_midday

# Mean midday PAR

# Extract Location and Depth from the combined column
all_data <- all_data %>%
  separate(Location, into = c("Location", "Depth"), sep = ", ") %>%
  mutate(Depth = trimws(Depth))

# Filter only July 2023
par_july <- all_data %>%
  filter(month(date) == 7 & year(date) == 2023)

# Compute mean PAR per Location and Depth
mean_par_july <- par_july %>%
  group_by(Location, Depth) %>%
  summarise(mean_PAR = mean(midday_mean_PAR, na.rm = TRUE))

# View result
print(mean_par_july)
# A tibble: 4 × 3
# Groups:   Location [2]
#Location         Depth mean_PAR
#<chr>            <chr>    <dbl>
#1 Coral City       10m      310 
#2 Coral City       40m      45.4
#3 Martha’s Finyard 10m      430 
#4 Martha’s Finyard 40m      137 

# Filter only nov 22
all_data2 <- all_data2 %>%
  separate(Location, into = c("Location", "Depth"), sep = ", ") %>%
  mutate(Depth = trimws(Depth))

par_nov <- all_data2 %>%
  filter(month(date) == 9 & year(date) == 2022)

# Compute mean PAR per Location and Depth
mean_par_nov <- par_nov %>%
  group_by(Location, Depth) %>%
  summarise(mean_PAR = mean(midday_mean_PAR, na.rm = TRUE))

# View result
print(mean_par_nov)
# Location         Depth mean_PAR
# <chr>            <chr>    <dbl>
#   1 Martha’s Finyard 10m       455.
# 2 Martha’s Finyard 40m       190.

#### joining data from winter and summer ####
all_data
all_data2
all_data3 <- bind_rows( all_data, all_data2)


gap_start <- as.Date("2022-09-22")
gap_end   <- as.Date("2023-03-13")

all_data3_gapfixed <- all_data3 %>%
  mutate(midday_mean_PAR = ifelse(date >= gap_start & date <= gap_end, NA, midday_mean_PAR))

PAR <- ggplot(all_data3_gapfixed, aes(x = date, y = midday_mean_PAR,
                             color = Location, linetype = Location)) +
  geom_line(size = 0.3) +
  # vertical lines to mark break
  geom_vline(xintercept = as.numeric(gap_start), linetype = "dashed", color = "gray") +
  geom_vline(xintercept = as.numeric(gap_end), linetype = "dashed", color = "gray") +
  scale_color_manual(values = setNames(all_data3$color, all_data3$Location)) +
  scale_linetype_manual(values = setNames(all_data3$line_type, all_data3$Location)) +
  scale_x_date(
    date_labels = "%m/%y", 
    date_breaks = "1 month", 
    expand = c(0,0)
  ) +
  scale_x_break(c(gap_start, gap_end)) +   # axis break
  labs(title = "Midday average PAR",
       x = "Date (MM/YY)",
       y = "PAR (µmol/(s·m²))") +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.x.top = element_blank(),
    axis.ticks.x.top = element_blank(),
    axis.title.x.top = element_blank()
  )
PAR


all_temp
all_temp2
all_temp3 <- bind_rows( all_temp, all_temp2)

temp <- ggplot(all_temp3, aes(x = date, y = mean_temp,
                             color = Location, linetype = Location)) +
  geom_line(size = 0.3) +
  # vertical lines to mark break
  geom_vline(xintercept = as.numeric(gap_start), linetype = "dashed", color = "gray") +
  geom_vline(xintercept = as.numeric(gap_end), linetype = "dashed", color = "gray") +
  scale_color_manual(values = setNames(all_temp3$color, all_temp3$Location)) +
  scale_linetype_manual(values = setNames(all_temp3$line_type, all_temp3$Location)) +
  scale_x_date(
    date_labels = "%m/%y", 
    date_breaks = "1 month", 
    expand = c(0,0)
  ) +
  scale_x_break(c(gap_start, gap_end)) +   # axis break
  labs(title = "Daily average temperature",
       x = "Date (MM/YY)",
       y = "Temperature (°C)",
       color = "Location",
       linetype = "Location") +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.x.top = element_blank(),
    axis.ticks.x.top = element_blank(),
    axis.title.x.top = element_blank()
  )
temp


#### combined plots  ####

combined <- (PAR / temp) 
combined <- combined + plot_layout(widths = c(1.4, 0.6)) +
 # plot_annotation(tag_levels = "B") +
  plot_layout(heights = c(1, 1)) +
  theme(plot.margin = margin(10, 10, 10, 10))
combined

ggsave("~/haifa/cayman/P.astreoides_physiology/github3/5/combined_loggers_midday.pdf", combined, width = 10, height = 8)

#### tables for daily average par ####
all_temp3
all_data3_gapfixed


# PAR 
par_summary <- all_data3_gapfixed %>%
  separate(Location, into = c("site", "depth"), sep = ", ") %>%
  mutate(depth = as.numeric(gsub("m", "", depth)),
         month = format(date, "%Y-%m")) %>%
  group_by(site, depth) %>%
  summarise(
    max_PAR_val  = max(midday_mean_PAR, na.rm = TRUE),
    month_max_PAR = month[which.max(midday_mean_PAR)],   # which.max on the original vector
    min_PAR_val  = min(midday_mean_PAR, na.rm = TRUE),
    month_min_PAR = month[which.min(midday_mean_PAR)],
    mean_PAR     = mean(midday_mean_PAR, na.rm = TRUE),
    sd_PAR = sd(midday_mean_PAR, na.rm = TRUE),
    n = sum(!is.na(midday_mean_PAR)),
    se_PAR = sd_PAR / sqrt(n),
    ci_lower = mean_PAR - 1.96 * se_PAR,
    ci_upper = mean_PAR + 1.96 * se_PAR,
    .groups = "drop"
  ) %>%
  rename(max_PAR = max_PAR_val, min_PAR = min_PAR_val)

# 2. Calculate monthly delta PAR
par_monthly_delta <- all_data3_gapfixed %>%
  separate(Location, into = c("site", "depth"), sep = ", ") %>%
  mutate(
    depth = as.numeric(gsub("m", "", depth)),
    date = as.Date(date),
    month = format(date, "%Y-%m")
  ) %>%
  group_by(site, depth, month) %>%
  summarise(
    monthly_delta_PAR = max(midday_mean_PAR, na.rm = TRUE) - min(midday_mean_PAR, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(site, depth) %>%
  summarise(
    mean_delta_PAR = mean(monthly_delta_PAR, na.rm = TRUE),
    sd_delta_PAR = sd(monthly_delta_PAR, na.rm = TRUE),
    n_delta = n(),
    se_delta_PAR = sd_delta_PAR / sqrt(n_delta),
    .groups = "drop"
  )

# 3. Join the two tables
par_overall <- left_join(par_summary, par_monthly_delta, by = c("site", "depth"))

write.csv2(par_overall, "~/haifa/cayman/P.astreoides_physiology/github3/5/par.midday.csv")
#  TEMPERATURE DATA 
temp_summary <- all_temp3 %>%
  separate(Location, into = c("site", "depth"), sep = ", ") %>%
  mutate(depth = as.numeric(gsub("m", "", depth)),
         month = format(date, "%Y-%m")) %>%
  group_by(site, depth) %>%
  summarise(
    min_temp = min(mean_temp, na.rm = TRUE),
    month_min_temp = month[which.min(mean_temp)],
    max_temp = max(mean_temp, na.rm = TRUE),
    month_max_temp = month[which.max(mean_temp)],
    mean_temp_val = mean(mean_temp, na.rm = TRUE),
    sd_temp = sd(mean_temp, na.rm = TRUE),
    n = sum(!is.na(mean_temp)),
    se_temp = sd_temp / sqrt(n),
    .groups = "drop"
  )

# 2. Calculate monthly delta temp
temp_monthly_delta <- all_temp3 %>%
  separate(Location, into = c("site", "depth"), sep = ", ") %>%
  mutate(depth = as.numeric(gsub("m", "", depth)),
         date = as.Date(date),
         month = format(date, "%Y-%m")
  ) %>%
  group_by(site, depth, month) %>%
  summarise(
    monthly_delta_temp = max(mean_temp, na.rm = TRUE) - min(mean_temp, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(site, depth) %>%
  summarise(
    mean_delta_temp = mean(monthly_delta_temp, na.rm = TRUE),
    sd_delta_temp = sd(monthly_delta_temp, na.rm = TRUE),
    n_delta = n(),
    se_delta_temp = sd_delta_temp / sqrt(n_delta),
    .groups = "drop"
  )

temp_overall <- left_join(temp_summary, temp_monthly_delta, by = c("site", "depth"))
write.csv2(temp_overall, "~/haifa/cayman/P.astreoides_physiology/github3/5/temp.csv")
