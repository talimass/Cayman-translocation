library(readxl)
library(dplyr)
library(lubridate)
library(ggplot2)
library(scales)

Sys.setlocale("LC_TIME", "en_US.UTF-8")

setwd("/home/gospozha/haifa/cayman/P.astreoides_physiology/loggers/")

## PAR logger

# Helper function
process_par_file <- function(filename, source_label, linetype_value) {
  read_excel(filename, skip = 5) %>%
    slice(-1) %>%
    mutate(
      datetime = ymd_hms(`Colombia Time`),
      PAR = as.numeric(PAR),
      date = as_date(datetime)
    ) %>%
    filter(!is.na(date), !is.na(PAR)) %>%
    group_by(date) %>%
    summarise(max_PAR = max(PAR, na.rm = TRUE)) %>%
    mutate(source = source_label, line_type = linetype_value)
}

# Process all four files
df1 <- process_par_file("CC_10m_PAR.xlsx", "Coral City, 10m", 'solid')
df2 <- process_par_file("CC_40m_PAR.xlsx", "Coral City, 40m", 'solid')
df3 <- process_par_file("MF_10m_PAR.xlsx", "Martha’s Finyard, 10m", 'dashed')
df4 <- process_par_file("MF_40m_PAR.xlsx", "Martha’s Finyard, 40m", 'dashed')

# Combine them
all_data <- bind_rows(df1, df2, df3, df4)

# Plot
PAR <- ggplot(all_data, aes(x = date, y = max_PAR, group = source)) +
  geom_line(aes(color = source, linetype = source), size = 0.3) +
  scale_linetype_manual(values = setNames(all_data$line_type, all_data$source)) +
  scale_x_date(date_labels = "%m/%y", date_breaks = "1 month") +
  labs(title = "Daily Maximum PAR",
       x = "Date (MM/YY)",
       y = "PAR (µmol/(s·m²))",
       color = "Location",
       linetype = "Location") +
  theme_minimal(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

PAR

# Save
ggsave("PAR.jpg", PAR, width = 6, height = 4)

## Temperature

# Function to process temperature data
process_temp_file <- function(filename, source_label, linetype_value) {
  read.csv(filename, skip = 3) %>%
    mutate(
      datetime = mdy_hms(.[[2]]),
      Temp = as.numeric(.[[3]]),
      date = as_date(datetime)
    ) %>%
    filter(!is.na(date), !is.na(Temp)) %>%
    group_by(date) %>%
    summarise(mean_temp = mean(Temp, na.rm = TRUE)) %>%
    mutate(source = source_label, line_type = linetype_value)
}

# Replace with your actual file names and location labels
df1 <- process_temp_file("CoralCity_10m_temp.csv",  "Coral City, 10m", 'solid')
df2 <- process_temp_file("CoralCity_45m_temp.csv",  "Coral City, 45m", 'solid')
df3 <- process_temp_file("Marthas_10m_temp.csv", "Martha’s Finyard, 10m", 'dashed')
df4 <- process_temp_file("Marthas40m_temp.csv", "Martha’s Finyard, 40m", 'dashed')

# Combine all
all_temp <- bind_rows(df1, df2, df3, df4)

# Plot
temp <- ggplot(all_temp, aes(x = date, y = mean_temp, group = source)) +
  geom_line(aes(color = source, linetype = source), size = 0.3) +
  scale_linetype_manual(values = setNames(all_temp$line_type, all_temp$source)) +
  scale_x_date(date_labels = "%m/%y", date_breaks = "1 month") +
  labs(title = "Daily Average Temperature",
       x = "Date (MM/YY)",
       y = "Temperature (°C)",
       color = "Location",
       linetype = "Location") +
  theme_minimal(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

temp
# Save
ggsave("temp.jpg", temp, width = 6, height = 4)
