# Load required libraries
library(ggplot2)
library(tidyr)
library(dplyr)
library(lubridate)
library(gridExtra)

# 1. Data Loading and Cleaning
# Reading the CSV, handling headers and missing values
raw_data <- read.csv("images_acquisition/BACs_sensors.dat", sep = ",", skip = 1, na.strings = c("NAN", "NA"))
# Removing the first 9 rows (likely sensor calibration/header noise)
clean_data <- raw_data[-c(1:9), ]

# Define custom color palette for publication
sensor_colors <- c("BAC 1" = "#c09048",  
                   "BAC 2" = "#a85400",  
                   "BAC 3" = "#8c3800",  
                   "Air"    = "#385438")

# 2. Data Formatting
# Convert TIMESTAMP to POSIXct date-time format
clean_data$TIMESTAMP <- as.POSIXct(clean_data$TIMESTAMP, format="%Y-%m-%d %H:%M")

# Ensure all character columns are converted to numeric for plotting
clean_data <- clean_data %>%
  mutate_if(is.character, as.numeric)

# Define the observation window
start_date <- as.POSIXct("2025-02-28 16:00:00")
end_date   <- as.POSIXct("2025-03-25 00:00:00")

# 3. Temperature (TC) Data Transformation
temp_long <- clean_data %>%
  select(TIMESTAMP, TC.1., TC.2., TC.3., TC.4.) %>%
  pivot_longer(cols = starts_with("TC"), names_to = "Source", values_to = "Temperature") %>%
  mutate(Source = recode(Source, 
                         "TC.1." = "BAC 1", 
                         "TC.2." = "BAC 2", 
                         "TC.3." = "BAC 3", 
                         "TC.4." = "Air")) %>%
  filter(TIMESTAMP >= start_date & TIMESTAMP <= end_date)

# 4. Volumetric Water Content (VWC) Data Transformation
hum_long <- clean_data %>%
  select(TIMESTAMP, VWC.1., VWC.2., VWC.3.) %>%
  pivot_longer(cols = starts_with("VWC"), names_to = "Source", values_to = "Moisture") %>%
  mutate(Source = recode(Source, 
                         "VWC.1." = "BAC 1", 
                         "VWC.2." = "BAC 2", 
                         "VWC.3." = "BAC 3")) %>%
  filter(TIMESTAMP >= start_date & TIMESTAMP <= end_date)

# 5. Visualization
# Temperature Plot
p_temp <- ggplot(temp_long, aes(x = TIMESTAMP, y = Temperature, color = Source)) +
  geom_line(linewidth = 1) + # Note: 'size' is deprecated in newer ggplot2, use 'linewidth'
  scale_color_manual(values = sensor_colors) +
  labs(title = "Temperature Evolution", x = "Date", y = "Temperature (°C)") +
  theme_bw()

# Moisture Plot
p_hum <- ggplot(hum_long, aes(x = TIMESTAMP, y = Moisture, color = Source)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = sensor_colors) +
  labs(title = "Volumetric Water Content (VWC)", x = "Date", y = "Humidity (%)") +
  theme_bw()

# Display both figures in a single column
grid.arrange(p_temp, p_hum, ncol = 1)