# Load required libraries
library(ggplot2)
library(tidyr)
library(dplyr)
library(lubridate)
library(gridExtra)

# 1. Data Loading and Cleaning
# Reading the CSV, handling headers and missing values
raw_data <- read.csv("images_acquisition/data/BACs_sensors.dat", sep = ",", skip = 1, na.strings = c("NAN", "NA"))
# Removing the first 9 rows (likely sensor calibration/header noise)
clean_data <- raw_data[-c(1:9), ]

# Define custom color palette for publication
sensor_colors <- c("Container 1" = "#c09048",  
                   "Container 2" = "#a85400",  
                   "Container 3" = "#8c3800",  
                   "Air"         = "#385438")

# Define the experimental cycles
cycles <- data.frame(
  name = c("Cycle 1", "Cycle 2", "Cycle 3", "Cycle 4"),
  start = as.POSIXct(c("2025-03-03", "2025-03-14", "2025-03-25", "2025-04-06")),
  end   = as.POSIXct(c("2025-03-13", "2025-03-24", "2025-04-05", "2025-04-16"))
)

# 2. Data Formatting
# Convert TIMESTAMP to POSIXct date-time format
clean_data$TIMESTAMP <- as.POSIXct(clean_data$TIMESTAMP, format="%Y-%m-%d %H:%M")

# Ensure all character columns are converted to numeric for plotting
clean_data <- clean_data %>%
  mutate_if(is.character, as.numeric)

# Define the observation window
start_date <- as.POSIXct("2025-02-28 16:00:00")
end_date   <- as.POSIXct("2025-04-16 00:00:00")

# 3. Temperature (TC) Data Transformation
temp_long <- clean_data %>%
  select(TIMESTAMP, TC.1., TC.2., TC.3., TC.4.) %>%
  pivot_longer(cols = starts_with("TC"), names_to = "Source", values_to = "Temperature") %>%
  mutate(Source = recode(Source, 
                         "TC.1." = "Container 1", 
                         "TC.2." = "Container 2", 
                         "TC.3." = "Container 3", 
                         "TC.4." = "Air")) %>%
  filter(TIMESTAMP >= start_date & TIMESTAMP <= end_date)

# 4. Volumetric Water Content (VWC) Data Transformation
hum_long <- clean_data %>%
  select(TIMESTAMP, VWC.1., VWC.2., VWC.3.) %>%
  pivot_longer(cols = starts_with("VWC"), names_to = "Source", values_to = "Humidity") %>%
  mutate(Source = recode(Source, 
                         "VWC.1." = "Container 1", 
                         "VWC.2." = "Container 2", 
                         "VWC.3." = "Container 3")) %>%
  filter(TIMESTAMP >= start_date & TIMESTAMP <= end_date) %>%
  group_by(Source) %>%
  mutate(Humidity_norm = (Humidity - min(Humidity, na.rm = TRUE)) / 
           (max(Humidity, na.rm = TRUE) - min(Humidity, na.rm = TRUE))) %>%
  ungroup()

# 5. Visualization
# Custom function to add cycles to any plot
add_cycle_layers <- function() {
  list(
    geom_rect(data = cycles, 
              aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = name), 
              inherit.aes = FALSE, alpha = 0.2),
    geom_text(data = cycles, 
              aes(x = start + (end - start)/2, y = Inf, label = name), 
              inherit.aes = FALSE, vjust = 2, size = 3, fontface = "bold", color = "black"),
    scale_fill_manual(values = c("Cycle 1" = "grey90", "Cycle 2" = "grey80", 
                                 "Cycle 3" = "grey90", "Cycle 4" = "grey80"), 
                      guide = "none")
  )
}

# Temperature Plot
p_temp <- ggplot(temp_long, aes(x = TIMESTAMP, y = Temperature, color = Source)) +
  add_cycle_layers() +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = sensor_colors) +
  labs(title = "Temperature Evolution across Cycles", x = "Date", y = "Temperature (°C)") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())


# Humidity Plot
p_hum <- ggplot(hum_long, aes(x = TIMESTAMP, y = Humidity_norm, color = Source)) +
  add_cycle_layers() + 
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = sensor_colors) +
  labs(title = "Relative Humidity Variation (Normalized)", x = "Date", y = "Normalized VWC (0-1)") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())

# Combine with shared legend
combined <- gridExtra::grid.arrange(
  p_hum + theme(axis.title.x = element_blank(), legend.position = "none"),
  p_temp + theme(legend.position = "bottom"),
  ncol = 1,
  heights = c(1, 1.2) # Give slightly more room to the bottom plot for the legend
)

# Save with a safe background
ggsave(
  filename = "images_acquisition/Figures/Figure_sensors_data.png", 
  plot = combined,
  bg = "white",
  width = 8, height = 10, dpi = 300
)
