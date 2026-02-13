# Chargement des bibliothèques nécessaires
library(ggplot2)
library(tidyr)
library(dplyr)
library(lubridate)

# Lecture des données CSV
data <- read.csv("images_acquisition/BACs_sensors.dat", sep = ",", skip = 1,  na.strings = c("NAN", "NA"))
data <- data[-c(1:9), ]

# Liste des couleurs personnalisées
bac_colors <- c("BAC 1" = "#c09048",   # Couleur pour BAC 1
                "BAC 2" = "#a85400",   # Couleur pour BAC 2
                "BAC 3" = "#8c3800",   # Couleur pour BAC 3
                "Air" = "#385438")     # Couleur pour Air

# Conversion de la colonne TIMESTAMP en format date-temps
data$TIMESTAMP <- as.POSIXct(data$TIMESTAMP, format="%Y-%m-%d %H:%M")

# Remplacement des valeurs "NAN" par NA
data <- data %>%
  mutate_if(is.character, as.numeric)   # Conversion des colonnes en numérique

# Définition de la période à visualiser
date_debut <- as.POSIXct("2025-02-28 16:00:00")
date_fin <- as.POSIXct("2025-03-25 00:00:00")

# Transformation des données pour le graphique des températures (TC)
temp_data <- data %>%
  select(TIMESTAMP, TC.1., TC.2., TC.3., TC.4.) %>%
  pivot_longer(cols = starts_with("TC"), names_to = "Bac", values_to = "Température") %>%
  mutate(Bac = recode(Bac, 
                      "TC.1." = "BAC 1", 
                      "TC.2." = "BAC 2", 
                      "TC.3." = "BAC 3", 
                      "TC.4." = "Air")) %>%
  filter(TIMESTAMP >= date_debut & TIMESTAMP <= date_fin)

# Transformation des données pour le graphique de l'humidité volumique (VWC)
hum_data <- data %>%
  select(TIMESTAMP, VWC.1., VWC.2., VWC.3.) %>%
  pivot_longer(cols = starts_with("VWC"), names_to = "Bac", values_to = "Humidité") %>%
  mutate(Bac = recode(Bac, 
                      "VWC.1." = "BAC 1", 
                      "VWC.2." = "BAC 2", 
                      "VWC.3." = "BAC 3")) %>%
  filter(TIMESTAMP >= date_debut & TIMESTAMP <= date_fin)

# Graphique des températures
temp <- ggplot(temp_data, aes(x = TIMESTAMP, y = Température, color = Bac)) +
  geom_line(size = 2) +
  scale_color_manual(values = bac_colors) +  # Applique la palette personnalisée
  labs(title = "Température", x = "Date", y = "Température (°C)") +
  theme_bw()

# Graphique de l'humidité
hum <- ggplot(hum_data, aes(x = TIMESTAMP, y = Humidité, color = Bac)) +
  geom_line(size = 2) +
  scale_color_manual(values = bac_colors) +  # Applique la palette personnalisée
  labs(title = "Humidité", x = "Date", y = "Humidité") +
  theme_bw()

# Affichage des graphiques
gridExtra::grid.arrange(temp, hum, ncol = 1)
