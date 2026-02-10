# ============================================================
# Table Sx — Ratios of detectability gains (taxon-specific models)
# Output: CSV ready for Supplementary Material
# ============================================================

library(dplyr)
library(emmeans)
library(readr)

# ---- helper: compute ratios from one model
# Ratios computed within each frequency level:
# - 1200/600
# - 2400/600
# (and optionally 2400/1200 if you want)
get_ratio_table <- function(model, taxon_label,
                            contrasts = list(
                              "1200/600" = c(-1, 1, 0),
                              "2400/600" = c(-1, 0, 1),
                              "2400/1200" = c(0, -1, 1) # optional
                            )) {
  
  # EMMs on response scale, by frequency
  emm <- emmeans(model, ~ resolution | frequence, type = "response")
  
  # Ratio contrasts (log scale internally; returned as ratios)
  rat <- contrast(
    emm,
    method = contrasts,
    by = "frequence",
    ratio = TRUE
  )
  
  # Convert to a clean table
  as.data.frame(summary(rat)) %>%
    transmute(
      taxon = taxon_label,
      frequence = as.character(frequence),
      contrast = as.character(contrast),
      ratio = ratio,
      SE = SE,
      df = df
    )
}

# ---- plug your models here
m_ench <- res_ench$model
m_acar <- res_acar$model
m_coll <- res_coll$model

tab_ench <- get_ratio_table(m_ench, "Enchytraeidae")
tab_acar <- get_ratio_table(m_acar, "Acari")
tab_coll <- get_ratio_table(m_coll, "Collembola")

Table_Sx <- bind_rows(tab_ench, tab_acar, tab_coll) %>%
  mutate(
    ratio = round(ratio, 1), 
    SE = round(SE, 1),
    frequence = factor(frequence, levels = c("15","60","360")),
    contrast  = factor(contrast, levels = c("1200/600","2400/600", "2400/1200"))
  ) %>%
  arrange(taxon, frequence, contrast) %>%
  select(-df) %>%
  filter(!is.na(ratio))

# ---- export
readr::write_csv(Table_Sx, "Table_Sx_resolution_ratios_by_taxon.csv")

