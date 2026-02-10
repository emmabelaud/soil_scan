# ============================================================
# Pipeline for scanner calibration
#  - Build run-level replicates (blocks of consecutive images)
#  - GLMM (negative binomial; optional zero-inflation)
#  - Full diagnostics with DHARMa
# ============================================================

# -----------------------
# 0) Packages
# -----------------------
pkgs <- c(
  "dplyr","tidyr","stringr","lubridate","forcats",
  "glmmTMB","DHARMa","emmeans","ggplot2","performance", "readr"
)
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install) > 0) install.packages(to_install)
invisible(lapply(pkgs, library, character.only = TRUE))

# -----------------------
# 1) Helper functions
# -----------------------

# Parse datetime like "09/03/2025  22:15:00" (note: possible double spaces)
parse_datetime <- function(x, tz = "Europe/Paris") {
  x <- stringr::str_squish(as.character(x))
  
  # 1) try ISO first: "2025-03-09 22:15:00"
  out <- suppressWarnings(lubridate::ymd_hms(x, tz = tz, quiet = TRUE))
  miss <- is.na(out)
  
  # 2) fallback French format: "09/03/2025  22:15:00"
  if (any(miss)) {
    out2 <- as.POSIXct(x[miss], format = "%d/%m/%Y %H:%M:%S", tz = tz)
    out[miss] <- out2
  }
  out
}

# Convert "1200dpi" -> 1200 numeric (keep factor later for modelling)
parse_resolution <- function(x) {
  as.integer(stringr::str_extract(as.character(x), "\\d+"))
}

# Convert "15min","1h","6h" -> minutes numeric
parse_frequency_min <- function(x) {
  x <- tolower(as.character(x))
  dplyr::case_when(
    stringr::str_detect(x, "min") ~ as.numeric(stringr::str_extract(x, "\\d+")),
    stringr::str_detect(x, "h")   ~ 60 * as.numeric(stringr::str_extract(x, "\\d+")),
    TRUE ~ NA_real_
  )
}

# Build run_id: blocks of consecutive images within bac×cycle×res×freq
assign_run_id <- function(df_img) {
  df_img %>%
    arrange(datetime) %>%
    mutate(
      dt_min = as.numeric(difftime(datetime, dplyr::lag(datetime), units = "mins")),
      # A new run starts if: first image OR gap bigger than 1.5x expected interval
      new_run = dplyr::if_else(is.na(dt_min) | dt_min > 1.5 * frequence_min, 1L, 0L),
      run_seq = cumsum(new_run),
      run_id  = paste(bac, cycle, resolution, frequence, sprintf("run%02d", run_seq), sep = "_")
    ) %>%
    select(-dt_min, -new_run)
}

# Quick DHARMa diagnostic wrapper
dharma_report <- function(model, name = deparse(substitute(model)), n = 2000) {
  cat("\n\n========================\nDHARMa diagnostics:", name, "\n========================\n")
  sim <- DHARMa::simulateResiduals(model, n = n, plot = FALSE)
  
  print(plot(sim))                       # general residual checks
  print(DHARMa::testDispersion(sim))     # over/under dispersion
  print(DHARMa::testZeroInflation(sim))  # excess zeros
  print(DHARMa::testOutliers(sim))       # outliers
  # Optional: temporal autocorrelation if you keep image-level models
  # print(DHARMa::testTemporalAutocorrelation(sim, time = df$datetime))
  invisible(sim)
}

# DHARMA figures
export_dharma_png <- function(model, prefix, n = 2000,
                              width = 1800, height = 1400, res = 300) {
  
  # 1) Simulated residuals object
  sim <- DHARMa::simulateResiduals(model, n = n, plot = FALSE)
  
  # 2) Main DHARMa diagnostic plot (4 panels)
  png(paste0(prefix, "_DHARMa_overview.png"), width = width, height = height, res = res)
  plot(sim)
  dev.off()
  
  # 3) Optional: specific tests as separate figures
  png(paste0(prefix, "_DHARMa_dispersion.png"), width = 1400, height = 1200, res = res)
  print(DHARMa::testDispersion(sim, plot = TRUE))
  dev.off()
  
  png(paste0(prefix, "_DHARMa_zeroinflation.png"), width = 1400, height = 1200, res = res)
  print(DHARMa::testZeroInflation(sim, plot = TRUE))
  dev.off()
  
  png(paste0(prefix, "_DHARMa_outliers.png"), width = 1400, height = 1200, res = res)
  print(DHARMa::testOutliers(sim, plot = TRUE))
  dev.off()
  
  # 4) Save the test outputs as a text log (nice for supplementary material)
  sink(paste0(prefix, "_DHARMa_tests.txt"))
  cat("Model:", prefix, "\n\n")
  print(DHARMa::testDispersion(sim))
  cat("\n")
  print(DHARMa::testZeroInflation(sim))
  cat("\n")
  print(DHARMa::testOutliers(sim))
  sink()
  
  return(invisible(sim))
}
# -----------------------
# 2) Read / prepare data
# -----------------------
df <- readr::read_csv("images_acquisition/database_calibration.csv")  # <-- set your real import

dat0 <- df %>%
  mutate(
    datetime = parse_datetime(datetime),
    resolution_num = parse_resolution(resolution),
    frequence_min  = parse_frequency_min(frequence),
    bac  = as.factor(bac),
    cycle = as.factor(cycle),
    resolution = factor(resolution_num, levels = sort(unique(resolution_num))),
    frequence = factor(frequence_min,  levels = sort(unique(frequence_min)))
  ) %>%
  filter(!is.na(datetime), !is.na(resolution_num), !is.na(frequence_min))

# Define "fauna" classes to keep
# Keep a clean fauna indicator + optional taxa mapping.
dat1 <- dat0 %>%
  mutate(
    class = as.character(class),
    class = dplyr::na_if(class, ""),
    class_clean = tolower(class),
    is_fauna = !is.na(class_clean) &
      !class_clean %in% c("na","<na>","root","som","unknown"),
    taxon = dplyr::case_when(
      class_clean %in% c("enchytraeidae","enchytraeid","enchytraeids") ~ "Enchytraeidae",
      class_clean %in% c("collembola","collembolan","collembolans")    ~ "Collembola",
      class_clean %in% c("acari","mesostigmata","oribatida","prostigmata","astigmata") ~ "Acari",
      TRUE ~ "Other"
    )
  )

# -----------------------
# 3) Build IMAGE-level table (unit = initial_image_name)
# -----------------------
# n_total: total fauna detections per image (all fauna classes)
img_total <- dat1 %>%
  filter(is_fauna) %>%
  group_by(bac, cycle, resolution, frequence, frequence_min, datetime, initial_image_name) %>%
  summarise(n_img = dplyr::n(), .groups = "drop") %>%
  group_by(bac, cycle, resolution, frequence) %>%
  assign_run_id() %>%
  ungroup()

# Taxon-specific image-level tables
img_taxa <- dat1 %>%
  filter(is_fauna, taxon %in% c("Enchytraeidae","Acari","Collembola")) %>%
  group_by(bac, cycle, resolution, frequence, frequence_min, taxon, datetime, initial_image_name) %>%
  summarise(n_img = dplyr::n(), .groups = "drop") %>%
  group_by(bac, cycle, resolution, frequence, taxon) %>%
  assign_run_id() %>%
  ungroup()

# -----------------------
# 4) Build RUN-level table
# -----------------------
# This is the key step to avoid pseudo-replication: aggregate the 12 images per run.
run_total <- img_total %>%
  group_by(bac, cycle, resolution, frequence, run_id) %>%
  summarise(
    sum_n   = sum(n_img),
    mean_n  = mean(n_img),
    n_images = n(),
    start_time = min(datetime),
    end_time   = max(datetime),
    .groups = "drop"
  ) %>%
  mutate(
    log_effort = log(n_images),
    # optional: coarse day/night covariate (if useful later)
    hour = lubridate::hour(start_time),
    period = dplyr::case_when(
      hour >= 6  & hour < 18 ~ "day",
      TRUE ~ "night"
    ),
    period = factor(period)
  )

run_taxa <- img_taxa %>%
  group_by(bac, cycle, resolution, frequence, taxon, run_id) %>%
  summarise(
    sum_n   = sum(n_img),
    mean_n  = mean(n_img),
    n_images = n(),
    start_time = min(datetime),
    end_time   = max(datetime),
    .groups = "drop"
  ) %>%
  mutate(
    log_effort = log(n_images),
    hour = lubridate::hour(start_time),
    period = dplyr::case_when(
      hour >= 6  & hour < 18 ~ "day",
      TRUE ~ "night"
    ),
    period = factor(period),
    taxon = factor(taxon)
  )

# Sanity checks
stopifnot(all(run_total$n_images > 0))
cat("\nRun-level sample sizes (total):\n")
print(run_total %>% count(resolution, frequence) %>% arrange(resolution, frequence))

# -----------------------
# 5) Fit GLMMs (Negative Binomial) + diagnostics
# -----------------------
# Recommended response: sum_n with offset(log(n_images)) OR mean_n without offset.
# Using sum_n + offset is robust if n_images varies.
# Random effects: bac + cycle (minimum).

m_total_nb <- glmmTMB(
  mean_n ~ resolution * frequence + (1|bac) + (1|cycle),
  family = nbinom2(),
  data = run_total
)

summary(m_total_nb)
performance::check_model(m_total_nb)  # quick visual check (non-DHARMa)
#sim_total <- dharma_report(m_total_nb, name = "m_total_nb")
sim_total <- export_dharma_png(m_total_nb, "images_acquisition/DHARMA/SuppFig_S1_total_fauna")

# Sensitivity to error distribution choice
m_total_pois <- update(m_total_nb, family = poisson)
AIC(m_total_nb, m_total_pois)

# If DHARMa indicates strong zero-inflation, refit as ZINB:
zi_test <- DHARMa::testZeroInflation(sim_total)
if (zi_test$p.value < 0.05) {
  message("\nDetected excess zeros -> refitting with zero-inflated NB...")
  m_total_zinb <- glmmTMB(
    sum_n ~ resolution * frequence + offset(log_effort) + (1|bac) + (1|cycle),
    ziformula = ~1,
    family = nbinom2(),
    data = run_total
  )
  summary(m_total_zinb)
  sim_total <- dharma_report(m_total_zinb, name = "m_total_zinb")
  m_total_final <- m_total_zinb
} else {
  m_total_final <- m_total_nb
}

# -----------------------
# 6) Post-hoc: estimated marginal means (ASE-friendly tables)
# -----------------------
emm_total <- emmeans::emmeans(m_total_final, ~ resolution * frequence, type = "response")
emm_total_tbl <- as.data.frame(emm_total)%>%
  dplyr::rename(lower.CL = asymp.LCL, upper.CL = asymp.UCL) %>%
  dplyr::filter(!is.na(lower.CL), !is.na(upper.CL)) %>%
  dplyr::mutate(
    resolution = as.factor(resolution),
    frequence  = as.factor(frequence)
  )
print(emm_total_tbl)

# Pairwise contrasts if needed:
pairs(emm_total)

# -----------------------
# 7) Plot (publication-oriented, simple)
# -----------------------
p_total <- emm_total_tbl %>%
  ggplot(aes(x = frequence, y = response, colour = resolution)) +
  geom_point(position = position_dodge(width = 0.2)) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                width = 0.1, position = position_dodge(width = 0.2)) +
  #facet_wrap(~ resolution, scales = "free_y") +
  labs(x = "Acquisition interval (min)",
       y = "Detected activity (expected count per image)",
       title = "Effects of resolution and acquisition frequency (total fauna)") +
  theme_bw()
print(p_total)

ggplot2::ggsave(
  filename = "images_acquisition/Figures/Fig_total_activity.png",
  plot     = p_total,
  width    = 180,      # largeur en mm (≈ une colonne ASE)
  height   = 120,      # hauteur en mm
  units    = "mm",
  dpi      = 300,
  bg       = "white"
)

# -----------------------
# 8) Taxon-specific models
# -----------------------
fit_taxon_model <- function(taxon_name, dat) {
  d <- dat %>% filter(taxon == taxon_name)
  
  m_nb <- glmmTMB(
    sum_n ~ resolution * frequence + offset(log_effort) + (1|bac) + (1|cycle),
    family = nbinom2(),
    data = d
  )
  
  sim <- dharma_report(m_nb, name = paste0("m_", taxon_name, "_nb"))
  zi <- DHARMa::testZeroInflation(sim)
  
  if (zi$p.value < 0.05) {
    message("\n[", taxon_name, "] excess zeros -> ZINB")
    m_zi <- glmmTMB(
      sum_n ~ resolution * frequence + offset(log_effort) + (1|bac) + (1|cycle),
      ziformula = ~1,
      family = nbinom2(),
      data = d
    )
    dharma_report(m_zi, name = paste0("m_", taxon_name, "_zinb"))
    m_final <- m_zi
  } else {
    m_final <- m_nb
  }
  
  emm <- emmeans::emmeans(m_final, ~ resolution * frequence, type = "response")
  list(model = m_final, emmeans = as.data.frame(emm))
}

res_ench <- fit_taxon_model("Enchytraeidae", run_taxa)
res_acar <- fit_taxon_model("Acari", run_taxa)
res_coll <- fit_taxon_model("Collembola", run_taxa)

# Combine emmeans into one table
emm_taxa_tbl <- bind_rows(
  res_ench$emmeans %>% mutate(taxon = "Enchytraeidae"),
  res_acar$emmeans %>% mutate(taxon = "Acari"),
  res_coll$emmeans %>% mutate(taxon = "Collembola")
)%>%
  dplyr::rename(lower.CL = asymp.LCL, upper.CL = asymp.UCL) %>%
  dplyr::filter(!is.na(lower.CL), !is.na(upper.CL)) %>%
  dplyr::mutate(
    resolution = as.factor(resolution),
    frequence  = as.factor(frequence))
print(emm_taxa_tbl)

# Taxon plots
p_taxa <- emm_taxa_tbl %>%
  ggplot(aes(x = frequence, y = response, colour = resolution)) +
  geom_point(position = position_dodge(width = 0.2)) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                width = 0.1, position = position_dodge(width = 0.2)) +
  facet_grid(taxon ~ ., scales = "free_y") +
  labs(x = "Acquisition interval (min)",
       y = "Detected activity (expected count per image)") +
  theme_bw()
print(p_taxa)

ggplot2::ggsave(
  filename = "images_acquisition/Figures/Fig_taxa_activity.png",
  plot     = p_taxa,
  width    = 180,      # largeur en mm (≈ une colonne ASE)
  height   = 120,      # hauteur en mm
  units    = "mm",
  dpi      = 300,
  bg       = "white"
)


sim_ench  <- export_dharma_png(res_ench$model, "images_acquisition/DHARMA/SuppFig_S2_Enchytraeidae")
sim_acar  <- export_dharma_png(res_acar$model, "images_acquisition/DHARMA/SuppFig_S3_Acari")
sim_coll  <- export_dharma_png(res_coll$model, "images_acquisition/DHARMA/SuppFig_S4_Collembola")

# -----------------------
# 9) Optional: sensitivity analysis at IMAGE-level
# -----------------------
# Sometimes reviewers ask "why aggregate?" You can show that inferences are consistent.
# This model is more fragile (autocorrelation), but defensible as sensitivity.
# Use only if needed.

m_img_nb <- glmmTMB(
  n_img ~ resolution * frequence + (1|bac) + (1|cycle) + (1|run_id),
  family = nbinom2(),
  data = img_total
)
summary(m_img_nb)
sim_img <- dharma_report(m_img_nb, name = "m_img_nb")
# NOTE: test temporal autocorrelation in DHARMa at image-level:
DHARMa::testTemporalAutocorrelation(sim_img, time = img_total$datetime)


# Resolution- and frequency-dependent shifts in apparent size distributions
# --- build "bbox size" table at detection level
parse_dpi <- function(x) {
  as.numeric(stringr::str_extract(tolower(as.character(x)), "\\d+"))
}

det_sizes <- dat1 %>%
  filter(taxon %in% c("Enchytraeidae","Acari","Collembola")) %>%
  mutate(
    dpi = parse_dpi(resolution),
    width_px  = as.numeric(xmax) - as.numeric(xmin),
    height_px = as.numeric(ymax) - as.numeric(ymin),
    area_px   = width_px * height_px,
    size_px   = sqrt(area_px),                      # proxy linéaire
    px_to_mm  = 2.54 / dpi,
    size_mm   = area_px * px_to_mm,
    width_mm  = width_px * px_to_mm,
    height_mm = height_px * px_to_mm) %>%
  filter(
    !is.na(dpi),
    is.finite(size_mm),
    width_px > 0, height_px > 0 ) 

det_sizes %>% count(resolution, dpi)
summary(det_sizes$size_mm)


det_q <- det_sizes %>%
  filter(
    taxon %in% c("Enchytraeidae","Acari","Collembola"),
    is.finite(size_mm),
    size_mm > 0
  ) %>%
  group_by(taxon, resolution, frequence) %>%
  summarise(
    q25 = quantile(size_mm, 0.25, na.rm = TRUE),
    q50 = quantile(size_mm, 0.50, na.rm = TRUE),
    q75 = quantile(size_mm, 0.75, na.rm = TRUE),
    n   = n(),
    .groups = "drop"
  )

# Sauvegarde table annexe
write_csv(det_q, "images_acquisition/Tables/Table_Sx_size_quantiles_q25_q50_q75.csv")

det_q_long <- det_q %>%
  pivot_longer(
    cols = c(q25, q50, q75),
    names_to = "quantile",
    values_to = "size_mm"
  ) %>%
  mutate(
    quantile = factor(quantile, levels = c("q25","q50","q75")),
    resolution = factor(resolution),
    frequence  = factor(frequence)
  )

p_quant_taxon <- det_q_long %>%
  ggplot(aes(x = resolution, y = size_mm, colour = quantile, group = quantile)) +
  geom_point(size = 2) +
  geom_line(linewidth = 0.9) +
  scale_y_log10() +
  facet_grid(taxon ~ frequence) +
  labs(
    x = "Spatial resolution (dpi)",
    y = "Apparent area (mm2, log scale)",
    colour = "Quantile")+ 
  theme_bw() +
  theme(legend.position = "bottom")+
  scale_colour_manual(
    values = c(
      q25 = "#deebf7",  # bleu très clair
      q50 = "#9ecae1",  # bleu moyen
      q75 = "#3182bd"   # bleu foncé
    ),
    labels = c(
      q25 = "25th percentile",
      q50 = "Median",
      q75 = "75th percentile"
    )
  )

ggsave(
  "images_acquisition/Figures/SuppFig_Sx_size_quantiles_by_taxon.png",
  p_quant_taxon,
  width = 220, height = 180, units = "mm", dpi = 300, bg = "white"
)

# ============================================================
# End   :-)
# ============================================================
