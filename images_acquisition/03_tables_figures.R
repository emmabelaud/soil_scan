
# Ensure output directories exist
table_dir <- "images_acquisition/Tables"
fig_dir   <- "images_acquisition/Figures"
if (!dir.exists(table_dir)) dir.create(table_dir, recursive = TRUE)
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

# -----------------------
# 1) Extract model outputs
# -----------------------

Table_terms <- bind_rows(
  extract_glmm_terms(res_ench$model,  "Enchytraeidae"),
  extract_glmm_terms(res_acar$model,  "Acari"),
  extract_glmm_terms(res_coll$model,  "Collembola")
) %>%
  tidyr::pivot_longer(
    cols = c(p_resolution, p_frequence, p_interaction),
    names_to = "term",
    values_to = "p_value"
  ) %>%
  mutate(
    p_value = round(p_value, 3)
  ) %>%
  arrange(label, term)

Table_quality <- dplyr::bind_rows(
  extract_model_quality(res_ench$model,  "Enchytraeidae"),
  extract_model_quality(res_acar$model,  "Acari"),
  extract_model_quality(res_coll$model,  "Collembola")
)

Table_main <- Table_terms %>%
  dplyr::mutate(
    term = dplyr::recode(term,
                         p_resolution  = "Resolution",
                         p_frequence   = "Frequency",
                         p_interaction = "Resolution × Frequency"),
    p = fmt_p(p_value),
    sig = p_stars(p_value),
    p_cell = paste0(p, " ", sig)
  ) %>%
  dplyr::select(label, term, p_cell) %>%
  tidyr::pivot_wider(names_from = term, values_from = p_cell) %>%
  dplyr::left_join(
    Table_quality %>% dplyr::select(label, R2_marginal, AIC, BIC),
    by = "label"
  ) %>%
  dplyr::arrange(factor(label, levels=c("Enchytraeidae","Acari","Collembola")))%>%
  dplyr::rename(
    `p (Resolution)` = Resolution,
    `p (Frequency)`  = Frequency,
    `p (Res × Freq)` = `Resolution × Frequency`,
    `R²m` = R2_marginal
  )

readr::write_csv(Table_main, "images_acquisition/Tables/Table_model_summary.csv")


# -----------------------
# 2) Ratio table (resolution gains)
# -----------------------

Table_ratios <- bind_rows(
  get_ratio_table(res_ench$model, "Enchytraeidae"),
  get_ratio_table(res_acar$model, "Acari"),
  get_ratio_table(res_coll$model, "Collembola")
) %>%
  arrange(taxon, frequence, contrast)

readr::write_csv(Table_ratios, "images_acquisition/Tables/Table_resolution_ratios_by_taxon.csv")

# -----------------------
# 3) Figures — GLMM effects
# -----------------------

# ---- Taxa: EMM plot
emm_taxa_tbl <- bind_rows(
  as.data.frame(emmeans::emmeans(res_acar$model, ~ resolution * frequence, type = "response")) %>% mutate(taxon = "Acari"),
  as.data.frame(emmeans::emmeans(res_coll$model, ~ resolution * frequence, type = "response")) %>% mutate(taxon = "Collembola"),
  as.data.frame(emmeans::emmeans(res_ench$model, ~ resolution * frequence, type = "response")) %>% mutate(taxon = "Enchytraeidae")
) %>%
  dplyr::rename(lower.CL = asymp.LCL, upper.CL = asymp.UCL) %>%
  dplyr::filter(!is.na(lower.CL), !is.na(upper.CL)) %>%
  dplyr::mutate(
    resolution = as.factor(resolution),
    frequence  = as.factor(frequence)
  )

p_taxa <- emm_taxa_tbl %>%
  ggplot(aes(x = resolution, y = response)) +
  geom_point(position = position_dodge(width = 0.2), size = 1.5) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                width = 0.1, position = position_dodge(width = 0.2)) +
  facet_grid(taxon ~ frequence) +
  scale_y_log10() +
  labs(
    x = "Acquisition interval (min)",
    y = "Detected individuals per image (log scale)"
  ) +
  theme_bw(base_size = 9)+
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "#F3F4F4", colour = NA),
        strip.text = element_text( face = "bold"),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))

print(p_taxa)

ggplot2::ggsave(
  filename = "images_acquisition/Figures/Figure2.png",
  plot     = p_taxa, width = 4, height = 4, dpi = 300,
  bg       = "white"
)

# -----------------------
# 4) shifts in apparent size distributions (quantiles)
# -----------------------
parse_dpi <- function(x) {
  # accepte "2400", 2400, "2400dpi", etc.
  as.numeric(stringr::str_extract(tolower(as.character(x)), "\\d+"))
}

bbox_cols <- c("xmin","xmax","ymin","ymax")
stopifnot(all(bbox_cols %in% names(dat1)))  # on force : tu dis qu'il y a toujours une bbox

det_sizes <- dat1 %>%
  dplyr::filter(taxon %in% c("Enchytraeidae","Acari","Collembola")) %>%  # ou taxon si tu l'as déjà
  dplyr::mutate(
    dpi = parse_dpi(resolution),
    
    xmin = as.numeric(xmin),
    xmax = as.numeric(xmax),
    ymin = as.numeric(ymin),
    ymax = as.numeric(ymax),
    
    width_px  = xmax - xmin,
    height_px = ymax - ymin,
    
    area_px2 = width_px * height_px,
    
    # conversion: 1 pixel = 25.4/dpi mm  =>  area_mm2 = area_px2 * (mm/px)^2
    mm_per_px = 25.4 / dpi,
    area_mm2  = area_px2 * (mm_per_px^2)
  ) %>%
  dplyr::filter(
    !is.na(dpi),
    is.finite(area_mm2),
    width_px > 0, height_px > 0,
    area_mm2 > 0
  )

det_q <- det_sizes %>%
  dplyr::group_by(taxon, resolution, frequence) %>%
  dplyr::summarise(
    q25 = stats::quantile(area_mm2, 0.25, na.rm = TRUE),
    q50 = stats::quantile(area_mm2, 0.50, na.rm = TRUE),
    q75 = stats::quantile(area_mm2, 0.75, na.rm = TRUE),
    n   = dplyr::n(),
    .groups = "drop"
  )

readr::write_csv(det_q, "images_acquisition/Tables/Table_Sx_size_quantiles_q25_q50_q75.csv")

det_q_long <- det_q %>%
  tidyr::pivot_longer(
    cols = c(q25, q50, q75),
    names_to = "quantile",
    values_to = "area_mm2"
  ) %>%
  dplyr::mutate(
    quantile = factor(quantile, levels = c("q25","q50","q75"),
                      labels = c("25th percentile","Median","75th percentile")),
    resolution = factor(resolution),
    frequence  = factor(frequence)
  )

p_quant_taxon <- ggplot(
  det_q_long,
  aes(x = resolution, y = area_mm2, colour = quantile, group = quantile)
) +
  geom_point(size = 1) +
  geom_line(linewidth = 0.5, alpha = 0.8) +
  scale_y_log10() +
  facet_grid(taxon ~ frequence) +
  labs(
    x = "Spatial resolution (dpi)",
    y = "Apparent area (mm², log scale)",
    colour = "Quantile"
  ) +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "#F3F4F4", colour = NA),
        strip.text = element_text( face = "bold"),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold", size=6),
        legend.text  = element_text(size = 6),
        legend.margin = margin(t = -5))+
  ggplot2::scale_colour_manual(
    values = c(
      "25th percentile" = "#95966eff",
      "Median"          = "#ccc04fff",
      "75th percentile" = "#DBCEA5"
    )
  )

print(p_quant_taxon)

ggplot2::ggsave(
  "images_acquisition/Figures/Figure3.png",
  p_quant_taxon,
  width = 4, height = 4, dpi = 300, bg = "white"
)


# -----------------------
# 5) Cumsum curves for each resolution × frequency
#    + summaries (mean/sd across runs)
#    + curve fitting per run and summary across runs
# -----------------------

# ---- 5.1 Build raw cumulative curves (per run_id) ----
cum_all <- make_cumsum_df(
  img_taxa,
  group_cols = c("taxon","resolution","frequence","run_id"),
  value_col  = "n_img"
) %>%
  transmute(
    taxon      = as.character(taxon),
    resolution = factor(resolution),
    frequence  = factor(frequence),
    datetime,
    run_id,
    x = effort_img,
    y = cum_obs
  ) %>%
  mutate(
    taxon = factor(taxon, levels = c("Acari","Collembola","Enchytraeidae"))
  )

readr::write_csv(cum_all, "images_acquisition/Tables/Table_cumsum_curves_raw.csv")

# ---- 5.2 Summary mean/sd across runs (for plotting mean curves) ----
cum_all_sum <- cum_all %>%
  group_by(taxon, resolution, frequence, x) %>%
  summarise(
    y_mean = mean(y, na.rm = TRUE),
    y_sd   = sd(y, na.rm = TRUE),
    n_runs = n_distinct(run_id),
    y_se   = y_sd / sqrt(n_runs),
    .groups = "drop"
  )

readr::write_csv(cum_all_sum, "images_acquisition/Tables/Table_cumsum_curves_mean_sd.csv")

# -----------------------
# 6) Fit curve types to each cumsum curve (per run), then summarise across runs
# -----------------------

# ---- 6.1 Fit candidate curves per run ----
Table_curves_all <- cum_all %>%
  group_by(taxon, resolution, frequence, run_id) %>%
  group_modify(~{
    fits <- fit_curve_candidates(.x %>% dplyr::select(x, y))
    if (nrow(fits) == 0) return(tibble())
    fits
  }) %>%
  ungroup() %>%
  arrange(taxon, resolution, frequence, run_id, AIC)

readr::write_csv(Table_curves_all, "images_acquisition/Tables/Table_cumsum_curve_fits_ALL.csv")

# ---- 6.2 Pick best curve per run (AIC-min) ----
Table_curves_best <- Table_curves_all %>%
  group_by(taxon, resolution, frequence, run_id) %>%
  slice_min(order_by = AIC, n = 1, with_ties = FALSE) %>%
  ungroup()

readr::write_csv(Table_curves_best, "images_acquisition/Tables/Table_cumsum_curve_fits_BEST_by_run.csv")

# ---- 6.3 Summarise best-fit parameters across runs (mean/sd) ----
# If curve_type varies across runs within a combination, this keeps curve_type explicit.
Table_curves_best_sum <- Table_curves_best %>%
  group_by(taxon, resolution, frequence, curve_type) %>%
  summarise(
    param_a    = mean(param_a, na.rm = TRUE),
    param_a_sd = sd(param_a,   na.rm = TRUE),
    param_b    = mean(param_b, na.rm = TRUE),
    param_b_sd = sd(param_b,   na.rm = TRUE),
    param_c    = mean(param_c, na.rm = TRUE),
    param_c_sd = sd(param_c,   na.rm = TRUE),
    n_runs     = n(),
    .groups = "drop"
  ) %>%
  arrange(taxon, resolution, frequence, desc(n_runs))

readr::write_csv(Table_curves_best_sum, "images_acquisition/Tables/Table_cumsum_curve_params_mean_sd.csv")

# Optional: if you want ONE curve_type per (taxon,resolution,frequence),
# keep the most frequent curve_type (mode) and its parameter summaries.
Table_curves_best_mode <- Table_curves_best_sum %>%
  group_by(taxon, resolution, frequence) %>%
  slice_max(order_by = n_runs, n = 1, with_ties = FALSE) %>%
  ungroup()

readr::write_csv(Table_curves_best_mode, "images_acquisition/Tables/Table_cumsum_curve_params_mode.csv")

# -----------------------
# 7) Graphical outputs (GLMM diagnostics + cumsum curves)
# -----------------------

# --- 7.1 GLMM diagnostics (DHARMa) ---
save_dharma <- function(model, name, n = 2000) {
  sim <- DHARMa::simulateResiduals(model, n = n, plot = FALSE)
  out_pdf <- file.path(fig_dir, paste0("DHARMa_", name, ".pdf"))
  grDevices::pdf(out_pdf, width = 8, height = 7)
  plot(sim)
  grDevices::dev.off()
  invisible(out_pdf)
}

invisible(save_dharma(res_ench$model,  "Enchytraeidae"))
invisible(save_dharma(res_acar$model,  "Acari"))
invisible(save_dharma(res_coll$model,  "Collembola"))


# --- 7.2 Best-fit curves overlay ---
predict_curve <- function(curve_type, a, b, c, x) {
  if (is.na(curve_type) || curve_type == "") return(rep(NA_real_, length(x)))
  if (curve_type == "linear") {
    a + b*x
  } else if (curve_type == "neg_exp") {
    a * (1 - exp(-b*x))
  } else if (curve_type == "michaelis_menten") {
    (a*x) / (b + x)
  } else if (curve_type == "logarithmic") {
    a + b*log(pmax(x, 1e-9))
  } else if (curve_type == "power") {
    a * (pmax(x, 1e-9)^b)
  } else if (curve_type == "logistic") {
    a / (1 + exp(-b*(x - c)))
  } else {
    rep(NA_real_, length(x))
  }
}

pred_best <- cum_all %>%
  filter(taxon != "Total") %>%
  group_by(taxon, resolution, frequence) %>%
  summarise(x_min = min(x), x_max = max(x), .groups = "drop") %>%
  left_join(Table_curves_best_mode, by = c("taxon","resolution","frequence")) %>%
  filter(!is.na(curve_type)) %>%
  mutate(x_grid = purrr::map2(x_min, x_max, ~seq(.x, .y, length.out = 200))) %>%
  mutate(
    y_hat = purrr::pmap(
      list(curve_type, param_a, param_b, param_c, x_grid),
      ~predict_curve(..1, ..2, ..3, ..4, ..5)
    )
  ) %>%
  tidyr::unnest(c(x_grid, y_hat)) %>%
  rename(x = x_grid, y = y_hat)

p_fit <- ggplot() +
  geom_line(
    data = subset(cum_all, taxon != "Total"),
    aes(x = x, y = log(y), group = interaction(run_id, resolution, frequence), colour = resolution),
    linewidth = 0.3, alpha = 0.25
  ) +
  geom_line(
    data = pred_best,
    aes(x = x, y = log(pmax(y, 1e-9)), group = interaction(resolution, frequence), colour = resolution),
    linewidth = 0.8, alpha = 0.8
  ) +
  facet_grid(taxon ~ frequence) +
  labs(x = "Number of events (images)", y = "Cumulative detections (log scale)", colour = "Resolution") +
  scale_x_continuous(limits = c(1, 11), breaks = c(3,6,9))+
  scale_color_manual(values = c('#FF4400', "#e89502ff", "#f0cb37ff"))+
  theme_bw(base_size = 9)+
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "#F3F4F4", colour = NA),
        strip.text = element_text( face = "bold"),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold", size=6),
        legend.text  = element_text(size = 6),
        legend.margin = margin(t = -5),
        panel.grid.minor = element_blank())

ggsave(
  filename = file.path(fig_dir, "Figure4.png"),
  plot = p_fit, width = 4, height = 4, dpi = 300
)
