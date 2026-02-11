# ============================================================
# ASE scanner calibration — merged pipeline
#
# What it does
#   1) Reads detection database (one row = one detection)
#   2) Builds IMAGE-level counts and RUN-level replicates
#   3) Fits GLMMs (NB or ZINB) for:
#        - Total fauna
#        - Acari, Collembola, Enchytraeidae
#   4) Extracts (i) marginal/conditional R2, and (ii) LRT p-values for:
#        - frequency (frequence)
#        - resolution
#        - interaction
#   5) Computes cumulative detections curves (cumsum) for each
#      resolution × frequency combination, for total + 3 taxa
#   6) Fits candidate curve types to each cumsum curve and exports
#      a table with the best curve type + parameter estimates
#
# Notes
#   - "pouvoir explicatif" is reported here as marginal/conditional R2
#     (Nakagawa) + likelihood-ratio chi-square for each term.
#   - Main-effect tests in presence of interaction are handled as:
#       * Interaction tested alone (drop interaction term)
#       * Resolution importance: drop resolution AND the interaction
#       * Frequency importance:  drop frequency AND the interaction
#     This yields interpretable "global" tests.
# ============================================================

# -----------------------
# 0) Packages
# -----------------------
pkgs <- c(
  "dplyr","tidyr","stringr","lubridate","forcats","purrr","tibble",
  "glmmTMB","DHARMa","emmeans","ggplot2","performance","readr"
)

to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install) > 0) install.packages(to_install)
invisible(lapply(pkgs, library, character.only = TRUE))

# -----------------------
# 1) Helper functions (parsing + run_id)
# -----------------------

parse_datetime <- function(x, tz = "Europe/Paris") {
  x <- stringr::str_squish(as.character(x))
  out <- suppressWarnings(lubridate::ymd_hms(x, tz = tz, quiet = TRUE))
  miss <- is.na(out)
  if (any(miss)) {
    out2 <- as.POSIXct(x[miss], format = "%d/%m/%Y %H:%M:%S", tz = tz)
    out[miss] <- out2
  }
  out
}

parse_resolution <- function(x) {
  as.integer(stringr::str_extract(as.character(x), "\\d+"))
}

parse_frequency_min <- function(x) {
  x <- tolower(as.character(x))
  dplyr::case_when(
    stringr::str_detect(x, "min") ~ as.numeric(stringr::str_extract(x, "\\d+")),
    stringr::str_detect(x, "h")   ~ 60 * as.numeric(stringr::str_extract(x, "\\d+")),
    TRUE ~ NA_real_
  )
}

assign_run_id <- function(df_img) {
  df_img %>%
    arrange(datetime) %>%
    mutate(
      dt_min = as.numeric(difftime(datetime, dplyr::lag(datetime), units = "mins")),
      new_run = dplyr::if_else(is.na(dt_min) | dt_min > 1.5 * frequence_min, 1L, 0L),
      run_seq = cumsum(new_run),
      run_id  = paste(bac, cycle, resolution, frequence, sprintf("run%02d", run_seq), sep = "_")
    ) %>%
    select(-dt_min, -new_run)
}

# -----------------------
# 2) GLMM utilities: ZI check + extraction of R2 & term tests
# -----------------------

fit_nb_or_zinb <- function(formula_nb, data, zi_alpha = 0.05, n_dharma = 2000,
                           quiet = TRUE) {
  m_nb <- glmmTMB::glmmTMB(formula_nb, family = glmmTMB::nbinom2(), data = data)
  sim <- DHARMa::simulateResiduals(m_nb, n = n_dharma, plot = FALSE)
  zi_test <- DHARMa::testZeroInflation(sim)

  if (!quiet) {
    message("Zero-inflation test p=", signif(zi_test$p.value, 3))
  }

  if (is.finite(zi_test$p.value) && zi_test$p.value < zi_alpha) {
    if (!quiet) message("Refitting with ZINB...")
    m_zi <- glmmTMB::glmmTMB(
      formula_nb,
      ziformula = ~ 1,
      family = glmmTMB::nbinom2(),
      data = data
    )
    return(list(model = m_zi, family = "ZINB", zi_p = zi_test$p.value))
  }

  list(model = m_nb, family = "NB", zi_p = zi_test$p.value)
}

# Likelihood-ratio test helper
lrt_models <- function(m_full, m_red) {
  a <- suppressWarnings(anova(m_full, m_red))
  # glmmTMB returns a data.frame-like with Chisq, Df, Pr(>Chisq)
  out <- tibble::tibble(
    chisq = as.numeric(a$Chisq[2]),
    df    = as.numeric(a$Df[2]),
    p     = as.numeric(a$`Pr(>Chisq)`[2])
  )
  out
}

# Extract explanatory power (R2) + term LRTs
extract_glmm_terms <- function(model, label = "model") {
  
  out <- list(label = label)
  
  # ---------- R2 ----------
  out$R2_marginal <- NA_real_
  out$R2_conditional <- NA_real_
  
  if (inherits(model, "merMod")) {
    r2 <- suppressWarnings(MuMIn::r.squaredGLMM(model))
    out$R2_marginal <- unname(r2[1, "R2m"])
    out$R2_conditional <- unname(r2[1, "R2c"])
  } else {
    # fallback (peut warning selon famille)
    r2 <- suppressWarnings(performance::r2_nakagawa(model))
    out$R2_marginal <- r2$R2_marginal
    out$R2_conditional <- r2$R2_conditional
  }
  
  # ---------- P-values type II ----------
  # (interaction testée explicitement, et effets principaux "globaux")
  a2 <- suppressWarnings(car::Anova(model, type = 2))
  a2_tab <- as.data.frame(a2)
  
  # helper
  grab_p <- function(term) {
    if (!term %in% rownames(a2_tab)) return(NA_real_)
    pv_col <- grep("^Pr\\(", colnames(a2_tab), value = TRUE)
    if (length(pv_col) == 0) return(NA_real_)
    as.numeric(a2_tab[term, pv_col[1]])
  }
  
  out$p_resolution <- grab_p("resolution")
  out$p_frequence  <- grab_p("frequence")
  out$p_interaction <- grab_p("resolution:frequence")
  
  # ---------- AIC/BIC ----------
  out$AIC <- AIC(model)
  out$BIC <- BIC(model)
  
  tibble::as_tibble(out)
}


# -----------------------
# 3) Cumsum curves + candidate curve fitting
# -----------------------

make_cumsum_df <- function(img_counts, group_cols, value_col = "n_img") {
  # img_counts must include datetime and the grouping columns
  img_counts %>%
    group_by(across(all_of(group_cols))) %>%
    arrange(datetime, .by_group = TRUE) %>%
    mutate(
      effort_img = dplyr::row_number(),
      cum_obs    = cumsum(.data[[value_col]])
    ) %>%
    ungroup()
}

safe_nls <- function(formula, data, start, lower = NULL, upper = NULL) {
  tryCatch({
    if (is.null(lower) && is.null(upper)) {
      nls(formula, data = data, start = start)
    } else {
      nls(formula, data = data, start = start, algorithm = "port",
          lower = lower %||% rep(-Inf, length(start)),
          upper = upper %||% rep(Inf, length(start)))
    }
  }, error = function(e) NULL)
}

`%||%` <- function(x, y) if (is.null(x)) y else x

fit_curve_candidates <- function(df_xy) {
  # df_xy must have x and y columns (x>0, y>=0)
  df_xy <- df_xy %>%
    mutate(x = as.numeric(x), y = as.numeric(y)) %>%
    filter(is.finite(x), is.finite(y), x > 0)

  if (nrow(df_xy) < 8) return(tibble())

  y_max <- max(df_xy$y, na.rm = TRUE)
  x_mid <- median(df_xy$x, na.rm = TRUE)

  fits <- list()

  # 1) Linear
  m_lin <- tryCatch(lm(y ~ x, data = df_xy), error = function(e) NULL)
  fits$linear <- m_lin

  # 2) Logarithmic: y = a + b*log(x)
  m_log <- tryCatch(lm(y ~ log(x), data = df_xy), error = function(e) NULL)
  fits$logarithmic <- m_log

  # 3) Power: y = a*x^b  (fit on log scale; requires y>0)
  if (all(df_xy$y > 0)) {
    m_pow <- tryCatch(lm(log(y) ~ log(x), data = df_xy), error = function(e) NULL)
    fits$power_loglog <- m_pow
  } else {
    fits$power_loglog <- NULL
  }

  # 4) Saturating exponential: y = a*(1 - exp(-b*x))
  m_exp <- safe_nls(
    y ~ a * (1 - exp(-b * x)),
    data = df_xy,
    start = list(a = y_max, b = 0.01),
    lower = c(a = 0, b = 1e-8)
  )
  fits$asymp_exponential <- m_exp

  # 5) Michaelis-Menten: y = (a*x)/(b + x)
  m_mm <- safe_nls(
    y ~ (a * x) / (b + x),
    data = df_xy,
    start = list(a = y_max, b = x_mid),
    lower = c(a = 0, b = 1e-8)
  )
  fits$michaelis_menten <- m_mm

  # 6) Logistic: y = a / (1 + exp(-(x - b)/c))
  m_logis <- safe_nls(
    y ~ a / (1 + exp(-(x - b)/c)),
    data = df_xy,
    start = list(a = y_max, b = x_mid, c = max(1, x_mid/5)),
    lower = c(a = 0, b = 1e-8, c = 1e-8)
  )
  fits$logistic <- m_logis

  # Build a tidy table
  out <- purrr::imap_dfr(fits, function(m, model_type) {
    if (is.null(m)) return(NULL)

    if (inherits(m, "lm")) {
      aic <- AIC(m)
      co  <- coef(m)

      if (model_type == "power_loglog") {
        # log(y)=alpha + beta*log(x) => y = exp(alpha) * x^beta
        tibble(
          curve_type = "power",
          AIC = aic,
          param_a = unname(exp(co[[1]])),
          param_b = unname(co[[2]]),
          param_c = NA_real_
        )
      } else if (model_type == "logarithmic") {
        tibble(
          curve_type = "logarithmic",
          AIC = aic,
          param_a = unname(co[[1]]),
          param_b = unname(co[[2]]),
          param_c = NA_real_
        )
      } else {
        tibble(
          curve_type = "linear",
          AIC = aic,
          param_a = unname(co[[1]]),
          param_b = unname(co[[2]]),
          param_c = NA_real_
        )
      }

    } else if (inherits(m, "nls")) {
      aic <- AIC(m)
      co  <- coef(m)

      # Harmonise parameter columns (a,b,c)
      param_a <- if ("a" %in% names(co)) unname(co[["a"]]) else NA_real_
      param_b <- if ("b" %in% names(co)) unname(co[["b"]]) else NA_real_
      param_c <- if ("c" %in% names(co)) unname(co[["c"]]) else NA_real_

      tibble(
        curve_type = model_type,
        AIC = aic,
        param_a = param_a,
        param_b = param_b,
        param_c = param_c
      )
    } else {
      NULL
    }
  })

  out %>% arrange(AIC)
}

# -----------------------
# 4) Read / prepare data
# -----------------------

# >>> Set your input file here
in_file <- "images_acquisition/database_calibration.csv"

df <- readr::read_csv(in_file, show_col_types = FALSE)

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

dat1 <- dat0 %>%
  mutate(
    class = as.character(class),
    class = dplyr::na_if(class, ""),
    class_clean = tolower(class),
    is_fauna = !is.na(class_clean) & !class_clean %in% c("na","<na>","root","som","unknown"),
    taxon = dplyr::case_when(
      class_clean %in% c("enchytraeidae","enchytraeid","enchytraeids") ~ "Enchytraeidae",
      class_clean %in% c("collembola","collembolan","collembolans")    ~ "Collembola",
      class_clean %in% c("acari","mesostigmata","oribatida","prostigmata","astigmata") ~ "Acari",
      TRUE ~ "Other"
    )
  )

# -----------------------
# 5) IMAGE-level counts
# -----------------------

img_total <- dat1 %>%
  filter(is_fauna) %>%
  group_by(bac, cycle, resolution, frequence, frequence_min, datetime, initial_image_name) %>%
  summarise(n_img = dplyr::n(), .groups = "drop") %>%
  group_by(bac, cycle, resolution, frequence) %>%
  assign_run_id() %>%
  ungroup()

img_taxa <- dat1 %>%
  filter(is_fauna, taxon %in% c("Enchytraeidae","Acari","Collembola")) %>%
  group_by(bac, cycle, resolution, frequence, frequence_min, taxon, datetime, initial_image_name) %>%
  summarise(n_img = dplyr::n(), .groups = "drop") %>%
  group_by(bac, cycle, resolution, frequence, taxon) %>%
  assign_run_id() %>%
  ungroup()

# -----------------------
# 6) RUN-level replicates
# -----------------------

run_total <- img_total %>%
  group_by(bac, cycle, resolution, frequence, run_id) %>%
  summarise(
    sum_n    = sum(n_img),
    mean_n   = mean(n_img),
    n_images = n(),
    start_time = min(datetime),
    end_time   = max(datetime),
    .groups = "drop"
  ) %>%
  mutate(
    log_effort = log(n_images),
    hour = lubridate::hour(start_time),
    period = dplyr::case_when(hour >= 6 & hour < 18 ~ "day", TRUE ~ "night"),
    period = factor(period)
  )

run_taxa <- img_taxa %>%
  group_by(bac, cycle, resolution, frequence, taxon, run_id) %>%
  summarise(
    sum_n    = sum(n_img),
    mean_n   = mean(n_img),
    n_images = n(),
    start_time = min(datetime),
    end_time   = max(datetime),
    .groups = "drop"
  ) %>%
  mutate(
    log_effort = log(n_images),
    hour = lubridate::hour(start_time),
    period = dplyr::case_when(hour >= 6 & hour < 18 ~ "day", TRUE ~ "night"),
    period = factor(period),
    taxon = factor(taxon)
  )

# -----------------------
# 7) GLMMs (total + taxa)
# -----------------------

# Use run-level sums with offset(log(n_images))
#form_run <- sum_n ~ resolution * frequence * cycle + offset(log_effort) + (1|bac) 
form_run <- sum_n ~ resolution * frequence + offset(log_effort) + (1|bac) + (1|cycle)

res_total <- fit_nb_or_zinb(form_run, run_total, quiet = TRUE)
res_ench  <- fit_nb_or_zinb(form_run, run_taxa %>% filter(taxon == "Enchytraeidae"), quiet = TRUE)
res_acar  <- fit_nb_or_zinb(form_run, run_taxa %>% filter(taxon == "Acari"), quiet = TRUE)
res_coll  <- fit_nb_or_zinb(form_run, run_taxa %>% filter(taxon == "Collembola"), quiet = TRUE)
# modèle collembole ne converge pas
# 1) Refit Collembola sans zero-inflation
res_coll$model <- glmmTMB::glmmTMB(
  form_run,
  family    = glmmTMB::nbinom2(link="log"),
  ziformula = ~0,
  data      = run_taxa %>% filter(taxon == "Collembola")
)
# 2) Vérifs
logLik(res_coll$model)
AIC(res_coll$model); BIC(res_coll$model)
summary(res_coll$model)


# -----------------------
# 8) Extract explanatory power + p-values (LRT)
# -----------------------

Table_terms <- bind_rows(
  extract_glmm_terms(res_total$model, "Total"),
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


safe_num <- function(x) ifelse(is.finite(x), x, NA_real_)

extract_model_quality <- function(model, label){
  # R2
  r2 <- suppressWarnings(performance::r2_nakagawa(model))
  R2m <- safe_num(r2$R2_marginal)
  R2c <- safe_num(r2$R2_conditional)
  
  # RE SD
  vc <- suppressWarnings(VarCorr(model))
  get_sd <- function(grp){
    if (!grp %in% names(vc$cond)) return(NA_real_)
    as.numeric(attr(vc$cond[[grp]], "stddev")[1])
  }
  sd_cycle <- get_sd("cycle")
  sd_bac   <- get_sd("bac")
  
  # likelihood + IC
  ll <- tryCatch(as.numeric(stats::logLik(model)), error=function(e) NA_real_)
  aic <- tryCatch(AIC(model), error=function(e) NA_real_)
  bic <- tryCatch(BIC(model), error=function(e) NA_real_)
  
  note <- ifelse(is.na(ll), "logLik undefined (model degeneracy) → AIC/BIC not reported", "")
  
  tibble::tibble(
    label = label,
    R2_marginal = round(R2m, 3),
    R2_conditional = round(R2c, 3),
    sd_cycle = round(sd_cycle, 3),
    sd_bac = round(sd_bac, 3),
    logLik = ifelse(is.na(ll), NA_character_, format(round(ll, 1), nsmall=1)),
    AIC = ifelse(is.na(aic), ".", format(round(aic, 1), nsmall=1)),
    BIC = ifelse(is.na(bic), ".", format(round(bic, 1), nsmall=1)),
    note = note
  )
}

p_stars <- function(p){
  ifelse(is.na(p), "",
         ifelse(p < 0.001, "***",
                ifelse(p < 0.01, "**",
                       ifelse(p < 0.05, "*", "ns"))))
}
fmt_p <- function(p){
  ifelse(is.na(p), NA_character_,
         ifelse(p < 1e-3, format(p, scientific = TRUE, digits = 2),
                format(round(p, 3), nsmall = 3)))
}

Table_fixed_effects <- Table_terms %>%
  mutate(
    term = recode(term,
                  p_resolution  = "Resolution",
                  p_frequence   = "Frequency",
                  p_interaction = "Resolution × Frequency"),
    cell = paste0(fmt_p(p_value), " ", p_stars(p_value))
  ) %>%
  select(label, term, cell) %>%
  tidyr::pivot_wider(names_from = term, values_from = cell) %>%
  arrange(factor(label, levels=c("Total","Enchytraeidae","Acari","Collembola")))

readr::write_csv(Table_terms, "images_acquisition/Tables/Table_fixed_effects.csv")


# -----------------------
# 9) Ratio table (resolution gains) — adapted from your second script
# -----------------------

get_ratio_table <- function(model, taxon_label,
                            contrasts = list(
                              "1200/600"   = c(-1, 1, 0),
                              "2400/600"   = c(-1, 0, 1),
                              "2400/1200"  = c(0, -1, 1)
                            )) {
  emm <- emmeans::emmeans(model, ~ resolution | frequence, type = "response")
  rat <- emmeans::contrast(emm, method = contrasts, by = "frequence", ratio = TRUE)

  as.data.frame(summary(rat)) %>%
    transmute(
      taxon = taxon_label,
      frequence = as.character(frequence),
      contrast = as.character(contrast),
      ratio = ratio,
      SE = SE
    )
}

Table_ratios <- bind_rows(
  get_ratio_table(res_ench$model, "Enchytraeidae"),
  get_ratio_table(res_acar$model, "Acari"),
  get_ratio_table(res_coll$model, "Collembola")
) %>%
  mutate(
    ratio = round(ratio, 2),
    SE    = round(SE, 2)
  ) %>%
  arrange(taxon, frequence, contrast)

readr::write_csv(Table_ratios, "images_acquisition/Tables/Table_resolution_ratios_by_taxon.csv")

# -----------------------
# 9b) Figures — GLMM effects (restoring original paper figures)
# -----------------------

dir.create("images_acquisition/Figures", recursive = TRUE, showWarnings = FALSE)

# ---- Total fauna: EMM plot (expected count per image)
emm_total <- emmeans::emmeans(res_total$model, ~ resolution * frequence, type = "response")
emm_total_tbl <- as.data.frame(emm_total) %>%
  dplyr::rename(lower.CL = asymp.LCL, upper.CL = asymp.UCL) %>%
  dplyr::filter(!is.na(lower.CL), !is.na(upper.CL)) %>%
  dplyr::mutate(
    resolution = as.factor(resolution),
    frequence  = as.factor(frequence)
  )

p_total <- emm_total_tbl %>%
  ggplot(aes(x = frequence, y = response, colour = resolution)) +
  geom_point(position = position_dodge(width = 0.2)) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                width = 0.1, position = position_dodge(width = 0.2)) +
  labs(
    x = "Acquisition interval (min)",
    y = "Detected activity (expected count per image)",
    title = "Effects of resolution and acquisition frequency (total fauna)"
  ) +
  theme_bw()

print(p_total)

ggplot2::ggsave(
  filename = "images_acquisition/Figures/Fig_total_activity.png",
  plot     = p_total,
  width    = 180,
  height   = 120,
  units    = "mm",
  dpi      = 300,
  bg       = "white"
)

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
  ggplot(aes(x = frequence, y = response, colour = resolution)) +
  geom_point(position = position_dodge(width = 0.2)) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                width = 0.1, position = position_dodge(width = 0.2)) +
  facet_grid(taxon ~ ., scales = "free_y") +
  labs(
    x = "Acquisition interval (min)",
    y = "Detected activity (expected count per image)"
  ) +
  theme_bw()

print(p_taxa)

ggplot2::ggsave(
  filename = "images_acquisition/Figures/Fig_taxa_activity.png",
  plot     = p_taxa,
  width    = 180,
  height   = 120,
  units    = "mm",
  dpi      = 300,
  bg       = "white"
)

# -----------------------
# 9c) shifts in apparent size distributions (quantiles)
# -----------------------

parse_dpi <- function(x) {
  as.numeric(stringr::str_extract(tolower(as.character(x)), "\\d+"))
}

# Safety: only run if bbox columns exist
bbox_cols <- c("xmin","xmax","ymin","ymax")
if (all(bbox_cols %in% names(dat1))) {

  det_sizes <- dat1 %>%
    filter(taxon %in% c("Enchytraeidae","Acari","Collembola")) %>%
    mutate(
      dpi = parse_dpi(resolution),
      width_px  = as.numeric(xmax) - as.numeric(xmin),
      height_px = as.numeric(ymax) - as.numeric(ymin),
      area_px   = width_px * height_px,
      size_px   = sqrt(area_px),
      px_to_mm  = 2.54 / dpi,
      size_mm   = area_px * px_to_mm
    ) %>%
    filter(!is.na(dpi), is.finite(size_mm), width_px > 0, height_px > 0, size_mm > 0)

  det_q <- det_sizes %>%
    group_by(taxon, resolution, frequence) %>%
    summarise(
      q25 = quantile(size_mm, 0.25, na.rm = TRUE),
      q50 = quantile(size_mm, 0.50, na.rm = TRUE),
      q75 = quantile(size_mm, 0.75, na.rm = TRUE),
      n   = n(),
      .groups = "drop"
    )

  readr::write_csv(det_q, "images_acquisition/Tables/Table_Sx_size_quantiles_q25_q50_q75.csv")

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
      colour = "Quantile"
    ) +
    theme_bw() +
    theme(legend.position = "bottom") +
    scale_colour_manual(
      values = c(q25 = "#deebf7", q50 = "#9ecae1", q75 = "#3182bd"),
      labels = c(q25 = "25th percentile", q50 = "Median", q75 = "75th percentile")
    )

  ggplot2::ggsave(
    "images_acquisition/Figures/SuppFig_Sx_size_quantiles_by_taxon.png",
    p_quant_taxon,
    width = 220, height = 180, units = "mm", dpi = 300, bg = "white"
  )

} else {
  message("Skipping size-quantiles figure: bbox columns not found in input data (xmin/xmax/ymin/ymax).")
}


# -----------------------
# 10) Cumsum curves for each resolution × frequency (total + 3 taxa)
# -----------------------

# Total
cum_total <- make_cumsum_df(
  img_total,
  group_cols = c("resolution","frequence"),
  value_col = "n_img"
) %>%
  transmute(
    taxon = "Total",
    resolution,
    frequence,
    datetime,
    x = effort_img,
    y = cum_obs
  )

# Taxa
cum_taxa <- make_cumsum_df(
  img_taxa,
  group_cols = c("taxon","resolution","frequence"),
  value_col = "n_img"
) %>%
  transmute(
    taxon = as.character(taxon),
    resolution,
    frequence,
    datetime,
    x = effort_img,
    y = cum_obs
  )

cum_all <- bind_rows(cum_total, cum_taxa) %>%
  mutate(
    taxon = factor(taxon, levels = c("Total","Acari","Collembola","Enchytraeidae")),
    resolution = factor(resolution),
    frequence  = factor(frequence)
  )

readr::write_csv(cum_all, "images_acquisition/Tables/Table_cumsum_curves_raw.csv")

# -----------------------
# 11) Fit curve types to each cumsum curve
# -----------------------

Table_curves_all <- cum_all %>%
  group_by(taxon, resolution, frequence) %>%
  group_modify(~{
    fits <- fit_curve_candidates(.x %>% select(x, y))
    if (nrow(fits) == 0) return(tibble())
    fits
  }) %>%
  ungroup() %>%
  arrange(taxon, resolution, frequence, AIC)

Table_curves_best <- Table_curves_all %>%
  group_by(taxon, resolution, frequence) %>%
  slice_min(order_by = AIC, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(taxon, resolution, frequence)

readr::write_csv(Table_curves_all,  "images_acquisition/Tables/Table_cumsum_curve_fits_ALL.csv")
readr::write_csv(Table_curves_best, "images_acquisition/Tables/Table_cumsum_curve_fits_BEST.csv")


# -----------------------
# 12) Graphical outputs (GLMM diagnostics + cumsum curves)
# -----------------------

fig_dir <- "images_acquisition/Figures"
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

# --- 12.1 GLMM diagnostics (DHARMa) ---
save_dharma <- function(model, name, n = 2000) {
  sim <- DHARMa::simulateResiduals(model, n = n, plot = FALSE)
  out_pdf <- file.path(fig_dir, paste0("DHARMa_", name, ".pdf"))
  grDevices::pdf(out_pdf, width = 8, height = 7)
  plot(sim)
  grDevices::dev.off()
  invisible(out_pdf)
}

invisible(save_dharma(res_total$model, "Total"))
invisible(save_dharma(res_ench$model,  "Enchytraeidae"))
invisible(save_dharma(res_acar$model,  "Acari"))
invisible(save_dharma(res_coll$model,  "Collembola"))

# --- 12.2 Raw cumsum curves ---

p_raw <- ggplot(cum_all, aes(x = x, y = log(y), group = interaction(resolution, frequence), colour = resolution)) +
  geom_line(linewidth = 0.6) +
  facet_grid(taxon ~ frequence, scales = "free_y") +
  labs(x = "Number of events", y = "Cumulative detections", colour = "Resolution") +
  theme_bw()

ggsave(filename = file.path(fig_dir, "Cumsum_raw_by_taxon_freq.png"), plot = p_raw, width = 11, height = 8, dpi = 300)

# --- 12.3 Best-fit curves overlay ---

predict_curve <- function(curve_type, a, b, c, x){
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


# Build predictions for each group using the best (AIC-min) curve
pred_best <- cum_all %>%
  group_by(taxon, resolution, frequence) %>%
  summarise(x_min = min(x), x_max = max(x), .groups = "drop") %>%
  left_join(Table_curves_best, by = c("taxon","resolution","frequence")) %>%
  mutate(x_grid = purrr::map2(x_min, x_max, ~seq(.x, .y, length.out = 200))) %>%
  filter(!is.na(curve_type)) %>%   # <- clé
  mutate(
    y_hat = purrr::pmap(
      list(curve_type, param_a, param_b, param_c, x_grid),
      ~predict_curve(..1, ..2, ..3, ..4, ..5)
    )
  ) %>%
  select(taxon, resolution, frequence, curve_type, param_a, param_b, param_c, x_grid, y_hat) %>%
  tidyr::unnest(c(x_grid, y_hat)) %>%
  rename(x = x_grid, y = y_hat)


p_fit <- ggplot() +
  geom_line(data = cum_all,
            aes(x = x, y = log(y), group = interaction(resolution, frequence), colour = resolution),
            linewidth = 0.5, alpha = 0.35) +
  geom_line(data = pred_best,
            aes(x = x, y = log(y), group = interaction(resolution, frequence), colour = resolution),
            linewidth = 0.9) +
  facet_grid(taxon ~ frequence)+ #, scales = "free_y") +
  labs(x = "Number of events", y = "Cumulative detections (log scale)",
       colour = "Resolution",
       title = "Cumulative detections with best-fit curve (AIC)") +
  theme_bw()

ggsave(filename = file.path(fig_dir, "Cumsum_bestfit_overlay_by_taxon_freq.png"), plot = p_fit, width = 11, height = 8, dpi = 300)


# 2) Calcul y_hat(100) = a * 100^b
best2 <- Table_curves_best  %>%
  mutate(
    x_ref = 100,
    y_hat_100 = param_a * (x_ref ^ param_b)
  )

# 3) Tableau synthétique (large) : une ligne par taxon × résolution, colonnes = fréquences
tab_yhat100 <- best2 %>%
  select(taxon, resolution, frequence, y_hat_100) %>%
  mutate(
    resolution = as.character(resolution),
    frequence  = as.character(frequence)
  ) %>%
  tidyr::pivot_wider(
    names_from  = frequence,
    values_from = y_hat_100,
    names_prefix = "yhat100_f"
  ) %>%
  arrange(taxon, as.numeric(resolution)) %>%
  mutate(across(starts_with("yhat100_"), ~round(.x, 1)))

# 4) Export (optionnel)
readr::write_csv(best2, "images_acquisition/Tables/Table_cumsum_yhat100_allcombos.csv")
readr::write_csv(tab_yhat100, "images_acquisition/Tables/Table_cumsum_yhat100_summary.csv")
