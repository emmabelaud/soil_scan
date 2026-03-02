# -----------------------
# 0) Packages
# -----------------------
pkgs <- c(
  "dplyr","tidyr","stringr","lubridate","forcats","purrr","tibble",
  "glmmTMB","DHARMa","emmeans","ggplot2","performance","readr",
  "car","MuMIn"
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


# Extract model ouputs
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
    r2 <- suppressWarnings(performance::r2_nakagawa(model))
    r2_df <- tryCatch(as.data.frame(r2), error = function(e) NULL)
    
    if (!is.null(r2_df) && all(c("R2_marginal","R2_conditional") %in% names(r2_df))) {
      out$R2_marginal <- as.numeric(r2_df$R2_marginal[1])
      out$R2_conditional <- as.numeric(r2_df$R2_conditional[1])
    } else if (is.numeric(r2) && length(r2) >= 2) {
      out$R2_marginal <- unname(r2[1])
      out$R2_conditional <- unname(r2[2])
    } else {
      # cas NA (singularité / variance RE=0)
      out$R2_marginal <- NA_real_
      out$R2_conditional <- NA_real_
    }
  }
  
  # ---------- P-values type II ----------
  a2 <- suppressWarnings(car::Anova(model, type = 2))
  a2_tab <- as.data.frame(a2)
  
  grab_p <- function(term) {
    if (!term %in% rownames(a2_tab)) return(NA_real_)
    pv_col <- grep("^Pr\\(", colnames(a2_tab), value = TRUE)
    if (length(pv_col) == 0) return(NA_real_)
    as.numeric(a2_tab[term, pv_col[1]])
  }
  
  out$p_resolution   <- grab_p("resolution")
  out$p_frequence    <- grab_p("frequence")
  out$p_interaction  <- grab_p("resolution:frequence")
  
  # ---------- AIC/BIC ----------
  out$AIC <- AIC(model)
  out$BIC <- BIC(model)
  
  tibble::as_tibble(out)
}

extract_model_quality <- function(model, label){
  
  safe_num <- function(x) ifelse(is.finite(as.numeric(x)), as.numeric(x), NA_real_)
  
  # R2 (robuste)
  r2 <- suppressWarnings(performance::r2_nakagawa(model))
  r2_df <- tryCatch(as.data.frame(r2), error = function(e) NULL)
  
  R2m <- NA_real_
  R2c <- NA_real_
  if (!is.null(r2_df) && all(c("R2_marginal","R2_conditional") %in% names(r2_df))) {
    R2m <- safe_num(r2_df$R2_marginal[1])
    R2c <- safe_num(r2_df$R2_conditional[1])
  } else if (is.numeric(r2) && length(r2) >= 2) {
    R2m <- safe_num(r2[1])
    R2c <- safe_num(r2[2])
  }
  
  # RE SD (glmmTMB)
  vc <- suppressWarnings(VarCorr(model))
  get_sd <- function(grp){
    if (!grp %in% names(vc$cond)) return(NA_real_)
    as.numeric(attr(vc$cond[[grp]], "stddev")[1])
  }
  sd_cycle <- get_sd("cycle")
  sd_bac   <- get_sd("bac")
  
  # IC
  ll  <- tryCatch(as.numeric(stats::logLik(model)), error=function(e) NA_real_)
  aic <- tryCatch(AIC(model), error=function(e) NA_real_)
  bic <- tryCatch(BIC(model), error=function(e) NA_real_)
  
  tibble::tibble(
    label = label,
    R2_marginal = round(R2m, 3),
    R2_conditional = round(R2c, 3),
    sd_cycle = round(sd_cycle, 3),
    sd_bac   = round(sd_bac, 3),
    AIC = round(aic, 1),
    BIC = round(bic, 1),
    converged = tryCatch(isTRUE(model$sdr$pdHess), error=function(e) NA),
    singular  = tryCatch(isTRUE(performance::check_singularity(model)), error=function(e) NA)
  )
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

# input file 
in_file <- "images_acquisition/data/database_calibration.csv"

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
      class_clean %in% c("acari","mesostigmata","oribatida","prostigmata","astigmata","trombidiformes","bdelloidea") ~ "Acari",
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
# 7) GLMMs (taxa)
# -----------------------

# Use run-level sums with offset(log(n_images))
form_run <- sum_n ~ resolution * frequence + offset(log_effort) + (1|bac) + (1|cycle)
form_run2 <- sum_n ~ resolution * frequence + offset(log_effort) + (1|cycle)

res_ench  <- fit_nb_or_zinb(form_run, run_taxa %>% filter(taxon == "Enchytraeidae"), quiet = TRUE)
res_acar  <- fit_nb_or_zinb(form_run2, run_taxa %>% filter(taxon == "Acari"), quiet = TRUE)
res_coll  <- fit_nb_or_zinb(form_run, run_taxa %>% filter(taxon == "Collembola"), quiet = TRUE)
# modèle collembole ne converge pas
# 1) Refit Collembola sans zero-inflation
res_coll$model <- glmmTMB::glmmTMB(
  form_run,
  family    = glmmTMB::nbinom2(link="log"),
  ziformula = ~0,
  data      = run_taxa %>% filter(taxon == "Collembola")
)



# -----------------------
# 8) Extract model outputs
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
# 9) Ratio table (resolution gains)
# -----------------------
get_ratio_table <- function(model, taxon_label){
  
  emm <- emmeans::emmeans(model, ~ resolution | frequence, type = "response")
  
  # ratio 2400/600, 1200/600, 2400/1200 (si niveaux présents dans ce frequence)
  rat <- emmeans::contrast(
    emm,
    method = list(
      "1200/600"  = c(-1, 1, 0),
      "2400/600"  = c(-1, 0, 1),
      "2400/1200" = c(0, -1, 1)
    ),
    by = "frequence",
    ratio = TRUE
  )
  
  as.data.frame(summary(rat, infer = c(TRUE, TRUE))) %>%  # CI
    transmute(
      taxon = taxon_label,
      frequence = as.character(frequence),
      contrast = as.character(contrast),
      ratio = ratio,
      low  = asymp.LCL,
      high = asymp.UCL
    ) %>%
    mutate(
      ratio_ci = sprintf("%.2f× [%.2f, %.2f]", ratio, low, high)
    )
}
Table_ratios <- bind_rows(
  get_ratio_table(res_ench$model, "Enchytraeidae"),
  get_ratio_table(res_acar$model, "Acari"),
  get_ratio_table(res_coll$model, "Collembola")
) %>%
  arrange(taxon, frequence, contrast)

readr::write_csv(Table_ratios, "images_acquisition/Tables/Table_resolution_ratios_by_taxon.csv")

# -----------------------
# 9b) Figures — GLMM effects
# -----------------------

dir.create("images_acquisition/Figures", recursive = TRUE, showWarnings = FALSE)

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
  geom_point(size = 2) +
  geom_line(linewidth = 0.9) +
  scale_y_log10() +
  facet_grid(taxon ~ frequence) +
  labs(
    x = "Spatial resolution (dpi)",
    y = "Apparent area (mm², log scale)",
    colour = "Quantile"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")+
  ggplot2::scale_colour_manual(
    values = c(
      "25th percentile" = "#deebf7",
      "Median"          = "#9ecae1",
      "75th percentile" = "#3182bd"
    )
  )

print(p_quant_taxon)

ggplot2::ggsave(
  "images_acquisition/Figures/SuppFig_Sx_size_quantiles_by_taxon.png",
  p_quant_taxon,
  width = 220, height = 180, units = "mm", dpi = 300, bg = "white"
)


# -----------------------
# 10) Cumsum curves for each resolution × frequency
#    + summaries (mean/sd across runs)
#    + curve fitting per run and summary across runs
# -----------------------

# ---- 10.1 Build raw cumulative curves (per run_id) ----
# Total
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

# ---- 10.2 Summary mean/sd across runs (for plotting mean curves) ----
# Note: averages at fixed x (= image count), across runs within taxon×resolution×frequence
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
# 11) Fit curve types to each cumsum curve (per run), then summarise across runs
# -----------------------

# ---- 11.1 Fit candidate curves per run ----
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

# ---- 11.2 Pick best curve per run (AIC-min) ----
Table_curves_best <- Table_curves_all %>%
  group_by(taxon, resolution, frequence, run_id) %>%
  slice_min(order_by = AIC, n = 1, with_ties = FALSE) %>%
  ungroup()

readr::write_csv(Table_curves_best, "images_acquisition/Tables/Table_cumsum_curve_fits_BEST_by_run.csv")

# ---- 11.3 Summarise best-fit parameters across runs (mean/sd) ----
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

invisible(save_dharma(res_ench$model,  "Enchytraeidae"))
invisible(save_dharma(res_acar$model,  "Acari"))
invisible(save_dharma(res_coll$model,  "Collembola"))


# --- 12.3 Best-fit curves overlay ---
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
    linewidth = 0.4, alpha = 0.25
  ) +
  geom_line(
    data = pred_best,
    aes(x = x, y = log(pmax(y, 1e-9)), group = interaction(resolution, frequence), colour = resolution),
    linewidth = 1.0
  ) +
  facet_grid(taxon ~ frequence) +
  labs(x = "Number of events (images)", y = "Cumulative detections (log scale)", colour = "Resolution") +
  theme_bw()+
  scale_x_continuous(limits = c(1, 11), breaks = c(3,6,9))

ggsave(
  filename = file.path(fig_dir, "Cumsum_bestfit_overlay_by_taxon_freq.png"),
  plot = p_fit, width = 11, height = 8, dpi = 300
)

ggsave(filename = file.path(fig_dir, "Cumsum_bestfit_overlay_by_taxon_freq.png"), plot = p_fit, width = 11, height = 8, dpi = 300)


