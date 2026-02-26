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

# Extract beta coeff
extract_glmm_betas <- function(model, label = "model") {
  
  coefs <- as.data.frame(summary(model)$coefficients$cond)
  coefs$term <- rownames(coefs)
  
  coefs <- dplyr::rename(coefs,
                         estimate  = Estimate,
                         std.error = `Std. Error`,
                         p.value   = `Pr(>|z|)`) %>%
    dplyr::mutate(
      label = label,
      conf.low  = estimate - 1.96 * std.error,
      conf.high = estimate + 1.96 * std.error
    ) %>%
    # ✅ garde tous les coefficients liés à resolution/frequence/interaction
    dplyr::filter(grepl("^(resolution|frequence)|resolution.*:frequence", term)) %>%
    # ✅ enlève les NA non-estimables (chez toi: resolution2400:frequence360)
    dplyr::filter(is.finite(estimate), is.finite(std.error)) %>%
    dplyr::mutate(
      estimate  = round(estimate, 3),
      std.error = round(std.error, 3),
      conf.low  = round(conf.low, 3),
      conf.high = round(conf.high, 3),
      p.value   = round(p.value, 3)
    ) %>%
    dplyr::select(label, term, estimate, std.error, conf.low, conf.high, p.value)
  
  tibble::as_tibble(coefs)
}

# Compute Incidence Rate Ratio
extract_betas_irr <- function(model, label = "model") {
  
  coefs <- as.data.frame(summary(model)$coefficients$cond)
  coefs$term_raw <- rownames(coefs)
  
  coefs <- dplyr::rename(coefs,
                         estimate  = Estimate,
                         std.error = `Std. Error`) %>%
    dplyr::mutate(
      label = label,
      conf.low  = estimate - 1.96 * std.error,
      conf.high = estimate + 1.96 * std.error,
      IRR      = exp(estimate),
      IRR_low  = exp(conf.low),
      IRR_high = exp(conf.high)
    ) %>%
    dplyr::filter(grepl("^(resolution|frequence)|resolution.*:frequence", term_raw)) %>%
    dplyr::filter(is.finite(estimate), is.finite(std.error)) %>%
    dplyr::mutate(
      term = dplyr::case_when(
        grepl("^resolution", term_raw) ~ "Resolution",
        grepl("^frequence",  term_raw) ~ "Frequency",
        grepl("resolution.*:frequence", term_raw) ~ "Resolution × Frequency",
        TRUE ~ NA_character_
      ),
      estimate  = round(estimate, 3),
      std.error = round(std.error, 3),
      conf.low  = round(conf.low, 3),
      conf.high = round(conf.high, 3),
      IRR       = round(IRR, 3),
      IRR_low   = round(IRR_low, 3),
      IRR_high  = round(IRR_high, 3)
    ) %>%
    dplyr::select(label, term, term_raw, estimate, std.error, conf.low, conf.high, IRR, IRR_low, IRR_high)
  
  tibble::as_tibble(coefs)
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
# 8) Extract model outputs
# -----------------------
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
fmt_beta_ci <- function(est, lo, hi, digits = 3){
  est <- round(est, digits); lo <- round(lo, digits); hi <- round(hi, digits)
  paste0(est, " [", lo, ", ", hi, "]")
}
fmt_beta <- function(est, lo, hi) paste0(est, " [", lo, ", ", hi, "]")
fmt_irr  <- function(irr, lo, hi) paste0(irr, " [", lo, ", ", hi, "]")

collapse_terms <- function(df){
  df %>%
    dplyr::mutate(
      beta_ci = fmt_beta(estimate, conf.low, conf.high),
      irr_ci  = fmt_irr(IRR, IRR_low, IRR_high)
    ) %>%
    dplyr::group_by(label, term) %>%
    dplyr::summarise(
      beta_ci = paste0(term_raw, ": ", beta_ci, collapse = "; "),
      irr_ci  = paste0(term_raw, ": ", irr_ci,  collapse = "; "),
      .groups = "drop"
    )
}


Table_betas_cells <- dplyr::bind_rows(
  extract_glmm_betas(res_total$model, "Total"),
  extract_glmm_betas(res_ench$model,  "Enchytraeidae"),
  extract_glmm_betas(res_acar$model,  "Acari"),
  extract_glmm_betas(res_coll$model,  "Collembola")
) %>%
  dplyr::mutate(
    term = dplyr::recode(term,
                         "resolution"="Resolution",
                         "frequence"="Frequency",
                         "resolution:frequence"="Resolution × Frequency"),
    beta_ci = fmt_beta_ci(estimate, conf.low, conf.high)
  ) %>%
  dplyr::select(label, term, beta_ci) %>%
  tidyr::pivot_wider(names_from = term, values_from = beta_ci)

Betas_total <- collapse_terms(extract_betas_irr(res_total$model, "Total"))
Betas_ench  <- collapse_terms(extract_betas_irr(res_ench$model,  "Enchytraeidae"))
Betas_acar  <- collapse_terms(extract_betas_irr(res_acar$model,  "Acari"))
Betas_coll  <- collapse_terms(extract_betas_irr(res_coll$model,  "Collembola"))

Table_betas_irr <- dplyr::bind_rows(Betas_total, Betas_ench, Betas_acar, Betas_coll) %>%
  tidyr::pivot_wider(
    names_from = term,
    values_from = c(beta_ci, irr_ci)
  )

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

Table_p_typeII <- Table_terms %>%
  dplyr::mutate(
    term = dplyr::recode(term,
                         p_resolution  = "Resolution",
                         p_frequence   = "Frequency",
                         p_interaction = "Resolution × Frequency"),
    p_cell = paste0(fmt_p(p_value), " ", p_stars(p_value))
  ) %>%
  dplyr::select(label, term, p_cell) %>%
  tidyr::pivot_wider(names_from = term, values_from = p_cell, names_prefix = "p_")


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

Table_quality <- dplyr::bind_rows(
  extract_model_quality(res_total$model, "Total"),
  extract_model_quality(res_ench$model,  "Enchytraeidae"),
  extract_model_quality(res_acar$model,  "Acari"),
  extract_model_quality(res_coll$model,  "Collembola")
)

Table_model_terms <- Table_betas_irr %>%
  dplyr::left_join(Table_p_typeII, by = "label") %>%
  dplyr::arrange(factor(label, levels=c("Total","Enchytraeidae","Acari","Collembola")))

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
  dplyr::arrange(factor(label, levels=c("Total","Enchytraeidae","Acari","Collembola")))

# Option: jolis noms de colonnes
Table_main <- Table_main %>%
  dplyr::rename(
    `p (Resolution)` = Resolution,
    `p (Frequency)`  = Frequency,
    `p (Res × Freq)` = `Resolution × Frequency`,
    `R²m` = R2_marginal
  )


readr::write_csv(Table_model_terms, "images_acquisition/Tables/Table_model_summary.csv")

# ---- SI table: IRR per contrast (readable) ----

extract_contrast_irr <- function(model, label){
  
  mf <- model.frame(model)
  ref_res <- levels(mf$resolution)[1]
  ref_frq <- levels(mf$frequence)[1]
  
  coefs <- as.data.frame(summary(model)$coefficients$cond)
  coefs$term_raw <- rownames(coefs)
  
  coefs %>%
    dplyr::rename(estimate = Estimate, se = `Std. Error`) %>%
    dplyr::filter(grepl("^(resolution|frequence)|resolution.*:frequence", term_raw)) %>%
    dplyr::filter(is.finite(estimate), is.finite(se)) %>%  # retire les NA non estimables
    dplyr::mutate(
      label = label,
      lo = estimate - 1.96 * se,
      hi = estimate + 1.96 * se,
      IRR = exp(estimate),
      IRR_lo = exp(lo),
      IRR_hi = exp(hi),
      Class = dplyr::case_when(
        grepl("^resolution", term_raw) ~ "Resolution",
        grepl("^frequence",  term_raw) ~ "Frequency",
        TRUE ~ "Resolution × Frequency"
      ),
      Contrast = dplyr::case_when(
        grepl("^resolution", term_raw) ~ paste0(sub("^resolution","", term_raw), " vs ", ref_res),
        grepl("^frequence", term_raw)  ~ paste0(sub("^frequence","", term_raw), " vs ", ref_frq),
        TRUE ~ term_raw
      ),
      `IRR (95% CI)` = paste0(round(IRR, 2), " (", round(IRR_lo, 2), "–", round(IRR_hi, 2), ")")
    ) %>%
    dplyr::select(label, Class, Contrast, `IRR (95% CI)`)
}

Table_SI_irr <- dplyr::bind_rows(
  extract_contrast_irr(res_total$model, "Total"),
  extract_contrast_irr(res_ench$model,  "Enchytraeidae"),
  extract_contrast_irr(res_acar$model,  "Acari"),
  extract_contrast_irr(res_coll$model,  "Collembola")
) %>%
  dplyr::arrange(factor(label, levels=c("Total","Enchytraeidae","Acari","Collembola")), Class, Contrast)

readr::write_csv(Table_SI_irr, "images_acquisition/Tables/Table_SI_IRR_contrasts.csv")

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
=======
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

Table_quality <- dplyr::bind_rows(
  extract_model_quality(res_total$model, "Total"),
  extract_model_quality(res_ench$model,  "Enchytraeidae"),
  extract_model_quality(res_acar$model,  "Acari"),
  extract_model_quality(res_coll$model,  "Collembola")
)

Table_model_summary <- Table_betas_cells %>%
  dplyr::left_join(Table_p_cells, by = "label", suffix = c("_beta", "_p")) %>%
  dplyr::left_join(Table_quality, by = "label")

# Option: construire des colonnes combinées β[CI] (p)
combine_beta_p <- function(beta, p){
  ifelse(is.na(beta), NA_character_, paste0(beta, " (p=", p, ")"))
}

Table_model_summary <- Table_model_summary %>%
  dplyr::mutate(
    `Resolution` = combine_beta_p(Resolution_beta, Resolution_p),
    `Frequency`  = combine_beta_p(Frequency_beta,  Frequency_p),
    `Resolution × Frequency` = combine_beta_p(`Resolution × Frequency_beta`, `Resolution × Frequency_p`)
  ) %>%
  dplyr::select(
    label,
    `Resolution`, `Frequency`, `Resolution × Frequency`,
    R2_marginal, R2_conditional, sd_bac, sd_cycle, AIC, BIC
  ) %>%
  dplyr::arrange(factor(label, levels=c("Total","Enchytraeidae","Acari","Collembola")))

readr::write_csv(Table_model_summary, "images_acquisition/Tables/Table_model_summary.csv")



### A enlever ######

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
##################################

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
  ggplot(aes(x = resolution, y = log(response))) +
  geom_point(position = position_dodge(width = 0.2)) +
  geom_errorbar(aes(ymin = log(lower.CL), ymax = log(upper.CL)),
                width = 0.1, position = position_dodge(width = 0.2)) +
  facet_grid(taxon ~ frequence) +
  labs(
    x = "Image resolution (dpi)",
    y = "Detected individuals per image\n(log scale)"
  ) +
  theme_bw()

print(p_taxa)

ggplot2::ggsave(
  filename = "images_acquisition/Figures/Fig_taxa_activity.png",
  plot     = p_taxa,
  width    = 180,
  height   = 180,
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
#    + summaries (mean/sd across runs)
#    + curve fitting per run and summary across runs
# -----------------------

# ---- 10.1 Build raw cumulative curves (per run_id) ----
# Total
cum_total <- make_cumsum_df(
  img_total,
  group_cols = c("resolution","frequence","run_id"),
  value_col  = "n_img"
) %>%
  transmute(
    taxon      = "Total",
    resolution = factor(resolution),
    frequence  = factor(frequence),
    datetime,
    run_id,
    x = effort_img,   # 1..12
    y = cum_obs
  )

# Taxa
cum_taxa <- make_cumsum_df(
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
  )

cum_all <- bind_rows(cum_total, cum_taxa) %>%
  mutate(
    taxon = factor(taxon, levels = c("Total","Acari","Collembola","Enchytraeidae"))
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

invisible(save_dharma(res_total$model, "Total"))
invisible(save_dharma(res_ench$model,  "Enchytraeidae"))
invisible(save_dharma(res_acar$model,  "Acari"))
invisible(save_dharma(res_coll$model,  "Collembola"))

# ---- 12.2 Raw curves (per run) + mean ± SD (exclude Total if desired) ----
p_raw_mean <- ggplot() +
  # raw run curves (thin)
  geom_line(
    data = subset(cum_all, taxon != "Total"),
    aes(x = x, y = log(y), group = interaction(run_id, resolution, frequence), colour = resolution),
    linewidth = 0.4, alpha = 0.25
  ) +
  # mean ± SE ribbon
  geom_ribbon(
    data = subset(cum_all_sum, taxon != "Total"),
    aes(
      x = x,
      ymin = log(pmax(y_mean - y_se, 1e-6)),
      ymax = log(pmax(y_mean + y_se, 1e-6)),
      group = interaction(resolution, frequence)
    ),
    alpha = 0.15
  )+
  # mean curve
  geom_line(
    data = subset(cum_all_sum, taxon != "Total"),
    aes(x = x, y = log(pmax(y_mean, 1e-9)), group = interaction(resolution, frequence), colour = resolution),
    linewidth = 0.9
  ) +
  facet_grid(taxon ~ frequence) +
  labs(x = "Number of events (images)", y = "Cumulative detections (log scale)", colour = "Resolution") +
  theme_bw()

ggsave(
  filename = file.path(fig_dir, "Cumsum_raw_plus_mean_sd_by_taxon_freq.png"),
  plot = p_raw_mean, width = 11, height = 8, dpi = 300
)

# ---- 12.2 Best-fit curves overlay (using mode curve_type per group) ----
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

