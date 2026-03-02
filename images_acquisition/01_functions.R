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
