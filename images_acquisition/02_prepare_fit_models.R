# -----------------------
# 1) Read / prepare data
# -----------------------

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
      class_clean %in% c("acari","mesostigmata","oribatida","prostigmata","astigmata") ~ "Acari",
      TRUE ~ "Other"
    )
  )

# -----------------------
# 2) IMAGE-level counts
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
# 3) RUN-level replicates
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
# 4) GLMMs 
# -----------------------
form_run <- sum_n ~ resolution * frequence + offset(log_effort) + (1|bac) + (1|cycle)
form_run2 <- sum_n ~ resolution * frequence + offset(log_effort) + (1|cycle)

res_total <- fit_nb_or_zinb(form_run, run_total, quiet = TRUE)
res_ench  <- fit_nb_or_zinb(form_run, run_taxa %>% filter(taxon == "Enchytraeidae"), quiet = TRUE)
res_acar  <- fit_nb_or_zinb(form_run, run_taxa %>% filter(taxon == "Acari"), quiet = TRUE)
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