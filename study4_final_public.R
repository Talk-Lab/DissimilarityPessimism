# Title: Dissimilarity Pessimism Study 4
# Author(s): Gus Cooney; Erica Boothby
# Description: Exploring mechanism for dissimilarity pessimism

# NOTES -------------------------------------------------------------------


# load libraries ----------------------------------------------------------

library(pacman)
p_load(tidyverse, glue, broom, lme4, nlme, sjPlot, lmerTest, emmeans, lavaan, car, RSA, tkrplot, qgraph, rgl, psych, rempsyc,
       ltm, rstatix,ggplot2, ggpubr, boot, boot.pval, JSmediation, mediation, bruceR, fwb, syuzhet, openai, jsonlite, purrr, 
       rmcorr, patchwork)

# set constants -----------------------------------------------------------

DATA_PATH <- "~/Dropbox/Manuscripts/Homophily/Data"
PLOT_PATH <- "~/Dropbox/Manuscripts/Homophily/Plots"

emm_options(pbkrtest.limit = 10000)

# load data ---------------------------------------------------------------

names <- read_csv(glue("{DATA_PATH}/study4_final.csv"), n_max = 0) %>% names()
raw <- read_csv(glue("{DATA_PATH}/study4_final.csv"), col_names = names, skip = 2)
df <- raw
df

# demographics ---------------------------------------------------------------

df <- df %>%
  filter(Finished == "1")

sex_map  <- c("1"="Male","2"="Female","3"="Other","4"="Prefer not to answer")

race_map <- c("1"="American Indian or Alaska Native", "2"="Asian",
  "3"="Black or African-American", "4"="Hispanic or Latino Origin",
  "5"="Hawaiian or Pacific Islander", "6"="White",
  "7"="Other", "8"="More than 1 of the above", "9"="Prefer not to answer"
)

df_completed <- df %>%
  filter(Finished %in% c(1, "1")) %>%
  mutate(
    age_num  = as.numeric(age),
    sex_chr  = dplyr::recode(as.character(sex),  !!!sex_map,  .default = NA_character_),
    race_chr = dplyr::recode(as.character(race), !!!race_map, .default = NA_character_)
  )

# 1) Session and person counts --------------------------
n_sessions <- nrow(df_completed)
n_people   <- n_distinct(df_completed$participant_id)

sessions_per_person <- df_completed %>%
  count(participant_id, name = "n_sessions")

n_repeat <- sum(sessions_per_person$n_sessions > 1)

# 2) person-level collapse if applicable ------------------
mode_chr <- function(x) {
  x <- x[!is.na(x)]
  if (!length(x)) return(NA_character_)
  tab <- sort(table(x), decreasing = TRUE)
  # deterministic tie-breaker: earliest value observed
  winners <- names(tab)[tab == max(tab)]
  # pick the first occurrence of any winner in the original (non-missing) order
  x_first <- x[match(winners, x)][1]
  x_first
}

first_nonmissing <- function(x) {
  y <- x[!is.na(x)]
  if (length(y)) y[1] else NA
}

# 3) One row per person --------------------------------
people <- df_completed %>%
  group_by(participant_id) %>%
  summarise(
    n_sessions     = n(),
    age            = first_nonmissing(age_num),
    sex            = mode_chr(sex_chr),
    race           = mode_chr(race_chr),
    sex_conflict   = n_distinct(sex_chr[!is.na(sex_chr)])  > 1,
    race_conflict  = n_distinct(race_chr[!is.na(race_chr)]) > 1,
    .groups = "drop"
  )

# 4) Person-level descriptive tables -------------------
age_desc <- people %>%
  summarise(
    n    = sum(!is.na(age)),
    mean = sprintf("%.2f", mean(age, na.rm = TRUE)),
    sd   = sd(age, na.rm = TRUE),
    med  = median(age, na.rm = TRUE),
    min  = min(age, na.rm = TRUE),
    max  = max(age, na.rm = TRUE)
  )

sex_tab <- people %>%
  count(sex, name = "count") %>%
  mutate(
    prop = count / sum(count),
    prop_formatted = sprintf("%.1f%%", prop * 100)  # "70.8%"
  )%>%
  arrange(desc(prop))

race_tab <- people %>%
  count(race, name = "count") %>%
  mutate(prop = count / sum(count),
  prop_formatted = sprintf("%.1f%%", prop * 100)
  ) %>%  
  arrange(desc(prop))

conflict_counts <- people %>%
  summarise(
    n_sex_conflict  = sum(sex_conflict,  na.rm = TRUE),
    n_race_conflict = sum(race_conflict, na.rm = TRUE)
  )

sessions_per_person %>% count(n_sessions, name = "n_people_with_this_many_sessions")

n_sessions 
n_people
age_desc
sex_tab 
race_tab

conflict_counts

# format data  -------------------------------------------------------------------------

# valence for all conditions

# For general and race conversation 1
df$valence_dissim_self_general_1 <- (df$thoughtvalence1_dissim_self_general_1 + df$thoughtvalence2_dissim_self_general_1) / 2
df$valence_dissim_other_general_1 <- (df$thoughtvalence1_dissim_other_general_1 + df$thoughtvalence2_dissim_other_general_1) / 2
df$valence_sim_self_general_1 <- (df$thoughtvalence1_sim_self_general_1 + df$thoughtvalence2_sim_self_general_1) / 2
df$valence_sim_other_general_1 <- (df$thoughtvalence1_sim_other_general_1 + df$thoughtvalence2_sim_other_general_1) / 2
df$valence_dissim_self_race_1 <- (df$thoughtvalence1_dissim_self_race_1 + df$thoughtvalence2_dissim_self_race_1) / 2
df$valence_dissim_other_race_1 <- (df$thoughtvalence1_dissim_other_race_1 + df$thoughtvalence2_dissim_other_race_1) / 2
df$valence_sim_self_race_1 <- (df$thoughtvalence1_sim_self_race_1 + df$thoughtvalence2_sim_self_race_1) / 2
df$valence_sim_other_race_1 <- (df$thoughtvalence1_sim_other_race_1 + df$thoughtvalence2_sim_other_race_1) / 2

# For general and race conversation 2
df$valence_dissim_self_general_2 <- (df$thoughtvalence1_dissim_self_general_2 + df$thoughtvalence2_dissim_self_general_2) / 2
df$valence_dissim_other_general_2 <- (df$thoughtvalence1_dissim_other_general_2 + df$thoughtvalence2_dissim_other_general_2) / 2
df$valence_sim_self_general_2 <- (df$thoughtvalence1_sim_self_general_2 + df$thoughtvalence2_sim_self_general_2) / 2
df$valence_sim_other_general_2 <- (df$thoughtvalence1_sim_other_general_2 + df$thoughtvalence2_sim_other_general_2) / 2
df$valence_dissim_self_race_2 <- (df$thoughtvalence1_dissim_self_race_2 + df$thoughtvalence2_dissim_self_race_2) / 2
df$valence_dissim_other_race_2 <- (df$thoughtvalence1_dissim_other_race_2 + df$thoughtvalence2_dissim_other_race_2) / 2
df$valence_sim_self_race_2 <- (df$thoughtvalence1_sim_self_race_2 + df$thoughtvalence2_sim_self_race_2) / 2
df$valence_sim_other_race_2 <- (df$thoughtvalence1_sim_other_race_2 + df$thoughtvalence2_sim_other_race_2) / 2

df <- df %>%
  mutate(
    thought_dissim_self_general_1 = paste(thought1_dissim_self_general_1, thought2_dissim_self_general_1, sep = " "),
    thought_dissim_other_general_1 = paste(thought1_dissim_other_general_1, thought2_dissim_other_general_1, sep = " "),
    thought_sim_self_general_1 = paste(thought1_sim_self_general_1, thought2_sim_self_general_1, sep = " "),
    thought_sim_other_general_1 = paste(thought1_sim_other_general_1, thought2_sim_other_general_1, sep = " "),
    thought_dissim_self_race_1 = paste(thought1_dissim_self_race_1, thought2_dissim_self_general_1, sep = " "),
    thought_dissim_other_race_1 = paste(thought1_dissim_other_race_1, thought2_dissim_other_general_1, sep = " "),
    thought_sim_self_race_1 = paste(thought1_sim_self_race_1, thought2_sim_self_general_1, sep = " "),
    thought_sim_other_race_1 = paste(thought1_sim_other_race_1, thought2_sim_other_general_1, sep = " "),
    
    thought_dissim_self_general_2 = paste(thought1_dissim_self_general_2, thought2_dissim_self_general_2, sep = " "),
    thought_dissim_other_general_2 = paste(thought1_dissim_other_general_2, thought2_dissim_other_general_2, sep = " "),
    thought_sim_self_general_2 = paste(thought1_sim_self_general_2, thought2_sim_self_general_2, sep = " "),
    thought_sim_other_general_2 = paste(thought1_sim_other_general_2, thought2_sim_other_general_2, sep = " "),
    thought_dissim_self_race_2 = paste(thought1_dissim_self_race_2, thought2_dissim_self_general_2, sep = " "),
    thought_dissim_other_race_2 = paste(thought1_dissim_other_race_2, thought2_dissim_other_general_2, sep = " "),
    thought_sim_self_race_2 = paste(thought1_sim_self_race_2, thought2_sim_self_general_2, sep = " "),
    thought_sim_other_race_2 = paste(thought1_sim_other_race_2, thought2_sim_other_general_2, sep = " ")
  )

# long format 
df_long <- df %>%
  pivot_longer(
    cols = starts_with("belief") | starts_with("conf") | starts_with("valence"),
    names_to = c(".value", "sim_dissim", "self_other", "general_race", "conversation"),
    names_sep = "_"
  )

df_long <- df_long %>%
  mutate(conversation = case_when(
    conversation == 1 ~ "one",
    conversation == 2 ~ "two"
  ))

head(df_long)

df_long

df_subset <- df_long[, c("participant_id", "sim_dissim", "self_other", "general_race", "conversation", "belief", "conf", "valence")]

#write.csv(df_subset, file = glue("{DATA_PATH}/z_study4_long_aggregate.csv"))

# primary effects  -------------------------------------------------------------------------

# Fit basic model
model <- lmer(belief ~ self_other * sim_dissim * general_race + (1|participant_id),
              data = df_subset, REML = FALSE)
summary(model)

# Total SD for standardization (between + within components)
vc <- as.data.frame(VarCorr(model))
sigma_total <- sqrt(sum(vc$vcov))

# emmeans contrast object into Cohen's d with CIs
contrast_to_d <- function(con_obj, sigma){
  s <- summary(con_obj, infer = c(TRUE, TRUE))  # get est, SE, df, CI
  s$d        <- s$estimate / sigma
  s$d.SE     <- s$SE / sigma
  s$d.lower  <- s$lower.CL / sigma
  s$d.upper  <- s$upper.CL / sigma
  s
}

# ----- H1: Main effect of partner type (similar vs dissimilar), collapsed over rating type
emm_H1  <- emmeans(model, ~ sim_dissim, pbkrtest.limit = 7856)                     # marginal means
con_H1  <- contrast(emm_H1, method = "pairwise", adjust="none")
d_H1   <- contrast_to_d(con_H1, sigma_total)
# report d_H1
emm_H1
con_H1
d_H1 

# ----- H2: Self vs other within dissimilar (and within similar, for completeness)
emm_H2  <- emmeans(model, ~ self_other | sim_dissim, pbkrtest.limit = 7856)        # simple effects by partner type
con_H2  <- contrast(emm_H2, method = "revpairwise")         # (other - self) within each sim_dissim
d_H2    <- contrast_to_d(con_H2, sigma_total)
# report row where sim_dissim == "dissim" as the preregistered H2
emm_H2
con_H2
d_H2

# ----- H3: Interaction (difference-in-differences)
# Take the self–other gap in dissimilar MINUS the self–other gap in similar
con_H3  <- contrast(con_H2, method = "pairwise", by = NULL) # dissimilar - similar
d_H3   <- contrast_to_d(con_H3, sigma_total)
# report d_H3
con_H3
d_H3

#note that this study now has a third factor, so our "overall 2x2 interaction with simple dummy coding, in the simple model output (summary(model)) 
#is not correct because that self_other:sim_dissim interaction is the value at the reference level of general_race
#And H3 is about the interaction in general, not just in one domain.
#so the contrast method below averages the two simple interactions (equal weights by default) 
#proof
#emm_by <- emmeans(model, ~ self_other | sim_dissim * general_race, pbkrtest.limit = 7856)
#gaps   <- contrast(emm_by, method = "revpairwise")  # (other - self), so mind the sign
#DID_by <- contrast(gaps, method = "pairwise", by = "general_race")
#summary(DID_by) 
#this shows the two interactions separately by domain

#more nuanced/correct model
# -----------------------------------------------
# Data preparation: effect-code the within-person factors
# -----------------------------------------------
# rating_type_c : +0.5 = self,    -0.5 = other
# partner_type_c: +0.5 = similar, -0.5 = dissimilar
# domain_type_c:  +0.5 = race,    -0.5 = general
# ("_c" denotes centered effect coding.)

df_subset <- df_subset |>
  dplyr::mutate(
    rating_type_c   = ifelse(self_other  == "self",  +0.5, -0.5),
    partner_type_c  = ifelse(sim_dissim  == "sim",   +0.5, -0.5),
    domain_type_c = ifelse(general_race == "race",  +0.5, -0.5)   # Domain
  )

# -----------------------------------------------
# Baseline mixed model: random intercepts only
# -----------------------------------------------
# Purpose: a compact specification that matches earlier reports
# and provides a reference when we compare with the preferred model.
fit_int <- lmer(
  belief ~ rating_type_c * partner_type_c * domain_type_c +
    (1 | participant_id),
  data  = df_subset,
  REML  = FALSE
)

summary(fit_int)

# -----------------------------------------------
# Preferred mixed model: random slopes (uncorrelated) by participant
# -----------------------------------------------
# We allow participants to vary in:
#   (i) the self-vs-other effect (rating_type_c),
#   (ii) the similar-vs-dissimilar effect (partner_type_c), and
#   (iii) their interaction (rating_type_c:partner_type_c).
# Using '||' (double-pipe) requests *uncorrelated* random effects,
# which reduces over-parameterization and helps avoid singular fits.
fit_rs <- lmer(
  belief ~ rating_type_c * partner_type_c * domain_type_c +
    (1 + rating_type_c + partner_type_c + rating_type_c:partner_type_c || participant_id),
  data    = df_subset,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)),
  REML    = FALSE
)
summary(fit_rs)
confint(fit_rs)

# Optionally compare (interpret LRTs cautiously when the richer model is near-boundary)
anova(fit_int, fit_rs)

# Use the preferred model for effect estimation
model <- fit_rs

vc <- as.data.frame(VarCorr(model))

get_v <- function(name) {
  v <- vc$vcov[vc$var1 == name & vc$grp != "Residual"]
  if (length(v) == 0) 0 else sum(v)
}

v_int <- get_v("(Intercept)")
v_S   <- get_v("rating_type_c")
v_D   <- get_v("partner_type_c")
v_SD  <- get_v("rating_type_c:partner_type_c")
v_res <- vc$vcov[vc$grp == "Residual"][1]

sigma_total <- sqrt(v_int + 0.25 * v_S + 0.25 * v_D + 0.0625 * v_SD + v_res)

# -----------------------------------------------
# Helper: convert an emmeans contrast to Cohen's d with CIs
# -----------------------------------------------
contrast_to_d <- function(con_obj, sigma) {
  s <- summary(con_obj, infer = c(TRUE, TRUE))  # estimate, SE, df, CI
  s$d        <- s$estimate / sigma
  s$d.SE     <- s$SE / sigma
  s$d.lower  <- s$lower.CL / sigma
  s$d.upper  <- s$upper.CL / sigma
  s
}

# -----------------------------------------------
# H1: Homophilous expectations (partner main effect)
# -----------------------------------------------
emm_H1 <- emmeans(
  model, ~ partner_type_c,
  pbkrtest.limit = 7856,
  at = list(partner_type_c = c(-0.5, +0.5))
)

# Define contrast explicitly as sim - dissim (positive = homophily)
con_H1 <- contrast(emm_H1, method = list("sim - dissim" = c(-1, 1)), adjust = "none")
d_H1   <- contrast_to_d(con_H1, sigma_total)

emm_H1
con_H1
d_H1

# -----------------------------------------------
# H2: Self vs. other *within dissimilar* (and within similar, for completeness)
# -----------------------------------------------
emm_H2 <- emmeans(
  model, ~ rating_type_c | partner_type_c,
  pbkrtest.limit = 7856,
  at = list(rating_type_c = c(-0.5, +0.5), partner_type_c = c(-0.5, +0.5))
)

# Contrast as self - other within each partner_type_c
con_H2 <- contrast(emm_H2, method = list("self - other" = c(-1, 1)), by = "partner_type_c")
d_H2   <- contrast_to_d(con_H2, sigma_total)

emm_H2
con_H2
d_H2   # Use the row with partner_type_c = -0.5 ("dissimilar") as the preregistered H2

# -----------------------------------------------
# H3: Interaction (difference-in-differences)
# -----------------------------------------------
# (self - other)_dissimilar  -  (self - other)_similar
# should be same as raw 2x2 interaction because now effects coded so raw not at reference level but grand mean
con_H3 <- contrast(con_H2, method = "pairwise", by = NULL)  # dissimilar - similar
d_H3   <- contrast_to_d(con_H3, sigma_total)

con_H3
d_H3

# confidence, valence, and projection  -------------------------------------------------------------------------

#confidence 
conf_model <- lmer(conf ~ self_other * sim_dissim * general_race + (1 | participant_id), data = df_subset)
summary(conf_model)
emmeans_conf_1 <- emmeans(conf_model, specs = ~ sim_dissim, pbkrtest.limit = 7856)
emmeans_conf_1
emmeans_conf_2 <- emmeans(conf_model, specs = ~ self_other * sim_dissim, pbkrtest.limit = 7856)
emmeans_conf_2
emmeans_conf_3 <- emmeans(conf_model, specs = ~ self_other * sim_dissim * general_race, pbkrtest.limit = 7856)
emmeans_conf_3

#valence
valence_model <- lmer(valence ~ self_other * sim_dissim + (1 | participant_id), data = df_subset)
summary(valence_model)
emmeans_valence_1 <- emmeans(valence_model, specs = ~ sim_dissim, pbkrtest.limit = 7856)
emmeans_valence_1
emmeans_valence_2 <- emmeans(valence_model, specs = ~ self_other * sim_dissim, pbkrtest.limit = 7856)
emmeans_valence_2

# mediation  -------------------------------------------------------------------------

set.seed(12345)

run_moderated_parallel_mediation <- function(data_subset, R = 1000, ci_type = "perc") {
  
  # Centered numeric contrasts (clear and stable)
  dat <- data_subset %>%
    mutate(
      rating_type_c  = ifelse(self_other == "self",  +0.5, -0.5),  # X
      partner_type_c = ifelse(sim_dissim == "sim",   +0.5, -0.5),  # Moderator
      M1 = conf,
      M2 = valence,
      Y  = belief
    )
  
  # Full random-slopes models (NO fallback; warnings will print if they occur)
  fm_m1 <- M1 ~ rating_type_c * partner_type_c +
    (1 + rating_type_c + partner_type_c + rating_type_c:partner_type_c || participant_id)
  
  fm_m2 <- M2 ~ rating_type_c * partner_type_c +
    (1 + rating_type_c + partner_type_c + rating_type_c:partner_type_c || participant_id)
  
  fm_y  <- Y  ~ rating_type_c * partner_type_c +
    M1 * partner_type_c + M2 * partner_type_c +
    (1 + rating_type_c + partner_type_c + rating_type_c:partner_type_c || participant_id)
  
  # Helper: compute the 15 path quantities from a data.frame
  compute_effects <- function(df) {
    # Fit full models; if any fitting error, propagate NA vector
    m1 <- try(lmer(fm_m1, data = df, REML = FALSE,
                   control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))),
              silent = FALSE)
    if (inherits(m1, "try-error")) return(rep(NA_real_, 15))
    
    m2 <- try(lmer(fm_m2, data = df, REML = FALSE,
                   control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))),
              silent = FALSE)
    if (inherits(m2, "try-error")) return(rep(NA_real_, 15))
    
    m3 <- try(lmer(fm_y,  data = df, REML = FALSE,
                   control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))),
              silent = FALSE)
    if (inherits(m3, "try-error")) return(rep(NA_real_, 15))
    
    fe1 <- fixef(m1); fe2 <- fixef(m2); fe3 <- fixef(m3)
    gb  <- function(fe, nm) if (nm %in% names(fe)) unname(fe[[nm]]) else 0
    
    # Moderator at dissimilar (-.5) vs similar (+.5)
    Wd <- -0.5; Ws <- +0.5
    
    # a-paths: X -> M1/M2
    a1_d <- gb(fe1,"rating_type_c") + Wd * gb(fe1,"rating_type_c:partner_type_c")
    a1_s <- gb(fe1,"rating_type_c") + Ws * gb(fe1,"rating_type_c:partner_type_c")
    a2_d <- gb(fe2,"rating_type_c") + Wd * gb(fe2,"rating_type_c:partner_type_c")
    a2_s <- gb(fe2,"rating_type_c") + Ws * gb(fe2,"rating_type_c:partner_type_c")
    
    # b-paths: M1/M2 -> Y
    b1_d <- gb(fe3,"M1") + Wd * gb(fe3,"M1:partner_type_c")
    b1_s <- gb(fe3,"M1") + Ws * gb(fe3,"M1:partner_type_c")
    b2_d <- gb(fe3,"M2") + Wd * gb(fe3,"M2:partner_type_c")
    b2_s <- gb(fe3,"M2") + Ws * gb(fe3,"M2:partner_type_c")
    
    # c' (direct): X -> Y
    cp_d <- gb(fe3,"rating_type_c") + Wd * gb(fe3,"rating_type_c:partner_type_c")
    cp_s <- gb(fe3,"rating_type_c") + Ws * gb(fe3,"rating_type_c:partner_type_c")
    
    # Indirects and totals
    ind1_d <- a1_d * b1_d; ind2_d <- a2_d * b2_d; tind_d <- ind1_d + ind2_d
    ind1_s <- a1_s * b1_s; ind2_s <- a2_s * b2_s; tind_s <- ind1_s + ind2_s
    
    c_d <- cp_d + tind_d
    c_s <- cp_s + tind_s
    
    # Differences (Sim − Dissim)
    d_tot  <- c_s  - c_d
    d_dir  <- cp_s - cp_d
    d_ind1 <- ind1_s - ind1_d
    d_ind2 <- ind2_s - ind2_d
    d_tind <- tind_s - tind_d
    
    c(c_d, cp_d, tind_d, ind1_d, ind2_d,
      c_s, cp_s, tind_s, ind1_s, ind2_s,
      d_tot, d_dir, d_ind1, d_ind2, d_tind)
  }
  
  # Point estimates on the observed data
  t0 <- compute_effects(dat)
  
  # Participant-level bootstrap (statistic returns 15-vector; failures return NAs)
  ids_all <- unique(dat$participant_id)
  stat_fun <- function(id_vec, i, original_data) {
    samp_ids <- id_vec[i]
    d2 <- do.call(rbind, lapply(samp_ids, function(pid)
      original_data[original_data$participant_id == pid, , drop = FALSE]))
    compute_effects(d2)
  }
  
  boot_obj <- boot(data = ids_all,
                   statistic = stat_fun,
                   R = R,
                   original_data = dat)
  
  # Replace t0 with our explicit computation
  boot_obj$t0 <- t0
  
  # Clean out any failed replicates (rows with any NA) and adjust R accordingly
  keep <- stats::complete.cases(boot_obj$t)
  R_eff <- sum(keep)
  if (R_eff < R) {
    message(sprintf("Removed %d/%d bootstrap replicates due to model failures; using R_eff = %d.",
                    R - R_eff, R, R_eff))
  }
  if (R_eff < 50L) {
    warning("Very few valid bootstrap replicates retained; CIs/p-values may be unstable.")
  }
  
  boot_obj$t <- boot_obj$t[keep, , drop = FALSE]
  boot_obj$R <- nrow(boot_obj$t)
  
  # Store the CI type used as an attribute (for downstream functions)
  attr(boot_obj, "ci_type") <- ci_type
  
  boot_obj
}

#-----------------------------------------------
# 2) Tidy table: percentile CIs + boot.pval p's
#-----------------------------------------------
create_tidy_boot_table_moderated <- function(boot_obj, ci_type = NULL) {
  
  if (is.null(ci_type)) ci_type <- attr(boot_obj, "ci_type") %||% "perc"
  
  # Map boot.ci 'type' to the element name in the return object
  ci_element <- switch(ci_type,
                       "perc" = "percent",
                       "basic" = "basic",
                       "norm" = "normal",
                       "bca"  = "bca",
                       "stud" = "student",
                       stop("Unknown ci_type: ", ci_type))
  
  terms <- c(
    "Total effect (c) – Dissimilar",
    "Direct effect (c′) – Dissimilar",
    "Total indirect – Dissimilar",
    "Indirect via Confidence – Dissimilar",
    "Indirect via Valence – Dissimilar",
    "Total effect (c) – Similar",
    "Direct effect (c′) – Similar",
    "Total indirect – Similar",
    "Indirect via Confidence – Similar",
    "Indirect via Valence – Similar",
    "Δ Total effect (Similar − Dissimilar)",
    "Δ Direct effect (Similar − Dissimilar)",
    "Δ Indirect via Confidence (Similar − Dissimilar)",
    "Δ Indirect via Valence (Similar − Dissimilar)",
    "Δ Total indirect (Similar − Dissimilar)"
  )
  
  est <- as.numeric(boot_obj$t0)
  J   <- length(est)
  
  ci_lo <- ci_hi <- rep(NA_real_, J)
  p_raw <- rep(NA_real_, J)
  
  for (j in seq_len(J)) {
    # Percentile CI (or other 'type'); extract the correct element name
    ci <- try(boot.ci(boot_obj, type = ci_type, index = j), silent = FALSE)
    if (!inherits(ci, "try-error")) {
      elem <- ci[[ci_element]]
      if (!is.null(elem) && length(elem) >= 5) {
        # boot.ci columns: conf, conf, point, lower, upper
        ci_lo[j] <- elem[4]
        ci_hi[j] <- elem[5]
      }
    }
    
    # CI-inversion p-value via boot.pval (no fallback)
    p_raw[j] <- boot.pval::boot.pval(boot_obj, type = ci_type, theta_null = 0, index = j)
  }
  
  p_fmt <- ifelse(is.na(p_raw), NA_character_,
                  ifelse(p_raw < 1e-4, sprintf("%.4e", p_raw), sprintf("%.4f", p_raw)))
  sig <- dplyr::case_when(
    is.na(p_raw) ~ NA_character_,
    p_raw < .001 ~ "< .001",
    p_raw < .01  ~ "< .01",
    p_raw < .05  ~ "< .05",
    TRUE ~ "ns"
  )
  
  tibble::tibble(
    Term = terms,
    Estimate = est,
    CI_lower = ci_lo,
    CI_upper = ci_hi,
    p_value_raw = p_fmt,
    p_value_formatted = sig
  )
}

# Overall
boot_overall <- run_moderated_parallel_mediation(df_subset, R = 200)  # e.g., 500 for testing
tbl_overall  <- create_tidy_boot_table_moderated(boot_overall)
print(tbl_overall)

# Split by general_race (same API)
boot_general <- run_moderated_parallel_mediation(dplyr::filter(df_subset, general_race == "general"), R = 200)
tbl_general  <- create_tidy_boot_table_moderated(boot_general)
print(tbl_general)

boot_race <- run_moderated_parallel_mediation(dplyr::filter(df_subset, general_race == "race"), R = 200)
tbl_race  <- create_tidy_boot_table_moderated(boot_race)
print(tbl_race)

#hand check a few numbers for sanity
idx     <- 10
theta0  <- boot_overall$t0[idx]           # your point estimate (~ 0.0197)
thetas  <- boot_overall$t[, idx]          # bootstrap replicates for this effect
R_eff   <- nrow(boot_overall$t)

# 95% percentile CI (matches your table):
ci_perc <- quantile(thetas, c(.025, .975), na.rm = TRUE)
ci_perc

k      <- sum(thetas <= 0, na.rm = TRUE)    # how many bootstrap effects ≤ 0
p_est  <- 2 * min(k, R_eff - k) / R_eff     # two-sided p (resolution 1/R_eff)
p_est
# With no draws ≤ 0, p_est = 0 -> reported as 1/R_eff = 0.001 via boot.pval

#confirm with BCa
boot.ci(boot_overall, type = c("perc", "bca"), index = idx)

#a and b paths
fit_mediation_models <- function(data_subset) {
  dat <- data_subset %>%
    mutate(
      rating_type_c  = ifelse(self_other == "self",  +0.5, -0.5),  # X
      partner_type_c = ifelse(sim_dissim == "sim",   +0.5, -0.5),  # W
      M1 = conf,                                                  # mediator 1
      M2 = valence,                                               # mediator 2
      Y  = belief                                                # outcome
    )
  
  fm_m1 <- M1 ~ rating_type_c * partner_type_c +
    (1 + rating_type_c + partner_type_c + rating_type_c:partner_type_c || participant_id)
  
  fm_m2 <- M2 ~ rating_type_c * partner_type_c +
    (1 + rating_type_c + partner_type_c + rating_type_c:partner_type_c || participant_id)
  
  fm_y  <- Y  ~ rating_type_c * partner_type_c +
    M1 * partner_type_c + M2 * partner_type_c +
    (1 + rating_type_c + partner_type_c + rating_type_c:partner_type_c || participant_id)
  
  m1 <- lmer(fm_m1, data = dat, REML = FALSE,
             control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
  m2 <- lmer(fm_m2, data = dat, REML = FALSE,
             control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
  m3 <- lmer(fm_y,  data = dat, REML = FALSE,
             control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
  
  list(dat = dat, m1 = m1, m2 = m2, m3 = m3)
}

## 2) wrapper around emtrends() that standardizes column names ----------
.emtr_to_tbl <- function(emm_grid, path_label, w_levels = c(-0.5, +0.5)) {
  
  s <- summary(
    emm_grid,
    infer = c(TRUE, TRUE)         # <-- ensures t & p are included
  ) %>% as.data.frame()
  
  # The estimate column is named "<var>.trend" (e.g., "rating_type_c.trend", "M1.trend").
  trend_col <- grep("\\.trend$", names(s), value = TRUE)
  if (length(trend_col) != 1L) stop("Couldn't locate the '.trend' column returned by emtrends().")
  
  # Standardize column names and compute t/p if not present.
  s <- s %>%
    mutate(
      estimate  = .data[[trend_col]],
      Condition = ifelse(partner_type_c == w_levels[1], "Dissimilar", "Similar"),
      t         = if ("t.ratio" %in% names(.)) `t.ratio` else estimate / SE,
      p.value   = if ("p.value" %in% names(.)) `p.value`
      else 2 * pt(abs(t), df = df, lower.tail = FALSE),
      Path      = path_label
    ) %>%
    select(Path, Condition, estimate, SE, df, t, p.value, lower.CL, upper.CL)
  
  s
}

## 3) Extract a1,a2,b1,b2 at W = −.5 / +.5 with tests & CIs ------------------
extract_paths_table <- function(mods) {
  atW <- list(partner_type_c = c(-0.5, +0.5))  # −.5 = Dissimilar, +.5 = Similar
  
  a1 <- emtrends(mods$m1, ~ partner_type_c, var = "rating_type_c", at = atW,
                 pbkrtest.limit = nrow(mods$dat)) %>%
    .emtr_to_tbl("a1: X \u2192 M1")
  
  a2 <- emtrends(mods$m2, ~ partner_type_c, var = "rating_type_c", at = atW,
                 pbkrtest.limit = nrow(mods$dat)) %>%
    .emtr_to_tbl("a2: X \u2192 M2")
  
  b1 <- emtrends(mods$m3, ~ partner_type_c, var = "M1", at = atW,
                 pbkrtest.limit = nrow(mods$dat)) %>%
    .emtr_to_tbl("b1: M1 \u2192 Y | X,M2")
  
  b2 <- emtrends(mods$m3, ~ partner_type_c, var = "M2", at = atW,
                 pbkrtest.limit = nrow(mods$dat)) %>%
    .emtr_to_tbl("b2: M2 \u2192 Y | X,M1")
  
  bind_rows(a1, a2, b1, b2) %>%
    arrange(Path, factor(Condition, levels = c("Dissimilar","Similar")))
}

##  Run --------------------------------------------------------
mods  <- fit_mediation_models(df_subset)   # will warn 'boundary (singular)' sometimes; OK for fixed effects
paths <- extract_paths_table(mods)
paths


## (Optional) Compute raw a*b products by condition ------------------------
compute_raw_indirects <- function(paths_tbl) {
  # reshape to a wide-by-path matrix of estimates
  wide_est <- paths_tbl %>%
    select(Path, Condition, estimate) %>%
    mutate(Path = case_when(
      str_detect(Path, "^a1") ~ "a1",
      str_detect(Path, "^a2") ~ "a2",
      str_detect(Path, "^b1") ~ "b1",
      str_detect(Path, "^b2") ~ "b2",
      TRUE ~ Path
    )) %>%
    pivot_wider(names_from = Path, values_from = estimate)
  
  wide_est %>%
    mutate(
      Indirect_M1 = a1 * b1,
      Indirect_M2 = a2 * b2
    ) %>%
    select(Condition, a1, b1, Indirect_M1, a2, b2, Indirect_M2)
}

raw_indirects <- compute_raw_indirects(paths)
raw_indirects

#don't match exacty - not a problem
#For inference on indirects, the bootstrap is the right tool because an indirect is a non‑linear function of multiple coefficients with non‑trivial covariance. 
#Bootstrap approach is doing exactly what it should. (The p‑values reported are obtained by CI‑inversion of the bootstrap distribution, which appears to be the recommended practice for such composites.

##################################

### AGGREGATE TABLE

##################################

d_long_1 <- read_csv(glue("{DATA_PATH}/study1_long_aggregate.csv"))
d_long_2 <- read_csv(glue("{DATA_PATH}/study2_long_aggregate.csv"))
d_long_3a <- read_csv(glue("{DATA_PATH}/study3a_long_aggregate.csv"))
d_long_3b <- read_csv(glue("{DATA_PATH}/study3b_long_aggregate.csv"))
d_long_4 <- read_csv(glue("{DATA_PATH}/study4_long_aggregate.csv"))
d_long_5 <- read_csv(glue("{DATA_PATH}/study5_long_aggregate.csv"))

#recodes
d_long_4 <- d_long_4 %>%
  mutate(
    belief = as.numeric(belief)
  ) %>%
  group_by(pid, conversation, general_race, sim_dissim, self_other) %>%
  mutate(
    rep = row_number()    # 1, 2, 3, ... for each repeat within pid × conversation × condition
  ) %>%
  ungroup() %>%
  mutate(
    pid_conv_rep = paste(pid, conversation, rep, sep = "_")  # unique per observation
  )

d_long_4 %>%
  filter(general_race == "race",
         sim_dissim == "dissim") %>%
  summarise(n = n(), .by = c(pid_conv_rep, self_other)) %>%
  filter(n > 1L)


classify_bias <- function(data,
                          dv,
                          study,
                          domain,
                          sample_type,
                          aggregate = FALSE,
                          sim_var      = "sim_dissim",
                          self_var     = "self_other",
                          dissim_code  = "dissim",   # adjust if needed
                          self_code    = "self",
                          other_code   = "other",
                          id_var       = "pid") {
  
  dv_sym   <- rlang::sym(dv)
  sim_sym  <- rlang::sym(sim_var)
  self_sym <- rlang::sym(self_var)
  id_sym   <- rlang::sym(id_var)
  
  # isolate dissimilar condition and the DV
  wide_df <- data %>%
    filter(!!sim_sym == dissim_code) %>%
    select(!!id_sym, !!self_sym, !!dv_sym) %>%
    pivot_wider(
      names_from  = !!self_sym,
      values_from = !!dv_sym
    ) %>%
    rename(
      self  = !!self_code,
      other = !!other_code
    )
  
  # bias score: Self_dissim - Other_dissim
  wide_df <- wide_df %>%
    mutate(
      bias  = self - other,
      thr   = if (aggregate) 1 else 0,  # threshold for "no bias"
      bias_cat = case_when(
        bias >  thr ~ "Negative bias",  # dissimilarity pessimism
        bias < -thr ~ "Positive bias",  # reverse
        TRUE        ~ "No bias"
      )
    )
  
  # percentages
  bias_summary <- wide_df %>%
    count(bias_cat) %>%
    mutate(
      prop = n / sum(n),
      study  = study,
      domain = domain,
      sample = sample_type
    ) %>%
    select(study, domain, sample, bias_cat, n, prop)
  
  bias_summary
}

bias_all <- bind_rows(

classify_bias(d_long_1,
              dv          = "rating",
              study       = "Study 1",
              domain      = "General",
              sample_type = "Lab",
              aggregate   = FALSE),

classify_bias(d_long_2,
              dv          = "rating",
              study       = "Study 2",
              domain      = "Race",
              sample_type = "Online",
              aggregate   = FALSE),

classify_bias(d_long_3a,
              dv          = "rating",
              study       = "Study 3a",
              domain      = "Personality",
              sample_type = "Online",
              aggregate   = FALSE),

classify_bias(d_long_3b,
              dv          = "rating",
              study       = "Study 3b",
              domain      = "Personality",
              sample_type = "Online",
              aggregate   = FALSE),


classify_bias(d_long_4 %>% filter(general_race == "race"),     # adjust 'domain_var' as needed
              dv          = "belief",
              study       = "Study 4",
              domain      = "Race",
              sample_type = "Lab",
              aggregate   = FALSE,
              id_var      = "pid_conv_rep"),   

classify_bias(d_long_4 %>% filter(general_race  == "general"),
              dv          = "belief",
              study       = "Study 4",
              domain      = "Personality",
              sample_type = "Lab",
              aggregate   = FALSE,
              id_var      = "pid_conv_rep"),

# Study 5: field; general and five organizational dimensions
classify_bias(d_long_5,
              dv          = "overall",        # "general"
              study       = "Study 5",
              domain      = "General",
              sample_type = "Field",
              aggregate   = FALSE)
)

classify_bias(d_long_5,
              dv          = "careerstage",
              study       = "Study 5",
              domain      = "Career stage",
              sample_type = "Field",
              aggregate   = FALSE)

classify_bias(d_long_5,
              dv          = "position",
              study       = "Study 5",
              domain      = "Job level",
              sample_type = "Field",
              aggregate   = FALSE)

classify_bias(d_long_5,
              dv          = "division",
              study       = "Study 5",
              domain      = "Division",
              sample_type = "Field",
              aggregate   = FALSE)

classify_bias(d_long_5,
              dv          = "industry",
              study       = "Study 5",
              domain      = "Industry",
              sample_type = "Field",
              aggregate   = FALSE)

classify_bias(d_long_5,
              dv          = "sociocultural",
              study       = "Study 5",
              domain      = "Sociocultural",
              sample_type = "Field",
              aggregate   = FALSE)

bias_plot_df <- bias_all %>%
  group_by(study, domain, sample) %>%
  mutate(total_n = sum(n)) %>%
  ungroup() %>%
  mutate(
    # Halve N for Study 4 only (both domains)
    total_n = ifelse(study == "Study 4", total_n / 2, total_n),
    bias_cat = factor(
      bias_cat,
      levels = c("Negative bias", "No bias", "Positive bias")
    ),
    StudyDomain = glue("{study}: {domain}\n({sample}, N = {total_n})"),
    StudyDomain = factor(StudyDomain, levels = unique(StudyDomain))
  ) %>%
  arrange(StudyDomain, bias_cat) # ensures within-study order: Neg, No, Pos

str(bias_plot_df)

fig5 <- ggplot(bias_plot_df,
       aes(x = StudyDomain, y = prop, fill = bias_cat)) +
  geom_col(
    width = 0.7,
    color = "black",
    linewidth = 0.5,
    position = position_stack(reverse = TRUE)  # <-- KEY FIX
  ) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_fill_manual(
    values = c(
      "Negative bias" = "grey15",   # darkest, should be bottom
      "No bias"       = "grey60",
      "Positive bias" = "grey90"    # lightest, should be top
    ),
    breaks = c("Negative bias", "No bias", "Positive bias"),
    labels = c("Negative bias (Self > Other)",
               "No bias",
               "Positive bias (Self < Other)")
  ) +
  labs(
    x    = NULL,
    y    = "Percentage of participants",
    fill = "Bias type"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x    = element_text(angle = 30, hjust = 1, vjust = 1),
    legend.position = "bottom",
    legend.title   = element_text(size = 12),
    legend.text    = element_text(size = 11),
    axis.line      = element_line(linewidth = 0.6),
    panel.spacing.x = unit(1, "lines")
  )


fig5

save_plot_multi <- function(plot,
                            file_stem,
                            exts   = c("png", "jpeg", "pdf"),
                            width  = 8,
                            height = 6,
                            dpi    = 300,
                            path   = PLOT_PATH) {
  # If the object is from grid.arrange() we convert it first
  if (inherits(plot, "gtable")) {
    plot <- gridExtra::arrangeGrob(plot)
  }
  
  # Iterate over requested extensions
  purrr::walk(exts, function(ext) {
    ggsave(
      filename = glue::glue("{file_stem}.{ext}"),
      plot     = plot,
      device   = ext,          # lets ggsave pick the right device
      path     = path,
      width    = width,
      height   = height,
      dpi      = dpi
    )
  })
  invisible(TRUE)
}

save_plot_multi(fig5,
                file_stem = "fig5",
                exts      = c("png", "pdf", "jpeg"),   # add/remove as needed
                width     = 10,
                height    = 8)


