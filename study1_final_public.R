# Title: Dissimilarity Pessimism Study 1
# Author(s): Gus Cooney; Erica Boothby
# Description: Initial demonstration of the effect

# load libraries ----------------------------------------------------------

library(pacman)
p_load(tidyverse, glue, broom, readr, lme4, lmerTest, emmeans, effsize, rempsyc, magrittr, psych)

# set constants -----------------------------------------------------------

DATA_PATH <- "~/Dropbox/Manuscripts/Homophily/Data"
PLOT_PATH <- "~/Dropbox/Manuscripts/Homophily/Plots"

# load data ---------------------------------------------------------------

d <- read_csv(glue("{DATA_PATH}/study1_final.csv"))
head(d)

# demographics ---------------------------------------------------------------

df <- d

sex_map  <- c("1"="Male","2"="Female","3"="Other")
race_map <- c(
  "1"="American Indian or Alaska Native", "2"="Asian",
  "3"="Black or African-American", "4"="Hispanic or Latino Origin",
  "5"="Hawaiian or Pacific Islander", "6"="White",
  "7"="Other", "8"="More than 1 of the above", "9"="Prefer not to answer"
)

df_completed <- df %>%
  mutate(
    age_num  = suppressWarnings(as.numeric(age)),
    sex_chr  = dplyr::recode(as.character(sex),  !!!sex_map,  .default = NA_character_),
    race_chr = dplyr::recode(as.character(race), !!!race_map, .default = NA_character_)
  )


# 1) Session and person counts --------------------------
n_sessions <- nrow(df_completed)
n_people   <- n_distinct(df_completed$pid)

sessions_per_person <- df_completed %>%
  dplyr::count(pid, name = "n_sessions")

n_repeat <- sum(sessions_per_person$n_sessions > 1)

# 2) helpers for person-level collapse ------------------
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
  group_by(pid) %>%
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
  dplyr::count(sex, name = "count") %>%
  mutate(
    prop = count / sum(count),
    prop_formatted = sprintf("%.1f%%", prop * 100)  # "70.8%"
  )%>%
  arrange(desc(prop))

race_tab <- people %>%
  dplyr::count(race, name = "count") %>%
  mutate(prop = count / sum(count),
         prop_formatted = sprintf("%.1f%%", prop * 100)
  ) %>%  
  arrange(desc(prop))

conflict_counts <- people %>%
  summarise(
    n_sex_conflict  = sum(sex_conflict,  na.rm = TRUE),
    n_race_conflict = sum(race_conflict, na.rm = TRUE)
  )

sessions_per_person %>% dplyr::count(n_sessions, name = "n_people_with_this_many_sessions")

n_sessions 
n_people
age_desc
sex_tab 
race_tab

conflict_counts

# format data ---------------------------------------------------------------

d_long <- d %>%
  pivot_longer(
    cols = c(sim_self, sim_other, dissim_self, dissim_other), 
    names_to = c("sim_dissim", "self_other"),
    names_sep = "_",
    values_to = "rating"
  )

d_long

d_long <- d_long |>
  dplyr::mutate(
    rating_type_c   = ifelse(self_other  == "self",  +0.5, -0.5),
    partner_type_c  = ifelse(sim_dissim  == "sim",   +0.5, -0.5),
  )

write.csv(d_long, file = glue("{DATA_PATH}/z_study1_long_aggregate.csv"))

# Modeling ---------------------------------------------------------------

# -----------------------------------------------
# Preferred mixed model: random slopes (uncorrelated) by participant
# -----------------------------------------------
# We allow participants to vary in:
#   (i) the self-vs-other effect (rating_type_c),
#   (ii) the similar-vs-dissimilar effect (partner_type_c), and
#   (iii) their interaction (rating_type_c:partner_type_c) - if can converge / relevant?
# Using '||' (double-pipe) requests *uncorrelated* random effects - reduces over-parameterization 

fit_rs_full <- lmer(
  rating ~ rating_type_c * partner_type_c +
    (1 + rating_type_c + partner_type_c + rating_type_c:partner_type_c || pid),
  data    = d_long,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)),
  REML    = FALSE
)
summary(fit_rs_full)
confint(fit_rs_full)

fit_rs_simple <- lmer(
  rating ~ rating_type_c * partner_type_c +
    (1|self_other:pid) + (1|sim_dissim:pid),
  data    = d_long,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)),
  REML    = FALSE
)
summary(fit_rs_simple)
confint(fit_rs_simple)

fit_rs <- lmer(
  rating ~ rating_type_c * partner_type_c +
    (1 + rating_type_c + partner_type_c || pid),
  data    = d_long,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)),
  REML    = FALSE
)
summary(fit_rs)
confint(fit_rs)

#Notes:

# Start Here: Effects coded, uncorrelated random effects
#(1 + rating_type_c + partner_type_c || pid)
#3 variance parameters
#Captures individual differences in main effects
#Stable convergence

# If You Care About Correlations - Try the single pipe version
#(1 + rating_type_c + partner_type_c | pid)
#Adds 3 correlation parameters between random effects
#If singular fit → revert to double pipe ||
  
#If You Think Interaction Varies:
#(1 + rating_type_c + partner_type_c + rating_type_c:partner_type_c || pid)
#Check if σ²_interaction > 0
#If ≈0, drop it

#Your Original (Also Fine!):
# Implicit random slopes through crossed grouping
#(1|self_other:pid) + (1|sim_dissim:pid)
#Captures main effect variation
#Doesn't include baseline person variance or interaction variance
#highly stable

#equivalent ANOVA
#anova.model <- aov(rating ~ self_other*sim_dissim + Error(pid/(self_other*sim_dissim)), data=d_long)
#summary(anova.model)
#Means.anova <- emmeans(anova.model, specs = ~ self_other*sim_dissim)
#Means.anova
#contr.anova <- contrast(Means.anova, method = "pairwise", adjust = "holm")
#contr.anova
#confint(contr.anova)

# Use the preferred model for effect estimation
model <- fit_rs

# -----------------------------------------------
# Total SD for standardization (Cohen's d)
# -----------------------------------------------
# With ±0.5 coding and *uncorrelated* random effects:
# Var_total = Var(b0) + .25 Var(b_rating) + .25 Var(b_partner) + .0625 Var(b_interaction) + Var(residual)
# If cell frequencies are not balanced, or you use different coding, the correct multipliers are the empirical second moments of the relevant columns:
# If use correlated random effects (single |), we'd also have covariance terms
# if interaction term is missing, it's just zero and code still works below

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

pretty_emm_H1 <- summary(emm_H1) %>%
  mutate(partner = dplyr::recode(as.character(partner_type_c),
                                 "-0.5" = "dissim",
                                 "0.5"  = "sim"),
  `95% CI` = sprintf("%.2f [%.2f, %.2f]", emmean, lower.CL, upper.CL)) %>%
  dplyr::select(partner, everything(), -partner_type_c)


# Define contrast explicitly as sim - dissim (positive = homophily)
con_H1 <- contrast(emm_H1, method = list("sim - dissim" = c(-1, 1)), adjust = "none")
con_d_H1   <- contrast_to_d(con_H1, sigma_total)
#note: df is slightly different than full model because emmeans re-computes the Satterthwaite/KR approximation, while summary(fit_rs) reports it for the individual coefficient.

emm_H1 
pretty_emm_H1
con_d_H1

# -----------------------------------------------
# H2: Self vs. other *within dissimilar* (and within similar, for completeness)
# -----------------------------------------------
emm_H2 <- emmeans(
  model, ~ rating_type_c | partner_type_c,
  pbkrtest.limit = 7856,
  at = list(rating_type_c = c(-0.5, +0.5), 
            partner_type_c = c(-0.5, +0.5))
)

pretty_emm_H2 <- summary(emm_H2) %>%
  dplyr::mutate(
    rating  = dplyr::recode(as.character(rating_type_c),
                            `-0.5` = "other", `0.5` = "self"),
    partner = dplyr::recode(as.character(partner_type_c),
                            `-0.5` = "dissim", `0.5` = "similar"),
    `95% CI` = sprintf("%.2f [%.2f, %.2f]", emmean, lower.CL, upper.CL)
  ) %>%
  dplyr::select(partner, rating, emmean, SE, df, `95% CI`) %>%
  dplyr::arrange(partner, rating)

# Contrast as self - other within each partner_type_c
con_H2 <- contrast(emm_H2, method = list("self - other" = c(-1, 1)), by = "partner_type_c")
con_d_H2   <- contrast_to_d(con_H2, sigma_total)

pretty_con_d_H2 <- con_d_H2 %>%
  dplyr::mutate(
    partner = dplyr::recode(as.character(partner_type_c),
                            `-0.5` = "dissim", `0.5` = "similar"),
    contrast = "self – other"
  ) %>%
  dplyr::transmute(
    partner, contrast,
    estimate, SE, df, t.ratio,
    `95% CI (raw)` = sprintf("%.2f [%.2f, %.2f]", estimate, lower.CL, upper.CL),
    p.value,
    d, d.SE,
    `95% CI (d)` = sprintf("%.2f [%.2f, %.2f]", d, d.lower, d.upper)
  ) %>%
  dplyr::arrange(partner)

pretty_emm_H2
pretty_con_d_H2

# -----------------------------------------------
# H3: Interaction (difference-in-differences)
# -----------------------------------------------
# (self - other)_dissimilar  -  (self - other)_similar
# should be same as raw 2x2 interaction because now effects coded so raw not at reference level but grand mean
con_H3 <- contrast(con_H2, method = "pairwise", by = NULL)  # dissimilar - similar
con_d_H3   <- contrast_to_d(con_H3, sigma_total)

con_H3
con_d_H3

pretty_con_d_H3 <- con_d_H3 %>%
  dplyr::mutate(
    contrast = "(self–other)_dissimilar – (self–other)_similar"
  ) %>%
  dplyr::transmute(
    contrast,
    estimate, SE, df,t.ratio,
    `95% CI (raw)` = sprintf("%.3f [%.3f, %.3f]", estimate, lower.CL, upper.CL),
    p.value,
    d, d.SE,
    `95% CI (d)` = sprintf("%.3f [%.3f, %.3f]", d, d.lower, d.upper)
  )

pretty_con_d_H3

#Social Desirability

fit_sd <- lmer(
  rating ~ rating_type_c * partner_type_c * social_desirability +
    (1 + rating_type_c + partner_type_c || pid),
  data = d_long,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)),
  REML = FALSE
)
summary(fit_sd)

emmip(fit_sd, rating_type_c ~ social_desirability, CIs = TRUE)


# plot primary effects ---------------------------------------------------------------------

plot_means <- summary(emm_H2) %>%
  dplyr::mutate(
    rating  = dplyr::recode(as.character(rating_type_c),
                            `-0.5` = "Other", `0.5` = "Self"),
    partner = dplyr::recode(as.character(partner_type_c),
                            `-0.5` = "Dissimilar Partner", `0.5` = "Similar Partner")
  ) %>%
  dplyr::select(partner, rating, emmean, SE, df, lower.CL, upper.CL) %>%
  dplyr::arrange(partner, rating)

plot_df <- plot_means %>%
  transmute(
    sim_dissim = factor(partner, levels = c("Similar Partner","Dissimilar Partner")),
    self_other = factor(rating,  levels = c("Self","Other")),
    emmean, lower.CL, upper.CL
  )

## --- 1) Compute y-positions for brackets from model CIs ---
tops <- plot_df %>%
  group_by(sim_dissim) %>%
  summarise(y_top = max(upper.CL), .groups = "drop")

h2_offset <- 0.12  # gap above the highest CI in each group
h1_offset <- 0.35  # main-effect bracket height
h3_offset <- 0.55  # interaction bracket height

y_global <- max(tops$y_top) + h1_offset

## --- 2) Build star-only bracket data frames ---

# H2: self vs other within each partner group (two rows)
ann_h2 <- tops %>%
  transmute(
    x          = sim_dissim,       # the x groups
    group1     = "Self",           # left bar in each group
    group2     = "Other",          # right bar in each group
    y.position = y_top + h2_offset,
    label      = c("***","***")    # put your stars here
  )

# H1: similar vs dissimilar (one bracket across the two x groups)
ann_h1 <- tibble(
  group1      = "Similar Partner",
  group2      = "Dissimilar Partner",
  y.position  = y_global,
  label       = "**"
)

# H3: interaction bracket (just above H1)
ann_h3 <- tibble(
  group1      = "Similar Partner",
  group2      = "Dissimilar Partner",
  y.position  = y_global + (h3_offset - h1_offset),
  label       = "**"
)

## --- 3) Plot using ggpubr (keeps your original look) ---

bp <- ggbarplot(
  data = plot_df,
  x = "sim_dissim", y = "emmean", fill = "self_other",
  add = "none", position = position_dodge(0.9), width = 0.7
) +
  geom_errorbar(
    data = plot_df,
    aes(x = sim_dissim, ymin = lower.CL, ymax = upper.CL, group = self_other),
    position = position_dodge(0.9), width = 0.15, linewidth = 0.5
  ) +
  scale_fill_grey(start = .5, end = .9,
                  labels = c("Self: Own Interest", "Other: Estimate of Partner's Interest")) +
  labs(y = "Interest in Having a Conversation",
       x = "Type of Conversation Partner") +
  theme(
    text = element_text(size = 16),        # <-- fix here
    legend.title = element_blank(),
    legend.text  = element_text(size = 16),
    legend.key   = element_blank()
  ) +
  coord_cartesian(ylim = c(3, 6.5))  # same vertical range

# H2: within-group brackets (need x + group1/group2)
bp <- bp + stat_pvalue_manual(
  ann_h2,
  label = "label", y.position = "y.position",
  x = "x", xmin = "group1", xmax = "group2",
  tip.length = 0.02, bracket.size = 0.6, size = 7
)

# H1 and H3: across-group brackets (no x column; provide group1/group2 only)
bp <- bp +
  #stat_pvalue_manual(
  #  ann_h1, label = "label", y.position = "y.position",
  #  xmin = "group1", xmax = "group2",
  #  tip.length = 0.02, bracket.size = 0.6, size = 7
  #) +
  stat_pvalue_manual(
    ann_h3, label = "label", y.position = "y.position",
    xmin = "group1", xmax = "group2",
    tip.length = 0.02, bracket.size = 0.6, size = 7
  )

bp

fig1 <- bp


#fig1 <- ggpar(bp, ylim = c(3,6.5))

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

save_plot_multi(fig1,
                file_stem = "fig1",
                exts      = c("png", "pdf", "jpeg"),   # add/remove as needed
                width     = 10,
                height    = 8)


