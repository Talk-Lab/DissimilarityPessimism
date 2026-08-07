# Title: Dissimilarity Pessimism Study 3a
# Author(s): Gus Cooney; Erica Boothby
# Description: Dissimilarity pessimism & homophilous choice

# load libraries ----------------------------------------------------------

library(pacman)
p_load(tidyverse, ggplot2, glue, corrr, broom, magrittr, ggExtra, lme4, ggeffects, 
       lmerTest, emmeans, pbkrtest, questionr, knitr, sjPlot, effects, report, ltm, 
       glm.predict, ggpubr, rstatix, margins)

# set constants -----------------------------------------------------------

DATA_PATH <- "~/Dropbox/Manuscripts/Homophily/Data"
PLOT_PATH <- "~/Dropbox/Manuscripts/Homophily/Plots"

# load data ---------------------------------------------------------------

d <- read_csv(glue("{DATA_PATH}/study3a_final.csv"))
head(d)

# exclusions ---------------------------------------------------------------
subset <- subset(d, exclude == 0)
d <- subset

# demographics ---------------------------------------------------------------

df <- d

sex_map  <- c("1"="Male","2"="Female","3"="Other")
race_map <- c(
  "1"="African- or Caribbean-American", "2"="Hispanic/Latino",
  "3"="Asian-American", "4"="Native American",
  "5"="White/Caucasian", "6"="Other"
)


df_completed <- df %>%
  mutate(
    age_num  = suppressWarnings(as.numeric(age)),
    sex_chr  = dplyr::recode(as.character(sex),  !!!sex_map,  .default = NA_character_),
    race_chr = dplyr::recode(as.character(race), !!!race_map, .default = NA_character_)
  )


# 1) Session and person counts --------------------------
n_sessions <- nrow(df_completed)
n_people   <- n_distinct(df_completed$participant_id)

sessions_per_person <- df_completed %>%
  dplyr::count(participant_id, name = "n_sessions")

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

# composite & manipulation check -------------------------------------------------

#composite
d.sim.self <- d[,c("sim_self1","sim_self2","sim_self3")]
cronbach.alpha(d.sim.self, na.rm=T)
d.sim.other <- d[,c("sim_other1","sim_other2","sim_other3")]
cronbach.alpha(d.sim.other, na.rm=T)
d.dif.self <- d[,c("dissim_self1","dissim_self2","dissim_self3")]
cronbach.alpha(d.dif.self, na.rm=T)
d.dif.other <- d[,c("dissim_other1","dissim_other2","dissim_other3")]
cronbach.alpha(d.dif.other, na.rm=T)

#manipulation
d.manipulation <- d %>%
  pivot_longer(cols = c(similar_same, similar_diff),
               names_to = "similarity",
               values_to = "rating")

t.test(d$similar_same, d$similar_diff, paired = TRUE)

d.manipulation %>%
  group_by(similarity) %>%
  summarise(
    n = n(),
    mean = mean(rating),
    sd = sd(rating),
    se = sd / sqrt(n),
    lower = mean - qt(0.975, n-1) * se,
    upper = mean + qt(0.975, n-1) * se
  )

# format data ---------------------------------------------------------------

d_long <- d %>%
  pivot_longer(
    cols = c(sim_self, sim_other, dissim_self, dissim_other), 
    names_to = c("sim_dissim", "self_other"),
    names_sep = "_",
    values_to = "rating"
  ) %>%
  mutate(
    pid = participant_id
  )

d_long

d_long <- d_long |>
  dplyr::mutate(
    rating_type_c   = ifelse(self_other  == "self",  +0.5, -0.5),
    partner_type_c  = ifelse(sim_dissim  == "sim",   +0.5, -0.5),
  )

d_long$rating <- as.numeric(d_long$rating)

#write.csv(d_long, file = glue("{DATA_PATH}/z_study3a_long_aggregate.csv"))


# modeling ---------------------------------------------------------------------

fit_rs <- lmer(
  rating ~ rating_type_c * partner_type_c +
    (1 + rating_type_c + partner_type_c || pid),
  data    = d_long,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)),
  REML    = FALSE
)
summary(fit_rs)
confint(fit_rs)


model <- fit_rs

# simr sensitivity analysis
#for h3
term <- "rating_type_c:partner_type_c"

sens <- function(model, term, beta, nsim = 200) {
  m <- model
  fixef(m)[term] <- beta
  ps <- powerSim(m, fixed(term, "t"), nsim = nsim, progress = FALSE)
  summary(ps)$mean
}

betas <- c(-0.06, -0.08, -0.10)
out <- data.frame(beta = betas,
                  power = sapply(betas, \(b) sens(model, term, b)))
out
#for h2
#use 2.8 x SE


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

h2_h3_sds <- d_long %>%
  group_by(sim_dissim, self_other) %>%
  summarise(M = mean(rating), SD = sd(rating), n_obs = n(), .groups = "drop")


pretty_emm_H2
pretty_con_d_H2
h2_h3_sds

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

# plot primary effects ---------------------------------------------------------------------

plot_means <- summary(emm_H2) %>%
  dplyr::mutate(
    rating  = dplyr::recode(as.character(rating_type_c),
                            `-0.5` = "Other", `0.5` = "Self"),
    partner = dplyr::recode(as.character(partner_type_c),
                            `-0.5` = "Dissimilar Personality", `0.5` = "Similar Personality")
  ) %>%
  dplyr::select(partner, rating, emmean, SE, df, lower.CL, upper.CL) %>%
  dplyr::arrange(partner, rating)

plot_df <- plot_means %>%
  transmute(
    sim_dissim = factor(partner, levels = c("Similar Personality","Dissimilar Personality")),
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
  group1      = "Similar Personality",
  group2      = "Dissimilar Personality",
  y.position  = y_global,
  label       = "***"
)

# H3: interaction bracket (just above H1)
ann_h3 <- tibble(
  group1      = "Similar Personality",
  group2      = "Dissimilar Personality",
  y.position  = y_global + (h3_offset - h1_offset),
  label       = "***"
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
  labs(y = "Expected Positivity of the Conversation",
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

fig3 <- bp


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

save_plot_multi(fig3,
                file_stem = "fig3",
                exts      = c("png", "pdf", "jpeg"),   # add/remove as needed
                width     = 10,
                height    = 8)


# choice modeling ---------------------------------------------------------------------

d$dissim_other <- as.numeric(d$dissim_other)
d$dissim_self <- as.numeric(d$dissim_self)
d$sim_other <- as.numeric(d$sim_other)
d$sim_self <- as.numeric(d$sim_self)

m0 <- glm(dissim_choice ~ 1, data = d, family = binomial)
summary(m0)

m1 <- glm(dissim_choice ~ dissim_other, data = d, family = binomial)
summary(m1)
confint(m1)

anova(m1, test="LRT")

# coefficient
beta_m1 <- coef(m1)["dissim_other"]

# 95% CI for beta (profile likelihood)
ci_beta_m1 <- confint(m1)["dissim_other", ]

# convert to OR
OR_m1 <- exp(beta_m1)
ci_OR_m1 <- exp(ci_beta_m1)

# Cohen's d analogue
d_m1 <- beta_m1 / 1.81

list(
  beta = beta_m1,
  beta_CI = ci_beta_m1,
  OR = OR_m1,
  OR_CI = ci_OR_m1,
  d = d_m1
)

#model 2: with self control
m2 <- glm(dissim_choice ~ dissim_self + dissim_other, data = d, family = binomial)
summary(m2)
confint(m2)
anova(m2, test = "LRT")
#anova(m1, m2, test = "LRT")

# betas
beta_m2_other <- coef(m2)["dissim_other"]
beta_m2_self  <- coef(m2)["dissim_self"]

# 95% CIs (profile likelihood)
ci_beta_m2_other <- confint(m2)["dissim_other", ]
ci_beta_m2_self  <- confint(m2)["dissim_self", ]

# ORs
OR_m2_other  <- exp(beta_m2_other)
OR_m2_self   <- exp(beta_m2_self)

ci_OR_m2_other <- exp(ci_beta_m2_other)
ci_OR_m2_self  <- exp(ci_beta_m2_self)

# effect size d
d_m2_other <- beta_m2_other / 1.81
d_m2_self  <- beta_m2_self  / 1.81

list(
  m2_other = list(
    beta = beta_m2_other,
    beta_CI = ci_beta_m2_other,
    OR = OR_m2_other,
    OR_CI = ci_OR_m2_other,
    d = d_m2_other
  ),
  m2_self = list(
    beta = beta_m2_self,
    beta_CI = ci_beta_m2_self,
    OR = OR_m2_self,
    OR_CI = ci_OR_m2_self,
    d = d_m2_self
  )
)


#reciprocity:

reciprocity <- lm(dissim_self ~ dissim_other, data = d)
summary(reciprocity)

cor(d$dissim_self, d$dissim_other)


#why is the average AME different than the 0.05-0.06 jump in plot 2???

# m1: dissim_choice ~ dissim_other
ame_m1 <- margins(m1)            # by default, AME for all predictors
summary(ame_m1)

# m2: dissim_choice ~ dissim_other + dissim_self
ame_m2 <- margins(m2)
summary(ame_m2)

# Separate AMEs for each predictor
ame_m2_other <- margins(m2, variables = "dissim_other")
ame_m2_self  <- margins(m2, variables = "dissim_self")

summary(ame_m2_other)
summary(ame_m2_self)

# predicted probabilities for each person at their actual values
p_hat <- predict(m2, type = "response")

# slope for dissim_other at each person's point (dp/dx_i)
beta_other <- coef(m2)["dissim_other"]
ME_i <- beta_other * p_hat * (1 - p_hat)

# AME is the mean of these marginal effects
AME <- mean(ME_i)
AME


df_me <- data.frame(
  p_hat = p_hat,
  ME_i  = ME_i
)

#ame average line in red
ggplot(df_me, aes(x = p_hat, y = ME_i)) +
  geom_point(alpha = 0.2) +
  geom_hline(yintercept = AME, color = "red", linetype = "dashed") +
  labs(
    x = "Predicted probability (p_hat)",
    y = "Marginal effect of dissim_other (dp/dx)",
    title = "Per-person marginal effects and the AME (red line)"
  ) +
  theme_bw()


# sequence of dissim_other from 1 to 7
xvals <- 1:7

# hold dissim_self at its mean
self_mean <- mean(d$dissim_self, na.rm = TRUE)

# predicted probabilities for an "average" participant
eta_avg <- coef(m2)["(Intercept)"] + beta_other * xvals + coef(m2)["dissim_self"] * self_mean
p_avg  <- plogis(eta_avg)  # logistic(eta)

data_plot <- data.frame(
  dissim_other = xvals,
  p_avg = p_avg
)

ggplot(data_plot, aes(x = dissim_other, y = p_avg)) +
  geom_point() +
  geom_line() +
  ylim(0, 1) +
  labs(
    x = "dissim_other",
    y = "Predicted probability for 'average' self-interest",
    title = "Predicted probability curve for an average participant"
  ) +
  theme_bw()


#force ggpredict to predict outside average! Show jumps are smaller at extremes

p_both <- ggpredict(
  m2,
  terms = c("dissim_other[1:7]", "dissim_self[1,7]")
)

p_both

df_both <- as.data.frame(p_both)

ggplot(df_both, aes(x = x, y = predicted, color = as.factor(group))) +
  geom_line() +
  geom_point() +
  ylim(0, 1) +
  labs(
    x = "dissim_other",
    y = "Predicted probability of choosing dissimilar partner",
    color = "dissim_self",
    title = "Logistic curve at low vs high self-interest"
  ) +
  theme_bw()

#clogit
d_alt <- d %>%
  transmute(participant_id,
            sim_choice,  dissim_choice,
            sim_other,   dissim_other,
            sim_self,    dissim_self) %>%    # perceived similarity 
  pivot_longer(
    cols = c(sim_choice, dissim_choice,
             sim_other, dissim_other,
             sim_self, dissim_self),
    names_to = c("alt", ".value"),             
    names_pattern = "(sim|dissim)_(.*)"
  ) %>%
  mutate(
    # alt-specific predictors: rename for clarity
    choice            = choice,
    other_interest    = other,
    self_interest     = self,
    alt               = factor(alt, levels = c("sim","dissim"))
  )

library(survival)
cl <- clogit(choice ~ other_interest + self_interest + strata(participant_id),
             data = d_alt)
summary(cl)

#summary
tibble(
  log_OR                 = coef(m1)) %>% 
  cbind(confint(model.dif.1)) %>% 
  dplyr::rename(lower_log_OR    = "2.5 %",
                upper_log_OR    = "97.5 %") %>%
  cbind(odds.ratio(model.dif.1)) %>% 
  dplyr::rename(lower_OR        = "2.5 %",
                upper_OR        = "97.5 %") %>% 
  mutate(percent_change  = ifelse(OR < 1, (1/OR - 1)*-100, (OR - 1)*100  ),
         lower_percent_change = ifelse(lower_OR < 1, (1/lower_OR - 1)*-100, (lower_OR - 1)*100  ),
         upper_percent_change = ifelse(upper_OR < 1, (1/upper_OR - 1)*-100, (upper_OR - 1)*100  )) %>%
  mutate_if(is.numeric, ~round(., 3)) %>% 
  t() %>% 
  kable()

tibble(
  log_OR                 = coef(m2)) %>% 
  cbind(confint(model.dif.2)) %>% 
  dplyr::rename(lower_log_OR    = "2.5 %",
                upper_log_OR    = "97.5 %") %>%
  cbind(odds.ratio(model.dif.2)) %>% 
  dplyr::rename(lower_OR        = "2.5 %",
                upper_OR        = "97.5 %") %>% 
  mutate(percent_change  = ifelse(OR < 1, (1/OR - 1)*-100, (OR - 1)*100  ),
         lower_percent_change = ifelse(lower_OR < 1, (1/lower_OR - 1)*-100, (lower_OR - 1)*100  ),
         upper_percent_change = ifelse(upper_OR < 1, (1/upper_OR - 1)*-100, (upper_OR - 1)*100  )) %>%
  mutate_if(is.numeric, ~round(., 3)) %>% 
  t() %>% 
  kable()


#plots
plot1 <- ggpredict(m1, terms = "dissim_other[1,2,3,4,5,6,7]") %>% 
  ggplot(aes(x = x, y = predicted, label = round(predicted, digit = 2))) +
  theme_bw() +
  geom_point() +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = .1) +
  geom_label(hjust = -.3) +
  #scale_x_continuous(expand = c(0.05,.15)) +
  scale_x_continuous(limits = c(.9,8)) +
  theme(
    text = element_text(size = 16)) +        
  ylim(0,1) +
  labs(x = "Perceived Positivity of Dissimilar Others' Expectations",
       y = "Probability of Choosing a Dissimilar Partner",
       title = "Expectations and Partner Choice")
plot1

### Model 2 - Model 1 + controlling for one's own preference


plot2 <- ggpredict(m2, terms = "dissim_other[1,2,3,4,5,6,7]") %>% 
  ggplot(aes(x = x, y = predicted, label = round(predicted, digit = 2))) +
  theme_bw() +
  geom_point() +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = .1) +
  geom_label(hjust = -.3) +
  #scale_x_continuous(expand = c(0.05,.50)) +
  ylim(0,1) +
  scale_x_continuous(limits = c(.9,8)) +
  theme(
    text = element_text(size = 16)) +    
  labs(x = "Perceived Positivity of Dissimilar Others' Expectations",
       y = "Probability of Choosing a Dissimilar Partner",
       title = "[...] Controlling for One's Own Expectations")
plot2
#predicted probabilities with self at mean
#at the mean level of self, if you go from the midpoint of the scale to the highest scale, 
#you go from 50% probability to 67% probaiilty of choosing the dissimlar partner, 


fig4 <- grid.arrange(plot1, plot2, ncol=2)
fig4

save_plot_multi(fig4,
                file_stem = "fig4",
                exts      = c("png", "pdf", "jpeg"),   # add/remove as needed
                width     = 12,
                height    = 6)
