# Title: Dissimilarity Pessimism Study 5
# Author(s): Gus Cooney; Erica Boothby
# Description: Dissimilarity pessimism among employees across organizationally relevant dimensions

# load libraries ----------------------------------------------------------

library(pacman)
p_load(tidyverse, ggplot, glue, corrr, broom, magrittr, ggExtra, lme4, ggeffects, lmerTest, 
       emmeans, pbkrtest, questionr, knitr, sjPlot, effects, report, brglm, scales, optimx, knitr,
       kableExtra)

# set constants -----------------------------------------------------------

DATA_PATH <- "~/Dropbox/Manuscripts/Homophily/Data"
PLOT_PATH <- "~/Dropbox/Manuscripts/Homophily/Plots"

# load data ---------------------------------------------------------------

d <- read_csv(glue("{DATA_PATH}/study5_final.csv"))
d

# demographics ---------------------------------------------------------------

#total
total <- length(unique(d$pid))

#age
age <- psych::describe(d$age)

#sex
sex <- table(d$sex)
names(sex)[names(sex) == "1"] = "Male"
names(sex)[names(sex) == "2"] = "Female"
names(sex)[names(sex) == "3"] = "Other"
prop.sex <- prop.table(sex)
df.sex <- data.frame(race=names(sex),count = as.numeric(sex), proportion = as.numeric(prop.sex))
df.sex <- df.sex[rev(order(df.sex$proportion)),]

#job level
job_level <- table(d$job_level)
names(job_level)[names(job_level) == "1"] = "Owner/Executive/C-Level"
names(job_level)[names(job_level) == "2"] = "Senior Management (Managing Director, VP, EVP, SVP, etc.)"
names(job_level)[names(job_level) == "3"] = "Middle Management (Supervisor, Manager, etc.)"
names(job_level)[names(job_level) == "4"] = "Intermediate (Associate, Senior Associate)"
names(job_level)[names(job_level) == "5"] = "Entry-Level (Assistant, Analyst, etc.)"
names(job_level)[names(job_level) == "6"] = "Other"
prop.job.level<- prop.table(job_level)
df.job.level <- data.frame(job_level=names(job_level),count = as.numeric(job_level), proportion = as.numeric(prop.job.level))
df.job.level<- df.job.level[rev(order(df.job.level$proportion)),]

#country
germany <- length(which(d$country == "65"))/total
usa <- length(which(d$country == "187"))/total

#demos
total
age
df.sex
df.job.level
germany
usa

# format data ------------------------------------------------------------------

d_long <- d %>%
  pivot_longer(
    cols = c(sim_self_overall, sim_other_overall, dissim_self_overall, dissim_other_overall,
             sim_self_careerstage, sim_other_careerstage, dissim_self_careerstage, dissim_other_careerstage,
             sim_self_position, sim_other_position, dissim_self_position, dissim_other_position,
             sim_self_division, sim_other_division, dissim_self_division, dissim_other_division,
             sim_self_industry, sim_other_industry, dissim_self_industry, dissim_other_industry,
             sim_self_sociocultural, sim_other_sociocultural, dissim_self_sociocultural, dissim_other_sociocultural
             ), 
    names_to = c("sim_dissim", "self_other",".value"),
    names_sep = "_",
    #values_to = NA
  )

#effects code for interpretation 
d_long <- d_long %>%
  mutate(
    rating_type_c = ifelse(self_other == "self", +0.5, -0.5),
    partner_type_c = ifelse(sim_dissim == "sim", +0.5, -0.5)
  )

d_long

write.csv(d_long, file = glue("{DATA_PATH}/z_study5_long_aggregate.csv"))

# modeling overall ----------------------------------------------------------------------

fit_rs <- lmer(
  overall ~ rating_type_c * partner_type_c +
    (1 + rating_type_c + partner_type_c || pid),
  data    = d_long,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)),
  REML    = FALSE
)
summary(fit_rs)
confint(fit_rs)


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


#-----------------------------
#rest of dimensions

dvs <- c("careerstage", "position", "division", "industry", "sociocultural")

dv_labels <- c(
  careerstage   = "Career Stage",
  position      = "Job Level",
  division      = "Division",
  industry      = "Industry",
  sociocultural = "Sociocultural"
)

# σ_total from model variance components (with random intercept + uncorrelated slopes)
sigma_total_from_model <- function(model) {
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
  
  sqrt(v_int + 0.25 * v_S + 0.25 * v_D + 0.0625 * v_SD + v_res)
}

contrast_to_d <- function(con_obj, sigma) {
  s <- summary(con_obj, infer = c(TRUE, TRUE))  # estimate, SE, df, CI
  s$d        <- s$estimate / sigma
  s$d.SE     <- s$SE / sigma
  s$d.lower  <- s$lower.CL / sigma
  s$d.upper  <- s$upper.CL / sigma
  s
}

analyze_dimension <- function(dv) {
  
  # Build formula DV ~ rating_type_c * partner_type_c
  form <- as.formula(
    glue("{dv} ~ rating_type_c * partner_type_c + (1 + rating_type_c + partner_type_c || pid)")
  )
  
  model <- lmer(
    form,
    data    = d_long,
    control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)),
    REML    = FALSE
  )
  
  sigma_total <- sigma_total_from_model(model)
  
  ## --- H1: Similar vs Dissim main effect (Table 3 analog) ---
  emm_H1 <- emmeans(
    model, ~ partner_type_c,
    pbkrtest.limit = 7856,
    at = list(partner_type_c = c(-0.5, 0.5))
  )
  
  emm_H1_sum <- summary(emm_H1)
  
  con_H1 <- contrast(emm_H1, list("sim - dissim" = c(-1, 1)), adjust = "none")
  con_H1_sum <- summary(con_H1, infer = c(TRUE, TRUE))
  con_H1_d   <- contrast_to_d(con_H1, sigma_total)
  
  # Extract means and H1 contrast
  sim_mean  <- emm_H1_sum$emmean[emm_H1_sum$partner_type_c ==  0.5]
  dis_mean  <- emm_H1_sum$emmean[emm_H1_sum$partner_type_c == -0.5]
  
  table3_row <- tibble(
    Dimension  = dv_labels[[dv]],
    Similar    = sim_mean,
    Dissimilar = dis_mean,
    Difference = con_H1_sum$estimate,
    CI         = sprintf("[%.2f, %.2f]", con_H1_sum$lower.CL, con_H1_sum$upper.CL),
    t          = con_H1_sum$t.ratio,
    df         = con_H1_sum$df,
    p          = con_H1_sum$p.value,
    d          = con_H1_d$d
  )
  
  ## --- H2: Self vs Other within dissimilar (Table 4 analog) ---
  emm_H2 <- emmeans(
    model, ~ rating_type_c | partner_type_c,
    pbkrtest.limit = 7856,
    at = list(rating_type_c = c(-0.5, 0.5),
              partner_type_c = c(-0.5, 0.5))
  )
  
  # Self vs Other contrast within each partner_type
  con_H2 <- contrast(emm_H2, list("self - other" = c(-1, 1)), by = "partner_type_c")
  con_H2_d <- contrast_to_d(con_H2, sigma_total)
  con_H2_sum <- summary(con_H2, infer = c(TRUE, TRUE))
  
  # We'll focus on the dissimilar condition (partner_type_c = -0.5) for H2
  dis_row   <- con_H2_sum[con_H2_sum$partner_type_c == -0.5, ]
  dis_row_d <- con_H2_d[con_H2_d$partner_type_c == -0.5, ]
  
  # Also need Self and Other emmeans in the dissimilar condition
  emm_dissim <- summary(emm_H2) %>%
    dplyr::filter(partner_type_c == -0.5) %>%
    dplyr::mutate(
      rating = dplyr::case_when(
        rating_type_c == -0.5 ~ "Other",
        rating_type_c ==  0.5 ~ "Self",
        TRUE                  ~ NA_character_
      )
    )  %>%
    arrange(rating)
  
  self_mean <- emm_dissim$emmean[emm_dissim$rating == "Self"]
  other_mean <- emm_dissim$emmean[emm_dissim$rating == "Other"]
  
  table4_row <- tibble(
    Dimension = dv_labels[[dv]],
    Self      = self_mean,
    Other     = other_mean,
    Difference = dis_row$estimate,
    CI        = sprintf("[%.2f, %.2f]", dis_row$lower.CL, dis_row$upper.CL),
    t         = dis_row$t.ratio,
    df        = dis_row$df,
    p         = dis_row$p.value,
    d         = dis_row_d$d
  )
  
  ## --- H3: Self–Other Gap Interaction (Table 5 analog) ---
  # difference-in-differences: (self-other)_dissim - (self-other)_similar
  con_H3 <- contrast(con_H2, method = "pairwise", by = NULL)
  con_H3_sum <- summary(con_H3, infer = c(TRUE, TRUE))
  con_H3_d   <- contrast_to_d(con_H3, sigma_total)
  
  table5_row <- tibble(
    Dimension = dv_labels[[dv]],
    Interaction = con_H3_sum$estimate,
    CI          = sprintf("[%.2f, %.2f]", con_H3_sum$lower.CL, con_H3_sum$upper.CL),
    t           = con_H3_sum$t.ratio,
    df          = con_H3_sum$df,
    p           = con_H3_sum$p.value,
    d           = con_H3_d$d
  )
  
  list(
    table3 = table3_row,
    table4 = table4_row,
    table5 = table5_row
  )
}

all_results <- map(dvs, analyze_dimension)

table3_all <- map_dfr(all_results, "table3")
table4_all <- map_dfr(all_results, "table4")
table5_all <- map_dfr(all_results, "table5")


# H1: Similar vs Dissimilar main effect
tab_H1 <- table3_all %>%
  transmute(
    Dimension,
    Hypothesis = "H1: Similar – Dissimilar",
    Estimate   = Difference,
    CI         = CI,
    t          = t,
    df         = df,
    p          = p,
    d          = d
  )

# H2: Self – Other in dissimilar condition
tab_H2 <- table4_all %>%
  transmute(
    Dimension,
    Hypothesis = "H2: Self – Other (Dissimilar Condition)",
    Estimate   = Difference,
    CI         = CI,
    t          = t,
    df         = df,
    p          = p,
    d          = d
  )

# H3: Interaction (gap_dissim – gap_sim)
tab_H3 <- table5_all %>%
  transmute(
    Dimension,
    Hypothesis = "H3: (Self–Other)_dissim – (Self–Other)_sim",
    Estimate   = Interaction,
    CI         = CI,
    t          = t,
    df         = df,
    p          = p,
    d          = d
  )

table_H123 <- bind_rows(tab_H1, tab_H2, tab_H3) %>%
  arrange(Dimension, Hypothesis)

table_for_print <- table_H123 %>%
  mutate(
    # Escape underscores
    Hypothesis = gsub("_", "\\\\_", Hypothesis),
    Estimate = sprintf("%.2f", Estimate),
    t        = sprintf("%.2f", t),
    df       = sprintf("%.0f", df),
    d        = sprintf("%.2f", d),
    p_val    = p,
    p        = case_when(
      p_val < .001 ~ "$\\textless .001$",
      p_val < .01  ~ "$\\textless .010$",
      p_val < .05  ~ sprintf("$%.3f$", p_val),
      TRUE         ~ "n.s."
    )
  ) %>%
  select(Dimension, Hypothesis, Estimate, CI, t, df, p, d)


# Start from table_H123 (numeric results for all domains × H1–H3)
table_for_print <- table_H123 %>%
  arrange(Dimension, Hypothesis) %>%
  mutate(
    # Domain will be printed only via group_rows; keep this column blank in the body
    Domain_print = "",
    
    # Indent and escape underscores in Hypothesis labels
    Hypothesis = paste0("\\hspace{2.0em}", gsub("_", "\\\\_", Hypothesis)),
    
    # Format numeric columns
    Estimate   = sprintf("%.2f", Estimate),
    t          = sprintf("%.2f", t),
    df         = sprintf("%.0f", df),
    d          = sprintf("%.2f", d),
    
    # Format p-values as LaTeX strings
    p_val      = p,
    p          = case_when(
      p_val < .001 ~ "$\\textless .001$",
      p_val < .01  ~ "$\\textless .010$",
      p_val < .05  ~ sprintf("$%.3f$", p_val),
      TRUE         ~ "n.s."
    )
  ) %>%
  select(Dimension, Domain_print, Hypothesis, Estimate, CI, t, df, p, d)

# indices for each domain's block (based on raw Dimension)
domain_groups <- split(seq_len(nrow(table_for_print)), table_for_print$Dimension)

kbl <- table_for_print %>%
  select(Domain_print, Hypothesis, Estimate, CI, t, df, p, d) %>%
  kable(
    format   = "latex",
    booktabs = TRUE,
    caption  = "Dissimilarity Pessimism Across Organizational Domains.",
    col.names = c("Domain", "Effect (Hypothesis)", "Estimate", "95\\% CI", "$t$", "df", "$p$", "$d$"),
    escape = FALSE,
    align  = c("l","l","r","c","r","r","c","r")
  ) %>%
  kable_styling(
    full_width = FALSE,
    position   = "center",
    latex_options = c("hold_position")
  )

# Add double midrule under the header row
kbl <- kbl %>%
  row_spec(0, extra_latex_after = "\\midrule")

# Add bold domain headings spanning the 3 rows for each domain
for (dom in names(domain_groups)) {
  rows <- domain_groups[[dom]]
  kbl <- kbl %>%
    group_rows(dom, min(rows), max(rows), bold = TRUE, italic = FALSE)
}

kbl


# plots ----------------------------------------------------------------------

# Mixed model for Sociocultural ratings
fit_sc <- lmer(
  sociocultural ~ rating_type_c * partner_type_c +
    (1 + rating_type_c + partner_type_c || pid),
  data    = d_long,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)),
  REML    = FALSE
)

emm_sc <- emmeans(
  fit_sc, ~ rating_type_c | partner_type_c,
  pbkrtest.limit = 7856,
  at = list(rating_type_c   = c(-0.5, +0.5),
            partner_type_c  = c(-0.5, +0.5))
)

plot_means_sc <- summary(emm_sc) %>%
  mutate(
    rating  = dplyr::recode(as.character(rating_type_c),
                     `-0.5` = "Other", `0.5` = "Self"),
    partner = dplyr::recode(as.character(partner_type_c),
                     `-0.5` = "Dissimilar Sociocultural Background",
                     `0.5`  = "Similar Sociocultural Background")
  ) %>%
  select(partner, rating, emmean, SE, df, lower.CL, upper.CL) %>%
  arrange(partner, rating)

plot_df_sc <- plot_means_sc %>%
  transmute(
    sim_dissim = factor(partner,
                        levels = c("Similar Sociocultural Background",
                                   "Dissimilar Sociocultural Background")),
    self_other = factor(rating, levels = c("Self", "Other")),
    emmean, lower.CL, upper.CL
  )

#new comparison for this graph - self to self dissimilar
emm_sc_full <- emmeans(
  fit_sc,
  ~ rating_type_c * partner_type_c,
  at = list(
    rating_type_c  = c(-0.5, 0.5),
    partner_type_c = c(-0.5, 0.5)
  )
)
emm_sc_full

con_self_sc <- contrast(
  emm_sc_full,
  list("Self: dissimilar - similar" = c(0, 0, -1, 1))
)
summary(con_self_sc, infer = TRUE)

# Brackets

#old places

tops_sc <- plot_df_sc %>%
  group_by(sim_dissim) %>%
  summarise(y_top = max(upper.CL), .groups = "drop")

h2_offset <- 0.12
h3_offset <- 0.55
y_global  <- max(tops_sc$y_top) + 0.35

# H2: Self vs Other within each partner type (both significant for Sociocultural)
ann_h2_sc <- tops_sc %>%
  transmute(
    x          = sim_dissim,
    group1     = "Self",
    group2     = "Other",
    y.position = y_top + h2_offset,
    label      = c("***", "***")   # you can use "**" or "n.s." if you want to be precise
  )

# H3: interaction (self–other gap_dissim – gap_sim); significant for Sociocultural
ann_h3_sc <- tibble(
  group1      = "Similar Sociocultural Background",
  group2      = "Dissimilar Sociocultural Background",
  y.position  = y_global + (h3_offset - 0.35),
  label       = "***"
)

#new bracket position

# Find the top of the self bars (using their CIs)
# max of Self bars (CIs) for a common bracket height
self_tops_sc <- plot_df_sc %>%
  filter(self_other == "Self") %>%
  summarise(y_top = max(upper.CL))

y_bracket <- self_tops_sc$y_top + 0.40  # same height for both brackets

# numeric x positions of the factor levels
x_sim  <- which(levels(plot_df_sc$sim_dissim) ==
                  "Similar Sociocultural Background")
x_dsim <- which(levels(plot_df_sc$sim_dissim) ==
                  "Dissimilar Sociocultural Background")

dodge_width <- 0.7
offset      <- 0  # 0.35

# bar centers
self_sim_center    <- x_sim  - .17
self_dsim_center   <- x_dsim - .12
other_dsim_center  <- x_dsim + .17

self_tops_sc <- plot_df_sc %>%
  dplyr::filter(self_other == "Self") %>%
  summarise(y_top = max(upper.CL))

y_bracket <- self_tops_sc$y_top + 0.25

# bracket 1: Self(similar) vs Self(dissimilar)
bracket_self_self <- tibble::tibble(
  xstart = self_sim_center,
  xend   = self_dsim_center,
  y      = y_bracket,
  label  = "***"
)

# bracket 2: Self vs Other in dissimilar condition
bracket_self_other_dsim <- tibble::tibble(
  xstart = self_dsim_center,
  xend   = other_dsim_center,
  y      = y_bracket,
  label  = "***"
)

gap1 <- 0.065   # tweak as needed (~1/20 of a bar width feels nice)
gap2 <- -0.04   # tweak as needed (~1/20 of a bar width feels nice)


bracket_self_self  <- bracket_self_self  %>%
  mutate(xend   = xend   - gap1)

bracket_self_other_dsim <- bracket_self_other_dsim %>%
  mutate(xstart = xstart + gap2)


## --- 3) Plot using ggpubr (keeps your original look) ---

bp_sc_simple <- ggbarplot(
  data = plot_df_sc,
  x = "sim_dissim", y = "emmean", fill = "self_other",
  add = "none", position = position_dodge(0.7), width = 0.7
) +
  geom_errorbar(
    data = plot_df_sc,
    aes(x = sim_dissim, ymin = lower.CL, ymax = upper.CL, group = self_other),
    position = position_dodge(0.7), width = 0.15, linewidth = 0.5
  ) +
  scale_fill_grey(start = .5, end = .9,
                  labels = c("Self: Own Interest", "Other: Estimate of Partner's Interest")) +
  labs(y = "Expected Interest in the Conversation",
       x = "Type of Conversation Partner") +
  theme(
    text = element_text(size = 16),
    legend.title = element_blank(),
    legend.text  = element_text(size = 16),
    legend.key   = element_blank()
  ) +
  coord_cartesian(ylim = c(4.0, 6.7))  # tweak as desired


bp_sc_simple <- bp_sc_simple +
  # Self(similar) vs Self(dissim)
  geom_segment(
    data = bracket_self_self,
    aes(x = xstart, xend = xend, y = y, yend = y),
    linewidth = 0.6
  ) +
  geom_segment(
    data = bracket_self_self,
    aes(x = xstart, xend = xstart, y = y, yend = y - 0.05),
    linewidth = 0.6
  ) +
  geom_segment(
    data = bracket_self_self,
    aes(x = xend, xend = xend, y = y, yend = y - 0.05),
    linewidth = 0.6
  ) +
  geom_text(
    data = bracket_self_self,
    aes(x = (xstart + xend) / 2, y = y + 0.05, label = label),
    size = 5
  ) +
  # Self vs Other within dissimilar
  geom_segment(
    data = bracket_self_other_dsim,
    aes(x = xstart, xend = xend, y = y, yend = y),
    linewidth = 0.6
  ) +
  geom_segment(
    data = bracket_self_other_dsim,
    aes(x = xstart, xend = xstart, y = y, yend = y - 0.05),
    linewidth = 0.6
  ) +
  geom_segment(
    data = bracket_self_other_dsim,
    aes(x = xend, xend = xend, y = y, yend = y - 0.05),
    linewidth = 0.6
  ) +
  geom_text(
    data = bracket_self_other_dsim,
    aes(x = (xstart + xend) / 2, y = y + 0.05, label = label),
    size = 5
  )

fig_sociocultural <- bp_sc_simple
fig_sociocultural 


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

save_plot_multi(fig_sociocultural,
                file_stem = "fig6",
                exts      = c("png", "pdf", "jpeg"),   # add/remove as needed
                width     = 10,
                height    = 8)







