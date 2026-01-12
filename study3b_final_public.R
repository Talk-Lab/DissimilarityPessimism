# Title: Dissimilarity Pessimism Study 3b
# Author(s): Gus Cooney; Erica Boothby
# Description: Dissimilarity pessimism & homophilous choice 2

# load libraries ----------------------------------------------------------

library(pacman)
p_load(tidyverse, glue, broom, readr, lme4, lmerTest, emmeans, effectsize, rempsyc, magrittr, 
       performance, RSA, lavaan, psych, flextable, ltm, rstatix, regclass, lm.beta)

# set constants -----------------------------------------------------------

DATA_PATH <- "~/Dropbox/Manuscripts/Homophily/Data"
PLOT_PATH <- "~/Dropbox/Manuscripts/Homophily/Plots"

# load data ---------------------------------------------------------------

names <- read_csv(glue("{DATA_PATH}/study3b_final.csv"), n_max = 0) %>% names()
raw <- read_csv(glue("{DATA_PATH}/study3b_final.csv"), col_names = names, skip = 2)
d <- raw
d
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
d.sim.self <- d[,c("sim_self1","sim_self2")]
cronbach.alpha(d.sim.self, na.rm=T)
d.sim.other <- d[,c("sim_other1","sim_other2")]
cronbach.alpha(d.sim.other, na.rm=T)
d.dif.self <- d[,c("dissim_self1","dissim_self2")]
cronbach.alpha(d.dif.self, na.rm=T)
d.dif.other <- d[,c("dissim_other1","dissim_other2")]
cronbach.alpha(d.dif.other, na.rm=T)

#manipulation
d.manipulation <- d[,c("similar_same", "similar_diff")]
d.manipulation <- d.manipulation %>%
  pivot_longer(cols = c(similar_same, similar_diff),
               names_to = "similarity",
               values_to = "rating")

t.test(d$similar_same, d$similar_diff, paired = TRUE)

d.manipulation %>%
  group_by(similarity) %>%
  get_summary_stats(rating, show = c("mean","sd","ci"))

# format data -------------------------------------------------------------------

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

write.csv(d_long, file = glue("{DATA_PATH}/z_study3b_long_aggregate.csv"))


# modeling ---------------------------------------------------------------------
#same as previous Study 3a - see that for additional notes

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


con_H1 <- contrast(emm_H1, method = list("sim - dissim" = c(-1, 1)), adjust = "none")
con_d_H1   <- contrast_to_d(con_H1, sigma_total)

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

# choice modeling ---------------------------------------------------------------------

d_regress <- d %>%
  mutate(self_minus_other_interest_gap = dissim_self - dissim_other,
         dissimilar_minus_similar_choice_gap = dissim_choice-sim_choice) %>%
  rename(selfs_interest_in_dissimilar_other = dissim_self, 
         perceived_dissimilar_others_interest_in_self = dissim_other,
         selfs_interest_in_similar_other = sim_self, 
         perceived_similar_others_interest_in_self = sim_other,
         desire_to_choose_dissimlar_partner = dissim_choice,
         desire_to_choose_similar_partner = sim_choice)
 
#1. other predicts choice
m1 <- lm(desire_to_choose_dissimlar_partner ~ perceived_dissimilar_others_interest_in_self, 
    data = d_regress)
summary(m1)
confint(m1)
effectsize(m1)
t_to_d(t = 19.078, df_error = 748)


#2. other predicts choice controlling for self
m2 <- lm(desire_to_choose_dissimlar_partner ~ selfs_interest_in_dissimilar_other + 
     perceived_dissimilar_others_interest_in_self, data = d_regress) 
summary(m2)
confint(m2)
effectsize(m2)
# Cohen's d self's interest
#t_to_d(t = 17.746, df_error = 747)

# Cohen's d for perceived other's interest
t_to_d(t = 3.491, df_error = 747)

#extra analysis
m.extra.1 <- lm(dissimilar_minus_similar_choice_gap ~ selfs_interest_in_dissimilar_other + 
       perceived_dissimilar_others_interest_in_self, data = d_regress)
summary(m.extra.1)
confint(m.extra.1)
  
m.extra.2 <- lm(dissimilar_minus_similar_choice_gap ~ selfs_interest_in_dissimilar_other + 
                  perceived_dissimilar_others_interest_in_self + selfs_interest_in_similar_other + 
                  perceived_similar_others_interest_in_self, data = d_regress)
summary(m.extra.2)
confint(m.extra.2)

# CRA ---------------------------------------------------------------------
# from Humberg et al., 2018

# simple difference score
m.dif.score <- lm(desire_to_choose_dissimlar_partner ~ self_minus_other_interest_gap, 
         data = d_regress)
summary(m.dif.score)
confint(m.dif.score)
t_to_d(t = 6.53, df_error = 748)


d_cra <- d %>%
  rename(S_dissim = dissim_self, R_dissim = dissim_other,S_sim = sim_self, R_sim = sim_other)%>%
  mutate(gap_self_minus_other = S_dissim - R_dissim,
         dissim_sim_choice_delta = dissim_choice-sim_choice)
# recommendation to standardize variables for CRA, but same scale already

vif.test <- lm(dissim_choice ~ S_dissim + R_dissim, data=d_cra)
VIF(vif.test)
# [...] We recommend a pre-analysis to test for potential multicollinearity problems 
# (e.g., by applying the variance inflation factor; VIF; Fox, 2016).

d_cra %>% 
  lm(dissim_choice ~ S_dissim + R_dissim, .) %>% 
  check_model()

model <- 'dissim_choice ~ 1 + c1*S_dissim + c2*R_dissim
a1 := c1+c2
a3 := c1-c2
c1_times_2 := 2*c1
c2_times_2 := 2*c2
minus_c1_times_2 := -2*c1
minus_c2_times_2 := -2*c2'

regression <- sem(model, data=d_cra)
summary(regression)
estimates <- parameterEstimates(regression)
estimates

subset(estimates, label == "c1")[c("est","se","pvalue", "ci.lower", "ci.upper")]

# Effect of R, indicated by c2:

subset(estimates, label == "c2")[c("est","se","pvalue", "ci.lower", "ci.upper")]

# Effect of SE

if (subset(estimates, label == "a3")["est"]>=0 & subset(estimates, label == "a1")["est"]>0){
  abs <- round(subset(estimates, label == "minus_c2_times_2")["est"],2)
  se <- round(subset(estimates, label == "minus_c2_times_2")["se"],2)
  a3_is <- "positive"
  
  # compute one-tailed p-value of abs for the hypothesis abs > 0, 
  # depending on the tail in which abs is positioned
  if (abs >= 0){pvalue <- round(subset(estimates, label == "minus_c2_times_2")["pvalue"]/2, 5)}
  if (abs < 0){pvalue <- round(1 - (subset(estimates, label == "minus_c2_times_2")["pvalue"]/2), 5)}
}


if (subset(estimates, label == "a3")["est"]<0 & subset(estimates, label == "a1")["est"]>0){
  abs <- round(subset(estimates, label == "minus_c1_times_2")["est"],2)
  se <- round(subset(estimates, label == "minus_c1_times_2")["se"],2)
  a3_is <- "negative"
  
  # compute one-tailed p-value of abs for the hypothesis abs > 0, 
  # depending on the tail in which abs is positioned
  if (abs >= 0){pvalue <- round(subset(estimates, label == "minus_c1_times_2")["pvalue"]/2, 5)}
  if (abs < 0){pvalue <- round(1 - (subset(estimates, label == "minus_c1_times_2")["pvalue"]/2), 5)}
}


if (subset(estimates, label == "a3")["est"]>=0 & subset(estimates, label == "a1")["est"]<0){
  abs <- round(subset(estimates, label == "c1_times_2")["est"],2)
  se <- round(subset(estimates, label == "c1_times_2")["se"],2)
  a3_is <- "positive"
 
  # compute one-tailed p-value of abs for the hypothesis abs > 0, 
  # depending on the tail in which abs is positioned
  if (abs >= 0){pvalue <- round(subset(estimates, label == "c1_times_2")["pvalue"]/2, 5)}
  if (abs < 0){pvalue <- round(1 - (subset(estimates, label == "c1_times_2")["pvalue"]/2), 5)}
  
}

if (subset(estimates, label == "a3")["est"]<0 & subset(estimates, label == "a1")["est"]<0){
  abs <- round(subset(estimates, label == "c2_times_2")["est"],2)
  se <- round(subset(estimates, label == "c2_times_2")["se"],2)
  a3_is <- "negative"
  
  # compute one-tailed p-value of abs for the hypothesis abs > 0, 
  # depending on the tail in which abs is positioned
  if (abs >= 0){pvalue <- round(subset(estimates, label == "c2_times_2")["pvalue"]/2, 5)}
  if (abs < 0){pvalue <- round(1 - (subset(estimates, label == "c2_times_2")["pvalue"]/2), 5)}
  
  }
#code is from original paper https://osf.io/hcsg6
#slightly different tests with no "minus_c1_times_2" in 2021 narcisism paper?

ci.lower <- round(abs -1.65*se, 3)

print(c("abs"=abs,se,pvalue,"a3_is"=a3_is, ci.lower = ci.lower))
#ci.upper is infinity
c1 <- as.tibble(subset(estimates, label == "c1")[c("est", "se", "pvalue", "ci.lower", "ci.upper")])
c2 <- as.tibble(subset(estimates, label == "c2")[c("est","se", "pvalue", "ci.lower", "ci.upper")])
abs <- as.tibble(c(abs, se, pvalue, ci.lower=ci.lower, ci.upper = "inf"))
abs <- abs %>%
  rename (ci.lower = ci.lower.est)

tidy_cra <- rbind(c1,c2,abs) %>%
  dplyr::select(-c("se")) %>%  
  rename(p = pvalue,
          Coefficient = est,
          CI_lower = ci.lower,
          CI_upper = ci.upper) %>%
    mutate(Criteria = c("c1", "c2","abs")) %>%
    mutate(Interpretation = c("Effect of Self's Interest in Talking to Dissimilar Other", 
                              "Effect of Perceptions of Dissimilar Other's Interest in Talking", 
                              "Effect of the Magnitude of the Difference Between Self and Other")) %>%
    relocate(Criteria, .before = Coefficient) %>%
    relocate(Interpretation, .before = Coefficient)
          
t_cra <- nice_table(tidy_cra, 
           title = c("Table 1. Dissimilarity Pessimism and Homophilous Choice.", "CRA Results"),
           note = c("Coefficients are unstandardized",
                    "The p values for c1 and c2 are two-tailed, and the p value for abs is one-tailed.",
                    "The 'c2' parameter reflects whether people's beliefs about dissimilar others' interest predicts
                    their choice of whether to talk to dissimilar others. People's mistakenly pessimistic beliefs about 
                    dissimilar others' interest will thus lead to fewer conversations with dissimlar others, and
                    this is one way dissimilarity pessimism is related to homophily (Hypothesis 5a).",
                    "The significance of 'abs' parameter reflects whether the magnitude of the self-other dissimilarity pessimism gap
                    (i.e., the difference between people's own interest in talking to a dissimilar other and people's
                    estimates of dissimlar others' interest in talking to them) predicts homophilous choice (Hypothesis 5b).",
                    "regression equation: 
                    willingness to talk to dissimilar other ~ co + c1*personal_interest + c2*perceived_others_interest + e"))

t_cra
theme_apa(t_cra)



# --- helpers ---------------------------------------------------------------

get_row <- function(estimates, lbl) {
  subset(estimates, label == lbl, select = c("est","se","pvalue","ci.lower","ci.upper"))
}

compute_abs <- function(estimates) {
  a1 <- get_row(estimates, "a1"); a3 <- get_row(estimates, "a3")
  a1_est <- as.numeric(a1$est);   a3_est <- as.numeric(a3$est)
  
  # Choose the linear combo per Humberg et al. (OSF Material G)
  if (a3_est >= 0 & a1_est > 0) {
    row  <- get_row(estimates, "minus_c2_times_2"); combo <- "-2 * c2"
  } else if (a3_est < 0 & a1_est > 0) {
    row  <- get_row(estimates, "minus_c1_times_2"); combo <- "-2 * c1"
  } else if (a3_est >= 0 & a1_est < 0) {
    row  <- get_row(estimates, "c1_times_2");       combo <- " 2 * c1"
  } else { # a3 < 0 & a1 < 0
    row  <- get_row(estimates, "c2_times_2");       combo <- " 2 * c2"
  }
  est <- as.numeric(row$est); se  <- as.numeric(row$se); p2 <- as.numeric(row$pvalue)
  p1 <- if (est > 0) p2/2 else 1 - p2/2  # one-tailed p for H1: abs > 0
  ci95 <- c(est - 1.96*se, est + 1.96*se)
  
  tibble::tibble(
    Criteria       = paste0("abs (Gap, |c1 − c2| − |c1 + c2|)"),
    Interpretation = "Effect of gap between Self and Other beyond positivity",
    Coefficient    = est, SE = se, p = p1,
    CI_lower       = ci95[1], CI_upper = ci95[2]
  )
}

build_table_short <- function(estimates) {
  c1 <- get_row(estimates, "c1") |> tibble::as_tibble()
  c2 <- get_row(estimates, "c2") |> tibble::as_tibble()
  a1 <- get_row(estimates, "a1") |> tibble::as_tibble()
  abs_row <- compute_abs(estimates)
  
  base <- tibble::tibble(
    Criteria = c("c1 (Self)", "c2 (Other)", "a1 (Positivity, c1+c2)"),
    Interpretation = c(
      "Effect of one's own interest in talking",
      "Effect of beliefs about partner's interest in talking",
      "Effect of positivity (Self and Other together)"
    ),
    Coefficient = c(c1$est, c2$est, a1$est),
    SE          = c(c1$se,  c2$se,  a1$se),
    p           = c(c1$pvalue, c2$pvalue, a1$pvalue),
    CI_lower    = c(c1$ci.lower, c2$ci.lower, a1$ci.lower),
    CI_upper    = c(c1$ci.upper, c2$ci.upper, a1$ci.upper)
  )
  
  short <- dplyr::bind_rows(base, abs_row)
  short
}

short_tbl <- build_table_short(estimates)

# Pretty print
t_short <- nice_table(
  short_tbl,
  title = c("CRA Summary"),
  note  = c(
    "A non-significant abs test suggests that people are not responding to the gap between self and other per se, beyond the positivity 
    of the beliefs (i.e., explicitly comparing how much they are interested in an interaction with a dissimilar partner to how much 
    they think a dissimilar partner is interested in them and using that algebraic difference to guide their willingness to choose a 
    dissimilar conversation partner). Instead, as revealed by the significant c2 pathway, people appear to be guided by the absolute 
    level of both their own interest, and also critically, their metaperceptions: pessimistic beliefs about dissimilar others’ interest 
    in talking predict less willingness to choose dissimilar partners, holding own interest constant.",
    "Coefficients are unstandardized.",
    "Model: Y ~ c0 + c1*Self + c2*Other + e, where Self = own interest, Other = believed partner’s interest.",
    "c1 = partial slope for Self (controlling for slope of Other).",
    "c2 = partial slope for Other (controlling for slope of Self).",
    "a1 = c1 + c2 (effect of positivity).",
    "For c1, c2, and a1, we report two-tailed p values and 95% CIs.",
    "abs = |c1 − c2| − |c1 + c2|, indexes whether the size of the self–other gap predicts the outcome, beyond positivity.",
    "for abs, the CRA gap-magnitude test, we report a one-tailed p for abs > 0, alongside a two-sided 95% CI for the linear contrast."
  )
)
t_short

# Fix the column range
t_short_final <- t_short |>
  width(j = "Criteria", width = 2.25) |>
  width(j = "Interpretation", width = 4.0) |>
  width(j = "b", width = 0.5) |>
  width(j = "SE", width = 0.5) |>
  width(j = "p", width = 0.9) |>
  width(j = "95% CI", width = 1.1)  # CI column might need to be a bit wider

# View it
t_short_final

# Or save as PNG
save_as_image(t_short_final, path = glue("{PLOT_PATH}/cra_summary_table.png"), 
              zoom = 5, expand = 10)

# Save word
#save_as_docx(t_short_final, path = glue("{PLOT_PATH}/cra_summary_table.docx"))

#bigger table-------------------------------------------------

# Compute |D| (distance) row from lm(Y ~ P + D + |D|)
compute_dabs_row <- function(m_mag) {
  co <- summary(m_mag)$coef
  if (!("Dabs" %in% rownames(co))) stop("Dabs not found in m_mag.")
  est <- co["Dabs","Estimate"]; se <- co["Dabs","Std. Error"]
  tval <- co["Dabs","t value"];   p2 <- co["Dabs","Pr(>|t|)"]
  ci95 <- c(est - 1.96*se, est + 1.96*se)
  tibble::tibble(
    Criteria       = "|D| (distance at fixed P; linear V-shape proxy)",
    Interpretation = "Does mismatch magnitude add beyond P and D? (two-tailed)",
    Coefficient = est, SE = se, p = p2,
    CI_lower = ci95[1], CI_upper = ci95[2]
  )
}

# Compute LOIC curvature a4 = b_S^2 - b_SR + b_R^2 from RSA model
compute_loic_curvature <- function(m_poly) {
  b    <- coef(m_poly); V <- vcov(m_poly)
  nms  <- names(b)
  iS2  <- which(nms == "I(Sc^2)")
  iSR  <- which(nms == "I(Sc * Rc)")
  iR2  <- which(nms == "I(Rc^2)")
  L    <- rep(0, length(b)); L[iS2] <- 1; L[iSR] <- -1; L[iR2] <- 1
  est  <- as.numeric(crossprod(L, b))
  se   <- sqrt(as.numeric(t(L) %*% V %*% L))
  z    <- est/se
  p2   <- 2*pnorm(-abs(z))
  ci95 <- c(est - 1.96*se, est + 1.96*se)
  tibble::tibble(
    Criteria       = "a4 (LOIC curvature; RSA)",
    Interpretation = "Distance-from-match (|S−R|) at fixed P; quadratic curvature (two-tailed)",
    Coefficient = est, SE = se, p = p2,
    CI_lower = ci95[1], CI_upper = ci95[2]
  )
}

# --- full table (CRA core + Q2 distance rows) -----------------------------

build_table_full <- function(estimates, m_mag, m_poly) {
  c1 <- get_row(estimates, "c1") |> tibble::as_tibble()
  c2 <- get_row(estimates, "c2") |> tibble::as_tibble()
  a1 <- get_row(estimates, "a1") |> tibble::as_tibble()
  a3 <- get_row(estimates, "a3") |> tibble::as_tibble()
  abs_row  <- compute_abs(estimates)
  dabs_row <- compute_dabs_row(m_mag)
  a4_row   <- compute_loic_curvature(m_poly)
  
  base <- tibble::tibble(
    Criteria = c("c1 (Self, S)",
                 "c2 (Other, R)",
                 "a1 (Positivity; LOC: move S & R together)",
                 "a3 (Directional tilt; LOIC: raise S, lower R equally)"),
    Interpretation = c(
      "Slope for S, controlling R (two-tailed)",
      "Slope for R, controlling S (two-tailed)",
      "Effect of overall positivity (two-tailed)",
      "Directional discrepancy at fixed positivity (two-tailed)"
    ),
    Coefficient = c(c1$est, c2$est, a1$est, a3$est),
    SE          = c(c1$se,  c2$se,  a1$se,  a3$se),
    p           = c(c1$pvalue, c2$pvalue, a1$pvalue, a3$pvalue),
    CI_lower    = c(c1$ci.lower, c2$ci.lower, a1$ci.lower, a3$ci.lower),
    CI_upper    = c(c1$ci.upper, c2$ci.upper, a1$ci.upper, a3$ci.upper)
  )
  
  full <- dplyr::bind_rows(base, abs_row, dabs_row, a4_row)
  full
}

full_tbl <- build_table_full(estimates, m_mag, m_poly)

# If you have your nice_table() util:
t_full <- nice_table(
  full_tbl,
  title = c("CRA/RSA Decomposition of Choice", "Comprehensive results"),
  note  = c(
    "c1 and c2 are slopes for S and R from Y ~ c1*S + c2*R.",
    "a1 = c1 + c2 (positivity; LOC). a3 = c1 − c2 (directional tilt; LOIC).",
    "abs tests Humberg’s path-independent difference rule (H1: abs>0, one-tailed); CI shown is the two-sided 95% CI for the tested linear contrast.",
    "|D| row is a linear V-shape proxy from lm(Y ~ P + D + |D|); a4 is LOIC curvature from the RSA polynomial (two-tailed).",
    "Interpretation: positivity and R-levels reflect overall expectations; a3 indexes directional tilt at fixed positivity; |D|/a4 index a mismatch (distance) effect at fixed positivity."
  )
)
t_full


#testing CRA "on both sides" - is it different for positive and negative! --------------------------------------
#which is actually distinct from Q2 - not necessary if Q2 is true that either side is also true (because might not be more than positivity)

library(lavaan)

# build H and interactions
d_int <- d_cra |>
  dplyr::mutate(H = as.numeric(S_dissim >= R_dissim),
                SxH = S_dissim * H,
                RxH = R_dissim * H)

model_regional <- '
  # --- Core regression with region interactions (H = 1 if S >= R; else 0)
  dissim_choice ~ c1*S_dissim + c2*R_dissim + g1*SxH + g2*RxH

  # --- Region-specific slopes
  c1pos := c1 + g1
  c2pos := c2 + g2
  # Negative-gap region uses c1, c2 (when H=0)

  # --- Region-specific a1/a3 (positivity and directional tilt)
  a1neg := c1 + c2
  a3neg := c1 - c2
  a1pos := c1pos + c2pos
  a3pos := c1pos - c2pos

  # --- Linear combos for Humberg’s abs test in each region
  #     (write in terms of the base coefficients to avoid re-referencing)
  c1_times_2           :=  2*c1
  c2_times_2           :=  2*c2
  minus_c1_times_2     := -2*c1
  minus_c2_times_2     := -2*c2

  c1pos_times_2        :=  2*(c1 + g1)
  c2pos_times_2        :=  2*(c2 + g2)
  minus_c1pos_times_2  := -2*(c1 + g1)
  minus_c2pos_times_2  := -2*(c2 + g2)
'

fit_regional <- lavaan::sem(model_regional, data = d_int, meanstructure = TRUE)
est <- lavaan::parameterEstimates(fit_regional)

get_lab <- function(est, lab) {
  out <- subset(est, label == lab, select = c("est","se","pvalue","ci.lower","ci.upper"))
  if (nrow(out) != 1) stop(paste("Label not found or not unique:", lab))
  out
}

compute_abs_region <- function(est, region = c("neg","pos")) {
  region <- match.arg(region)
  
  # Pull a1, a3 and the needed contrast rows per region
  if (region == "neg") {
    a1 <- get_lab(est, "a1neg"); a3 <- get_lab(est, "a3neg")
    # choose the proper linear combo label based on the signs of a1 and a3
    if (a3$est >= 0 & a1$est > 0)      combo_lab <- "minus_c2_times_2"
    else if (a3$est < 0 & a1$est > 0)  combo_lab <- "minus_c1_times_2"
    else if (a3$est >= 0 & a1$est < 0) combo_lab <- "c1_times_2"
    else                                combo_lab <- "c2_times_2"
  } else {
    a1 <- get_lab(est, "a1pos"); a3 <- get_lab(est, "a3pos")
    if (a3$est >= 0 & a1$est > 0)      combo_lab <- "minus_c2pos_times_2"
    else if (a3$est < 0 & a1$est > 0)  combo_lab <- "minus_c1pos_times_2"
    else if (a3$est >= 0 & a1$est < 0) combo_lab <- "c1pos_times_2"
    else                                combo_lab <- "c2pos_times_2"
  }
  
  combo <- get_lab(est, combo_lab)
  est_lin <- as.numeric(combo$est)
  se_lin  <- as.numeric(combo$se)
  p2      <- as.numeric(combo$pvalue)   # two-tailed p reported by lavaan
  
  # One-tailed p for H1: abs > 0
  p1 <- if (est_lin > 0) p2/2 else 1 - p2/2
  ci95 <- c(est_lin - 1.96*se_lin, est_lin + 1.96*se_lin)
  
  data.frame(
    Region = region,
    a1 = a1$est, a3 = a3$est,
    Combo = combo_lab,  # which linear combination implements abs in this sign case
    abs_est = est_lin,  SE = se_lin,
    p_one_tailed = p1,
    CI_lower = ci95[1], CI_upper = ci95[2],
    stringsAsFactors = FALSE
  )
}

# Build a tidy table with both regions
abs_neg <- compute_abs_region(est, "neg")
abs_pos <- compute_abs_region(est, "pos")
abs_tbl <- rbind(abs_neg, abs_pos)
abs_tbl




#Q1
# Create P and D
d_cra <- d_cra |>
  dplyr::mutate(P = 0.5*(S_dissim + R_dissim),
                D = 0.5*(S_dissim - R_dissim))

m_PD <- lm(dissim_choice ~ P + D, data = d_cra)
summary(m_PD)   # the coefficient on D equals a3


#Q2

#Minimal proxy (linear + absolute value):
d_cra <- d_cra |>
  dplyr::mutate(Dabs = abs(0.5*(S_dissim - R_dissim)))
m_mag <- lm(dissim_choice ~ P + D + Dabs, data = d_cra)
summary(m_mag)  # the coefficient on Dabs tests |gap| while controlling positivity and direction

#Principled RSA (second order):
# Center S and R before squaring/interacting
d_center <- d_cra |>
  dplyr::mutate(Sc = scale(S_dissim, center = TRUE, scale = FALSE)[,1],
                Rc = scale(R_dissim, center = TRUE, scale = FALSE)[,1])

m_poly <- lm(dissim_choice ~ Sc + Rc + I(Sc^2) + I(Sc*Rc) + I(Rc^2), data = d_center)

# Derived LOIC curvature a4 = b3 - b4 + b5
b <- coef(m_poly)
a4 <- b["I(Sc^2)"] - b["I(Sc * Rc)"] + b["I(Rc^2)"]
# For a CI/p-value for a4, use car::linearHypothesis or fit the same in lavaan with 'a4 :='

# Coefficients are in order: Intercept, Sc, Rc, I(Sc^2), I(Sc*Rc), I(Rc^2)
# So positions 4, 5, 6 for the quadratic terms
# a4 = b4 - b5 + b6

a4_result <- deltaMethod(m_poly, 
                         "b4 - b5 + b6",
                         parameterNames = paste0("b", 1:6))

print(a4_result)

# Calculate p-value
a4_pvalue <- 2 * (1 - pnorm(abs(a4_result$Estimate / a4_result$SE)))

cat(sprintf("a4 = %.4f, SE = %.4f, 95%% CI [%.4f, %.4f], p = %.4f\n",
            a4_result$Estimate, 
            a4_result$SE, 
            a4_result$`2.5 %`, 
            a4_result$`97.5 %`, 
            a4_pvalue))

# Test hypothesis: coefficient 4 - coefficient 5 + coefficient 6 = 0
# Create hypothesis matrix [0, 0, 0, 1, -1, 1]
hyp_test <- linearHypothesis(m_poly, 
                             hypothesis.matrix = c(0, 0, 0, 1, -1, 1))
print(hyp_test)  # This gives you the F-test and p-value





plotmodel <- function(estimates, 
                      filename="result_plot",
                      type="3d", model="full", main="",
                      xlim=NULL, ylim=NULL, zlim=NULL, 
                      xlab="Self-view", ylab="Criterion", zlab="Narcissism",
                      axes=c("LOC","LOIC"), 
                      project=c("LOC","LOIC"), 
                      param=F,
                      legend=FALSE){
  
  # extract parameter estimates of c1 and c2
  c1 <- estimates[estimates$label=="c1","est"]
  c2 <- estimates[estimates$label=="c2","est"]
  
  # build plot of specified model
  plot <- plotRSA(x=c1, y=c2, x2=0, xy=0, y2=0,
                  type=type, main = main, 
                  xlim=xlim, ylim=ylim, zlim=zlim, 
                  xlab=xlab, ylab=ylab, zlab=zlab,
                  cex.main = 1.3, 
                  axes=axes, 
                  project=project, 
                  param=param,
                  legend=legend
  )
  
  # show plot
  print(plot)
  
  # prepare saving of plots 
  dir.create("Result_plots", showWarnings = FALSE)
  filename <- paste0("Result_plots/plot_",filename,".png")
  
  # save plot
  png(filename=filename,units="in",width=5,height=5,res=300,pointsize=18)
  print(plot)
  dev.off()
}


### RSA -------------------------------------------------------------

d_rsa <-  d %>%
  rename(perceptions_of_dissimilar_others_interest = dissim_other,
         self_interest_in_dissimilar_others = dissim_self, 
         choice_same = sim_choice,
         choice_different = dissim_choice) %>%
  mutate(interest_interaction = perceptions_of_dissimilar_others_interest * self_interest_in_dissimilar_others,
         others_interest2 = perceptions_of_dissimilar_others_interest * perceptions_of_dissimilar_others_interest,
         self_interest2 = self_interest_in_dissimilar_others * self_interest_in_dissimilar_others,
         dissim_sim_choice_delta = choice_different-choice_same)

# significance test of the full polynomial model before proceeding
d_rsa %>%
  lm(choice_different~ perceptions_of_dissimilar_others_interest + self_interest_in_dissimilar_others + interest_interaction + others_interest2 + self_interest2, .) %>%
  summary()

grand_mean <- mean(c(d_rsa$perceptions_of_dissimilar_others_interest, d_rsa$self_interest_in_dissimilar_others), na.rm = TRUE)
grand_mean

d_rsa <-d_rsa %>%
  mutate(perceptions_of_dissimilar_others_interest_centered = perceptions_of_dissimilar_others_interest - grand_mean,
         self_interest_in_dissimilar_others_centered = self_interest_in_dissimilar_others - grand_mean)

rsa.model_same <-  RSA(formula = choice_different ~ perceptions_of_dissimilar_others_interest_centered*self_interest_in_dissimilar_others_centered, 
                       data    = d_rsa,   
                       scale   = FALSE,           
                       na.rm   = TRUE,       
                       out.rm  = TRUE,       
                       models  = c("onlyx","onlyx2","onlyy","onlyy2","additive", "absdiff",
                                   "diff","SQD","SSQD","RR","SRR","SRRR", "IA","null","full"), 
                       missing = "listwise") 
aictab(rsa.model_same)

rsa.model_same.reduced <-  RSA(formula = choice_different ~ perceptions_of_dissimilar_others_interest_centered*self_interest_in_dissimilar_others_centered, 
                               data    = d_rsa,     
                               scale   = FALSE,           
                               na.rm   = TRUE,       
                               out.rm  = TRUE,       
                               models  = c("SRR","additive", "null"), 
                               missing = "listwise") 
aictab(rsa.model_same.reduced)
summary(rsa.model_same.reduced, model = "additive" )

plot(rsa.model_same.reduced, model = "SRR", legend=F)
plot(rsa.model_same.reduced, model = "additive", legend=F)


plot(rsa.model_same.reduced, model = "additive", legend=F, bw=TRUE, points=TRUE, border = TRUE, project = c("contour"), hull = FALSE,
     xlab = "Perceptions of Dissimilar Other's Interest",  ylab = "Interest in Dissimilar Other", zlab = "Willingness to Talk to Dissimilar Person")


plot(rsa.model_same.reduced, model = "additive", legend=F, bw=FALSE, points=TRUE, border = TRUE, project = c("contour"), hull = FALSE,
     xlab = "Perceptions of Dissimilar Other's Interest",  ylab = "Interest in Dissimilar Other", zlab = "Willingness to Talk to Dissimilar Person")

getPar(rsa.model_same.reduced, model = "additive", type = "coef")

plotRSA(rsa.model_same.reduced, model = "additive")


d_rsa %>%
  lm(choice_different ~ perceptions_of_dissimilar_others_interest + self_interest_in_dissimilar_others + interest_interaction + others_interest2 + self_interest2, .) %>%
  summary()

plotRSA(x  = 0.30634,  # Enter main effect of predictor 1 from mlm model
        y  = 0.77656,  # Enter main effect of predictor 2 from mlm model
        x2 = -0.06509, # Enter squared effect of predictor 1 from mlm model
        y2 = -0.02721 , # Enter squared effect of predictor 2 from mlm model
        xy = 0.08160, # Enter interaction effect from mlm model
        b0 = -0.37778, # Enter intercept from mlm model
        type = "3d", surface = "predict", model = "additive", bw=TRUE, points = TRUE)

