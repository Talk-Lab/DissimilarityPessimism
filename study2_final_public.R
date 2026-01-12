# Title: Dissimilarity Pessimism Study 2
# Author(s): Gus Cooney; Erica Boothby
# Description: Dissimilarity pessimism & conversations across an intergroup divide

# load libraries ----------------------------------------------------------

library(pacman)
p_load(tidyverse, ggplot2, glue, corrr, broom, magrittr, ggExtra, lme4, ggeffects, lmerTest, 
       emmeans, pbkrtest, questionr, knitr, sjPlot, effects, report, brglm, scales, optimx)

# set constants -----------------------------------------------------------

DATA_PATH <- "~/Dropbox/Manuscripts/Homophily/Data"
PLOT_PATH <- "~/Dropbox/Manuscripts/Homophily/Plots"

# load data ---------------------------------------------------------------

raw <- read_csv(glue("{DATA_PATH}/study2_final.csv"))
d <- raw
d

# exclusions ---------------------------------------------------------------

# exclusions
subset <- subset(d, check5 == 0)
d <- subset

#by matched race
subset <- subset(d, race == 1 | race == 5)
d <- subset

#by matched gender
subset <- subset(d, sex == 1 | sex == 2)
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

# manipulation check ------------------------------------------------------------

#manipulation
d.manipulation <- d[,c("sim_similar", "dissim_similar")]
d.manipulation <- d.manipulation %>%
  pivot_longer(cols = c(sim_similar, dissim_similar),
               names_to = "similarity",
               values_to = "rating")

t.test(d$sim_similar, d$dissim_similar, paired = TRUE)

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

# format data -------------------------------------------------------------------

test <- d %>%
  pivot_longer(
    cols = c(sim_self, sim_other, dissim_self, dissim_other), 
    names_to = c("sim_diff", "self_other"),
    names_sep = "_",
    values_to = "rating"
  )

test

d_long <- d %>%
  pivot_longer(
    cols = c(sim_self, sim_other, dissim_self, dissim_other),
    names_to   = c("sim_dissim", "self_other"),
    names_sep  = "_",
    values_to  = "rating"
  ) %>%
  mutate(
    pid = participant_id
  )

d_long <- d_long %>%
  mutate(
    # within-person effect coding (as in Study 1)
    rating_type_c  = ifelse(self_other  == "self", +0.5, -0.5),
    partner_type_c = ifelse(sim_dissim  == "sim",  +0.5, -0.5),
    
    # between-person race factor (assumes d$black_white is "white"/"black")
    race_f = factor(black_white, levels = c("white","black")),
    
    # effect-coded race (±0.5); choose sign convention and be consistent in the text
    race_c = ifelse(race_f == "white", +0.5, -0.5)
  )

d_long

write.csv(d_long, file = glue("{DATA_PATH}/z_study2_long_aggregate.csv"))

# modeling -----------------------------------------------------------------------

fit_s2 <- lmer(
  rating ~ rating_type_c * partner_type_c * race_c +
    (1 + rating_type_c + partner_type_c || pid),
  data    = d_long,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)),
  REML    = FALSE
)

summary(fit_s2)
confint(fit_s2)

model <- fit_s2

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
#note: the extra factor compared to Study 1 (race, a between-person fixed effect) does not change the denominator formula, given the final random-effects structure.

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

emm_H1 <- emmeans(model, ~ partner_type_c, 
                       at = list(partner_type_c = c(-0.5, +0.5)))
con_H1 <- contrast(emm_H1, list("similar - dissimilar" = c(-1, 1)))
con_d_H1   <- contrast_to_d(con_H1, sigma_total)

emm_H1

pretty_emm_H1 <- summary(emm_H1) %>%
  dplyr::mutate(
    partner = dplyr::recode(as.character(partner_type_c),
                            `-0.5` = "Dissimilar Partner",
                            `0.5`  = "Similar Partner"),
    `95% CI` = sprintf("%.2f [%.2f, %.2f]", emmean, lower.CL, upper.CL)
  ) %>%
  dplyr::select(partner, emmean, SE, df, `95% CI`) %>%
  dplyr::arrange(partner)

pretty_emm_H1
con_d_H1 

#by race
emm_H1_race <- emmeans(model, ~ partner_type_c | race_c,
                       at = list(partner_type_c = c(-0.5, +0.5)))
con_H1_race <- contrast(emm_H1_race, list("similar - dissimilar" = c(-1, 1)), by = "race_c")
con_d_H1_race   <- contrast_to_d(con_H1_race, sigma_total)

emm_H1_race
con_d_H1_race

pretty_emm_H1_race <- summary(emm_H1_race) %>%
  dplyr::mutate(
    race = dplyr::recode(as.character(race_c),
                         `-0.5` = "Black participants",
                         `0.5`  = "White participants"),
    partner = dplyr::recode(as.character(partner_type_c),
                            `-0.5` = "Dissimilar Partner",
                            `0.5`  = "Similar Partner"),
    `95% CI` = sprintf("%.2f [%.2f, %.2f]", emmean, lower.CL, upper.CL)
  ) %>%
  dplyr::select(race, partner, emmean, SE, df, `95% CI`) %>%
  dplyr::arrange(race, partner)

# Pretty contrast (similar - dissimilar) by race, with raw and d
pretty_con_d_H1_race <- con_d_H1_race %>%
  dplyr::mutate(
    race = dplyr::recode(as.character(race_c),
                         `-0.5` = "Black participants",
                         `0.5`  = "White participants")
  ) %>%
  dplyr::transmute(
    race, contrast,
    estimate, SE, df, t.ratio,
    `95% CI (raw)` = sprintf("%.2f [%.2f, %.2f]", estimate, lower.CL, upper.CL),
    p.value,
    d, d.SE,
    `95% CI (d)`   = sprintf("%.2f [%.2f, %.2f]", d, d.lower, d.upper)
  ) %>%
  dplyr::arrange(race)

pretty_emm_H1_race
pretty_con_d_H1_race

# -----------------------------------------------
# H2: Self vs. other *within dissimilar* (by race)
# -----------------------------------------------

emm_H2 <- emmeans(model, ~ rating_type_c | partner_type_c * race_c,
                         at = list(rating_type_c = c(-0.5, +0.5),
                                   partner_type_c = c(-0.5, +0.5)))
con_H2 <- contrast(emm_H2, list("self - other" = c(-1, 1)), by = c("partner_type_c","race_c"))
con_d_H2   <- contrast_to_d(con_H2, sigma_total)

emm_H2
con_d_H2


pretty_emm_H2 <- summary(emm_H2) %>%
  dplyr::mutate(
    # label the coded factors inline (no external maps)
    race    = dplyr::recode(as.character(race_c),
                            `-0.5` = "Black",
                            `0.5`  = "White"),
    partner = dplyr::recode(as.character(partner_type_c),
                            `-0.5` = "Dissimilar Partner",
                            `0.5`  = "Similar Partner"),
    rating  = dplyr::recode(as.character(rating_type_c),
                            `-0.5` = "Other",
                            `0.5`  = "Self"),
    `95% CI` = sprintf("%.2f [%.2f, %.2f]", emmean, lower.CL, upper.CL)
  ) %>%
  dplyr::select(race, partner, rating, emmean, SE, df, `95% CI`) %>%
  dplyr::arrange(race, partner, rating)

pretty_con_d_H2 <- con_d_H2 %>%
  dplyr::mutate(
    race    = dplyr::recode(as.character(race_c),
                            `-0.5` = "Black",
                            `0.5`  = "White"),
    partner = dplyr::recode(as.character(partner_type_c),
                            `-0.5` = "Dissimilar Partner",
                            `0.5`  = "Similar Partner")
  ) %>%
  dplyr::transmute(
    race, partner, contrast,        # "self - other"
    estimate, SE, df, t.ratio,
    `95% CI (raw)` = sprintf("%.2f [%.2f, %.2f]", estimate, lower.CL, upper.CL),
    p.value,
    d, d.SE,
    `95% CI (d)` = sprintf("%.2f [%.2f, %.2f]", d, d.lower, d.upper)
  ) %>%
  dplyr::arrange(race, partner)

pretty_emm_H2
pretty_con_d_H2

#actual gap

d_acc <- d_long %>%
  filter(partner_type_c == -0.5) %>%  # cross-race condition
  mutate(
    measure_type = ifelse(self_other == "other", "perceived_other",
                          ifelse(self_other == "self", "actual_self", NA))
  )

fit_accuracy <- lmer(
  rating ~ measure_type * race_c +
    (1 | pid),
  data = d_acc,
  REML = FALSE
)
summary(fit_accuracy)

#testing normal effects
emm_acc <- emmeans(fit_accuracy, ~ measure_type * race_c)
emm_acc

#the correct perceived versus actual contrasts
emm_acc@grid

summary(emm_acc)$emmean

# Contrast 1: White participants' perceived-other (row 4)
#              vs Black participants' actual self (row 1)
cross1 <- contrast(
  emm_acc,
  method = list("White belief − Black actual" = c(-1, 0, 0, 1)),
  adjust = "none"
)

# Contrast 2: Black participants' perceived-other (row 2)
#              vs White participants' actual self (row 3)
cross2 <- contrast(
  emm_acc,
  method = list("Black belief − White actual" = c(0, 1, -1, 0)),
  adjust = "none"
)

summary(cross1)
summary(cross2)

contrast_to_d(cross1, sigma_total)
contrast_to_d(cross2, sigma_total)


#interaction test

# All 2×2 cells in the dissimilar condition only
emm_H2_race <- emmeans(
  model,
  ~ race_c * rating_type_c,
  at = list(
    partner_type_c = -0.5,          # dissimilar only
    rating_type_c  = c(-0.5, 0.5)
  ),
  pbkrtest.limit = 7856
)

summary(emm_H2_race)

## Check the order of rows in summary(emm_H2_race)$emmean.
## With the usual ordering, you should see something like:
## race_c  rating_type_c
## -0.5    -0.5   (Black, Other)
##  0.5    -0.5   (White, Other)
## -0.5     0.5   (Black, Self)
##  0.5     0.5   (White, Self)

# Now define self–other contrasts within each race:
H2_race <- contrast(
  emm_H2_race,
  list(
    "Black_H2" = c(-1,  0,  1,  0),   # (Self - Other) for Black, dissimilar
    "White_H2" = c( 0,  -1, 0,  1)    # (Self - Other) for White, dissimilar
  )
)

summary(H2_race, infer = c(TRUE, TRUE))

# Race difference in the H2 effect: White_H2 - Black_H2
H2_race_diff <- contrast(
  H2_race,
  list("White_minus_Black_H2" = c(-1, 1))
)

summary(H2_race_diff, infer = c(TRUE, TRUE))


# -----------------------------------------------
# H3:  Interaction (difference-in-differences)
# -----------------------------------------------
con_H3_race <- contrast(con_H2 , method = "pairwise", by = "race_c")  # (dissim - similar) within race
con_d_H3_race   <- contrast_to_d(con_H3_race, sigma_total)

pretty_con_d_H3_race <- con_d_H3_race %>%
  dplyr::mutate(
    race = dplyr::recode(as.character(race_c),
                         `-0.5` = "Black",
                         `0.5`  = "White"),
    contrast_label = "(self–other)_dissimilar – (self–other)_similar"
  ) %>%
  dplyr::transmute(
    race, contrast = contrast_label,
    estimate, SE, df, t.ratio,
    `95% CI (raw)` = sprintf("%.2f [%.2f, %.2f]", estimate, lower.CL, upper.CL),
    p.value,
    d, d.SE,
    `95% CI (d)` = sprintf("%.2f [%.2f, %.2f]", d, d.lower, d.upper)
  ) %>%
  dplyr::arrange(race)

pretty_con_d_H3_race

#interaction

# full 2×2×2 grid as you already had
emm_full <- emmeans(
  model,
  ~ race_c * partner_type_c * rating_type_c,
  at = list(
    rating_type_c  = c(-0.5, 0.5),
    partner_type_c = c(-0.5, 0.5)
  ),
  pbkrtest.limit = 7856
)

# H3 for Black and White
L_Black_H3 <- c(-1, 0,  1, 0,  1, 0, -1, 0)
L_White_H3 <- c( 0,-1,  0, 1,  0, 1,  0,-1)

H3_race <- contrast(
  emm_full,
  list(
    "Black_H3" = L_Black_H3,
    "White_H3" = L_White_H3
  )
)

summary(H3_race, infer = c(TRUE, TRUE))

# Race difference in H3: White_H3 – Black_H3
H3_race_diff <- contrast(
  H3_race,
  list("White_minus_Black_H3" = c(-1, 1))
)

summary(H3_race_diff, infer = c(TRUE, TRUE))




# plot --------------------------------------------------------------------------

## 1) Build plotting data from emmeans — dissimilar (cross-race) only
black_lab <- "BLACK participant thinking they will\n talk to a WHITE participant"
white_lab <- "WHITE participant thinking they will\n talk to a BLACK participant"

plot_means_dissim <- summary(emm_H2) %>%
  dplyr::filter(partner_type_c == -0.5) %>%                 # dissimilar/cross-race only
  dplyr::mutate(
    race = dplyr::recode(as.character(race_c),
                         `-0.5` = "Black", `0.5` = "White"),
    x_label = dplyr::recode(race, Black = black_lab, White = white_lab),
    rating  = dplyr::recode(as.character(rating_type_c),
                            `-0.5` = "Other", `0.5` = "Self")
  ) %>%
  dplyr::select(x_label, rating, emmean, lower.CL, upper.CL)

# The two x positions in the order you want
plot_df <- plot_means_dissim %>%
  dplyr::mutate(
    x_label  = factor(x_label, levels = c(black_lab, white_lab)),
    self_other = factor(rating, levels = c("Self","Other"))
  ) %>%
  dplyr::rename(sim_dissim = x_label)

## 2) Compute star positions (no brackets), using model CIs to place them
tops <- plot_df %>%
  group_by(sim_dissim) %>%
  summarise(y_top = max(upper.CL), .groups = "drop")

# derive stars from your H2 p-values for dissimilar/cross-race
stars_tbl <- pretty_con_d_H2 %>%
  dplyr::filter(partner == "Dissimilar Partner") %>%
  dplyr::mutate(
    sim_dissim = dplyr::recode(race,
                               "Black" = black_lab,
                               "White" = white_lab),
    stars = dplyr::case_when(
      p.value < .001 ~ "***",
      p.value < .01  ~ "**",
      p.value < .05  ~ "*",
      TRUE           ~ "ns"
    )
  ) %>%
  dplyr::select(sim_dissim, stars)

stars_pos <- tops %>%
  left_join(stars_tbl, by = "sim_dissim") %>%
  mutate(y = y_top + 0.12)  # small gap above the taller CI in each group

## 3) Plot (same ggpubr look as Study 1), add stars (no brackets)
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
                  labels = c("Self: Own Interest",
                             "Other: Estimate of Partner's Interest")) +
  labs(y = "Interest in Having a Conversation", x = NULL) +
  theme(
    text = element_text(size = 16),
    legend.title = element_blank(),
    legend.text  = element_text(size = 16),
    legend.key   = element_blank(),
    axis.text.x  = element_text(size = 14, vjust = 1)  # long labels
  ) +
  coord_cartesian(ylim = c(3, 6.5)) +
  # stars only (no brackets)
  geom_text(
    data = stars_pos,
    aes(x = sim_dissim, y = y, label = stars),
    size = 6, vjust = 0 #fontface = "bold"
  )

fig2 <- bp

fig2


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

save_plot_multi(fig2,
                file_stem = "fig2",
                exts      = c("png", "pdf", "jpeg"),   # add/remove as needed
                width     = 10,
                height    = 8)


