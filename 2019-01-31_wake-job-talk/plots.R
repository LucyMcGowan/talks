library(tidyverse)
library(tipr)
library(survival)
library(survey)
groups <- list(
  "Accessories" = 
    c("paco21", "ph1", "meanbp1"))

rhc <- read_csv("http://biostat.mc.vanderbilt.edu/wiki/pub/Main/DataSets/rhc.csv")
rhc <- rhc %>%
  mutate(
    exposure = case_when(
      swang1 == "RHC" ~ 1,
      TRUE ~ 0
    ),
    time = t3d30,
    event = case_when(
      dth30 == "Yes" ~ 1,
      TRUE ~ 0
    )
  )

vars <- c("meanbp1", "paco21", "dnr1", "ph1", "resp1",
          "aps1", "surv2md1", "resp1", "pafi1")

ps_frm <- as.formula(
  paste("exposure ~ ", paste(vars, collapse = "+"))
)

ps_mod <- glm(ps_frm,
              family = binomial(link = "logit"),
              data = rhc)

rhc <- rhc %>%
  mutate(
    p = predict(ps_mod, type = "response"),
    weight = case_when(
      exposure == 1 ~ 1 - p,
      TRUE ~ p
    )
  )

weight_des <- svydesign(ids = ~ 1,
                        data = rhc,
                        weights = ~ weight)

run_outcome_mod <- function(formula, wt) {
  des <- svydesign(ids = ~ 1,
                   data = rhc, 
                   weights = wt)
  svycoxph(formula = formula,
           design = des)
}

outcome_frm <- as.formula(
  paste("Surv(time, event) ~ exposure + ",
        paste(vars, collapse = "+"))
)

outcome_mod <- svycoxph(formula = outcome_frm,
                        design = weight_des)


full_model <- data.frame(
  point_estimate = exp(coef(outcome_mod)["exposure"]),
  lb = exp(confint(outcome_mod, "exposure")[, 1]),
  ub = exp(confint(outcome_mod, "exposure")[, 2])
)

o <- observed_bias_tbl(ps_mod, outcome_mod, groups = groups)

o <- o %>%
  mutate(
    weight = 
      map(p, ~ case_when(rhc$exposure == 1 ~ 1 - .x, TRUE ~ .x))
  )

o <- o %>%
  mutate(
    outcome_model =
      map2(outcome_formula, weight, run_outcome_mod)
  )

o <- o %>%
  mutate(
    point_estimate = 
      map_dbl(outcome_model, ~ exp(coef(.x)["exposure"])),
    lb = 
      map_dbl(outcome_model, ~ exp(confint(.x, "exposure")[, 1])),
    ub = 
      map_dbl(outcome_model, ~ exp(confint(.x, "exposure")[, 2]))
  )

o <- bind_rows(
  o,
  observed_bias_tip(
    tip = full_model$lb,
    full_model$point_estimate,
    full_model$lb,
    full_model$ub,
    "Hypothetical unmeasured confounder (Tip LB)"),
  observed_bias_tip(
    tip = full_model$point_estimate,
    full_model$point_estimate,
    full_model$lb,
    full_model$ub,
    "Hypothetical unmeasured confounder (Tip Point Est)")
)

o <- o %>%
  mutate(
    e_value = map2_dbl(
      lb,
      ub,
      ~ observed_covariate_e_value(
        lb = full_model$lb,
        ub = full_model$ub,
        lb_adj = .x,
        ub_adj = .y,
        transform = "HR")
    )
  )

o <- o %>%
  left_join(
    read_csv("var_labels.csv"),
    by = c("dropped" = "var")
  ) %>%
  mutate(
    dropped = case_when(
      type == "covariate" ~ var_label,
      TRUE ~ dropped
    )
  )

o <- observed_bias_order(o, "lb")
o_ <- o[-c(1:2), ]

## no e-value ----
g <- ggplot(data = o_)
g <- g + 
  geom_hline(yintercept = full_model$point_estimate,
             lwd = 2, 
             color = "light blue") + 
  geom_hline(yintercept = 1, lty = 2) +
  geom_rect(ymin = full_model$lb,
            ymax = full_model$ub,
            xmin = 0,
            xmax = 25,
            alpha = 0.01,
            fill = "light blue")

g <- g + 
  geom_pointrange(aes(x = dropped,
                      y = point_estimate,
                      ymin = lb,
                      ymax = ub)) + 
  coord_flip()


g <- g + 
  ylab("Effect modeled without variable(s) of interest") +
  xlab("") +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = "bottom")

g
## yes e-value ----
g <- ggplot(data = o)
g <- g + 
  geom_hline(yintercept = full_model$point_estimate,
             lwd = 2, 
             color = "light blue") + 
  geom_hline(yintercept = 1, lty = 2) +
  geom_rect(ymin = full_model$lb,
            ymax = full_model$ub,
            xmin = 0,
            xmax = 25,
            alpha = 0.01,
            fill = "light blue")

g <- g + 
  geom_pointrange(aes(x = dropped,
                      y = point_estimate,
                      ymin = lb,
                      ymax = ub)) + 
  coord_flip()

g <- g + 
  geom_point(aes(x = dropped, y = e_value, color = type),
             pch = 8) + 
  scale_color_manual(name = "",
                     values = c("covariate" = "purple",
                                "group" = "orange",
                                "tip" = "red"),
                     labels = c("Observed Covariate E-value",
                                "Observed Covariate E-value (group)",
                                "E-value")
  )

g <- g + 
  ylab("Effect modeled without variable(s) of interest") +
  xlab("") +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = "bottom")

g
