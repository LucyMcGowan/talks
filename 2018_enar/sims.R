#cor = cov/sqrt(v1*v2)
library(tidyverse)
library(sandwich)
library(lmtest)

dr <- function(weight, y, m_1, m_0, z) {
  (sum(weight * (m_1 - m_0)) / sum(weight)) + 
    (sum(weight * z * (y - m_1)) / sum(weight * z)) -
    (sum(weight * (1 - z) * (y - m_0)) / sum(weight * (1 - z)))
}

ps <- function(weight, y, z) {
  sum(weight * z * y) / sum(weight * z) - sum(weight * (1 - z) * y) / sum(weight * (1 - z))
}

sim <- function(form_ps = "z ~ x_1",
                form_out = "y ~ z + x_1",
                ps_family = binomial("probit"),
                outcome_family = "gaussian",
                n = 1000,
                mu_x_1 = 0.5,
                mu_x_2 = 1,
                cov_z = 1,
                var_x_1 = 2,
                var_x_2 = 1,
                a = 1, 
                b = 1,
                c_1 = 1,
                c_2 = 2,
                d = 1,
                e = 0.5,
                f_1 = 0.25,
                f_2 = 0.75) {
  X <- mvtnorm::rmvnorm(n,
                        mean = c(mu_x_1, mu_x_2),
                        sigma = matrix(c(var_x_1, cov_z, cov_z, var_x_2), ncol = 2)
  )
  x_1 <- X[, 1]
  x_2 <- X[, 2]
  V <- rnorm(n, 0, 1)
  z <- as.numeric(e + f_1 * x_1 + f_2 * x_2 + V > 0)
  if (outcome_family == "gaussian") {
    U <- rnorm(n, 0, 1)
    y <- a + b * z + c_1 * x_1 + c_2 * x_2 + d * U
  } else {
    U <- rlogis(n, 0, 1)
    y_lin <- a + b * z + c_1 * x_1 + c_2 * x_2 + d * U
    y <- (y_lin > 0)
  }
  
  dat <- data.frame(
    y = y,
    z = z,
    x_1 = x_1,
    x_2 = x_2
  )
  
  form_out2 <- gsub("z", "", form_out)
  dat <- dat %>%
    mutate(
      p = predict(glm(as.formula(form_ps), data = dat, family = ps_family), type = "response"),
      m1 = predict(
        glm(as.formula(form_out2), data = dat[dat$z == 1, ], family = outcome_family),
        newdata = dat
      ),
      m0 = predict(
        glm(as.formula(form_out2), data = dat[dat$z == 0, ], family = outcome_family),
        newdata = dat
      ),
      p_0 = 1 - p,
      p_assign = case_when(
        z == 1 ~ p,
        z == 0 ~ p_0
      ),
      p_min = pmin(p, p_0),
      atm = p_min / p_assign,
      ato = 1 - p_assign,
      ate = 1 / p_assign
    )
  
  ato <- dr(dat$ato, dat$y, dat$m1, dat$m0, dat$z)
  atm <- dr(dat$atm, dat$y, dat$m1, dat$m0, dat$z)
  ate <- dr(dat$ate, dat$y, dat$m1, dat$m0, dat$z)
  
  ato_ps <- ps(dat$ato, dat$y, dat$z)
  atm_ps <- ps(dat$atm, dat$y, dat$z)
  ate_ps <- ps(dat$ate, dat$y, dat$z)

  naive <- glm(as.formula(form_out), data = dat, family = outcome_family)
  
  t <- tibble(
    method = c("pscore_only_ato", "pscore_only_atm", "pscore_only_ate", "outcome_and_pscore_ato", "outcome_and_pscore_atm", "outcome_and_pscore_ate","outcome_model_only"),
    est = c(
      ato_ps, atm_ps, ate_ps,
      ato, atm, ate,
      naive$coefficients[["z"]])
  )
  t
}

set.seed(924)

corr_0_200 <- purrr::map_df(1:1000, ~ sim(form_ps = "z ~ x_1",
                                          form_out = "y ~ z + x_1", cov_z = 0, n = 200)
)
corr_0_200 %>%
  group_by(method) %>%
  summarise(mean(est))

corr_0_1000 <- purrr::map_df(1:1000, ~ sim(form_ps = "z ~ x_1",
                                           form_out = "y ~ z + x_1", cov_z = 0)
)
corr_0_1000 %>%
  group_by(method) %>%
  summarise(mean(est))
corr_0_5000 <- purrr::map_df(1:1000, ~ sim(form_ps = "z ~ x_1",
                                           form_out = "y ~ z + x_1", cov_z = 0, n = 5000)
)
corr_0_5000 %>%
  group_by(method) %>%
  summarise(mean(est))
# corr_neg85 <- purrr::map_df(1:1000, ~sim(form_ps = "z ~ x_1",
#     form_out = "y ~ z + x_1", cov_z = -1.2)
# )
# 
# corr_neg85 %>%
#   group_by(method) %>%
#   summarise(mean(est))

corr_85_200 <- purrr::map_df(1:1000, ~sim(form_ps = "z ~ x_1",
                                          form_out = "y ~ z + x_1", cov_z = 1.2, n = 200)
)

corr_85_200 %>%
  group_by(method) %>%
  summarise(mean(est))


corr_85_1000 <- purrr::map_df(1:1000, ~sim(form_ps = "z ~ x_1",
                                           form_out = "y ~ z + x_1", cov_z = 1.2)
)

corr_85_1000 %>%
  group_by(method) %>%
  summarise(mean(est))

corr_85_5000 <- purrr::map_df(1:1000, ~sim(form_ps = "z ~ x_1",
                                           form_out = "y ~ z + x_1", cov_z = 1.2, n = 5000)
)

corr_85_5000 %>%
  group_by(method) %>%
  summarise(mean(est))

corr_1_200 <- purrr::map_df(1:1000, ~sim(form_ps = "z ~ x_1",
                                         form_out = "y ~ z + x_1", cov_z = 1.38, n = 200)
)
corr_1_200 %>%
  group_by(method) %>%
  summarise(mean(est), sd(est)/mean(sd), mean(sd))
corr_1_1000 <- purrr::map_df(1:1000, ~sim(form_ps = "z ~ x_1",
                                          form_out = "y ~ z + x_1", cov_z = 1.38)
)
corr_1_1000 %>%
  group_by(method) %>%
  summarise(mean(est), sd(est)/mean(sd), mean(sd))

corr_1_5000 <- purrr::map_df(1:1000, ~sim(form_ps = "z ~ x_1",
                                          form_out = "y ~ z + x_1", cov_z = 1.38, n = 5000)
)
corr_1_5000 %>%
  group_by(method) %>%
  summarise(mean(est), sd(est)/mean(sd), mean(sd))

d <- bind_rows("corr_0_200" = corr_0_200, 
               "corr_0_1000" = corr_0_1000,
               "corr_0_5000" = corr_0_5000, 
               "corr_85_200" = corr_85_200,
               "corr_85_1000" = corr_85_1000, 
               "corr_85_5000" = corr_85_5000,
               "corr_1_200" = corr_1_200, 
               "corr_1_1000" = corr_1_1000, 
               "corr_1_5000" = corr_1_5000, .id = "corr")

d %>%
  group_by(corr, method) %>%
  summarise(bias = mean(est) - 1) %>%
  mutate(
    n = case_when(
      grepl("200", corr) ~ 200,
      grepl("1000", corr) ~ 1000,
      grepl("5000", corr) ~ 5000
    ),
    model = case_when(
      method == "outcome_and_pscore" ~ "propensity score + outcome",
      method == "outcome_model_only" ~ "outcome only"
    ),
    correlation = case_when(
      grepl("corr_0", corr) ~ "~0",
      grepl("corr_85", corr) ~ "~0.85",
      grepl("corr_1", corr) ~ "~1"
    ),
    group = case_when(
      correlation == "~0" & model == "propensity score + outcome" ~ 1,
      correlation == "~0.85" & model == "propensity score + outcome" ~ 2,
      correlation == "~1" & model == "propensity score + outcome" ~ 3, 
      correlation == "~0" & model == "outcome only" ~ 4,
      correlation == "~0.85" & model == "outcome only" ~ 5,
      correlation == "~1" & model == "outcome only" ~ 6
    )) %>%
  ggplot(aes(x = n, y = bias, group = group, color = correlation, linetype = model)) +
  geom_point() + 
  geom_line() +     
  scale_x_continuous(breaks = c(200, 1000, 5000)) 

# d <- purrr::map_df(1:1000, ~sim(form_ps = "z ~ x_1",
#                                 form_out = "y ~ z + x_1", cov_z = 1.35, outcome_family = "binomial")
# )
# 
# d %>%
#   group_by(method) %>%
#   summarise(mean(est), sd(est)/mean(sd))
# 
# d <- purrr::map_df(1:100, ~sim(form_ps = "z ~ x_1",
#                                 form_out = "y ~ z + x_1", cov_z = 0, outcome_family = "binomial")
# )
# 
# d %>%
#   group_by(method) %>%
#   summarise(mean(est), sd(est)/mean(sd))
