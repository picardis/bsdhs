# # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # Behavioral State-Dependent Habitat Selection # # #
# # # And Implications For Animal Translocations # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # Simona Picardi # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# This code reproduces the analysis published in Picardi et al., 2021 JApplEcol

# Part II, Integrated Step Selection Analysis

# Load packages ####

library(amt)
library(tidyverse)
library(patchwork)

# Load iSSA data ####

nd_amt1 <- readRDS("input/Picardi-et-al_JApplEcol_iSSA-Data-State1.rds")
nd_amt2 <- readRDS("input/Picardi-et-al_JApplEcol_iSSA-Data-State2.rds")

# Seasonal models ####

nd_amt1 <- nd_amt1 %>% 
  dplyr::select(hen_name, rsteps) %>% 
  unnest(cols = rsteps) %>% 
  mutate(hen_season = case_when(
    lubridate::month(t1_) %in% c(4:5) ~ paste0(hen_name, "_spring"),
    lubridate::month(t1_) %in% c(6:8) ~ paste0(hen_name, "_summer"),
    lubridate::month(t1_) %in% c(1:3, 9:12) ~ paste0(hen_name, "_winter")
  )) %>% 
  nest(rsteps = -hen_season)

nd_amt2 <- nd_amt2 %>% 
  dplyr::select(hen_name, rsteps) %>% 
  unnest(cols = rsteps) %>% 
  mutate(hen_season = case_when(
    lubridate::month(t1_) %in% c(4:5) ~ paste0(hen_name, "_spring"),
    lubridate::month(t1_) %in% c(6:8) ~ paste0(hen_name, "_summer"),
    lubridate::month(t1_) %in% c(1:3, 9:12) ~ paste0(hen_name, "_winter")
  )) %>% 
  nest(rsteps = -hen_season)

# Restricted state, without brood
nd_amt1_n <- nd_amt1 %>% 
  mutate(rsteps = lapply(rsteps, FUN = function(x) {
    d <- x %>% 
      filter(with_brood == 0,)
  }))

nd_amt1_n <- nd_amt1_n %>% 
  unnest(rsteps) %>% 
  nest(rsteps = -hen_season)

# Exploratory state, without brood
nd_amt2_n <- nd_amt2 %>% 
  mutate(rsteps = lapply(rsteps, FUN = function(x) {
    d <- x %>% 
      filter(with_brood == 0,)
  }))

nd_amt2_n <- nd_amt2_n %>% 
  unnest(rsteps) %>% 
  nest(rsteps = -hen_season)

# Restricted state, with brood
nd_amt1_b <- nd_amt1 %>% 
  mutate(rsteps = lapply(rsteps, FUN = function(x) {
    d <- x %>% 
      filter(with_brood == 1,)
  }))

nd_amt1_b <- nd_amt1_b %>% 
  unnest(rsteps) %>% 
  nest(rsteps = -hen_season)

# Exploratory state, with brood
nd_amt2_b <- nd_amt2 %>% 
  mutate(rsteps = lapply(rsteps, FUN = function(x) {
    d <- x %>% 
      filter(with_brood == 1,)
  }))

nd_amt2_b <- nd_amt2_b %>% 
  unnest(rsteps) %>% 
  nest(rsteps = -hen_season)

# How many steps in final dataset?
nd_amt1_b %>% 
  unnest(rsteps) %>% 
  bind_rows(unnest(nd_amt1_n, rsteps)) %>% 
  bind_rows(unnest(nd_amt2_b, rsteps)) %>% 
  bind_rows(unnest(nd_amt2_n, rsteps)) %>% 
  filter(case_) %>% 
  nrow()

# Discard individuals with < 16 used steps

toss <- nd_amt1_n %>% 
  unnest(rsteps) %>% 
  bind_rows(unnest(nd_amt1_b, rsteps)) %>% 
  bind_rows(unnest(nd_amt2_n, rsteps)) %>% 
  bind_rows(unnest(nd_amt2_b, rsteps)) %>% 
  filter(case_) %>% 
  group_by(hen_season, state) %>% 
  tally() %>% 
  arrange(n) %>% 
  filter(n < 16) %>% 
  pull(hen_season)

nd_amt1_n <- nd_amt1_n %>%
  filter(!hen_season %in% toss) 

nd_amt1_b <- nd_amt1_b %>%
  filter(!hen_season %in% toss) 

nd_amt2_n <- nd_amt2_n %>%
  filter(!hen_season %in% toss)

nd_amt2_b <- nd_amt2_b %>%
  filter(!hen_season %in% toss)

# No brood
# State 1
issa1_n_s <- nd_amt1_n %>% 
  mutate(issa = lapply(rsteps, FUN = function (x) {
    res <- NA
    res <- try(fit_issf(data = x, 
                        formula = case_ ~  
                          dist_to_mesic_log_scaled +
                          dist_to_roads_log_scaled +
                          dist_to_well_pads_log_scaled +
                          perennial_herb_scaled +
                          sagebrush_scaled +
                          sage_contig_scaled +
                          asp_sin +
                          asp_cos +
                          slope_scaled +
                          sl_ + 
                          log_sl_ + 
                          cos_ta_ +
                          strata(step_id_),
                        model = TRUE))
  }))

# State 2
issa2_n_s <- nd_amt2_n %>% 
  mutate(issa = lapply(rsteps, FUN = function (x) {
    res <- try(fit_issf(data = x, 
                        formula = case_ ~ 
                          dist_to_mesic_log_scaled +
                          dist_to_roads_log_scaled +
                          dist_to_well_pads_log_scaled +
                          perennial_herb_scaled +
                          sagebrush_scaled +
                          sage_contig_scaled +
                          asp_sin +
                          asp_cos +
                          slope_scaled +
                          sl_ + 
                          log_sl_ + 
                          cos_ta_ +
                          strata(step_id_), 
                        model = TRUE))
  }))

# Brood
# State 1
issa1_b_s <- nd_amt1_b %>% 
  mutate(issa = lapply(rsteps, FUN = function (x) {
    res <- NA
    res <- try(fit_issf(data = x, 
                        formula = case_ ~  
                          dist_to_mesic_log_scaled +
                          dist_to_roads_log_scaled +
                          dist_to_well_pads_log_scaled +
                          perennial_herb_scaled +
                          sagebrush_scaled +
                          sage_contig_scaled +
                          asp_sin +
                          asp_cos +
                          slope_scaled +
                          sl_ + 
                          log_sl_ + 
                          cos_ta_ +
                          strata(step_id_),
                        model = TRUE))
  }))

# State 2
issa2_b_s <- nd_amt2_b %>% 
  mutate(issa = lapply(rsteps, FUN = function (x) {
    res <- try(fit_issf(data = x, 
                        formula = case_ ~ 
                          dist_to_mesic_log_scaled +
                          dist_to_roads_log_scaled +
                          dist_to_well_pads_log_scaled +
                          perennial_herb_scaled +
                          sagebrush_scaled +
                          sage_contig_scaled +
                          asp_sin +
                          asp_cos +
                          slope_scaled +
                          sl_ + 
                          log_sl_ + 
                          cos_ta_ +
                          strata(step_id_), 
                        model = TRUE))
  }))
# Not enough information here. 

# Keep model output only

issa1_n_seas <- issa1_n_s %>% 
  dplyr::select(hen_season, issa)

issa2_n_seas <- issa2_n_s %>% 
  dplyr::select(hen_season, issa)

issa1_b_seas <- issa1_b_s %>% 
  dplyr::select(hen_season, issa)

# Get model output in tidy format ####

issa1_b_tidy <- issa1_b_seas %>% 
  mutate(tidy_res = lapply(issa, function(x) {
    if(inherits(x, "fit_clogit")) {
      res <- broom::tidy(x$model)
      return(res)
    } else {
      return(NA)
    }})) %>% 
  filter(!is.na(tidy_res))

issa1_n_tidy <- issa1_n_seas %>% 
  mutate(tidy_res = lapply(issa, function(x) {
    if(inherits(x, "fit_clogit")) {
      res <- broom::tidy(x$model)
      return(res)
    } else {
      return(NA)
    }})) %>% 
  filter(!is.na(tidy_res))

issa2_n_tidy <- issa2_n_seas %>% 
  mutate(tidy_res = lapply(issa, function(x) {
    if(inherits(x, "fit_clogit")) {
      res <- broom::tidy(x$model)
      return(res)
    } else {
      return(NA)
    }})) %>% 
  filter(!is.na(tidy_res))

# Unnest

issa1_b_df <- issa1_b_tidy %>% 
  dplyr::select(hen_season, tidy_res) %>% 
  unnest(cols = tidy_res)

issa1_n_df <- issa1_n_tidy %>% 
  dplyr::select(hen_season, tidy_res) %>% 
  unnest(cols = tidy_res)

issa2_n_df <- issa2_n_tidy %>% 
  dplyr::select(hen_season, tidy_res) %>% 
  unnest(cols = tidy_res) 

# This gives us estimates and CIs for the beta coefficients (log-RSS for a
# 1-unit change in each covariate in isolation)

issa1_b_iw <- issa1_b_df %>% 
  nest(data = c(-term))

issa1_n_iw <- issa1_n_df %>% 
  nest(data = c(-term))

issa2_n_iw <- issa2_n_df %>% 
  nest(data = c(-term)) 

issa1_b_iw <- issa1_b_iw %>% 
  mutate(iw = lapply(data, function(x) {
    
    mod <- lm(estimate ~ 1, data = x, weights = 1/(std.error)^2)
    
    return(mod)
  })) %>% 
  mutate(pred = lapply(iw, function(x) {
    
    pred <- predict(x, 
                    newdata = data.frame(dummy = NA),
                    se.fit = TRUE)
    
    est <- data.frame(mean = pred$fit,
                      lwr = pred$fit - 1.96 * pred$se.fit,
                      upr = pred$fit + 1.96 * pred$se.fit)
    
    return(est)
  }))

issa1_n_iw <- issa1_n_iw %>% 
  mutate(iw = lapply(data, function(x) {
    x <- x %>% 
      mutate(season = stringr::word(hen_season, 2, 2, "_"))
    mod <- lm(estimate ~ season, data = x, weights = 1/(std.error)^2)
    return(mod)
  })) %>% 
  mutate(pred = lapply(iw, function(x) {
    pred <- predict(x, 
                    newdata = data.frame(season = c("spring", 
                                                    "summer", 
                                                    "winter")),
                    se.fit = TRUE)
    est <- data.frame(season = c("Spring", 
                                 "Summer", 
                                 "Winter"),
                      mean = pred$fit,
                      lwr = pred$fit - 1.96 * pred$se.fit,
                      upr = pred$fit + 1.96 * pred$se.fit)
    return(est)
  }))

issa2_n_iw <- issa2_n_iw %>% 
  mutate(iw = lapply(data, function(x) {
    x <- x %>% 
      mutate(season = stringr::word(hen_season, 2, 2, "_"))
    mod <- lm(estimate ~ season, data = x, weights = 1/(std.error)^2)
    return(mod)
  })) %>% 
  mutate(pred = lapply(iw, function(x) {
    pred <- predict(x, 
                    newdata = data.frame(season = c("spring", 
                                                    "summer", 
                                                    "winter")),
                    se.fit = TRUE)
    est <- data.frame(season = c("Spring", 
                                 "Summer", 
                                 "Winter"),
                      mean = pred$fit,
                      lwr = pred$fit - 1.96 * pred$se.fit,
                      upr = pred$fit + 1.96 * pred$se.fit)
    return(est)
  }))

iw_1b <- issa1_b_iw %>% 
  dplyr::select(term, pred) %>% 
  unnest(cols = pred)

iw_1n <- issa1_n_iw %>% 
  dplyr::select(term, pred) %>% 
  unnest(cols = pred)

iw_2n <- issa2_n_iw %>% 
  dplyr::select(term, pred) %>% 
  unnest(cols = pred) 

# Plot parameter estimates ####

covs <- c("Aspect (cosine)",
          "Aspect (sine)",
          "Distance to mesic habitat",
          "Distance to roads",
          "Distance to well pads",
          "Perennial herbaceous cover",
          "Sagebrush contiguity",
          "Sagebrush cover",
          "Slope")

iw_1b_p <- iw_1b %>% 
  filter(!term %in% c("cos_ta_", "log_sl_", "sl_")) %>% 
  ggplot(aes(x = term, y = mean, color = term)) +
  geom_point() +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "Restricted state, with brood", x = " ", y = "log-RSS") +
  scale_color_discrete(name = "Covariate", labels = covs) +
  theme_bw() +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
  coord_cartesian(ylim = c(-0.75, 0.65))

iw_1n_p <- iw_1n %>% 
  filter(!term %in% c("cos_ta_", "log_sl_", "sl_")) %>% 
  ggplot(aes(x = term, y = mean, color = term)) +
  geom_point() +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.2) +
  facet_wrap(~ season) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "Restricted state, without brood", x = " ", y = "log-RSS") +
  scale_color_discrete(name = "Covariate", labels = covs) +
  theme_bw() +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
  coord_cartesian(ylim = c(-0.75, 0.65))

iw_2n_p <- iw_2n %>% 
  filter(!term %in% c("cos_ta_", "log_sl_", "sl_")) %>% 
  ggplot(aes(x = term, y = mean, color = term)) +
  geom_point() +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.2) +
  facet_wrap(~ season) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "Exploratory state, without brood", x = " ", y = "log-RSS") +
  scale_color_discrete(name = "Covariate", labels = covs) +
  theme_bw() +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
  coord_cartesian(ylim = c(-0.75, 0.65))

iw_2n_p / iw_1n_p / iw_1b_p + plot_layout(guides = "collect") & 
  theme(legend.position = "bottom") 
