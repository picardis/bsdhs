# # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # Behavioral State-Dependent Habitat Selection # # #
# # # And Implications For Animal Translocations # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # Simona Picardi # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# This code reproduces the analysis published in Picardi et al., 2021 JApplEcol

# Part I, Hidden Markov Model 

# Load packages ####

library(momentuHMM)
library(tidyverse)
library(lubridate)
library(patchwork)

# Load HMM data ####

hmm_data <- readRDS("input/Picardi-et-al_JApplEcol_HMM-Data.rds")

# Run HMM ####

form <- ~ dst_scaled
formD <- ~ status

hmm <- fitHMM(data = hmm_data,
               nbStates = 2,
               formula = form, 
               formulaDelta = formD, 
               dist = list(step = "gamma", angle = "vm"),
               estAngleMean = list(angle = TRUE),
               Par0 = list(step = c(mean_1 = 100, 
                                    mean_2 = 1500, 
                                    sd_1 = 100, 
                                    sd_2 = 10000,  
                                    zeromass_1 = 0, 
                                    zeromass_2 = 0),
                           angle = c(mean_1 = 3, 
                                     mean_2 = 0.001,
                                     concentration_1 = 0.1,
                                     concentration_2 = 0.99)))

# Assign behavioral state using Viterbi algorithm ####

hmm_data <- hmm_data %>%
  as.data.frame() %>% 
  mutate(state = viterbi(hmm)) 

# Plot sequential phases ####

hmm_data %>% 
  group_by(hen_name) %>% 
  mutate(release_date = as_date(first(timestamp)),
         state_words = case_when(
           state == 1 ~ "Restricted",
           state == 2 ~ "Exploratory"
         )) %>% 
  ungroup() %>% 
  # To plot the full time extent, comment out line below 
  filter(timestamp <= ymd(paste0(year(release_date), "12-31"))) %>% 
  mutate(julian = as.numeric(format(timestamp, "%j")),
         tyear = year(timestamp) - year(release_date),
         tjulian = julian + (365 * tyear)) %>%
  ggplot(aes(y = hen_name, x = tjulian, color = state_words)) + 
  geom_line(size = 2) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(x = "Julian day", y = " ", color = "State")

# Plot model predictions ####

mean_dst <- mean(hmm_data$dst)
sd_dst <- sd(hmm_data$dst)

# Transition probabilities

b <- hmm$mle$beta

tp <- data.frame(dst = 1:200) %>% 
  mutate(dst_stand = (dst - mean_dst)/sd_dst) %>% 
  mutate(en_to_ex = b[1,1] + 
           b[2,1] *  dst_stand,         
         ex_to_en = b[1,2] + 
           b[2,2] * dst_stand) %>% 
  mutate(en_to_en = 1 - en_to_ex,
         ex_to_ex = 1 - ex_to_en) %>% 
  mutate_at(vars(en_to_ex:ex_to_ex), plogis) %>% 
  pivot_longer(cols = en_to_ex:ex_to_ex, names_to = "tp", values_to = "value")

facet_labs <- c("en_to_en" = "Staying in restricted",
                "en_to_ex" = "Switching to exploratory",
                "ex_to_en" = "Switching to restricted",
                "ex_to_ex" = "Staying in exploratory")

tp_switch <- tp %>% 
  filter(tp %in% c("en_to_ex", "ex_to_en"))

ggplot(tp_switch, aes(x = dst, y = value)) +
  geom_line(size = 1) +
  facet_wrap(~ tp, labeller = as_labeller(facet_labs)) +
  labs(x = "Days since translocation", y = "Transition probability") +
  theme_bw()

# Initial state probabilities

ind_status <- hmm_data %>% 
  dplyr::select(ID, status) %>% 
  distinct()

isp <- data.frame(ID = word(rownames(hmm$CIreal$delta$lower), 2, 2, ":"),
                  lwr = c(hmm$CIreal$delta$lower[,1],
                          hmm$CIreal$delta$lower[,2]),
                  mle = c(hmm$CIreal$delta$est[,1],
                          hmm$CIreal$delta$est[,2]),
                  upr = c(hmm$CIreal$delta$upper[,1],
                          hmm$CIreal$delta$upper[,2]),
                  state = rep(c(1, 2), each = nrow(hmm$CIreal$delta$lower)),
                  row.names = NULL) %>% 
  left_join(ind_status) %>% 
  dplyr::select(-ID) %>% 
  distinct() %>% 
  mutate(state_descript = case_when(
    state == 1 ~ "Restricted",
    state == 2 ~ "Exploratory"
  ),
  status_descript = case_when(
    status == "NBH" ~ "Non-brood",
    status == "PNH" ~ "Pre-nesting",
    status == "TBH" ~ "Brood"
  ))

ggplot(isp, aes(x = state_descript, y = mle, color = factor(state_descript))) +
  geom_point(size = 2) + 
  geom_errorbar(mapping = aes(ymin = lwr, ymax = upr), width = 0.2) +
  facet_wrap(~ status_descript) + 
  theme_bw() +
  theme(axis.text.x = element_blank(), legend.position = "bottom") +
  labs(x = element_blank(), 
       y = "Initial state probability",
       color = "State")

# Distributions by state ####

# Step length

# Get canonical parameters of gamma distr from mean and variance

# State 1
shape1 <- hmm$mle$step[1,1]^2/hmm$mle$step[2,1]^2
scale1 <- hmm$mle$step[2,1]^2/hmm$mle$step[1,1]

# State 2
shape2 <- hmm$mle$step[1,2]^2/hmm$mle$step[2,2]^2
scale2 <- hmm$mle$step[2,2]^2/hmm$mle$step[1,2]

# Simulate and plot
sim_step <- data.frame(x = 0:10000,
                       sim = c(dgamma(x = 0:10000, shape = shape1, scale = scale1),
                               dgamma(x = 0:10000, shape = shape2, scale = scale2)),
                       state = rep(c("Restricted", "Exploratory"), each = 10001)) %>% 
  mutate(state = factor(state, levels = c("Restricted", "Exploratory")))

# Turning angles

# Simulate and plot
sim_angle <- data.frame(x = seq(-pi, pi, length.out = 100),
                        sim = c(CircStats::dvm(theta = seq(-pi, pi, length.out = 100), 
                                               mu = hmm$mle$angle[1,1], 
                                               kappa = hmm$mle$angle[2,1]),
                                CircStats::dvm(theta = seq(-pi, pi, length.out = 100), 
                                               mu = hmm$mle$angle[1,2], 
                                               kappa = hmm$mle$angle[2,2])),
                        state = rep(c("Restricted", "Exploratory"), each = 100)) %>% 
  mutate(state = factor(state, levels = c("Restricted", "Exploratory")))

# Overlay observed data with fitted distribution in each state
hmm_data <- hmm_data %>% 
  mutate(state_words = case_when(
    state == 1 ~ "Restricted",
    state == 2 ~ "Exploratory"
  ))

sl_w_hist <- hmm_data %>% 
  mutate(state = case_when(
    state == 1 ~ "Restricted",
    state == 2 ~ "Exploratory"
  )) %>% 
  ggplot(aes(x = step)) +
  geom_histogram(aes(y = ..density.., fill = factor(state)), bins = 80, color = "black") +
  geom_density(data = sim_step, mapping = aes(x = x, y = sim), stat = "identity",
               lwd = 0.8) +
  facet_wrap(~ state) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "Step length (m)", y = "Frequency") +
  coord_cartesian(ylim = c(0, 0.004), xlim = c(0, 6500))

ta_w_hist <- hmm_data %>% 
  mutate(state = case_when(
    state == 1 ~ "Restricted",
    state == 2 ~ "Exploratory"
  )) %>% 
  ggplot(aes(x = angle)) +
  geom_histogram(aes(y = ..density.., fill = factor(state)), bins = 16, color = "black") +
  geom_line(data = sim_angle, mapping = aes(x = x, y = sim), stat = "identity",
            lwd = 0.8) +
  facet_wrap(~ state) +
  theme_bw() +
  theme(legend.position = "none") +
  scale_x_continuous(breaks = c(0, pi/2, pi, -pi/2), 
                     labels = c("0", "90", "180", "270")) +
  labs(y = "Frequency", x = "Turning angle (degrees)") +
  coord_polar(start = pi) 

sl_w_hist / ta_w_hist
