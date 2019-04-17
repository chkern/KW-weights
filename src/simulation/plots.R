# Setup
library(tidyverse)

est  <- read.table("est.txt")
wt_m <- read.table("wt_m.txt")
wt_v <- read.table("wt_v.txt")
p_scores <- read.table("p_scores.txt")
m.pop_y <- 0.18383

est <-
  est %>%
  rename(method = V2, model = V3, est = V4) %>%
  group_by(method, model) %>% 
  mutate(simu = row_number()) %>%
  ungroup %>%
  mutate(method = recode(method,
                         "A" = "Cht (naïve)",
                         "B" = "Cht (wtd)",
                         "C" = "Svy (wtd)",
                         "D" = "IPSW (true)",
                         "E" = "IPSW (main)", 
                         "F" = "IPSW (pairwise)",
                         "G" = "KW (true)", 
                         "H" = "KW (main)",
                         "I" = "KW (pairwise)",
                         "J" = "KW (MOB)",
                         "K" = "KW (RF)",
                         "L" = "KW (CRF)",
                         "M" = "KW (XTree)",
                         "N" = "KW (GBM)",
                         "O" = "KW (MBoost)")) %>%
  mutate(model = recode(model,
                        "A" = "Additive + Linear",
                        "B" = "Additive + Mild non-linear",
                        "C" = "Additive + Moderate non-linear",
                        "D" = "Mild non-additive + Linear",
                        "E" = "Mild non-additive + Mild non-linear",
                        "F" = "Mod non-additive + Linear",
                        "G" = "Mod non-additive + Mod non-linear",
                        "H" = "Strong non-additive + Strong non-linear",
                        "I" = "Transformation"))

wt <-
  wt_m %>%
  bind_cols(wt_v) %>%
  select(V1, V2, V3, V4, V41) %>%
  rename(method = V2, model = V3, wt_m = V4, wt_v = V41) %>% 
  group_by(method, model) %>%
  mutate(simu = row_number()) %>%
  ungroup %>%
  filter(method != "A") %>%
  mutate(method = recode(method,
                         "B" = "Cht (wtd)",
                         "C" = "Svy (wtd)",
                         "D" = "IPSW (true)",
                         "E" = "IPSW (main)", 
                         "F" = "IPSW (pairwise)",
                         "G" = "KW (true)", 
                         "H" = "KW (main)",
                         "I" = "KW (pairwise)",
                         "J" = "KW (MOB)",
                         "K" = "KW (RF)",
                         "L" = "KW (CRF)",
                         "M" = "KW (XTree)",
                         "N" = "KW (GBM)",
                         "O" = "KW (MBoost)")) %>%
  mutate(model = recode(model,
                        "A" = "Additive + Linear",
                        "B" = "Additive + Mild non-linear",
                        "C" = "Additive + Moderate non-linear",
                        "D" = "Mild non-additive + Linear",
                        "E" = "Mild non-additive + Mild non-linear",
                        "F" = "Mod non-additive + Linear",
                        "G" = "Mod non-additive + Mod non-linear",
                        "H" = "Strong non-additive + Strong non-linear",
                        "I" = "Transformation"))

# Estimate by method and model (boxplots over simu)

ggplot(est) +
  geom_boxplot(aes(x = method, y = est), outlier.size = 0.1) +
  facet_wrap(~ model, ncol = 3) +
  geom_hline(yintercept = m.pop_y) +
  coord_flip() +
  labs(x = "", y = "") +
  scale_x_discrete(limits = rev(levels(est$method))) +
  theme_bw() +
  theme(legend.position = "none")

ggsave("simu_p01.png", width = 11, height = 11)

# Relative bias by method and model (boxplots over simu)

est %>%
  mutate(bias = (est - m.pop_y)/m.pop_y*100) %>%
  ggplot() +
  geom_boxplot(aes(x = method, y = bias), outlier.size = 0.1) +
  facet_wrap(~ model, ncol = 3, scales = "free_x") +
  geom_hline(yintercept = 0) +
  coord_flip() +
  labs(x = "", y = "") +
  scale_x_discrete(limits = rev(levels(est$method))) +
  theme_bw() +
  theme(legend.position = "none")

ggsave("simu_p02.png", width = 11, height = 11)

# Squared error by model and method (mean and std over simu)

est %>%
  mutate(err = (est - m.pop_y)^2) %>%
  ggplot(aes(x = model, y = err)) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), 
               geom = "pointrange") +
  facet_wrap(~ method, ncol = 3, scales = "free_x") +
  geom_hline(yintercept = 0) +
  coord_flip() +
  labs(x = "", y = "") +
  scale_x_discrete(limits = rev(levels(est$model))) +
  theme_bw() +
  theme(legend.position = "none")

ggsave("simu_p03.png", width = 11, height = 11)

# CV of weights by model and method (mean and std over simu)

wt %>%
  mutate(cv = sqrt(wt_v)/wt_m) %>%
  ggplot(aes(x = model, y = cv)) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), 
               geom = "pointrange") +
  facet_wrap(~ method, ncol = 3) +
  geom_hline(yintercept = 0) +
  coord_flip() +
  labs(x = "", y = "") +
  scale_x_discrete(limits = rev(levels(wt$model))) +
  theme_bw() +
  theme(legend.position = "none")

ggsave("simu_p04.png", width = 11, height = 11)
