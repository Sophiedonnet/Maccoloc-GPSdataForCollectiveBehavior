#===============================================================================
# This script does the following
# - Import the reorder and filted latent representation dataset of cvae3
# - Import the trajectory dataset of XGV
# - Plot the trajectories, with colour depending on the parameters of their
# latent variational distribution
#===============================================================================

#===============================================================================
# Import function and library
#===============================================================================
library(tidyverse)
library(lubridate)
library(matrixStats)
library(patchwork)
source("./99_helper_functions.R")

df_traj <- read_csv(file = "./Rdata/xgv_100_14-18.csv")
df_traj <- df_traj %>% filter(wday(time) %in% c(2,3,4,5,6)) |>
  filter(maree != 12097796)
df_cvae <- read_csv(file = "./Rdata/xgv_cvae_3_reordered_filtered.csv")

df_traj <- df_traj |> 
  mutate(x_std = coord2std(x, -4.28),
         y_std = coord2std(y, 47.78)) 


#===============================================================================
# Plotting the trajectory data in function of the color in the latent space
#===============================================================================
mu_names <- c(
  "mu[1]",
  "mu[2]",
  "mu[3]"
)


df_traj %>% left_join(df_cvae %>% mutate(`...1` = `...1` + 1), by = "maree") %>% 
  #filter(sig_2 >= .2) %>% 
  mutate(h = rep(seq(24), times = 492264/24)) %>% 
  #  mutate(Group = factor(LETTERS[Group], levels = LETTERS[1:5])) %>% 
  na.omit() %>% 
  pivot_longer(cols = starts_with("mu"), names_to = "Dim", values_to = "mu") %>% 
  mutate(Dim = factor(Dim, labels = c("mu_1", "mu_2", "mu_3"))) %>% 
  ggplot(aes(x = x, y = y, col = mu, group = maree)) + 
  annotation_map(map = fortify(maps::map("france", xlim = c(-5.3, -3.6), ylim = c(47, 48.3)))) +
  geom_path(alpha = .1) + 
  scale_colour_viridis_c(direction =  - 1) + 
  annotate("point", x = -4.28, y = 47.78, shape = 23, size = 6, fill = "red") +
  facet_wrap(~ Dim, labeller = label_parsed) + 
  coord_equal() + xlab("Longitude") + ylab("Latitude") +
  labs(color=expression(mu)) +
  theme_bw(base_size = 15)
# ylim(-5, 1) + 
# xlim(-7, 5)



#ggsave("../figures/plot_traj_mu.png")


df_traj %>% left_join(df_cvae %>% mutate(`...1` = `...1` + 1), by = "maree") %>% 
  #filter(sig_2 >= .2) %>% 
  mutate(h = rep(seq(24), times = 492264/24)) %>% 
  #  mutate(Group = factor(LETTERS[Group], levels = LETTERS[1:5])) %>% 
  na.omit() %>% 
  pivot_longer(cols = starts_with("sig"), names_to = "Dim", values_to = "sig") %>% 
  mutate(Dim = factor(Dim, labels = c("gamma^2_1", "gamma^2_2", "gamma^2_3"))) %>% 
  ggplot(aes(x = x, y = y, col = sig, group = maree)) + 
  annotation_map(map = fortify(maps::map("france", xlim = c(-5.3, -3.6), ylim = c(47, 48.3)))) +
  geom_path(alpha = .1) + 
  annotate("point", x = -4.28, y = 47.78, shape = 23, size = 6, fill = "red") +
  scale_colour_viridis_c(direction =  - 1) + 
  facet_wrap(~ Dim) + 
  coord_equal() + xlab("Longitude") + ylab("Latitude") +
  labs(color=expression(gamma^2)) +
  theme_bw(base_size = 15) +
  theme(text=element_text(family="Helvetica"))
# ylim(-5, 1) + 
# xlim(-7, 5)
# ylim(-5, 1) + 
# xlim(-7, 5)



ggsave("../figures/plot_traj_sig.png", width = 12.3, height = 3.79)
