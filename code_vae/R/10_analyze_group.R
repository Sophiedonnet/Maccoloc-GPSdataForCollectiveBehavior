#===============================================================================
# This script analyze the trajectories, latent space and income of the fleet per
# group (last part and plot of the manuscript).
# (Figure 11 of the paper, files latent_space_group.png and
#  income_per_trip_noguide.png)
#
# It does the following thing:
# - Import fitted colSBM model, and the groups of each vessel in each quarter
# - Import latent space of cvae3, trajectories and catch data
# - Assign to each trip the group of its vessel during the quarter
# - Figure 11 (b): income per trip, by quarter and group
# - Figure 11 (a): trajectories (top) and latent means mu_1, mu_2 (bottom)
#   of the trips of each group
#
# Inputs:
# - Rdata/sacrois_xgv_capt.csv (09_extract_xgv_capture.R): catch data
# - Rdata/xgv_cvae_3_reordered_filtered.csv (03_computing_bc.R)
# - Rdata/trim_peche.rds and Rdata/save_sbmfit/iidcolsbm_trim_round8.rds
#   (07b_fit_colsbm_proxy_graphs_replicate_results.R)
# - Rdata/xgv_100_14-18.csv (02_merge_xgv.R)
#
# Outputs:
# - images2/Ammount_per_trip_per_group.png: Figure 11 (b)
#   (income_per_trip_noguide.png in the paper)
# - images2/latent_space_group.png: Figure 11 (a)
# - Rdata/sacrois_xgv_capt_group.csv: positions of the trips with their group
#
# The income figure for the earlier fits is in archive/images2/.
#===============================================================================

library(tidyverse)
library(lubridate)
library(ggplotFL)
library(patchwork)
source("./99_helper_functions.R")
source('./98_functions_resampling.R')
#library(devtools)
#install_github("flr/ggplotFL")
#install_github("flr/FLCore")


#===============================================================================
# Data
#===============================================================================
df_xgv_capt <- read_csv("Rdata/sacrois_xgv_capt.csv",show_col_types = FALSE)
df_cvae <- read_csv(file = "./Rdata/xgv_cvae_3_reordered_filtered.csv",show_col_types = FALSE)
trim_peche <- readRDS("./Rdata/trim_peche.rds")
df_traj <- read_csv("Rdata/xgv_100_14-18.csv",show_col_types = FALSE)
df_traj <- df_traj %>% filter(wday(time) %in% c(2,3,4,5,6)) |>
  filter(maree != 12097796)

# Selected model of the colSBM fit of 07b_...R (6 groups). Earlier fits are
# in archive/ (see archive/README.md).
iidcolsbm_trim <- readRDS("./Rdata/save_sbmfit/iidcolsbm_trim_round8.rds")$best_fit

Q <-iidcolsbm_trim$nb_clusters

# Groups of each vessel in each quarter, relabelled in the plotting order
# (as in 08_analyse_behaviour.R)
best_ord = get_order(iidcolsbm_trim)
Z_trim <- get_Ztrim_relabel(iidcolsbm_trim)
Z_trim <- relabel_Ztrim(Z_trim,best_ord)
Z_trim_map <- get_Ztrim_map(Z_trim)

navire_id <- df_cvae |> pull(navire) |> unique()

# Group of each trip = group of its vessel during the quarter
# (quarteryear = 1..20 from Q1 2014)
df_cvae_group <- df_cvae %>%
  mutate(quarteryear = quarter(time) + 4*(year(time) - 2014)) %>% 
  mutate(navire = factor(navire, levels = navire_id)) %>% 
  group_by(maree) %>% 
  mutate(Group = Z_trim[navire, quarteryear])%>%mutate(quarter=factor(quarter(time)))


#===============================================================================
# Income per trip and per group (Figure 11 (b))
#===============================================================================
# Gross value per trip (sum over the landed species, catch lines above
# 10 000 euros removed) by quarter and group: black line = median, shades =
# [10%, 90%] and [25%, 75%] quantiles, dotted = 1% and 99%.
# Each trip is dated at the last month of its quarter.
df_xgv_capt2 <- df_xgv_capt %>%
  #  mutate(MONTANT_EUROS_MOYENNE = pmin(MONTANT_EUROS_MOYENNE, 8000)) %>% 
  filter(MONTANT_EUROS_MOYENNE <= 10000) %>% 
  rename(maree = MAREE_ID) %>% 
  left_join(df_cvae_group) %>%
  mutate(quarter = factor(quarter(time)),
         Group = LETTERS[Group]) %>% 
  mutate(date = as_date(paste(year(time), 3*quarter(time),"01", sep = "-"))) %>% 
  group_by(maree, Group, quarter, date) %>% 
  summarise(MONTANT = sum(MONTANT_EUROS_MOYENNE, na.rm = TRUE))

# Mean income per trip of each group
rev_moyen <- df_xgv_capt2%>% group_by(Group) %>% summarise(MONTANT = mean(MONTANT, na.rm = TRUE))
print(rev_moyen)

ggplot_incomes <- ggplot(df_xgv_capt2,aes(y = MONTANT, x = date, fill = Group, col = Group)) +
  #  geom_smooth() +
  #       stat_ecdf(linewidth = 1.1) +
  stat_flquantiles(probs=c(0.10, 0.90), geom = "ribbon",
                   alpha=0.30) +
  stat_flquantiles(probs=c(0.01), geom = "line",
                   linetype=3) +
  stat_flquantiles(probs=c(0.99), geom = "line",
                   linetype=3) +
  stat_flquantiles(probs=c(0.25, 0.75), geom = "ribbon",
                   alpha=0.30) +
  stat_flquantiles(probs=c(0.50), geom = "line", size=1.5) +
  stat_flquantiles(probs=c(0.50), geom = "line", col = "black") +
  #    geom_boxplot(position = "dodge") +
  ylim(c(0, 5000)) + 
  xlab("Date") +
  ylab("Income per trip (Euros)") +
  facet_wrap(~ Group, ncol = Q) +
  theme_bw(base_size = 20) #+

ggplot_incomes + ylab('') + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave("./images2/Ammount_per_trip_per_group.png", height = 8, width = 12)


#===============================================================================
# Trajectories and latent space per group (Figure 11 (a))
#===============================================================================

# Bottom: 2d histogram of the latent means (mu_1, mu_2) of the trips of each
# group
p_lat <- df_cvae_group %>%
  mutate(Group = factor(LETTERS[Group], levels = LETTERS[1:Q])) %>% 
  na.omit() %>% 
  #  filter(Group %in% c(2,4)) %>% 
  ggplot(aes(x = mu_1, y = mu_2,  fill = Group)) + 
  #  geom_abline() +
  stat_bin_2d(aes(alpha = after_stat(ndensity)), show.legend = FALSE) +
  coord_fixed() + ylim(-2.5, 2.5) +
  facet_wrap(~ Group, ncol = Q) +
  theme_bw(base_size = 20)  +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )
# Same, by quarter of the year (exploratory, not in the paper)
#%filter(quarteryear %in%   c(1:4)) %>
p_lat2 <- df_cvae_group %>%
  mutate(Group = factor(LETTERS[Group], levels = LETTERS[1:Q])) %>% 
  na.omit() %>% 
  #  filter(Group %in% c(2,4)) %>% 
  ggplot(aes(x = mu_1, y = mu_2,  fill = Group)) + 
  #  geom_abline() +
  stat_bin_2d(aes(alpha = after_stat(ndensity)), show.legend = FALSE) +
  coord_fixed() + ylim(-2.5, 2.5) +
  facet_grid(quarter~ Group) +
  theme_bw(base_size = 20)  +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )
p_lat2


# Hourly positions of each trip with its group (20511 trips x 24 hours)
df_traj_group <- df_traj %>%  mutate(navire = as_factor(navire)) |>
  left_join(df_cvae_group %>% mutate(`...1` = `...1` + 1), 
            by = c("month", "maree", "time", "navire")) %>% 
  select(x, y, Group, time, maree) %>% 
  mutate(h = rep(seq(24), times = 20511)) %>% 
  mutate(Group = factor(LETTERS[Group], levels = LETTERS[1:Q])) %>% 
  na.omit()

write_csv(df_traj_group, "Rdata/sacrois_xgv_capt_group.csv")


# Top: trajectories of each group, over all the trajectories in grey
p_traj <- df_traj_group %>%
  ggplot(aes(x = x, y = y, col = Group, group = maree)) + #, size = 1 + .01*mu_3)) +
  #    annotation_map(map = fortify(maps::map("france", xlim = c(-6, -2.5), ylim = c(47, 48.3)))) +
  annotate("path", x = df_traj_group$x, y = df_traj_group$y, 
           col = "gray70", alpha = .5) + #, group = df_traj_std_group$maree) + 
  geom_path(show.legend = FALSE, alpha = .1) +
  xlab("y_1") + ylab("y_2") +
  # gghighlight::gghighlight(,
  #                          calculate_per_facet = TRUE) + 
  #  geom_point(aes(alpha = h/24), show.legend = FALSE) +
  coord_equal() +
  facet_wrap(~ Group, ncol = Q) +
  theme_bw(base_size = 20)

p_all <- (p_traj/p_lat ) + plot_layout(heights = unit(c(-1,2),c("null","null")))
p_all <- cowplot::align_plots(p_traj, p_lat)
p_all
ggsave("./images2/latent_space_group.png", plot = p_all, height = 7.5, width = 15)
