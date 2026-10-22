#===============================================================================
# Latent representations of 3 selected vessels
# (Figures 5 and 6 of the paper, files plot_ts_mu.png and
#  plot_traj_latent_2-6-19_new.png)
#
# This script does the following
# - Import the reorder and filted latent representation dataset of vae4 and
# cvae3
# - Import the trajectory dataset of XGV
# - Standardized the geographical coordinates
# - Selects 3 vessels (2, 6 and 19)
# - Figure 5: time series of the latent distributions of the 3 vessels
# - Figure 6: trajectories and latent space of the 3 vessels
#
# Inputs:
# - Rdata/xgv_100_14-18.csv (02_merge_xgv.R): hourly positions of each trip
# - Rdata/xgv_{cvae_3,vae_4}_reordered_filtered.csv (03_computing_bc.R)
#
# Outputs: images2/plot_ts_mu.png (Figure 5) and
# images2/plot_traj_latent_2-6-19_new.png (Figure 6)
#
# Note: the vessels are numbered 1..33 in their order of appearance in df_cvae;
# this numbering is used in all the figures of the paper.
#===============================================================================

#===============================================================================
# Import function and library
#===============================================================================
library(tidyverse)
library(lubridate)
library(matrixStats)
library(patchwork)
source("./99_helper_functions.R")

#===============================================================================
# Import and filter dataset
#===============================================================================
# Trips on week days only (wday 2..6 = Monday..Friday), and one trip removed
# (12097796), as for the training of the VAEs
df_traj <- read_csv(file = "./Rdata/xgv_100_14-18.csv")
df_traj <- df_traj %>% filter(wday(time) %in% c(2,3,4,5,6)) |>
  filter(maree != 12097796)
df_cvae <- read_csv(file = "./Rdata/xgv_cvae_3_reordered_filtered.csv")
df_vae <- read_csv(file = "./Rdata/xgv_vae_4_reordered_filtered.csv")

navire_id <- unique(df_cvae$navire)

# The 3 vessels shown in the paper
b = 2
bp = 6
bpp = 19

#===============================================================================
# Plotting the latent distribution time series (Figure 5)
#===============================================================================
# One panel per vessel (b, bp, bpp), one row per latent dimension: mean mu of
# the variational distribution of each trip, and band mu +/- z * gamma, with
# gamma = sqrt(sig) the standard deviation.
# Note: z = qnorm(.95) ~ 1.645 here, while the caption of the paper states 1.96.
p_ts_mu <- df_vae %>% filter(navire == navire_id[b]) %>% 
  mutate(model = "VAE") %>%
  bind_rows(df_cvae %>% filter(navire == navire_id[b]) %>% 
              mutate(model = "CVAE")) %>% 
  pivot_longer(cols = starts_with("mu_") | starts_with("sig"), names_sep = "_",
               names_to = c("moment", "dim")) %>%
  pivot_wider(names_from = moment, values_from = value) %>% 
  ggplot(aes(x = time, y = mu, col = model, fill = model,lty = model)) + 
  geom_ribbon(aes(ymin = mu - qnorm(.95)*sqrt(sig), 
                  ymax = mu + qnorm(.95)*sqrt(sig)), 
              linetype = "dotted", alpha = .3, linewidth = .01) +
  scale_color_manual(values = c("tomato", "steelblue")) +
  scale_fill_manual(values = c("tomato", "steelblue")) +
  geom_line(alpha = .8) + 
  ylim(-2.5, 2.5) +
  xlim(as.Date("2014-01-01"), as.Date("2018-12-31")) +
  xlab("") + 
  ylab("Z") + 
  facet_grid( dim ~ .) +
  ggtitle(paste0(b)) +
  theme_bw(base_size = 20) +
  (df_vae %>% filter(navire == navire_id[bp]) %>% 
     mutate(model = "VAE") %>%
     bind_rows(df_cvae %>% filter(navire == navire_id[bp]) %>% 
                 mutate(model = "CVAE")) %>% 
     pivot_longer(cols = starts_with("mu_") | starts_with("sig"), names_sep = "_",
                  names_to = c("moment", "dim")) %>%
     pivot_wider(names_from = moment, values_from = value) %>% 
     ggplot(aes(x = time, y = mu, col = model, fill = model,lty = model)) + 
     geom_ribbon(aes(ymin = mu - qnorm(.95)*sqrt(sig), 
                     ymax = mu + qnorm(.95)*sqrt(sig)), 
                 linetype = "dotted", alpha = .3, linewidth = .01) +
     scale_color_manual(values = c("tomato", "steelblue")) +
     scale_fill_manual(values = c("tomato", "steelblue")) +
     ylab("") +
     geom_line(alpha = .8) + 
     ylim(-2.5, 2.5) +
     xlim(as.Date("2014-01-01"), as.Date("2018-12-31")) +
     facet_grid( dim ~ .) +
     ggtitle(paste0(bp)) +
     theme_bw(base_size = 20)) + 
  (df_vae %>% filter(navire == navire_id[bpp]) %>% 
     mutate(model = "VAE") %>%
     bind_rows(df_cvae %>% filter(navire == navire_id[bpp]) %>% 
                 mutate(model = "CVAE")) %>% 
     pivot_longer(cols = starts_with("mu_") | starts_with("sig"), names_sep = "_",
                  names_to = c("moment", "dim")) %>%
     pivot_wider(names_from = moment, values_from = value) %>% 
     ggplot(aes(x = time, y = mu, col = model, fill = model,lty = model)) + 
     geom_ribbon(aes(ymin = mu - qnorm(.95)*sqrt(sig), ymax = mu + qnorm(.95)*sqrt(sig)), 
                 linetype = "dotted", alpha = .3, linewidth = .01) +
     scale_color_manual(values = c("tomato", "steelblue")) +
     scale_fill_manual(values = c("tomato", "steelblue")) +
     geom_line(alpha = .8) +
     ylab("") + xlab("") + 
     ylim(-2.5, 2.5) +
     xlim(as.Date("2014-01-01"), as.Date("2018-12-31")) +
     facet_grid( dim ~ .) +
     ggtitle(paste0(bpp)) +
     theme_bw(base_size = 20)) +   
  plot_layout(guides = "collect")
p_ts_mu
ggsave("./images2/plot_ts_mu.png", p_ts_mu, width = 16, height = 9, dpi = 300)




#===============================================================================
# Plotting the latent mean and the standardized trajectories (Figure 6)
#===============================================================================
# Left: trajectories of the 3 vessels, in standardized coordinates (offset in
# longitude and latitude from the harbour at (-4.28, 47.78), converted to
# units of 10 km, see coord2std).
# Right: latent means of their trips, x = mu_1, y = mu_2, size = mu_3.
# Trips of June 2015 and December 2015 are highlighted.

df_traj <- df_traj |>
  mutate(x_std = coord2std(x, -4.28),
        y_std = coord2std(y, 47.78)) 

boat <- c(2, 6, 19)
p_traj_latent <- (df_traj   %>% 
    filter(navire %in% navire_id[boat]) %>% 
    mutate(maree = as_factor(maree),
           #navire = as_factor(navire),
           month = as_factor(month(time)),
           year = year(time)) %>% 
    mutate(h = rep(seq(24), nrow(.)/24)) %>% # hour of the trip (24 positions per trip)
    #mutate(metier = df_cvae$metier[df_cvae$maree == maree]) %>% 
    mutate(maree = as_factor(maree)) %>% 
    mutate(navire = case_match(navire,
                               navire_id[2] ~ 2,
                               navire_id[6] ~ 6,
                               navire_id[19] ~ 19#,
                               #navire_id[20] ~ 20
    )) %>% 
    mutate(navire = as_factor(navire)) %>% 
    #filter(time < "2014-01-29" & wday(time) !=1) %>% 
    ggplot(aes(x = x_std, y = y_std , col = month)) +
    geom_path(linewidth = .5, alpha = .5, show.legend = FALSE) +
    geom_point(aes(alpha = month), size = 1, show.legend = FALSE) + 
    facet_grid(navire ~ .) + coord_fixed() +
    xlim(-7.5, 5) + ylim(-6,2) +
    gghighlight::gghighlight((month %in% c("6","12") & year == "2015"),
                             calculate_per_facet = TRUE) +
    
    scale_color_manual(values = c("tomato", "steelblue")) +
    theme_bw(base_size = 20)) |
  (df_vae %>%   filter(navire %in% navire_id[boat]) %>% 
     mutate(maree = as_factor(maree),
            #navire = as_factor(navire),
            month = as_factor(month(time)),
            year = year(time),
            Model = "VAE") %>%
     bind_rows(df_cvae %>%   filter(navire %in% navire_id[boat]) %>% 
                 mutate(maree = as_factor(maree),
                        #navire = as_factor(navire),
                        month = as_factor(month(time)),
                        year = year(time),
                        Model = "CVAE")) %>% 
     mutate(navire = case_match(navire,
                                navire_id[2] ~ 2,
                                navire_id[6] ~ 6,
                                navire_id[19] ~ 19#,
                                #navire_id[20] ~ 20
     )) %>% 
     mutate(navire = as_factor(navire)) %>% 
     ggplot(aes(x = mu_1, y = mu_2, size = mu_3,
                col = month)) +
     geom_point(show.legend = FALSE) +   facet_grid(navire ~ Model) +
     coord_fixed() +
     gghighlight::gghighlight((month %in% c("12","6") & year == "2015"), 
                              calculate_per_facet = TRUE) +
     xlim(-3, 3) + 
     ylim(-3,3) + 
     scale_color_manual(values = c("tomato", "steelblue")) +
     labs(x = expression(mu[1]), y = expression(mu[2])) +
     theme_bw(base_size = 20))
p_traj_latent
ggsave("./images2/plot_traj_latent_2-6-19_new.png", p_traj_latent,
       width = 14, height = 12, dpi = 300)
