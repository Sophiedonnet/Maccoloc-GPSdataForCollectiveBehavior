#===============================================================================
# This script does the following
# - Import the reorder and filted latent representation dataset of vae4 and 
# cvae3
# - Import the trajectory dataset of XGV
# - Standardized the geographical coordinates
# - Selects 3 vessels
# - Plot the latent space, trajectory and times series of those 3 vessels
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
df_traj <- read_csv(file = "./Rdata/xgv_100_14-18.csv")
df_traj <- df_traj %>% filter(wday(time) %in% c(2,3,4,5,6)) |>
  filter(maree != 12097796)
df_cvae <- read_csv(file = "./Rdata/xgv_cvae_3_reordered_filtered.csv")
df_vae <- read_csv(file = "./Rdata/xgv_vae_4_reordered_filtered.csv")

navire_id <- unique(df_cvae$navire)

b = 2
bp = 6
bpp = 19

#===============================================================================
# Plotting the latent distribution time series
#===============================================================================
df_vae %>% filter(navire == navire_id[b]) %>% 
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




#===============================================================================
# Plotting the latent mean and the standardized trajectories
#===============================================================================

df_traj <- df_traj |> 
  mutate(x_std = coord2std(x, -4.28),
        y_std = coord2std(y, 47.78)) 

boat <- c(2, 6, 19)
(df_traj   %>% 
    filter(navire %in% navire_id[boat]) %>% 
    mutate(maree = as_factor(maree),
           #navire = as_factor(navire),
           month = as_factor(month(time)),
           year = year(time)) %>% 
    mutate(h = rep(seq(24), nrow(.)/24)) %>% 
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
     theme_bw(base_size = 20))