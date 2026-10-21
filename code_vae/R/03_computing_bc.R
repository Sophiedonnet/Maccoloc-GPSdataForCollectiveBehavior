#===============================================================================
# This script does the following:
# - Import the latent representations learned by the cvae/vae (3 and 4 dims)
# - Import the xgv dataset
# - Filter the week-end and stange maree trajectories
# - Compute the BC by dimensions
# - Order the dimension by importance
# - Export the concatenated and ordered dataset in pydata into Rdata 
# 
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
# df_traj_std <- read_csv(file = "./Rdata/xgv_traj_std.csv")
# df_traj_std <- df_traj_std[df_traj_std$...1 != 5390,]
# df_traj <- read_csv(file = "./Rdata/xgv_100_14-18.csv")
# df_traj <- df_traj %>% filter(wday(time) %in% c(2,3,4,5,6)) |> 
#   filter(maree != 12097796)
df_cvae <- read_csv(file = "./pydata/latentspace_cvae_10k_3.csv")
df_cvae <- df_cvae %>% group_by(navire) %>% filter(! duplicated(time)) %>% 
  ungroup() 
df_vae <- read_csv(file = "./pydata/latentspace_vae_10k_4.csv")
df_vae <- df_vae %>% group_by(navire) %>% filter(! duplicated(time)) %>% 
  ungroup()


#===============================================================================
# Rename the mean and variance
#===============================================================================
df_vae <- df_vae %>% 
  rename(mu_1 = `0...7`, mu_2 = `1...8`, mu_3 = `2...9`, mu_4 = `3...10`,
         logsig_1 = `0...11`, logsig_2 = `1...12`, logsig_3 = `2...13`, logsig_4 = `3...14`)
df_cvae <- df_cvae %>% rename(mu_1 = `0...7`, mu_2 = `1...8`, mu_3 = `2...9`,
                              logsig_1 = `0...10`, logsig_2 = `1...11`, logsig_3 = `2...12`)
#===============================================================================
# Sample 1e5 couple of trajectories 
#===============================================================================
n_traj <- nrow(df_cvae)
ord1 <- sample(n_traj, 2e5, replace = TRUE)
ord2 <- sample(n_traj, 2e5, replace = TRUE)
ind <- ord1 != ord2
ord1 <- ord1[ind][1:1e5]
ord2 <- ord2[ind][1:1e5]

mu_cvae <- as.matrix(df_cvae[,7:9])
sig2_cvae <- as.matrix(exp(df_cvae[,7:9+3]))
mu_vae <- as.matrix(df_vae[,1:4+6])
sig2_vae <- as.matrix(exp(df_vae[,1:4+10]))


#===============================================================================
# PLots of the reorder dimensions
#===============================================================================


(bhattacharyya_1d(mu_vae[ord1,], mu_vae[ord2,], 
                  sig2_vae[ord1,], sig2_vae[ord2,]) |>
    as_tibble() %>% 
    rename("3" = logsig_1, "4" = logsig_2, "2" = logsig_3, "1" = logsig_4) %>% 
    pivot_longer(cols = everything(), names_to = "Dim", values_to = "BC") %>% 
    mutate(Model = "VAE") %>% 
    bind_rows(bhattacharyya_1d(mu_cvae[ord1,], mu_cvae[ord2,], 
                               sig2_cvae[ord1,], sig2_cvae[ord2,]) |>
                as_tibble() %>% 
                rename("1" = logsig_1, "3" = logsig_2, "2" = logsig_3) %>% 
                pivot_longer(cols = everything(), names_to = "Dim", values_to = "BC") %>% 
                mutate(Model = "CVAE")) %>% 
    ggplot(aes(y = BC, x = Dim)) +
    geom_boxplot() +
    geom_violin(col = "tomato", alpha = .6) +
    facet_wrap(~ Model, scale = "free_x") +
    xlab("Dimension") +
    theme_bw(base_size = 20)) +
  (kullback_1d(mu_vae, 0,  sig2_vae, 1) %>%   as_tibble() %>% 
     rename("3" = logsig_1, "4" = logsig_2, "2" = logsig_3, "1" = logsig_4) %>% 
     pivot_longer(cols = everything(), names_to = "Dim", values_to = "DKL") %>% 
     mutate(Model = "VAE") %>% 
     bind_rows(kullback_1d(mu_cvae, 0, sig2_cvae, 1) |>
                 as_tibble() %>% 
                 rename("1" = logsig_1, "3" = logsig_2, "2" = logsig_3) %>% 
                 pivot_longer(cols = everything(), names_to = "Dim", values_to = "DKL") %>% 
                 mutate(Model = "CVAE")) %>% 
     ggplot(aes(y = DKL, x = Dim)) +
     geom_boxplot() +
     geom_violin(col = "tomato", alpha = .6) +
     #  geom_jitter() +
     facet_wrap(~ Model, scale = "free_x") +
     xlab("Dimension") +
     theme_bw(base_size = 20))


#===============================================================================
# Rename and reorder the dimensions
#===============================================================================


df_vae <- df_vae %>% rename(mu_1 = mu_4, mu_2 = mu_3, mu_3 = mu_1, mu_4 = mu_2, 
                            sig_1 = logsig_4, sig_2 = logsig_3, 
                            sig_3 = logsig_1, sig_4 = logsig_2) %>%
  mutate(sig_1 = exp(sig_1), sig_2 = exp(sig_2), 
         sig_3 = exp(sig_3), sig_4 = exp(sig_4))
df_cvae <- df_cvae %>% rename(mu_2 = mu_3,  mu_3 = mu_2, sig_1 = logsig_1,
                              sig_2 = logsig_3, sig_3 = logsig_2) %>% 
  mutate(sig_1 = exp(sig_1), sig_2 = exp(sig_2), sig_3 = exp(sig_3))

write_csv(df_cvae, file = "Rdata/xgv_cvae_3_reordered_filtered.csv")
write_csv(df_vae, file = "Rdata/xgv_vae_4_reordered_filtered.csv")
