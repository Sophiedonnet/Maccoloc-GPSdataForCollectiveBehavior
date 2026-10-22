#===============================================================================
# Ordering of the latent dimensions
# (Figure 4 of the paper, file plot_dimension.png)
#
# This script does the following:
# - Import the latent representations learned by the CVAE (d_z = 3) and the
#   VAE (d_z = 4): for each trip, mean mu_k and log-variance of the Gaussian
#   variational distribution on each dimension k
# - For 1e5 random couples of trips, compute the Bhattacharyya coefficient
#   (BC) on each dimension, and the Kullback-Leibler divergence (DKL) between
#   each variational distribution and the N(0,1) prior (Figure 4)
# - Order the dimensions by importance (from the figure: the less informative
#   a dimension, the closer its BC to 1 and its DKL to 0)
# - Export the latent representations with ordered dimensions into Rdata/
#
# Inputs (exported by py_workflow/cvae_analysis.ipynb):
# - pydata/latentspace_cvae_10k_3.csv, pydata/latentspace_vae_10k_4.csv
#
# Outputs (used by all the following scripts):
# - Rdata/xgv_cvae_3_reordered_filtered.csv, Rdata/xgv_vae_4_reordered_filtered.csv
#   columns mu_k (means) and sig_k (VARIANCES gamma^2_k, despite the name)
#
# Notes:
# - Figure 4 is not saved by the script (no ggsave).
# - The couples of trips are sampled without a seed: the figure changes
#   slightly from one run to the next.
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
# One trip per vessel and per day (duplicated days are removed)
df_cvae <- read_csv(file = "./pydata/latentspace_cvae_10k_3.csv")
df_cvae <- df_cvae %>% group_by(navire) %>% filter(! duplicated(time)) %>% 
  ungroup() 
df_vae <- read_csv(file = "./pydata/latentspace_vae_10k_4.csv")
df_vae <- df_vae %>% group_by(navire) %>% filter(! duplicated(time)) %>% 
  ungroup()


#===============================================================================
# Rename the mean and variance
#===============================================================================
# In the Python export, the means and the log-variances are both in columns
# named 0, 1, 2(, 3); read_csv makes the names unique by appending the column
# position (`0...7`, ...). The means come first.
df_vae <- df_vae %>% 
  rename(mu_1 = `0...7`, mu_2 = `1...8`, mu_3 = `2...9`, mu_4 = `3...10`,
         logsig_1 = `0...11`, logsig_2 = `1...12`, logsig_3 = `2...13`, logsig_4 = `3...14`)
df_cvae <- df_cvae %>% rename(mu_1 = `0...7`, mu_2 = `1...8`, mu_3 = `2...9`,
                              logsig_1 = `0...10`, logsig_2 = `1...11`, logsig_3 = `2...12`)
#===============================================================================
# Sample 1e5 couple of trajectories
#===============================================================================
# 2e5 couples are drawn, couples of a trip with itself are removed and the
# first 1e5 are kept
n_traj <- nrow(df_cvae)
ord1 <- sample(n_traj, 2e5, replace = TRUE)
ord2 <- sample(n_traj, 2e5, replace = TRUE)
ind <- ord1 != ord2
ord1 <- ord1[ind][1:1e5]
ord2 <- ord2[ind][1:1e5]

# Means and variances (exp of the log-variances) of each dimension
mu_cvae <- as.matrix(df_cvae[,7:9])
sig2_cvae <- as.matrix(exp(df_cvae[,7:9+3]))
mu_vae <- as.matrix(df_vae[,1:4+6])
sig2_vae <- as.matrix(exp(df_vae[,1:4+10]))


#===============================================================================
# PLots of the reorder dimensions (Figure 4)
#===============================================================================
# Left: BC between the two trips of each couple, dimension by dimension.
# Right: DKL between the variational distribution of each trip and N(0,1).
# The columns are renamed with the rank of each dimension in the final order:
# - VAE : original dimensions 4, 3, 1, 2 become 1, 2, 3, 4
# - CVAE: original dimensions 1, 3, 2 become 1, 2, 3

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
# Same order as in the figure above. The log-variances are replaced by the
# variances sig_k = exp(logsig_k).

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
