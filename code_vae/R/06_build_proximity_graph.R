#===============================================================================
# This script does the following
# - Import the reorder and filted latent representation dataset of cvae3
# - Compute BC between trajectories from the same day
# - Compute BC between trajectories from adjacent days
# - Export an array of dim TxBxB and an array of dim (T-1)xBxB
# - Export ldates, a vector of dates to make a link between the 
# dimension of the matrices and the fishing day
#===============================================================================


library(tidyverse)
source("./99_helper_functions.R")

df_cvae <- read_csv(file = "./Rdata/xgv_cvae_3_reordered_filtered.csv")


#===============================================================================
# Reshaping the dataframe to have one maree each day for each vessel
# Extract the variational parameters as an array
#===============================================================================


df_cvae_c <- tibble(time = rep(seq(min(df_cvae$time), max(df_cvae$time), by = "days"), times = 33),
                    navire = rep(unique(df_cvae$navire), each = length(time)/33)) %>% 
  filter(wday(time) %in% c(2,3,4,5,6)) %>% 
  full_join(df_cvae, by = c("time", "navire")) %>% 
  arrange(time)


ar_cvae_c <- df_cvae_c %>% select(-`...1`) %>% 
  pivot_longer(cols = starts_with("mu") | starts_with("sig"), names_to = "param", values_to = "mu_sigma2") %>%
  select("mu_sigma2") %>% as_vector() %>% 
  array(dim = c(6, 33, length(.)/(33*6))) %>% aperm(c(1,3,2))

ar_cvae_mu <-   ar_cvae_c[1:3,,]
ar_cvae_sig2 <-   ar_cvae_c[4:6,,]




#===============================================================================
# Compute all the possible BC form same day trajectories
# Reshape them as an array where the first dimension is the time and 
# Dim2 and 3 is a matrix between vessels
#===============================================================================

d <- expand.grid(seq(33), seq(33))
bc_jour_cvae <- lapply(seq(nrow(d)),
                       function(i) {
                         apply(
                           X = bhattacharyya_1d(
                             mu1 = ar_cvae_mu[,,d[i,1]], 
                             mu2 = ar_cvae_mu[,,d[i,2]], 
                             var1 = ar_cvae_sig2[,,d[i,1]], 
                             var2 = ar_cvae_sig2[,,d[i,2]]),
                           MARGIN = 2, FUN = "prod")
                       })
bc_jour_cvae_ar <- aperm(array(unlist(bc_jour_cvae), 
                               dim = c(length(bc_jour_cvae[[1]]), 33, 33)), 
                         perm = c(1,2,3))

saveRDS(bc_jour_cvae_ar, "./Rdata/bc_jour_cvae.rds")


#===============================================================================
# Compute all the possible BC form same day trajectories lagged by one day
# Reshape them as an array where the first dimension is the time and 
# Dim2 and 3 is a matrix between vessels
#===============================================================================


bc_lag_cvae <- lapply(seq(nrow(d)),
                      function(i) {
                        apply(
                          X = bhattacharyya_1d(
                            mu1 = ar_cvae_mu[,1:(dim(ar_cvae_mu)[2]-1),d[i,1]], 
                            mu2 = ar_cvae_mu[,2:dim(ar_cvae_mu)[2],d[i,2]], 
                            var1 = ar_cvae_sig2[,1:(dim(ar_cvae_mu)[2]-1),d[i,1]], 
                            var2 = ar_cvae_sig2[,2:dim(ar_cvae_mu)[2],d[i,2]]),
                          MARGIN = 2, FUN = "prod")
                      })
bc_lag_cvae_ar <- aperm(array(unlist(bc_lag_cvae), 
                              dim = c(length(bc_lag_cvae[[1]]), 33, 33)), 
                        perm = c(1,2,3))

saveRDS(bc_lag_cvae_ar, "./Rdata/bc_lag_cvae.rds")

ldates <- df_cvae_c %>% pull(time) %>% unique()

saveRDS(ldates, "./Rdata/ldates.rds")
