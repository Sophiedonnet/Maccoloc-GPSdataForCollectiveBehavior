#===============================================================================
# colSBM fit on the quarterly proximity graphs
# (fit used for all the results of the paper: Figures 9-11, A2, A3)
#
# This script does the following
# - Build the 20 quarterly proximity graphs from the imported proximity array:
#   vessels b and b' are connected in a quarter if the `quantile`-quantile of
#   their daily BC over the days of the quarter is above `threshold`
# - Fit an iid-colSBM on the collection of the 20 graphs
# - Export the fitted model
#
# Inputs (produced by 06_build_proximity_graph.R):
# - Rdata/bc_jour_cvae.rds : array (days x 33 x 33) of the BC between the
#                            trajectories of two vessels on the same day
# - Rdata/ldates.rds       : dates corresponding to the first dimension
#
# Outputs:
# - Rdata/trim_peche.rds   : for each quarter, vessels (1..33) that fished
# - Rdata/save_sbmfit/iidcolsbm_trim_round8.rds : colSBM fit (object of class
#   bmpop; the selected model is $best_fit, with 6 blocks)
#
# Every step is skipped when its output file already exists.
#
# Note: the colSBM options below are the ones stored in the saved fit
# ($global_opts). The number of runs (nb_run) is not stored in the fit.
#
# Earlier fits (first version of the paper, intermediate rounds) and the
# script comparing them with this fit are in archive/ (see archive/README.md).
#===============================================================================


library(tidyverse)
library(lubridate)
#devtools::install_github("GrossSBM/colSBM",force=TRUE)
library(colSBM)

threshold <-  0.8
quantile   <- 0.8


bc_jour_cvae_ar <- readRDS("Rdata/bc_jour_cvae.rds")
ldates <- readRDS("Rdata/ldates.rds")
# The 20 quarters, ordered by year then quarter: (1, 2014), (2, 2014), ...
trim_year <- expand.grid(seq(4), seq(2014, 2018)) %>%
  rename(Trim = Var1, Year = Var2)

#===============================================================================
# Build a list of vector of boolean. TRUE if vessels [i] has fishing days in
# quarter [[k]]
#===============================================================================
file_trim_peche = "./Rdata/trim_peche.rds"
if (! file.exists(file_trim_peche)){
trim_peche <-
   lapply(seq(20),
          function(t) {
            id <- which(
              year(ldates) == trim_year$Year[t] &
              (((month(ldates) + 3) %%4) + 1) == trim_year$Trim[t])
            net <- 1*(
              apply(
                bc_jour_cvae_ar[id,,],
                c(2,3),
                "quantile",
                na.rm = TRUE, probs = quantile) > threshold)
            trim_peche_t <- ! is.na(diag(net))
            return(trim_peche_t)
          })
saveRDS(trim_peche, file_trim_peche)
}

#===============================================================================
# Construct collection of full network
#===============================================================================
# Adjacency matrix of each quarter, restricted to the vessels that fished
# during the quarter (rows/columns named by the vessel number 1..33). NA if
# two vessels never fished the same day.
net_BCq8_trim <-
  lapply(seq(20),
         function(t) {
           id <- which(
             year(ldates) == trim_year$Year[t] &
               (((month(ldates) + 3) %%4) + 1) == trim_year$Trim[t])

           net <- 1*(
             apply(
               bc_jour_cvae_ar[id,,],
               c(2,3),
               "quantile",
               na.rm = TRUE, probs = quantile) > threshold)
           rownames(net) <- seq(33)
           colnames(net) <- seq(33)
           net <- net[! is.na(diag(net)), ! is.na(diag(net))]
           diag(net) <- 0
           net
         })

#===============================================================================
# Fit the iid-colSBM
#===============================================================================
my_global_opts = list(nb_cores = 1,
                      nb_models = 3,
                      nb_init = 10,
                      depth = 4,
                      verbosity = 1,
                      spectral_init = TRUE,
                      Q_max = 7,
                      plot_details = 1)

file_res_full = paste0("Rdata/save_sbmfit/iidcolsbm_trim_round8.rds")
if (! file.exists(file_res_full)){
  iidcol_trim_clust <- colSBM::estimate_colSBM(netlist = net_BCq8_trim,
                                             colsbm_model = "iid",
                                             net_id =  trim_year %>% transmute(trimyear = paste0(Trim, "_", Year)) %>% pull(trimyear),
                                             nb_run = 20,
                                             directed = FALSE, global_opts=my_global_opts)
  saveRDS(iidcol_trim_clust, file_res_full)
  rm(iidcol_trim_clust)
}

# Selected model
bestfit <- readRDS(file_res_full)$best_fit
bestfit$nb_clusters
