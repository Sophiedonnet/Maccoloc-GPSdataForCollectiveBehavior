#===============================================================================
# Bootstrap of the collective-behaviour analysis
# (Appendix "Uncertainty quantification" of the paper)
#
# This script does the following:
# 1. Build the 20 quarterly proximity graphs from the full data (reference)
# 2. Generate N_boot bootstrap collections: for each quarter, keep a random
#    fraction `frac` of the fishing days (sampling without replacement) and
#    rebuild the proximity graph from these days only
# 3. Fit an iid-colSBM on each bootstrap collection
# 4. Check that the bootstrap collections (and fitted networks) are distinct
# 5. Compare the bootstrap fits with the full-data fit:
#    - distribution of the estimated number of blocks R_hat
#    - mesoscale structure of one bootstrap sample per value of R_hat
#    - ARI between bootstrap and full-data clusterings, per quarter
#
# Inputs (produced by 06_build_proximity_graph.R and 07b_...R):
# - Rdata/bc_jour_cvae.rds : array (days x 33 x 33) of the Bhattacharyya
#                            coefficients between the trajectories of two
#                            vessels on the same day (NA if a vessel did not
#                            fish that day)
# - Rdata/ldates.rds       : dates corresponding to the first dimension
# - Rdata/save_sbmfit/iidcolsbm_trim_round8.rds : colSBM fit on the full data
#
# Outputs (in `dir_boot`, one file per bootstrap replicate k):
# - subnets_<frac>/sub_net_<k>.rds        : list(days, netlist) of replicate k
# - fits_<frac>/iidcolsbm_trim_sub_net_<k>.rds : colSBM fit of replicate k
#
# Every step is skipped when its output file already exists, so the script
# can be interrupted and re-run. Delete the files of a step to recompute it.
#
# Computing time: one colSBM fit takes from several minutes to more than an
# hour depending on `nb_run`, `my_global_opts` and the number of cores.
#
# Note on the results of the first version of the paper: the files in
# Rdata/save_colsbmfit_subsample/subnets_90/sub_net_1..34.Rdata contain an
# object named `net_BCq8_trim` (the full-data networks) instead of
# `sub_net_BCq8_trim`. The former estimation loop used `load()` and then
# `sub_net_BCq8_trim`, i.e. an object left in the R session: the fits
# res_iidcolsbm_subnets_90/iidcolsbm_trim_sub_net_1..30.rds were all
# estimated on the same network collection. This script avoids the problem by
# saving and reading each replicate with saveRDS()/readRDS() and by checking
# (step 4) that the replicates differ.
#===============================================================================

rm(list = ls())
library(tidyverse)
library(lubridate)
library(mclust)    # adjustedRandIndex()
library(patchwork)
library(colSBM)    # devtools::install_github("GrossSBM/colSBM")
source(file = "99_helper_functions.R")
source(file = "98_functions_resampling.R") # get_order*(), get_Ztrim(); needs graphclust


#===============================================================================
# 0. Parameters
#===============================================================================

# Construction of the proximity graphs (same values as for the full data,
# see 07_fit_colsbm_proxy_graphs.R): vessels b and b' are connected in a
# quarter if the `quantile`-quantile of their daily BC is above `threshold`
threshold <- 0.8
quantile  <- 0.8

# Bootstrap
frac      <- 0.9   # fraction of the fishing days kept in each quarter
N_boot    <- 30    # number of bootstrap replicates
seed_boot <- 1234  # replicate k uses seed_boot + k (sub-sampling and fit)

# colSBM options: same as for the reference fit on the full data
# (Rdata/save_sbmfit/iidcolsbm_trim_round8.rds)
nb_run <- 8
my_global_opts <- list(nb_cores = 4,
                       nb_models = 5,
                       nb_init = 10,
                       depth = 4,
                       verbosity = 1,
                       spectral_init = TRUE,
                       Q_max = 7,
                       plot_details = 1)

# Output folders
dir_boot    <- "Rdata/save_colsbm_bootstrap"
dir_subnets <- file.path(dir_boot, paste0("subnets_", 100 * frac))
dir_fits    <- file.path(dir_boot, paste0("fits_", 100 * frac))
dir.create(dir_subnets, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_fits, recursive = TRUE, showWarnings = FALSE)

file_subnet <- function(k) file.path(dir_subnets, paste0("sub_net_", k, ".rds"))
file_fit    <- function(k) file.path(dir_fits, paste0("iidcolsbm_trim_sub_net_", k, ".rds"))


#===============================================================================
# 1. Data and construction of the quarterly proximity graphs
#===============================================================================

bc_jour_cvae_ar <- readRDS("Rdata/bc_jour_cvae.rds")
ldates <- readRDS("Rdata/ldates.rds")

# The 20 quarters, ordered by year then quarter: (1, 2014), (2, 2014), ...
trim_year <- expand.grid(seq(4), seq(2014, 2018)) %>%
  rename(Trim = Var1, Year = Var2)
net_id <- paste0(trim_year$Trim, "_", trim_year$Year)

# Indices (in ldates) of the days of quarter t
days_of_quarter <- lapply(seq(20), function(t) {
  which(year(ldates) == trim_year$Year[t] &
          (((month(ldates) + 3) %% 4) + 1) == trim_year$Trim[t])
})

# Proximity graph built from the days `days`:
# - A[b, b'] = 1 if the `quantile`-quantile of BC(b, b') over the days is
#   larger than `threshold`. NA if b and b' never fished the same day.
# - vessels that did not fish during these days (NA diagonal) are removed,
#   rows/columns are named by the vessel number (1..33)
build_network <- function(days) {
  net <- 1 * (apply(bc_jour_cvae_ar[days, , ], c(2, 3), "quantile",
                    na.rm = TRUE, probs = quantile) > threshold)
  rownames(net) <- colnames(net) <- seq(33)
  net <- net[!is.na(diag(net)), !is.na(diag(net))]
  diag(net) <- 0
  net
}


#===============================================================================
# 2. Reference: networks of the full data and colSBM fit
#===============================================================================

net_BCq8_trim <- lapply(days_of_quarter, build_network)

bestfit_full <- readRDS("Rdata/save_sbmfit/iidcolsbm_trim_round8.rds")$best_fit
# Sanity check: the reference fit was estimated on these networks
stopifnot(isTRUE(all.equal(
  lapply(net_BCq8_trim, unname),
  lapply(bestfit_full$A, function(A) { A <- unname(A); A[A == -1] <- NA; A }))))


#===============================================================================
# 3. Generate the bootstrap collections of networks
#===============================================================================
# For replicate k and each quarter, floor(frac * n_days) days are sampled
# without replacement among the n_days fishing days of the quarter. Removing
# days removes opportunities for two vessels to fish the same day, hence
# more NA in the bootstrap networks.

for (k in seq(N_boot)) {
  if (!file.exists(file_subnet(k))) {
    set.seed(seed_boot + k)
    days <- lapply(days_of_quarter, function(id) sort(sample(id, floor(frac * length(id)))))
    netlist <- lapply(days, build_network)
    saveRDS(list(days = days, netlist = netlist), file_subnet(k))
  }
}


#===============================================================================
# 4. Fit an iid-colSBM on each bootstrap collection
#===============================================================================

for (k in seq(N_boot)) {
  if (!file.exists(file_fit(k))) {
    message("Bootstrap replicate ", k, "/", N_boot)
    sub_net <- readRDS(file_subnet(k))$netlist
    set.seed(seed_boot + k)
    fit_k <- colSBM::estimate_colSBM(netlist = sub_net,
                                     colsbm_model = "iid",
                                     net_id = net_id,
                                     nb_run = nb_run,
                                     directed = FALSE,
                                     global_opts = my_global_opts)
    saveRDS(fit_k, file_fit(k))
    rm(fit_k)
    gc()
  }
}


#===============================================================================
# 5. Checks
#===============================================================================

boot_fits <- lapply(seq(N_boot), function(k) readRDS(file_fit(k))$best_fit)

# Each fit was estimated on the networks of its own replicate
fit_networks <- lapply(boot_fits, function(f) {
  lapply(f$A, function(A) { A <- unname(A); A[A == -1] <- NA; A })
})
for (k in seq(N_boot)) {
  stopifnot(isTRUE(all.equal(lapply(readRDS(file_subnet(k))$netlist, unname),
                             fit_networks[[k]])))
}
# The replicates are all different (and different from the full data)
stopifnot(length(unique(fit_networks)) == N_boot)
stopifnot(!any(sapply(fit_networks, identical, lapply(net_BCq8_trim, unname))))


#===============================================================================
# 6. Results
#===============================================================================

#------------- Number of blocks selected on each bootstrap replicate
R_hat <- sapply(boot_fits, function(f) f$nb_clusters)
table(R_hat) / N_boot * 100


#------------- Mesoscale structures: one replicate per value of R_hat
# The blocks are ordered as for the full-data fit (get_order), and the blocks
# of each replicate are permuted to minimise the L2 distance between its
# graphon and the one of the full-data fit.
o_full <- get_order(bestfit_full)

plot_boot <- function(m) {
  f <- boot_fits[[m]]
  o <- get_order_withref(fit = f, fit_ref = bestfit_full, order_ref = o_full)
  plot(f, type = "graphon", ord = o) +
    ggtitle(bquote(m == .(m) ~ "," ~ hat(R) == .(f$nb_clusters))) +
    theme(plot.title = element_text(size = 12), legend.position = "none")
}
# first replicate for each value of R_hat, most frequent value first
R_freq <- sort(table(R_hat), decreasing = TRUE)
m_show <- sapply(as.integer(names(R_freq)), function(r) which(R_hat == r)[1])

m_show=c(3,4)

#------------- ARI between bootstrap and full-data clusterings
# For replicate k and quarter t, the ARI is computed on the vessels present
# in both the full-data and the bootstrap networks of quarter t.
Z_trim_full <- get_Ztrim(bestfit_full) # 33 x 20 matrix of blocks (NA = absent)

ARI <- sapply(boot_fits, function(f) {
  Z_trim <- get_Ztrim(f)
  sapply(seq(20), function(t) {
    keep <- !is.na(Z_trim[, t]) & !is.na(Z_trim_full[, t])
    adjustedRandIndex(Z_trim_full[keep, t], Z_trim[keep, t])
  })
})
summary(c(ARI))

p_hist_ARI <- ggplot(tibble(ARI = c(ARI)), aes(x = ARI)) +
  # boundary = 1: the last bin ends at 1, so that ARI = 1 is not dropped
  geom_histogram(aes(y = after_stat(count / sum(count))), binwidth = 0.03,
                 boundary = 1, fill = "#4C72B0", color = "white") +
  geom_vline(xintercept = mean(ARI, na.rm = TRUE), color = "#DD8452",
             linewidth = 1.1) +
  annotate("text", x = mean(ARI, na.rm = TRUE) - 0.02, y = Inf, hjust = 1,
           vjust = 2, label = "Mean of the ARI", color = "#DD8452") +
  coord_cartesian(xlim = c(0, 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) + # room for the label
  ylab("Frequency") +
  theme_minimal()


#------------- Figure of the appendix (bootstrap_results.png)
p_struct <- lapply(m_show, plot_boot)
final_plot <- wrap_plots(p_struct, nrow = 1) /
  (plot_spacer() + p_hist_ARI + plot_spacer() + plot_layout(widths = c(1, 2, 1)))
final_plot
ggsave(paste0("./images2/bootstrap_results_", 100 * frac, ".png"),
       final_plot, width = 9, height = 8, dpi = 300)
