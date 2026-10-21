#===============================================================================
# This script does the following
# - Build quarterly graph from imported proximity array
# - Fit a colSBM model on the data
# - Export the fitted model
#===============================================================================


library(tidyverse)
library(lubridate)
# install.packages("devtools")
#devtools::install_github("Chabert-Liddell/colSBM")
library(colSBM)

threshold <-  0.8
quantile   <- 0.8


bc_jour_cvae_ar <- readRDS("Rdata/bc_jour_cvae.rds")
ldates <- readRDS("Rdata/ldates.rds")



trim_year <- expand.grid(seq(4), seq(2014, 2018)) %>%
  rename(Trim = Var1, Year = Var2)

#===============================================================================
# Build a list of vector of boolean. TRUE if vessels [i] has fishing day in 
# quarter [[k]]
#===============================================================================
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
           ! is.na(diag(net))
         })
saveRDS(trim_peche, "./Rdata/trim_peche.rds")


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



if (! file.exists("Rdata/save_sbmfit/iidcolsbm_trim.rds")) {
  iidcol_trim_clust <- colSBM::clusterize_unipartite_networks(netlist = net_BCq8_trim,
                                                   colsbm_model = "iid",
                                                   net_id =  trim_year %>% transmute(trimyear = paste0(Trim, "_", Year)) %>% pull(trimyear),
                                                   nb_run = 5,
                                                   directed = FALSE)
  saveRDS(iidcol_trim_clust, "Rdata/save_sbmfit/iidcolsbm_trim.rds")
  rm(iidcol_trim_clust)
  gc()
}

