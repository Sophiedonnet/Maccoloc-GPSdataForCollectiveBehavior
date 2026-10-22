#===============================================================================
# This script does the following
# - Build quarterly graph from a fraction of all the imported proximity array 
# - Fit a colSBM model on the data for each sub sample
# - Export the fitted model for each sub sample
#===============================================================================


library(tidyverse)
library(lubridate)
library(mclust)
#devtools::install_github("GrossSBM/colSBM",force=TRUE)
library(colSBM)
library(patchwork)
source(file='99_helper_functions.R')
source(file='98_functions_resampling.R')

threshold <-  0.8
quantile   <- 0.8


bc_jour_cvae_ar <- readRDS("Rdata/bc_jour_cvae.rds")
ldates <- readRDS("Rdata/ldates.rds")
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
# Construct collection of full network
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
# Replicate saint clair'e results 
#===============================================================================
my_global_opts = list(nb_cores = 4,
                      nb_models = 1,
                      nb_init = 10,
                      depth = 4,
                      verbosity = 1,
                      spectral_init = TRUE,
                      Q_max = 7,
                      plot_details = 1)


#-------------- de Zinit 
iidcol_trim_clust_full_replicate7 <- readRDS(paste0("Rdata/save_sbmfit/iidcolsbm_trim_round7.rds"))
iidcol_trim_clust_full_saintClair <- readRDS(paste0("Rdata/save_sbmfit/iidcolsbm_trim.rds"))



Z_init <- lapply(1:7,function(Q){
    iidcol_trim_clust_full_replicate7$model_list[[1]][[Q]][[1]]$Z})
Z_init[[5]] <- iidcol_trim_clust_full_saintClair$Z
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

  

#===============================================================================
#  Exploitation résultats. 
#===============================================================================
file_res_full_replicate8 = paste0("Rdata/save_sbmfit/iidcolsbm_trim_round8.rds")
file_res_full_saintClair = paste0("Rdata/save_sbmfit/iidcolsbm_trim.rds")
iidcol_trim_clust_full_saintClair <- readRDS(file_res_full_saintClair)
iidcol_trim_clust_full_replicate8 <- readRDS(file_res_full_replicate8)

bestfit_replicate <- iidcol_trim_clust_full_replicate8$best_fit
bestfit_saintClair <- iidcol_trim_clust_full_saintClair

#file_res_full_replicate7 = paste0("Rdata/save_sbmfit/iidcolsbm_trim_round7.rds")
#file_res_full_replicate6= paste0("Rdata/save_sbmfit/iidcolsbm_trim_round6.rds")
#file_res_full_replicate5= paste0("Rdata/save_sbmfit/iidcolsbm_trim_round5.rds")












#==================================================== 
#-------------- Reorder by minimizing the graphon L2 distance
#====================================================

o_saintClair = get_order(bestfit_saintClair)
o_replicate  = get_order_withref(fit = bestfit_replicate, fit_ref=bestfit_saintClair,order_ref=o_saintClair)

p_replicate <- plot(bestfit_replicate, type = "graphon",ord=o_replicate)  + theme(legend.position = "none") + ggtitle("Replicate")
p_saintClair <- plot(bestfit_saintClair, type = "graphon",ord=o_saintClair) + theme(legend.position = "none")  +  ggtitle("Saint Clair")

(p_saintClair | p_replicate ) 


#==================================================== 
#--------------------- compar final clustering bateau 
#==================================================== 

Ztrim_replicate <- get_Ztrim_relabel(bestfit_replicate)
Ztrim_replicate_2  <- relabel_Ztrim(Ztrim_replicate,o_replicate_withref)
Ztrim_map_replicate <- get_Ztrim_map(Ztrim_replicate_2)

Ztrim_saintClair <- get_Ztrim(bestfit_saintClair)
Ztrim_saintClair_2  <- relabel_Ztrim(Ztrim_saintClair,o_saintClair)
Ztrim_map_saintClair <- get_Ztrim_map(Ztrim_saintClair_2)

table(Ztrim_map_replicate ,Ztrim_map_saintClair)
adjustedRandIndex(Ztrim_map_replicate ,Ztrim_map_saintClair)




