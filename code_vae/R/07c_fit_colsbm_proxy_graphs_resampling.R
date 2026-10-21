#===============================================================================
# This script does the following
# - Build quarterly graph from a fraction of all the imported proximity array 
# - Fit a colSBM model on the data for each sub sample
# - Export the fitted model for each sub sample
#===============================================================================
rm(list=ls())
library(tidyverse)
library(lubridate)
library(mclust)
library(patchwork)
library(colSBM)
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
# Create bootstrap data  
#===============================================================================

N_sub_sampl <- 30

#---------- Create subnet 
# for (k in 1:(N_sub_sampl)) {
#   print(k)
#   file_sim=paste0("Rdata/save_colsbmfit_subsample/subnets_90/sub_net_",k,".Rdata")
#   sub_net_BCq8_trim  <- list()
#   
#   if (! file.exists(file_sim)) {
#     for (t in 1:20){
#       id.t <- which(year(ldates) == trim_year$Year[t] &
#                  (((month(ldates) + 3) %%4) + 1) == trim_year$Trim[t])
#       l_id.t <- length(id.t)
#       sub_id.t <- sample(id.t,floor(l_id.t*0.9))
#       sub_net.t <- 1*(apply(
#                  bc_jour_cvae_ar[sub_id.t,,],
#                  c(2,3),
#                  "quantile",
#                  na.rm = TRUE, probs = quantile) > threshold)
#       #trim_peche.t <- ! is.na(diag(sub_net.t))
#       rownames(sub_net.t) <- seq(33)
#       colnames(sub_net.t) <- seq(33)
#       sub_net.t <-sub_net.t[! is.na(diag(sub_net.t)), ! is.na(diag(sub_net.t))]
#       diag(sub_net.t) <- 0
#       sub_net_BCq8_trim[[t]] <- sub_net.t
#       #sub_trim_peche[[t]] <- trim_peche.t 
#       }
#    save(sub_net_BCq8_trim,file=file_sim)
#   }
#  
# }
#===============================================================================
#---------- estim iid col sbm on bootstrap data
#===============================================================================

# my_global_opts = list(nb_cores = 4,
#                       nb_models = 5,
#                       nb_init = 10,
#                       depth = 4,
#                       verbosity = 1,
#                       spectral_init = TRUE,
#                       Q_max = 7,
#                       plot_details = 1)
#                   
# for (k in 1:N_sub_sampl){
#    
#    print(paste("sub_sampl=",k))
#    load(file=paste0("Rdata/save_colsbmfit_subsample/subnets_90/sub_net_",k,".Rdata"))
#    file_res = paste0("Rdata/save_colsbmfit_subsample/res_iidcolsbm_subnets_90/iidcolsbm_trim_sub_net_",k,".rds")
# #   
# # 
#    if (! file.exists(file_res)) {
#      iidcol_trim_clust <- colSBM::estimate_colSBM(netlist = sub_net_BCq8_trim,
#                                                                  colsbm_model = "iid",
#                                                                  net_id =  trim_year %>% transmute(trimyear = paste0(Trim, "_", Year)) %>% pull(trimyear),
#                                                                  nb_run = 8,
#                                                                  directed = FALSE, global_opts=my_global_opts)
#      saveRDS(iidcol_trim_clust, file_res)
#      rm(iidcol_trim_clust)
#      gc()  
#    }
# # 
# }


#===============================================================================
#  Exploitation résultats. 
#===============================================================================

#================ Kestim
Kestim = rep(0,N_sub_sampl)
for (k in 1:N_sub_sampl){
  print(k)
  res_k <- readRDS(paste0("Rdata/save_colsbmfit_subsample/res_iidcolsbm_subnets_90/iidcolsbm_trim_sub_net_",k,".rds"))
  Kestim[k] <- res_k$best_fit$nb_clusters
  }
table(Kestim)/30*100



#================ Plot resulting structure

#=========  Full data estimation: used as a reference 
file_res_full = paste0("Rdata/save_sbmfit/iidcolsbm_trim_round8.rds")
iidcol_trim_clust_full <- readRDS(file_res_full)
bestfit_full <- iidcol_trim_clust_full$best_fit


file_res_full_saintClair = paste0("Rdata/save_sbmfit/iidcolsbm_trim.rds")
iidcol_trim_clust_full_saintClair <- readRDS(file_res_full_saintClair)
bestfit_full_saintClair <- iidcol_trim_clust_full_saintClair

o_saintClair = get_order(bestfit_full_saintClair )
o_full  = get_order_withref(fit = bestfit_full, fit_ref=bestfit_full_saintClair,order_ref=o_saintClair)

#--- prepare plot for k6 
k6 = sample((1:N_sub_sampl)[Kestim==6],1)
file_res_k6 = paste0("Rdata/save_colsbmfit_subsample/res_iidcolsbm_subnets_90/iidcolsbm_trim_sub_net_",k6,".rds")
res_k6 <- readRDS(file_res_k6)
o_k6 <- get_order_withref(res_k6$best_fit,bestfit_full,o_full)


#---  prepare plot for k7
k7 = sample((1:N_sub_sampl)[Kestim==7],1)
file_res_k7 = paste0("Rdata/save_colsbmfit_subsample/res_iidcolsbm_subnets_90/iidcolsbm_trim_sub_net_",k7,".rds")
res_k7 <- readRDS(file_res_k7)
o_k7 <- get_order_withref(res_k7$best_fit,bestfit_full,o_full)



p_k6 <- plot(res_k6$best_fit, type = "graphon",ord=o_k6)  + theme(plot.title = element_text(size = 12),legend.position = "none") + ggtitle(bquote(m == .(k6) ~ "," ~  hat(R) == 6))
p_k7 <- plot(res_k7$best_fit, type = "graphon",ord=o_k7)  + theme(plot.title = element_text(size = 12),legend.position = "none") + ggtitle(bquote(m == .(k7)~ "," ~  hat(R) == 7))
 




############## recup des Zfull pour tous les neouds avec des NA si le neoud n'est pas présent. 
Z_trim_full <- get_Ztrim(bestfit_full)


########### for all subsample

N_sub_sampl = 30
Ari_estim <- matrix(0,20,N_sub_sampl) 
for (k in 1:N_sub_sampl){
   file_res = paste0("Rdata/save_colsbmfit_subsample/res_iidcolsbm_subnets_90/iidcolsbm_trim_sub_net_",k,".rds")
   res_k <- readRDS(file_res)
   res_k <- res_k$best_fit
   Z_trim <- get_Ztrim(res_k)
  
   Ari_estim[,k] <- sapply(1:20,function(t){
     Z.t.full <- Z_trim_full[,t]
     Z.t <- Z_trim[,t]
     w.t <- which(is.na(Z.t))
     return(adjustedRandIndex(Z.t.full[-w.t],Z.t[-w.t]))
   })
}

ARI_df <- as.data.frame(c(Ari_estim))
names(ARI_df)=c('ARI')

mu  <- mean(ARI_df$ARI, na.rm = TRUE)

p_hist_ARI <- ggplot(ARI_df,aes(x=ARI))
p_hist_ARI <- p_hist_ARI +geom_histogram(aes(y = after_stat(count / sum(count))),binwidth = 0.03,fill = "#4C72B0",
                                       color = "white")
p_hist_ARI <- p_hist_ARI +scale_x_continuous(limits = c(0, 1))+ylab("Frequency")
p_hist_ARI <- p_hist_ARI +  geom_vline(
  xintercept = mu,
  color = "#DD8452",
  linewidth = 1.1,
  linetype = "solid"
) +theme_minimal() + annotate("text", x = mu-0.15, y = Inf, label = "Mean of the ARI",
                              vjust = 2, color = "#DD8452")




final_plot <- (p_k6 + p_k7) / 
  (plot_spacer() + p_hist_ARI + plot_spacer()) +
  plot_layout(heights = c(1, 1), widths = c(1, 2, 1))

final_plot


summary(ARI_df$ARI)
