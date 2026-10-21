#===============================================================================
# This script analyze the income of the fleet per
# group (last part and plot of the manuscript) from our results (estim="nous").
#===============================================================================

library(tidyverse)
library(lubridate)
library(ggplotFL)
library(patchwork)
source("./99_helper_functions.R")
source('./98_functions_resampling.R')

estim="nous"

#---------------- Data
df_xgv_capt <- read_csv("Rdata/sacrois_xgv_capt.csv",show_col_types = FALSE)
df_cvae <- read_csv(file = "./Rdata/xgv_cvae_3_reordered_filtered.csv",show_col_types = FALSE)
trim_peche <- readRDS("./Rdata/trim_peche.rds")
df_traj <- read_csv("Rdata/xgv_100_14-18.csv",show_col_types = FALSE)
df_traj <- df_traj %>% filter(wday(time) %in% c(2,3,4,5,6)) |>
  filter(maree != 12097796)

#----------------- Groups
if(estim=='saintclair'){
  iidcolsbm_trim <- readRDS("./Rdata/save_sbmfit/iidcolsbm_trim.rds")}
if(estim =='nous'){
  iidcolsbm_trim <- readRDS("./Rdata/save_sbmfit/iidcolsbm_trim_round8.rds")
  iidcolsbm_trim <- iidcolsbm_trim$best_fit
}
if(estim =='nousQ5'){
  iidcolsbm_trim <- readRDS("./Rdata/save_sbmfit/iidcolsbm_trim_round7.rds")
  iidcolsbm_trim <- iidcolsbm_trim$model_list[[1]][[5]][[1]] 
}

Q <-iidcolsbm_trim$nb_clusters

best_ord = get_order(iidcolsbm_trim) 
Z_trim <- get_Ztrim_relabel(iidcolsbm_trim)
Z_trim <- relabel_Ztrim(Z_trim,best_ord) # groups pas trim
Z_trim_map <- get_Ztrim_map(Z_trim)  

#-------------------------------------- 


navire_id <- df_cvae |> pull(navire) |> unique()
df_cvae_group <- df_cvae %>% 
  mutate(quarteryear = quarter(time) + 4*(year(time) - 2014)) %>% 
  mutate(navire = factor(navire, levels = navire_id)) %>% 
  group_by(maree) %>% 
  mutate(Group = Z_trim[navire, quarteryear])%>%mutate(quarter=factor(quarter(time)))


df_xgv_capt2 <- df_xgv_capt %>% 
  #  mutate(MONTANT_EUROS_MOYENNE = pmin(MONTANT_EUROS_MOYENNE, 8000)) %>% 
  filter(MONTANT_EUROS_MOYENNE <= 10000) %>% 
  rename(maree = MAREE_ID) %>% 
  left_join(df_cvae_group) %>%
  mutate(quarter = factor(quarter(time)),
         Group = LETTERS[Group]) %>% 
  mutate(date = as_date(paste(year(time), 3*quarter(time),"01", sep = "-")))%>%
  group_by(maree, Group, quarter, date) %>% 
  summarise(MONTANT = sum(MONTANT_EUROS_MOYENNE, na.rm = TRUE))

nb_lines = matrix(0,Q*5*4,4)
r = 0
for (g in 1:Q){
  for (y in 2014:2018){
    for (q in 1:4){
      r = r+1
    w_g_y_q = which(df_xgv_capt2$Group==LETTERS[g]&year(df_xgv_capt2$date)==y & df_xgv_capt2$quarter==q)
    nb_lines[r,]=c(g,y,q, length(w_g_y_q))
    }
  }
}

summary(nb_lines[,4])

nb_lines[which(nb_lines[,4]<=10 & (nb_lines[,4]>0)),]
length(w_E_1)
length(w_E_2)
length(w_E_3)
length(w_E_4)


w_F_1 = which(df_xgv_capt2$Group=='F'&year(df_xgv_capt2$date)==2017 & df_xgv_capt2$quarter==1)
w_F_2 = which(df_xgv_capt2$Group=='F'&year(df_xgv_capt2$date)==2017 & df_xgv_capt2$quarter==2)
w_F_3 = which(df_xgv_capt2$Group=='F'&year(df_xgv_capt2$date)==2017 & df_xgv_capt2$quarter==3)
w_F_4 = which(df_xgv_capt2$Group=='F'&year(df_xgv_capt2$date)==2017 & df_xgv_capt2$quarter==4)



length(w_F_1)
length(w_F_2)
length(w_F_3)
length(w_F_4)

