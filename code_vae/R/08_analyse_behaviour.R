rm(list=ls())
#===============================================================================
# Individual and collective behaviours
# (Figures 9, 10 and A2 of the paper)
#
# This script does the following
# - Import fitted colSBM model
# - Import latent space of cvae3
# - Plot the 20 quarterly proximity graphs (Figure A2)
# - Compute DSI and GSI for individual
# - Plot DSI and GSI for individual (Figure 9)
# - Plot mesoscale structure and groups (Figure 10)
#
# Inputs:
# - Rdata/xgv_cvae_3_reordered_filtered.csv (03_computing_bc.R)
# - Rdata/trim_peche.rds (07b_...R): for each quarter,
#   the vessels (1..33) that fished during the quarter
# - Rdata/bc_lag_cvae.rds (06_build_proximity_graph.R): BC between the trip
#   of vessel b on day t and the trip of vessel b' on day t+1, for all couples
#   (b, b'), array (days x 33 x 33)
# - Rdata/save_sbmfit/iidcolsbm_trim_round8.rds (07b_...R): colSBM fit on the
#   20 quarterly proximity graphs
#
# Outputs (in images2/):
# - proxy_graph_trim_nous.png      : Figure A2 (proxy_graph_trim.png in the paper)
# - bc_group_colsbm_trim_cvae.png  : Figure 9
# - iidcolsbm_group_trim.png       : Figure 10
#
# The groups are the blocks of the colSBM, relabelled A, B, ... with the
# ordering functions of 98_functions_resampling.R.
#===============================================================================


library(tidyverse)
library(lubridate)
library(matrixStats)
library(patchwork)
source("./99_helper_functions.R")
source('./98_functions_resampling.R')
#-------------------------------- 
whereImages <-"./images2" 

################################################
########   data
################################################


#---------- cvae 
df_cvae <- read.csv(file = "./Rdata/xgv_cvae_3_reordered_filtered.csv")
df_cvae$time <- as.Date(df_cvae$time)

trim_peche <- readRDS("./Rdata/trim_peche.rds")
bc_lag_cvae_ar <- readRDS("./Rdata/bc_lag_cvae.rds")

#------------ colsbm
# Selected model of the colSBM fit of 07b_...R (6 groups). Earlier fits are
# in archive/ (see archive/README.md).
iidcolsbm_trim <- readRDS("./Rdata/save_sbmfit/iidcolsbm_trim_round8.rds")$best_fit

Q <-iidcolsbm_trim$nb_clusters # number of groups


#===============================================================================
# Export the matrices from the colSBM
#===============================================================================
# 33 x 33 adjacency matrix of each quarter, with NA for the vessels that did
# not fish during the quarter and for the couples of vessels that never
# fished the same day (coded -1 in the colSBM object)
trim_matrices <- lapply(
  seq(20),
  function(l) {
    A <- matrix(NA, 33, 33)
    A[trim_peche[[l]], trim_peche[[l]]] <- matrix(iidcolsbm_trim$A[[l]])
    A[A==-1] <- NA
    diag(A) <- NA
    return(A)
  }
)


#===============================================================================
# Compute the map of memberships for each vessel
#===============================================================================
# Z_trim: 33 x 20 matrix, group of each vessel in each quarter (NA if absent),
# with groups relabelled so that 1, 2, ... correspond to the plotting order
# best_ord (A, B, ... in the figures).
# Note: get_Ztrim_relabel also swaps the labels of the two largest groups
# (and of their most connected groups) in some quarters, so that vessel 1 is
# never in the one of the two largest groups it belongs to less often.
# Z_trim_map: most frequent group of each vessel over the quarters.
best_ord = get_order(iidcolsbm_trim) 
Z_trim <- get_Ztrim_relabel(iidcolsbm_trim)
Z_trim <- relabel_Ztrim(Z_trim,best_ord)
Z_trim_map <- get_Ztrim_map(Z_trim)  


#===============================================================================
# Plot list of adjacency matrices reordered by Z_map with diagonal in color
#===============================================================================
# Exploratory, not in the paper. The diagonal is coloured with the raw labels
# of the colSBM (not relabelled); note that hue_pal()(5) gives 5 colours while
# the fit has Q = 6 groups (vessels of group 6 get no colour).
p_list_trim <- lapply(
  seq(20),
  function(l) {
    Zl <- rep(NA,33)
    Zl[trim_peche[[l]]] <- iidcolsbm_trim$Z[[l]]
    trim_matrices[[l]][order(Z_trim_map), order(Z_trim_map)] |> 
      t() |> 
      reshape2::melt() |> 
      ggplot(aes(x = Var1, y = Var2, fill = value)) +
      geom_tile(show.legend = FALSE) +
      scale_y_reverse() +
      coord_fixed(expand = FALSE) +
      scale_fill_gradient2(low = "white", high = "black", limits = c(0,1), na.value = "grey80") +
      annotate(geom = "tile", x = seq(33), y = seq(33), fill = scales::hue_pal()(5)[Zl[order(Z_trim_map)]]) +
      geom_hline(yintercept = seq(0,33) + 0.5, linewidth = .1,color = "grey20") +
      geom_vline(xintercept = seq(0,33)+ 0.5, linewidth = .1, color = "grey20") +
      theme_void() 
  })

patchwork::wrap_plots(p_list_trim, ncol = 4, tag_level = 'keep', byrow = TRUE) + 
  plot_annotation(tag_levels = '1')


#===============================================================================
# Plot list of adjacency matrices unordered (see plot from appendix)
#===============================================================================
# Figure A2: the 20 quarterly proximity graphs (rows: quarters, columns:
# years). Black: edge; pink: NA (vessel not active during the quarter, or
# couple of vessels that never fished the same day).
array(unlist(trim_matrices), dim = c(33,33,20)) |>
  reshape2::melt() |> 
  mutate(Quarter = (Var3-1)%%4 + 1) |> 
  mutate(Year = (Var3-1)%/%4 + 1) |> 
  mutate(Year = c("2014","2015","2016","2017","2018")[Year]) |> 
  ggplot(aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(show.legend = FALSE) +
  scale_y_reverse() +
  coord_fixed(expand = FALSE) +
  scale_fill_gradient2(low = "white", high = "black", limits = c(0,1), na.value = "lightpink1") +
#       annotate(geom = "tile", x = seq(33), y = seq(33), fill = scales::hue_pal()(5)[Zl]) +
  geom_hline(yintercept = seq(0,33) + 0.5, linewidth = .1,color = "grey20") +
  geom_vline(xintercept = seq(0,33)+ 0.5, linewidth = .1, color = "grey20") +
  facet_grid(Quarter ~ Year) +
  xlab("") + ylab("")+
  theme_bw(base_size = 15) +
  theme(axis.text = element_blank(), axis.ticks = element_blank())

ggsave("./images2/proxy_graph_trim_nous.png", scale = 2)

#trim_matrices[[1]]

#===============================================================================
# Compute the BC coeefficient for DSI and GSI
#===============================================================================
# DSI (daily stability index) of a vessel: proportion of its couples of trips
# on adjacent days with BC >= threshold.
# GSI (global stability index): same for random couples of its trips.
#
# Latent parameters in an array (param x day x vessel): one row per vessel and
# per week day of the period (full_join), NA when the vessel did not fish, so
# that adjacent indices are adjacent week days. Rows 1-3: mu_k, rows 4-6: sig_k.
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
# Compute DSI
#===============================================================================
# BC between the trips of day t and day t+1 of each vessel (product over the
# 3 dimensions): (days - 1) x 33 matrix, NA if one of the two days is missing
bclag_navire_cvae <- bhattacharyya_1d(
  mu1 = ar_cvae_mu[,1:(dim(ar_cvae_mu)[2]-1),], 
  mu2 = ar_cvae_mu[,2:dim(ar_cvae_mu)[2],], 
  var1 = ar_cvae_sig2[,1:(dim(ar_cvae_mu)[2]-1),], 
  var2 = ar_cvae_sig2[,2:dim(ar_cvae_mu)[2],]) %>% 
  apply(c(2,3), "prod")


#===============================================================================
# Compute GSI
#===============================================================================
# For each vessel, BC between 2e4 random couples of distinct trips
# (2e4 x 33 matrix). No seed: the GSI vary slightly from one run to the next.
mu_cvae <- as.matrix(df_cvae |> select(starts_with("mu")))
sig2_cvae <- as.matrix(df_cvae |> select(starts_with("sig")))


bc_navire <- sapply(seq(33),
                       function(b) {
                         ind <- df_cvae$navire == unique(df_cvae$navire)[b]
                         ord1 <- sample(nrow(mu_cvae[ind,]), size = 3e4, replace = TRUE)
                         ord2 <- sample(nrow(mu_cvae[ind,]), size = 3e4, replace = TRUE)
                         diff <- ord1 != ord2
                         ord1 <- ord1[diff][1:2e4]
                         ord2 <- ord2[diff][1:2e4]
                         bc_cvae <- rowProds(bhattacharyya_1d(mu_cvae[ind,][ord1,], 
                                                     mu_cvae[ind,][ord2,], 
                                                     sig2_cvae[ind,][ord1,], 
                                                     sig2_cvae[ind,][ord2,]))
                         return(bc_cvae)
                       }
)

#===============================================================================
# Treshold DSI and GSI
#===============================================================================
threshold = 0.8

bclag_navire_th <- colMeans(bclag_navire_cvae >= threshold, na.rm = TRUE)
bc_navire_th <- colMeans(bc_navire >= threshold, na.rm = TRUE)

# lag = DSI, dist = GSI
bc_seuil_navire <- tibble(navire = seq(33),
                          lag = bclag_navire_th,
                          dist = bc_navire_th,
                          Model = 'c-vae')


#===============================================================================
# Computing fleet and group average
#===============================================================================
# Z_trim_mean: 33 x Q matrix, proportion of its quarters each vessel spends
# in each group
Z_trim_mean <- outer(Z_trim, seq(Q), "==")  |> 
  apply(c(1,3), "mean", na.rm = TRUE) 


Mat_navire =cbind(seq(33),Z_trim_mean)
colnames(Mat_navire) <- c('navire',LETTERS[1:Q])
Mat_navire <- as_tibble(Mat_navire)


# Group averages: mean of the indices of the vessels, weighted by the
# proportion of time they spend in the group (squares in Figure 9)
mean_group <- bc_seuil_navire%>%
  left_join(Mat_navire) %>% 
  filter(Model == "c-vae") %>% 
  pivot_longer(cols = LETTERS[1:Q],
               names_to = "Group", values_to = "Prop") |> 
  group_by(Group) |> 
  summarise(lag = sum(as.numeric(Prop)*as.numeric(lag))/sum(Prop),
            dist = sum(as.numeric(Prop)*as.numeric(dist))/sum(Prop))

# Average of the vessels (triangle in Figure 9)
mean_bc_boat <- bc_seuil_navire %>%
  filter(Model == "c-vae") %>% 
  select(lag, dist) %>% 
  colMeans()

# Fleet (diamond in Figure 9): same indices computed on all the trips,
# whatever the vessel. GSI on 1e5 random couples of trips; DSI on the BC
# between trips of adjacent days of any two vessels (bc_lag_cvae.rds).
n_traj <- nrow(df_cvae)
ord1 <- sample(n_traj, 2e5, replace = TRUE)
ord2 <- sample(n_traj, 2e5, replace = TRUE)
ind <- ord1 != ord2
ord1 <- ord1[ind][1:1e5]
ord2 <- ord2[ind][1:1e5]

bc_cvae <- rowProds(bhattacharyya_1d(mu_cvae[ord1,],
                                     mu_cvae[ord2,],
                                     sig2_cvae[ord1,],
                                     sig2_cvae[ord2,]))

anno <- tibble(dist = mean(bc_cvae > .8),
               lag = mean(bc_lag_cvae_ar >.8, na.rm = TRUE),
               Model = c("c-vae"))


#===============================================================================
# Plot DSI and GSI for individual (Figure 9)
#===============================================================================
# Each vessel is a pie at (DSI, GSI) showing the proportions of its group
# memberships. Vessels 2, 6 and 19 (Figures 5-6) are labelled in red.
bc_seuil_navire%>%
  left_join(Mat_navire) %>%
  filter(Model == "c-vae") %>%
  ggplot() + xlab("DSI") + ylab("GSI") + xlim(0.15,NA) + ylim(0,0.2) +
  annotate("point", x = mean_bc_boat[1], y = mean_bc_boat[2], 
           col = "black", size = 8, shape = 17) +  
  ggrepel::geom_text_repel(data = bc_seuil_navire%>% 
                             left_join(Mat_navire) %>% 
                             filter(Model == "c-vae") %>% 
                             filter(! navire %in% c(2, 6, 19)),
                           mapping = aes(x = lag, y = dist, 
                                         col = navire, label = navire), 
                           show.legend = FALSE, size = 5, #min.segment.length = 0.5, 
                           max.overlaps = 20,
                           force = 5, force_pull =5,
                           # direction = "y", #nudge_y = -0.03, nudge_x =.005,
                           box.padding = 1#,
                           #col = "black"
  ) +
  scatterpie::geom_scatterpie(aes(y = dist, x = lag), pie_scale = 1, alpha = .9,
                              cols = LETTERS[1:Q], color = "grey20",
                              linewidth = .5,
                              data = (bc_seuil_navire%>%  
                                        left_join(Mat_navire)) %>% 
                                filter(Model == "c-vae"), 
                              legend_name = "Group") + coord_equal() +
  ggrepel::geom_label_repel(data = bc_seuil_navire%>% 
                              left_join(Mat_navire) %>% 
                              filter(Model == "c-vae") %>% 
                              filter(navire %in% c(2, 6, 19)),
                            mapping = aes(x = lag, y = dist, label = navire), 
                            show.legend = FALSE , size = 7, 
                            min.segment.length = 0.1, 
                            direction = "y", 
                            nudge_y = -0.03, 
                            nudge_x =.005,
                            point.padding = 4,
                            col = "red") +
  geom_point(data = mean_group, 
             mapping = aes(x = lag, y = dist, fill = Group), 
             size = 8, shape = 22, stroke = 2, show.legend = FALSE) + 
  theme_bw(base_size = 20) +
  geom_point(data = anno %>% filter(Model == "c-vae"), 
             inherit.aes = FALSE, aes(x = lag, y = dist), 
             col = "black", size = 8, shape = 18)


ggsave("./images2/bc_group_colsbm_trim_cvae.png", dpi = 300)


#===============================================================================
# Plot the matrix view of the mesoscale structure and the block memberhips.
#===============================================================================
# Figure 10. Left: proportion of the quarters each vessel spends in each
# group, vessels 2, 6 and 19 framed in red. Right: connectivity matrix of the
# colSBM, blocks in the order best_ord with widths = group proportions.

#------------------- matrix
# Axis breaks at the block boundaries, group letters on the secondary axes
p <-   plot(iidcolsbm_trim, ord = best_ord)
breaks1 <-  c(p$data$xmax[1:Q], 1)
breaks2 <- (breaks1[-1] + breaks1[-(Q+1)])/2
labels1 <- round(breaks1, 2)
labels2 <- LETTERS[1:Q]
p <- p +  scale_y_continuous(labels = labels1, breaks = breaks1,trans = "reverse", 
                      sec.axis = sec_axis(trans = "identity", 
                                          labels = labels2, 
                                          breaks = breaks2)) +
   scale_x_continuous(labels = labels1, breaks = breaks1, 
                      sec.axis = sec_axis(trans = "identity", 
                                          labels = labels2, 
                                          breaks = breaks2)) + 
   theme_bw(base_size = 20) + 
   guides(fill = guide_colorsteps()) + 
   theme(axis.ticks = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1)) 


#--------------- groups 
Mat_groups <- (outer(Z_trim, seq(Q), "==") |> apply(c(1,3), "mean", na.rm = TRUE)) |> t() |> 
    reshape2::melt() %>% 
    rename(Group = Var1, Vessel = Var2, Prop = value) %>% 
    mutate(Group = LETTERS[Group]) 

gg_Mat_groups <- ggplot(Mat_groups ,aes(x = Group, y = Vessel, alpha = Prop, fill=Group)) +
    geom_tile() + #scale_fill_gradient2() +
    scale_fill_discrete() +
    geom_hline(yintercept = seq(0,33) + .5) +
    geom_vline(xintercept = seq(0,5) + .5) + # assumes 5 groups: the line after group F is missing when Q = 6
    scale_alpha_continuous(range = c(0,1)) +
    geom_hline(yintercept = c(5, 6, 1, 2, 18, 19) + 0.5, 
               col = "red", linewidth = 1.1) +
    annotate("segment", x = rep( c(0.5, Q+.5), 3), xend = rep(c(0.5, Q+.5), 3), 
             y = rep(c(5.5, 1.5, 18.5), 2), 
             yend = rep(c(6.5, 2.5, 19.5), 2), col = "red", linewidth = 1.1) +
    scale_y_continuous(labels = c("2", "6", "19"), breaks = c(2, 6, 19)) +
    coord_fixed(ratio = 2/4, expand = FALSE, clip = 'off') +
    guides(fill = guide_none(),
           alpha = guide_bins(reverse = TRUE)) +
    theme_minimal(base_size = 20, base_line_size = 0)

gg_Mat_groups + guide_area() +  p +    plot_layout(guides = 'collect', design = c("CBA"), widths = c(6.5,1,2)) 

ggsave("./images2/iidcolsbm_group_trim.png", dpi = 300)

 
 