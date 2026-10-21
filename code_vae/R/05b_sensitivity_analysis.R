#===============================================================================
# This script does the following:
# - Import generated trajectories from cvae3
#   The imported files contains the following:
#   - computed sensitivity indices 
#   - mu_1, mu_2 and eta 
#   - 1200 trajectoires (10 mu1, 10 mu_2, et 12 mois)
#   - Compute dmax and theta[dmax]
#     - With the article euclidean distance (uncomment to use)
#     - With a corrected projected version (comment if needed)
#  - Compute the sobol MDA indices 
#
#===============================================================================

library(tidyverse)
library(ranger)
library(sobolMDA) # available here: https://gitlab.com/drti/sobolmda
# To download and install from the .tar.gz using
# install.packages("sobolmda-sobolMDA.tar.gz")
#devtools::install_git("https://gitlab.com/drti/sobolmda")
rm(list = ls()) # Nettoyage
source("./99_helper_functions.R")
# Chargement données

donnees <- read.table("./pydata/generated_traj_cvae3.csv",
                      sep = ",",
                      header = TRUE,
                      row.names = 1) %>% 
  filter(mu1 <= 1.2) 


donnees

#===============================================================================
# Article version
#===============================================================================
# donnees <- donnees |> 
#   mutate(dist = sqrt((x*10)**2 + (y*10)**2)) |> 
#   mutate(theta = atan(y/x)) |> 
#   group_by(mu1, mu2, cos_eta, sin_eta) |> 
#   mutate(max_dist = max(dist),
#          max_theta = theta) |> 
#   filter(dist == max_dist)


#===============================================================================
# New and more exact version
#===============================================================================

donnees <- donnees |> 
  mutate(lon = std2coord(x, -4.28),
         lat = std2coord(y, 47.78)) |> 
  mutate(x = (lon+4.28)*60*1.8*cos(lat/180*pi),
         y = (lat-47.78)*60*1.8) |> 
  mutate(dist = sqrt(x**2 + y**2)) |> 
  mutate(theta = atan(y/x)) |> 
  group_by(mu1, mu2, cos_eta, sin_eta) |> 
  mutate(max_dist = max(dist),
         max_theta = theta) |> 
  filter(dist == max_dist)

#===============================================================================
# New and more exact version using just lon and lat
#===============================================================================

# donnees <- donnees |> 
#   mutate(lon = std2coord(x, -4.28),
#          lat = std2coord(y, 47.78)) |> 
#   mutate(dist = acos(sin(47.78)*sin(lat)+cos(47.78)*cos(lat)*cos(lon-(-4.28)))*63.71) |> 
#   mutate(theta = atan2( 
#     sin(lon + 4.28) * cos(lat) , 
#     cos(47.78) * sin(lat) - sin(47.78)*cos(lat)*cos(lon+4.28))) |> 
#   group_by(mu1, mu2, cos_eta, sin_eta) |> 
#   mutate(max_dist = max(dist),
#          max_theta = theta) |> 
#   filter(dist == max_dist)


best_model_max_theta <- sobolMDA::ranger(formula = max_theta ~ mu1 + mu2 + cos_eta + sin_eta,
                                         num.trees = 500, data = donnees,
                                         importance = "sobolMDA", mtry = 3, 
                                         sample.fraction = 0.1)
best_model_max_dist <- sobolMDA::ranger(formula = max_dist ~ mu1 + mu2 + cos_eta + sin_eta,
                                        num.trees = 500, data = donnees,
                                        importance = "sobolMDA", mtry = 3, 
                                        sample.fraction = 0.1)

best_model_max_dist
best_model_max_theta
# Graphiques --------------------------------------------------------------

importances <- sobolMDA::importance(best_model_max_dist) %>% 
  enframe(name = "rowname", value = "Importance") %>% 
  mutate(metric = "D[max]") %>% 
  bind_rows(sobolMDA::importance(best_model_max_theta) %>% 
              enframe(name = "rowname", value = "Importance") %>%
              mutate(metric = "theta[max]")) %>% 
  mutate(rowname = factor(rowname,
                          levels = c("mu1", "mu2", "cos_eta", "sin_eta") %>% 
                            rev())) 


importances

my_theme <- function(base_size = 12){
  theme_bw()  %+replace%
    theme(
      panel.border = element_rect(colour = "black", 
                                  fill = rgb(0, 0, 0, 0)),
      # plot.background = element_rect(fill = "white"),# bg around panel
      legend.background = element_blank(), 
      text = element_text(family = "LM Roman 10", size = base_size),
      axis.title = element_text(size = rel(1)),
      legend.text = element_text(size = rel(1)),
      legend.title = element_text(size = rel(1)),
      axis.text = element_text(size = rel(1)),
      strip.background = element_rect(fill = "lightgoldenrod1",
                                      color = "black"),
      plot.subtitle = element_text(hjust = 0.5, size = rel(1)),
      plot.title = element_text(face = "bold", size = rel(1.2), hjust = 0.5))
}
theme_set(my_theme(20))
ggplot(importances) + 
  aes(y = Importance, x = rowname, col = metric) +
  geom_point(size = 5) +
  ylim(c(-0.01,1)) +
  scale_x_discrete(labels = c(
    "mu1" = expression(z[1]),
    "mu2" = expression(z[2]),
    "sin_eta" = expression(sin(m)),
    "cos_eta" = expression(cos(m))
  )) +
  scale_color_discrete("Trajectory \n Charateristic",
                       labels = c(latex2exp::TeX("$D_{max}$"),
                                  latex2exp::TeX("$\\theta_{max}$"))) +
  labs(x = "Variable",
       y = "Importance") #+ 
#  facet_wrap(~metric, labeller = label_parsed)
