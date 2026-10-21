#===============================================================================
# This script merge the yearly trajectories
# Filters:
# - Vessels with less than 100 fishing days
# Does not filter (if needed just uncomment the lines):
# - Week-end trajectories
# - A strange maree with duplicated day where the vessels is staying at the 
# harbour
# Extract the merged data to a new .csv
# This is the .csv that must be used in the python script.
# 
# It also export another .csv with just the standardized trajectories
# (centered at 0) 
#=============================================================================== 

library(tidyverse)
library(lubridate)

df_xgv <- map_dfr(.x = c("2014", "2015", "2016", "2017", "2018"),
                  .f = function(y) {
                    df <-  read_csv(file = paste0("Rdata/le_guilvinec_", y, ".csv"))
                    return (df)
                  })


#===============================================================================
# Filters vessels with less than 100 trajectories
# Uncomment to Filter week-end and a strange maree
#===============================================================================
navire100_id <- df_xgv %>% 
  group_by(navire) %>% 
  summarise(nb_maree = n_distinct(maree)) %>% 
  filter(nb_maree >= 100) %>% 
  select(navire) %>% as_vector()

 df_xgv100 <- df_xgv %>% filter(navire %in% navire100_id) #|> 
#   filter(wday(time) %in% c(2,3,4,5,6)) |> 
#   filter(maree != 12097796)

write_csv(df_xgv100, file = "Rdata/xgv_100_14-18.csv")

