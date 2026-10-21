#===============================================================================
# This script do the following:
# - Read and import SACROIS and SIH-GEOLOC files
# - Select Le Guilvinnec trawlers doing daily maree during more than 
# 2 timestamps
# - Match the maree of both tables
# - Smooth the trajectories to match time stamps
# - Export the new datasets
# - Codes can be uncommented for plotting and analysis purposes
#===============================================================================


library(tidyverse)
library(moveVis)
library(sf)
library(lubridate)

#setwd("~/Documents/maccoloc/maccoRoc/")

for (year in c("2014", "2015", "2016", "2017", "2018")) {
  df_traj <- read_delim(paste0("rawdata/Geoloc/SIH-GEOLOC-positions-macco-", year, ".csv"),
                        delim = ";",
                        col_types = cols(
                          ID = col_integer(),
                          NAVIRE = col_integer(),
                          PAVILLON = col_skip(),#col_factor(),
                          DATE_HEURE = col_datetime("%d/%m/%Y %H:%M:%S"),
                          LONGITUDE = col_double(),
                          LATITUDE = col_double(),
                          VITESSE_MOYENNE = col_double(),
                          VITESSE_INSTANTANEE = col_double(),
                          DIRECTION = col_double(),
                          EN_PECHE = col_logical(),
                          ENGF_COD_PRINC = col_factor(),
                          ENGF_LIB_PRINC = col_character(),
                          S_FLOTTILLE_IFREMER_COD = col_factor(),
                          S_FLOTTILLE_IFREMER_LIB = col_character(),
                          TAILLE = col_skip(),#col_double(),
                          PUISSANCE = col_skip(),#col_double(),
                          MAREE_ID = col_integer()
                        ))
  
  df_sacrois <- read_delim(paste0("data/refNavs/SACROIS-macco_", year, ".csv"),
                           delim = "\t", 
                           locale = locale(decimal_mark = ",", 
                                           tz = "UTC",
                                           encoding = "UTF-8"),
                           col_types = cols(
                             ORIGINE = col_factor(),
                             CROISEMENT_DEBARQS_CAPTURES = col_skip(),#col_character(),
                             CROISEMENT_CAPTURES_VENTES = col_skip(),#col_character(),
                             CFR_COD = col_skip(),#col_integer(),
                             NAVS_COD = col_integer(),
                             PAVILLON = col_skip(),#col_factor(),
                             MAREE_ID = col_integer(),
                             MAREE_DATE_DEP = col_datetime("%d-%m-%Y %H:%M:%S"),
                             MAREE_DATE_RET = col_datetime("%d-%m-%Y %H:%M:%S"),
                             QUAL_GEOLOC_DEP_COD = col_skip(),#col_factor(),
                             QUAL_GEOLOC_RET_COD = col_skip(),#col_factor(),
                             MAREE_DATE_DEP_GEOLOC = col_skip(),#col_datetime("%d-%m-%Y %H:%M:%S"),
                             MAREE_DATE_RET_GEOLOC = col_skip(),#col_datetime("%d-%m-%Y %H:%M:%S"),
                             TPS_MER = col_double(),
                             LIEU_COD_DEP_SACAPT = col_skip(),#col_factor(),
                             LIEU_COD_RET_SACAPT = col_skip(),#col_factor(),
                             LIEU_COD_DEP_SACROIS = col_factor(),
                             LIEU_COD_RET_SACROIS = col_factor(),
                             SEQ_ID = col_integer(),
                             DATE_SEQ = col_datetime("%d-%m-%Y %H:%M:%S"),
                             DATE_SEQ_ORIGINE = col_factor(),
                             METIER_COD_SACROIS = col_character(),
                             METIER_QUALITE_SACROIS = col_factor(),
                             METIER_DCF_5_COD = col_character(),
                             METIER_DCF_6_COD = col_character(),
                             ENGIN_COD_SACROIS = col_factor(),
                             ENGIN_COD = col_factor(),
                             ENGIN_ORIGINE = col_factor(),
                             MAILLAGE = col_skip(),#col_factor(),
                             MAILLAGE_ORIGINE = col_skip(),#col_character(),
                             DIMENSION = col_skip(),#col_character(),
                             DIMENSION_ORIGINE = col_skip(),#col_factor(),
                             SECT_COD_SACROIS_NIV1 = col_skip(),#col_character(),
                             SECT_COD_SACROIS_NIV2 = col_skip(),#col_character(),
                             SECT_COD_SACROIS_NIV3 = col_skip(),#col_character(),
                             SECT_COD_SACROIS_NIV4 = col_skip(),#col_character(),
                             SECT_COD_SACROIS_NIV5 = col_skip(),#col_character(),
                             SECT_COD_SACROIS_NIV6 = col_skip(),#col_character(),
                             SECT_COD_SACROIS_NIV7 = col_skip(),#col_character(),
                             SECT_COD_SACROIS = col_character(),
                             TSECT_COD_SACROIS = col_character(),
                             SECT_COD_ORIGINE = col_factor(),
                             GRADIENT_COD = col_factor(),
                             ZEE_COD = col_skip(),#col_factor(),
                             ZPT_COD = col_skip(),#col_factor(),
                             ZEE_ORIGINE = col_skip(),#col_factor(),
                             UGA = col_skip(),#col_character(),
                             TP_NAVIRE_GEOLOC = col_skip(),#col_double(),
                             TP_NAVIRE_SACROIS = col_skip(),#col_double(),
                             TP_NAVIRE_SIPA = col_skip(),#col_double(),
                             TP_NAVIRE_MOYENNE = col_skip(),#col_double(),
                             TP_ENGIN_SIPA = col_skip(),#col_double(),
                             CAPT_ID = col_integer(),
                             ORIGINE_ESP_COD_FAO = col_skip(),#col_factor(),
                             ESP_COD_FAO = col_character(),
                             STOCK_ORGP = col_character(),
                             ORIGINE_QUANT_POIDS_VIF = col_skip(),#col_character(),
                             QUANT_POIDS_VIF_SACROIS = col_skip(),#col_double(),
                             QUANT_POIDS_VIF_SIPA = col_skip(),#col_double(),
                             QUANT_POIDS_VIF_MOYENNE = col_double(),
                             ORIGINE_MONTANT_EUROS = col_skip(),#col_factor(),
                             MONTANT_EUROS_SACROIS = col_skip(),#col_double(),
                             MONTANT_EUROS_SIPA = col_skip(),#col_double(),
                             MONTANT_EUROS_MOYENNE = col_double()
                           ))
  
  #===============================================================================
  #  Obtain the otter trailer from Le Guilvinec doing daily maree
  #===============================================================================
  
  maree_id <- df_sacrois %>% 
    filter(LIEU_COD_DEP_SACROIS %in% c("XGV", "FGV") &
           LIEU_COD_RET_SACROIS %in% c("XGV", "FGV")) %>% 
    filter(str_starts(METIER_COD_SACROIS, pattern = "OT")) %>% 
    filter(TPS_MER <= 24) %>% 
    select(MAREE_ID) %>% 
    as_vector() %>% 
    unique()
  
  
  traj_guilvinnec <- df_traj %>% 
    filter(MAREE_ID %in% maree_id)
  
  
  traj_guilvinnec %>% 
    group_by(NAVIRE) %>% 
    summarise(n_distinct(MAREE_ID)) %>% 
    print(n = Inf)
  
  
  #===============================================================================
  #  Removing Maree with just one time stamp
  #===============================================================================
  
  traj_guilvinnec <- traj_guilvinnec %>% 
    group_by(MAREE_ID) %>% 
    filter(n() >= 2)
  
  
  #===============================================================================
  # Filtering for maree starting and ending on the same day
  #===============================================================================
  #sf_guilvinnec_2014 <- st_set_geometry()
  same_day <- traj_guilvinnec %>% 
    group_by(MAREE_ID) %>% 
    summarise(DEPARTURE = min(DATE_HEURE), 
              RETURN = max(DATE_HEURE)) %>%
    mutate(SAME_DAY = day(DEPARTURE) == day(RETURN)) %>% 
    filter(SAME_DAY == TRUE) %>% 
    select(MAREE_ID) %>% 
    as_vector()
  
  traj_guilvinnec <- traj_guilvinnec %>% filter(MAREE_ID %in% same_day)
  navire_maree <- traj_guilvinnec %>% select(NAVIRE, MAREE_ID) %>% distinct()
  
  #===============================================================================
  # Plotting length of maree
  #===============================================================================
  
  # traj_guilvinnec %>%
  #   group_by(NAVIRE, MAREE_ID) %>%
  #   summarise(nb_step = n()) %>%
  #   ggplot(aes(x = nb_step)) +
  #   geom_histogram()
  # 
  # 
  # traj_guilvinnec %>%
  #   group_by(NAVIRE, MAREE_ID) %>%
  #   summarise(nb_step = n()) %>%
  #   group_by(nb_step) %>%
  #   summarise( n()) %>%
  #   print(n = Inf)
  # 
  # traj_guilvinnec %>%
  #   group_by(NAVIRE) %>%
  #   summarise(nb_marees = n_distinct(MAREE_ID)) %>%
  #   ggplot(aes(x = nb_marees)) +
  #   geom_histogram()
  
  
  #===============================================================================
  # Smoothing trajectories. 1 point per hour at h:00 minute
  #===============================================================================
  
  move_data  <- df2move(df = traj_guilvinnec,
                        proj =  "+init=epsg:4326 +proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0",#sf::st_crs(4326),
                        x = "LONGITUDE",
                        y = "LATITUDE",
                        time = "DATE_HEURE",
                        track_id = "MAREE_ID")
  
  move_data <- align_move(move_data, res = 60, digit = 0, unit = "mins")
  
  # str(move_data)  
  
  
  #===============================================================================
  # Exporting the modified dataset
  #===============================================================================
  df_guilv <- move_data@data %>% 
    mutate(month = month(time))
  df_guilv <- as_tibble(df_guilv, rownames = "maree") %>% 
    mutate(maree = as.integer(str_sub(maree, 2, 9)))
  
  time2pos <- function(t) {
    1 + hour(t)#+ minute(t)/15 + 4*hour(t)
  }
  
  #time2pos(df_guilv$time[2])
  
  df_guilv_full <- map_dfr(
    .x = seq_along(same_day),
    .f = function(id) {
      data <- df_guilv %>% filter(maree == same_day[id])
      dates <- data$time[1]
      month <- data$month[1]
      tb <- tibble(maree = rep(data$maree[1], 24*1),
                   x = rep(-4.28, 24*1),
                   y = rep(47.78, 24*1),
                   month = rep(data$month[1], 24*1),
                   time = date(data$time[1])
      )
      tb[time2pos(data$time), c("x","y")] <- data[, c("x", "y")]
      return(tb)
    })
  
  df_guilv_full <- df_guilv_full %>% 
    left_join(navire_maree %>% rename(maree = MAREE_ID, navire = NAVIRE), by = "maree")
  
  write_csv(df_guilv_full, file = paste0("Rdata/le_guilvinec_", year, ".csv"))
}


