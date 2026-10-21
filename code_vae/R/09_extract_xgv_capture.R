#===============================================================================
# This script exports the capture data concerning the xgv fleet
#===============================================================================

library(tidyverse)
library(lubridate)

df_cvae <- read_csv(file = "./Rdata/xgv_cvae_3_reordered_filtered.csv")

# df_sacrois_xgv_capt <- map_dfr(.x = c("2014", "2015", "2016", "2017", "2018"),
#                                .f = function(y) {
#                                  df <- read_delim(paste0("data/refNavs/SACROIS-macco_", y, ".csv"),
#                                                   delim = "\t", 
#                                                   locale = locale(decimal_mark = ",", 
#                                                                   tz = "UTC",
#                                                                   encoding = "UTF-8"),
#                                                   col_types = cols(
#                                                     ORIGINE = col_factor(),
#                                                     CROISEMENT_DEBARQS_CAPTURES = col_skip(),#col_character(),
#                                                     CROISEMENT_CAPTURES_VENTES = col_skip(),#col_character(),
#                                                     CFR_COD = col_skip(),#col_integer(),
#                                                     NAVS_COD = col_integer(),
#                                                     PAVILLON = col_skip(),#col_factor(),
#                                                     MAREE_ID = col_integer(),
#                                                     MAREE_DATE_DEP = col_datetime("%d-%m-%Y %H:%M:%S"),
#                                                     MAREE_DATE_RET = col_datetime("%d-%m-%Y %H:%M:%S"),
#                                                     QUAL_GEOLOC_DEP_COD = col_skip(),#col_factor(),
#                                                     QUAL_GEOLOC_RET_COD = col_skip(),#col_factor(),
#                                                     MAREE_DATE_DEP_GEOLOC = col_skip(),#col_datetime("%d-%m-%Y %H:%M:%S"),
#                                                     MAREE_DATE_RET_GEOLOC = col_skip(),#col_datetime("%d-%m-%Y %H:%M:%S"),
#                                                     TPS_MER = col_double(),
#                                                     LIEU_COD_DEP_SACAPT = col_skip(),#col_factor(),
#                                                     LIEU_COD_RET_SACAPT = col_skip(),#col_factor(),
#                                                     LIEU_COD_DEP_SACROIS = col_factor(),
#                                                     LIEU_COD_RET_SACROIS = col_factor(),
#                                                     SEQ_ID = col_integer(),
#                                                     DATE_SEQ = col_datetime("%d-%m-%Y %H:%M:%S"),
#                                                     DATE_SEQ_ORIGINE = col_factor(),
#                                                     METIER_COD_SACROIS = col_character(),
#                                                     METIER_QUALITE_SACROIS = col_factor(),
#                                                     METIER_DCF_5_COD = col_character(),
#                                                     METIER_DCF_6_COD = col_character(),
#                                                     ENGIN_COD_SACROIS = col_factor(),
#                                                     ENGIN_COD = col_factor(),
#                                                     ENGIN_ORIGINE = col_factor(),
#                                                     MAILLAGE = col_skip(),#col_factor(),
#                                                     MAILLAGE_ORIGINE = col_skip(),#col_character(),
#                                                     DIMENSION = col_skip(),#col_character(),
#                                                     DIMENSION_ORIGINE = col_skip(),#col_factor(),
#                                                     SECT_COD_SACROIS_NIV1 = col_skip(),#col_character(),
#                                                     SECT_COD_SACROIS_NIV2 = col_skip(),#col_character(),
#                                                     SECT_COD_SACROIS_NIV3 = col_skip(),#col_character(),
#                                                     SECT_COD_SACROIS_NIV4 = col_skip(),#col_character(),
#                                                     SECT_COD_SACROIS_NIV5 = col_skip(),#col_character(),
#                                                     SECT_COD_SACROIS_NIV6 = col_skip(),#col_character(),
#                                                     SECT_COD_SACROIS_NIV7 = col_skip(),#col_character(),
#                                                     SECT_COD_SACROIS = col_character(),
#                                                     TSECT_COD_SACROIS = col_character(),
#                                                     SECT_COD_ORIGINE = col_factor(),
#                                                     GRADIENT_COD = col_factor(),
#                                                     ZEE_COD = col_skip(),#col_factor(),
#                                                     ZPT_COD = col_skip(),#col_factor(),
#                                                     ZEE_ORIGINE = col_skip(),#col_factor(),
#                                                     UGA = col_skip(),#col_character(),
#                                                     TP_NAVIRE_GEOLOC = col_skip(),#col_double(),
#                                                     TP_NAVIRE_SACROIS = col_skip(),#col_double(),
#                                                     TP_NAVIRE_SIPA = col_skip(),#col_double(),
#                                                     TP_NAVIRE_MOYENNE = col_skip(),#col_double(),
#                                                     TP_ENGIN_SIPA = col_skip(),#col_double(),
#                                                     CAPT_ID = col_integer(),
#                                                     ORIGINE_ESP_COD_FAO = col_skip(),#col_factor(),
#                                                     ESP_COD_FAO = col_character(),
#                                                     STOCK_ORGP = col_character(),
#                                                     ORIGINE_QUANT_POIDS_VIF = col_skip(),#col_character(),
#                                                     QUANT_POIDS_VIF_SACROIS = col_skip(),#col_double(),
#                                                     QUANT_POIDS_VIF_SIPA = col_skip(),#col_double(),
#                                                     QUANT_POIDS_VIF_MOYENNE = col_double(),
#                                                     ORIGINE_MONTANT_EUROS = col_skip(),#col_factor(),
#                                                     MONTANT_EUROS_SACROIS = col_skip(),#col_double(),
#                                                     MONTANT_EUROS_SIPA = col_skip(),#col_double(),
#                                                     MONTANT_EUROS_MOYENNE = col_double()
#                                                   ))
#                                  df <-   df %>%  filter(MAREE_ID %in% df_cvae$maree) %>% 
#                                    mutate(ENGIN = str_sub(METIER_COD_SACROIS, 1, 3)) %>%  
#                                    select(MAREE_ID, ENGIN, METIER_DCF_5_COD, ESP_COD_FAO, QUANT_POIDS_VIF_MOYENNE, MONTANT_EUROS_MOYENNE) #%>% 
#                                  # group_by(MAREE_ID) %>%
#                                  # summarise(montant = sum(MONTANT_EUROS_MOYENNE),
#                                  #           engin = str_sub(METIER_COD_SACROIS, 1, 3),
#                                  #           metier = METIER_DCF_5_COD) %>% 
#                                  # distinct()
#                                  df
#                                })
# 
# write_csv(df_sacrois_xgv_capt, "Rdata/sacrois_xgv_capt.csv")
