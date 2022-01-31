# R 3.6.2

# Compare aquametbio metrics vs DEQ generated metrics (from pinned data on server)

# libraries
library(tidyverse)
library(sf)
library(pool)
library(DBI)
library(pins)
library(config)
library(lubridate)
library(dbplyr)
library(readxl)


# get configuration settings
conn <- config::get("connectionSettings")

# use API key to register board
board_register_rsconnect(key = conn$CONNECT_API_KEY,  #Sys.getenv("CONNECT_API_KEY"),
                         server = conn$CONNECT_SERVER)#Sys.getenv("CONNECT_SERVER"))

# pull pinned VSCI metrics
VSCI <- pin_get('ejones/VSCIresults', board = 'rsconnect')

# bring in aquametbio metrics
aquamet <- left_join(
  read_excel('data/aquametBioDataOut/aquametBio_run200_01192022.xlsx', sheet = 'Taxa Metrics'),
  read_excel('data/aquametBioDataOut/aquametBio_run200_01192022.xlsx', sheet = 'Tolerance Metrics'), by = 'BENSAMPID') %>% 
  left_join(read_excel('data/aquametBioDataOut/aquametBio_run200_01192022.xlsx', sheet = 'FFG Metrics'), by = 'BENSAMPID') %>% 
  left_join(read_excel('data/aquametBioDataOut/aquametBio_run200_01192022.xlsx', sheet = 'Habit Metrics'), by = 'BENSAMPID') %>% 
  left_join(read_excel('data/aquametBioDataOut/aquametBio_run200_01192022.xlsx', sheet = 'Dominance Metrics'), by = 'BENSAMPID')
  
# compare datasets
compare <- left_join(aquamet, VSCI, by = c('BENSAMPID' = 'BenSampID')) %>% 
  group_by(BENSAMPID) %>% 
  rowwise() %>% 
  mutate(
    #`Comp_Family Total Taxa` =   
    `Comp_Family EPT Taxa`  = abs(sum(`Family EPT Taxa`, -  EPT_NTAX, na.rm = T)),
    `Comp_%Ephem` = abs(sum(`%Ephem`, - EPHEPIND, na.rm= T)),        
    #`Comp_%PT - Hydropsychidae` not in aquametbio
#    `Comp_%FamilyScraper` = abs(sum(`%FamilyScraper`, - SCRPPTAX, na.rm=T)), # habit will change from family to genus level, not good comparison
    `Comp_%Chiro` = abs(sum(`%Chiro`, - CHIRPTAX, na.rm = T))#,
   # `Comp_Family %2 Dominant`   # not in aquamet
   # `Comp_Family HBI`  = abs(sum(`Family HBI`,- WTD_TV, na.rm= T)) not good comparison
   ) %>% 
  dplyr::select(StationID, BENSAMPID, `Comp_Family EPT Taxa`, `Family EPT Taxa`, EPT_NTAX,
                `Comp_%Ephem`, EPHEPIND, `%Ephem`,
        #        `Comp_%FamilyScraper`, `%FamilyScraper`, SCRPPTAX,
                `Comp_%Chiro`, `%Chiro`,  CHIRPTAX)#,
               # `Comp_Family HBI`, `Family HBI`, WTD_TV          )


write.csv(compare, 'data/aquametBioDataOut/compareMetrics.csv')
