# This script builds out the removeRareTaxa() function
# it is intended to work with any taxa list (bugs, fish) and be flexible enough to return
# allow users to different thresholds to remove taxa

# Jason spreadsheet with math
#benthics <- read_excel('data/Jason/benthics_Emma_June11.xlsx', sheet = 'Ref_Bugs_All')
# 
# refstress <- read_csv('data/dataOut/stationsSpatial_withAllData.csv') %>% 
#   dplyr::select(JRH_Final_Ref_Cod, StationID, `Collection Date`)
# benthics <- read_csv('data/dataOut/benthics.csv') %>% 
#   left_join(refstress, by = c('StationID', 'Collection Date')) %>% 
#   dplyr::select(JRH_Final_Ref_Cod, everything()) %>% 
#   filter(JRH_Final_Ref_Cod %in% c('Ref', 'Ref?','Ref-Pied?'))
# threshold <- 5
# groupingFields <- c("JRH_Final_Ref_Cod","StationID", "BenSampID" ,"Collection Date", "Collection DateTime", "RepNum" )


removeRareTaxa <- function(benthics, # wide tibble of taxa count by station
                           threshold, # numeric 
                           groupingFields # vector of all non taxa fields
){

  taxaMath <- benthics %>% 
    dplyr::select(-c( !!groupingFields) ) %>% 
    pivot_longer(cols = everything(), names_to = 'Taxa', values_to = 'Count') %>% 
    filter(!is.na(Count)) %>% 
    mutate(Count = ifelse(Count > 0, 1, 0)) %>% 
    group_by(Taxa) %>% 
    summarise(Count = sum(Count)) %>% 
    mutate(Percent = Count / nrow(benthics) * 100)
  # Separate these for easy reporting to user
  taxaToKeep <- filter(taxaMath, Percent >= !! threshold) 
  taxaToDrop <- filter(taxaMath, Percent < !! threshold) 
  
  benthicsSlim <- benthics %>% 
    group_by_at( groupingFields) %>% 
    pivot_longer(cols = -groupingFields, names_to = 'Taxa', values_to = 'Count') %>% 
    filter(Taxa %in% taxaToKeep$Taxa) %>% 
    replace(is.na(.), 0) %>% 
    mutate(Count = log1p(Count)) %>% 
    pivot_wider(names_from = Taxa, values_from = Count)
    
  return(
    list(`Taxa Math` = taxaMath,
         `Taxa Kept` = taxaToKeep,
         `Taxa Dropped` = taxaToDrop,
         `Benthics Slim` = benthicsSlim) )
}

dataList <- removeRareTaxa(benthics, # dataset to analyze, can filter on whatever you want beforehand (e.g. ref, stress, play with ecoregions)
                           threshold = 5, # drop taxa below 5% occurrence in dataset 
                           groupingFields = c("JRH_Final_Ref_Cod", "StationID", "BenSampID" ,
                                              "Collection Date", "Collection DateTime", "RepNum")) # fields that you don't want thrown into the math
# how to view data that comes out of function
View(dataList$`Taxa Math`) # taxa occurrence and percent breakdown
View(dataList$`Taxa Kept`) # taxa higher or equal to input threshold
View(dataList$`Taxa Dropped`) # taxa lower than threshold, dropped from further analyses
View(dataList$`Benthics Slim`) # Taxa below threshold removed, individual counts logged
