# organize WVA bug data for analyses
# tested them running through metrics and that went well. For efficiency, I'll add these in to benthics run to avoid doing everything twice
benthicsWVA <- read_csv('data/dataOut/benthicsVA_WVA.csv') %>% 
  filter(DataSource == 'WVA_BUG') %>% 
  mutate(BenSampID = UID) %>% 
  dplyr::select(-c(UID, `Collection DateTime`, RepNum, Season, Year, JulianDate, SampleN)) %>% 
  group_by(DataSource, StationID, BenSampID, `Collection Date`) %>% 
  pivot_longer(-c('DataSource', 'StationID', 'BenSampID', 'Collection Date'), names_to = 'FinalID', values_to = 'Individuals') %>% 
  filter(Individuals > 0) %>% 
  rename(BENSAMPID = BenSampID) %>%  # assignDistinct() changes names to caps so adjust here to avoid issues later
  mutate(FinalID = as.character(FinalID)) %>% 
  mutate_if(is.character, str_to_upper) %>%    # convert names to all caps
  # convert names to newest taxa list
  mutate(FinalID = case_when(FinalID == 'NIXE' ~ 'AFGHANURUS',
                             FinalID == 'ANCHYTARSUS' ~ 'ANCHYTARSUS BICOLOR',
                             FinalID == 'STENONEMA' ~ 'STENONEMA FEMORATUM',
                             TRUE ~ FinalID))



# join to master taxa list to see what FinalID's need to be converted
benthicsWVAjoin <- left_join(benthicsWVA, masterTaxaList, by = 'FinalID')
benthicsWVAjoinIssues <- filter(benthicsWVAjoin, is.na(PHYLUM))



# Now bring in master taxa list and join for metrics

masterTaxaList <- read_excel('data/masterTaxaGenus_01132022_JRH.xlsx',  sheet = 'masterTaxaGenusJRHFeb10') %>%  #change with genus fixes
  #read_excel('data/masterTaxaGenus_01132022.xlsx', sheet = 'masterTaxaGenus') %>% # Emma update
  mutate_at(c("FamTolVal", "TolVal", "GVSCI_TolVal"), as.numeric) %>% 
  rename(FinalID_FFG = FFG) %>% # fix this before things get confusing
  
  # Optioservus fix
  mutate(Genus = case_when(FinalID == 'Promoresia' ~ 'Promoresia',
                           TRUE ~ as.character(Genus))) %>% 
  # Stenonema fix
  mutate(Genus = case_when(FinalID %in% c("Maccaffertium/Stenonema", "Maccaffertium") ~ "Maccaffertium",
                           TRUE ~ as.character(Genus)),
         GVSCI_TolVal = case_when(FinalID %in% c("Maccaffertium/Stenonema", "Maccaffertium") ~ 5.4,
                                  TRUE ~ as.numeric(GVSCI_TolVal)) ) %>% #,
  # GregVSCI_TolVal2 = case_when(FinalID %in% c("Maccaffertium/Stenonema", "Maccaffertium") ~ 5.4,
  #                   TRUE ~ as.numeric(GregVSCI_TolVal2))) %>% 
  
  
  # rename to desired format
  rename(PHYLUM = Phylum, 
         CLASS = Class,
         ORDER = Order, 
         FAMILY = Family,
         SUBFAMILY = Subfamily, 
         TRIBE = Tribe, 
         GENUS = Genus,
         # use genus level bugs for metrics
         #TARGET_TAXON = GVSCI_FinalID,
         FFG = GVSCI_FFG,
         HABIT = GVSCI_Habit, 
         PTV = GVSCI_TolVal) %>% 
  #mutate(TAXA_ID = 1:n()) %>%  # give numeric unique identifier
  mutate_if(is.character, str_to_upper) %>% # package needs uppercase taxa names
  
  # recode FFG to desired format: CF=collector-filterer, CG=collector-gatherer,PA=parasite,PI=Piercer,PR=predator,SC=scraper,SH=shredder
  mutate(FFG = case_when(FFG == "COLLECTOR" ~ "CG", 
                         FFG == "SCRAPER" ~ "SC",
                         FFG == "SHREDDER" ~ "SH",
                         FFG == "PREDATOR" ~ "PR",
                         FFG == "FILTERER" ~ "CF",
                         TRUE ~ as.character(NA)),
         #  recode HABIT to desired format:  AT=attached,BU=burrower,CB=climber, CN=clinger,DV=diver,PK=planktonic,SK=skater,SP=sprawler,SW=swimmer    
         HABIT = case_when(HABIT == "SWIMMER" ~ "SW", 
                           HABIT == "CLINGER" ~ "CN",
                           HABIT ==  "CLIMBER" ~ "CB",
                           HABIT == "BURROWER" ~ "BU",
                           HABIT == "SPRAWLER" ~ "SP",
                           TRUE ~ as.character(NA))) 

masterTaxaList[ masterTaxaList == "NA" ] <- NA

masterTaxaListTarget <- masterTaxaList %>% 
  #filter(Target == 1) # this will drop all higher level taxa from joining, not a good idea
  arrange(GVSCI_FinalID, Target) %>% 
  
  filter(!is.na(GVSCI_FinalID)) %>% 
  
  
  distinct(GVSCI_FinalID, .keep_all = T) %>% 
  rename(TARGET_TAXON = GVSCI_FinalID) %>% 
  mutate(TAXA_ID = 1:n()) %>% 
  as.data.frame()
#rowwise() %>% 
#mutate(TAXA_ID = sum(3000, TAXA_ID))

benthicsWVA <- benthicsWVA %>%
  #rename(BENSAMPID = BenSampID) %>%  # assignDistinct() changes names to caps so adjust here to avoid issues later
  #mutate(FinalID = as.character(FinalID)) %>% 
  #mutate_if(is.character, str_to_upper) %>%   # convert names to all caps
  
  ## join to full master taxa list
  #left_join(dplyr::select(masterTaxaList, FinalID, TARGET_TAXON, TAXA_ID),  by = "FinalID")
  left_join(dplyr::select(masterTaxaList, FinalID, GVSCI_FinalID#TARGET_TAXON),  by = "FinalID")
  ),  by = "FinalID") %>% 
  # join to target level master taxa list for unique TAXA_ID
  left_join(dplyr::select(masterTaxaListTarget, TARGET_TAXON, TAXA_ID),  by = c("GVSCI_FinalID" = 'TARGET_TAXON')) #"FinalID"
# left_join(dplyr::select(masterTaxaListTarget, FinalID, TARGET_TAXON, TAXA_ID),  by = "FinalID")



benthicsWVAPrep <- prepBentCts_WSA_EVJ(
  inCts = benthicsWVA,
  inTaxa = masterTaxaListTarget,
  sampID = "BENSAMPID",
  ct = "Individuals",
  taxa_id ="TAXA_ID"# "FinalID"#"TAXA_ID"
)

outTaxWVA <- calcBentTaxMetsEVJ(benthicsWVAPrep,
                          masterTaxaListTarget,
                          sampID = "BENSAMPID",
                          dist='IS_DISTINCT',
                          ct='TOTAL')
outTolWVA <- calcBentTolMets(benthicsWVAPrep,
                          masterTaxaListTarget,
                          sampID="BENSAMPID",
                          dist='IS_DISTINCT',
                          ct='TOTAL',
                          ptv='PTV')

ffgMetWVA <- calcBentFFGmets(benthicsWVAPrep,
                          masterTaxaListTarget,
                          sampID="BENSAMPID",
                          dist='IS_DISTINCT',
                          ct='TOTAL')

habitMetWVA <- calcBentHabitMets(benthicsWVAPrep,
                              masterTaxaListTarget,
                              sampID="BENSAMPID",
                              dist='IS_DISTINCT',
                              ct='TOTAL',
                              habit = "HABIT")

domMetWVA <- calcBentDominMets(benthicsWVAPrep,
                            masterTaxaListTarget,
                            sampID="BENSAMPID",
                            dist='IS_DISTINCT',
                            ct='TOTAL')
# Bonus VA metrics
source('aquametBio_specialFunctions/IBIspecialRequestFunctions_IsDistinct.R')

BCG <- read_excel('data/MASTER_ATTRIBUTES_BUGS_06062019.xlsx', sheet = 'Sheet1') %>% 
  dplyr::select(FinalID, GenusFinal, BCGattribute = `Virginia General BCG Attribute   First Choice`, 
                `BCG DisOxy`:`BCG pctIMP`) %>% 
  mutate_at(c("BCGattribute", "BCG DisOxy", "BCG acidity", "BCG Alkalinity", "BCG spCond", "BCG Chloride", 
              "BCG Sulfate", "BCG TN.TP",  "BCG totHab", "BCG RBS", "BCG pctIMP"), as.numeric) %>% 
  mutate(BCGattribute = case_when(BCGattribute == 0 ~ as.numeric(NA),
                                  TRUE ~ BCGattribute))

masterTaxaListBCG <- read_excel('data/masterTaxaGenus_01132022.xlsx', sheet = 'masterTaxaGenus') %>% # Emma update
  mutate_at(c("FamTolVal", "TolVal", "GVSCI_TolVal"), as.numeric) %>% 
  rename(FinalID_FFG = FFG) %>% # fix this before things get confusing
  
  
  # Optioservus fix
  mutate(Genus = case_when(FinalID == 'Promoresia' ~ 'Promoresia',
                           TRUE ~ as.character(Genus))) %>% 
  
  # Stenonema fix
  mutate(Genus = case_when(FinalID %in% c("Maccaffertium/Stenonema", "Maccaffertium") ~ "Maccaffertium",
                           TRUE ~ as.character(Genus)),
         GVSCI_TolVal = case_when(FinalID %in% c("Maccaffertium/Stenonema", "Maccaffertium") ~ 5.4,
                                  TRUE ~ as.numeric(GVSCI_TolVal)) ) %>% #,
  # GregVSCI_TolVal2 = case_when(FinalID %in% c("Maccaffertium/Stenonema", "Maccaffertium") ~ 5.4,
  #                   TRUE ~ as.numeric(GregVSCI_TolVal2))) %>% 
  # 
  
  
  
  # add BCG information
  left_join(BCG, by = c('OldFinalID' = 'FinalID', 'GVSCI_FinalID' = 'GenusFinal')) %>% 
  
  
  # rename to desired format
  rename(PHYLUM = Phylum, 
         CLASS = Class,
         ORDER = Order, 
         FAMILY = Family,
         SUBFAMILY = Subfamily, 
         TRIBE = Tribe, 
         GENUS = Genus,
         # use genus level bugs for metrics
         #TARGET_TAXON = GVSCI_FinalID,
         FFG = GVSCI_FFG,
         HABIT = GVSCI_Habit, 
         PTV = GVSCI_TolVal) %>% 
  #mutate(TAXA_ID = 1:n()) %>%  # give numeric unique identifier
  mutate_if(is.character, str_to_upper) %>% # package needs uppercase taxa names
  
  # recode FFG to desired format: CF=collector-filterer, CG=collector-gatherer,PA=parasite,PI=Piercer,PR=predator,SC=scraper,SH=shredder
  mutate(FFG = case_when(FFG == "COLLECTOR" ~ "CG", 
                         FFG == "SCRAPER" ~ "SC",
                         FFG == "SHREDDER" ~ "SH",
                         FFG == "PREDATOR" ~ "PR",
                         FFG == "FILTERER" ~ "CF",
                         TRUE ~ as.character(NA)),
         #  recode HABIT to desired format:  AT=attached,BU=burrower,CB=climber, CN=clinger,DV=diver,PK=planktonic,SK=skater,SP=sprawler,SW=swimmer    
         HABIT = case_when(HABIT == "SWIMMER" ~ "SW", 
                           HABIT == "CLINGER" ~ "CN",
                           HABIT ==  "CLIMBER" ~ "CB",
                           HABIT == "BURROWER" ~ "BU",
                           HABIT == "SPRAWLER" ~ "SP",
                           TRUE ~ as.character(NA))) 


masterTaxaListBCG[ masterTaxaListBCG == "NA" ] <- NA

masterTaxaListBCGTarget <- masterTaxaListBCG %>% 
  #filter(Target == 1) # this will drop all higher level taxa from joining, not a good idea
  arrange(GVSCI_FinalID, Target) %>% 
  distinct(GVSCI_FinalID, .keep_all = T) %>% 
  rename(TARGET_TAXON = GVSCI_FinalID) %>% 
  mutate(TAXA_ID = 1:n()) %>% 
  as.data.frame()

bonusMetricsWVA <- IBImetrics(benthicsWVAPrep, masterTaxaListBCGTarget)
