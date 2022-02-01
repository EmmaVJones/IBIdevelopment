
# Tolerance and BCG attribute info
BCG <- read_excel('data/MASTER_ATTRIBUTES_BUGS_06062019.xlsx', sheet = 'Sheet1') %>% 
  dplyr::select(FinalID, GenusFinal, BCGattribute = `Virginia General BCG Attribute   First Choice`, 
                `BCG DisOxy`:`BCG pctIMP`)

masterTaxaGenus <- read_excel('data/masterTaxaGenus_01132022.xlsx', sheet = 'masterTaxaGenus') %>% # Emma update
  mutate_at(c("FamTolVal", "TolVal", "GVSCI_TolVal"), as.numeric) %>% 
  left_join(BCG, by = c('OldFinalID' = 'FinalID', 'GVSCI_FinalID' = 'GenusFinal'))

### USE FROM OTHER SIDE OF APP
vmast <- masterTaxaGenus %>%
  mutate(e=ifelse(Order=="Ephemeroptera", 1, 0),
         p=ifelse(Order=="Plecoptera",1,0),
         t=ifelse(Order=="Trichoptera", 1, 0),
         tmin=ifelse((Order=="Trichoptera" & Family != "Hydropsychidae") | 
                       (Order=="Trichoptera" & is.na(Family)) , 1, 0), 
         ept=ifelse(e+p+t>0,1,0), 
         scraper = ifelse(GVSCI_FFG=="Scraper", 1, 0),
         chiro = ifelse(Family=="Chironomidae",1, 0),
         ptmin = ifelse(p + tmin > 0,1,0),
         `clinger-HS` = ifelse(GVSCI_Habit == 'Clinger' & ! Family %in% c("Hydropsychidae","Simuliidae"), 1, 0),
         
         # new for IBI development
         c = ifelse(Genus == 'Cheumatopsyche', 1, 0),
         `ept-h` = ifelse(ept == 1 & Family != 'Hydropsychidae', 1, 0),
         `ept-h+c` = ifelse(`ept-h` > 0 | c > 0, 1 ,0),
         
         elmid = ifelse(Family == 'Elmidae', 1, 0),
         
         `e-b` = ifelse(Order=="Ephemeroptera" & Family != "Baetidae", 1 , 0),
         
         ept4 = ifelse(ept == 1 & GVSCI_TolVal <= 4, 1, 0),
         
         coll = ifelse(GVSCI_FFG %in% c("Collector", "Filterer"), 1, 0) ,
         
         # BCG metrics
         att23 = ifelse(BCGattribute %in% c(2,3), 1, 0),
         att5 = ifelse(BCGattribute %in% c(5), 1, 0)
         
         # stressor specific metrics? Jason request: what's the threshold?
         
         ) %>%
  # Then put it in long format so it can be merged to and input taxa list
  dplyr::select(GVSCI_FinalID, GVSCI_TolVal, e,p,t, ept,ptmin, scraper, chiro,`clinger-HS`,
         `ept-h+c`, elmid, `e-b`, ept4, coll,
         att23, att5) %>% 
  distinct(GVSCI_FinalID, .keep_all = T) %>% # drop multiple rows bc working back to family level data from genus
  filter(!is.na(GVSCI_FinalID)) %>%
  pivot_longer(-GVSCI_FinalID, names_to = 'metric', values_to = 'metric_val') %>%
  #  pivot_longer(-`Final VA Family ID`, names_to = 'metric', values_to = 'metric_val') %>%
  filter(!is.na(metric_val))



edas_options_genus <- select(masterTaxaGenus, Class, Subclass, Order, Suborder, Superfamily, Family, GVSCI_FinalID, FinalID) %>%
  mutate(across(where(is.factor), as.character))
edas_list_genus <- select(edas_options_genus, GVSCI_FinalID, FinalID)
# for Excluding taxa, need a list of all Family level designations that may end up as a FinalID
# these are all unique Family names and the few taxa that are the only 
GenusNames <- c(unique(edas_options_genus$FinalID)[!is.na(unique(edas_options_genus$FinalID))])
FamilyNames <- unique(edas_options_genus$Family)[!is.na(unique(edas_options_genus$Family))]
SuperfamilyNames <- unique(edas_options_genus$Superfamily)[!is.na(unique(edas_options_genus$Superfamily))]
SuborderNames <- unique(edas_options_genus$Suborder)[!is.na(unique(edas_options_genus$Suborder))]
OrderNames <- unique(edas_options_genus$Order)[!is.na(unique(edas_options_genus$Order))]
SubclassNames <- unique(edas_options_genus$Subclass)[!is.na(unique(edas_options_genus$Subclass))]
ClassNames <- unique(edas_options_genus$Class)[!is.na(unique(edas_options_genus$Class))]


EDASrare <- read.csv('data/aquametBioDataOut/benthics01132022.csv') %>%
  rename(`Excluded Taxa` = Excluded.Taxa) %>% 
  
  ########## #filter(str_detect(BenSampID, 'R110') & RepNum == 1) %>% # keep only rarified data and Rep1's
  mutate(Count = Individuals) %>% # Rename to match formatting of functions
  ######`Excluded Taxa` = ifelse(`Excluded Taxa` == T, -1, 0)) %>% 
  select(BenSampID, FinalID, Count, `Excluded Taxa`) %>%
  mutate(GenusTaxaLevel = ifelse(FinalID %in% GenusNames, T, F),
         FamilyTaxaLevel = ifelse(FinalID %in% FamilyNames, T, F),
         SuperfamilyTaxaLevel = ifelse(FinalID %in% SuperfamilyNames, T, F),
         SuborderTaxaLevel = ifelse(FinalID %in% SuborderNames, T, F),
         OrderTaxaLevel = ifelse(FinalID %in% OrderNames, T, F),
         SubclassTaxaLevel = ifelse(FinalID %in% SubclassNames, T, F),
         ClassTaxaLevel = ifelse(FinalID %in% ClassNames, T, F))

# Work FinalID back up to Genus Level
EDASrare2 <- left_join(EDASrare,edas_list_genus, by="FinalID") %>%
  #filter(!is.na(`Final VA Family ID`)) %>%
  rename( `Genus Level Excluded Taxa` = `Excluded Taxa`)


# We also need to do a little data manipulation to incorporate biologist exclusion information appropriately.
exclusionMath  <- EDASrare2 %>%
  mutate(`Genus Level Excluded Taxa` = 
           ifelse(`Genus Level Excluded Taxa` == -1, 
                  ifelse(`SuperfamilyTaxaLevel` == TRUE | `SuborderTaxaLevel` == TRUE | `OrderTaxaLevel` == TRUE | 
                           `SubclassTaxaLevel` == TRUE | `ClassTaxaLevel` == TRUE |
                           `FamilyTaxaLevel` == TRUE , -1, 0), 0 )) %>%
  # had to get super ifelse nesty here to make this logic work, ugly but works
  group_by(BenSampID, GVSCI_FinalID) %>%
  summarise(`Genus Level Count` = sum(Count), 
            #`Genus Level Excluded Taxa` = sum(`Genus Level Excluded Taxa`),
            `Genus Level Taxa` = n(),
            `Genus Level Excluded Taxa` = sum(`Genus Level Excluded Taxa`),
            `Final Genus Level Taxa` = `Genus Level Taxa` + sum(`Genus Level Excluded Taxa`) )

# Join bug traits
bugTraits <- left_join(exclusionMath,vmast,by='GVSCI_FinalID') 


bugdatrare <- exclusionMath


IBImetricCalculation <- function(bugTraits,bugdatrare,vmast){
  bugtotals <- bugtotals_function(bugdatrare)
  metprop <- metprop_function(bugTraits, bugtotals)
  hiltax <- filter(vmast, metric == 'GVSCI_TolVal') %>%
    group_by(GVSCI_FinalID) %>%
    summarise(TolVal = mean(metric_val))
  
  # IBI Calculations
  IBImetrics <- suppressMessages(
    rich(bugdatrare, 'Genus Total Taxa') %>%
      left_join(summaryStress(metprop,'ept', percent = F, 'Genus EPT Taxa')) %>%
      left_join(summaryStress(metprop,'e', percent = T, 'Ephem')) %>%
      left_join(summaryStress(metprop, 'ptmin', T, 'PT - Hydropsychidae')) %>% 
      left_join(summaryStress(metprop,'scraper', T, 'GenusScraper')) %>% 
      left_join(summaryStress(metprop,'chiro', T, 'Chiro')) %>%
      left_join(pDom2(bugdatrare, bugtotals,'Genus %2 Dominant')) %>%
      left_join(hilsindex(bugdatrare,hiltax,'Genus HBI')) %>% 
      # new metrics
      left_join(summaryStress(metprop,'ept-h+c', percent = T, 'EPT-H+C')) %>%
      left_join(summaryStress(metprop,'elmid', percent = F, 'Elmid')) %>%
      left_join(summaryStress(metprop,'e-b', percent = T, 'Ephem-B')) %>%
  #    left_join(summaryStress(metprop,'ept4', percent = T, 'EPT-H+C')) %>%  # this won't work # / % EPT Taxa (with tol 0-4) what???
      left_join(summaryStress(metprop,'att23', percent = F, 'BCGatt2&3')) %>%
      left_join(summaryStress(metprop,'att23', percent = T, 'BCGatt2&3')) %>%
      left_join(summaryStress(metprop,'att5', percent = F, 'BCGatt5')) %>%
      left_join(summaryStress(metprop,'att5', percent = T, 'BCGatt5')) %>%
      
      
      replace(is.na(.), 0) )
  return(IBImetrics)
}


