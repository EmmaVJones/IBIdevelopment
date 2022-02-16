# IBI special request functions, new starting point with IS_DISTINCT from aquametBio
# R 4.1.2
# Tolerance and BCG attribute info
# BCG <- read_excel('data/MASTER_ATTRIBUTES_BUGS_06062019.xlsx', sheet = 'Sheet1') %>% 
#   dplyr::select(FinalID, GenusFinal, BCGattribute = `Virginia General BCG Attribute   First Choice`, 
#                 `BCG DisOxy`:`BCG pctIMP`) %>% 
#   mutate_at(c("BCGattribute", "BCG DisOxy", "BCG acidity", "BCG Alkalinity", "BCG spCond", "BCG Chloride", 
#               "BCG Sulfate", "BCG TN.TP",  "BCG totHab", "BCG RBS", "BCG pctIMP"), as.numeric) %>% 
#   mutate(BCGattribute = case_when(BCGattribute == 0 ~ as.numeric(NA),
#                                   TRUE ~ BCGattribute))
# 
# masterTaxaListBCG <- read_excel('data/masterTaxaGenus_01132022.xlsx', sheet = 'masterTaxaGenus') %>% # Emma update
#   mutate_at(c("FamTolVal", "TolVal", "GVSCI_TolVal"), as.numeric) %>% 
#   rename(FinalID_FFG = FFG) %>% # fix this before things get confusing
#   
#   # add BCG information
#   left_join(BCG, by = c('OldFinalID' = 'FinalID', 'GVSCI_FinalID' = 'GenusFinal')) %>% 
#   
#   
#   # rename to desired format
#   rename(PHYLUM = Phylum, 
#          CLASS = Class,
#          ORDER = Order, 
#          FAMILY = Family,
#          SUBFAMILY = Subfamily, 
#          TRIBE = Tribe, 
#          GENUS = Genus,
#          # use genus level bugs for metrics
#          #TARGET_TAXON = GVSCI_FinalID,
#          FFG = GVSCI_FFG,
#          HABIT = GVSCI_Habit, 
#          PTV = GVSCI_TolVal) %>% 
#   #mutate(TAXA_ID = 1:n()) %>%  # give numeric unique identifier
#   mutate_if(is.character, str_to_upper) %>% # package needs uppercase taxa names
#   
#   # recode FFG to desired format: CF=collector-filterer, CG=collector-gatherer,PA=parasite,PI=Piercer,PR=predator,SC=scraper,SH=shredder
#   mutate(FFG = case_when(FFG == "COLLECTOR" ~ "CG", 
#                          FFG == "SCRAPER" ~ "SC",
#                          FFG == "SHREDDER" ~ "SH",
#                          FFG == "PREDATOR" ~ "PR",
#                          FFG == "FILTERER" ~ "CF",
#                          TRUE ~ as.character(NA)),
#          #  recode HABIT to desired format:  AT=attached,BU=burrower,CB=climber, CN=clinger,DV=diver,PK=planktonic,SK=skater,SP=sprawler,SW=swimmer    
#          HABIT = case_when(HABIT == "SWIMMER" ~ "SW", 
#                            HABIT == "CLINGER" ~ "CN",
#                            HABIT ==  "CLIMBER" ~ "CB",
#                            HABIT == "BURROWER" ~ "BU",
#                            HABIT == "SPRAWLER" ~ "SP",
#                            TRUE ~ as.character(NA))) 
#   
# 
# masterTaxaListBCG[ masterTaxaListBCG == "NA" ] <- NA
# 
# masterTaxaListBCGTarget <- masterTaxaList %>% 
#   #filter(Target == 1) # this will drop all higher level taxa from joining, not a good idea
#   arrange(GVSCI_FinalID, Target) %>% 
#   distinct(GVSCI_FinalID, .keep_all = T) %>% 
#   rename(TARGET_TAXON = GVSCI_FinalID) %>% 
#   mutate(TAXA_ID = 1:n()) %>% 
#   as.data.frame()




# fieldName <- 'BCGattribute'
# metricName <- 'att23'
# includedAttributes <- c(2,3)

BCGmath <- function(BCGorg, fieldName, metricName, includedAttributes){
  BCGorg %>% as_tibble() %>% 
    group_by(BENSAMPID) %>% 
    filter(IS_DISTINCT > 0) %>% # drop non-distinct taxa
    dplyr::select(BENSAMPID:TARGET_TAXON, filterField = !! fieldName) %>% 
    filter(!is.na(filterField)) %>% # drop rows without attribution
    mutate(total = n()) %>% 
    filter(filterField %in% includedAttributes) %>% 
    summarise(metric = unique(metricName),
              metrich = n(),
              total = unique(total),
              metprop = metrich / total)
}
#BCGmath(BCGorg, 'BCGattribute', 'att5', c(5))


#x <- benthicsPrepTaxa
bugtotals_function_NEW <- function(x){ x %>% 
    group_by(BENSAMPID) %>% 
    distinct(TAXA_ID, .keep_all = T) %>% 
    summarise(bugtotal=sum(TOTAL)) }
#bugtotals <- bugtotals_function_NEW(benthicsPrepTaxa)

#metrich = taxa richness of each metric by sample
# Drop exclued taxa for richness metrics
metrichNEW <- function(x){ x %>%
    filter(metric_val > 0 & IS_DISTINCT > 0) %>%
    group_by(BENSAMPID,metric)%>%
    summarise(metrich=n_distinct(TAXA_ID))}


#metsum = number of individuals per metric by sample 
# Keep Excluded Taxa for this metric
#x <- benthicsPrepTaxa
metsumNEW <- function(x){ x %>%
    filter(metric_val > 0) %>%
    group_by(BENSAMPID,metric)%>%
    summarise(metsum=sum(TOTAL))}


#metprop=proprotion of individuals in the sample with that metric  
metprop_function_NEW <- function(x,bugtotals){
  mets <- full_join(metrichNEW(x),metsumNEW(x),by=c('BENSAMPID','metric'))
  left_join(mets, bugtotals, by="BENSAMPID") %>% 
    mutate(metprop=metsum/bugtotal)}
#metprop <- metprop_function_NEW(benthicsPrepTaxa, bugtotals)


# y <- benthicsPrepTaxa

#The first 5 metrics decrease with stress
#**Total Richness----
richNEW <- function(y, specialName){ 
  y %>%
    group_by(BENSAMPID) %>%
    filter(IS_DISTINCT >= 1) %>%
    summarise(totTaxa = n_distinct(TAXA_ID)) %>% 
    rename(!! quo_name(specialName) := totTaxa ) }# Richness
#richNEW(benthicsPrepTaxa, 'Genus Total Taxa')


#metricName <- 'e'
#percent <- FALSE
#specialName <- 'Genus EPT Taxa'

summaryStressNEW <- function(metpropNEW, metricName, percent, specialName){
  z <- filter_(metpropNEW, lazyeval::interp(~metric == metricName)) %>%
    group_by(BENSAMPID) 
  if(percent == FALSE){
    justSummary <- select(z, BENSAMPID, metrich)
    names(justSummary) <- c('BENSAMPID',ifelse(is.na(specialName),metricName,specialName))
    return(justSummary)
  }else{
    percent <- mutate(z, percent = metprop*100) %>%
      select(BENSAMPID, percent)
    names(percent) <- c('BENSAMPID',paste('%',specialName,sep = ''))
    return(percent)
  }
}
#summaryStressNEW(metpropNEW,'ept', percent = F, NA)

# ** % Dominant 2 Taxa---- *** note this is programmed differently than % Dom 5 metric (on purpose)
pickTopTwoDomNEW <- function(y){
  #y %>%
  filter(y, IS_DISTINCT >= 1 ) %>%
    group_by(BENSAMPID, TAXA_ID) %>% 
    summarise(twodtot=sum(TOTAL)) %>%
    top_n(n=2, wt=twodtot) %>% #returns top n rows.  but if 2nd taxa is tied, returns all rows
    ungroup() %>% 
    group_by(BENSAMPID) %>%
    arrange(BENSAMPID, desc(twodtot)) %>% #These two arrange lines put the family in descending order
    slice(1:2)} #incase there is tie drop extra taxa, dom 5 metrics keep everyone!

pDom2NEW <- function(y, bugtotalsNEW, specialName){
  top2 <- pickTopTwoDomNEW(y)
  left_join(top2, bugtotalsNEW, by ="BENSAMPID") %>% 
    mutate(pdom=(twodtot/bugtotal)*100) %>%   
    select(BENSAMPID,pdom) %>%
    group_by(BENSAMPID) %>% 
    summarise(pdom2=sum(pdom)) %>% 
    rename(!! quo_name(specialName) := pdom2 ) }
#pDom2NEW(benthicsPrep, bugtotalsNEW,'Genus %2 Dominant')


# **Hisenhoff Index----

#Hilsenhoff Index taxa info
#hiltax <- filter(vmast, metric == 'TolVal') %>%
#  group_by(GVSCI_FinalID) %>%
#  summarise(TolVal = mean(metric_val))

#y <-  benthicsPrep

#Raw Calculations for Hilsenhoff Index and Shannon Diversity  
hilsindexNEW <- function(y, hiltaxNEW, specialName){
  left_join(y, hiltaxNEW, by="TAXA_ID") %>%
    group_by(BENSAMPID) %>% 
    mutate(nxa = TOTAL * TolVal, 
           sumn = sum(TOTAL)) %>% 
    filter(!is.na(TolVal) & IS_DISTINCT > 0) %>% 
    group_by(BENSAMPID) %>% 
    summarise(hilsindex = sum(nxa)/sum(TOTAL)) %>% 
    rename(!! quo_name(specialName) := hilsindex ) 
}
#hilsindexNEW(y,hiltax,'Genus HBI')


IBImetrics <- function(benthicsPrep, masterTaxaListBCGTarget){
  # taxa math
  vmastNEW <- masterTaxaListBCGTarget %>%
    mutate(e=ifelse(ORDER== "EPHEMEROPTERA", 1, 0),
           p=ifelse(ORDER=="PLECOPTERA",1,0),
           t=ifelse(ORDER=="TRICHOPTERA", 1, 0),
           tmin=ifelse((ORDER=="TRICHOPTERA" & TARGET_TAXON != 'HYDROPSYCHE') , 1, 0),  #View(filter(masterTaxaListBCGTarget, str_detect(TARGET_TAXON, 'HYDROPS')))
           ept=ifelse(e+p+t>0,1,0), 
           scraper = ifelse(FFG=="SC", 1, 0),
           chiro = ifelse(FAMILY=="CHIRONOMIDAE",1, 0),
           ptmin = ifelse(p + tmin > 0,1,0),
           `clinger-HS` = ifelse(HABIT == 'CN' & TARGET_TAXON != 'HYDROPSYCHE' & FAMILY != "SIMULIIDAE", 1, 0),
           
           # new for IBI development
           `ept-h+c` = ifelse(ept == 1 & TARGET_TAXON != "HYDROPSYCHE" & GENUS != 'CHEUMATOPSYCHE', 1, 0), # removing genus == Cheumatopsyche doesn't matter since it is already filtered with family hydropsychidae & Genus != 'Cheumatopsyche', 1, 0),
           elmid = ifelse(FAMILY == "ELMIDAE", 1, 0),
           `e-b` = ifelse(ORDER=="EPHEMEROPTERA" & FAMILY != "BAETIDAE", 1 , 0),
           `ept-h-b` = ifelse(ept == 1 & TARGET_TAXON != 'HYDROPSYCHE' & FAMILY != "BAETIDAE", 1, 0),
           coll = ifelse(FFG %in% c("CG", "CF"), 1, 0) ) %>%
    
    # BCG metrics are calculated outside this bc the n taxa attributed (denominator for percent metrics) depends on each individual sample
          
    # Then put it in long format so it can be merged to and input taxa list
    dplyr::select(TAXA_ID, TARGET_TAXON, PTV, e,p,t, ept,ptmin, scraper, chiro,`clinger-HS`,
                  `ept-h+c`, elmid, `e-b`, `ept-h-b`, coll) %>%
    distinct(TARGET_TAXON, .keep_all = T) %>% # drop multiple rows bc working back to family level data from genus
    filter(!is.na(TARGET_TAXON)) %>%
    pivot_longer(-c(TARGET_TAXON, TAXA_ID), names_to = 'metric', values_to = 'metric_val') %>%
    #  pivot_longer(-`Final VA Family ID`, names_to = 'metric', values_to = 'metric_val') %>%
    filter(!is.na(metric_val))
  
  
  benthicsPrepTaxa <- left_join(benthicsPrep, vmastNEW, by = 'TAXA_ID')
  
  # BCG Math
  BCGorg <- left_join(benthicsPrep, 
                      dplyr::select(masterTaxaListBCGTarget, TAXA_ID, TARGET_TAXON, contains('BCG')),
                      by = c('TAXA_ID'))
  
  
  
  bugtotalsNEW <- bugtotals_function_NEW( benthicsPrepTaxa )
  metpropNEW <- metprop_function_NEW(benthicsPrepTaxa, bugtotalsNEW)
  hiltaxNEW <- filter(vmastNEW, metric == 'PTV') %>%
    group_by(TAXA_ID, TARGET_TAXON) %>%
    summarise(TolVal = mean(metric_val))
  
  IBImetrics <- suppressMessages(
    richNEW(benthicsPrepTaxa, 'Genus Total Taxa') %>%
      left_join(summaryStressNEW(metpropNEW,'ept', percent = F, 'Genus EPT Taxa')) %>%
      left_join(summaryStressNEW(metpropNEW,'e', percent = T, 'Ephem')) %>%
      left_join(summaryStressNEW(metpropNEW, 'ptmin', T, 'PT - Hydropsychidae')) %>% 
      left_join(summaryStressNEW(metpropNEW,'scraper', T, 'GenusScraper')) %>% 
      left_join(summaryStressNEW(metpropNEW,'chiro', T, 'Chiro')) %>%
      left_join(pDom2NEW(benthicsPrep, bugtotalsNEW,'Genus %2 Dominant')) %>%
      left_join(hilsindexNEW(benthicsPrep,hiltaxNEW,'Genus HBI')) %>% 
      # new metrics
      left_join(summaryStressNEW(metpropNEW,'ept-h+c', percent = T, 'EPT-H+C')) %>%
      left_join(summaryStressNEW(metpropNEW,'elmid', percent = F, 'Elmid')) %>%
      left_join(summaryStressNEW(metpropNEW,'e-b', percent = T, 'Ephem-B')) %>%
      # long form here since there is a variable denominator
      left_join(summaryStressNEW(left_join(benthicsPrep, masterTaxaListBCGTarget, by = 'TAXA_ID') %>% 
                                   group_by(BENSAMPID) %>% 
                                   filter(!is.na(PTV)) %>% 
                                   mutate(total = n()) %>% 
                                   filter(ORDER %in% c("EPHEMEROPTERA", "PLECOPTERA", "TRICHOPTERA")) %>% 
                                   filter(PTV <= 4.5) %>% 
                                   summarise(metric = 'ept4.5',
                                             metrich = n(),
                                             total = total,
                                             metprop = metrich / total) %>% # why isn't this summarizing??
                                   distinct(BENSAMPID, .keep_all = T), 'ept4.5', percent = T, 'EPT4.5')) %>% 
      # long form here since there is a variable denominator
      left_join(summaryStressNEW(left_join(benthicsPrep, masterTaxaListBCGTarget, by = 'TAXA_ID') %>% 
                                   group_by(BENSAMPID) %>% 
                                   filter(!is.na(PTV)) %>% 
                                   mutate(total = n()) %>% 
                                   filter(ORDER %in% c("EPHEMEROPTERA", "PLECOPTERA", "TRICHOPTERA")) %>% 
                                   filter(PTV <= 6.5) %>% 
                                   summarise(metric = 'ept6.5',
                                             metrich = n(),
                                             total = total,
                                             metprop = metrich / total) %>% # why isn't this summarizing??
                                   distinct(BENSAMPID, .keep_all = T), 'ept6.5', percent = T, 'EPT6.5')) %>% 
      left_join(summaryStressNEW(BCGmath(BCGorg, 'BCGattribute', 'BCGatt2&3', c(2,3)), 'BCGatt2&3', percent = F, 'BCGatt2&3')) %>%
      left_join(summaryStressNEW(BCGmath(BCGorg, 'BCGattribute', 'BCGatt2&3', c(2,3)), 'BCGatt2&3', percent = T, 'BCGatt2&3')) %>%
      left_join(summaryStressNEW(BCGmath(BCGorg, 'BCGattribute', 'att5', c(5)), 'att5', percent = F, 'BCGatt5')) %>%
      left_join(summaryStressNEW(BCGmath(BCGorg, 'BCGattribute', 'att5', c(5)), 'att5', percent = T, 'BCGatt5')) %>%
      left_join(summaryStressNEW(BCGmath(BCGorg, "BCG DisOxy", 'BCGatt2&3', c(2,3)), 'BCGatt2&3', percent = T, 'BCG_DO_att2&3')) %>%
      left_join(summaryStressNEW(BCGmath(BCGorg, "BCG DisOxy", 'BCGatt5', c(5)), 'BCGatt5', percent = T, 'BCG_DO_att5')) %>%
      left_join(summaryStressNEW(BCGmath(BCGorg, "BCG acidity", 'BCGatt2&3', c(2,3)), 'BCGatt2&3', percent = T, 'BCG_Acidity_att2&3')) %>%
      left_join(summaryStressNEW(BCGmath(BCGorg, "BCG acidity", 'BCGatt5', c(5)), 'BCGatt5', percent = T, 'BCG_Acidity_att5')) %>%
      left_join(summaryStressNEW(BCGmath(BCGorg, "BCG Alkalinity", 'BCGatt2&3', c(2,3)), 'BCGatt2&3', percent = T, 'BCG_Alkalinity_att2&3')) %>%
      left_join(summaryStressNEW(BCGmath(BCGorg, "BCG Alkalinity", 'BCGatt5', c(5)), 'BCGatt5', percent = T, 'BCG_Alkalinity_att5')) %>%
      left_join(summaryStressNEW(BCGmath(BCGorg, "BCG spCond", 'BCGatt2&3', c(2,3)), 'BCGatt2&3', percent = T, 'BCG_spCond_att2&3')) %>%
      left_join(summaryStressNEW(BCGmath(BCGorg, "BCG spCond", 'BCGatt5', c(5)), 'BCGatt5', percent = T, 'BCG_spCond_att5')) %>%
      left_join(summaryStressNEW(BCGmath(BCGorg, "BCG Chloride", 'BCGatt2&3', c(2,3)), 'BCGatt2&3', percent = T, 'BCG_Chloride_att2&3')) %>%
      left_join(summaryStressNEW(BCGmath(BCGorg, "BCG Chloride", 'BCGatt5', c(5)), 'BCGatt5', percent = T, 'BCG_Chloride_att5')) %>%
      left_join(summaryStressNEW(BCGmath(BCGorg, "BCG Sulfate" , 'BCGatt2&3', c(2,3)), 'BCGatt2&3', percent = T, 'BCG_Sulfate_att2&3')) %>%
      left_join(summaryStressNEW(BCGmath(BCGorg, "BCG Sulfate" , 'BCGatt5', c(5)), 'BCGatt5', percent = T, 'BCG_Sulfate_att5')) %>%
      left_join(summaryStressNEW(BCGmath(BCGorg, "BCG TN.TP" , 'BCGatt2&3', c(2,3)), 'BCGatt2&3', percent = T, 'BCG_TN.TP_att2&3')) %>%
      left_join(summaryStressNEW(BCGmath(BCGorg, "BCG TN.TP" , 'BCGatt5', c(5)), 'BCGatt5', percent = T, 'BCG_TN.TP_att5')) %>%
      left_join(summaryStressNEW(BCGmath(BCGorg, "BCG totHab" , 'BCGatt2&3', c(2,3)), 'BCGatt2&3', percent = T, 'BCG_totHab_att2&3')) %>%
      left_join(summaryStressNEW(BCGmath(BCGorg, "BCG totHab" , 'BCGatt5', c(5)), 'BCGatt5', percent = T, 'BCG_totHab_att5')) %>%
      left_join(summaryStressNEW(BCGmath(BCGorg, "BCG RBS" , 'BCGatt2&3', c(2,3)), 'BCGatt2&3', percent = T, 'BCG_RBS_att2&3')) %>%
      left_join(summaryStressNEW(BCGmath(BCGorg, "BCG RBS" , 'BCGatt5', c(5)), 'BCGatt5', percent = T, 'BCG_RBS_att5')) %>%
      left_join(summaryStressNEW(BCGmath(BCGorg, "BCG pctIMP"  , 'BCGatt2&3', c(2,3)), 'BCGatt2&3', percent = T, 'BCG_pctIMP_att2&3')) %>%
      left_join(summaryStressNEW(BCGmath(BCGorg, "BCG pctIMP"  , 'BCGatt5', c(5)), 'BCGatt5', percent = T, 'BCG_pctIMP_att5')) )
   
  return(IBImetrics)
}



