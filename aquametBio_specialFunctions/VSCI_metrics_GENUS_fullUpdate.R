# These are the metric functions that build the VSCI, converted to genus level for IBI work

# R version 3.6.2
# Adapted from code by Lou Reynolds and Leah Ettema (USEPA)

#x <- bugTraits

#metrich = taxa richness of each metric by sample
# Drop exclued taxa for richness metrics
metrich <- function(x){ x %>%
    #filter(metric_val > 0) %>%
    filter(metric_val > 0 & `Final Genus Level Taxa` > 0) %>%
    group_by(BenSampID,metric)%>%
    summarise(metrich=n_distinct(GVSCI_FinalID))}


#metsum = number of individuals per metric by sample 
# Keep Excluded Taxa for this metric
#x <- bugTraits
metsum <- function(x){ x %>%
    filter(metric_val > 0) %>%
    #filter(metric_val > 0 & `Excluded Taxa` == 0) %>%
    group_by(BenSampID,metric)%>%
    summarise(metsum=sum(`Genus Level Count`))}

#x <- bugdatrare
bugtotals_function <- function(x){ x %>% 
    group_by(BenSampID) %>% 
    summarise(bugtotal=sum(`Genus Level Count`)) }


#bugtotals <- bugtotals_function(bugdatrare)

#metprop=proprotion of individuals in the sample with that metric  
metprop_function <- function(x,bugtotals){
  mets <- full_join(metrich(x),metsum(x),by=c('BenSampID','metric'))
  left_join(mets, bugtotals, by="BenSampID") %>% 
    mutate(metprop=metsum/bugtotal)}
#metprop <- metprop_function(bugTraits, bugtotals)


#IBI Metric Calculations----
#Because these are Family based, the bugdata must be joined to the EDAS list and converted to Family before calculations

#y <- bugdatrare

#The first 5 metrics decrease with stress
#**Total Richness----
rich <- function(y, specialName){ 
  z <- y %>%
    group_by(BenSampID) %>%
    filter(`Genus Level Excluded Taxa` >= 0) %>%
    summarise(totTaxa = n_distinct(GVSCI_FinalID))
  names(z) <- c('BenSampID',specialName)
  return(z)} # Richness
#rich(bugdatrare, 'Genus Total Taxa')

#metricName <- 'e'
#percent <- FALSE
#specialName <- 'Genus EPT Taxa'

summaryStress <- function(metprop, metricName, percent, specialName){
  z <- filter_(metprop, lazyeval::interp(~metric == metricName)) %>%
    group_by(BenSampID) 
  if(percent == FALSE){
    justSummary <- select(z, BenSampID, metrich)
    names(justSummary) <- c('BenSampID',ifelse(is.na(specialName),metricName,specialName))
    return(justSummary)
  }else{
    percent <- mutate(z, percent = metprop*100) %>%
      select(BenSampID, percent)
    names(percent) <- c('BenSampID',paste('%',specialName,sep = ''))
    return(percent)
  }
}
#summaryStress(metprop,'ept', percent = F, NA)
#summaryStress(metprop,'ept', percent = F, 'Genus EPT Taxa')
#summaryStress(metprop,'e', percent = T, 'Ephem')
#summaryStress(metprop, 'ptmin', T, 'PT - Hydro')
#summaryStress(metprop,'scraper', T, 'FamilyScraper')
#summaryStress(metprop,'chiro', T, 'Chironomid')



# ** % Dominant 2 Taxa---- *** note this is programmed differently than % Dom 5 metric (on purpose)
pickTopTwoDom <- function(y){
  #y %>%
  filter(y, `Genus Level Excluded Taxa` >= 0 ) %>%
    group_by(BenSampID,GVSCI_FinalID) %>% 
    summarise(twodtot=sum(`Genus Level Count`)) %>%
    top_n(n=2, wt=twodtot) %>% #returns top n rows.  but if 2nd taxa is tied, returns all rows
    ungroup() %>% 
    group_by(BenSampID) %>%
    arrange(BenSampID, desc(twodtot)) %>% #These two arrange lines put the family in descending order
    slice(1:2)} #incase there is tie drop extra taxa, dom 5 metrics keep everyone!

pDom2 <- function(y, bugtotals, specialName){
  top2 <- pickTopTwoDom(y)
  z <- left_join(top2, bugtotals, by ="BenSampID") %>% 
    mutate(pdom=(twodtot/bugtotal)*100) %>%   
    select(BenSampID,pdom) %>%
    group_by(BenSampID) %>% 
    summarise(pdom2=sum(pdom))
  names(z) <- c('BenSampID',specialName)
  return(z)}
#pDom2(bugdatrare, bugtotals,'Family %2 Dominant')

# **Hisenhoff Index----

#Hilsenhoff Index taxa info
#hiltax <- filter(vmast, metric == 'TolVal') %>%
#  group_by(GVSCI_FinalID) %>%
#  summarise(TolVal = mean(metric_val))


#Raw Calculations for Hilsenhoff Index and Shannon Diversity  
hilsindex <- function(y, hiltax, specialName){
  z <- left_join(y, hiltax, by="GVSCI_FinalID") %>%
    group_by(BenSampID)%>% 
    mutate(nxa=`Genus Level Count`*TolVal, sumn=sum(`Genus Level Count`)) %>% 
    filter(!is.na(TolVal) & `Final Genus Level Taxa` > 0) %>% 
    group_by(BenSampID) %>% 
    summarise(hilsindex=sum(nxa)/sum(`Genus Level Count`))#summarise(hilsindex=sum(nxa)/mean(sumn))
  names(z) <- c('BenSampID',specialName)
  return(z)
}

#hilsindex(y,hiltax,'Genus HBI')

#This function truncates scores to fall between 0-100
truncfxn <- function(x){ifelse(x>=100,100,x)}
