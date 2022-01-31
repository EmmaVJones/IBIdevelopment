# rewrite of the calcBentTaxMets() to skip ORTHOCLADIINAE and TANYTARSINI issues

# outTax <- calcBentTaxMets(benthicsPrep,
#                           masterTaxaListTarget,
#                           sampID = "BENSAMPID", 
#                           dist='IS_DISTINCT',
#                           ct='TOTAL')


# inCts = benthicsPrep
# inTaxa =masterTaxaListTarget#masterTaxaList %>% mutate(TAXA_ID = FinalID)
# sampID = "BENSAMPID" #c('UID','SAMPLE_TYPE','SAMPLE_CAT') option to concat more
# dist = "IS_DISTINCT" 
# ct =  "TOTAL"
# taxa_id = "TAXA_ID" 

calcBentTaxMetsEVJ <- function (inCts, inTaxa, sampID = "UID", dist = "IS_DISTINCT", 
                                ct = "TOTAL", taxa_id = "TAXA_ID") 
{
  ctVars <- c(sampID, dist, ct, taxa_id)
  if (any(ctVars %nin% names(inCts))) {
    msgTraits <- which(ctVars %nin% names(inCts))
    print(paste("Missing variables in input data frame:", 
                paste(names(inCts)[msgTraits], collapse = ",")))
    return(NULL)
  }
  inCts <- subset(inCts, select = c(sampID, ct, dist, taxa_id))
  names(inCts)[names(inCts) == ct] <- "TOTAL"
  names(inCts)[names(inCts) == dist] <- "IS_DISTINCT"
  names(inCts)[names(inCts) == taxa_id] <- "TAXA_ID"
  for (i in 1:length(sampID)) {
    if (i == 1) 
      inCts$SAMPID <- inCts[, sampID[i]]
    else inCts$SAMPID <- paste(inCts$SAMPID, inCts[, sampID[i]], 
                               sep = ".")
  }
  checkTaxa <- inCts[!(inCts$TAXA_ID %in% inTaxa$TAXA_ID), 
  ]
  if (nrow(checkTaxa) > 0) {
    return(print("Taxa in counts that do not have matches in taxalist! Cannot continue."))
  }
  if ("NON_TARGET" %in% names(inTaxa)) {
    inTaxa <- subset(inTaxa, is.na(NON_TARGET) | NON_TARGET == 
                       "" | NON_TARGET == "N")
  }
  necTraits <- c("PHYLUM", "CLASS", "ORDER", 
                 "FAMILY", "TRIBE", "SUBFAMILY", "GENUS")
  if (any(necTraits %nin% names(inTaxa))) {
    msgTraits <- which(necTraits %nin% names(inTaxa))
    return(paste("Some of the traits are missing from the taxa list. The following are \nrequired for metric calculations to run:\n", 
                 necTraits[msgTraits], "\n"))
  }
  samples <- unique(subset(inCts, select = c(sampID, "SAMPID")))
  inCts.1 <- inCts[inCts$TAXA_ID %in% inTaxa$TAXA_ID, c("SAMPID", 
                                                        "TAXA_ID", "TOTAL", "IS_DISTINCT")]
  inCts.1 <- inCts.1[!is.na(inCts.1$TOTAL) & inCts.1$TOTAL > 
                       0, ]
  inTaxa.1 <- inTaxa
  inTaxa.1$EPT_ <- with(inTaxa.1, ifelse(ORDER %in% c("PLECOPTERA", 
                                                      "EPHEMEROPTERA", "TRICHOPTERA"), 1, NA))
  inTaxa.1$EPHE <- with(inTaxa.1, ifelse(ORDER %in% c("EPHEMEROPTERA"), 
                                         1, NA))
  inTaxa.1$PLEC <- with(inTaxa.1, ifelse(ORDER %in% c("PLECOPTERA"), 
                                         1, NA))
  inTaxa.1$TRIC <- with(inTaxa.1, ifelse(ORDER %in% c("TRICHOPTERA"), 
                                         1, NA))
  inTaxa.1$CHIR <- with(inTaxa.1, ifelse(FAMILY %in% c("CHIRONOMIDAE"), 
                                         1, NA))
  inTaxa.1$CRUS <- with(inTaxa.1, ifelse(CLASS %in% c("MALACOSTRACA", 
                                                      "MAXILLOPODA", "BRANCHIOPODA", "CEPHALOCARIDA", 
                                                      "OSTRACODA", "REMIPEDIA"), 1, NA))
  inTaxa.1$NOIN <- with(inTaxa.1, ifelse(CLASS %nin% c("INSECTA"), 
                                         1, NA))
  inTaxa.1$DIPT <- with(inTaxa.1, ifelse(ORDER %in% c("DIPTERA"), 
                                         1, NA))
  inTaxa.1$MOLL <- with(inTaxa.1, ifelse(PHYLUM %in% c("MOLLUSCA"), 
                                         1, NA))
  inTaxa.1$AMPH <- with(inTaxa.1, ifelse(ORDER %in% c("AMPHIPODA"), 
                                         1, NA))
  inTaxa.1$EPOT <- with(inTaxa.1, ifelse(ORDER %in% c("ODONATA", 
                                                      "PLECOPTERA", "EPHEMEROPTERA", "TRICHOPTERA"), 
                                         1, NA))
  inTaxa.1$HEMI <- with(inTaxa.1, ifelse(ORDER %in% c("HEMIPTERA"), 
                                         1, NA))
  inTaxa.1$MITE <- with(inTaxa.1, ifelse(ORDER %in% c("TROMBIDIFORMES", 
                                                      "SARCOPTIFORMES"), 1, NA))
  inTaxa.1$ODON <- with(inTaxa.1, ifelse(ORDER %in% c("ODONATA"), 
                                         1, NA))
  inTaxa.1$OLLE <- with(inTaxa.1, ifelse(CLASS %in% c("OLIGOCHAETA", 
                                                      "HIRUDINEA", "CLITELLATA"), 1, NA))
  inTaxa.1$ORTH <- with(inTaxa.1, ifelse(SUBFAMILY %in% c("ORTHOCLADIINAE"),
                                         1, NA))
  inTaxa.1$TANY <- with(inTaxa.1, ifelse(TRIBE %in% c("TANYTARSINI"), 
                                         1, NA))
  inTaxa.1$TUBINAID <- with(inTaxa.1, ifelse(FAMILY %in% c("TUBIFICIDAE", 
                                                           "NAIDIDAE"), 1, NA))
  if (length(grep("NON_TARGET", names(inTaxa.1))) > 0) {
    inTaxa.1 <- subset(inTaxa.1, is.na(NON_TARGET) | NON_TARGET == 
                         "")
  }
  params <- c("EPT_", "EPHE", "PLEC", "TRIC", 
              "CHIR", "CRUS", "DIPT", "MOLL", 
              "NOIN", "TUBINAID", "AMPH", "EPOT", 
              "HEMI", "MITE", "ODON", "OLLE", 
              "ORTH", # EVJ edit
              "TANY")
  taxalong <- reshape(inTaxa.1[, c("TAXA_ID", params)], 
                      idvar = "TAXA_ID", direction = "long", varying = params, 
                      timevar = "TRAIT", v.names = "value", times = params)
  taxalong <- taxalong[!is.na(taxalong$value), ]
  taxalong$TRAIT <- as.character(taxalong$TRAIT)
  totals <- aggregate(x = list(TOTLNIND = inCts.1$TOTAL, TOTLNTAX = inCts.1$IS_DISTINCT), 
                      by = inCts.1[c("SAMPID")], FUN = sum)
  
  
  
  
  
  
  inCts.1 <- merge(inCts.1, totals, by = "SAMPID")
  inCts.1$CALCPIND <- with(inCts.1, TOTAL/TOTLNIND)
  inCts.1$CALCPTAX <- with(inCts.1, IS_DISTINCT/TOTLNTAX)
  traitDF <- merge(inCts.1, taxalong, by = "TAXA_ID")
  outMet.1 <- aggregate(x = list(NIND = traitDF$TOTAL, NTAX = traitDF$IS_DISTINCT), 
                        by = traitDF[c("SAMPID", "TRAIT", "TOTLNTAX", 
                                       # EVJ edit add total individual for final output
                                       'TOTLNIND')], 
                        FUN = sum)
  outMet.2 <- aggregate(x = list(PIND = traitDF$CALCPIND, PTAX = traitDF$CALCPTAX), 
                        by = traitDF[c("SAMPID", "TRAIT", "TOTLNTAX",
                                       # EVJ edit add total individual for final output
                                       'TOTLNIND'
                                       )], 
                        FUN = function(x) {
                          round(sum(x) * 100, 2)
                        })
  outMet <- merge(outMet.1, outMet.2, by = c("SAMPID", 
                                             "TRAIT", "TOTLNTAX",
                                             # EVJ edit add total individual for final output
                                             'TOTLNIND'
                                             ))
  print("Done calculating taxonomy metrics.")
  outLong <- reshape(outMet, idvar = c("SAMPID", "TOTLNTAX",
                                       # EVJ edit add total individual for final output
                                       'TOTLNIND',
                                       "TRAIT"), direction = "long", varying = names(outMet)[!names(outMet) %in% 
                                                                                               c("SAMPID", "TOTLNTAX", "TRAIT",
                                                                                                 # EVJ edit add total individual for final output
                                                                                                 'TOTLNIND')], 
                     timevar = "variable", v.names = "value", 
                     times = names(outMet)[!names(outMet) %in% c("SAMPID", 
                                                                 "TOTLNTAX", "TRAIT",
                                                                 # EVJ edit add total individual for final output
                                                                 'TOTLNIND')])
  outLong$variable <- paste(outLong$TRAIT, outLong$variable, 
                            sep = "")
  outLong$TRAIT <- NULL
  outWide <- reshape(outLong, idvar = c("SAMPID", "TOTLNTAX",
                                        # EVJ edit add total individual for final output
                                        'TOTLNIND'), 
                     direction = "wide", timevar = "variable", 
                     v.names = "value")
  names(outWide) <- gsub("value\\.", "", names(outWide))
  outWide <- merge(outWide, samples, by = "SAMPID")
  # EVJ Edit
  # outWide$ORTHCHIRPIND <- with(outWide, round(ORTHNIND/CHIRNIND * 
  #                                               100, 2))
  # outWide$ORTHCHIRPIND <- with(outWide, ifelse(is.na(ORTHCHIRPIND) | 
  #                                                is.nan(ORTHCHIRPIND), 0, ORTHCHIRPIND))
  empty_tax <- data.frame(t(rep(NA, 56)), stringsAsFactors = F)
  names(empty_tax) <- c("TOTLNTAX", 
                        
                        # EVJ edit add total individual
                        "TOTLNIND", 
                        
                        "AMPHNTAX", 
                        "AMPHPIND", "AMPHPTAX", "CHIRNTAX", 
                        "CHIRPIND", "CHIRPTAX", "CRUSNTAX", 
                        "CRUSPIND", "CRUSPTAX", "DIPTNTAX", 
                        "DIPTPIND", "DIPTPTAX", "EPHENTAX", 
                        "EPHEPIND", "EPHEPTAX", "EPOTNTAX", 
                        "EPOTPIND", "EPOTPTAX", "EPT_NTAX", 
                        "EPT_PIND", "EPT_PTAX", "HEMINTAX", 
                        "HEMIPIND", "HEMIPTAX", "MITENTAX", 
                        "MITEPIND", "MITEPTAX", "MOLLNTAX", 
                        "MOLLPIND", "MOLLPTAX", "NOINNTAX", 
                        "NOINPIND", "NOINPTAX", "ODONNTAX", 
                        "ODONPIND", "ODONPTAX", "OLLENTAX", 
                        "OLLEPIND", "OLLEPTAX", "ORTHNTAX", 
                        "ORTHPIND", "ORTHPTAX", "PLECNTAX", 
                        "PLECPIND", "PLECPTAX", "TANYNTAX", 
                        "TANYPIND", "TANYPTAX", "TRICNTAX", 
                        "TRICPIND", "TRICPTAX", "TUBINAIDNTAX", 
                        "TUBINAIDPIND", "TUBINAIDPTAX")#, "ORTHCHIRPIND")# EVJ Edit
  outWide.all <- merge(outWide, empty_tax, all = TRUE)
  outWide.all <- subset(outWide.all, !is.na(SAMPID))
  outWide.all[is.na(outWide.all)] <- 0
  outWide.fin <- outWide.all
  outWide.fin$SAMPID <- NULL
  outWide.fin <- outWide.fin[, c(sampID, 
                                 
                                 # EVJ edit, add in total taxa
                                 "TOTLNTAX","TOTLNIND",
                                 
                                 "AMPHNTAX", 
                                 "AMPHPIND", "AMPHPTAX", "CHIRNTAX", 
                                 "CHIRPIND", "CHIRPTAX", "CRUSNTAX", 
                                 "CRUSPIND", "CRUSPTAX", "DIPTNTAX", 
                                 "DIPTPIND", "DIPTPTAX", "EPHENTAX", 
                                 "EPHEPIND", "EPHEPTAX", "EPOTNTAX", 
                                 "EPOTPIND", "EPOTPTAX", "EPT_NTAX", 
                                 "EPT_PIND", "EPT_PTAX", "HEMINTAX", 
                                 "HEMIPIND", "HEMIPTAX", "MITENTAX", 
                                 "MITEPIND", "MITEPTAX", "MOLLNTAX", 
                                 "MOLLPIND", "MOLLPTAX", "NOINNTAX", 
                                 "NOINPIND", "NOINPTAX", "ODONNTAX", 
                                 "ODONPIND", "ODONPTAX", "OLLENTAX", 
                                 "OLLEPIND", "OLLEPTAX", "ORTHNTAX", 
                                 "ORTHPIND", "ORTHPTAX", "PLECNTAX", 
                                 "PLECPIND", "PLECPTAX", "TANYNTAX", 
                                 "TANYPIND", "TANYPTAX", "TRICNTAX", 
                                 "TRICPIND", "TRICPTAX", "TUBINAIDNTAX", 
                                 "TUBINAIDPIND", "TUBINAIDPTAX")]# , "ORTHCHIRPIND")]# EVJ Edit
  return(outWide.fin)
}
  