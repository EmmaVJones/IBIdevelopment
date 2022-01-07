# rewrite of the calcAllBentMets() to use calcBentTaxMetsEVJ() instead of default method that doesn't work with VA master taxa

# calcAllBentMets(indf = benthicsPrep,
#                 inTaxa = masterTaxaListTarget,#masterTaxaList %>% mutate(TAXA_ID = FinalID), 
#                 sampID = "BENSAMPID", #c('UID','SAMPLE_TYPE','SAMPLE_CAT') option to concat more
#                 dist = "IS_DISTINCT", 
#                 ct = "TOTAL", 
#                 taxa_id = "TAXA_ID", 
#                 ffg = "FFG",
#                 habit = "HABIT", 
#                 ptv = "PTV")


# indf = benthicsPrep
# inTaxa = masterTaxaListTarget#masterTaxaList %>% mutate(TAXA_ID = FinalID)
# sampID = "BENSAMPID" #c('UID','SAMPLE_TYPE','SAMPLE_CAT') option to concat more
# dist = "IS_DISTINCT" 
# ct = "TOTAL"
# taxa_id = "TAXA_ID" 
# ffg = "FFG"
# habit = "HABIT"
# ptv = "PTV"


calcAllBentMetsEVJ <- function (indf, inTaxa, sampID = "UID", dist = "IS_DISTINCT",
                                ct = "TOTAL", taxa_id = "TAXA_ID", ffg = "FFG",
                                habit = "HABIT", ptv = "PTV"){
  checkTaxa <- indf[indf$TAXA_ID %in% setdiff(indf$TAXA_ID, 
                                              inTaxa$TAXA_ID), ]
  if (nrow(checkTaxa) > 0) {
    return(print("Taxa in counts that do not have matches in taxalist! Cannot continue."))
  }
  if ("NON_TARGET" %in% names(inTaxa)) {
    inTaxa <- subset(inTaxa, is.na(NON_TARGET) | NON_TARGET == 
                       "" | NON_TARGET == "N")
  }
  ctVars <- c(sampID, dist, ct, taxa_id)
  if (any(ctVars %nin% names(indf))) {
    msgTraits <- which(ctVars %nin% names(indf))
    print(paste("Missing variables in input count data frame:", 
                paste(names(indf)[msgTraits], collapse = ",")))
    return(NULL)
  }
  necTraits <- c("PHYLUM", "CLASS", "ORDER", 
                 "FAMILY", "TRIBE", "SUBFAMILY", "GENUS", 
                 ffg, habit, ptv)
  if (any(necTraits %nin% names(inTaxa))) {
    msgTraits <- which(necTraits %nin% names(inTaxa))
    return(paste("Some of the traits are missing from the taxa list. The following are \nrequired for metric calculations to run:\n", 
                 necTraits[msgTraits], "\n"))
  }
  inTaxa <- subset(inTaxa, select = names(inTaxa) %in% c("TAXA_ID", 
                                                         "PHYLUM", "CLASS", "ORDER", "FAMILY", 
                                                         "TRIBE", "SUBFAMILY", "GENUS", ffg, 
                                                         habit, ptv))
  indf[, c(ct, taxa_id, dist)] <- lapply(indf[, c(ct, taxa_id, 
                                                  dist)], as.numeric)
  inTaxa[, c(ptv, taxa_id)] <- lapply(inTaxa[, c(ptv, taxa_id)], 
                                      as.numeric)
  # EVJ Edit
  taxMet <- calcBentTaxMetsEVJ(indf, inTaxa, sampID, dist, ct, 
                            taxa_id)
  # taxMet <- calcBentTaxMets(indf, inTaxa, sampID, dist, ct, 
  #                           taxa_id)
  tax.1 <- reshape(taxMet, idvar = sampID, direction = "long", 
                   varying = c("AMPHNTAX", "AMPHPIND", "AMPHPTAX", 
                               "CHIRNTAX", "CHIRPIND", "CHIRPTAX", 
                               "CRUSNTAX", "CRUSPIND", "CRUSPTAX", 
                               "DIPTNTAX", "DIPTPIND", "DIPTPTAX", 
                               "EPHENTAX", "EPHEPIND", "EPHEPTAX", 
                               "EPOTNTAX", "EPOTPIND", "EPOTPTAX", 
                               "EPT_NTAX", "EPT_PIND", "EPT_PTAX", 
                               "HEMINTAX", "HEMIPIND", "HEMIPTAX", 
                               "MITENTAX", "MITEPIND", "MITEPTAX", 
                               "MOLLNTAX", "MOLLPIND", "MOLLPTAX", 
                               "NOINNTAX", "NOINPIND", "NOINPTAX", 
                               "ODONNTAX", "ODONPIND", "ODONPTAX", 
                               "OLLENTAX", "OLLEPIND", "OLLEPTAX", 
                               "ORTHNTAX", "ORTHPIND", "ORTHPTAX", 
                               "PLECNTAX", "PLECPIND", "PLECPTAX", 
                               "TANYNTAX", "TANYPIND", "TANYPTAX", 
                               "TRICNTAX", "TRICPIND", "TRICPTAX", 
                               "TUBINAIDNTAX", "TUBINAIDPIND", "TUBINAIDPTAX"#, 
                   ), v.names = "value", timevar = "variable", 
                   # EVJ edit
                              # "ORTHCHIRPIND"), v.names = "value", timevar = "variable", 
                   times = c("AMPHNTAX", "AMPHPIND", "AMPHPTAX", 
                             "CHIRNTAX", "CHIRPIND", "CHIRPTAX", 
                             "CRUSNTAX", "CRUSPIND", "CRUSPTAX", 
                             "DIPTNTAX", "DIPTPIND", "DIPTPTAX", 
                             "EPHENTAX", "EPHEPIND", "EPHEPTAX", 
                             "EPOTNTAX", "EPOTPIND", "EPOTPTAX", 
                             "EPT_NTAX", "EPT_PIND", "EPT_PTAX", 
                             "HEMINTAX", "HEMIPIND", "HEMIPTAX", 
                             "MITENTAX", "MITEPIND", "MITEPTAX", 
                             "MOLLNTAX", "MOLLPIND", "MOLLPTAX", 
                             "NOINNTAX", "NOINPIND", "NOINPTAX", 
                             "ODONNTAX", "ODONPIND", "ODONPTAX", 
                             "OLLENTAX", "OLLEPIND", "OLLEPTAX", 
                             "ORTHNTAX", "ORTHPIND", "ORTHPTAX", 
                             "PLECNTAX", "PLECPIND", "PLECPTAX", 
                             "TANYNTAX", "TANYPIND", "TANYPTAX", 
                             "TRICNTAX", "TRICPIND", "TRICPTAX", 
                             "TUBINAIDNTAX", "TUBINAIDPIND", "TUBINAIDPTAX"))#, 
                            # "ORTHCHIRPIND")) # EVJ edit
  ffgMet <- calcBentFFGmets(indf, inTaxa, sampID, dist, ct, 
                            taxa_id, ffg)
  ffg.1 <- reshape(ffgMet, idvar = sampID, direction = "long", 
                   varying = c("COFINTAX", "COFIPIND", "COFIPTAX", 
                               "COFITRICNTAX", "COFITRICPIND", "COFITRICPTAX", 
                               "COGANTAX", "COGAPIND", "COGAPTAX", 
                               "PREDNTAX", "PREDPIND", "PREDPTAX", 
                               "SCRPNTAX", "SCRPPIND", "SCRPPTAX", 
                               "SHRDNTAX", "SHRDPIND", "SHRDPTAX"), 
                   v.names = "value", timevar = "variable", 
                   times = c("COFINTAX", "COFIPIND", "COFIPTAX", 
                             "COFITRICNTAX", "COFITRICPIND", "COFITRICPTAX", 
                             "COGANTAX", "COGAPIND", "COGAPTAX", 
                             "PREDNTAX", "PREDPIND", "PREDPTAX", 
                             "SCRPNTAX", "SCRPPIND", "SCRPPTAX", 
                             "SHRDNTAX", "SHRDPIND", "SHRDPTAX"))
  habitMet <- calcBentHabitMets(indf, inTaxa, sampID, dist, 
                                ct, taxa_id, habit)
  habit.1 <- reshape(habitMet, idvar = sampID, direction = "long", 
                     varying = c("BURRNTAX", "BURRPIND", "BURRPTAX", 
                                 "CLMBNTAX", "CLMBPIND", "CLMBPTAX", 
                                 "CLNGNTAX", "CLNGPIND", "CLNGPTAX", 
                                 "SPWLNTAX", "SPWLPIND", "SPWLPTAX", 
                                 "SWIMNTAX", "SWIMPIND", "SWIMPTAX"), 
                     v.names = "value", timevar = "variable", 
                     times = c("BURRNTAX", "BURRPIND", "BURRPTAX", 
                               "CLMBNTAX", "CLMBPIND", "CLMBPTAX", 
                               "CLNGNTAX", "CLNGPIND", "CLNGPTAX", 
                               "SPWLNTAX", "SPWLPIND", "SPWLPTAX", 
                               "SWIMNTAX", "SWIMPIND", "SWIMPTAX"))
  tolMet <- calcBentTolMets(indf, inTaxa, sampID, dist, ct, 
                            taxa_id, ptv)
  tol.1 <- reshape(tolMet, idvar = sampID, direction = "long", 
                   varying = c("FACLNTAX", "FACLPIND", "FACLPTAX", 
                               "INTLNTAX", "INTLPIND", "INTLPTAX", 
                               "NTOLNTAX", "NTOLPIND", "NTOLPTAX", 
                               "STOLNTAX", "STOLPIND", "STOLPTAX", 
                               "TL01NTAX", "TL01PIND", "TL01PTAX", 
                               "TL23NTAX", "TL23PIND", "TL23PTAX", 
                               "TL45NTAX", "TL45PIND", "TL45PTAX", 
                               "TL67NTAX", "TL67PIND", "TL67PTAX", 
                               "TOLRNTAX", "TOLRPIND", "TOLRPTAX", 
                               "WTD_TV"), v.names = "value", timevar = "variable", 
                   times = c("FACLNTAX", "FACLPIND", "FACLPTAX", 
                             "INTLNTAX", "INTLPIND", "INTLPTAX", 
                             "NTOLNTAX", "NTOLPIND", "NTOLPTAX", 
                             "STOLNTAX", "STOLPIND", "STOLPTAX", 
                             "TL01NTAX", "TL01PIND", "TL01PTAX", 
                             "TL23NTAX", "TL23PIND", "TL23PTAX", 
                             "TL45NTAX", "TL45PIND", "TL45PTAX", 
                             "TL67NTAX", "TL67PIND", "TL67PTAX", 
                             "TOLRNTAX", "TOLRPIND", "TOLRPTAX", 
                             "WTD_TV"))
  domMet <- calcBentDominMets(indf, inTaxa, sampID, dist, ct, 
                              taxa_id)
  dom.1 <- reshape(domMet, idvar = sampID, direction = "long", 
                   varying = c("HPRIME", "DOM1PIND", "DOM3PIND", 
                               "DOM5PIND", "CHIRDOM1PIND", "CHIRDOM3PIND", 
                               "CHIRDOM5PIND"), v.names = "value", timevar = "variable", 
                   times = c("HPRIME", "DOM1PIND", "DOM3PIND", 
                             "DOM5PIND", "CHIRDOM1PIND", "CHIRDOM3PIND", 
                             "CHIRDOM5PIND"))
  names(indf)[names(indf) == ct] <- "FINAL_CT"
  names(indf)[names(indf) == dist] <- "IS_DISTINCT"
  rhs <- paste(sampID, collapse = "+")
  form <- paste("cbind(FINAL_CT, IS_DISTINCT)", rhs, 
                sep = "~")
  totals <- aggregate(formula(form), data = indf, FUN = function(x) sum = sum(x))
  names(totals)[names(totals) == "FINAL_CT"] <- "TOTLNIND"
  names(totals)[names(totals) == "IS_DISTINCT"] <- "TOTLNTAX"
  totals.long <- reshape(totals, idvar = sampID, direction = "long", 
                         varying = c("TOTLNIND", "TOTLNTAX"), v.names = "value", 
                         timevar = "variable", times = c("TOTLNIND", 
                                                         "TOTLNTAX"))
  mets <- rbind(tax.1, ffg.1, habit.1, tol.1, dom.1, totals.long)
  metOut <- reshape(mets, direction = "wide", idvar = sampID, 
                    timevar = "variable")
  names(metOut) <- gsub("value\\.", "", names(metOut))
  return(metOut)
}


