prepBentCts_WSA_EVJ <- function (inCts, inTaxa = bentTaxa_nrsa, sampID = "UID", ct = "TOTAL", 
          taxa_id = "TAXA_ID") 
{
  ctVars <- c(sampID, ct, taxa_id)
  if (any(ctVars %nin% names(inCts))) {
    msgTraits <- which(ctVars %nin% names(inCts))
    print(paste("Missing variables in input data frame:", 
                paste(names(inCts)[msgTraits], collapse = ",")))
    return(NULL)
  }
  inCts <- subset(inCts, select = c(sampID, ct, taxa_id))
  if (is.null(inTaxa)) {
    inTaxa <- bentTaxa
    inTaxa <- subset(inTaxa, is.na(NON_TARGET) | NON_TARGET == 
                       "")
  }
  necTraits <- c("PHYLUM", "CLASS", "ORDER", "FAMILY", "GENUS", 
                 "TARGET_TAXON", taxa_id)
  if (any(necTraits %nin% names(inTaxa))) {
    msgTraits <- which(necTraits %nin% names(inTaxa))
    return(paste("Some of the traits are missing from the taxa list. The following are \nrequired for metric calculations to run:\n", 
                 necTraits[msgTraits], "\n"))
  }
  names(inCts)[names(inCts) == ct] <- "TOTAL"
  names(inCts)[names(inCts) == taxa_id] <- "TAXA_ID"
  names(inTaxa)[names(inTaxa) == taxa_id] <- "TAXA_ID"
  inTaxa.1 <- inTaxa[, c("TAXA_ID", "TARGET_TAXON", "PHYLUM", 
                         "CLASS", "ORDER", "FAMILY", "GENUS")]
  inCts.1 <- merge(inCts, inTaxa.1, by = c("TAXA_ID"))
  inCts.1$TOTAL <- as.numeric(inCts.1$TOTAL)
  inCts.1 <- inCts.1[inCts.1$TOTAL > 0, ]
  fixTaxa <- with(inCts.1, which(CLASS %in% c("ARACHNIDA", 
                                              "POLYCHAETA", "OLIGOCHAETA") & !is.na(FAMILY) & FAMILY != 
                                   ""))
  inCts.1$TARGET_TAXON[fixTaxa] <- inCts.1$FAMILY[fixTaxa]
  inCts.2 <- merge(inCts.1, subset(inTaxa.1, select = c("TAXA_ID", 
                                                        "TARGET_TAXON")), by = "TARGET_TAXON", all.x = TRUE)
  names(inCts.2)[names(inCts.2) == "TAXA_ID.y"] <- "TAXA_ID"
  # inCts.2$TAXA_ID <- with(inCts.2, ifelse(TARGET_TAXON %in% 
  #                                           c("CRICOTOPUS/ORTHOCLADIUS", "THIENEMANNIMYIA GENUS GR."), 
  #                                         3581, ifelse(TARGET_TAXON %in% c("CERATOPOGONINAE"), 
  #                                                      3566, TAXA_ID)))
  # inCts.2$TARGET_TAXON <- with(inCts.2, ifelse(TAXA_ID == 
  #                                                3581, "CHIRONOMIDAE", ifelse(TAXA_ID == 3566, "CERATOPOGONIDAE", 
  #                                                                             TARGET_TAXON)))
  totals <- aggregate(x = list(TOTAL = inCts.2$TOTAL), by = inCts.2[c(sampID, 
                                                                      "TAXA_ID")], FUN = sum)
  inCts.3 <- merge(totals, inTaxa.1, by = "TAXA_ID")
  inCts.4 <- assignDistinct(inCts.3, c(sampID), taxlevels = c("PHYLUM", 
                                                              "CLASS", "ORDER", "FAMILY", "GENUS"), final.name = "TARGET_TAXON")
  inCts.4$IS_DISTINCT <- with(inCts.4, ifelse(is.na(IS_DISTINCT), 
                                              0, IS_DISTINCT))
  outCts <- subset(inCts.4, select = c(sampID, "TAXA_ID", 
                                       "TOTAL", "IS_DISTINCT"))
  return(outCts)
}
