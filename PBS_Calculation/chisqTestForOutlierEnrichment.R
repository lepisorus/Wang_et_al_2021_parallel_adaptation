df <- matrix(c(13017, 157641, 96840308, 1962932621), 2, 2, dimnames=list(c("coding", "non-coding"), c("selected", "non-selected")))
df
c <- chisq.test(df)
c
c$expected
c$observd
c$observed
savehistory("chisqTestForOutlierEnrichment.R")
