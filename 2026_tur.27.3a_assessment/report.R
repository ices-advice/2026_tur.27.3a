## Prepare plots and tables for report

## Before:
## After:

library(icesTAF)
library(spict)
library(TMB)
library(dplyr)
library(rmarkdown)

mkdir("report")
cp("bootstrap/initial/report/*", "report/")

outdir <- "report/"

output_format <- NULL # "all"
quiet <- FALSE

icesTAF::msg("Report: Making catch working document")
render("report/tur.27.3a_catch_WD.Rmd", output_dir = outdir,
       #output_format = output_format,
       clean = TRUE, quiet = quiet,  encoding = 'UTF-8')


  icesTAF::msg("Report: Making assessment working document")
render("report/tur.27.3a_assessment_WD.Rmd", output_dir = outdir,
       #output_format = output_format,
       clean = TRUE, quiet = quiet,  encoding = 'UTF-8')

)icesTAF::msg("Report: Making catch and assessment presentation")
render("report/tur.27.3a_assessment_Presentation.Rmd",
       output_dir = outdir,
       #output_format = output_format,
       clean = TRUE, quiet = quiet,  encoding = 'UTF-8')

# -------------------------------------------------------------------------------
# extra Foptions too big to run in Rmarkdown
fitExtraFoptions <- readRDS("model/tur.27.3a_fit_extraFoptions.Rds")
assessmentsummary <- read.csv("report/tur.27.3a_assessment_summary.csv")
ly <- 2025
lastTAC <- 395
meandis <- mean((assessmentsummary$Discards / assessmentsummary$Catches * 100) [assessmentsummary$Year %in% seq(ly-2, ly)])

fit3b <- readRDS("model/tur.27.3a_fit.Rds")
bind <- which(fit3b$inp$time == ly+2)

bf <- readRDS("report/tur.27.3a_stockstatus.Rds")
b_for <- exp(bf$BBmsy[bind, 2])

getmanline <- function(fit) {
  cfy <- fit$inp$maninterval[1]
  by <- fit$inp$maneval
  find <- which(fit$inp$time == cfy) 
  bind <- which(fit$inp$time == by)
  cind <- which(fit$inp$timeCpred == cfy)
  
  ct <- get.par("logCpred", fit, TRUE)[cind,2]
  dis <- meandis * ct / 100
  lan <- (1 - meandis  / 100) * ct
  f <- get.par("logFFmsy", fit, TRUE)[find,2]
  if (f < 0.000001) f <- 0
  b <- get.par("logBBmsy", fit, TRUE)[bind,2]
  bchange <- (b-b_for)/b_for*100
  
  nms <- gsub("YYY", by,
              gsub("XXX", cfy, c("Total catch (XXX)", "Proj. land (XXX)", 
                                 "Proj. disc (XXX)", "Fmort (FXXX/ FMSY)",
                                 "Stock size (BYYY/ BMSY)", "% B ch.", "% Advice ch.")))
  setNames(data.frame(round(ct), round(lan), round(dis), 
                      icesRound(f), icesRound(b), icesRound(bchange),
                      icesRound((ct - lastTAC) / lastTAC * 100), row.names = FALSE),
           nms)
  
}

library(dplyr)
library(spict)
library(icesAdvice)
catchscenariotableExtra <- t(sapply(fitExtraFoptions$man, getmanline)) %>% 
  write.csv("report/tur.27.3a_ExtraFoptions.csv", row.names = FALSE)


