## Run analysis, write model results

## Before: tur.27.3a_catch.csv, tur.27.3a_index.csv
## After: tur.27.3a_fit.Rds

library(icesTAF)
library(spict)
library(TMB)

mkdir("model")

icesTAF::msg("Model: read in data")

## Load catch and index data
indexQ1 <- read.taf("data/tur_27_3a_indexQ1.csv")
indexQ1 <- subset(indexQ1,!(is.na(indexQ1$Index)))
indexQ3 <- read.taf("data/tur_27_3a_indexQ3.csv")
indexQ3 <- subset(indexQ3,!(is.na(indexQ3$Index))) #rm years with missing data
indexQ4 <- read.taf("data/tur_27_3a_indexQ4.csv")
indexQ4 <- subset(indexQ4,!(is.na(indexQ4$Index))) #rm years with missing data
catch <- read.taf("data/tur_27_3a_catch.csv")

## What is the last year in the data?
lastYear <- max(catch$Year)

## Management period and evaluation time
maninterval <- seq(lastYear + 2, lastYear + 3)
maneval <- lastYear + 3

## Priors
logsdi <- list(c(log(mean(indexQ1$sdlogI)),0.5,1),
               c(log(mean(indexQ3$sdlogI)),0.5,1),
               c(log(mean(indexQ4$sdlogI)),0.5,1)) # estimate of index noise

priors <- list(
  logbeta  = c(0,0,0),
  logalpha = c(0,0,0),
  logsdi = logsdi, # index noise 
  logsdc = c(log(0.1), 0.5, 1), # usually catch noise is low
  logsdb = c(log(0.15), 2, 1), # value based on meta study of Berg et al. (in prep.)
  logr = c(log(0.5),0.5,1), # based on Turbot estimate in Fishbase 
  logbkfrac = c(log(0.75),1, 1)
)

baseinp <- inp <- list(
  timeC = catch$Year,
  obsC  = catch$Total,
  timeI = list(indexQ1$Year,indexQ3$Year+0.5,indexQ4$Year+0.75),
  obsI  = list(indexQ1$Index,indexQ3$Index,indexQ4$Index),
  stdevfacI = list(indexQ1$sdlogI / mean(indexQ1$sdlogI),
                   indexQ3$sdlogI / mean(indexQ3$sdlogI),
                   indexQ4$sdlogI / mean(indexQ4$sdlogI)),
  stdevfacC = c(rep(3, length(1960:2001)), rep(1, length(2002:lastYear))),
  priors = priors,
  maninterval = maninterval,
  maneval = maneval,
  optimiser.control = list(iter.max = 1e4, eval.max = 1e4),
  nspinup = 160
)

add.thorson.gamma.prior<-function(inp){
  n.est <- 1.478
  sdn <- 0.849
  x90 <- qnorm(0.9,n.est,sdn)
  sr <- modefrac2shaperate(log(n.est),log(x90))
  
  inp$priors$logn <- c(log(2), 1, 0)
  inp$priors$logngamma <- c(sr[1], sr[2], 1)
  inp
}

inp <- add.thorson.gamma.prior(inp)
baseinp <- inp

inp <- check.inp(inp)

icesTAF::msg("Model: model fit")
fit <- fit.spict(inp)
fit <- calc.osa.resid(fit)
fit <- calc.process.resid(fit)

### #################################

fitExtraFoptions <- fit

icesTAF::msg("Model: management scenarios")
fit <- add.man.scenario(fit, "F=Fmsy_C_fractile", fractiles = list(catch = 0.35), breakpointB = c(1/2))
fit <- add.man.scenario(fit, "F=Fmsy", breakpointB = c(1/2))
fit <- add.man.scenario(fit, "F=Fsq", ffac = 1)
fit <- add.man.scenario(fit, "F=0", ffac = 0)
fit <- add.man.scenario(fit, "F=Fmsy_All_fractiles", fractiles = list(catch = 0.35, bbmsy = 0.35, ffmsy = 0.35), breakpointB = c(1/2))

icesTAF::msg("Model: retro")
fit <- retro(fit)

icesTAF::msg("Model: saving the results")
saveRDS(fit, "model/tur.27.3a_fit.Rds")

icesTAF::msg("Model: extra Foption management scenarios")
## F options from 0.01 to upper 95% bound of Fmsy estimate
for (fopt in seq(0.01, get.par("Fmsy", fit)[3], 0.01 )) {
  fitExtraFoptions <- add.man.scenario(fitExtraFoptions, paste0("F=", fopt), fractiles = list(catch = 0.35), fabs = fopt)
}

saveRDS(fitExtraFoptions, "model/tur.27.3a_fit_extraFoptions.Rds")


icesTAF::msg("Model: correct retro")
correctRetro <- fit 

# Define quarters and input files
RI1 <- readRDS("bootstrap/data/retroindexQ1.Rds")
RI3 <- readRDS("bootstrap/data/retroindexQ3.Rds")
RI4 <- readRDS("bootstrap/data/retroindexQ4.Rds")

# Use lapply to generate retro inputs from the 3 retroindex lists
retroinps <- lapply(seq_along(RI1), function(j) {
  retroinp <- baseinp
  
  # Add 3 indices (Q1, Q3, Q4)
  retroinp$timeI <- list(RI1[[j]]$Year, RI3[[j]]$Year, RI4[[j]]$Year)
  retroinp$obsI  <- list(RI1[[j]]$Index, RI3[[j]]$Index, RI4[[j]]$Index)
  retroinp$stdevfacI <- list(
    RI1[[j]]$sdlogI / mean(RI1[[j]]$sdlogI),
    RI3[[j]]$sdlogI / mean(RI3[[j]]$sdlogI),
    RI4[[j]]$sdlogI / mean(RI4[[j]]$sdlogI)
  )
  
  # Restrict catch data to retrospective period
  keep <- retroinp$timeC %in% c(1960:1982, RI1[[j]]$Year)
  retroinp$timeC <- retroinp$timeC[keep]
  retroinp$obsC  <- retroinp$obsC[keep]
  
  # stdev for catch
  lstyr <- max(RI1[[j]]$Year)
  retroinp$stdevfacC <- c(rep(3, length(1960:2001)), rep(1, length(2002:lstyr)))
  
  # Clean previous manual inputs
  retroinp$maneval <- NULL
  retroinp$maninterval <- NULL
  
  check.inp(retroinp)
  retroinp
})

# Add the base model as the first item
inpretro <- c(list(baseinp), retroinps)

# Fit all models
correctRetro$retro <- lapply(inpretro, fit.spict)
saveRDS(correctRetro, file = "model/tur.27.3a_fit_correctRetro.Rds")

# icesTAF::msg("Model: extra run using sdI and sdB priors")
# from1975 <- catch$Year >= 1975
# inp <- list(
#   timeC = catch$Year[from1975],
#   obsC  = catch$Total[from1975],
#   timeI = index$Year,
#   obsI  = index$Index,
#   stdevfacI = index$sdlogI,
#   stdevfacC = c(rep(3, length(1975:2001)), rep(1, length(2002:lastYear))),
#   priors = list(
#     logn = c(0,0,0),
#     logalpha = c(0,0,0),
#     logsdi = c(log(0.5), 0.5, 1),
#     logsdb = c(log(0.5), 0.5, 1),
#     logbkfrac = c(logbkfrac, 0.5, 1)
#   ),
#   maninterval = maninterval,
#   maneval = maneval,
#   ini = list(logn = log(2)),
#   phases = list(logn = -1),
#   optimiser.control = list(iter.max = 1e4, eval.max = 1e4)
# )
#
# inp3c <- check.inp(inp)
# fit3c <- fit.spict(inp3c)
# fit3c <- calc.osa.resid(fit3c)
# fit3c <- retro(fit3c)
# fit3c <- add.man.scenario(fit3c, "F=Fmsy_C_fractile", fractiles = list(catch = 0.35), breakpointB = c(1/2))
# saveRDS(fit3c, "model/tur.27.3a_fit_2021_priors_for_sdb_sdi.Rds")

