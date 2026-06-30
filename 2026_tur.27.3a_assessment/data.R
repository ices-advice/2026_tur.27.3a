## Preprocess data, write TAF data tables

## Before: tur.27.3a.InterCatch_raised_discards.Rds, ICES_1950-2010.csv,
##         ICES_2006-2023.csv, Preliminary landings 2024&2025, indexQ1,3,4.csv
## After: tur.27.3a_catch.csv, indexQ1.csv

library(icesTAF)
library(dplyr)

## Helper functions ----
vsum <- function(x, y) {
  sum(x, y, na.rm = TRUE )
}

getBMS <- function(fn) {
  read.csv(fn,
    stringsAsFactors = FALSE,
    header = TRUE,
    fileEncoding = "UTF-8-BOM") %>% 
    as_tibble() %>%
    filter(Species.Latin.Name %in% c("Psetta maxima", "Scophthalmus maximus"),
      Area %in% c("27_3_A","27_3_A_20","27_3_A_21",
                  "27.3.a.20","27.3.a.21","27.3.a")) %>%
    mutate(Country = ifelse(Country == "GB", "UK", Country)) %>%
    {
      if ("ConfidentialityFlag" %in% names(.)) {
        mutate(., ConfidentialityFlag = ifelse(ConfidentialityFlag == "Y", "Y", "N"))
      } else {
        mutate(., ConfidentialityFlag = "N")
      }
    } %>%
    {
      if ("AMS_Catch" %in% names(.)) {
        . %>%
          mutate(
            AMS.Catch.TLW. = AMS_Catch,
            BMS.Catch.TLW. = suppressWarnings(as.numeric(BMS_Catch))
          )
      } else {
        .
      }
    } %>%
    group_by(Year, ConfidentialityFlag) %>%
    summarise(BMS = sum(BMS.Catch.TLW., na.rm = TRUE))
}
    
get_catchcat_percent <- function(dat) {
  remove_subdivision <- function(x) paste0(head(strsplit(x, "[.]")[[1]], -1), collapse = ".")
  perarea <- dat %>%
    filter(! Catch.Cat. %in% c("Logbook Registered Discard", "BMS landing")) %>%
    group_by(Catch.Cat., Area) %>%
    summarise(catchpercatchcat = sum(Catch..kg), .groups = "drop") %>%
    ungroup() %>%
    reshape2::dcast(Area ~ Catch.Cat., value.var = "catchpercatchcat") %>%
    mutate(Total = Discards + Landings,
           `Discard ratio` = ifelse(Total > 0, Discards / (Landings + Discards) * 100, 0),
           DLratio = ifelse(Total > 0, Discards / Landings, 0),
           Year = first(dat$Year))
  total <- perarea %>%
    summarise(Area = "Total", ## remove_subdivision(first(Area)),
              Discards = sum(Discards, na.rm = TRUE),
              Landings = sum(Landings, na.rm = TRUE),
              Total = sum(Total, na.rm = TRUE),
              Year = first(Year),
              .groups = "drop") %>%
    mutate(#`Raising factor` = Discards / Landings * 100,
      DLratio = ifelse(Total > 0, Discards / Landings, 0),
      `Discard ratio` = Discards / (Landings + Discards) * 100)
  bind_rows(perarea, total)
}


#' Reads ICES preliminary catch statistics and returns only
#'
#' @param fn file name, it should contain the
#' @param latin string or vector, latin names of species to select
#' @param areas vector with ICES areas, e.g. "27_3_A_20"
#' @param areaout string, the name of the combined area
#' @param speciesout string, species name in the returned data.frame
#'
#' @return data.frame of selected species by country
#' @export
#'
#' @examples
readPrel <- function(fn,
                     latin = c("Psetta maxima", "Scophthalmus maximus"),
                     areas = c("27_3_A", "27_3_A_20", "27_3_A_21",
                               "27.3.a.20","27.3.a.21","27.3.a"),
                     areaout = "27.3.a",
                     speciesout = "TUR") {
  res <- read.csv(fn, stringsAsFactors = FALSE, header = TRUE,fileEncoding="UTF-8-BOM") %>% as_tibble()
  yr <- unique(res$Year)
  stopifnot(length(yr) == 1)
  if("ConfidentialityFlag" %in% names(res)) {
    res$ConfidentialityFlag <- ifelse(res$ConfidentialityFlag == "Y", "Y", "N")
  } else {
    res$ConfidentialityFlag <- "N"
  }
  if ("AMS_Catch" %in% names(res)) {
    res$AMS.Catch.TLW. <- res$AMS_Catch
    suppressWarnings(res$BMS.Catch.TLW. <- as.numeric(res$BMS_Catch))
  }
  res %>%
    filter(Species.Latin.Name %in% latin,
           Area %in% areas) %>%
    mutate(Country = ifelse(Country == "GB", "UK", Country)) %>%
    group_by(Year, Country,ConfidentialityFlag) %>%
    summarise(Species = speciesout,
              Area = areaout,
              Units = "TLW",
              Landings = vsum(AMS.Catch.TLW., BMS.Catch.TLW.)) %>%
    ungroup() %>%
    transmute(Species, Area, Units, Country = factor(Country, levels = unique(labels), ordered = TRUE), 
              Year = as.integer(yr), Landings = Landings, ConfidentialityFlag = ConfidentialityFlag)
}

## Loading data and settings ----
mkdir("data")

levels <- c("Denmark", "Sweden", "Netherlands", "Norway", "UK (England)",
            "UK(Scotland)", "Germany", "UKS", "UKE", "Belgium",
            "Germany, Fed. Rep. of", "UK - Eng+Wales+N.Irl.",
            "UK - England & Wales", "GB")
labels <- c("DK",      "SE",     "NL",          "NO",     "UK",
            "UK",           "DE",      "UK",  "UK",  "BE",
            "DE",                    "UK",
            "UK",                   "UK")

icd <- readRDS(file = "bootstrap/data/tur.27.3a.InterCatch_raised_discards.Rds")
lastyr <- max(icd$Year)

## Official landings ----
olhist <-
  read.csv("bootstrap/data/ICES_1950-2010.csv",
           stringsAsFactors = FALSE, header = TRUE) %>% as_tibble() %>%
  filter(Species == "Turbot", Division == "III a") %>%
  reshape2:::melt.data.frame(id.vars = c("Species", "Division",  "Country"), variable.name = "Year") %>%
  mutate(Year = gsub("X","",Year))  %>%
  as_tibble() %>%
  filter(as.integer(Year) >= 1950) %>%
  mutate(Year = as.integer(Year),
         Country = ifelse(Country %in% levels, Country, "Other"),
         Country = factor(Country, levels = levels, labels = labels, ordered = TRUE)) %>%
  filter(! value %in% c("-", ".", "<0.5")) %>%
  group_by(Year, Country) %>%
  summarise(value = sum(as.numeric(value), na.rm = TRUE)) %>%
  ungroup %>%
  transmute(Species = "TUR", Area = "27.3.a", Units = "TLW", Country, 
            Year, Landings = value,  ConfidentialityFlag = "N") %>%
  filter(Year < 2006, Landings > 0)
ol <- read.csv("bootstrap/data/ICESCatchDataset2006-2023.csv",
               stringsAsFactors = FALSE, header = TRUE, na.strings = "c",fileEncoding="UTF-8-BOM") %>% as_tibble()

prels <- bind_rows(readPrel("bootstrap/data/Preliminary_landings_allSpecies_2024.csv"),
                   readPrel("bootstrap/data/Preliminary_landings_allSpecies_2025.csv"))

oltur <- ol %>% filter(Species == "TUR", Area %in% c("27.3.a")) %>%
  reshape2:::melt.data.frame(id.vars = c("Species", "Area", "Units", "Country"), variable.name = "Year") %>%
  mutate(Year = gsub("X","",Year))  %>%
  as_tibble() %>%
  mutate(Country = ifelse(Country == "GB", "UK", Country),
         Year = as.integer(Year),
         Country = factor(Country, levels = unique(labels), ordered = TRUE)) %>%
  group_by(Species, Area, Units, Country, Year) %>%
  summarise(Landings = sum(as.numeric(value))) %>%
  mutate(ConfidentialityFlag = "N") %>%
  bind_rows(olhist, prels)

saveRDS(oltur, file = "data/official_landings_tur.27.3a.Rds")

bms <- bind_rows(
  getBMS("bootstrap/data/Preliminary_landings_allSpecies_2024.csv"),
  getBMS("bootstrap/data/Preliminary_landings_allSpecies_2025.csv"))
saveRDS(bms, file = "data/tur.27.3a.BMS.Rds")


## Catches ----

## Discard rate and estimation
percent_discards <- lapply(split(icd, icd$Year), get_catchcat_percent)%>%
  bind_rows() %>%
  filter(Area == "Total") %>%
  reshape2:::melt.data.frame(id.vars = c("Area", "Year"))

## Mean discard rate
dlrate <- mean(filter(percent_discards, Year <= 2018, variable == "DLratio")$value, na.rm = TRUE)


turcatch_mean <- oltur %>% group_by(Year) %>%
  summarise(Landings = sum(Landings)) %>%
  filter(Year <= 2001) %>%
  transmute(Year, Landings, Discards = dlrate * Landings, Total = Landings + Discards) %>%
  reshape2:::melt.data.frame(id.vars = "Year", factorsAsStrings = TRUE) %>%
  as_tibble() %>%
  bind_rows(
    percent_discards %>%
      filter(variable != "Discard ratio") %>%
      mutate(value = value / 1000)
    ##%>% reshape2::dcast(Year~variable) %>%
    ##transmute(Year, Landings = Landings / 1000, Catch = Total / 1000, Discards = Discards / 1000)
  )

save(turcatch_mean, file = "data/tur.27.3a.catch_data.Rdata")


tur_27_3a_catch <- turcatch_mean  %>% reshape2::dcast(Year ~ variable)
from1960 <- tur_27_3a_catch$Year >= 1960
tur_27_3a_catch <- tur_27_3a_catch[from1960, ]

write.taf(tur_27_3a_catch, dir = "data")

tur_27_3a_indexQ1 <- read.taf("bootstrap/data/indexQ1.csv")
tur_27_3a_indexQ3 <- read.taf("bootstrap/data/indexQ3.csv")
tur_27_3a_indexQ4 <- read.taf("bootstrap/data/indexQ4.csv")

write.taf(c("tur_27_3a_indexQ1"), dir = "data")
write.taf(c("tur_27_3a_indexQ3"), dir = "data")
write.taf(c("tur_27_3a_indexQ4"), dir = "data")
