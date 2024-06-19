library(here)
library(data.table)
library(glue)
library(googlesheets4)


## NOTE new sheet
## create an ID to access the googlesheets results sheet
yourl <- "https://docs.google.com/spreadsheets/d/1FU7OH9HXDwFS48qYqmoRTbQpGr29Kvx6AeC7FTbK_tk/edit#gid=0"
shid <- as.character(as_sheets_id(yourl))


## utility function
upload.to.sheets <- function(basename,filename,sheetid
                             ){
  filename <- gsub("\\.csv$","",filename) #safety in case csv included at and
  fn <- glue(basename) + filename + ".csv"
  tmp <- fread(file=fn)
  write_sheet(tmp,sheetid,sheet=filename)
}

## upload relevant table data
upload.to.sheets(here('outdata//'),'cascadetab',shid) #first will need to re-authenticate

# ## rest can be run as block

## need article tables
yurl <- "https://docs.google.com/spreadsheets/d/1jLlXsVAPUhwNP2gMhBVy74wfthREhMdTxT8oPtjCNSU/edit#gid=1897282126"
shidneat <- as.character(as_sheets_id(yurl))

## ---- Table 1 -------
## build 1st
# load(here('outdata/Table1ATT.Rdata'))
# 
# Table1 <- Table1ATT
# 
# write_sheet(Table1,shidneat,sheet="Tab1RAW")

flz1 <- c(
  "Table1ATT.csv", "Table1ATTscreen.csv", "Table1ATTscreenint.csv")
for( fn in flz1)
  upload.to.sheets(here('outdata//'),fn,shidneat)

## ---- Table 2 -------

## ATT part
load(here('outdata/Table2ATT_formatted.Rdata'))
write_sheet(Table2ATT_formatted,shidneat,sheet="Tab2ATT")

load(here('outdata/Table2ATT_formattedscreen.Rdata'))
write_sheet(Table2ATT_formatted,shidneat,sheet="Tab2ATTscreen")

load(here('outdata/Table2ATT_formattedscreenINT.Rdata'))
write_sheet(Table2ATT_formatted,shidneat,sheet="Tab2ATTscreenINT")

# Appendix tables 
# cascades, costs etc
flz1 <- c(
  'cascade_table.csv', 'country_unit_costs.csv')
for( fn in flz1)
  upload.to.sheets(here('costs//outdata//'),fn,shidneat)

## Table S3
TableS2 <- fread(here('outdata/hiv_prevalence.csv'))
write_sheet(TableS2,shidneat,sheet="TableS2")

TableS3 <- fread(here('outdata/TableS3.csv'))
write_sheet(TableS3,shidneat,sheet="TableS3")

load(here('outdata/PD.Rdata'))
write_sheet(PD,shidneat,sheet="TableS4")

# age sub-groups
load(here('outdata/Table2ATTY_formatted.Rdata'))
write_sheet(Table2ATTY_formatted,shidneat,sheet="Table2ATTY")

load(here('outdata/Table2ATTO_formatted.Rdata'))
write_sheet(Table2ATTO_formatted,shidneat,sheet="Table2ATTO")

## --- SA table ---
sa.base <- fread(here('outdata/ICERall.csv'))
sa.hi <- fread(here('outdata/ICERallhi.csv'))
sa.lo <- fread(here('outdata/ICERalllo.csv'))
sa.effects <- fread(here('outdata/ICERallpooled_eff.csv'))
# sa.cfr <- fread(here('outdata/ICERallcfr.csv'))
sa.cfrltfu <- fread(here('outdata/ICERallcfrltfu.csv'))
sa.cfrwho <- fread(here('outdata/ICERallwhocfr.csv'))
sa.cfrsystRev <- fread(here('outdata/ICERallsystRev.csv'))
sa.cfrsystRevtxd <- fread(here('outdata/ICERallsystRevtxd.csv'))
sa.screen <- fread(here('outdata/ICERallscreen.csv'))
sa.screenINT <- fread(here('outdata/ICERallscreenINT.csv'))

names(sa.base)[3] <- 'Base case'
names(sa.hi)[3] <- '5% discount rate'
names(sa.lo)[3] <- '0% discount rate'
names(sa.effects)[3] <- 'Pooled intervention effects'
names(sa.cfrltfu)[3] <- 'CFR including LTFU'
names(sa.cfrwho)[3] <- 'WHO CFR'
names(sa.cfrsystRev)[3] <- 'Systematic review CFR'
names(sa.cfrsystRevtxd)[3] <- 'Systematic review CFR with improvements in treatment outcome'

# names(sa.screen)[3] <- 'Tuberculosis screening rates observed in TIPPI study under SoC'
# names(sa.screenINT)[3] <- 'Tuberculosis screening for all children'
names(sa.screen)[3] <- 'Lower tuberculosis screening rates (-20%)'
names(sa.screenINT)[3] <- 'Higher tuberculosis screening rates (+20%)'

SAll <- Reduce(merge,list(sa.base,sa.lo,sa.hi,sa.cfrltfu,sa.cfrwho,sa.cfrsystRev,sa.cfrsystRevtxd,sa.effects, sa.screen,sa.screenINT))
SAll <- SAll[,iso3:=NULL]
SAll <- SAll |>
  dplyr::mutate_if(is.integer, as.character)
SAll <- melt(SAll, id.vars = 'country', variable.name = 'Assumption', value.name = 'ICER (US$/DALY averted')
setkey(SAll, country)
write_sheet(SAll,shidneat,sheet="SAll")
