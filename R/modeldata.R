## data preparation and cascade modelling
rm(list=ls())

## flags for sensitivity analyses
shell <- FALSE # whether running from shell script or not
if(shell){
  ## running from shell
  args <- commandArgs(trailingOnly=TRUE)
  print(args)
  SA <- args[1]                  # none,base/lo/hi,cdr,txd
  if(SA == 'none'){
    SA <- ''
  } 
} else { #set by hand
  rm(list=ls()) #clear all
  shell <- FALSE #whether running from shell script or not
  
  ##sensitivity analyses (mostly for PT):
  ## '' = basecase assumes screening observed in TIPPI
  ## 'screen' = # assumes screening observed in TIPPI for SOC and all children presenting are screened for INT
  ## 'screenINT' = # assumes all children presenting are screened 
  sacases <- c('', 'screen', 'screenINT')
  SA <- sacases[1]
}

library(here)
library(glue)
library(ggplot2)
library(data.table)
library(ggthemes)
library(scales)
library(readxl)
library(discly)
library(HEdtree)
library(tidyverse)
rexel <- function(x,...) as.data.table(read_excel(x,...))
ssum <- function(x) sqrt(sum(x^2))
gh <- function(x) glue(here(x))

set.seed(123)
cns <- c("Cameroon","Kenya")
cnisos <- c('CMR','KEN')

## number of reps
nreps <- 1e4  # for generation of random samples

## country key
CK <- data.table(iso3=c('CMR', 'KEN'), country=c('Cameroon', 'Kenya')) #NOTE this is where the countries involved are coded
CK

countrykey <- data.table(iso3=cnisos,country=cns)
setkey(countrykey,iso3)
countrykey

save(countrykey,file=here('outdata/countrykey.Rdata'))


## ===== DISCOUNTED LIFE-YEARS TABLES
## NOTE discount rate set here
discount.rate <- 0.03
DRL <- c(0.03,0,0.05) #for sensitivity analyses
## make life-years

## calculate discounted life-years
## template:
LYT <- data.table(age=0:4,
                  age_group=c(rep('0-1',2),rep('2-4',3)),
                  LYS=0.0,LYS0=0.0)
for(DR in DRL){
  fn <- here('outdata/LYK.Rdata')
  if(DR!=0.03)
    fn <- gh('outdata/LYK_{round(1e2*DR)}.Rdata')
  ## make country/age key
  LYK <- list()
  for(iso in cnisos){
    ## iso <- cn
    tmp <- copy(LYT)
    tmp[,iso3:=iso]
    for(ag in tmp$age) tmp[age==ag,LYS:=discly(iso,ag,2020,
                                               dr=DR)]
    for(ag in tmp$age) tmp[age==ag,LYS0:=discly(iso,ag,2020,
                                                dr=0)]
    LYK[[iso]] <- tmp
  }
  LYK <- rbindlist(LYK)
  ## assume unweighted & collapse
  LYK <- LYK[,.(LYS=mean(LYS),LYS0=mean(LYS0)),
             by=.(iso3,age=age_group)]
  setkey(LYK,age)
  save(LYK,file=fn)
}

## ===== CASCADE DATA

# complete cascade 
# CASCP <- fread(here('indata/CASCP.csv')) # not disaggregated
CASCP <- fread(here('indata/CASCPAGE.csv')) # disaggregated by age
# CASCP <- fread(here('indata/CASCPOAGE.csv')) # disaggregated by age & now including treatment success

# CASCP[ ,arm:=ifelse(arm==0, 'SOC', 'INT')]
CASCP[country=='CMR' & age=='u2' & arm=='SOC']

CASCP <- CASCP[, c('metric', 'dx'):=NULL]
CASCP <- CASCP[,variable:=gsub('frac.', '', variable)]

tmp <- CASCP[variable=='presumptive']
tmp <- tmp[,variable:='screened1']
tmp <- tmp[,cascade:='Screened for TB symptoms']

tmp[country=='CMR']
tmp[ ,value:=ifelse(country=='CMR' & arm=='SOC', value*48.12, value)]
tmp[ ,value:=ifelse(country=='CMR' & arm=='INT', value*50.72, value)]
tmp[ ,value:=ifelse(country!='CMR' & arm=='SOC', value*66.79, value)]
tmp[ ,value:=ifelse(country!='CMR' & arm=='INT', value*73.15, value)]

tmp[ ,frac:=ifelse(country=='CMR' & arm=='SOC', frac*48.12, frac)]
tmp[ ,frac:=ifelse(country=='CMR' & arm=='INT', frac*50.72, frac)]
tmp[ ,frac:=ifelse(country!='CMR' & arm=='SOC', frac*66.79, frac)]
tmp[ ,frac:=ifelse(country!='CMR' & arm=='INT', frac*73.15, frac)]

tmp[country=='CMR']
tmp[country!='CMR']

CASCP <- rbind(CASCP, tmp)
unique(CASCP$cascade)
unique(CASCP$variable)


CASCP[,cascade:=factor(variable, 
                       levels = c('presented', 'screened1', 'screened', 'presumptive', 'evaluated', 'bac.assess', 'Xpert', 'diagnosed', 'att'),
                       labels = c('Initial care seeking', 'Screened for TB symptoms', 'TB symptoms','Presumptive TB', 'Evaluated for TB', 
                                  'Bacteriologically assessed', 'Tested on Xpert', 'TB diagnosis', 'TB treatment'),
                       ordered = TRUE)]

setkey(CASCP, country, age, arm, cascade, variable)
CASCP[country=='CMR' & age=='u2' & arm=='SOC']

ATT <- copy(CASCP)
ATT <- ATT[,arm:=ifelse(arm=='SOC', 'Standard of care', 'Intervention')]
ATT <- ATT[,arm:=factor(arm, levels = c('Standard of care', 'Intervention'))]
ATT <- ATT[cascade!='Screened for TB symptoms',]
ATT <- ATT[,age:=ifelse(age=='u2', '<2 years',  '2-5 years')]

## absolute
ggplot(ATT,aes(cascade,value,fill=arm)) +
  geom_bar(stat = 'identity',position='dodge') +
  facet_grid(country~age,scales='free')+
  scale_fill_colorblind() +
  xlab('') + scale_y_log10(label=comma)+ ylab('Number (log scale)')+
  theme(axis.text.x = element_text(angle = 45, vjust = 1.0, hjust=1)) +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 12))

ggsave(filename=here('plots/cascade_intervention_abs_age.pdf'),
       w=11,h=10)
ggsave(filename=here('plots/cascade_intervention_abs_age.png'),
       w=11,h=10)

## intervention cascade
## relative
ATR <- copy(CASCP)

setkey(ATR, country, age, arm, cascade, variable)
ATR[country=='CMR' & age=='u2' & arm!='SOC']
unique(ATR$cascade)

# ATR[cascade=='DS-TB treatment',tx:=value,by=country]
# ATR[,frac:=value/tx[9],by=country]
# ATR <- ATR[cascade!='DS-TB treatment',]
# ATR <- ATR[cascade!='Screened for TB symptoms',]
# ATR <- ATR[,aged:=ifelse(age=='u2', '<2 years',  '2-5 years')]

ggplot(ATR,aes(cascade,frac,color=arm,group=arm)) +
  geom_point() + geom_line()+
  facet_grid(country~age,scales='free')+
  xlab('') +
  scale_y_log10()+
  ylab('Number per child successfully treated for TB (log scale)')+
  theme_classic() + ggpubr::grids()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1.0, hjust=1)) +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 12))

ggsave(filename=here('plots/intervention_cascade_age.pdf'),
       w=6,h=5)
ggsave(filename=here('plots/intervention_cascade_age.png'),
       w=6,h=5)


# ATR <- dcast(ATR, country+cascade+variable~arm, value.var = 'frac')
# setorder(ATR, country, variable)
# 
# ATR <- ATR[,.(metric=variable, iso3=country, Baseline=SOC, Intervention=INT, ratio=INT/SOC)]
# ATR <- melt(ATR, id.vars = c(metric))
ATR <- ATR[arm=='INT',.(metric=cascade,iso3=country,age, frac)]

save(ATR,file=here('outdata/ATRage.Rdata')) #cascade per ATT initiation

load(file=here('outdata/ATRage.Rdata')) #cascade

# tmp <- CASCP

## merge against cascade/cost data
# unit.costs <- fread(gh('outdata/model_mean_total_costs.csv'))    #read old cost data
# unit.costs <- fread(gh('indata/model_mean_total_costs.csv'))    #read new cost data
unit.costs <- fread(gh('indata/model_mean_total_costs_age.csv'))    #read new cost data
# unit.costs <- unit.costs[facility=='DH',]
unit.costs <- unit.costs[, cascade:=factor(cascade, levels = c('screening', 'diagnosis', 'att'), labels = c('presented', 'presumptive', 'att'))]
unit.costs <- unit.costs[study_phase!=2,.(iso3, age=age_cat, arm=study_phase, metric=cascade, cost.m, cost.sd)]


tmp2 <- copy(unit.costs)
CD <- merge(tmp2[arm==1,.(iso3,age,metric, uc.soc=cost.m, uc.soc.sd=cost.sd)],
            tmp2[arm==3,.(iso3,age,metric, uc.int=cost.m, uc.int.sd=cost.sd)],
            by=c('iso3','age','metric'))

CD <- CD[,metric:=factor(metric, 
                         levels = c('presented', 'presumptive', 'att'),
                         labels = c('Initial care seeking', 'Presumptive TB', 'DS-TB treatment'))]

## merge into cascade
ATR <- merge(ATR,countrykey,by='iso3')
ATR[metric=='TB diagnosis',frac:=round(frac,0)]

# ART2 <- ATR[,.(iso3,age,metric,frac)]
ART2 <- rbind(ATR[,.(iso3,age,metric,frac)],
              ATR[metric=='TB diagnosis',.(iso3,age,metric='DS-TB treatment', frac)]) # TB treatment added to cascade

## adding TB symptoms & dx, which isn't in itself associated with resources
CD2 <- rbindlist(list(
  CD[,.(iso3,age, metric,uc.soc,uc.soc.sd,uc.int, uc.int.sd)],
  CD[metric=='DS-TB treatment',.(iso3,age,metric='TB diagnosis', uc.soc=0,uc.soc.sd=0,uc.int=0, uc.int.sd=0)]))

CD2 <- rbindlist(list(
  CD2[,.(iso3,age,metric,uc.soc,uc.soc.sd,uc.int, uc.int.sd)],
  CD2[metric=='Initial care seeking',.(iso3,age,metric='Screened for TB symptoms', uc.soc,uc.soc.sd,uc.int, uc.int.sd)]))

CD2 <- rbindlist(list(
  CD2[,.(iso3,age,metric,uc.soc,uc.soc.sd,uc.int, uc.int.sd)],
  CD2[metric=='Presumptive TB',.(iso3,age,metric='TB symptoms', uc.soc=0,uc.soc.sd=0,uc.int=0, uc.int.sd=0)]))

CD2 <- rbindlist(list(
  CD2[,.(iso3,age,metric,uc.soc,uc.soc.sd,uc.int, uc.int.sd)],
  CD2[metric=='Presumptive TB',.(iso3,age,metric='Evaluated for TB', uc.soc,uc.soc.sd,uc.int, uc.int.sd)]))

CD2 <- rbindlist(list(
  CD2[,.(iso3,age,metric,uc.soc,uc.soc.sd,uc.int, uc.int.sd)],
  CD2[metric=='Presumptive TB',.(iso3,age,metric='Bacteriologically assessed', uc.soc=0,uc.soc.sd=0,uc.int=0, uc.int.sd=0)]))

CD2 <- rbindlist(list(
  CD2[,.(iso3,age,metric,uc.soc,uc.soc.sd,uc.int, uc.int.sd)],
  CD2[metric=='Presumptive TB',.(iso3,age,metric='Tested on Xpert', uc.soc=0,uc.soc.sd=0,uc.int=0, uc.int.sd=0)]))

CD2[metric=='Initial care seeking',c('uc.soc','uc.soc.sd','uc.int', 'uc.int.sd'):=.(0,0,0, 0)]
CD2[metric=='Presumptive TB',c('uc.soc','uc.soc.sd','uc.int', 'uc.int.sd'):=.(0,0,0, 0)]

# merge into cascade
ART2 <- merge(ART2,CD2,by=c('iso3','age','metric'))
ART2[iso3!='CMR']

save(ART2,file=here('outdata/ART2age.Rdata')) #cascade + costs

load(file=here('outdata/ART2age.Rdata')) #cascade + costs

## ==================== baseline data ===========
B1 <- CASCP[arm=='SOC',]
B1

## --- compare cascade

## join
DB <- CASCP[,.(iso3=country,age,arm, metric=cascade,frac)]
DB <- DB[order(iso3,metric)]

vrs <- DB[,unique(metric)]
DB[,metric:=factor(metric,levels=vrs,ordered = TRUE)]
DB[,arm:=factor(arm,levels=c('SOC', 'INT'),ordered = TRUE)]
DB[,age:=factor(age,levels=c('u2', 'o2'),ordered = TRUE)]
DB[metric=='TB diagnosis',frac:=round(frac,0)]

DBR <- DB

ggplot(DBR,aes(metric,frac,col=iso3,group=paste(arm,iso3),
               lty=arm)) +
  geom_point() + geom_line() +
  facet_grid(~age) +
  scale_y_log10(label=comma)+ xlab('')+
  ylab('Number per children diagnosed for TB (log scale)')+
  theme_classic() + ggpubr::grids()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1.0, hjust=1))

## ggsave(filename=here('plots/cascade_compare.pdf'),w=7,h=8)
ggsave(filename=here('plots/cascade_compare_age.png'),w=7,h=8)

DBR <- DBR[metric %in% vrs[1:7]]

ggplot(DBR,aes(arm,frac,col=iso3,lty=metric,
               group=paste(iso3,metric))) +
  geom_point() + geom_line()+
  facet_grid(metric~age, scales = 'free') +
  ylab('Number per child diagnosed with TB') + xlab('Period')## +
## theme_classic() + ggpubr::grids()

## ggsave(filename=here('plots/cascade_compare2.pdf'),w=8,h=5)
ggsave(filename=here('plots/cascade_compare2.png'),w=8,h=5)

## comparison of increases in resources vs effect
## edatm <- edat[,.(RR=mean(RR)),by=.(quant,age,country)]

DBC <- dcast(DBR,iso3+age+metric ~ arm,value.var = 'frac')
DBC[,ratio:=INT/SOC]
## compute average for each metric to fill in gaps
# tmp <- DBC[,.(rat=mean(ratio,na.rm=TRUE)),by=metric]
# DBC <- merge(DBC,tmp,by='metric')
# DBC[is.na(Baseline),ratio:=rat] #replace missing with average
# DBC[,rat:=NULL]                 #remove additional data
# tmp <- DBC[metric=='Presumptive TB identified']
# tmp[,metric:="Screened for symptoms"]
# DBC <- rbind(DBC,tmp) #NOTE apply presumptive to screened also

## this gives the int/soc ratio of presumptive and testing per treatment
save(DBC,file=here('outdata/DBCage.Rdata')) #ratios int v bl
load(file=here('outdata/DBCage.Rdata')) #ratios int v bl


# redoing to capture screening sensitivity analyses
DBC <- merge(DBC,CK,by='iso3', all.x = TRUE)
ART2 <- merge(ART2[,.(iso3, age, metric, uc.soc, uc.soc.sd, uc.int, uc.int.sd)],
              DBC[,.(iso3, age, metric, frac=INT, ratio)],by=c('iso3','age','metric'),all.x=TRUE)
# ART2[,c('SOC','INT'):=NULL]

ART2[,.(iso3,age,metric,ratio)]
ART2 <- ART2[is.na(ratio), ratio:=1]
ART2 <- ART2[is.na(frac), frac:=1]

ART2[is.na(uc.soc.sd)|uc.soc.sd<0,uc.soc.sd:=uc.soc/40]
ART2[is.na(uc.int.sd),uc.int.sd:=uc.int/40]

## merge and compute both sets of costs
xtra <- ART2[,.(cost.soc=sum(uc.soc*frac/ratio),
                cost.soc.sd=ssum(frac*uc.soc.sd/ratio),
                cost.int=sum((uc.int)*frac),
                cost.int.sd=ssum(frac*uc.int.sd)
),
by=.(iso3, age)]

xtra[,frac.soc:=paste0(round(cost.soc),' (',
                       round(cost.soc.sd),')')]
xtra[,frac.int:=paste0(round(cost.int),' (',
                       round(cost.int.sd),')')]
xtra <- melt(xtra[,.(iso3,age,frac.soc,frac.int)],id=c('iso3','age'))
xtra[,metric:='cost']

cascadetab.soc <- rbind(ART2[,.(iso3,age,metric,
                                frac=paste0(
                                  round(frac/ratio,2)))],
                        xtra[variable=='frac.soc',
                             .(iso3,age,metric,frac=value)])
cascadetab.int <- rbind(ART2[,.(iso3,age,metric,
                                frac=paste0(round(frac,2)))],
                        xtra[variable=='frac.int',
                             .(iso3,age,metric,frac=value)])


cascadetab.soc <- dcast(cascadetab.soc,iso3+age~metric,
                        value.var = 'frac')
cascadetab.int <- dcast(cascadetab.int,iso3+age~metric,
                        value.var = 'frac')

## this block is for table 1 below
ctsm <- melt(cascadetab.soc,id=c('iso3','age'))
ctsm[,condition:='SOC']
ctim <- melt(cascadetab.int,id=c('iso3','age'))
ctim[,condition:='INT']
ctm <- rbind(ctsm,ctim)

unique(ctm$variable)
ccs <- c("iso3","age","Initial care seeking",
         "Screened for TB symptoms",
         "TB symptoms",
         "Presumptive TB",
         "Evaluated for TB",
         "Bacteriologically assessed",
         "Tested on Xpert",
         "TB diagnosis",
         "DS-TB treatment","cost")
setcolorder(cascadetab.soc,ccs)
setcolorder(cascadetab.int,ccs)
cascadetab.int[,iso3:=NULL]
cascadetab <- cbind(cascadetab.soc,cascadetab.int)
cascadetab

fn1 <- glue(here('outdata/cascadetab_age')) + SA + '.csv'
fwrite(cascadetab,file=fn1)

## ------ reformatting for upper part of Table 1 -----
ctm <- merge(ctm,CK,by='iso3')
ctm[,iso3:=NULL]
ctm <- ctm[order(country,age, condition,variable)]
ctmc <- dcast(ctm,age + condition + variable ~ country,value.var = 'value')
ctmc$condition <- factor(ctmc$condition,levels=c('SOC','INT'),ordered=TRUE)
ctmc$variable <- factor(ctmc$variable,levels=ccs[-1],ordered=TRUE)
setkey(ctmc,age,condition,variable)

levels(ctmc$variable)
levels(ART2$metric)
uc.soc <- ART2[,.(iso3, age, metric, uc.soc=paste0(uc.soc, ' (', uc.soc.sd, ')'))]
uc.int <- ART2[,.(iso3, age, metric, uc.int=paste0(uc.int, ' (', uc.int.sd, ')'))]

cmr_soc <- merge(ctmc[condition=='SOC',.(age, condition, variable, Cameroon)], 
                 uc.soc[iso3=='CMR',.(age, variable=metric, uc.soc)],
                 by=c('age','variable'))

cmr_soc <- cmr_soc[order(age,condition,variable, Cameroon, uc.soc)]
cmr_soc <- rbind(cmr_soc, ctmc[condition=='SOC' & variable=='cost',.(age,condition, variable, Cameroon=1, uc.soc=Cameroon)])

cmr_int <- merge(ctmc[condition=='INT',.(age, condition, variable, Cameroon)], 
                 uc.int[iso3=='CMR',.(age, variable=metric, uc.int)],
                 by=c('age','variable'))
cmr_int <- cmr_int[order(age, condition,variable, Cameroon, uc.int)]
cmr_int <- rbind(cmr_int, ctmc[condition=='INT' & variable=='cost',.(age, condition, variable, Cameroon=1, uc.int=Cameroon)])

ken_soc <- merge(ctmc[condition=='SOC',.(age, condition, variable, Kenya)], 
                 uc.soc[iso3=='KEN',.(age, variable=metric, uc.soc)],
                 by=c('age','variable'))
ken_soc <- ken_soc[order(age, condition,variable, Kenya, uc.soc)]
ken_soc <- rbind(ken_soc, ctmc[condition=='SOC' & variable=='cost',.(age, condition, variable, Kenya=1, uc.soc=Kenya)])

ken_int <- merge(ctmc[condition=='INT',.(age, condition, variable, Kenya)], 
                 uc.int[iso3=='KEN',.(age, variable=metric, uc.int)],
                 by=c('age','variable'))
ken_int <- ken_int[order(age, condition,variable, Kenya, uc.int)]
ken_int <- rbind(ken_int, ctmc[condition=='INT' & variable=='cost',.(age, condition, variable, Kenya=1, uc.int=Kenya)])

ctmc_soc_csts <- cbind(cmr_soc[,.(age, condition, variable, Cameroon, Cameroon.uc=uc.soc)],ken_soc[,.(Kenya, Kenya.uc=uc.soc)])
ctmc_int_csts <- cbind(cmr_int[,.(age, condition, variable, Cameroon, Cameroon.uc=uc.int)], ken_int[,.(Kenya, Kenya.uc=uc.int)])
ctmc_csts <- rbind(ctmc_soc_csts, ctmc_int_csts)
Table1ATT <- copy(ctmc_csts)

setkey(Table1ATT, age, condition, variable)
Table1ATT

fn1 <- glue(here('outdata/Table1ATTage')) + SA + '.csv'
fwrite(Table1ATT,file=fn1)

# save top cascade for CE table in modeloutcomes.R
cascadetop <- Table1ATT[condition=='SOC' & variable %in% c('TB symptoms', 'Evaluated for TB'), .(age,variable, Cameroon, Kenya)]
cascadetop <- cascadetop[,variable:=ifelse(variable=='TB symptoms', 'symptoms.soc', 'evaluated.soc')]
cascadetop <- cascadetop[,age:=ifelse(age=='u2', '0-1', '2-4')]
cascadetop <- melt(cascadetop, id.vars = c('age', 'variable'), variable.name = 'country')
cascadetop <- cascadetop[,value:=as.numeric(value)]
cascadetop <- dcast(cascadetop, country+age~variable)

fn1 <- glue(here('outdata/cascadetop')) + SA + '.Rdata'
save(cascadetop,file=fn1)

# Plotting
K <- merge(ART2[,.(iso3, age, metric, uc.soc, uc.soc.sd, uc.int, uc.int.sd)],
           DBC[,.(iso3, age, metric, frac=INT, ratio)],by=c('iso3','age','metric'),all.x=TRUE)

# tmp <- DBC1[metric=='Screened for TB symptoms',.(iso3, age, metric, frac=INT)]
# tmp <- cbind(tmp,
#              K[metric=='Screened for TB symptoms',.(uc.soc, uc.soc.sd, uc.int, uc.int.sd, ratio)])
# K <- rbind(K[metric!='Screened for TB symptoms'],
#              tmp)

K[is.na(ratio),ratio:=1] #TB tx or dx
K[is.na(frac),frac:=1]

K[iso3=='CMR'] #check
K[iso3!='CMR'] 

##NOTE this same structure (K) is developed in model data
## would be best to ensure consistency and remove from this file

## computing costs by activity: soc -> int cascade change + increment
KA <- copy(K)
KA[,cost.soc:=frac*uc.soc/ratio] #cost assuming resource use same as soc
KA[,cost.int:=frac*(uc.int)] #  scale-up already included
KA[iso3=='CMR']                     #inspect
KA[iso3!='CMR'] 

# if(SA=='screen'){                        #sensitivity analysis
#   KA <- KA[metric!='Initial care seeking',]
#   KA <- KA[metric=='Screened for TB symptoms',metric:='Initial care seeking']
# } 
# else {
#   KA <- KA[metric!='Screened for TB symptoms',]
# }


## formatting and plotting
KAM <- melt(KA[metric %in% c('Screened for TB symptoms',
                             'Evaluated for TB',
                             'DS-TB treatment'),
               .(metric,iso3,age,cost.soc,cost.int)],
            id.vars = c('metric','age','iso3'))
KAM <- merge(KAM,CK,by='iso3')
KAM$metric <- factor(KAM$metric,
                     levels=c('Screened for TB symptoms',
                              'Evaluated for TB',
                              'DS-TB treatment'),
                     labels=c('TB symptom screening',
                              'TB investigations',
                              'DS-TB treatment'),
                     ordered = TRUE)
KAM[grepl('soc',variable),variable:='Standard of care']
KAM[grepl('int',variable),variable:='Intervention']
KAM[,iso3:=ifelse(iso3=='CMR', 'Cameroon', 'Kenya')]
KAM[,age:=ifelse(age=='u2', '< 2 years', '2-5 years')]

GP <- ggplot(KAM,aes(age,value,fill=metric)) +
  geom_bar(stat='identity',position = position_stack(reverse = TRUE)) +
  facet_grid(country~variable) +
  scale_fill_colorblind() +
  scale_y_continuous(label=comma)+
  xlab('') + ylab('Cost per child treated (USD)')  +
  theme(axis.text.x = element_text(angle = 0, vjust = 1.0, hjust=1),
        legend.position = 'top',legend.title = element_blank())
if(!shell) GP

fn1 <- glue(here('plots/cost_cascade_age')) + SA + '.png'
ggsave(GP,file=fn1,w=6,h=5); ## ggsave(GP,file=fn2,w=10,h=10)

GP <- ggplot(KAM,
             aes(age,value,fill=metric)) +
  geom_bar(stat='identity',position = position_fill(reverse = TRUE)) +
  facet_grid(country~variable) +
  scale_fill_colorblind() +
  scale_y_continuous(label=percent)+
  xlab('') + ylab('Fraction of cost per child treated by stage') +
  theme(axis.text.x = element_text(angle = 0, vjust = 1.0, hjust=1),
        legend.position = 'top',legend.title = element_blank())
if(!shell) GP

fn1 <- glue(here('plots/cost_cascade_age2')) + SA + '.png'
ggsave(GP,file=fn1,w=6,h=5); ## ggsave(GP,file=fn2,w=10,h=10)

## ASM age splits
age_split <- fread(here('indata/cascade_summary_age_arm.csv'))
age_split <- melt(age_split, id.vars = c('Arm', 'cascade', 'age'), variable.name = 'iso3')
age_split <- age_split[,.(n=sum(value)),
                       by=.(iso3,arm=Arm, age, cascade)]

age_split <- age_split[, age:=ifelse(age=='u2', '0-1', '2-4')]
age_split <- age_split[, arm:=ifelse(arm==1, 'SOC', 'INT')]
age_split <- dcast(age_split, iso3+arm+cascade~age, value.var = 'n')

present <- age_split[cascade=='presented',.(iso3, arm, `0-1`,`2-4`)]
symptoms <- age_split[cascade=='presumptive TB (screened)',.(iso3, arm, `0-1`,`2-4`)]
presumptive <- age_split[cascade=='presumptive TB (enrolled)',.(iso3, arm, `0-1`,`2-4`)]
diagnosis <- age_split[cascade=='TB diagnosed',.(iso3, arm, `0-1`,`2-4`)]
# treatment <- age_split[cascade=='presented',]

present[,qty:='present']
symptoms[,qty:='symptoms']
presumptive[,qty:='presumptive']
diagnosis[,qty:='diagnosis']
# treatment[,qty:='diagnosis']

ASM <- rbindlist(list(present,symptoms,presumptive,diagnosis))

# ASM <- ASM[qty=='present']
# ASM <- ASM[,arm:=ifelse(arm=='SOC', 'Baseline', 'Intervention')]

save(ASM,file=here('outdata/ASM.Rdata'))

# HIV prevalence from INPUT study data
load(here('outdata/cohortdata.Rdata'))
popn_hiv <- setDT(popn_hiv)

popn_hiv <- popn_hiv[arm!=2,.(country, arm=as.character(arm), age=age_cat,hiv_status, n)]
# popn_hiv <- popn_hiv |> ungroup()
#   group_by(country, arm, hiv_status) |> 
#   summarise(n())

popn_hiv <- popn_hiv[, n:=sum(n),
                     by=.(country, arm, hiv_status)]

popn_hiv <- popn_hiv[hiv_status==0,value:=n]
popn_hiv[,total:=n+value[1],by=.(country, arm, age)]
popn_hiv <- popn_hiv[hiv_status==1, .(iso3=country, arm, age, hiva=n, total)]
popn_hiv <- popn_hiv[, c('hiva', 'total'):=.(sum(hiva), sum(total)),
                     by=.(iso3, arm, age)]
popn_hiv <- unique(popn_hiv)

BL <- popn_hiv[arm==1,]
BL <- BL[rep(1:nrow(BL),each=nreps)]
BL[,id:=rep(1:nreps,nrow(CK)*2)] 
# BL <- rbind(BL,BL)
# BL[,age:=rep(c('0-1','2-4'),each=nrow(BL)/2)]
BL[,age:=ifelse(age=='u2', '0-1', '2-4')]
BL[,hiv:=rbeta(nrow(BL),hiva,total)]
BL <- BL[,.(id,iso3,age,hiv)]

## HIV prevalence from INT data
INT <- popn_hiv[arm==3,]
INT <- INT[rep(1:nrow(INT),each=nreps)]
# INT[,id:=rep(1:max(PSA$id),nrow(CK))] 
# INT <- rbind(INT,INT)
# INT[,age:=rep(c('0-1','2-4'),each=nrow(INT)/2)]
INT[,id:=rep(1:nreps,nrow(CK)*2)] 
INT[,age:=ifelse(age=='u2', '0-1', '2-4')]
INT[,hiv:=rbeta(nrow(INT),hiva,total)]
INT[,hiv:=rbeta(nrow(INT),hiva,total)]
INT <- INT[,.(id,iso3,age,hiv.int=hiv)]

HIVtrial <- merge(BL,INT,by=c('id','iso3', 'age'))

HIVtrial <- merge(HIVtrial,CK,by='iso3')

HIVtrial[,.(hiv.soc=mean(hiv),hiv.int=mean(hiv.int)), by=country]

fn <- here('outdata/HIVtrial.Rdata')
save(HIVtrial,file=fn)

# save for Appendix
x <- popn_hiv[,.(hiva=sum(hiva), total=sum(total)), by = .(iso3)]
x[,prop:=round(hiva/total*100,1)]
hiv <- popn_hiv[,.(iso3, arm, age, hiva=round(hiva/total*100,1))]
hiv <- hiv[, arm:=ifelse(arm==1, 'Standard of care', 'Intervention')]
hiv <- hiv[, Country:=ifelse(iso3=='CMR', 'Cameroon', 'Kenya')]
hiv[,age:=ifelse(age=='u2', '0-1', '2-4')]
hiv <- dcast(hiv, Country+age~arm, value.var = 'hiva')
hiv
fn <- here('outdata/hiv_prevalence.csv')
fwrite(hiv,file=fn)

## RRs
# country-specific
effects <- fread(here('indata/int_effects_details.csv'), header = TRUE)    #read effect data

tmp <- effects |>
  extract(`MEDIAN (IQR)`, into = c("mid", "lo"), "([^(]+)\\s*[^0-9]+([0-9].*).") |>
  separate(lo,c("lo","hi"),"-") |>
  mutate_at(c("mid", "lo","hi"), as.numeric)
tmp <- setDT(tmp)
tmp <- tmp[c(1:9,40:42)]
tmp1 <- getLNparms(tmp[,mid],(tmp[,hi]-tmp[,lo])^2/3.92^2,med=FALSE) #NOTE changed to mean vs median
curve(dlnorm(x,tmp1$mu[3],tmp1$sig[3]),from=0,to=10)
curve(dlnorm(x,tmp1$mu[6],tmp1$sig[6]),from=0,to=10)
curve(dlnorm(x,tmp1$mu[9],tmp1$sig[9]),from=0,to=10)

curve(dlnorm(x,tmp1$mu[10],tmp1$sig[10]),from=0,to=10)
curve(dlnorm(x,tmp1$mu[11],tmp1$sig[11]),from=0,to=10)
curve(dlnorm(x,tmp1$mu[12],tmp1$sig[12]),from=0,to=10)
# 
# curve(dgamma(x,tmp1$mu[3],scale=tmp1$sig[3]),from=0,to=100)
# curve(dgamma(x,tmp1$mu[6],scale=tmp1$sig[6]),from=0,to=5)
# 
# curve(dgamma(x,tmp1$mu[8],scale=tmp1$sig[8]),from=0,to=100)

curve(dgamma(x,tmp[,mid][3]^2/((tmp[,hi][3]-tmp[,lo][3])^2/3.92^2),scale=((tmp[,hi][3]-tmp[,lo][3])^2/3.92^2)/tmp[,mid][3]),from=0,to=10)
curve(dgamma(x,tmp[,mid][6]^2/((tmp[,hi][6]-tmp[,lo][6])^2/3.92^2),scale=((tmp[,hi][6]-tmp[,lo][6])^2/3.92^2)/tmp[,mid][6]),from=0,to=20)
curve(dgamma(x,tmp[,mid][9]^2/((tmp[,hi][9]-tmp[,lo][9])^2/3.92^2),scale=((tmp[,hi][9]-tmp[,lo][9])^2/3.92^2)/tmp[,mid][9]),from=0,to=10)
curve(dgamma(x,tmp[,mid][10]^2/((tmp[,hi][10]-tmp[,lo][10])^2/3.92^2),scale=((tmp[,hi][10]-tmp[,lo][10])^2/3.92^2)/tmp[,mid][10]),from=0,to=10)
curve(dgamma(x,tmp[,mid][11]^2/((tmp[,hi][11]-tmp[,lo][11])^2/3.92^2),scale=((tmp[,hi][11]-tmp[,lo][11])^2/3.92^2)/tmp[,mid][11]),from=0,to=10)
curve(dgamma(x,tmp[,mid][12]^2/((tmp[,hi][12]-tmp[,lo][12])^2/3.92^2),scale=((tmp[,hi][12]-tmp[,lo][12])^2/3.92^2)/tmp[,mid][12]),from=0,to=10)
## curve(dlnorm(x,tmp1$mu[6],tmp1$sig[6]),from=0,to=100)
tmp[,DISTRIBUTION:=paste0("LN(",tmp1$mu,",",tmp1$sig,")")] #LN distributions

# dig_effects <- tmp[1:6]
# dig_effects[,DISTRIBUTION:=paste0("G(",dig_effects[,mid]^2/((dig_effects[,hi]-dig_effects[,lo])^2/3.92^2),",",((dig_effects[,hi]-dig_effects[,lo])^2/3.92^2)/dig_effects[,mid],")")] # Beta distributions
# 
# txsx <- tmp[10:12]
# txsx1 <- getAB(txsx[,mid],(txsx[,hi]-txsx[,lo])^2/3.92^2) #NOTE changed to mean vs median
# 
# curve(dbeta(x,txsx1$a[1],txsx1$b[1]))
# txsx[,DISTRIBUTION:=paste0("LN(",txsx1$mu,",",txsx1$sig,")")] #LN distributions
# 
# tmp <- rbind(dig_effects,txsx)
PE <- parse.parmtable(data.frame(tmp[,.(NAME, DISTRIBUTION)]))
PSAE <- makePSA(nreps,PE) # PSAE = PSA effects
PSAE[,id:=1:nrow(PSAE)]
names(PSAE)
names(PSAE)[grepl('tb.symptomsRR', names(PSAE))] <- gsub('tb.symptomsRR', 'RRpresTB', names(PSAE)[grepl('tb.symptomsRR', names(PSAE))])
names(PSAE)[grepl('tb.investigationRR', names(PSAE))] <- gsub('tb.investigationRR', 'RRevalTB', names(PSAE)[grepl('tb.investigationRR', names(PSAE))])
names(PSAE)[grepl('tb.diagnosis', names(PSAE))] <- gsub('tb.diagnosis', '', names(PSAE)[grepl('tb.diagnosis', names(PSAE))])
names(PSAE)
PSAEL <- melt(PSAE, id.vars = 'id')  # PSAE = PSA effects long
PSAEL <- PSAEL[,parm:=gsub('ken.|cmr.', '', variable)]
PSAEL[,country:='Both']; PSAEL[grepl('cmr.',variable),country:='Cameroon']; PSAEL[grepl('ken',variable),country:='Kenya']
PSAEL[,variable:=gsub("ken.|cmr.","",variable)]
PSAEL <- dcast(PSAEL,id+country~variable,value.var = 'value')
PSAEL

Parms <- copy(tmp)
Parms[,DISTRIBUTION:=tmp$DISTRIBUTION]
Parms[DISTRIBUTION=='LN(NA,NA)', DISTRIBUTION:=NA]
Parms[,DISTRIBUTION:=paste0("LN(",round(tmp1$mu,5),",",round(tmp1$sig,5),")")] #LN distributions

Parms <- Parms |> 
  mutate(`MEDIAN (IQR)`=paste0(round(mid,5)," (",round(lo,5),"-",round(hi,5),")")) |> 
  select(NAME, DISTRIBUTION, `MEDIAN (IQR)`, DESCRIPTION, SOURCE)
fn1 <- glue(here('outdata/TableS3')) + SA + '.csv'
fwrite(Parms,file=fn1)

summary(PSAEL[country=='Cameroon',.(RR)])
summary(PSAEL[country=='Kenya',.(RR)])

INTE <- copy(PSAEL)

# Pooled intervention effects
tmp <- effects |>
  extract(`MEDIAN (IQR)`, into = c("mid", "lo"), "([^(]+)\\s*[^0-9]+([0-9].*).") |>
  separate(lo,c("lo","hi"),"-") |>
  mutate_at(c("mid", "lo","hi"), as.numeric)
tmp <- setDT(tmp)
tmp <- tmp[c(1:3,41),]
tmp1 <- getLNparms(tmp[,mid],(tmp[,hi]-tmp[,lo])^2/3.92^2,med=FALSE) #NOTE changed to mean vs median
## curve(dlnorm(x,tmp1$mu[3],tmp1$sig[3]),from=0,to=3)
# tmp[,DISTRIBUTION:=paste0("LN(",tmp1$mu,",",tmp1$sig,")")] #LN distributions
tmp[,DISTRIBUTION:=paste0("G(",tmp[,mid]^2/((tmp[,hi]-tmp[,lo])^2/3.92^2),",",((tmp[,hi]-tmp[,lo])^2/3.92^2)/tmp[,mid],")")] # Beta distributions

# dig_effects <- tmp[1:3]
# dig_effects[,DISTRIBUTION:=paste0("G(",dig_effects[,mid]^2/((dig_effects[,hi]-dig_effects[,lo])^2/3.92^2),",",((dig_effects[,hi]-dig_effects[,lo])^2/3.92^2)/dig_effects[,mid],")")] # Beta distributions
# 
# txsx <- tmp[4]
# txsx1 <- getLNparms(txsx[,mid],(txsx[,hi]-txsx[,lo])^2/3.92^2,med=FALSE) #NOTE changed to mean vs median
# txsx[,DISTRIBUTION:=paste0("LN(",txsx1$mu,",",txsx1$sig,")")] #LN distributions
# tmp <- rbind(dig_effects,txsx)

PE <- parse.parmtable(data.frame(tmp[,.(NAME, DISTRIBUTION)]))
PSAE <- makePSA(nreps,PE) # PSAE = PSA effects
PSAE[,id:=1:nrow(PSAE)]
names(PSAE)
names(PSAE)[grepl('tb.symptomsRR', names(PSAE))] <- gsub('tb.symptomsRR', 'RRpresTB', names(PSAE)[grepl('tb.symptomsRR', names(PSAE))])
names(PSAE)[grepl('tb.investigationRR', names(PSAE))] <- gsub('tb.investigationRR', 'RRevalTB', names(PSAE)[grepl('tb.investigationRR', names(PSAE))])
names(PSAE)[grepl('tb.diagnosis', names(PSAE))] <- gsub('tb.diagnosis', '', names(PSAE)[grepl('tb.diagnosis', names(PSAE))])
names(PSAE)
PSAEL <- melt(PSAE, id.vars = 'id')  # PSAE = PSA effects long
PSAEL <- PSAEL[,parm:=gsub('ken.|cmr.', '', variable)]
PSAEL[,country:='Both']; PSAEL[grepl('cmr.',variable),country:='Cameroon']; PSAEL[grepl('ken',variable),country:='Kenya']
PSAEL[,variable:=gsub("ken.|cmr.","",variable)]
PSAEL <- dcast(PSAEL,id+country~variable,value.var = 'value')
PSAEL

ParmsP <- copy(tmp)
ParmsP[,DISTRIBUTION:=tmp$DISTRIBUTION]
ParmsP[DISTRIBUTION=='LN(NA,NA)', DISTRIBUTION:=NA]
ParmsP[,DISTRIBUTION:=paste0("LN(",round(tmp1$mu,5),",",round(tmp1$sig,5),")")] #LN distributions

# ParmsTab3 <- rbind(ParmsP, Parms)     # save for later use in creating parameter tables for Appendix
# fn1 <- glue(here('outdata/TableS3')) + SA + '.csv'
# fwrite(ParmsTab3,file=fn1)

INTEP <- copy(PSAEL)
fn <- here('outdata/INTE.Rdata')
save(INTE, INTEP,file=fn)

## trial treatment outcomes
trt_outcomes <- data.frame(
  iso3 = rep(c('CMR', 'KEN'), each=2),  
  arm =  rep(c('SOC', 'INT'), 2),    
  cases = c(2, 35, 77, 41),
  trt_completed = c(2, 24, 57, 40),
  cured = c(0, 4, 6, 0),
  ltfu = c(0, 2, 4, 1),
  not_evaluated = c(0, 1, 4, 0),
  who_favourable = c(2, 28, 63, 40),
  deaths = c(0, 2, 6, 0))
trt_outcomes <- setDT(trt_outcomes)

## trial=based CFR
TCFR <- trt_outcomes[,.(iso3, arm, deaths, cases=cases-not_evaluated)]

# WHO CFR
who_cfr <- data.table(
  iso3 = c('CMR', 'KEN'),
  deaths = c(7400+3100, 17000+8900),
  cases = c(44000, 128000)
)


pooled <- who_cfr[,.(alpha_who=sum(deaths)-1, beta_who=sum(cases)-sum(deaths)-1), by = .(iso3)]

curve(dbeta(x,10499,33499),n=1e3)
fac <- 10000 #how much weaker (same mean)? I want a>1 so there is a peak within interval
curve(dbeta(x,10499/fac,33499/fac),n=1e3)
pooled <- pooled[,.(iso3, alpha_who=alpha_who/fac, beta_who=beta_who/fac)]

TCFR <- merge(TCFR,
              pooled,
              by = 'iso3')

# TCFR <- trt_outcomes[,.(iso3, arm, deaths, cases)]
LTFU <- trt_outcomes[,.(iso3, arm, ltfu, ltfu2=ltfu+not_evaluated)]

TCFR <- merge(TCFR,
              LTFU,
              by = c('iso3', 'arm'))

# set priors using pooled SOC data from the two countries
pooled <- rbind(TCFR[arm=='SOC',.(alpha=sum(deaths)-1, beta=sum(cases)-sum(deaths)-1), by = .(arm)],
                TCFR[arm=='SOC',.(arm='INT',alpha=sum(deaths)-1, beta=sum(cases)-sum(deaths)-1)])

fac <- 4 # how much weaker (same mean)? I want a>1 so there is a peak within interval
pooled <- pooled[,.(arm, alpha=alpha/fac, beta=beta/fac)]

TCFR <- merge(TCFR,
              pooled,
              by = 'arm')


TCFR <- TCFR[rep(1:nrow(TCFR),each=nreps)]
TCFR <- TCFR[,cfr:=rbeta(nrow(TCFR),alpha+deaths, beta+cases-deaths)] # approx flat conjugate prior
TCFR[,who_cfr:=rbeta(nrow(TCFR),alpha_who+deaths, beta_who+cases-deaths)] # who prior
TCFR[,id:=rep(1:nreps,nrow(TCFR)/nreps)]


TCFR[,.(cfr=mean(cfr), ltfu = mean(ltfu), who_cfr=mean(who_cfr)), by = .(iso3, arm)]

TCFR[,.(cfr=mean(cfr), ltfu = mean(ltfu), cfrlftu=mean(cfr+(cfr*ltfu/cases)), who_cfr=mean(who_cfr)), by = .(iso3, arm)]

# extract prior and mean CFRs for each arm
CFR <- TCFR[,.(DISTRIBUTION=paste0("B(",round(mean(alpha),5),",",round(mean(beta),5),")"),
               cfr=round(median(cfr),5), cfr_lo=round(quantile(cfr, 0.25),5), cfr_hi=round(quantile(cfr, 0.95),5)), by = .(iso3, arm)]

CFR <- CFR |> 
  mutate(`MEDIAN (IQR)`=paste0(round(cfr,5)," (",round(cfr_lo,5),"-",round(cfr_hi,5),")")) |> 
  select(-cfr, -cfr_lo, -cfr_hi)
fn1 <- glue(here('outdata/CFR')) + '.csv'
fwrite(CFR,file=fn1)

# 
# TCFR <- dcast(TCFR[,.(iso3,id,arm,cfr,who_cfr, ltfu, ltfu2)],
#               iso3+id+ltfu~arm,value.var = c('cfr', 'who_cfr'), drop = FALSE)
soc <- TCFR[arm=='SOC',.(iso3,id,arm,cfr_soc=cfr,who_cfr_soc=who_cfr, cfr_ltfu_soc=cfr+(cfr*ltfu/cases), ltfu_soc=ltfu/cases, ltfu2_soc=ltfu2/cases)]
int <- TCFR[arm!='SOC',.(iso3,id,arm,cfr_int=cfr,who_cfr_int=who_cfr, cfr_ltfu_int=cfr+(cfr*ltfu/cases), ltfu_int=ltfu/cases, ltfu2_int=ltfu2/cases)]

TCFR <- cbind(soc[,.(iso3,id,cfr_soc,who_cfr_soc, cfr_ltfu_soc, ltfu_soc, ltfu2_soc)],
              int[,.(cfr_int,who_cfr_int, cfr_ltfu_int, ltfu_int, ltfu2_int)])


TCFR[,.(cfr_soc=mean(cfr_soc), cfr_ltfu_soc=mean(cfr_ltfu_soc), cfr_int=mean(cfr_int), cfr_ltfu_int=mean(cfr_ltfu_int)), by = .(iso3)]

## trial=based treatment success
TSXS <- trt_outcomes[,.(iso3, arm, who_favourable, cases)]
LTFU <- trt_outcomes[,.(iso3, arm, ltfu, ltfu2=ltfu+not_evaluated)]

TSXS <- merge(TSXS,
              LTFU,
              by = c('iso3', 'arm'))

TSXS[,.(sum(who_favourable)/sum(cases)), by = .(arm, iso3)]

# set priors using pooled SOC data from the two countries
pooled <- rbind(TSXS[arm=='SOC',.(alpha=sum(who_favourable)-1, beta=sum(cases)-sum(who_favourable)-1), by = .(arm)],
                TSXS[arm=='SOC',.(arm='INT',alpha=sum(who_favourable)-1, beta=sum(cases)-sum(who_favourable)-1)])

curve(dbeta(x,64, 13),n=1e3)
fac <- 4 # how much weaker (same mean)? I want a>1 so there is a peak within interval
curve(dbeta(x,64/fac,13/fac),n=1e3)
pooled <- pooled[,.(arm, alpha=alpha/fac, beta=beta/fac)]

TSXS <- merge(TSXS,
              pooled,
              by = 'arm')
TSXS <- TSXS[rep(1:nrow(TSXS),each=nreps)]
TSXS <- TSXS[,TSXS:=rbeta(nrow(TSXS),alpha+who_favourable, beta+cases-who_favourable)] # approx flat conjugate prior
TSXS[,id:=rep(1:nreps,times=4)]


treatment_success <- TSXS[,.(DISTRIBUTION=paste0("B(",round(mean(alpha),5),",",round(mean(beta),5),")"),
                             tsxs=round(mean(TSXS),5), tsxs_lo=round(quantile(TSXS, 0.25),5), tsxs_hi=round(quantile(TSXS, 0.95),5)), by = .(iso3, arm)]

treatment_success <- treatment_success |> 
  mutate(`MEDIAN (IQR)`=paste0(round(tsxs,5)," (",round(tsxs_lo,5),"-",round(tsxs_hi,5),")")) |> 
  select(-tsxs, -tsxs_lo, -tsxs_hi)

fn1 <- glue(here('outdata/treatment_success')) + '.csv'
fwrite(treatment_success,file=fn1)

soc <- TSXS[arm=='SOC',.(iso3,id,arm,SOC=TSXS,SOC_LTFU=TSXS + TSXS*ltfu/cases)]
int <- TSXS[arm=='INT',.(iso3,id,arm,INT=TSXS,INT_LTFU=TSXS + TSXS*ltfu/cases)]

TSXS <-cbind(soc[,.(iso3,id,arm,SOC,SOC_LTFU)], 
             int[,.(INT,INT_LTFU)])

TSXS[,.(TSXS=mean(SOC), TSXSL=mean(SOC_LTFU), TSXSI=mean(INT), TSXSIL=mean(INT_LTFU)), by = .(iso3)]
TSXS[,.(mean(INT/SOC))]
TSXS[,.(mean(INT/SOC)), by = .(iso3)]

fn <- here('outdata/TCFR.Rdata')
save(TCFR, TSXS, file=fn)

## NOTE below here is WHO data and doesn't relate to costing
## === making CDRs and getting some relevant HIV data ===
whodir <- glue('~/Dropbox/WHO_TBreports/data2020/')

## === load WHO notifications
## http://www.who.int/tb/country/data/download/en/
fn <- here('indata/TB_notifications_2020-10-15.csv')
N <- fread(fn)

nmz <- paste0('newrel_',c('m04','f04',
                          'm59','f59',
                          'm1014','f1014',
                          'm014','f014'))
nmz <- c('iso3','year',nmz)

## reduce to relevant data: 2018 has no NA
NP <- N[year==2018,..nmz]
NP <- melt(NP,id.vars = c('iso3','year'))
NP[,sex:=ifelse(grepl("f",variable),'F','M')]
NP[,age:=gsub("[a-z]|_","",variable)]
NP <- NP[iso3 %in% cnisos]
NP
NP <- NP[age %in% c('04','014')]
NP

## === load WHO age-specific incidence estimates
fn <- here('indata/TB_burden_age_sex_2020-10-15.csv')
A <- fread(fn)
## keep only relevant categories
A <- A[year==2019]
A <- A[sex!='a']
A <- A[age_group %in% c('0-4','0-14','5-14')]
A <- A[risk_factor=='all']
A[,age:=gsub('-','',age_group)]
## harmonize namings
A[sex=='f',sex:='F']
A[sex=='m',sex:='M']
unique(A[,.(sex,age)])                  #check
A[,best.sd:=(hi-lo)/3.92]
A <- A[iso3 %in% cnisos]
A

## HIV
fn <- here('indata/TB_burden_countries_2020-10-15.csv')
H <- fread(fn)
H <- H[year==2019,.(iso3,e_tbhiv_prct,e_tbhiv_prct_lo,e_tbhiv_prct_hi)]
H <- H[iso3 %in% cnisos]
H[,hiv:=e_tbhiv_prct/100]
H[,hiv.sd:=(e_tbhiv_prct_hi-e_tbhiv_prct_lo)/392] #adults of course
H <- H[,.(iso3,hiv,hiv.sd)]

hfn <- here('outdata/H.Rdata')
save(H,file=hfn)

## === merge data
AN <- merge(NP[,.(iso3,sex,age,notes=value)],
            A[,.(iso3,sex,age,inc=best,lo,hi)],
            by=c('iso3','sex','age'),all.x=TRUE,all.y=FALSE)

ANO <- AN[age=='014']
ANY <- AN[age=='04']
ANB <- merge(ANY[,.(iso3,sex,notes.04=notes,inc.04=inc,
                    inc.sd.04=(hi-lo)/3.92)],
             ANO[,.(iso3,sex,notes.514=notes,inc.514=inc,
                    inc.sd.514=(hi-lo)/3.92)],
             by=c('iso3','sex')
)
CDRu <- ANB[,.(rel.sd=mean(inc.sd.514/inc.514)),by=iso3]
ANB[,c('notes.514','inc.514'):=.(notes.514-notes.04,inc.514-notes.04)]
CDR <- ANB[,.(notes.04=sum(notes.04),inc.04=sum(inc.04),
              notes.514=sum(notes.514),inc.514=sum(inc.514)),by=iso3]
CDR[,c('cdr04','cdr514'):=.(notes.04/inc.04,notes.514/inc.514)]
CDR[,totnotes:=notes.04+notes.514]
CDR[,c('frac04','frac514'):=.(notes.04/totnotes,notes.514/totnotes)]
CDR <- melt(CDR[,.(iso3,frac04,frac514,cdr04,cdr514)],id='iso3')
CDR[,qty:='cdr']
CDR[grepl('frac',variable),qty:='frac']
CDR[,age:='0-4']
CDR[grepl('14',variable),age:='5-14']
CDR <- CDR[,.(iso3,qty,age,value)]
CDR <- merge(CDR,CDRu,by='iso3',all.x=TRUE)
CDR[,cdr.v:=NA_real_]
CDR[,cdr.v:=(rel.sd*value)^2]
CDR[,rel.sd:=NULL]
CDR <- CDR[iso3 %in% cnisos]

cdrfn <- here('outdata/CDR.Rdata')
save(CDR,file=cdrfn)

## Figure 1

## ===== SIMPLE CONVERSIONS FROM CSV TO RDATA FOR CONSISTENCY
## extract HIV
# already available as cohortdata.Rdata
# 
# ## INT cascade data
# INPUT <- fread(here('indata/CASCP.csv'))
# save(INPUT,file=here('outdata/INPUT.Rdata'))

## CE thresholds
CET <- fread(here('indata/CEthresholds.csv'))
CET[,`1x GDP`:=as.numeric(gsub(",","",`1x GDP`))] ## format threshold data
CET[,`3x GDP`:=as.numeric(gsub(",","",`3x GDP`))]
CETM <- CET[,.(iso3=Country,`1x GDP`,`3x GDP`,Y1a,Y1b,Y2a,Y2b)]
CETM <- melt(CETM,id='iso3')
CETM[,value:=as.numeric(value)]
CETM <- merge(CETM,countrykey,by='iso3')
names(CETM)[2] <- 'threshold'
tmp <- CETM[threshold=='1x GDP'] #add in 1/2 GDP as in Chi Vassall
tmp[,value:=value/2]
tmp[,threshold:='0.5x GDP']
CETM <- rbind(CETM,tmp)
save(CETM,file=here('outdata/CETM.Rdata'))

## modelling parmeters
PD <- read.csv(here('indata/INPUTparms.csv'))
save(PD,file=here('outdata/PD.Rdata'))

