# rm(list=ls())

set.seed(1234)

## flags for sensitivity analyses
shell <- FALSE #whether running from shell script or not
if(shell){
  ## running from shell
  args <- commandArgs(trailingOnly=TRUE)
  print(args)
  SA <- args[1]                  #none,base/lo/hi,cdr,txd
  if(SA == 'none'){
    SA <- ''
  } 
} else { #set by hand
  rm(list=ls()) #clear all
  shell <- FALSE #whether running from shell script or not
  ##sensitivity analyses (mostly for PT):
  ## '' = basecase
  ## 'discr'='base'/'lo'/'hi'
  ## 'screen' = # assumes screening observed in TIPPI for SOC and all children presenting are screened for INT
  ## 'screenINT' = # assumes all children presenting are screened 
  ## 'txd' = making the completion influence tx outcome
  ## 'pooled_eff' = pooled intervention effect (not stratified by country)
  sacases <- c('','lo','hi', 'screen', 'screenINT', 'cfrltfu', 'whocfr', 'systRev', 'txd', 'pooled_eff')
  SA <- sacases[1]
}

## libraries
library(here)
library(tidyverse)
library(data.table)
library(HEdtree)
library(lattice)
library(glue)

## for CEAC plotting
source(here('R/inputfunctions.R')) #CEAC & plotting utils
gh <- function(x) glue(here(x))

## ===== INPUT DATA

load(file=here('outdata/LYK.Rdata'))  # LYs discounted
load(file=here('outdata/ART2age.Rdata')) # ATT cascade & costs NOTE 2 update w/costs
# load(file=here('outdata/CDR.Rdata')) # CDR
load(file=here('outdata/DBCage.Rdata')) # cascade ratios for int v soc
load(file=here('outdata/ASM.Rdata')) # age splits from pre/post data
load(file=here('outdata/H.Rdata')) # HIV by country from baseline data
# load(file=here('outdata/HIVtrial.Rdata')) # HIV by country from study data
load(file=here('outdata/cascadetop.Rdata')) # top cascade data
load(file=here('outdata/INTE.Rdata')) # trial-based INT effects
load(file=here('outdata/TCFR.Rdata')) # trial-based CFR
load(file=here('outdata/CETM.Rdata')) # CE thresholds
CK <- data.table(iso3=c('CMR', 'KEN'), country=c('Cameroon', 'Kenya'))

PD <- read.csv(here('indata/INPUTparms.csv'))
PZ <- parse.parmtable(PD)              #make into parm object

set.seed(1234)

## number of reps
nreps <- 1e4

#top to plot in CEAC curves
ceactop <- 3e3 

## data brought in on:
## tx success - txsuccess
## parameters - CFRdatam
## trial effects - PAEL
## unit costs & cascade for INT - ART2
## ratio of cascades INT v SOC - DCM
## HIV data - in this file from H
## life expectancy - LYK

## ================= PREPARE DATA SETS ===================

## ====== make PSA df for outcome data
PSA <- makePSA(nreps,PZ,
               dbls = list(c('hivartOR:mn','hivartOR:sg')))
PSA[,id:=1:nrow(PSA)]
## make ontx HIV+ outcomes
PSA[,ontxHAY:=ilogit(logit(ontxY) + `hivartOR:mn` + `hivartOR:sg`)] #on ART
PSA[,ontxHAO:=ilogit(logit(ontxO) + `hivartOR:mn` + `hivartOR:sg`)] #on ART
PSA[,ontxHAY:=pmax(ontxHAY,ontxY)]; PSA[,ontxHAO:=pmax(ontxHAO,ontxO)]

CFRdata <- PSA[,.(id,notxY,ontxY,
                  notxO,ontxO,
                  notxHAY,notxHAO,
                  ontxHAY,ontxHAO)] #,ontxHAY,ontxHAO
CFRdatam <- melt(CFRdata,id='id')
CFRdatam[,age:='2-4']; CFRdatam[grepl('Y',variable),age:='0-1']
CFRdatam[,variable:=gsub("Y$|O$","",variable)]
CFRdatam <- dcast(CFRdatam,id+age~variable,value.var = 'value')
CFRdatam[,c('dN','dA'):=.(notx-ontx,notxHA-ontxHA)] #delta-CFR

## RRs
if(SA=='pooled_eff'){
  INTE <- INTEP
} else {
  INTE <- INTE
}


## check
INTE[country!='Both',median(RR),by=country]
INTE[country!='Both',mean(RR),by=country]
summary(INTE[country=='Cameroon',.(RR)])
summary(INTE[country=='Kenya',.(RR)])
##NOTE means and medians very different important, in particular mean>1 in CMR if fitted to median

## changes in ATT success
nmz <- c('Cameroon', 'Kenya')
txsuccess <- data.table(
  country=nmz,
  BL=c(1, 0.818),
  INT=c(0.848, 0.976)
)


## ==== cascades and costs
DBC1 <- copy(DBC)  # assumes screening observed in TIPPI

# screening for TB symptoms sensitivity analysis
if(SA=='screen'){                        # assumes all children presenting are screened
  tmp <- cbind(DBC[metric=='Screened for TB symptoms',.(iso3, age, metric)],
               DBC[metric=='Initial care seeking',.(SOC,INT, ratio)])
  # tmp[,ratio:=INT/SOC]
  DBC1 <- rbind(DBC[metric!='Screened for TB symptoms'],
                tmp)
}

if(SA=='screenINT'){                     # assumes screening observed in TIPPI for SOC and all children presenting are screened for INT
  # DBC <- DBC
  tmp <- cbind(DBC[metric=='Screened for TB symptoms',.(iso3, age, metric, SOC)],
               DBC[metric=='Initial care seeking',.(INT)])
  tmp[,ratio:=INT/SOC]
  DBC1 <- rbind(DBC[metric!='Screened for TB symptoms'],
                tmp)
} 

# DBC1 <- copy(DBC)
K <- merge(ART2[,.(iso3, age, metric, uc.soc, uc.soc.sd, uc.int, uc.int.sd)],
           DBC1[,.(iso3, age, metric, frac=INT, ratio)],by=c('iso3','age','metric'),all.x=TRUE)
K[is.na(uc.soc.sd)|uc.soc.sd<0,uc.soc.sd:=uc.soc/40]
K[is.na(uc.int.sd),uc.int.sd:=uc.int/40]
K[is.na(ratio),ratio:=1] # TB tx or dx
K[is.na(frac),frac:=1]

## compute total costs & SD NOTE int not incremental in this analysis
K2 <- K[uc.soc>0,.(cost.soc=sum(frac*uc.soc/ratio),
                   cost.soc.sd=ssum(frac*uc.soc.sd/ratio),
                   cost.int=sum(frac*(uc.int)),
                   cost.int.sd=ssum(frac*uc.int.sd)),
        by=.(iso3,age)]
K2 <- merge(K2,CK,by='iso3')
K2 <- K2[,age:=ifelse(age=='u2', '0-1', '2-4')]
K2 #per child treated


## ================= MERGE TO CREATE MAIN DATA ===================
## main data for results
T <- data.table(
  id=rep(1:nrow(PSA), 4),
  age=rep(c('0-1','2-4'),each=nrow(PSA)),
  country=rep(c('Cameroon', 'Kenya'), each=nrow(PSA)*2)
)

T <- merge(T,CK,by=c('country'))

## ## merge in u2/o2 age split 
S <- ASM[qty=='present']
S <- S[,arm:=ifelse(arm=='SOC', 'Baseline', 'Intervention')]
S <- S[rep(1:nrow(S),each=max(T$id))]
S <- S[,frac:=rbeta(nrow(S),`0-1`,`2-4`)] #approx flat conjugate prior
S[,id:=rep(1:max(T$id),nrow(S)/max(T$id))]
S <- dcast(S[,.(iso3,id,arm,frac)],
           iso3+id~arm,value.var = 'frac')
T <- merge(T,S,by=c('iso3','id'),all.x=TRUE)
T[age=='2-4',Baseline:=1-Baseline] #defined as prop u5
T[age=='2-4',Intervention:=1-Intervention] #defined as prop u5
T[iso3=='CMR' & id==1] #check
names(T)[names(T)=='Baseline'] <- 'frac'
names(T)[names(T)=='Intervention'] <- 'fracI'

## merge in effects
if(length(unique(INTE$country))==1){                        
  T <- merge(T,
             INTE[,country:=NULL],
             by=c('id'))
} else {
  T <- merge(T,
             INTE[country!='Both',],
             by=c('id', 'country'))
}

## merge in costs
T <- merge(T,K2,by=c('iso3', 'country', 'age'))
## gamma samples for cost uncertainty (with safeties for 0)
T[,costv.soc:=rgamma(nrow(T),
                     shape=(cost.soc/(cost.soc.sd+1e-6))^2,
                     scale = cost.soc.sd^2/(cost.soc+1e-6))]
T[,costv.int:=rgamma(nrow(T),
                     shape=(cost.int/(cost.int.sd+1e-6))^2,
                     scale = cost.int.sd^2/(cost.int+1e-6))]

T[,c('costt.soc','costt.int'):=.(costv.soc,costv.int*RR)] # TODO: Check if this is not double counting 
T[,Dcost:=costt.int-costt.soc]

## merge in systematic review based CFR
T <- merge(T,CFRdatam[age=='0-1',.(id,notx, notxHA, ontx, ontxHA, dN, dA)],by=c('id')) #merge in CFR

## merge in trial-based CFR
T <- merge(T,TCFR,by=c('iso3','id'),all.x=TRUE)

## merge in tx success
# T <- merge(T,txsuccess,by='country')    #change in tx success
T <- merge(T,TSXS,by=c('iso3','id'),all.x=TRUE)

## merge in HIV
H
H[abs(hiv.sd)==0, hiv.sd:=hiv/40]
T <- merge(T,H[,.(iso3,hiv,hiv.sd)],by='iso3',all.x=TRUE)
T[,hiv:=rgamma(nrow(T),
               shape=(hiv/(hiv.sd+1e-6))^2,
               scale = hiv.sd^2/(hiv+1e-6))]
T[,.(hiv=mean(hiv), hiv.sd=mean(hiv.sd)), by = .(country)]

# # hiv from study
# T <- merge(T,HIVtrial[,.(id,iso3,age,hiv, hiv.int)],
#            by=c('id','age','iso3'),all.x = TRUE)

## merge in LYS
if(SA %in% c('hi','lo')){
  sa.drn <- ifelse(SA=='lo',0,5)
  fn <- gh('outdata/LYK_{sa.drn}.Rdata')
  load(fn)
}
T <- merge(T,LYK[,.(iso3,age,LYS,LYS0)],by=c('iso3','age'),all.x=TRUE)

## ================= CALCULATIONS ===================
T[,CFRnotx:=notx*(1-hiv) + notxHA*hiv]

# trial-based CFR
if(SA=='cfrltfu'){   
  T[,CFRtx:=cfr_int + cfr_int*ltfu_int]
  T[,CFRtx.soc:=cfr_soc + cfr_soc*ltfu_soc]
} else {
  T[,CFRtx:=cfr_int] 
  T[,CFRtx.soc:=cfr_soc]
}

# WHO-based CFR
if(SA=='whocfr'){   
  T[,CFRtx:=who_cfr_int]
  T[,CFRtx.soc:=who_cfr_soc]
} else {
  T[,CFRtx:=cfr_int] 
  T[,CFRtx.soc:=cfr_soc]
}

# literature-based CFR - original setup which had poor results for Kenya
if(SA=='systRev'){   
  T[,CFRtx:=ontx*(1-hiv) + ontxHA*hiv]
  T[,CFRtx.soc:=CFRtx]
} else {
  T[,CFRtx:=cfr_int] 
  T[,CFRtx.soc:=cfr_soc]
}

## NOTE this is important - if not tx mortality effect, KEN effect is -ve; with it +ve
## T[,CFRtx:=ontx*(1-hiv) + ontxHA*hiv]
## T[,CFRtx.soc:=CFRtx] 

## incremental lives saved - change in implied by RR
## before - 1 ontx:(RR-1) notx
## after - RR ontx: 0     notx
T[,deaths.int:=CFRtx*RR]         #deaths in intervention
T[,deaths.soc:=CFRtx.soc + (RR-1)*CFRnotx]         #deaths in soc
T[,LS:=deaths.soc-deaths.int]
T[,c('LYL.soc','LYL0.soc','LYL.int','LYL0.int'):=
    .(deaths.soc*LYS,deaths.soc*LYS0,deaths.int*LYS,deaths.int*LYS0)]
T[,c('dDALY','dDALY0'):=.(LYS*LS,LYS0*LS)]
T[,dCost:= costt.int-costt.soc]

T[,mean(RR),by=country]
T[,.(cfr_soc=mean(cfr_soc), cfr_int=mean(cfr_int),who_cfr_soc=mean(who_cfr_soc),who_cfr_int=mean(who_cfr_int)), by = .(iso3)]
T[,.(CFRnotx=mean(CFRnotx), CFRtx.soc=mean(CFRtx.soc), CFRtx=mean(CFRtx), deaths.soc=mean(deaths.soc), deaths.int=mean(deaths.int), LS=mean(LS)), by = .(country)]

(cfr <- T[,.(
  cfr.soc=mean(CFRtx.soc),
  cfr.soc.lo=lof(CFRtx.soc),
  cfr.soc.hi=hif(CFRtx.soc),
  cfr.int=mean(CFRtx),
  cfr.int.lo=lof(CFRtx),
  cfr.int.hi=hif(CFRtx)),
  by=.(iso3)])

(RR <- T[,.(
  RR=mean(RR),
  RR.lo=lof(RR),
  RR.hi=hif(RR)),
  by=.(iso3)])

txsuccess

## ================= COST-EFFECTIVENESS OUTPUT ===================
T[,.(dCost=mean(dCost),dDALY=mean(dDALY)),by=country]
T[,.(ICER=mean(dCost)/mean(dDALY)),by=country]

## CE plot
xyplot(data=T,dCost ~ dDALY | country,
       group=country,
       panel=function(...) {
         panel.xyplot(...)
         panel.abline(h=0)
         panel.abline(v=0)
       })

# consider putting this into another file
T <- merge(T,
           cascadetop,
           by=c('country', 'age'))

T1 <- T[,.(cost.soc=sum(costt.soc*frac),
           cost.int=sum(costt.int*frac),
           Dcost=sum(dCost*frac),
           ## deaths.soc=sum((deaths.int+LS)*frac),
           symptoms.soc = symptoms.soc,
           evaluated.soc = evaluated.soc,
           symptoms.int=sum(symptoms.soc*RRpresTB),
           evaluated.int=sum(evaluated.soc*RRevalTB),
           tx=sum(RR*frac),
           deaths.soc=sum(deaths.soc*frac),
           deaths.int=sum(deaths.int*frac),
           LYL.soc=sum(LYL.soc*frac),
           LYL.int=sum(LYL.int*frac),
           LYL0.soc=sum(LYL0.soc*frac),
           LYL0.int=sum(LYL0.int*frac),
           LS=sum(LS*frac),
           dDALY0=sum(dDALY0*frac),
           dDALY=sum(dDALY*frac)),
        by=.(country,id)]
## T[,mean(RR-1),by=.(country,age)]
T1 <- merge(T1,CK,by='country') #country iso3 merged on

## same with age - basically just consistent renaming
T2 <- T[,.(cost.soc=(costt.soc),
           cost.int=(costt.int),
           Dcost=(Dcost),
           symptoms.soc = symptoms.soc,
           evaluated.soc = evaluated.soc,
           symptoms.int=sum(symptoms.soc*RRpresTB),
           evaluated.int=sum(evaluated.soc*RRevalTB),
           tx = RR,
           ## deaths.soc=((deaths.int+LS)),
           deaths.soc=(deaths.soc),
           deaths.int=(deaths.int),
           LYL.soc=(LYL.soc),
           LYL.int=(LYL.int),
           LYL0.soc=(LYL0.soc),
           LYL0.int=(LYL0.int),
           LS=(LS),
           dDALY0=(dDALY0),
           dDALY=(dDALY)),
        by=.(country,id,age)]
T2 <- merge(T2,CK,by='country') #country iso3 merged on

## CEA plot
icer <- T1 %>% 
  group_by(country) %>% 
  summarise(cst = mean(Dcost),
            eff = mean(dDALY),
            icer = format(round(mean(Dcost)/mean(dDALY)), big.mark = ',', scientific=FALSE))

GP <- ggplot(T1,aes(dDALY,Dcost)) +
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  geom_point(alpha=0.1,shape=1) +
  geom_point(data = icer, aes(y = cst, x = eff, group=country), 
             size = 2, color = 'red', shape = 20) +
  geom_text(data = icer, aes(y = ifelse(cst>10000, 5e5, 2e3), x = ifelse(eff>10, 2000, 20), label=paste0('ICER=', 'US$',icer, '/DALY averted'), group=country), size = 2)+
  geom_abline(data=CETM[threshold %in% c("0.5x GDP","1x GDP")],
              aes(intercept=0,slope=value,col=threshold)) +
  facet_wrap(~country, scales = 'free') +
  scale_y_continuous(label=comma) +
  xlab('Discounted disability adjusted life-years averted')+
  ylab('Incremental cost (USD)')+
  theme(legend.position = "top")
if(!shell) GP

fn1 <- glue(here('plots/CEall')) + SA + '.png'
## fn2 <- glue(here('plots/CEall')) + SAT + '.pdf' 
ggsave(GP,file=fn1,w=6,h=5); ## ggsave(GP,file=fn2,w=10,h=10)
## PDF versions too big

## make CEAC data
lz <- seq(from = 0,to=ceactop*4,length.out = 1000)
ceacd <- list()
for(iso in unique(T1$iso3)){
  tmp <- T1[iso3==iso,.(Q=dDALY,P=Dcost)]
  ceacd[[iso]] <- data.table(
    iso3=iso,x=lz,
    y=make.ceac(tmp,lz))
}
ceacd <- rbindlist(ceacd)
ceacd[, iso3:=ifelse(iso3=='CMR', 'Cameroon', 'Kenya')]
## make CEAC plot
CEAC <- make.ceac.plot(ceacd,xpad=50) +
  theme(legend.position="top", legend.title = element_blank()) +
  xlab('Cost-effectiveness threshold (US$/DALY averted') +
  xlim(0,5000)
if(!shell) CEAC


fn1 <- glue(here('plots/CEAC')) + SA + '.png'
## fn2 <- glue(here('plots/CEAC')) + SAT + '.pdf'
ggsave(CEAC,file=fn1,w=7,h=7); ## ggsave(CEAC,file=fn2,w=7,h=7)

## ICERs by country
ice <- T1[,
          .(cost.soc=mean(cost.soc), #costs
            cost.soc.lo=lof(cost.soc),
            cost.soc.hi=hif(cost.soc),
            cost.int=mean(cost.int),
            cost.int.lo=lof(cost.int),
            cost.int.hi=hif(cost.int),
            Dcost=mean(Dcost),
            Dcost.lo=lof(Dcost),
            Dcost.hi=hif(Dcost),
            symptoms.soc =mean(symptoms.soc),
            evaluated.soc = mean(evaluated.soc),
            symptoms.int=mean(symptoms.int),
            symptoms.int.lo=lof(symptoms.int),
            symptoms.int.hi=hif(symptoms.int),
            evaluated.int=mean(evaluated.int),
            evaluated.int.lo=lof(evaluated.int),
            evaluated.int.hi=hif(evaluated.int),
            tx=mean(tx),
            tx.lo=lof(tx),
            tx.hi=hif(tx),
            ## deaths
            deaths.soc=mean(deaths.soc),
            deaths.soc.lo=lof(deaths.soc),
            deaths.soc.hi=hif(deaths.soc),
            LYL.soc=mean(LYL.soc),
            LYL.soc.lo=lof(LYL.soc),
            LYL.soc.hi=hif(LYL.soc),
            LYL0.soc=mean(LYL0.soc),
            LYL0.soc.lo=lof(LYL0.soc),
            LYL0.soc.hi=hif(LYL0.soc),
            deaths.int=mean(deaths.int),
            deaths.int.lo=lof(deaths.int),
            deaths.int.hi=hif(deaths.int),
            LYL.int=mean(LYL.int),
            LYL.int.lo=lof(LYL.int),
            LYL.int.hi=hif(LYL.int),
            LYL0.int=mean(LYL0.int),
            LYL0.int.lo=lof(LYL0.int),
            LYL0.int.hi=hif(LYL0.int),
            Ddeaths=mean(deaths.int-deaths.soc),
            Ddeaths.lo=lof(deaths.int-deaths.soc),
            Ddeaths.hi=hif(deaths.int-deaths.soc),
            LS=mean(LS),LS.lo=lof(LS),LS.hi=hif(LS),
            ## dalys
            dDALY0=mean(dDALY0),
            dDALY=mean(dDALY),
            # dDALY.nohiv=mean(dDALY.nohiv),
            dDALY.hi=hif(dDALY),
            dDALY.lo=lof(dDALY),
            dDALY0.hi=hif(dDALY0),
            dDALY0.lo=lof(dDALY0),
            ## ICER
            ICER=mean(Dcost)/mean(dDALY)
          ),
          by=country]

## multiply most by 100
vec <- names(ice)
vec <- vec[!vec %in% c('country','ICER')] #all but country and ICER
ice[,(vec):=lapply(.SD, function(x) 100*x), .SDcols = vec]

# Table1ATT
icer <- ice[,.(country=country,
               symptoms.soc =symptoms.soc,
               evaluated.soc = evaluated.soc,
               treated.soc=100,         #NOTE normalized to 100
               cost.soc = bracket(cost.soc,cost.soc.lo,cost.soc.hi),
               symptoms.int = bracket(symptoms.int, symptoms.int.lo, symptoms.int.hi),
               evaluated.int = bracket(evaluated.int, evaluated.int.lo, evaluated.int.hi),
               treated.int=bracket(tx,tx.lo,tx.hi),
               cost.int = bracket(cost.int,cost.int.lo,cost.int.hi),
               symptoms.dif=bracket(symptoms.int-symptoms.soc,symptoms.int.lo-symptoms.soc,symptoms.int.hi-symptoms.soc), #NOTE also 100
               evaluated.dif=bracket(evaluated.int-evaluated.soc,evaluated.int.lo-evaluated.soc,evaluated.int.hi-evaluated.soc), #NOTE also 100
               treated.dif=bracket(tx-1e2,tx.lo-1e2,tx.hi-1e2), #NOTE also 100
               cost.dif=bracket(Dcost,Dcost.lo,Dcost.hi),
               deaths.soc = bracket(deaths.soc,deaths.soc.lo,deaths.soc.hi),
               LYL.soc = bracket(LYL.soc,LYL.soc.lo,LYL.soc.hi),
               LYL0.soc = bracket(LYL0.soc,LYL0.soc.lo,LYL0.soc.hi),
               deaths.int = bracket(deaths.int,deaths.int.lo,deaths.int.hi),
               LYL.int = bracket(LYL.int,LYL.int.lo,LYL.int.hi),
               LYL0.int = bracket(LYL0.int,LYL0.int.lo,LYL0.int.hi),
               deaths.dif=bracket(-LS,-LS.hi,-LS.lo),
               LY0.dif=bracket(dDALY0,dDALY0.lo,dDALY0.hi),
               LY.dif=bracket(dDALY,dDALY.lo,dDALY.hi),
               ICER = format(round(ICER),big.mark = ',')
)]
icer

fn1 <- glue(here('outdata/ICERatt')) + SA + '.csv'
fwrite(icer,file=fn1)

##  export ICERS
iceb <- T1[,.(ICER=mean(Dcost)/mean(dDALY)),
           by=.(country,iso3)]

icebrr <- iceb[,.(country,iso3,
                  ICER = format(round(ICER),big.mark = ',') )]
icebrr

fn <- gh('outdata/ICERall') + SA + '.csv'
fwrite(icebrr,file=fn)

## --table 2 format
Table2ATT <- icer[,.(country,
                     symptoms.soc, evaluated.soc, treated.soc,cost.soc, 
                     symptoms.int,evaluated.int, treated.int, cost.int,
                     diff.symptoms=symptoms.dif, diff.evaluated=evaluated.dif,diff.treated=treated.dif,
                     deaths.soc, deaths.int,
                     diff.deaths=deaths.dif,
                     LYL.soc, LYL.int,
                     LYL0.soc, LYL0.int,
                     diff.LYL0=LY0.dif,diff.LYL=LY.dif,
                     diff.cost=cost.dif,ICER)]

fn1 <- glue(here('outdata/Table2ATT')) + SA + '.Rdata'
save(Table2ATT,file=fn1)

tmp <- melt(Table2ATT, id.vars = 'country')
soc <- tmp[grepl('.soc', variable),]
soc <- soc[,variable:=gsub('.soc', '', variable)]
names(soc)[3] <- 'Control'

int <- tmp[grepl('.int', variable),]
# int <- int[!grepl('ICER', variable),]
int <- int[,variable:=gsub('.int', '', variable)]
names(int)[3] <- 'Intervention'

d <- tmp[grepl('diff.|ICER', variable),]
d <- d[,variable:=gsub('diff.', '', variable)]
names(d)[3] <- 'Increment'

merged <- merge(soc, int, by=c('country', 'variable'), all.x = TRUE, all.y = TRUE)
merged <- merge(merged, d,by=c('country', 'variable'), all.x = TRUE, all.y = TRUE)

vars <- c('symptoms', 'evaluated', "treated", "deaths", "LYL0","LYL","cost","ICER")
var_labs <- c('symptoms', 'evaluated', "treated", "deaths", "LYL0","LYL","cost","ICER")

merged$variable <- factor(merged$variable, 
                          levels = vars)   
setorder(merged, country, variable)

Table2ATT_formatted <- cbind(merged[country=='Cameroon',.(variable, Control, Intervention, Increment)], merged[country=='Kenya', .(Control, Intervention, Increment)])
Table2ATT_formatted$variable <- var_labs
Table2ATT_formatted
fn1 <- glue(here('outdata/Table2ATT_formatted')) + SA + '.Rdata'
save(Table2ATT_formatted,file=fn1)