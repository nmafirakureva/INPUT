library(readxl)
library(data.table)
library(ggplot2)
library(scales)
library(glue)

## functions used in the other analyses

## utilities
rotx <- theme(axis.text.x = element_text(angle = 90, hjust = 1))
shhs <- data.table(sheet=1:3,age=c('014','04','514'),aged=c('0-14','0-4','5-14'))
setkey(shhs,sheet) ## sheet hash
gh <- function(x) glue(here(x))


## utility functions
hif <- function(x) quantile(x,probs = 0.975)
lof <- function(x) quantile(x,probs = 0.025)
logit <- function(x) log(x/(1-x))
ilogit <- function(x) 1/(1+exp(-x))
ssum <- function(x) sqrt(sum(x^2))

bracket <- function(x,y,z,ns=0,dp=0) paste0(format(round(x,dp),big.mark = ',',nsmall=ns),
                                            ' (',
                                            format(round(y,dp),big.mark = ',',nsmall=ns),
                                            ' - ',
                                            format(round(z,dp),big.mark = ',',nsmall=ns),
                                            ')')

bracket2 <- function(x,y,z,dp=0) paste0(formatC(round(x,dp),big.mark = ','),
                                            ' (',
                                            formatC(round(y,dp),big.mark = ','),
                                            ' - ',
                                            formatC(round(z,dp),big.mark = ','),
                                            ')')

rexel <- function(x,...) as.data.table(read_excel(x,...))

## drop & rename

droprename <- function(D,
                       nmz = c('Country','Facility',
                               'Baseline.Rate','Intervention.Rate',
                               'Baseline.Num','Intervention.Num')){
    ## drop % improvements
    (drop <- names(D)[5:6])
    D[,c(drop):=NULL]
    ## nmz <- c('Country','Facility','Baseline.Rate','Intervention.Rate',
    ##          'Baseline.Num','Intervention.Num')
    names(D)[1:6] <- nmz
    D <- D[Country!='India']
}



## NOTE this renames in the original data
dataprep <- function(D){
    DM <- melt(D[,.(Country,Facility,
                    Baseline.Rate,Intervention.Rate,
                    Baseline.Num,Intervention.Num)],
               id=c('Country','Facility'))
    DM[,gp:=paste0(Country,Facility)]
    DM[,qty:='Rate']
    DM[grepl('Num',variable),qty:='Number']
    DM[,variable:=gsub('\\.Rate|\\.Num','',variable)]
    DM <- dcast(DM,
                Country + Facility + variable + gp ~ qty,
                value.var = "value")
    DM[,value:=Rate]
    return(DM)
}

ptdataprep <- function(D){
  DM <- melt(D[,.(Country,Facility,Baseline,Intervention)],
             id=c('Country','Facility'))
  DM[,gp:=paste0(Country,Facility)]
  return(DM)
}


## plot function
intplot <- function(D,ylb='Monthly rate'){
    ggplot(D,aes(variable,value,group=gp,col=gp)) +
        geom_point(aes(size=Number),shape=1) +
        geom_line() +
        xlab('') +
        ylab(ylb) +
        facet_wrap(~Country,scales='free') +
        theme_classic()+
        guides(colour="none") + 
        theme(legend.position = 'bottom',
              strip.background = element_blank()) +
        ggpubr::grids()
}

## summary table
improvesummary <- function(D04,D514,D014){
    tmp <- data.table(age=c('04','514','014'),
                    Improving=c(D04[ ,1e2*mean(Intervention.Rate >
                                               Baseline.Rate)],
                                D514[ ,1e2*mean(Intervention.Rate >
                                                Baseline.Rate)],
                                D014[ ,1e2*mean(Intervention.Rate >
                                                Baseline.Rate)]),
                    Ratio=c(D04[Baseline.Rate>0,
                                median(Intervention.Rate/
                                               Baseline.Rate)],
                            D514[Baseline.Rate>0,
                                 median(Intervention.Rate/
                                      Baseline.Rate)],
                            D014[Baseline.Rate>0,
                                 median(Intervention.Rate/
                                      Baseline.Rate)]),
                    Difference=c(D04[,2*median((Intervention.Rate-Baseline.Rate)/(Intervention.Rate+Baseline.Rate))],
                                 D514[,2*median((Intervention.Rate-Baseline.Rate)/(Intervention.Rate+Baseline.Rate))],
                                 D014[,2*median((Intervention.Rate-Baseline.Rate)/(Intervention.Rate+Baseline.Rate))]))
    tmp[,c('Improving','Ratio','Difference'):=.(round(Improving),
                                                  round(Ratio,1),
                                                  round(Difference,1))]
    return(tmp)
}

## ----------------- functions for second set -------
droprename2 <- function(D,
                        nmz = c('Country','Facility',
                                'Baseline.Num','Baseline.FT',
                                'Intervention.Num','Intervention.FT')
                        ){
    keep <- names(D)[c(1:4,6:7)]
    D <- D[,..keep]
    names(D) <- nmz
    D <- D[Country!='India']
    D[,Baseline.Rate:=Baseline.Num/Baseline.FT]
    D[,Intervention.Rate:=Intervention.Num/Intervention.FT]
}


## NOTE this renames in the original data
dataprep <- function(D){
    DM <- melt(D[,.(Country,Facility,
                    Baseline.Rate,Intervention.Rate,
                    Baseline.Num,Intervention.Num)],
               id=c('Country','Facility'))
    DM[,gp:=paste0(Country,Facility)]
    DM[,qty:='Rate']
    DM[grepl('Num',variable),qty:='Number']
    DM[,variable:=gsub('\\.Rate|\\.Num','',variable)]
    DM <- dcast(DM,
                Country + Facility + variable + gp ~ qty,
                value.var = "value")
    DM[,value:=Rate]
    return(DM)
}


make.ceac <- function(CEA,lamz){
    crv <- lamz
    for(i in 1:length(crv)) crv[i] <- CEA[,mean(lamz[i]*Q-P>0)]
    crv
}

make.ceac.plot <- function(D,thresholds=NULL,xpad=0){
    ## cols
    cbPalette0 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
                    "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
                    "darkorchid1")
    cbPalette <- c("#999999", "#E69F00", "#56B4E9","#009E73",
                   "#F0E442", "#0072B2","#D55E00", "#CC79A7",
                   "darkorchid1")
    ## plot
    GP <- ggplot(D,aes(x,y,col=iso3)) +
        geom_line() +
        theme_classic() +
        theme(legend.position = 'top')+
        scale_y_continuous(label=percent,
                           limits=c(0,1),
                           expand = c(0,0))+
        scale_x_continuous(label=comma,
                           limits=c(0,max(D$x)+1e3),
                           expand = c(0,xpad),
                           breaks = seq(0,max(D$x),2000))+
        ## scale_color_colorblind()+
        scale_colour_manual(values=cbPalette)+
        ylab('Probability cost-effective')+
        xlab('Cost-effectiveness threshold (USD/DALY averted)')+
        ggpubr::grids()
        ## refs
        if(!is.null(thresholds))
            GP <- GP + geom_vline(data=thresholds,aes(xintercept=x,
                                                      col=iso3),lty=2)
    GP
}
