#Dependancies------
## Packages-------

# Install Packages if necessary
# install.packages("dplyr")
# install.packages("lme4")
# install.packages("glmulti")
# install.packages("lmerTest")
# install.packages("stringr")

library(dplyr)
library(visreg)
library(lme4)
library(glmulti)
library(lmerTest)
library(stringr)


## Functions-------

lmer.glmulti<-function(formula,data,random="",...) {
  newf <- formula
  newf[[3]] <- substitute(f+r,
                          list(f=newf[[3]],
                               r=reformulate(random)[[2]]))
  lmer(newf,data=data,
       REML=FALSE,...)
}
glmulti_topfnct<-function(res=res){
  top <- weightable(res)
  top <- top[top$aic <= min(top$aic) + 2,]
  print(top$model[1])
  return(top)
  remove(res,df,response,m1)
}

sigfig <- function(vec, n=3){ 
    formatC(signif(vec,digits=n), digits=n,format="fg", flag="#") 
}     
N<-function(x){length(x)}
se<-function(x){sd(x)/sqrt(length(x))}

# Load Data and directories-------

wdir<-"~/Dropbox/KU/Dimensions/RootBiomass"
wdname<-"~/Dropbox/KU/Dimensions/SoilHarvest/2022/Tables"

Biomass<-read.csv(file=paste(wdir,"AbovegroundBiomass2020-2022.csv",sep="/"))
Roots<-read.csv(file=paste(wdir,"Roots2020-2022.csv",sep="/"))
PlotTRMT<-read.csv(file=paste(wdir,"ExpDesign.csv",sep="/"))

PlotTRMT$RAINTRT<-as.factor(PlotTRMT$RAINTRT)
PlotTRMT$DISPERSION<-as.factor(PlotTRMT$DISPERSION)
PlotTRMT$BLOCK<-as.factor(PlotTRMT$BLOCK)
PlotTRMT$SUBBLOCK<-as.factor(PlotTRMT$SUBBLOCK)
PlotTRMT$DISPERSION.FOURLEVELS<-as.factor(PlotTRMT$DISPERSION.FOURLEVELS)

summary(PlotTRMT)
summary(Biomass)
summary(Roots)

#Data Joining-----

Biomass<-left_join(Biomass,PlotTRMT,by="PLOTNO")

Roots<-left_join(Roots,PlotTRMT,by="PLOTNO")

RootsTime<-subset(Roots,SPPNO!=0)

Roots20<-subset(RootsTime,YEAR==2020)
Roots22<-subset(RootsTime,YEAR==2022)

Biomass20<-subset(Biomass,YEAR==2020)
Biomass22<-subset(Biomass,YEAR==2022)

# Exploratory Data Analysis-----

par(mfrow=c(2,2))

hist((Roots20$DryWt.g.m2), las=1,font=2,xlab="",ylab="",main="")
box()
mtext(expression(paste(bold("Root Biomass 2020 (g m")^bold("-2"),bold(")"))),
      side=1,line=3,cex=0.95)
mtext(expression(paste(bold("Frequency"))),
      side=2,line=2.5,cex=0.95)
hist((Roots22$DryWt.g.m2), las=1,font=2,xlab="",ylab="",main="")
box()
mtext(expression(paste(bold("Root Biomass 2022 (g m")^bold("-2"),bold(")"))),
      side=1,line=3,cex=0.95)
mtext(expression(paste(bold("Frequency"))),
      side=2,line=2.5,cex=0.95)
hist(log10(Roots20$DryWt.g.m2), las=1,font=2,xlab="",ylab="",main="")
box()
mtext(expression(paste(bold("log10 Root Biomass 2020 (g m")^bold("-2"),bold(")"))),
      side=1,line=3,cex=0.95)
mtext(expression(paste(bold("Frequency"))),
      side=2,line=2.5,cex=0.95)
hist(log10(Roots22$DryWt.g.m2), las=1,font=2,xlab="",ylab="",main="")
box()
mtext(expression(paste(bold("log10 Root Biomass 2022 (g m")^bold("-2"),bold(")"))),
      side=1,line=3,cex=0.95)
mtext(expression(paste(bold("Frequency"))),
      side=2,line=2.5,cex=0.95)


# Root Glmulti--------

df=RootsTime
response=df$log10RootBiomass

m1 <- lmer(response ~ SPPNO*RAINTRT*DISPERSION*YEAR
           +(1|RAINTRT:SUBBLOCK)+(1|SUBBLOCK)+(1|PLOTNO)
           ,data=df)
res <- glmulti(formula(m1,fixed.only=TRUE),
               random="+ (1|RAINTRT:SUBBLOCK) + (1|SUBBLOCK) + (1|PLOTNO)",
               data=df,method="g",
               fitfunc=lmer.glmulti,
               intercept=TRUE,marginality=FALSE,level=2)
topRB.mod.D2<-glmulti_topfnct(res=res)
write.csv(topRB.mod.D2,file=paste(wdname,"topRB.mod.D2.csv",sep="/"))
remove(res,m1)

m1 <- lmer(response ~ SPPNO*RAINTRT*DISPERSION.FOURLEVELS*YEAR
           +(1|RAINTRT:SUBBLOCK)+(1|SUBBLOCK)+(1|PLOTNO)
           ,data=df)
res <- glmulti(formula(m1,fixed.only=TRUE),
               random="+ (1|RAINTRT:SUBBLOCK) + (1|SUBBLOCK) + (1|PLOTNO)",
               data=df,method="g",
               fitfunc=lmer.glmulti,
               intercept=TRUE,marginality=FALSE,level=2)
topRB.mod.D4<-glmulti_topfnct(res=res)
write.csv(topRB.mod.D4,file=paste(wdname,"topRB.mod.D4.csv",sep="/"))
remove(res,m1)

topRB.mod.D2<-read.csv(file=paste(wdname,"topRB.mod.D2.csv",sep="/"))
topRB.mod.D2<-topRB.mod.D2[,-1]
topRB.mod.D2[1,1]

topRB.mod.D4<-read.csv(file=paste(wdname,"topRB.mod.D4.csv",sep="/"))
topRB.mod.D4<-topRB.mod.D4[,-1]
topRB.mod.D4[1,1]

# Root Models--------

df=RootsTime
response=df$log10RootBiomass

formula.1<-paste(topRB.mod.D4$model[1],"(1|RAINTRT:SUBBLOCK) + (1|SUBBLOCK) + (1|PLOTNO)",sep=" +" )
formula.1 <- str_replace(formula.1, "response ~ ", "")
lmerrb.1<-lmer(paste("log10RootBiomass ~ ",formula.1,sep=""),data=df)

formula.2<-paste(topRB.mod.D4$model[2],"(1|RAINTRT:SUBBLOCK) + (1|SUBBLOCK) + (1|PLOTNO)",sep=" + ")
formula.2 <- str_replace(formula.2, "response ~ ", "")
lmerrb.2<-lmer(paste("log10RootBiomass ~ ",formula.2,sep=""),data=df)

formula.3<-paste(topRB.mod.D4$model[3],"(1|RAINTRT:SUBBLOCK) + (1|SUBBLOCK) + (1|PLOTNO)",sep=" + ")
formula.3 <- str_replace(formula.3, "response ~ ", "")
lmerrb.3<-lmer(paste("log10RootBiomass ~ ",formula.3,sep=""),data=df)

formula.4<-paste(topRB.mod.D4$model[4],"(1|RAINTRT:SUBBLOCK) + (1|SUBBLOCK) + (1|PLOTNO)",sep=" + ")
formula.4 <- str_replace(formula.4, "response ~ ", "")
lmerrb.4<-lmer(paste("log10RootBiomass ~ ",formula.4,sep=""),data=df)

formula.5<-paste(topRB.mod.D4$model[5],"(1|RAINTRT:SUBBLOCK) + (1|SUBBLOCK) + (1|PLOTNO)",sep=" + ")
formula.5 <- str_replace(formula.5, "response ~ ", "")
lmerrb.5<-lmer(paste("log10RootBiomass ~ ",formula.5,sep=""),data=df)

formula.1<-paste(topRB.mod.D2$model[1],"(1|RAINTRT:SUBBLOCK) + (1|SUBBLOCK) + (1|PLOTNO)",sep=" +" )
formula.1 <- str_replace(formula.1, "response ~ ", "")
lmerrb.1.D2<-lmer(paste("log10RootBiomass ~ ",formula.1,sep=""),data=df)

formula.2<-paste(topRB.mod.D2$model[2],"(1|RAINTRT:SUBBLOCK) + (1|SUBBLOCK) + (1|PLOTNO)",sep=" + ")
formula.2 <- str_replace(formula.2, "response ~ ", "")
lmerrb.2.D2<-lmer(paste("log10RootBiomass ~ ",formula.2,sep=""),data=df)

formula.3<-paste(topRB.mod.D2$model[3],"(1|RAINTRT:SUBBLOCK) + (1|SUBBLOCK) + (1|PLOTNO)",sep=" + ")
formula.3 <- str_replace(formula.3, "response ~ ", "")
lmerrb.3.D2<-lmer(paste("log10RootBiomass ~ ",formula.3,sep=""),data=df)

formula.4<-paste(topRB.mod.D2$model[4],"(1|RAINTRT:SUBBLOCK) + (1|SUBBLOCK) + (1|PLOTNO)",sep=" + ")
formula.4 <- str_replace(formula.4, "response ~ ", "")
lmerrb.4.D2<-lmer(paste("log10RootBiomass ~ ",formula.4,sep=""),data=df)

formula.5<-paste(topRB.mod.D2$model[5],"(1|RAINTRT:SUBBLOCK) + (1|SUBBLOCK) + (1|PLOTNO)",sep=" + ")
formula.5 <- str_replace(formula.5, "response ~ ", "")
lmerrb.5.D2<-lmer(paste("log10RootBiomass ~ ",formula.5,sep=""),data=df)

anova(lmerrb.1)
anova(lmerrb.2)
anova(lmerrb.3)
anova(lmerrb.4)
anova(lmerrb.5)

summary(lmerrb.1)
summary(lmerrb.2)
summary(lmerrb.3)
summary(lmerrb.4)
summary(lmerrb.5)

rand(lmerrb.1)
rand(lmerrb.2)
rand(lmerrb.3)
rand(lmerrb.4)
rand(lmerrb.5)


RootAICglmulti<-data.frame(c("lmerrb.1","lmerrb.2","lmerrb.3",
                             "lmerrb.4","lmerrb.5",
                             "lmerrb.1.D2","lmerrb.2.D2","lmerrb.3.D2",
                             "lmerrb.4.D2","lmerrb.5.D2"),
                           c(topRB.mod.D4$aic,topRB.mod.D2$aic))
names(RootAICglmulti)<-paste(c("Model","aic.glmulti"))

RootAICglmulti

#Aboveground Glmulti-------

hist(log10(Biomass$BMTOTPLOT.g))
Biomass$log10BMTOTPLOT.g<-log10(Biomass$BMTOTPLOT.g)

df=Biomass
response=df$log10BMTOTPLOT.g


m1 <- lmer(response ~ SPPNO*RAINTRT*DISPERSION.FOURLEVELS*YEAR
           +(1|RAINTRT:SUBBLOCK)+(1|SUBBLOCK)+(1|PLOTNO)
           ,data=df)
res <- glmulti(formula(m1,fixed.only=TRUE),
               random="+ (1|RAINTRT:SUBBLOCK) + (1|SUBBLOCK) + (1|PLOTNO)",
               data=df,method="g",
               fitfunc=lmer.glmulti,
               intercept=TRUE,marginality=FALSE,level=2)
topAB.mod.D4<-glmulti_topfnct(res=res)
write.csv(topAB.mod.D4,file=paste(wdname,"topAB.mod.D4.csv",sep="/"))
remove(res,m1)


m1 <- lmer(response ~ SPPNO*RAINTRT*DISPERSION*YEAR
           +(1|RAINTRT:SUBBLOCK)+(1|SUBBLOCK)+(1|PLOTNO)
           ,data=df)
res <- glmulti(formula(m1,fixed.only=TRUE),
               random="+ (1|RAINTRT:SUBBLOCK) + (1|SUBBLOCK) + (1|PLOTNO)",
               data=df,method="g",
               fitfunc=lmer.glmulti,
               intercept=TRUE,marginality=FALSE,level=2)
topAB.mod.D2<-glmulti_topfnct(res=res)
write.csv(topAB.mod.D2,file=paste(wdname,"topAB.mod.D2.csv",sep="/"))
remove(res,m1)

topAB.mod.D4<-read.csv(file=paste(wdname,"topAB.mod.D4.csv",sep="/"))
topAB.mod.D4<-topAB.mod.D4[,-1]
topAB.mod.D4[1,1]

topAB.mod.D2<-read.csv(file=paste(wdname,"topAB.mod.D2.csv",sep="/"))
topAB.mod.D2<-topAB.mod.D2[,-1]
topAB.mod.D2[1,1]


#Aboveground Models-------

df=Biomass
response=df$log10BMTOTPLOT.g

formula.1<-paste(topAB.mod.D4$model[1],"(1|RAINTRT:SUBBLOCK) + (1|SUBBLOCK) + (1|PLOTNO)",sep=" +" )
formula.1 <- str_replace(formula.1, "response ~ ", "")
lmerab.1<-lmer(paste("log10BMTOTPLOT.g ~ ",formula.1,sep=""),data=df)

formula.2<-paste(topAB.mod.D4$model[2],"(1|RAINTRT:SUBBLOCK) + (1|SUBBLOCK) + (1|PLOTNO)",sep=" +" )
formula.2 <- str_replace(formula.2, "response ~ ", "")
lmerab.2<-lmer(paste("log10BMTOTPLOT.g ~ ",formula.2,sep=""),data=df)

formula.3<-paste(topAB.mod.D4$model[3],"(1|RAINTRT:SUBBLOCK) + (1|SUBBLOCK) + (1|PLOTNO)",sep=" +" )
formula.3 <- str_replace(formula.3, "response ~ ", "")
lmerab.3<-lmer(paste("log10BMTOTPLOT.g ~ ",formula.3,sep=""),data=df)

formula.1<-paste(topAB.mod.D2$model[1],"(1|RAINTRT:SUBBLOCK) + (1|SUBBLOCK) + (1|PLOTNO)",sep=" +" )
formula.1 <- str_replace(formula.1, "response ~ ", "")
lmerab.1.D2<-lmer(paste("log10BMTOTPLOT.g ~ ",formula.1,sep=""),data=df)

formula.2<-paste(topAB.mod.D2$model[2],"(1|RAINTRT:SUBBLOCK) + (1|SUBBLOCK) + (1|PLOTNO)",sep=" +" )
formula.2 <- str_replace(formula.2, "response ~ ", "")
lmerab.2.D2<-lmer(paste("log10BMTOTPLOT.g ~ ",formula.2,sep=""),data=df)

formula.3<-paste(topAB.mod.D2$model[3],"(1|RAINTRT:SUBBLOCK) + (1|SUBBLOCK) + (1|PLOTNO)",sep=" +" )
formula.3 <- str_replace(formula.3, "response ~ ", "")
lmerab.3.D2<-lmer(paste("log10BMTOTPLOT.g ~ ",formula.3,sep=""),data=df)

formula.4<-paste(topAB.mod.D2$model[4],"(1|RAINTRT:SUBBLOCK) + (1|SUBBLOCK) + (1|PLOTNO)",sep=" +" )
formula.4 <- str_replace(formula.4, "response ~ ", "")
lmerab.4.D2<-lmer(paste("log10BMTOTPLOT.g ~ ",formula.4,sep=""),data=df)

formula.5<-paste(topAB.mod.D2$model[5],"(1|RAINTRT:SUBBLOCK) + (1|SUBBLOCK) + (1|PLOTNO)",sep=" +" )
formula.5 <- str_replace(formula.5, "response ~ ", "")
lmerab.5.D2<-lmer(paste("log10BMTOTPLOT.g ~ ",formula.5,sep=""),data=df)

anova(lmerab.1)
anova(lmerab.2)
anova(lmerab.3)

summary(lmerab.1)
summary(lmerab.2)
summary(lmerab.3)

rand(lmerab.1)
rand(lmerab.2)
rand(lmerab.3)

ShootAICglmulti<-data.frame(c("lmerab.1","lmerab.2","lmerab.3",
                              "lmerab.1.D2","lmerab.2.D2","lmerab.3.D2","lmerab.4.D2","lmerab.5.D2"),
                            c(topAB.mod.D4$aic,topAB.mod.D2$aic))
names(ShootAICglmulti)<-paste(c("Model","aic.glmulti"))

ShootAICglmulti

#Total Plot Prod--------------

RootsTime$PLOTYEAR<-as.factor(paste(RootsTime$PLOTNO,RootsTime$YEAR))
Biomass$PLOTYEAR<-as.factor(paste(Biomass$PLOTNO,Biomass$YEAR))

names(RootsTime)
names(Biomass)

kp1<-which(colnames(RootsTime)=="PLOTNO") 
kp2<-which(colnames(RootsTime)=="YEAR")
kp3<-which(colnames(RootsTime)=="PLOTYEAR")
kp4<-which(colnames(RootsTime)=="DryWt.g.m2") 
kp5<-which(colnames(RootsTime)=="log10RootBiomass") 
kp6<-which(colnames(Biomass)=="PLOTYEAR")
kp7<-which(colnames(Biomass)=="BMTOTPLOT.g") 
kp8<-which(colnames(Biomass)=="log10BMTOTPLOT.g") 

TotalBiomass<-left_join(RootsTime[,c(kp1,kp2,kp3,kp4,kp5)],Biomass[,c(kp6,kp7,kp8)],by="PLOTYEAR")
remove(kp1,kp2,kp3,kp4,kp5,kp6,kp7,kp8)

TotalBiomass<-left_join(TotalBiomass,PlotTRMT)
TotalBiomass$TotalBiomass<-TotalBiomass$BMTOTPLOT.g+TotalBiomass$DryWt.g.m2

hist(log10(TotalBiomass$TotalBiomass))
TotalBiomass$log10TotalBiomass<-log10(TotalBiomass$TotalBiomass)

#Total Plot Prod Glmulti--------------

df=TotalBiomass
response=df$log10TotalBiomass


m1 <- lmer(response ~ SPPNO*RAINTRT*DISPERSION*YEAR
           +(1|RAINTRT:SUBBLOCK)+(1|SUBBLOCK)+(1|PLOTNO)
           ,data=df)
res <- glmulti(formula(m1,fixed.only=TRUE),
               random="+ (1|RAINTRT:SUBBLOCK) + (1|SUBBLOCK) + (1|PLOTNO)",
               data=df,method="g",
               fitfunc=lmer.glmulti,
               intercept=TRUE,marginality=FALSE,level=2)
topTB.mod.D2<-glmulti_topfnct(res=res)
write.csv(topTB.mod.D2,file=paste(wdname,"topTB.mod.D2.csv",sep="/"))
remove(res,m1)

m1 <- lmer(response ~ SPPNO*RAINTRT*DISPERSION.FOURLEVELS*YEAR
           +(1|RAINTRT:SUBBLOCK)+(1|SUBBLOCK)+(1|PLOTNO)
           ,data=df)
res <- glmulti(formula(m1,fixed.only=TRUE),
               random="+ (1|RAINTRT:SUBBLOCK) + (1|SUBBLOCK) + (1|PLOTNO)",
               data=df,method="g",
               fitfunc=lmer.glmulti,
               intercept=TRUE,marginality=FALSE,level=2)
topTB.mod.D4<-glmulti_topfnct(res=res)
write.csv(topTB.mod.D4,file=paste(wdname,"topTB.mod.D4.csv",sep="/"))
remove(res,m1)

topTB.mod.D2<-read.csv(file=paste(wdname,"topTB.mod.D2.csv",sep="/"))
topTB.mod.D2<-topTB.mod.D2[,-1]
topTB.mod.D2[1,1]

topTB.mod.D4<-read.csv(file=paste(wdname,"topTB.mod.D4.csv",sep="/"))
topTB.mod.D4<-topTB.mod.D4[,-1]
topTB.mod.D4[1,1]

#Total Plot Prod Models--------------

df=TotalBiomass
response=df$log10TotalBiomass

formula.1<-paste(topTB.mod.D4$model[1],"(1|RAINTRT:SUBBLOCK) + (1|SUBBLOCK) + (1|PLOTNO)",sep=" +" )
formula.1 <- str_replace(formula.1, "response ~ ", "")
lmertb.1<-lmer(paste("log10TotalBiomass ~ ",formula.1,sep=""),data=df)

formula.2<-paste(topTB.mod.D4$model[2],"(1|RAINTRT:SUBBLOCK) + (1|SUBBLOCK) + (1|PLOTNO)",sep=" + ")
formula.2 <- str_replace(formula.2, "response ~ ", "")
lmertb.2<-lmer(paste("log10TotalBiomass ~ ",formula.2,sep=""),data=df)

formula.3<-paste(topTB.mod.D4$model[3],"(1|RAINTRT:SUBBLOCK) + (1|SUBBLOCK) + (1|PLOTNO)",sep=" + ")
formula.3 <- str_replace(formula.3, "response ~ ", "")
lmertb.3<-lmer(paste("log10TotalBiomass ~ ",formula.3,sep=""),data=df)

formula.4<-paste(topTB.mod.D4$model[4],"(1|RAINTRT:SUBBLOCK) + (1|SUBBLOCK) + (1|PLOTNO)",sep=" + ")
formula.4 <- str_replace(formula.4, "response ~ ", "")
lmertb.4<-lmer(paste("log10TotalBiomass ~ ",formula.4,sep=""),data=df)

formula.1<-paste(topTB.mod.D2$model[1],"(1|RAINTRT:SUBBLOCK) + (1|SUBBLOCK) + (1|PLOTNO)",sep=" + ")
formula.1 <- str_replace(formula.1, "response ~ ", "")
lmertb.1.D2<-lmer(paste("log10TotalBiomass ~ ",formula.1,sep=""),data=df)

formula.2<-paste(topTB.mod.D2$model[2],"(1|RAINTRT:SUBBLOCK) + (1|SUBBLOCK) + (1|PLOTNO)",sep=" + ")
formula.2 <- str_replace(formula.2, "response ~ ", "")
lmertb.2.D2<-lmer(paste("log10TotalBiomass ~ ",formula.2,sep=""),data=df)

formula.3<-paste(topTB.mod.D2$model[3],"(1|RAINTRT:SUBBLOCK) + (1|SUBBLOCK) + (1|PLOTNO)",sep=" + ")
formula.3 <- str_replace(formula.3, "response ~ ", "")
lmertb.3.D2<-lmer(paste("log10TotalBiomass ~ ",formula.3,sep=""),data=df)

formula.4<-paste(topTB.mod.D2$model[4],"(1|RAINTRT:SUBBLOCK) + (1|SUBBLOCK) + (1|PLOTNO)",sep=" + ")
formula.4 <- str_replace(formula.4, "response ~ ", "")
lmertb.4.D2<-lmer(paste("log10TotalBiomass ~ ",formula.4,sep=""),data=df)

formula.5<-paste(topTB.mod.D2$model[5],"(1|RAINTRT:SUBBLOCK) + (1|SUBBLOCK) + (1|PLOTNO)",sep=" + ")
formula.5 <- str_replace(formula.5, "response ~ ", "")
lmertb.5.D2<-lmer(paste("log10TotalBiomass ~ ",formula.5,sep=""),data=df)

formula.6<-paste(topTB.mod.D2$model[6],"(1|RAINTRT:SUBBLOCK) + (1|SUBBLOCK) + (1|PLOTNO)",sep=" + ")
formula.6 <- str_replace(formula.6, "response ~ ", "")
lmertb.6.D2<-lmer(paste("log10TotalBiomass ~ ",formula.6,sep=""),data=df)

formula.7<-paste(topTB.mod.D2$model[7],"(1|RAINTRT:SUBBLOCK) + (1|SUBBLOCK) + (1|PLOTNO)",sep=" + ")
formula.7 <- str_replace(formula.7, "response ~ ", "")
lmertb.7.D2<-lmer(paste("log10TotalBiomass ~ ",formula.7,sep=""),data=df)

anova(lmertb.1.D2)
anova(lmertb.2.D2)
anova(lmertb.3.D2)
anova(lmertb.4.D2)
anova(lmertb.5.D2)
anova(lmertb.6.D2)
anova(lmertb.7.D2)

summary(lmertb.1.D2)
summary(lmertb.2.D2)
summary(lmertb.3.D2)
summary(lmertb.4.D2)
summary(lmertb.5.D2)
summary(lmertb.6.D2)
summary(lmertb.7.D2)

rand(lmertb.1.D2)
rand(lmertb.2.D2)
rand(lmertb.3.D2)
rand(lmertb.4.D2)
rand(lmertb.5.D2)
rand(lmertb.6.D2)
rand(lmertb.7.D2)

anova(lmerrb.1)

TotalAICglmulti<-data.frame(c("lmertb.1","lmertb.2","lmertb.3","lmertb.4",
                              "lmertb.1.D2","lmertb.2.D2","lmertb.3.D2","lmertb.4.D2","lmertb.5.D2","lmertb.6.D2","lmertb.7.D2"),
                            c(topTB.mod.D4$aic,topTB.mod.D2$aic))
names(TotalAICglmulti)<-paste(c("Model","aic.glmulti"))

TotalAICglmulti

