#Load Packages-------
library(lme4)

#Load Data-------
FileName<-"SITable_k.csv"
WorkDir<-"~/"
kdf<-read.csv(file=paste(WorkDir,FileName,sep="/"))
remove(FileName,WorkDir)

str(kdf)
kdf$Rain<-as.factor(kdf$Rain)
kdf$Subblock<-as.factor(kdf$Subblock)
kdf$CommComp<-as.factor(kdf$CommComp)
kdf$Ast10<-as.factor(kdf$Ast10)
kdf$Fab10<-as.factor(kdf$Fab10)
kdf$Poa10<-as.factor(kdf$Poa10)
kdf$TwoSpp<-as.factor(kdf$TwoSpp)
kdf$TwoSppFab<-as.factor(kdf$TwoSppFab)

levels(kdf$CommComp)
kdf$CommComp <- 
  relevel(kdf$CommComp, ref = "OVER")

kdf2<-subset(kdf,Richness==2)

kdfmix<-subset(kdf,Richness!=1)

#LME 1: Richness, CommComp, Rain, Prop Panvir-------

anova(kd4lmermod<-lmer(k~Richness*CommComp*Rain
                       +PropPanvir
                       +(1|Subblock)
                       ,data=kdf))
summary(kd4lmermod)
rand(kd4lmermod)

#LME 2: Richness, Rain, Phylo Dist from Panvir-------

anova(kd4lmermod.Phylo<-lmer(k~Richness*Rain
                         +PhyloDistPv
                         +(1|Subblock)
                         ,data=kdf))

summary(kd4lmermod.Phylo)
rand(kd4lmermod.Phylo)

#Two Species Mixture ANOVA--------

kdf2aov2spp<-aov(k~TwoSpp,data=kdf2)

par(mfrow=c(1,2))
plot(kdf2aov2spp,which=c(1,2))
shapiro.test(resid(kdf2aov2spp))

anova(kdf2aov2spp)
summary(kdf2aov2spp)
TukeyHSD(kdf2aov2spp)


#Multiple Species Mixtures three-way ANOVA--------

kdfpresabsmod<-aov(k~Fab10+Poa10+Ast10,data=kdfmix)

par(mfrow=c(1,2))
plot(kdfpresabsmod,which=c(1,2))
shapiro.test(resid(kdfpresabsmod))

anova(kdfpresabsmod)
summary(kdfpresabsmod)
TukeyHSD(kdfpresabsmod)

