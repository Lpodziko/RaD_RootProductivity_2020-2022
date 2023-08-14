#Load Packages-------
library(lme4)
library(lmerTest)

#Load Data-------
FileName<-"SITable_k.csv"
WorkDir<-"~/Dropbox/KU/Dimensions/aManuscriptPreparations/LitterDecomp-Dim4a-2020"
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
kdf$CommComp.rl <- 
  relevel(kdf$CommComp, ref = "OVER")
str(kdf)

kdf2<-subset(kdf,Richness==2)

kdfmix<-subset(kdf,Richness!=1)

kdfmixPOA<-subset(kdfmix,CommComp=="POA")

kdfmixFAB<-subset(kdfmix,CommComp=="FAB")

#Experimental Design models-------
##LME 1: Richness, CommComp, Rain, Prop Panvir-------

anova(lme1<-lmer(k~Richness*CommComp.rl*Rain
                       +PropPanvir.Rlzd
                       +(1|Subblock)
                       ,data=kdf))
summary(lme1)
rand(lme1)

##LME 2: Richness, Rain, Phylo Dist from Panvir-------

anova(lme2<-lmer(k~Richness*Rain
                         +PhyloDistPv
                         +(1|Subblock)
                         ,data=kdf))

summary(lme2)
rand(lme2)

##AOV 1: Two Species Mixture ANOVA--------

aov1<-aov(k~TwoSpp,data=kdf2)

par(mfrow=c(1,2))
plot(aov1,which=c(1,2))
shapiro.test(resid(aov1))

anova(aov1)
summary(aov1)
TukeyHSD(aov1)


##AOV 2: Multiple Species Mixtures three-way ANOVA--------

aov2<-aov(k~Fab10+Poa10+Ast10,data=kdfmix)

par(mfrow=c(1,2))
plot(aov2,which=c(1,2))
shapiro.test(resid(aov2))

anova(aov2)
summary(aov2)
TukeyHSD(aov2)

#Microbial Community-------
##LME 3-6: Microbial Diversity-----

lme3<-lmer(k~ShannonsH.FungSapr+(1|Subblock),data=kdf)
anova(lme3)
summary(lme3)

lme4<-lmer(k~Richness.FungSapr+(1|Subblock),data=kdf)
anova(lme4)
summary(lme4)

lme5<-lmer(k~ShannonsH.Bacteria+(1|Subblock),data=kdf)
anova(lme5)
summary(lme5)

lme6<-lmer(k~Richness.Bacteria+(1|Subblock),data=kdf)
anova(lme6)
summary(lme6)

##Piecewise Reg: Microbial Dissimilarity from switchgrass-----

piecewise1 <- lm(k ~ AvgDistFung*(PropPanvir.Pntd > 0) + AvgDistFung*(PropPanvir.Pntd == 0),data=kdf)
anova(piecewise1)
summary(piecewise1)

piecewise2 <- lm(k ~ AvgDistBact*(PropPanvir.Pntd > 0) + AvgDistBact*(PropPanvir.Pntd == 0),data=kdf)
anova(piecewise2)
summary(piecewise2)

#Litter Chemistry-------
##LME 7: N difference--------

trnsno<-(-1*floor(min(kdf$Ndiff,na.rm=T)))
lme7<-lmer(log(Ndiff+trnsno)~scale(k)*CommComp+(1|Subblock),data=kdf)

anova(lme7)
summary(lme7)
rand(lme7)
remove(trnsno)

##LME 8: C difference--------

trnsno<-(-1*floor(min(kdf$Cdiff,na.rm=T)))
lme8<-lmer(log(Cdiff+trnsno)~scale(k)*CommComp+(1|Subblock),data=kdf)

anova(lme8)
summary(lme8)
rand(lme8)
remove(trnsno)

##LME 9: CN difference--------

lme9<-lmer(CNdiff~scale(k)*CommComp+(1|Subblock),data=kdf)

anova(lme9)
summary(lme9)
rand(lme9)

#Overyielding by k-------
##LME 10: Grasses----

lme10<-lmer(CE~scale(k)*PropPanvir.Pntd+(1|Subblock),data=kdfmixPOA)

anova(lme10)
summary(lme10)
rand(lme10)

##LME 11: Legumes----

lme11<-lmer(CE~scale(k)+(1|Subblock),data=kdfmixFAB)

anova(lme11)
summary(lme11)
rand(lme11)

