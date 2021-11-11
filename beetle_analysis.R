#clear workspace
rm(list=ls())

#open packages
library(readr)       #to import and read .csv file
library(dplyr)       #data manipulation (inc %>% function)
library(performance) # for checking model quality
#from Tom's notes -->
library(car)         #needed to use "Anova"
library(glmmTMB)     #to run zero-inflated GLMMs
library(DHARMa)      #to check model residuals
library(ggplot2)     #for graphing results
library(multcomp)    #to run post-hoc tests (pairwise comparisons)

#import data
beetles <- read_csv("~/University/ASAB beetle project/beetle_data.csv")
View(beetles)

#create new columns for direct care
#female
beetles<- beetles %>% 
  mutate(f_direct_care = f_freq_f)
#male
beetles<- beetles %>%
  mutate(m_direct_care = m_freq_f)

#create new columns for indirect care
#female
beetles<- beetles %>%
  mutate(f_indirect_care = f_freq_c + f_freq_m)
#male
beetles<- beetles %>%
  mutate(m_indirect_care = m_freq_c + m_freq_m)

#does week of observation correlate with level of care?
#female direct
plot(f_direct_care~block_id, data=beetles)    
abline(lm(beetles$f_direct_care~beetles$block_id),col="red",lwd=3)  
cor.test(beetles$f_direct_care,beetles$block_id, method="pearson")
#females provided more direct care in later weeks

#female indirect
plot(f_indirect_care~block_id, data=beetles)
abline(lm(beetles$f_indirect_care~beetles$block_id),col="red",lwd=3)
cor.test(beetles$f_indirect_care,beetles$block_id, method="pearson")

#female food provisioning
plot(f_freq_f~block_id, data=beetles)
abline(lm(beetles$f_freq_f~beetles$block_id),col="red",lwd=3)
cor.test(beetles$f_freq_f,beetles$block_id, method="pearson")

#female carrion consumption
plot(f_freq_c~block_id, data=beetles)
abline(lm(beetles$f_freq_c~beetles$block_id),col="red",lwd=3)
cor.test(beetles$f_freq_c,beetles$block_id, method="pearson")

#male direct
plot(m_direct_care~block_id, data=beetles)
abline(lm(beetles$m_direct_care~beetles$block_id),col="red",lwd=3)
cor.test(beetles$m_direct_care,beetles$block_id, method="pearson")

#male indirect
plot(m_indirect_care~block_id, data=beetles)
abline(lm(beetles$m_indirect_care~beetles$block_id),col="red",lwd=3)
cor.test(beetles$m_indirect_care,beetles$block_id, method="pearson")

#testing basic correlations
#female direct care and handicap
plot(f_direct_care~female_treatment, data=beetles)
abline(lm(beetles$f_direct_care~beetles$female_treatment),col="red",lwd=3)
cor.test(beetles$f_direct_care,beetles$female_treatment, method="pearson")

#male care and female condition
plot(m_direct_care~female_treatment, data=beetles)
abline(lm(beetles$m_direct_care~beetles$female_treatment),col="red",lwd=3)
cor.test(beetles$m_direct_care,beetles$female_treatment, method="pearson")

#male presence and larvae number (at obs)
plot(obs_larvae_no~male_treatment, data=beetles)
abline(lm(beetles$obs_larvae_no~beetles$male_treatment),col="red",lwd=3)
cor.test(beetles$obs_larvae_no,beetles$male_treatment, method="pearson")

#female treatment and larvae number (obs)
plot(obs_larvae_no~female_treatment, data=beetles)
abline(lm(beetles$obs_larvae_no~beetles$female_treatment),col="red",lwd=3)
cor.test(beetles$obs_larvae_no,beetles$female_treatment, method="pearson")

#male presence and larvae number (at disp)
plot(disp_larvae_no~male_treatment, data=beetles)
abline(lm(beetles$disp_larvae_no~beetles$male_treatment),col="red",lwd=3)
cor.test(beetles$disp_larvae_no,beetles$male_treatment, method="pearson")

#female treatment and larvae number (at disp)
plot(disp_larvae_no~female_treatment, data=beetles)
abline(lm(beetles$disp_larvae_no~beetles$female_treatment),col="red",lwd=3)
cor.test(beetles$disp_larvae_no,beetles$female_treatment, method="pearson")


# ---------- working through Tom's notes ------------- #

#check number of zeros (female direct care)
table(beetles$f_direct_care==0)
#doesn't seem to be zero inflated (but what is the expectation?)

#view this as histogram
hist(beetles$f_direct_care)

#check number of zeros (male direct care)
table(beetles$m_direct_care==0)
#seems very zero inflated

#view as histogram
hist(beetles$m_direct_care)

# MODEL 1 - binomial model
#add column for number of counts females were NOT provisioning
beetles<- beetles %>% 
  mutate(f_freq_not_f = 30 - f_freq_f)

#run simple model
mod1<- glmmTMB(cbind(f_freq_f, f_freq_not_f)~female_treatment,
               data=beetles,family="binomial")
summary(mod1)

#check residuals
plot(simulateResiduals(mod1))
#left hand plot indicates likely problem with overdispersion

# MODEL 2 - overdispersed binomial model
#deal with overdispersion by includig individual level random effect
beetles$brood_id<- as.factor(1:nrow(beetles))
View(beetles$brood_id) #this has numbered brood_ID rather than using letter codes

mod1.2<- glmmTMB(cbind(f_freq_f, f_freq_not_f)~female_treatment+(1|brood_id),
                 data=beetles,family="binomial")
plot(simulateResiduals(mod1.2))
#this model does a better job at dealing with overdispersion (left hand plot)

#how well does model fit the data? - simulate data using the model and compare to original data
sim_max<- apply(simulate(mod1.2, nsim = 1000), 2, max)
hist(sim_max, breaks=10) 
abline(v=max(beetles$f_freq_f, na.rm=T),col="red",lwd=2)
#this seems okay (line falls within predicted values)

#check how good model is at predicting number of zeros in data
simy<- simulate(mod1.2, 1000)
nz<- c()  
for(i in seq(1, length(simy),1)){nz[i]<- colSums(simy[,i]==0)[1]}  
hist(nz, main="Number of zeros")
abline(v=sum(beetles$f_freq_f==0,na.rm=T),col="red",lwd=2)
#red line is actual number of zeros
#seems to estimate number of zeros okay (what does spike at 30 mean?)

# MODEL 3 - zero-inflated, overdispersed binomial model
mod1.3<- glmmTMB(cbind(f_freq_f, f_freq_not_f)~female_treatment+(1|brood_id),
                 data=beetles, ziformula=~1, family="binomial")
summary(mod1.3)
plot(simulateResiduals(mod1.3))  #deals well with overdispersion

sim_max<- apply(simulate(mod1.3, nsim = 1000),2,max)
hist(sim_max, breaks=10)
abline(v=max(beetles$f_freq_f,na.rm=T),col="red",lwd=2)
#is this too far to one side?

simy<- simulate(mod1.3, 1000)
nz<- c()  
for(i in seq(1, length(simy), 1)) {nz[i]<- colSums(simy[,i]==0)[1]} 
hist(nz, main="Number of zeros")
abline(v=sum(beetles$f_freq_f==0, na.rm=T),col="red",lwd=2)  
#this predicts number of zeros very well

# ANOVA - used to test overall effect of treatment
Anova(mod1.3)
summary(mod1.3)  

# PAIRWISE COMPARISONS
#female treatment has to be formatted as a factor
beetles$female_treatment<- as.factor(beetles$female_treatment)
#re-run model 1.3 with treatment as a factor 

summary(glht(mod1.3, mcp(female_treatment="Tukey")), test=adjusted("bonferroni"))


# --------- PRACTICE MODELS -----------  

# MODEL 2 including female and male treatment INCLUDING BLOCK
# focusing on female response
mod2<- glmmTMB(cbind(f_freq_f, f_freq_not_f)~female_treatment+male_treatment+
                 female_treatment*male_treatment+(1|block_id),data=beetles,family="binomial")
summary(mod2)
plot(simulateResiduals(mod2))  #overdispersed

#deal with overdispersal by including (1|brood_id)
mod2.2<- glmmTMB(cbind(f_freq_f,f_freq_not_f)~female_treatment+male_treatment+
                   female_treatment*male_treatment+(1|block_id)+(1|brood_id),
                 data=beetles,family="binomial")
summary(mod2.2)  
plot(simulateResiduals(mod2.2))  #deals with overdispersal

#how well does it fit the data
sim_max<- apply(simulate(mod2.2, nsim = 1000), 2, max)
hist(sim_max, breaks=10) 
abline(v=max(beetles$f_freq_f, na.rm=T),col="red",lwd=2)
#this seems okay (line falls within predicted values)

#check how good model is at predicting number of zeros in data
simy<- simulate(mod2.2, 1000)
nz<- c()  
for(i in seq(1, length(simy),1)){nz[i]<- colSums(simy[,i]==0)[1]}  
hist(nz, main="Number of zeros")
abline(v=sum(beetles$f_freq_f==0,na.rm=T),col="red",lwd=2)
#very not good at predicting zeros

mod2.3<- glmmTMB(cbind(f_freq_f,f_freq_not_f)~female_treatment+male_treatment+
                   female_treatment*male_treatment+(1|block_id)+(1|brood_id),
                 data=beetles,ziformula=~1,family="binomial")

#how well does it deal with overdispersal?  
plot(simulateResiduals(mod2.3)) #good

#how well does it fit the data?
sim_max<- apply(simulate(mod2.3, nsim = 1000), 2, max)
hist(sim_max, breaks=10) 
abline(v=max(beetles$f_freq_f, na.rm=T),col="red",lwd=2)
#good

#how well does it predict number of zeros
simy<- simulate(mod2.3, 1000)
nz<- c()  
for(i in seq(1, length(simy),1)){nz[i]<- colSums(simy[,i]==0)[1]}  
hist(nz, main="Number of zeros")
abline(v=sum(beetles$f_freq_f==0,na.rm=T),col="red",lwd=2)
#better than previous model but not great

#model outputs
Anova(mod2.3)
summary(mod2.3)  

# MODEL 3 including female and male treatment NOT INCLUDING BLOCK
# focusing on female response
mod3<- glmmTMB(cbind(f_freq_f,f_freq_not_f)~female_treatment+male_treatment+
                 female_treatment*male_treatment,data=beetles,family="binomial")
summary(mod3)
plot(simulateResiduals(mod3))  #overdispersed

#deal with overdispersal by including (1|brood_id)
mod3.2<- glmmTMB(cbind(f_freq_f,f_freq_not_f)~female_treatment+male_treatment+
                   female_treatment*male_treatment+(1|brood_id),data=beetles,family="binomial")
summary(mod3.2)  
plot(simulateResiduals(mod3.2))  #deals with overdispersal

#how well does it fit the data
sim_max<- apply(simulate(mod3.2, nsim = 1000), 2, max)
hist(sim_max, breaks=10) 
abline(v=max(beetles$f_freq_f, na.rm=T),col="red",lwd=2)
#this seems okay (line falls within predicted values)

#how well does it predict number of zeros
simy<- simulate(mod3.2, 1000)
nz<- c()  
for(i in seq(1, length(simy),1)){nz[i]<- colSums(simy[,i]==0)[1]}  
hist(nz, main="Number of zeros")
abline(v=sum(beetles$f_freq_f==0,na.rm=T),col="red",lwd=2)
#very bad at predicting zeros

#include zero inflation model
mod3.3<- glmmTMB(cbind(f_freq_f,f_freq_not_f)~female_treatment+male_treatment+
                   female_treatment*male_treatment+(1|brood_id),
                 data=beetles,ziformula=~1,family="binomial")

#how well does it deal with overdispersal?  
plot(simulateResiduals(mod3.3)) #good

#how well does it fit the data?
sim_max<- apply(simulate(mod3.3, nsim = 1000), 2, max)
hist(sim_max, breaks=10) 
abline(v=max(beetles$f_freq_f, na.rm=T),col="red",lwd=2)
#good

#how well does it predict number of zeros
simy<- simulate(mod3.3, 1000)
nz<- c()  
for(i in seq(1, length(simy),1)){nz[i]<- colSums(simy[,i]==0)[1]}  
hist(nz, main="Number of zeros")
abline(v=sum(beetles$f_freq_f==0,na.rm=T),col="red",lwd=2)
#much better - almost perfect

#model outputs
Anova(mod3.3)
summary(mod3.3)

# not including (1|block_id) in subsequent models --> randomization means block just adds noise to data, no systematic bias

#------- prep for final models ---------

#add column for number of counts females were not providing direct care
beetles<- beetles %>% 
  mutate(ffreq_no_direct = 30 - f_direct_care)

#add column for number of counts females were not providing indirect care
beetles<- beetles %>%
  mutate(ffreq_no_indirect = 30 - f_indirect_care)

#add column for number of counts males were not providing direct care
beetles<- beetles %>%
  mutate(mfreq_no_direct = 30 - m_direct_care)

#add column for number of counts males were not providing indirect care
beetles<- beetles %>%
  mutate(mfreq_no_indirect = 30 - m_indirect_care)


#---------- FINAL MODELS -------------

#---- MODEL 4 - FEMALE DIRECT CARE -------

#plot data
hist(beetles$f_direct_care) #lots of zeros

mod4.1<- glmmTMB(cbind(f_direct_care,ffreq_no_direct)~female_treatment+male_treatment+
                   female_treatment*male_treatment,data=beetles,family="binomial")

#check overdispersal
plot(simulateResiduals(mod4.1)) #overdispersed

#deal with overdispersion by including individual level random effect
mod4.2<- glmmTMB(cbind(f_direct_care,ffreq_no_direct)~female_treatment+male_treatment+
                   female_treatment*male_treatment+(1|brood_id),data=beetles,family="binomial")

#check overdispersal
plot(simulateResiduals(mod4.2)) #better at dealing with overdispersal

#how well does model fit the data?
sim_max<- apply(simulate(mod4.2,nsim=1000),2,max)
hist(sim_max,breaks=10)
abline(v=max(beetles$f_direct_care,na.rm=T),col="red",lwd=2)
#seems okay - line falls within predicted values

#how well does model predict zeros
simy<- simulate(mod4.2,1000)
nz<- c()
for(i in seq(1,length(simy),1)){nz[i]<- colSums(simy[,i]==0)[1]}
hist(nz,main="Number of zeros")
abline(v=sum(beetles$f_direct_care==0,na.rm=T),col="red",lwd=2)
#very bad at predicting zeros

#include zero inflation model
mod4.3<- glmmTMB(cbind(f_direct_care,ffreq_no_direct)~female_treatment+male_treatment+
                   female_treatment*male_treatment+(1|brood_id),
                 ziformula=~1,data=beetles,family=binomial)

#ceck overdispersion
plot(simulateResiduals(mod4.3))

#how well does model fit the data?
sim_max<- apply(simulate(mod4.3,nsim=1000),2,max)
hist(sim_max,breaks=10)
abline(v=max(beetles$f_direct_care,na.rm=T),col="red",lwd=2)

#how well does it predict zeros?
simy<- simulate(mod4.3,1000)
nz<- c()
for(i in seq(1,length(simy),1)){nz[i]<- colSums(simy[,i]==0)[1]}
hist(nz,main="Number of zeros")
abline(v=sum(beetles$f_direct_care==0,na.rm=T),col="red",lwd=2)
#predicts zeros well

#model outputs
Anova(mod4.3)
summary(mod4.3)

#pairwise comparisons
summary(glht(mod4.3,mcp(female_treatment="Tukey")),test=adjusted("bonferroni"))
#no effect of treatment on female direct care?

#plot results (female treatment)
plot_mod4_f<- ggplot(subset(beetles, !is.na(f_direct_care)),aes(x=female_treatment,y=f_direct_care))+
  geom_dotplot(binaxis="y",stackdir="center",alpha=0.6,dotsize=0.5,
               position=position_dodge(0.8))+
  stat_summary(fun.data=mean_se,geom="errorbar",fun.args=list(mult=1),
               colour="black",width=0.2,position=position_dodge(0.8))+
  stat_summary(fun=mean,geom="point",shape=21,fill="grey",
               size=3,position=position_dodge(0.8))+
  ylab("Direct care (number of scans)")+
  xlab("Female Treatment")+
  theme_bw()
  
plot_mod4_f

#format male treatment so the plot works
str(beetles$female_treatment)
str(beetles$male_treatment)
beetles$male_treatment<- as.factor(beetles$male_treatment)
str(beetles$male_treatment)

#plot results (male presence)
plot_mod4_m<- ggplot(subset(beetles, !is.na(f_direct_care)),aes(x=male_treatment,y=f_direct_care))+
  geom_dotplot(binaxis="y",stackdir="center",alpha=0.6,dotsize=0.5,
               position=position_dodge(0.8))+
  stat_summary(fun.data=mean_se,geom="errorbar",fun.args=list(mult=1),
               colour="black",width=0.2,position=position_dodge(0.8))+
  stat_summary(fun=mean,geom="point",shape=21,fill="grey",
               size=3,position=position_dodge(0.8))+
  ylab("Direct care (number of scans)")+
  xlab("Male presence")+
  theme_bw()

plot_mod4_m  

#------ MODEL 5 - FEMALE INDIRECT CARE ---------
hist(beetles$f_indirect_care) #lots of zeros but otherwise normal distribution

mod5.1<- glmmTMB(cbind(f_indirect_care,ffreq_no_indirect)~female_treatment+male_treatment+
                   female_treatment*male_treatment,data=beetles,family="binomial")

summary(mod5.1) #looks like male treatment may affect level of indirect care female provides

#check for overdispersal
plot(simulateResiduals(mod5.1)) #looks overdispersed

#include individual level random effect to deal with overdispersal
mod5.2<- glmmTMB(cbind(f_indirect_care,ffreq_no_indirect)~female_treatment+male_treatment+
                   female_treatment*male_treatment+(1|brood_id),data=beetles,family="binomial")

plot(simulateResiduals(mod5.2)) #better

#how does it fit the data
sim_max<- apply(simulate(mod5.2,nsim=1000),2,max)
hist(sim_max,breaks=10)
abline(v=max(beetles$f_indirect_care,na.rm=T),col="red",lwd=2)
#seems okay - line falls within predicted values (but only just)

#how well does model predict zeros
simy<- simulate(mod5.2,1000)
nz<- c()
for(i in seq(1,length(simy),1)){nz[i]<- colSums(simy[,i]==0)[1]}
hist(nz,main="Number of zeros")
abline(v=sum(beetles$f_indirect_care==0,na.rm=T),col="red",lwd=2)
#underestimates number of zeros

#include ziformula
mod5.3<- glmmTMB(cbind(f_indirect_care,ffreq_no_indirect)~female_treatment+male_treatment+
                   female_treatment*male_treatment+(1|brood_id),ziformula=~1,
                 data=beetles,family="binomial")

plot(simulateResiduals(mod5.3))

#how well does model fit the data?
sim_max<- apply(simulate(mod5.3,nsim=1000),2,max)
hist(sim_max,breaks=10)
abline(v=max(beetles$f_indirect_care,na.rm=T),col="red",lwd=2)
#good

#how well does it predict zeros?
simy<- simulate(mod5.3,1000)
nz<- c()
for(i in seq(1,length(simy),1)){nz[i]<- colSums(simy[,i]==0)[1]}
hist(nz,main="Number of zeros")
abline(v=sum(beetles$f_indirect_care==0,na.rm=T),col="red",lwd=2)
#predicts zeros well

#model outputs
Anova(mod5.3)
summary(mod5.3)

#pairwise comparison
summary(glht(mod5.3,mcp(female_treatment="Tukey")),test=adjusted("bonferroni"))

#plot results (female treatment)
plot_mod5_f<- ggplot(subset(beetles, !is.na(f_indirect_care)),aes(x=female_treatment,y=f_indirect_care))+
  geom_dotplot(binaxis="y",stackdir="center",alpha=0.6,dotsize=0.5,
               position=position_dodge(0.8))+
  stat_summary(fun.data=mean_se,geom="errorbar",fun.args=list(mult=1),
               colour="black",width=0.2,position=position_dodge(0.8))+
  stat_summary(fun=mean,geom="point",shape=21,fill="grey",
               size=3,position=position_dodge(0.8))+
  ylab("Indirect care (number of scans)")+
  xlab("Female Treatment")+
  theme_bw()

plot_mod5_f

#plot results (male presence)
plot_mod5_m<- ggplot(subset(beetles, !is.na(f_indirect_care)),aes(x=male_treatment,y=f_indirect_care))+
  geom_dotplot(binaxis="y",stackdir="center",alpha=0.6,dotsize=0.5,
               position=position_dodge(0.8))+
  stat_summary(fun.data=mean_se,geom="errorbar",fun.args=list(mult=1),
               colour="black",width=0.2,position=position_dodge(0.8))+
  stat_summary(fun=mean,geom="point",shape=21,fill="grey",
               size=3,position=position_dodge(0.8))+
  ylab("Female indirect care (number of scans)")+
  xlab("Male Presence")+
  theme_bw()

plot_mod5_m  


#----- MODEL 6 - MALE DIRECT CARE -----

#came across errors but cant look at effect of male presence on male care (no males)
#so -> create subset df for only when male was included
beetleMales<- beetles[beetles$male_treatment == 1,]
View(beetleMales)

hist(beetleMales$m_direct_care)

mod6.1<- glmmTMB(cbind(m_direct_care,mfreq_no_direct)~female_treatment,
                 data=beetleMales,family="binomial")

summary(mod6.1)

plot(simulateResiduals(mod6.1)) #does not seem overdispersed

#try overdispersal model to see if it works better?
mod6.2<- glmmTMB(cbind(m_direct_care,mfreq_no_direct)~female_treatment+
                   (1|brood_id),data=beetleMales,family="binomial")

plot(simulateResiduals(mod6.2)) #unsure which is better

#how does it fit the data (mod6.1)
sim_max<- apply(simulate(mod6.1,nsim=1000),2,max)
hist(sim_max,breaks=100)
abline(v=max(beetleMales$m_direct_care,na.rm=T),col="red",lwd=2)
#?????

#how well does it fit the data (mod6.2)
sim_max<- apply(simulate(mod6.2,nsim=1000),2,max)
hist(sim_max,breaks=100)
abline(v=max(beetleMales$m_direct_care,na.rm=T),col="red",lwd=2)

#how well does model predict zeros (mod6.1)
simy<- simulate(mod6.1,1000)
nz<- c()
for(i in seq(1,length(simy),1)){nz[i]<- colSums(simy[,i]==0)[1]}
hist(nz,main="Number of zeros")
abline(v=sum(beetleMales$m_direct_care==0,na.rm=T),col="red",lwd=2)
#seems okay 

#how well does model predict zeros (mod6.2)
simy<- simulate(mod6.2,1000)
nz<- c()
for(i in seq(1,length(simy),1)){nz[i]<- colSums(simy[,i]==0)[1]}
hist(nz,main="Number of zeros")
abline(v=sum(beetleMales$m_direct_care==0,na.rm=T),col="red",lwd=2)
#also seems okay

summary(mod6.1)
summary(mod6.2) 

#try including zero inflation to see if it does better (no overdispersal)
mod6.3<- glmmTMB(cbind(m_direct_care,mfreq_no_direct)~female_treatment,ziformula=~1,
                 data=beetleMales,family="binomial")

plot(simulateResiduals(mod6.3)) 

#how does it fit the data
sim_max<- apply(simulate(mod6.3,nsim=1000),2,max)
hist(sim_max,breaks=100)
abline(v=max(beetleMales$m_direct_care,na.rm=T),col="red",lwd=2)

#how well does model predict zeros
simy<- simulate(mod6.3,1000)
nz<- c()
for(i in seq(1,length(simy),1)){nz[i]<- colSums(simy[,i]==0)[1]}
hist(nz,main="Number of zeros")
abline(v=sum(beetleMales$m_direct_care==0,na.rm=T),col="red",lwd=2)
#better

#try including zero inflation to see if it does better (including overdispersal)
mod6.4<- glmmTMB(cbind(m_direct_care,mfreq_no_direct)~female_treatment+
                   (1|brood_id),ziformula=~1,
                 data=beetleMales,family="binomial")

plot(simulateResiduals(mod6.4)) 

#how does it fit the data
sim_max<- apply(simulate(mod6.4,nsim=1000),2,max)
hist(sim_max,breaks=100)
abline(v=max(beetleMales$m_direct_care,na.rm=T),col="red",lwd=2)

#how well does model predict zeros
simy<- simulate(mod6.4,1000)
nz<- c()
for(i in seq(1,length(simy),1)){nz[i]<- colSums(simy[,i]==0)[1]}
hist(nz,main="Number of zeros")
abline(v=sum(beetleMales$m_direct_care==0,na.rm=T),col="red",lwd=2)
#better

#how to choose whether (1|brood_id) should be included or not 
#(use simplest mode? so don't include)

# USE THIS MODEL?
mod6.3<- glmmTMB(cbind(m_direct_care,mfreq_no_direct)~female_treatment,
                 ziformula=~1,data=beetleMales,family="binomial")

#model outputs
Anova(mod6.3)
summary(mod6.3)

#pairwise comparisons
summary(glht(mod6.3,mcp(female_treatment="Tukey")),test=adjusted("bonferroni"))
#no effect of female treatment on male direct care?

#plot results
plot_mod6<- ggplot(subset(beetleMales, !is.na(m_direct_care)),aes(x=female_treatment,y=m_direct_care))+
  geom_dotplot(binaxis="y",stackdir="center",alpha=0.6,dotsize=0.5,
               position=position_dodge(0.8))+
  stat_summary(fun.data=mean_se,geom="errorbar",fun.args=list(mult=1),
               colour="black",width=0.2,position=position_dodge(0.8))+
  stat_summary(fun=mean,geom="point",shape=21,fill="grey",
               size=3,position=position_dodge(0.8))+
  ylab("Direct care (number of scans)")+
  xlab("Female Treatment")+
  theme_bw()

plot_mod6     
      

#----- MODEL 7 - MALE INDIRECT CARE -----
mod7.1<- glmmTMB(cbind(m_indirect_care,mfreq_no_indirect)~female_treatment,
                 data=beetleMales,family="binomial")

summary(mod7.1)

plot(simulateResiduals(mod7.1))  #overdispersed?

mod7.2<- glmmTMB(cbind(m_indirect_care,mfreq_no_indirect)~female_treatment+
                   (1|brood_id),data=beetleMales,family="binomial")

plot(simulateResiduals(mod7.2))  #looks much better

#chek how well it fits the data
sim_max<- apply(simulate(mod7.2,nsim=1000),2,max)
hist(sim_max,breaks=100)
abline(v=max(beetleMales$m_indirect_care,na.rm=T),col="red",lwd=2)
#this seems bad!!!!!

#how well does model predict zeros
simy<- simulate(mod7.2,1000)
nz<- c()
for(i in seq(1,length(simy),1)){nz[i]<- colSums(simy[,i]==0)[1]}
hist(nz,main="Number of zeros")
abline(v=sum(beetleMales$m_indirect_care==0,na.rm=T),col="red",lwd=2)
#seems good

#try zi
mod7.3<- glmmTMB(cbind(m_indirect_care,mfreq_no_indirect)~female_treatment+
                   (1|brood_id),ziformula=~1,data=beetleMales,family="binomial")

plot(simulateResiduals(mod7.3))

#chek how well it fits the data
sim_max<- apply(simulate(mod7.3,nsim=1000),2,max)
hist(sim_max,breaks=10)
abline(v=max(beetleMales$m_indirect_care,na.rm=T),col="red",lwd=2)
#this seems better

#how well does model predict zeros
simy<- simulate(mod7.3,1000)
nz<- c()
for(i in seq(1,length(simy),1)){nz[i]<- colSums(simy[,i]==0)[1]}
hist(nz,main="Number of zeros")
abline(v=sum(beetleMales$m_indirect_care==0,na.rm=T),col="red",lwd=2)

#NOTE: this model was not zero inflated but including ziformula ...
#allowed it to fit the data much better 

#model outputs
Anova(mod7.3)
summary(mod7.3)

#pairwise comparisons
summary(glht(mod7.3,mcp(female_treatment="Tukey")),test=adjusted("bonferroni"))
#no effect of female treatment on male indirect care?

#plot results
plot_mod7<- ggplot(subset(beetleMales, !is.na(m_indirect_care)),aes(x=female_treatment,y=m_indirect_care))+
  geom_dotplot(binaxis="y",stackdir="center",alpha=0.6,dotsize=0.5,
               position=position_dodge(0.8))+
  stat_summary(fun.data=mean_se,geom="errorbar",fun.args=list(mult=1),
               colour="black",width=0.2,position=position_dodge(0.8))+
  stat_summary(fun=mean,geom="point",shape=21,fill="grey",
               size=3,position=position_dodge(0.8))+
  ylab("Inirect care (number of scans)")+
  xlab("Female Treatment")+
  theme_bw()

plot_mod7  


# ----- EFFECT ON LARVAL PERFORMANCE -----

# DISPERSAL LARVAE NUMBER
hist(beetles$disp_larvae_no) #poisson distributed

#initial model
mod8.1<- glm(disp_larvae_no~female_treatment+male_treatment+
             female_treatment*male_treatment,data=beetles,family="poisson")

summary(mod8.1)

#check overdispersal
plot(simulateResiduals(mod8.1)) #not overdispersed

#check how well it fits the data
sim_max<- apply(simulate(mod8.1,nsim=1000),2,max)
hist(sim_max,breaks=10)
abline(v=max(beetles$disp_larvae_no,na.rm=T),col="red",lwd=2)
#within predicted values

Anova(mod8.1)
summary(mod8.1)

# DISPERSAL LARVAE WEIGHT (BROOD)
hist(beetles$disp_larvae_weight) #poisson distributed

#initial model
mod9.1<- glm(disp_larvae_weight~female_treatment+male_treatment+
               female_treatment*male_treatment,data=beetles,family="poisson")

summary(mod8.1)

#check overdispersal
plot(simulateResiduals(mod8.1)) #not overdispersed

#check how well it fits the data
sim_max<- apply(simulate(mod8.1,nsim=1000),2,max)
hist(sim_max,breaks=10)
abline(v=max(beetles$disp_larvae_no,na.rm=T),col="red",lwd=2)
#within predicted values

Anova(mod8.1)
summary(mod8.1)


# dispersal larvae weight (average individual)
