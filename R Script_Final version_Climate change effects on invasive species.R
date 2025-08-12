#Author: Nicola Smith
#Date: 16 July 2025
#Project: Climate change effects on marine invasive species global meta-analysis
#Version: Final


#### Step 1 - Retrieve and Clean Data ####
getwd()
setwd("/Users/nicolas.smith/Desktop/Climate change and invasive species/Climate change-invasive species interactions")
investigate<-read.csv("Coded data.csv")
investigate
head(investigate)

str(investigate)

#load packages
install.packages("tidyverse")
install.packages("metafor")
library(tidyverse)
library(metafor)
# Check R version
R.Version()$version.string


#Get R documentation for metafor package
help(metafor)


#Remove Stachowicz study b/c it uses the correlation coefficient and we are
  #dealing only with means for the LnRR

investigate<-(subset(investigate, authors!="Stachowicz et al."))
investigate
nrow(investigate)
str(investigate)



#Convert means into LnRR for investigate data 
outcomes<-escalc(measure="ROM",m1i=after_mean, m2i=before_mean, n1i=dummy_n, n2i=dummy_n,
                 sd1i=sd_after, sd2i=sd_before,data=investigate)

outcomes %>%
  View


# Remove studies that estimate more than 1 time frame so there are no duplicate studies

# I need to remove studies that give results for more than one time frame (e.g., models that 
# predict an effect 50 vs 100 years later, use only one, i.e. 50 years)
outcomes

no.duplicates<-subset(outcomes, duration !="100 years")

no.duplicates

no.duplicates %>%
  View

complete_data_no_duplicates<-subset(no.duplicates, duration !="140 years (1960 to 2100)")

complete_data_no_duplicates %>%
  View

nrow(complete_data_no_duplicates)
str(complete_data_no_duplicates)

# There are n = 69-2 = 67 studies in the complete_data_no_duplicates dataset b/c two of the studies are NA b/c effect size = Infinity

#### Step 2: Run multi-level meta-regression models ####

#Note: I shouldn't run a random effects model because the exps, obs, and models are done on 
#completely different time scales so the overall effect size for all studies is not meaningful

# ...........................................

# Mod 1: study type as a fixed effect and study/author as a random effect

mixed_no_duplicates_study.type<- rma.mv(yi,vi, mods = ~ study.type - 1, random = ~1 | authors/code, 
                                        method = "REML", digits = 3, data = complete_data_no_duplicates)
mixed_no_duplicates_study.type

# _____ Back-transformation LnRR to risk ratios and then to percentages for model outputs

experiments<-exp(0.562) -1
experiments*100
#75.41773

models<-exp(3.984) -1
models*100
#5273.153

observational<-exp(4.159) - 1
observational*100
#6300.748


#Rosenthal Fail Safe N
fsn(yi, vi, data=complete_data_no_duplicates)

#Fail-safe N Calculation Using the Rosenthal Approach

#Observed Significance Level: <.0001
#Target Significance Level:   0.05

#Fail-safe N: 47694

# .............................................

# Mod 2: Mixed-effects model for exps only where study is a random effect and functional group is a fixed effect

# First, collect only exp data

exp_only<-subset(complete_data_no_duplicates, study.type=="experiment")
exp_only
str(exp_only)

exp_only %>%
  View

# However, there are only n=2 studies for primary produces in exp, so we should probably remove it b/c of the small n
#and then run the mixed effects model again
exp_only_no_seaweed<-subset(exp_only,taxon!="seaweed")
exp_only_no_seaweed

exp_only_no_seaweed %>%
    View

mix_effect_functional.grp_exp_no_seaweed<-rma.mv(yi,vi, mods = ~ functional.group - 1, random = ~1 | authors/code,
                                         method = "REML", digits = 4, data = exp_only_no_seaweed)
mix_effect_functional.grp_exp_no_seaweed


# There are n=14 effect sizes that are just from experiments without primary producers 
  #but n = 12 b/c 2 NAs in the dataframe (i.e., 2 studies from Sorte et al. 2010)


#................................................

# Mod 3: # Mixed-effects model for models only where study is a random effect and functional group is a fixed effect

# First, collect only model data

models_only<-subset(complete_data_no_duplicates, study.type=="model")
models_only

#Remove omnivore group b/c n=1
models_only_no_omnivore<-subset(models_only,functional.group!="omniovre")
models_only_no_omnivore

# Remove primary producers because n=2
models_only_no_omnivore_no_pp<-subset(models_only_no_omnivore,functional.group!="primary producer")
models_only_no_omnivore_no_pp


mix_effect_functional.grp_models<-rma.mv(yi,vi, mods = ~ functional.group - 1, random = ~1 | authors/code,
                                         method = "REML", digits = 4, data = models_only_no_omnivore_no_pp)

mix_effect_functional.grp_models

# .....................................................

# Model 4:Mixed-effects model for obs only where study is a random effect and functional group is a fixed effect

# First, collect only obs data

obs_only<-subset(complete_data_no_duplicates, study.type=="observational")
obs_only
str(obs_only)

mix_effect_functional.grp_obs<-rma.mv(yi,vi, mods = ~ functional.group - 1, random = ~1 | authors/code,
                                      method = "REML", digits = 4, data = obs_only)

mix_effect_functional.grp_obs

# ......................................................

# Model 5: Mixed-effects model for exps only where study is a random effect and invasion stage is a fixed effect

# Remove introduction stage from model analysis b/c n=2, which is too small

intro<-subset(exp_only, invasion.stage=="introduction")
intro
nrow(intro)
#There are n=2 exps where invasion stage is introduction

mod1<-subset(exp_only,invasion.stage !="introduction")
mod1
nrow(mod1)
str(mod1)
View(mod1)

# Remove spread stage from model analysis b/c n=1, which is too small
mod1<-subset(mod1,invasion.stage!="spread")
mod1
View(mod1)

#I've used an optimizer in the below model because the default did not converge

mix_effect_invasion.stage_exp<-rma.mv(yi,vi, mods = ~ invasion.stage - 1, random = ~1 | authors/code,
                                      method = "REML", digits = 4, data = mod1, 
                                      verbose=TRUE, control=list(rel.tol=1e-8))
mix_effect_invasion.stage_exp

# Backtransform LnRR to risk ratio and then to % for significant model output

establishment<-exp(0.6905)-1
establishment*100

# .....................................................

# Model 6: Mixed-effects model for models only where study is a random effect and invasion stage is a fixed effect

# Remove establishment stage from model analysis b/c n=2, which is too small
view(models_only)
est<-subset(models_only, invasion.stage=="establishment")
est
nrow(est)
#There are only  n=2 exps where invasion stage is establishment

mod4<-subset(models_only,invasion.stage!="establishment")

# I cannot run mixed models for invasion stage in models b/c all models only looked at the spread stage of invasion

# ......................................................

# Model 7: Mixed-effects model for obs only where study is a random effect and invasion stage is a fixed effect

# Remove introduction stage from the model b/c n=2 for intros
view(obs_only)
mod3<-subset(obs_only,invasion.stage !="introduction")
view(mod3)

mix_effect_invasion.stage_obs<-rma.mv(yi,vi, mods = ~ invasion.stage - 1, random = ~1 | authors/code,
                                      method = "REML", digits = 4, data = mod3)

mix_effect_invasion.stage_obs


# Model 8: Mixed-effects model for exp only where study is a random effect and invader response is a fixed effect

# Remove species richness b/c n=2, which is too small

# n = 4 for abundance
#n = 4 for impact
# n = 4 for performance
#n = 2 species richness

mod7<-subset(exp_only, biological.response !="invasive species richness")
mod7
view(mod7)


mix_effect_invader.response_exp<-rma.mv(yi,vi, mods = ~ biological.response - 1, random = ~1 | authors/code,
                                        method = "REML", digits = 4, data = mod7)

mix_effect_invader.response_exp

# Backtransform LnRR to risk ratio then to % for significant model output

abundance<-exp(1.1175)-1
abundance*100

# Model 9: Mixed-effects model for models only where study is a random effect and invader response is a fixed effect

#n = 2 for fecundity so remove it from the model

mod5<-subset(models_only, biological.response != "invader fecundity")
mod5

View(mod5)


mix_effect_invader.response_models<-rma.mv(yi,vi, mods = ~ biological.response - 1, random = ~1 | authors/code,
                                           method = "REML", digits = 4, data = mod5)

mix_effect_invader.response_models


# Model 10: Mixed-effects model for obs only where study is a random effect and invader response is a fixed effect

mix_effect_invader.response_obs<-rma.mv(yi,vi, mods = ~ biological.response - 1, random = ~1 | authors/code,
                                        method = "REML", digits = 4, data = obs_only)

mix_effect_invader.response_obs


# Backtransform LnRR to risk ratio then to % for significant model output

abundance<-exp(6.0949)-1
abundance*100

# Model 11: Mixed-effects model for exps where study is a random effect and stressor type is a fixed effect

# Remove seawater freshening and acidification b/c in both cases, n = 2

stressor_exp<-subset(exp_only, stressor.type != "seawater freshening")
stressor_exp

stressor_exp2<-subset(stressor_exp, stressor.type != "acidification")
stressor_exp2

view(stressor_exp2)

mix_effect_stressor_exp <-rma.mv(yi,vi, mods = ~ stressor.type - 1, random = ~1 | authors/code,
                                 method = "REML", digits = 4, data = stressor_exp2)

mix_effect_stressor_exp 

# Backtransform LnRR to risk ratio then to % for significant model output

extreme<-exp(0.6876)-1
extreme*100



# Model 12: Mixed-effects model for obs where study is a random effect and stressor type is a fixed effect

mix_effect_stressor_obs <-rma.mv(yi,vi, mods = ~ stressor.type - 1, random = ~1 | authors/code,
                                 method = "REML", digits = 4, data = obs_only)

mix_effect_stressor_obs 

# Backtransform LnRR to risk ratio then to % for significant model output

warming<-exp( 6.0949)-1
warming*100

#### ..................................... Figures .................................####

# Fig.2 = LnRR by study type (i.e., Mod1 mixed effects model)

#First make a dataframe of results in excel and call it results
results<-read.csv("results1.csv") 
results

A<-ggplot(data=results, aes(x=Study.type, y=LnRR,  shape=Temporal.scale, size=Temporal.scale, color = Temporal.scale)) +
  geom_point() +
  geom_pointrange(aes(ymax=ub, ymin=lb)) +
  geom_hline(yintercept=0, linetype="dashed") +
  scale_shape_manual(values=c(16, 16,1 )) +
  scale_color_manual(values=c('orange','blue','purple'))+
  scale_size_manual(values=c(0.85,0.85,0.85))+
  coord_flip() +
  theme_bw() +
  theme(legend.text=element_text(size=15), legend.title=element_text(size=17)) +
  theme(axis.title.y = element_blank()) 
A

A1<-A + theme(axis.text = element_text(size = 20)) 
A1          

A2<- A1 + theme(axis.title.x = element_text(vjust = -0.8))
A2

A3<- A2 + theme(axis.title = element_text(size = 23))
A3

# Fig 3. LnRR for all hypotheses/moderators

# First create fig for mixed model for exp where invasion stage is a fixed effect

results3<-read.csv("results3.csv")
results3

B<-ggplot(data=results3, aes(x=Invasion.stage, y=LnRR, shape=Invasion.stage, color=Temporal.scale, size=Temporal.scale)) +
  geom_point() +
  geom_pointrange(aes(ymax=ub, ymin=lb)) +
  geom_hline(yintercept=0, linetype="dashed") +
  scale_shape_manual(values=c(1,1)) +
  scale_color_manual(values = c('purple')) +
  scale_size_manual(values=c(0.85,0.85))+
  scale_y_continuous(name="LnRR", limits=c(-6,15) )+
  coord_flip() +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.title.y = element_blank()) +
  theme(axis.title.x = element_blank())
B

B3<- B + theme(text = element_text(size = 20))
B3

# Fig for mixed model for obs where invasion stage is a fixed effect

results_obs_stage<-read.csv("results_obs_stage.csv")
results_obs_stage

C<-ggplot(data=results_obs_stage, aes(x=Invasion.stage, y=LnRR, color=Temporal.scale, shape=Invasion.stage, size=Temporal.scale)) +
  geom_point() +
  geom_pointrange(aes(ymax=ub, ymin=lb)) +
  geom_hline(yintercept=0, linetype="dashed") +
  scale_shape_manual(values=c(1,1)) +
  scale_color_manual(values = c('blue')) +
  scale_size_manual(values=c(0.85,0.85))+
  scale_y_continuous(name="LnRR", limits=c(-6,15) )+
  coord_flip() +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.title.y = element_blank()) +
  theme(axis.title.x = element_blank())
C

C3<-C + theme(text = element_text(size = 20))
C3


# Fig for mixed model for models where invasion stage is a fixed effect cannot be done
#because all models examined the spread stage except for two that examined establishment



# Fig for mixed model for exps where invasion measure is a fixed effect

results4<-read.csv("results4.csv")
results4

E<-ggplot(data=results4, aes(x=Invader.response, y=LnRR, color=code, shape=Invader.response, size=code)) +
  geom_point() +
  geom_pointrange(aes(ymax=ub, ymin=lb)) +
  geom_hline(yintercept=0, linetype="dashed") +
  scale_shape_manual(values=c(16,1,1)) +
  scale_color_manual(values = c('purple')) +
  scale_size_manual(values=c(0.85,0.85,0.85))+
  scale_y_continuous(name="LnRR", limits=c(-6,15) )+
  coord_flip() +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.title.y = element_blank()) +
  theme(axis.title.x = element_blank())
E

E3<- E + theme(text = element_text(size = 20))
E3


# Fig for mixed effects model for obs where invasion measure is a fixed effect
results5<-read.csv("results5.csv")
results5

F<-ggplot(data=results5, aes(x=Invader.response, y=LnRR, color=code, shape=Invader.response, size=code)) +
  geom_point() +
  geom_pointrange(aes(ymax=ub, ymin=lb)) +
  geom_hline(yintercept=0, linetype="dashed") +
  scale_shape_manual(values=c(16,1)) +
  scale_color_manual(values = c('blue')) +
  scale_size_manual(values=c(0.85,0.85))+
  scale_y_continuous(name="LnRR", limits=c(-6,15) )+
  coord_flip() +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.title.y = element_blank()) +
  theme(axis.title.x = element_blank())
F

F3<- F + theme(text = element_text(size = 20))
F3

# Fig for mixed effects model for models where invasion measure is a fixed effect

response_mods_only<-read.csv("response_mods_only.csv")
response_mods_only

G<-ggplot(data=response_mods_only, aes(x=Response, y=LnRR, color=code, shape=Response, size=code)) +
  geom_point() +
  geom_pointrange(aes(ymax=ub, ymin=lb)) +
  geom_hline(yintercept=0, linetype="dashed") +
  scale_shape_manual(values=c(1,1)) +
  scale_color_manual(values = c('orange')) +
  scale_size_manual(values=c(0.85,0.85))+
  scale_y_continuous(name="LnRR", limits=c(-6,15) )+
  coord_flip() +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.title.y = element_blank()) +
  theme(axis.title.x = element_blank())
G

G3<- G + theme(text = element_text(size = 20))
G3

# Fig for mixed model for exps where functional group is a fixed effect

results_exp_functional_group<-read.csv("results_exp_functional_group.csv")
results_exp_functional_group

H<-ggplot(data=results_exp_functional_group, aes(x=Functional.group, y=LnRR, color=code, shape=Functional.group, size=code)) +
  geom_point() +
  geom_pointrange(aes(ymax=ub, ymin=lb)) +
  geom_hline(yintercept=0, linetype="dashed") +
  scale_shape_manual(values=c(1,1,1)) +
  scale_color_manual(values = c('purple')) +
  scale_size_manual(values=c(0.85,0.85,0.85))+
  scale_y_continuous(name="LnRR", limits=c(-6,15) )+
  coord_flip() +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.title.y = element_blank()) +
  theme(axis.title.x = element_blank())
H

H3<- H + theme(text = element_text(size = 20))
H3


# Fig for mixed model for obs where functional group is a fixed effect
results_functional_group_obs<-read.csv("results_functional_group_obs.csv")
results_functional_group_obs

W<-ggplot(data=results_functional_group_obs, aes(x=Functional_group, y=LnRR, color=code, shape=Functional_group, size=code)) +
  geom_point() +
  geom_pointrange(aes(ymax=ub, ymin=lb)) +
  geom_hline(yintercept=0, linetype="dashed") +
  scale_shape_manual(values=c(1,1)) +
  scale_color_manual(values = c('blue')) +
  scale_size_manual(values=c(0.85,0.85))+
  scale_y_continuous(name="LnRR", limits=c(-6,15) )+
  coord_flip() +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.title.y = element_blank()) +
  theme(axis.title.x = element_blank())
W

W3<- W + theme(text = element_text(size = 20))
W3

# Fig for mixed model for models where functional group is a fixed effect
results_functional_group_models<-read.csv("results_functional_group_models.csv")
results_functional_group_models


# Fig for mixed model for exp where climate stressor is a fixed effect
Results_Stressor_exp<-read.csv("Results_Stressor_exp.csv")
Results_Stressor_exp

J<-ggplot(data=Results_Stressor_exp, aes(x=Stressor, y=LnRR, color=code, shape=Stressor, size=code)) +
  geom_point() +
  geom_pointrange(aes(ymax=ub, ymin=lb)) +
  geom_hline(yintercept=0, linetype="dashed") +
  scale_shape_manual(values=c(1,1)) +
  scale_color_manual(values = c('purple')) +
  scale_size_manual(values=c(0.85,0.85))+
  scale_y_continuous(name="LnRR", limits=c(-6,15) )+
  coord_flip() +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.title.y = element_blank()) +
  theme(axis.title.x = element_blank())
J

J3<- J + theme(text = element_text(size = 20))
J3

# Fig for mixed model for obs where climate stressor is a fixed effect
Results_Stressor_obs<-read.csv("Results_Stressor_obs.csv")
Results_Stressor_obs

K<-ggplot(data=Results_Stressor_obs, aes(x=Stressor, y=LnRR, color=code, shape=Stressor, size=code)) +
  geom_point() +
  geom_pointrange(aes(ymax=ub, ymin=lb)) +
  geom_hline(yintercept=0, linetype="dashed") +
  scale_shape_manual(values=c(1,16)) +
  scale_color_manual(values = c('blue')) +
  scale_size_manual(values=c(0.85,0.85))+
  scale_y_continuous(name="LnRR", limits=c(-6,15) )+
  coord_flip() +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.title.y = element_blank()) +
  theme(axis.title.x = element_blank())
K

K3<- K + theme(text = element_text(size = 20))
K3


#### Combine all plots on one figure ####

library("ggplot2")
library("ggpubr")
library("gridExtra")

figure<-ggarrange(B3,C3, E3, F3, G3,H3, W3,I3,J3,K3,
                  labels = c("A","B","C","D","E","F","G","H","I","J"),
                  ncol = 3, nrow = 4)
figure


figureA<-ggarrange(B3,C3,
                   labels = c("A","B"),
                   ncol = 2, nrow = 1)
figureA


figureB<-ggarrange( E3, F3, G3,
                  labels = c("C","D","E"),
                  ncol = 3, nrow = 1)
figureB


figureC<-ggarrange( H3, W3,I3,
                    labels = c("F","G","H"),
                    ncol = 3, nrow = 1)
figureC

figureD<-ggarrange(J3,K3,
                   labels = c("I","J"),
                   ncol = 2, nrow = 1)
figureD

?par()





#### Sengrid#### Sensitivity analyses ####

# Analysis 1 - Remove Pacific oyster from the dataset

complete_data_no_duplicates
view(complete_data_no_duplicates)
nrow(complete_data_no_duplicates)

   # Remove the Pacific oyster from the main dataset
no.oyster<-subset(complete_data_no_duplicates,taxon!="oyster")
no.oyster
nrow(no.oyster)
view(no.oyster)
# N = 38 for main dataset with Pacific oyster removed, but remove 2 studies b/c of NAs for LnRR
 # so N = 38-2 = 36

#Run mixed-effects model with study type as a fixed effect and study as a random effect

mixed_no.oyster_study.type<- rma.mv(yi,vi, mods = ~ study.type - 1, random = ~1 | authors/code, 
                                        method = "REML", digits = 3, data = no.oyster)
mixed_no.oyster_study.type

# Plot results from sensitivity analysis
    #First make a dataframe of results in excel and call it results
results.no.oyster<-read.csv("results.no.oyster.csv") 
results.no.oyster

A<-ggplot(data=results.no.oyster, aes(x=Study.type, y=LnRR,  shape=Temporal.scale, size=Temporal.scale, color = Temporal.scale)) +
  geom_point() +
  geom_pointrange(aes(ymax=ub, ymin=lb)) +
  geom_hline(yintercept=0, linetype="dashed") +
  scale_shape_manual(values=c(16, 16,1 )) +
  scale_color_manual(values=c('orange','blue','purple'))+
  scale_size_manual(values=c(0.85,0.85,0.85))+
  coord_flip() +
  theme_bw() +
  theme(legend.text=element_text(size=15), legend.title=element_text(size=17)) +
  theme(axis.title.y = element_blank()) 
A

A1<-A + theme(axis.text = element_text(size = 20)) 
A1          

A2<- A1 + theme(axis.title.x = element_text(vjust = -0.8))
A2

A3<- A2 + theme(axis.title = element_text(size = 23))
A3


#What is the sample size by study type?
exp.no.oyster<-subset(no.oyster,study.type=="experiment")
exp.no.oyster
nrow(exp.no.oyster)
# N = 16 but 2 NAs for LnRR, so N = 16-2=14 for exps
model.no.oyster<-subset(no.oyster,study.type=="model")
model.no.oyster
nrow(model.no.oyster)
#N=17 for models
obs.no.oyster<-subset(no.oyster,study.type=="observational")
obs.no.oyster
nrow(obs.no.oyster)
#N=5 for observational studies

