####################################
####################################
#
# Key questions to address:
#
# 1. explore the relationship between CFQ and each of the set size effect, cuing cost, cuing benefit and flanker interference effects.
#
# 2. explore the relationship between DAT1 and CFQ/ADHD scores and the cuing data (probably nothing)
#    and competition data (without regard to target position).
#
####################################
####################################
#
#####  IF RUNNING THIS SCRIPT FOR THE FIRST TIME:
#####  Install some handy packages:
#####  Remove the "# "in front of the line below and run the code. Replace the # after installing the packages, otherwise the R 
#  will reinstall the packates every time you run the script
#  install.packages(c("MASS", "robustbase", "robust", "mgcv", "scatterplot3d", "quantreg", "rrcov", "lars", "pwr", "mc2d", "psych", 
#                     "Rfit","MBESS", "BayesFactor", "PoweR", "ggplot2", "reshape2", "plyr",
#                      "rmarkdown", "car", "gridExtra", "bootES", "BEST","foreign","nlme","pastecs","multcomp","ggplot2","compute.es"
#                     ,"ez","lattice","lme4","effects","diagram","png", "grid"))
## load relevant packages 
library(foreign)
library(car)
library(ggplot2)
library(pastecs)
library(psych)
library(plyr)
library(multcomp)
library(reshape2)
library(compute.es)
library(ez)
library(lattice)
library(lme4)
library(png)
library(grid)
#
##############################################################################
##############################################################################
####                                                                      ####
####        ~~~~~~         Competition data           ~~~~~~~             ####
####                                                                      ####
##############################################################################
##############################################################################
#
##### 1. Import single trial Competition data
compdata <- read.csv("S:/R-MNHS-SPP/Bellgrove-data/CogBattery_UQ&Monash/AllParticipantData_comp_CFQCAARS.csv", header=T, na.strings="#N/A")
#na.strings="#N/A" converts all #N/As in the document to NA.

# 2. Competition Data: Look at Column names:
colnames(compdata)
#Rename some of the columns:
compdata<-rename(compdata, c("ResponseTime"="RT", "ResponseStatus"="Accuracy"))

# 3. Remove any row that has a NAs in RT data:
compdata<-compdata[complete.cases(compdata$RT),] 

# 4. Check number of Trials for each participant by running the function 'length', 
#Have a look at this
#Competition: everybody has 320 trials to start with
num_trials1 <- ddply(compdata, c("ID"), summarise,
                     Trials    = length(RT))
summary(num_trials1$Trials)

# 5.Calculate each participant's accuracy as a %
Comp_accuracy_checker <- ddply(compdata, c("ID"), summarise,
                          Hits  = sum(Accuracy=="Correct"),
                          Misses = sum(Accuracy=="Incorrect"))
Comp_accuracy_checker$Total=Comp_accuracy_checker$Hits+Comp_accuracy_checker$Misses
Comp_accuracy_checker$Accuracy=(Comp_accuracy_checker$Hits/Comp_accuracy_checker$Total)*100
summary(Comp_accuracy_checker$Accuracy)
#Plot the accuracy distribution 
hist(Comp_accuracy_checker$Accuracy)

# 6. Remove any one who is guessing
#
#This task has 50% chance of guessing each trial correctly
#a binom.test shows that participants must have at least 179 correct out of
#the 320 trials in order for us to be 95% confident that they were not guessing
#the answer
binom.test(x=179, n=320, p=0.5)
#Most participant's have greater than 178 trials correct, check who doesn't:
ThoseGuessing<-Comp_accuracy_checker[Comp_accuracy_checker$Hits<178,]
ThoseGuessing
#make a new logical column in the main dataframe for whether each
#participant was guessing or not (TRUE/FALSE):
for (i in 1:length(compdata$ID)) {compdata$guessing[i]<-any(ThoseGuessing$ID==compdata$ID[i])}


# 7. Make the required factors:
compdata$ID <- factor(compdata$ID)
compdata$SetSize <- factor(compdata$SetSize)
compdata$TargetHem <- factor(compdata$TargetHem)
compdata$TargetLoc <- factor(compdata$TargetLoc)
compdata$DisplayType <- factor(compdata$DisplayType)
compdata$TargetType <- factor(compdata$TargetType)
compdata$Trial <- factor(compdata$Trial)
compdata$Site <- factor(compdata$Site)
compdata$dat1utr <- factor(compdata$dat1utr)
compdata$dat1i8 <- factor(compdata$dat1i8)
compdata$dat1haplo <- factor(compdata$dat1haplo)
compdata$drd4 <- factor(compdata$drd4)

#Rename factor Levels:
compdata$TargetHem <- revalue(compdata$TargetHem, c("1"="Left", "2"="Right")) #Hi Beth, please just check this is correct i.e. that "1" does indeed stand for "Left" here
compdata$DisplayType <- revalue(compdata$DisplayType, c("1"="Unilateral", "2"="Bilateral")) #Beth also better double check that "1" was actually "Unilateral
compdata$SetSize <- revalue(compdata$SetSize , c("4"="Four", "8"="Eight"))
compdata$dat1utr<- revalue(compdata$dat1utr, c("0"="zero","1"="one","2"="two")) 

# 8. Tidy data
#
# Look at RT distribution for each DAT1 group
ggplot(compdata, aes(RT))  + geom_histogram(aes(y=..count..), colour="black", fill="white") + facet_wrap(~ dat1utr)

#Now kick out those participants who were guessing 
compdata<-compdata[compdata$guessing!=TRUE,] 

#Kick out trials with stim. timing problems
compdata<-compdata[compdata$Timing.Status=="OK",]

#Remove incorrect trials
compdata<-compdata[compdata$Accuracy=="Correct",]

#Remove trials where RT too fast 
compdata<-compdata[!compdata$RT<200,]
#....or too slow to be real RTs
compdata<-compdata[!compdata$RT>2000,]
#Now have a look at the RT distribution:
hist(compdata$RT)

compdata$log_RT<-log(compdata$RT) #log
hist(compdata$log_RT)

##### Z-score each participant's log(RT) data inside the task Conditions ####
compdata$IDbySetSizebyDisplayType<-interaction(compdata$ID, compdata$SetSize, compdata$DisplayType) #ID by SetSize by DisplayType
#calculate mean and sd 
m <- tapply(compdata$log_RT,compdata$IDbySetSizebyDisplayType,mean)
s <- sqrt(tapply(compdata$log_RT,compdata$IDbySetSizebyDisplayType,var))
#calculate RT.Z and save it inside data.frame
compdata$RT.Z <- (compdata$log_RT-m[compdata$IDbySetSizebyDisplayType])/s[compdata$IDbySetSizebyDisplayType]
#check that Z scores have mean=0 and std=1 
RT.Z_checker <- ddply(compdata, c("ID", "SetSize", "DisplayType"), summarise,
                      N    = length(RT.Z ),
                      mean = round(mean(RT.Z )),
                      sd   = sd(RT.Z ),
                      se   = sd / sqrt(N))
#Check to make sure the z-transformation worked (i.e. should see Z mean=0 and SD=1)
summary(RT.Z_checker$mean)
summary(RT.Z_checker$sd)


##Remove trials where absolute RT.Z>3 (i.e. remove outlier RTs)
compdata<-compdata[!abs(compdata$RT.Z)>3,]

#Plot log(RT) historgam again after outlier removal 
hist(compdata$log_RT)

#Check number of Trials for each participant by running the function 'length', 
#on "compdata$RT" for each group, broken down by ID
num_trials2 <- ddply(compdata, c("ID"), summarise,
                     Trials    = length(RT))

summary(num_trials2$Trials)

#Plot the distribution of each participant's number of trials left after kicking out Incorrect and too fast/slow trials:
hist(num_trials2$Trials)

## Kick out Subjs based on CFQ and CAARS scores 
## (need to actually import these scores, but for now just kick out those who 
##  I know have bad CFQ and CAARS scores from previous J.Neuro paper)
toBeRemoved<-which(compdata$Subj=="subject240" | compdata$Subj=="subject341" | compdata$Subj=="subject358" | compdata$Subj=="subject241")

compdata<-compdata[-toBeRemoved,]

compdata<-compdata[!compdata$CAARS_H<1,]

#Calculate Distracter Hemifield (DistHem) as used in the 
#Newman et al. (2014) JoN paper, rather then using TargetHem:
compdata$TargetHembyDisplayType<-interaction(compdata$TargetHem, compdata$DisplayType)
compdata$DistHem <- revalue(compdata$TargetHembyDisplayType, c("Left.Unilateral"="Left", "Left.Bilateral"="Right", "Right.Unilateral"="Right", "Right.Bilateral"="Left"))


####################################
###  Crossed vs Nested Factors:  ###
####################################


#Test if ID is nested within Site:
with(compdata, isNested(ID,Site)) #SO yes, ID is nested within Site...
#Print a table of this nesting:
xtabs(~ ID + Site, compdata, drop = TRUE, sparse = TRUE)
# Because each Participant ("ID")" occurs within one and only one level of 
# Site we say that ID is nested within Site

# The task factors (SetSize, TargetHem, DisplayType, etc) are crossed or 
# partially crossed factors, "crossed" means that they have an observation for 
# each combination level each of the other factors (some factors might be partially 
# crossed in some participants if there is a lot of data missing due to bad accuracy)

# For Multilevel modeling with the lme4 package, nested factors are represended 
# like this "(1 | Site/ID) - which means ID is nested inside 
# Site", while crossed and partially cross factors are represented like 
# this (1 |SetSize) +(1|TargetHem) +(1|DisplayType), etc.

####################################

#Remove the DAT1 rows that are NA:
compdata2<-compdata[complete.cases(compdata$dat1utr),] 
#Remove the CFQ rows that are NA:
compdata2<-compdata2[complete.cases(compdata2$CFQ),] 
#Remove the CAARS (total ADHD) rows that are NA:
compdata2<-compdata2[complete.cases(compdata2$CAARS_H),] 

summary(compdata2$Site) #Summarize numbers of observations still left after NAs have been excluded from the dataset

############################################################################################################
#######    Aggregate to mean instead of trial-by-trial and check for partcicipant-level outliers     #######
############################################################################################################
#  This will alow us to calculate set-size effect which can not be measured on the single trial level      #
#
DF_collapsed<-ddply(compdata2, .(ID, Site, dat1utr, SetSize,TargetLoc,DisplayType,TargetHem, DistHem, CFQ, CAARS_H), summarise, RT=mean(RT))
ggplot(DF_collapsed, aes(RT))  + geom_histogram(aes(y=..count..), colour="black", fill="white") + facet_wrap(~ dat1utr)

# Collapse to each participant so I can Z-score individual participants' overall log(RT)s:
DF_collapsed_IDs<-ddply(data2, .(ID, dat1utr), summarise, RT=mean(RT))
DF_collapsed_IDs$RT.Z<-scale(DF_collapsed_IDs$RT)
ggplot(DF_collapsed_IDs, aes(RT.Z))  + geom_histogram(aes(y=..count..), colour="black", fill="white")
#Split by DAT1 group
ggplot(DF_collapsed_IDs, aes(RT, fill=dat1utr, colour=dat1utr)) + geom_density(alpha=.05) 

#######    Calculate Set-size effect by DistHem    #######
require(reshape2)
#Bring SetSize Factor up into wide format and then calculate SetSize effect
DF_DistHem_wide <- dcast(DF_collapsed, ID + Site + dat1utr + DisplayType + TargetHem + DistHem + CFQ + CAARS_H ~ SetSize, value.var="RT", fun.aggregate=mean)
#Calculate SetSize effect for each DistHem:
DF_DistHem_wide$SetSizeEffect<-DF_DistHem_wide$Eight - DF_DistHem_wide$Four ########### BETH, this is an example of calculating the derived "Setsize effect" (mean RTs from trials with Eight minus those from from trials with Eight). The other tasks have similar derived measures you can calculate, such as the valid/invalid cueing effect 
summary(DF_DistHem_wide$DistHem)

###Now try model DAT1 x DistHem in the collapsed setsize-effect data
random_intercepts_only_SetSizeDistHem<-lmer(SetSizeEffect ~ 1 + (1 | Site/ID) +(1|DisplayType) + (1|DistHem), data = DF_DistHem_wide, na.action = na.omit, REML=FALSE)
DistHem<-update(random_intercepts_only_SetSizeDistHem, .~. + DistHem)
dat1utr<-update(DistHem, .~. + dat1utr)
DistHem_by_dat1utr<-update(dat1utr, .~. + dat1utr*DistHem)
CFQ<-update(DistHem_by_dat1utr, .~. + CFQ)
CAARS_H<-update(DistHem_by_dat1utr, .~. + CAARS_H)
DistHem_by_CFQ<-update(CAARS_H, .~. + CFQ*DistHem)
DistHem_by_CAARS_H<-update(DistHem_by_CFQ, .~. + CAARS_H*DistHem)
dat1utr_by_CFQ<-update(DistHem_by_CAARS_H, .~. + CFQ*dat1utr)
dat1utr_by_CAARS_H<-update(dat1utr_by_CFQ, .~. + CAARS_H*dat1utr)
DistHem_by_dat1utr_by_CFQ<-update(dat1utr_by_CAARS_H, .~. + CFQ*dat1utr*DistHem)
DistHem_by_dat1utr_by_CAARS_H<-update(DistHem_by_dat1utr_by_CFQ, .~. + CAARS_H*dat1utr*DistHem)

anova(random_intercepts_only_SetSizeDistHem, DistHem, dat1utr, DistHem_by_dat1utr,CFQ,CAARS_H,DistHem_by_CFQ,DistHem_by_CAARS_H,dat1utr_by_CFQ,dat1utr_by_CAARS_H,DistHem_by_dat1utr_by_CFQ,DistHem_by_dat1utr_by_CAARS_H)

# DistHem                                7 29475 29516 -14731    29461 5.5418      1    0.01857 *  
# DistHem_by_dat1utr                    11 29477 29541 -14727    29455 5.4835      2    0.06446 .  
# CAARS_H                               12 29467 29538 -14722    29443 8.6401      0    < 2e-16 ***
# dat1utr_by_CFQ                        17 29468 29567 -14717    29434 8.3232      2    0.01558 *  


############################################################################################################
#######                               Plot significant effects                                       #######
############################################################################################################

require(effects)
# Plot CAARS_H  p< 2e-16 ***
#...with effects package:
eff.CAARS_H_SetSizeEffect <- Effect(c("CAARS_H"), CAARS_H)
plot(eff.CAARS_H_SetSizeEffect, main=NULL, rug=T, ticks.x=NULL)
#now plot the same effect with ggplot: 
ggplot(DF_DistHem_wide, aes(x=DF_DistHem_wide$CAARS_H, y=DF_DistHem_wide$SetSizeEffect)) +
  stat_smooth(level = 0.95,size=1) # - good they both look the same - I just wanted to check to make sure the fitted "effects" plot looked the same as plotting the raw data using ggplot

# Plot dat1utr_by_CFQ p=0.01558 * 
eff.dat1utr_by_CFQ_SetSizeEffect <- Effect(c("CFQ", "dat1utr"), dat1utr_by_CFQ)
plot(eff.dat1utr_by_CFQ_SetSizeEffect, layout=c(3,1), main=NULL, rug=T, multiline =F, alternating=F, ticks.x=NULL) 

# Plot DistHem_by_dat1utr p=0.06446  
eff.DistHem_by_dat1utr_SetSizeEffect <- Effect(c("DistHem", "dat1utr"), DistHem_by_dat1utr)
plot(eff.DistHem_by_dat1utr_SetSizeEffect, layout=c(3,1), main=NULL, rug=F, multiline =F, alternating=F, ticks.x=NULL) #(Figure 2A of Newman et al. (2014) JoN paper, the trend still there)
#Sweet, this effect looks the same as the effect in 
#Figure 2A of Newman et al. (2014) JoN paper  Which is nice to see. 
#plot with ggplot:
source("summarySE.R") 
source("summarySEwithin.R") #function to calculate Std.Er of mean
source("normDataWithin.R")
plotdata <- summarySEwithin(DF_DistHem_wide, measurevar="SetSizeEffect", withinvars=c("DistHem"), betweenvars="dat1utr", idvar="ID")
ggplot(plotdata, aes(x=DistHem, y=SetSizeEffect, fill=dat1utr)) +
  geom_bar(position=position_dodge(.9), colour="Black", stat="identity") + 
  geom_errorbar(position=position_dodge(.9), width=.3, aes(ymin=SetSizeEffect-se, ymax=SetSizeEffect+se)) + #can change "se" to "ci" if I want to use 95%ci instead
  geom_hline(yintercept=0) +  coord_cartesian(ylim = c(-15, 65)) +
  xlab("Distractor Hemifield") + ylab("Set Size Effect (ms)") +
  theme(axis.title.x = element_text(face="bold", size=12),
        axis.text.x  = element_text(face="bold", angle=0,  size=12)) +
  theme(axis.title.y = element_text(face="bold", size=12),
        axis.text.y  = element_text(angle=0, vjust=0.5, size=12)) +
  theme(legend.title = element_text(size=11, face="bold")) +
  theme(legend.text = element_text(size = 11, face = "bold")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"))  # - good they both look the same - I just wanted to check to make sure the fitted "effects" plot looked the same as plotting the raw data using ggplot


##############################################################################
##############################################################################
####                                                                      ####
####            ~~~~~~~         Cueing data           ~~~~~~~             ####
####                                                                      ####
##############################################################################
##############################################################################

##### 1. Import single trial cueing data
cuedata <- read.csv("S:/R-MNHS-SPP/Bellgrove-data/CogBattery_UQ&Monash/AllParticipantData_cueing_CFQCAARS.csv", header=T, na.strings="#N/A")
#na.strings="#N/A" converts all #N/As in the document to NA.

# 2. cueing Data: Look at Column names:
colnames(cuedata)
#Rename some of the columns:
cuedata<-rename(cuedata, c("ResponseTime.ms."="RT", "ResponseStatus"="Accuracy"))
colnames(cuedata)

# 3. Remove any row that has a NAs in RT data:
cuedata<-cuedata[complete.cases(cuedata$RT),] 

# 4. Check number of Trials for each participant by running the function 'length', 
#Have a look at this
#cueing: everybody has 288 trials to start with
num_trials1 <- ddply(cuedata, c("ID"), summarise,
                     Trials    = length(RT))
summary(num_trials1$Trials)

# 5.Calculate each participant's accuracy as a %
cue_accuracy_checker <- ddply(cuedata, c("ID"), summarise,
                              Hits  = sum(Accuracy=="Correct"),
                              Misses = sum(Accuracy=="Incorrect"))
cue_accuracy_checker$Total=cue_accuracy_checker$Hits+cue_accuracy_checker$Misses
cue_accuracy_checker$Accuracy=(cue_accuracy_checker$Hits/cue_accuracy_checker$Total)*100
summary(cue_accuracy_checker$Accuracy)
#Plot the accuracy distribution 
hist(cue_accuracy_checker$Accuracy)

# 6. Remove any one who is guessing
#
#This task has 50% chance of guessing each trial correctly
#a binom.test shows that participants must have at least 161 correct out of
#the 288 trials in order for us to be 95% confident that they were not guessing
#the answer
binom.test(x=161, n=288, p=0.5)
#Most participant's have greater than 178 trials correct, check who doesn't:
ThoseGuessing<-cue_accuracy_checker[cue_accuracy_checker$Hits<161,]
ThoseGuessing
#make a new logical column in the main dataframe for whether each
#participant was guessing or not (TRUE/FALSE):
for (i in 1:length(cuedata$ID)) {cuedata$guessing[i]<-any(ThoseGuessing$ID==cuedata$ID[i])}

# 7. Make the required factors:
cuedata$ID <- factor(cuedata$ID)
#Valid/Invalid/Bilateral (neutral)
cuedata$ExptCond <- factor(cuedata$ExptCond)
#Left vs right
cuedata$ConditionTag1 <- factor(cuedata$ConditionTag1)
cuedata$Trial <- factor(cuedata$Trial)
cuedata$Site <- factor(cuedata$Site)
cuedata$dat1utr <- factor(cuedata$dat1utr)
cuedata$dat1i8 <- factor(cuedata$dat1i8)
cuedata$dat1haplo <- factor(cuedata$dat1haplo)
cuedata$drd4 <- factor(cuedata$drd4)

#Rename factor Levels:
cuedata$dat1utr<- revalue(cuedata$dat1utr, c("0"="zero","1"="one","2"="two")) 

# 8. Tidy data
#
# Look at RT distribution for each DAT1 group
ggplot(cuedata, aes(RT))  + geom_histogram(aes(y=..count..), colour="black", fill="white") + facet_wrap(~ dat1utr)

#Remove incorrect trials
cuedata<-cuedata[cuedata$Accuracy=="Correct",]

#Remove trials where RT too fast 
cuedata<-cuedata[!cuedata$RT<200,]
#....or too slow to be real RTs
cuedata<-cuedata[!cuedata$RT>2000,]

#Now have a look at the RT distribution:
hist(cuedata$RT)
cuedata$log_RT<-log(cuedata$RT) #log
hist(cuedata$log_RT)

#Remove CAARS data that is <1
cuedata<-cuedata[!cuedata$CAARS_H<1,]
summary(cuedata$CAARS_H)

##### Z-score each participant's log(RT) data inside the task Conditions ####
#ID by Condition (valid/invalid/neutral) by location (left/right)
cuedata$IDbyExptCondbyConditionTag1<-interaction(cuedata$ID, cuedata$ExptCond, cuedata$ConditionTag1) 

#calculate mean and sd 
m <- tapply(cuedata$log_RT,cuedata$IDbyExptCondbyConditionTag1,mean)
s <- sqrt(tapply(cuedata$log_RT,cuedata$IDbyExptCondbyConditionTag1,var))
#calculate RT.Z and save it inside data.frame
cuedata$RT.Z <- (cuedata$log_RT-m[cuedata$IDbyExptCondbyConditionTag1])/s[cuedata$IDbyExptCondbyConditionTag1]
#check that Z scores have mean=0 and std=1 
RT.Z_checker <- ddply(cuedata, c("ID", "ExptCond", "ConditionTag1"), summarise,
                      N    = length(RT.Z ),
                      mean = round(mean(RT.Z )),
                      sd   = sd(RT.Z ),
                      se   = sd / sqrt(N))
#Check to make sure the z-transformation worked (i.e. should see Z mean=0 and SD=1)
summary(RT.Z_checker$mean)
summary(RT.Z_checker$sd)

##Remove trials where absolute RT.Z>3 (i.e. remove outlier RTs)
cuedata<-cuedata[!abs(cuedata$RT.Z)>3,]

#Plot log(RT) historgam again after outlier removal 
hist(cuedata$log_RT)

#Check number of Trials for each participant by running the function 'length', 
#on "cuedata$RT" for each group, broken down by ID
cue_num_trials2 <- ddply(cuedata, c("ID"), summarise,
                     Trials    = length(RT))

#Plot the distribution of each participant's number of trials left after kicking out Incorrect and too fast/slow trials:
hist(num_trials2$Trials)

## Kick out Subjs based on CFQ and CAARS scores 
## (need to actually import these scores, but for now just kick out those who 
## I know have bad CFQ and CAARS scores from previous J.Neuro paper)
## Also added in Subj86 who has only 24 trials and no other data.

toBeRemoved<-which(cuedata$Subj=="subject240" | cuedata$Subj=="subject341" | cuedata$Subj=="subject358" | cuedata$Subj=="subject241")
#**For some reason this removes all of my data
# cuedata<-cuedata[-toBeRemoved,]

####################################
###  Crossed vs Nested Factors:  ###
####################################


#Test if ID is nested within Site:
with(cuedata, isNested(ID,Site)) #SO yes, ID is nested within Site...
#Print a table of this nesting:
xtabs(~ ID + Site, cuedata, drop = TRUE, sparse = TRUE)
# Because each Participant ("ID")" occurs within one and only one level of 
# Site we say that ID is nested within Site

# The task factors (SetSize, TargetHem, DisplayType, etc) are crossed or 
# partially crossed factors, "crossed" means that they have an observation for 
# each combination level each of the other factors (some factors might be partially 
# crossed in some participants if there is a lot of data missing due to bad accuracy)

# For Multilevel modeling with the lme4 package, nested factors are represended 
# like this "(1 | Site/ID) - which means ID is nested inside 
# Site", while crossed and partially cross factors are represented like 
# this (1 |SetSize) +(1|TargetHem) +(1|DisplayType), etc.

####################################

#Remove the DAT1 rows that are NA:
cuedata2<-cuedata[complete.cases(cuedata$dat1utr),] 
#Remove the CFQ rows that are NA:
cuedata2<-cuedata2[complete.cases(cuedata2$CFQ),] 
#Remove the CAARS (total ADHD) rows that are NA:
cuedata2<-cuedata2[complete.cases(cuedata2$CAARS_H),] 

summary(cuedata2$Site) #Summarize numbers of observations still left after NAs have been excluded from the dataset

############################################################################################################
#######    Aggregate to mean instead of trial-by-trial and check for partcicipant-level outliers     #######
############################################################################################################
#  This will alow us to calculate set-size effect which can not be measured on the single trial level      #
#
DF_collapsed<-ddply(cuedata2, .(ID, Site, dat1utr, ExptCond, ConditionTag1, CFQ, CAARS_H), summarise, RT=mean(RT))
ggplot(DF_collapsed, aes(RT))  + geom_histogram(aes(y=..count..), colour="black", fill="white") + facet_wrap(~ dat1utr)

# Collapse to each participant so I can Z-score individual participants' overall log(RT)s:
DF_collapsed_IDs<-ddply(cuedata2, .(ID, dat1utr), summarise, RT=mean(RT))
DF_collapsed_IDs$RT.Z<-scale(DF_collapsed_IDs$RT)
ggplot(DF_collapsed_IDs, aes(RT.Z))  + geom_histogram(aes(y=..count..), colour="black", fill="white")
#Split by DAT1 group
ggplot(DF_collapsed_IDs, aes(RT, fill=dat1utr, colour=dat1utr)) + geom_density(alpha=.05) 

#######    Calculate Cue effects   #######
require(reshape2)
#Bring Cueing Factors up into wide format and then calculate Cueing effects
DF_CueEffect_wide <- dcast(DF_collapsed, ID + Site + dat1utr + ConditionTag1 + CFQ + CAARS_H ~ ExptCond, value.var="RT", fun.aggregate=mean)
#Calculate CueCost (invalid - bilateral):
DF_CueEffect_wide$CueCost<-DF_CueEffect_wide$Invalid - DF_CueEffect_wide$Bilateral ########### BETH, this is an example of calculating the derived "Setsize effect" (mean RTs from trials with Eight minus those from from trials with Eight). The other tasks have similar derived measures you can calculate, such as the valid/invalid cueing effect 
summary(DF_CueEffect_wide$CueCost)
#Calculate CueCBenefit (bilateral - valid):
DF_CueEffect_wide$CueBenefit<-DF_CueEffect_wide$Bilateral - DF_CueEffect_wide$Valid ########### BETH, this is an example of calculating the derived "Setsize effect" (mean RTs from trials with Eight minus those from from trials with Eight). The other tasks have similar derived measures you can calculate, such as the valid/invalid cueing effect 
summary(DF_CueEffect_wide$CueBenefit)

###Now try model the collapsed CueCost data
random_intercepts_only_CueCost<-lmer(CueCost ~ 1 + (1 | Site/ID) +(1|ConditionTag1), data = DF_CueEffect_wide, na.action = na.omit, REML=FALSE)

ConditionTag1<-update(random_intercepts_only_CueCost, .~. + ConditionTag1)
dat1utr<-update(ConditionTag1, .~. + dat1utr)
ConditionTag1_by_dat1utr<-update(dat1utr, .~. + dat1utr*ConditionTag1)
CFQ<-update(ConditionTag1_by_dat1utr, .~. + CFQ)
CAARS_H<-update(ConditionTag1_by_dat1utr, .~. + CAARS_H)
ConditionTag1_by_CFQ<-update(CAARS_H, .~. + CFQ*ConditionTag1)
ConditionTag1_by_CAARS_H<-update(ConditionTag1_by_CFQ, .~. + CAARS_H*ConditionTag1)
dat1utr_by_CFQ<-update(ConditionTag1_by_CAARS_H, .~. + CFQ*dat1utr)
dat1utr_by_CAARS_H<-update(dat1utr_by_CFQ, .~. + CAARS_H*dat1utr)
ConditionTag1_by_dat1utr_by_CFQ<-update(dat1utr_by_CAARS_H, .~. + CFQ*dat1utr*ConditionTag1)
ConditionTag1_by_dat1utr_by_CAARS_H<-update(ConditionTag1_by_dat1utr_by_CFQ, .~. + CAARS_H*dat1utr*ConditionTag1)

anova(random_intercepts_only_CueCost, ConditionTag1, dat1utr, ConditionTag1_by_dat1utr,CFQ,CAARS_H,ConditionTag1_by_CFQ,ConditionTag1_by_CAARS_H,dat1utr_by_CFQ,dat1utr_by_CAARS_H,ConditionTag1_by_dat1utr_by_CFQ,ConditionTag1_by_dat1utr_by_CAARS_H)

# CFQ                                    11 14048 14104 -7012.8    14026 0.2534      1     0.6147  
# CAARS_H                                12 29467 29538 -14722     29443 8.6401      0    < 2e-16 ***
# dat1utr_by_CFQ                         17 29468 29567 -14717     29434 8.3232      2    0.01558 *  
# dat1utr_by_CFQ                         16 14054 14136 -7010.9    14022 1.3733      2     0.5033    
# dat1utr_by_CAARS_H                     18 14050 14143 -7006.9    14014 7.8743      2     0.0195 *  
# ConditionTag1_by_dat1utr_by_CFQ        20 14054 14157 -7006.9    14014 0.1252      2     0.9393    
# ConditionTag1_by_dat1utr_by_CAARS_H    22 14056 14170 -7006.0    14012 1.8114      2     0.4043  


###Now try model the collapsed CueBenefit data
random_intercepts_only_CueBenefit<-lmer(CueBenefit ~ 1 + (1 | Site/ID) +(1|ConditionTag1), data = DF_CueEffect_wide, na.action = na.omit, REML=FALSE)

ConditionTag1<-update(random_intercepts_only_CueBenefit, .~. + ConditionTag1)
dat1utr<-update(ConditionTag1, .~. + dat1utr)
ConditionTag1_by_dat1utr<-update(dat1utr, .~. + dat1utr*ConditionTag1)
CFQ<-update(ConditionTag1_by_dat1utr, .~. + CFQ)
CAARS_H<-update(ConditionTag1_by_dat1utr, .~. + CAARS_H)
ConditionTag1_by_CFQ<-update(CAARS_H, .~. + CFQ*ConditionTag1)
ConditionTag1_by_CAARS_H<-update(ConditionTag1_by_CFQ, .~. + CAARS_H*ConditionTag1)
dat1utr_by_CFQ<-update(ConditionTag1_by_CAARS_H, .~. + CFQ*dat1utr)
dat1utr_by_CAARS_H<-update(dat1utr_by_CFQ, .~. + CAARS_H*dat1utr)
ConditionTag1_by_dat1utr_by_CFQ<-update(dat1utr_by_CAARS_H, .~. + CFQ*dat1utr*ConditionTag1)
ConditionTag1_by_dat1utr_by_CAARS_H<-update(ConditionTag1_by_dat1utr_by_CFQ, .~. + CAARS_H*dat1utr*ConditionTag1)

anova(random_intercepts_only_CueBenefit, ConditionTag1, dat1utr, ConditionTag1_by_dat1utr,CFQ,CAARS_H,ConditionTag1_by_CFQ,ConditionTag1_by_CAARS_H,dat1utr_by_CFQ,dat1utr_by_CAARS_H,ConditionTag1_by_dat1utr_by_CFQ,ConditionTag1_by_dat1utr_by_CAARS_H)

#CFQ                                 11 13793 13850 -6885.7    13771 0.0001      1     0.9908    
#CAARS_H                             11 13793 13850 -6885.7    13771 0.0641      0     <2e-16 ***
#ConditionTag1_by_CFQ                13 13796 13864 -6885.2    13770 0.9625      2     0.6180    
#ConditionTag1_by_CAARS_H            14 13798 13871 -6885.1    13770 0.0588      1     0.8084    
#dat1utr_by_CFQ                      16 13801 13883 -6884.4    13769 1.5656      2     0.4571    
#dat1utr_by_CAARS_H                  18 13802 13895 -6883.1    13766 2.6033      2     0.2721    
#ConditionTag1_by_dat1utr_by_CFQ     20 13804 13907 -6882.0    13764 2.1934      2     0.3340    
#ConditionTag1_by_dat1utr_by_CAARS_H 22 13808 13921 -6881.9    13764 0.1803      2     0.9138  


############################################################################################################
#######                               Plot significant effects                                       #######
############################################################################################################

require(effects)
# Plot CueCost CAARS_H  p< 2e-16 ***
#...with effects package:
eff.CAARS_H_CueCost <- Effect(c("CAARS_H"), CAARS_H)
plot(eff.CAARS_H_CueCost, main=NULL, rug=T, ticks.x=NULL)
#now plot the same effect with ggplot: 
ggplot(DF_CueEffect_wide, aes(x=DF_CueEffect_wide$CAARS_H, y=DF_CueEffect_wide$CueCost)) +
  stat_smooth(level = 0.95,size=1) 

# Plot CueCost dat1utr_by_CFQ p=0.0195 * 
eff.dat1utr_by_CAARS_H_CueCost <- Effect(c("dat1utr", "CAARS_H"), dat1utr_by_CAARS_H)
plot(eff.dat1utr_by_CAARS_H_CueCost, layout=c(3,1), main=NULL, rug=T, multiline =F, alternating=F, ticks.x=NULL) 

# Plot CueBenefit CAARS_H  p< 2e-16 ***
#...with effects package:
eff.CAARS_H_CueBenefit <- Effect(c("CAARS_H"), CAARS_H)
plot(eff.CAARS_H_CueBenefit, main=NULL, rug=T, ticks.x=NULL)
#now plot the same effect with ggplot: 
ggplot(DF_CueEffect_wide, aes(x=DF_CueEffect_wide$CAARS_H, y=DF_CueEffect_wide$CueBenefit)) +
  stat_smooth(level = 0.95,size=1)




##############################################################################
##############################################################################
####                                                                      ####
####            ~~~~~~~         Flanker data           ~~~~~~~            ####
####                                                                      ####
##############################################################################
##############################################################################

##### 1. Import single trial flanker data
flankerdata <- read.csv("S:/R-MNHS-SPP/Bellgrove-data/CogBattery_UQ&Monash/AllParticipantData_flanker_CFQCAARS.csv", header=T, na.strings="#N/A")
#na.strings="#N/A" converts all #N/As in the document to NA.

# 2. flanker data: Look at Column names:
colnames(flankerdata)
#Rename some of the columns:
flankerdata<-rename(flankerdata, c("ResponseType"="Accuracy", "TrialNo"="Trial"))
colnames(flankerdata)

# 3. Remove any row that has a NAs in RT data:
flankerdata<-flankerdata[complete.cases(flankerdata$RT),] 

# 4. Check number of Trials for each participant by running the function 'length', 
#Have a look at this
#cueing: everybody has 288 trials to start with
flanker_num_trials1 <- ddply(flankerdata, c("ID"), summarise,
                     Trials    = length(RT))
summary(flanker_num_trials1$Trials)

#Calculate each participant's accuracy as a %
flanker_accuracy_checker <- ddply(flankerdata, c("ID"), summarise,
                          Hits  = sum(Accuracy=="Correct"),
                          Misses = sum(Accuracy=="Incorrect"))
flanker_accuracy_checker$Total=flanker_accuracy_checker$Hits+flanker_accuracy_checker$Misses
flanker_accuracy_checker$Accuracy=(flanker_accuracy_checker$Hits/flanker_accuracy_checker$Total)*100
summary(flanker_accuracy_checker$Accuracy)
#Plot the accuracy distribution 
hist(flanker_accuracy_checker$Accuracy)

# 6. Remove any one who is guessing
#
#This task has 50% chance of guessing each trial correctly
#a binom.test shows that participants must have at least 161 correct out of
#the 288 trials in order for us to be 95% confident that they were not guessing
#the answer
binom.test(x=161, n=288, p=0.5)
#Most participant's have greater than 178 trials correct, check who doesn't:
ThoseGuessing<-flanker_accuracy_checker[flanker_accuracy_checker$Hits<161,]
ThoseGuessing
#make a new logical column in the main dataframe for whether each
#participant was guessing or not (TRUE/FALSE):
for (i in 1:length(flankerdata$ID)) {flankerdata$guessing[i]<-any(ThoseGuessing$ID==flankerdata$ID[i])}

hist(flanker_accuracy_checker$Accuracy)

# 7. Make the required factors:
flankerdata$ID <- factor(flankerdata$ID)
flankerdata$Trial <- factor(flankerdata$Trial)
flankerdata$Site <- factor(flankerdata$Site)
flankerdata$PostGoLabel <- factor(flankerdata$PostGoLabel)
flankerdata$dat1utr <- factor(flankerdata$dat1utr)
flankerdata$dat1i8 <- factor(flankerdata$dat1i8)
flankerdata$dat1haplo <- factor(flankerdata$dat1haplo)
flankerdata$drd4 <- factor(flankerdata$drd4)

#Rename factor Levels:
flankerdata$dat1utr<- revalue(flankerdata$dat1utr, c("0"="zero","1"="one","2"="two")) 

# 8. Tidy data
#
# Look at RT distribution for each DAT1 group
ggplot(flankerdata, aes(RT))  + geom_histogram(aes(y=..count..), colour="black", fill="white") + facet_wrap(~ dat1utr)

#Now kick out those participants who were guessing 
flankerdata<-flankerdata[flankerdata$guessing!=TRUE,] 

#Kick out trials with stim. timing problems
flankerdata<-flankerdata[flankerdata$TimingStatus=="OK",]

#Remove incorrect trials
flankerdata<-flankerdata[flankerdata$Accuracy=="Correct",]

#Remove trials where RT too fast 
flankerdata<-flankerdata[!flankerdata$RT<200,]
#....or too slow to be real RTs
flankerdata<-flankerdata[!flankerdata$RT>2000,]
#Now have a look at the RT distribution:
hist(flankerdata$RT)



flankerdata$log_RT<-log(flankerdata$RT) #log

# 3. Remove any row that has a NAs in RT data:
flankerdata<-flankerdata[complete.cases(flankerdata$log_RT),] 
hist(flankerdata$log_RT)

##### Z-score each participant's log(RT) data inside the task Conditions ####
flankerdata$IDbyPostGoLabel<-interaction(flankerdata$ID, flankerdata$PostGoLabel) #ID by PostGoLabel
#calculate mean and sd 
m <- tapply(flankerdata$log_RT,flankerdata$IDbyPostGoLabel,mean)
s <- sqrt(tapply(flankerdata$log_RT,flankerdata$IDbyPostGoLabel,var))
#calculate RT.Z and save it inside data.frame
flankerdata$RT.Z <- (flankerdata$log_RT-m[flankerdata$IDbyPostGoLabel])/s[flankerdata$IDbyPostGoLabel]
#check that Z scores have mean=0 and std=1 
RT.Z_checker <- ddply(flankerdata, c("ID", "PostGoLabel"), summarise,
                      N    = length(RT.Z ),
                      mean = round(mean(RT.Z )),
                      sd   = sd(RT.Z ),
                      se   = sd / sqrt(N))
#Check to make sure the z-transformation worked (i.e. should see Z mean=0 and SD=1)
summary(RT.Z_checker$mean)
summary(RT.Z_checker$sd)


##Remove trials where absolute RT.Z>3 (i.e. remove outlier RTs)
flankerdata<-flankerdata[!abs(flankerdata$RT.Z)>3,]

#Plot log(RT) historgam again after outlier removal 
hist(flankerdata$log_RT)

#Check number of Trials for each participant by running the function 'length', 
#on "flankerdata$RT" for each group, broken down by ID
flanker_num_trials2 <- ddply(flankerdata, c("ID"), summarise,
                     Trials    = length(RT))

summary(flanker_num_trials2$Trials)

#Plot the distribution of each participant's number of trials left after kicking out Incorrect and too fast/slow trials:
hist(flanker_num_trials2$Trials)

flankerdata<-flankerdata[!flankerdata$CAARS_H<1,]
summary(flankerdata$CAARS_H)

####################################
###  Crossed vs Nested Factors:  ###
####################################


#Test if ID is nested within Site:
with(flankerdata, isNested(ID,Site)) #SO yes, ID is nested within Site...
#Print a table of this nesting:
xtabs(~ ID + Site, flankerdata, drop = TRUE, sparse = TRUE)
# Because each Participant ("ID")" occurs within one and only one level of 
# Site we say that ID is nested within Site

# The task factors (SetSize, TargetHem, DisplayType, etc) are crossed or 
# partially crossed factors, "crossed" means that they have an observation for 
# each combination level each of the other factors (some factors might be partially 
# crossed in some participants if there is a lot of data missing due to bad accuracy)

# For Multilevel modeling with the lme4 package, nested factors are represended 
# like this "(1 | Site/ID) - which means ID is nested inside 
# Site", while crossed and partially cross factors are represented like 
# this (1 |SetSize) +(1|TargetHem) +(1|DisplayType), etc.

####################################

#Remove the DAT1 rows that are NA:
flankerdata2<-flankerdata[complete.cases(flankerdata$dat1utr),] 
#Remove the CFQ rows that are NA:
flankerdata2<-flankerdata2[complete.cases(flankerdata2$CFQ),] 
#Remove the CAARS (total ADHD) rows that are NA:
flankerdata2<-flankerdata2[complete.cases(flankerdata2$CAARS_H),] 

summary(flankerdata2$Site) #Summarize numbers of observations still left after NAs have been excluded from the dataset

############################################################################################################
#######    Aggregate to mean instead of trial-by-trial and check for partcicipant-level outliers     #######
############################################################################################################
#  This will alow us to calculate set-size effect which can not be measured on the single trial level      #
#
DF_collapsed<-ddply(flankerdata2, .(ID, Site, dat1utr, PreGoLabel, CFQ, CAARS_H), summarise, RT=mean(RT))
ggplot(DF_collapsed, aes(RT))  + geom_histogram(aes(y=..count..), colour="black", fill="white") + facet_wrap(~ dat1utr)

# Collapse to each participant so I can Z-score individual participants' overall log(RT)s:
DF_collapsed_IDs<-ddply(flankerdata2, .(ID, dat1utr), summarise, RT=mean(RT))
DF_collapsed_IDs$RT.Z<-scale(DF_collapsed_IDs$RT)
ggplot(DF_collapsed_IDs, aes(RT.Z))  + geom_histogram(aes(y=..count..), colour="black", fill="white")
#Split by DAT1 group
ggplot(DF_collapsed_IDs, aes(RT, fill=dat1utr, colour=dat1utr)) + geom_density(alpha=.05) 

#######    Calculate incongruence cost   #######
require(reshape2)
#Bring Cueing Factors up into wide format and then calculate Flanker effect
DF_IncongCost_wide <- dcast(DF_collapsed, ID + Site + dat1utr + CFQ + CAARS_H ~ PreGoLabel, value.var="RT", fun.aggregate=mean)
#Collapse left and right for 1) congruent, 2) incongruent, 3) neutral using reshape2 melt function

#Collapse left & right congruent
DF_IncongCost_wide<-melt(DF_IncongCost_wide, id=c("ID", "Site", "dat1utr", "CFQ", "CAARS_H","Left Incongruent","Left Neutral","Right Incongruent","Right Neutral"), measured=c("Left Congruent","Right Congruent"))
DF_IncongCost_wide<-rename(DF_IncongCost_wide, c("variable"="CongruentLR", "value"="CongruentRT"))

DF_IncongCost_wide<-melt(DF_IncongCost_wide, id=c("ID", "Site", "dat1utr", "CFQ", "CAARS_H","Left Neutral","Right Neutral","CongruentLR","CongruentRT"), measured=c("Left Incongruent","Right Incongruent"))
DF_IncongCost_wide<-rename(DF_IncongCost_wide, c("variable"="IncongruentLR", "value"="IncongruentRT"))

DF_IncongCost_wide<-melt(DF_IncongCost_wide, id=c("ID", "Site", "dat1utr", "CFQ", "CAARS_H","CongruentLR","CongruentRT","IncongruentLR","IncongruentRT"), measured=c("Left Neutral","Right Neutral"))
DF_IncongCost_wide<-rename(DF_IncongCost_wide, c("variable"="NeutralLR", "value"="NeutralRT"))

#Calculate Flanker Cost (incongruent - neutral):
DF_IncongCost_wide$FlankerCost<-DF_IncongCost_wide$IncongruentRT - DF_IncongCost_wide$NeutralRT  
summary(DF_IncongCost_wide$FlankerCost)

###Now try model the collapsed FlankerCost data
#First make a baseline model with random intercepts and no fixed effects:
random_intercepts_only_FlankerCost<-lmer(FlankerCost ~ 1 + (1 | Site/ID) + (1|IncongruentLR) + (1|NeutralLR), data = DF_IncongCost_wide, na.action = na.omit, REML=FALSE)

#Then we add in the fixed effect of Flanker Cost:
IncongruentLR<-update(random_intercepts_only_FlankerCost, .~. + IncongruentLR)
NeutralLR<-update(IncongruentLR, .~. + NeutralLR)
dat1utr<-update(NeutralLR, .~. + dat1utr)
CFQ<-update(dat1utr, .~. + CFQ)
CAARS_H<-update(CFQ, .~. + CAARS_H)

##Testing for interactions
IncongruentLR_by_NeutralLR<-update(dat1utr, .~. + IncongruentLR*NeutralLR)

dat1utr_by_IncongruentLR<-update(IncongruentLR_by_NeutralLR, .~. + dat1utr*IncongruentLR)
dat1utr_by_NeutralLR<-update(dat1utr_by_IncongruentLR, .~. + dat1utr*NeutralLR)
NeutralLR_by_IncongruentLR_by_dat1utr<-update(dat1utr_by_IncongruentLR, .~. + NeutralLR*IncongruentLR*dat1utr)

anova(random_intercepts_only_FlankerCost, IncongruentLR, NeutralLR, dat1utr, CFQ, CAARS_H, IncongruentLR_by_NeutralLR, dat1utr_by_IncongruentLR,
      dat1utr_by_IncongruentLR, dat1utr_by_NeutralLR, NeutralLR_by_IncongruentLR_by_dat1utr)


# MODEL6: IncongruentLR + NeutralLR + dat1utr + CFQ + CAARS_H
#         Df   AIC   BIC logLik deviance  Chisq Chi Df   Pr(>Chisq)
# MODEL6  12 46101 46180 -23039    46077 4.1805      1    0.04089 *

