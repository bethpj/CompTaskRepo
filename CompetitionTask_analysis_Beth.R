#Set working directory
# setwd(("//ad.monash.edu/home/User063/bethj/Documents/GitHub/CompTaskRepo"))
setwd(("C:/Users/Dan/Documents/GitHub/CompTask"))
# setwd(("C:/GitHub/CompTask"))

####################################
####################################
#
# Key questions to address:
#
# 1. explore the relationship between CFQ and each of the set size effect, cuing cost, cuing benefit and flanker interference effects.
#
# 2. explore the relationship between DAT1 and CFQ/ADHD scores and the cuing data (probably nothing)
#   and competition data (without regard to target position).
#
####################################
####################################
#
######  FIRST TIME ONLY -
#####  Install some handy packages:
####  Remove the "# "in front of the line below and run the code. Replace the # after installing the packages, otherwise the R 
#will reinstall the packates every time you run the script
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

###### Import single trial data:
# data <- read.csv("C:/subj/AllParticipantData_comp_CFQCAARS.csv", header=T)
# data <- read.csv("//ad.monash.edu/home/User088/newmand/Desktop/AllParticipantData_comp_CFQCAARS/AllParticipantData_comp_CFQCAARS.csv", header=T, na.strings="#N/A")
data <- read.csv("S:/R-MNHS-SPP/Bellgrove-data/CogBattery_UQ&Monash/AllParticipantData_comp_CFQCAARS.csv", header=T, na.strings="#N/A")
#na.strings="#N/A" converts all #N/As in the document to NA.

# Look at Column names:
colnames(data)
#Rename some of the columns how I like them:
data<-rename(data, c("ResponseTime"="RT", "ResponseStatus"="Accuracy"))

#Remove the one row that has a NAs:
data<-data[complete.cases(data$RT),] 
#Check number of Trials for each participant by running the function 'length', 
#on "data$RT", broken down by ID
num_trials1 <- ddply(data, c("ID"), summarise,
                     Trials    = length(RT))
#Have a look at this
summary(num_trials1$Trials)#so everybody has 320 trials to start with

#Calculate each participant's accuracy as a %
Accuracy_checker <- ddply(data, c("ID"), summarise,
                          Hits  = sum(Accuracy=="Correct"),
                          Misses = sum(Accuracy=="Incorrect"))
Accuracy_checker$Total=Accuracy_checker$Hits+Accuracy_checker$Misses
Accuracy_checker$Accuracy=(Accuracy_checker$Hits/Accuracy_checker$Total)*100
summary(Accuracy_checker$Accuracy)
#Plot the accuracy distribution 
hist(Accuracy_checker$Accuracy)

#This task has 50% chance of guessing each trial correctly
#a binom.test shows that participants must have at least 179 correct out of
#the 320 trials in order for us to be 95% confident that they were not guessing
#the answer
binom.test(x=179, n=320, p=0.5)
#Most participant's have greater than 178 trials correct, check who doesn't:
ThoseGuessing<-Accuracy_checker[Accuracy_checker$Hits<178,]
ThoseGuessing
#make a new logical column in the main dataframe for whether each
#participant was guessing or not (TRUE/FALSE):
for (i in 1:length(data$ID)) {data$guessing[i]<-any(ThoseGuessing$ID==data$ID[i])}


#Make the required factors:
data$ID <- factor(data$ID)
data$SetSize <- factor(data$SetSize)
data$TargetHem <- factor(data$TargetHem)
data$TargetLoc <- factor(data$TargetLoc)
data$DisplayType <- factor(data$DisplayType)
data$TargetType <- factor(data$TargetType)
data$Trial <- factor(data$Trial)
data$Site <- factor(data$Site)
data$dat1utr <- factor(data$dat1utr)
data$dat1i8 <- factor(data$dat1i8)
data$dat1haplo <- factor(data$dat1haplo)
data$drd4 <- factor(data$drd4)

#these are not really factors, they are things we measured from the participant:
# #Cognitive failures  
# data$CFQ <-factor(data$CFQ)
# #CAARS Inattention
# data$CAARS_A <-factor(data$CAARS_A)
# #CAARS hyperactivity
# data$CAARS_B <-factor(data$CAARS_B)
# #CAARS ADHD total
# data$CAARS_H <-factor(data$CAARS_H)

#Rename factor Levels:
data$TargetHem <- revalue(data$TargetHem, c("1"="Left", "2"="Right")) #Hi Beth, please just check this is correct i.e. that "1" does indeed stand for "Left" here
data$DisplayType <- revalue(data$DisplayType, c("1"="Unilateral", "2"="Bilateral")) #Beth also better double check that "1" was actually "Unilateral
data$SetSize <- revalue(data$SetSize , c("4"="Four", "8"="Eight"))

data$dat1utr<- revalue(data$dat1utr, c("0"="zero","1"="one","2"="two")) 

summary(data$CFQ)
summary(data$CAARS_H)

#Look at RT distribution for each DAT1 group
ggplot(data, aes(RT))  + geom_histogram(aes(y=..count..), colour="black", fill="white") + facet_wrap(~ dat1utr)

#Now kick out those participants who were guessing 
data<-data[data$guessing!=TRUE,] 

#Kick out trials with stim. timing problems
data<-data[data$Timing.Status=="OK",]

#Remove incorrect trials
data<-data[data$Accuracy=="Correct",]

#Remove trials where RT too fast 
data<-data[!data$RT<200,]
#....or too slow to be real RTs
data<-data[!data$RT>2000,]
#Now have a look at the RT distribution:
hist(data$RT)

data$log_RT<-log(data$RT) #log
hist(data$log_RT)

#####Z-score each participant's log(RT) data inside the task Conditions####
data$IDbySetSizebyDisplayType<-interaction(data$ID, data$SetSize, data$DisplayType) #ID by SetSize by DisplayType
#calculate mean and sd 
m <- tapply(data$log_RT,data$IDbySetSizebyDisplayType,mean)
s <- sqrt(tapply(data$log_RT,data$IDbySetSizebyDisplayType,var))
#calculate RT.Z and save it inside data.frame
data$RT.Z <- (data$log_RT-m[data$IDbySetSizebyDisplayType])/s[data$IDbySetSizebyDisplayType]
#check that Z scores have mean=0 and std=1 
RT.Z_checker <- ddply(data, c("ID", "SetSize", "DisplayType"), summarise,
                      N    = length(RT.Z ),
                      mean = round(mean(RT.Z )),
                      sd   = sd(RT.Z ),
                      se   = sd / sqrt(N))
#Check to make sure the z-transformation worked (i.e. should see Z mean=0 and SD=1)
summary(RT.Z_checker$mean)
summary(RT.Z_checker$sd)

##Remove trials where absolute RT.Z>3 (i.e. remove outlier RTs)
data<-data[!abs(data$RT.Z)>3,]

#Plot log(RT) historgam again after outlier removal 
hist(data$log_RT)

#Check number of Trials for each participant by running the function 'length', 
#on "data$RT" for each group, broken down by ID
num_trials2 <- ddply(data, c("ID"), summarise,
                     Trials    = length(RT))

summary(num_trials2$Trials)
#Plot the distribution of each participant's number of trials
#left after kicking out Incorrect and too fast/slow trials:
hist(num_trials2$Trials)


## Kick out Subjs based on CFQ and CAARS scores 
## (need to actually import these scores, but for now just kick out those who 
##  I know have bad CFQ and CAARS scores from previous J.Neuro paper)
toBeRemoved<-which(data$Subj=="subject240" | data$Subj=="subject341" | data$Subj=="subject358" | data$Subj=="subject241")

data<-data[-toBeRemoved,]

#Calculate Distracter Hemifield (DistHem) as used in the 
#Newman et al. (2014) JoN paper, rather then using TargetHem:
data$TargetHembyDisplayType<-interaction(data$TargetHem, data$DisplayType)
data$DistHem <- revalue(data$TargetHembyDisplayType, c("Left.Unilateral"="Left", "Left.Bilateral"="Right", "Right.Unilateral"="Right", "Right.Bilateral"="Left"))


###########################################################################
#########################################################################
###Crossed vs Nested Factors:

#Test if ID is nested within Site:
with(data, isNested(ID,Site)) #SO yes, ID is nested within Site...
#Print a table of this nesting:
xtabs(~ ID + Site, data, drop = TRUE, sparse = TRUE)
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


########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################

#Remove the DAT1 rows that are NA:
data2<-data[complete.cases(data$dat1utr),] 
#Remove the CFQ rows that are NA:
data2<-data2[complete.cases(data2$CFQ),] 
#Remove the CAARS (total ADHD) rows that are NA:
data2<-data2[complete.cases(data2$CAARS_H),] 

summary(data2$Site) #Summarize numbers of observations still left after NAs have been excluded from the dataset

################################################################################################################
################Aggregate to mean instead of trial-by-trial and check for partcicipant-level outliers###########
############This will also alow me to calculate Set-size effect 
############which can not be measured on the single trial level
############for obvious reasons.
DF_collapsed<-ddply(data2, .(ID, Site, dat1utr, SetSize,TargetLoc,DisplayType,TargetHem, DistHem, CFQ, CAARS_H), summarise, RT=mean(RT))

ggplot(DF_collapsed, aes(RT))  + geom_histogram(aes(y=..count..), colour="black", fill="white") + facet_wrap(~ dat1utr)

#####collapse to each participant so I can Z-score individual participants' overall log(RT)s:
DF_collapsed_IDs<-ddply(data2, .(ID, dat1utr), summarise, RT=mean(RT))
DF_collapsed_IDs$RT.Z<-scale(DF_collapsed_IDs$RT)
ggplot(DF_collapsed_IDs, aes(RT.Z))  + geom_histogram(aes(y=..count..), colour="black", fill="white")
#Split by DAT1 group
ggplot(DF_collapsed_IDs, aes(RT, fill=dat1utr, colour=dat1utr)) + geom_density(alpha=.05) 

################Calculate Set-size effect by DistHem##############
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

########## Plot the significant effects:
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




########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
























# Hi Beth, so Mark is more interested in looking at the effects of Genetics on
# the derived measures from the task such as "SetSize 
# effect" for Comp task like I have done above here, 
#or "cueing validity effects" for cueing task, etc, rather 
# than looking at single trial RTs. But if you are interested, below is an example
# of single trial RT analysis of Comp task.


# lets look at modeling the single trial RT data as a function of TargetLoc:
#IDs are nested inside Site (1 | Site/ID)
#Also, TargetLocs (1-16) are nested inside TargetHem (left vs. right):
with(data2, isNested(TargetLoc,TargetHem))

#First make a baseline model with random intercepts and no fixed effects:
random_intercepts_only<-lmer(log(RT) ~ 1 + (1 | Site/ID) +(1|TargetHem/TargetLoc) + (1 |SetSize) 
                             +(1|DisplayType), data = data2, na.action = na.omit, REML=FALSE)
#Then we add in the fixed effect of SetSize, etc.:
SetSize<-update(random_intercepts_only, .~. + SetSize)
DisplayType<-update(SetSize, .~. + DisplayType)
TargetLoc<-update(DisplayType, .~. + TargetLoc)
dat1utr<-update(TargetLoc, .~. + dat1utr)
CFQ<-update(dat1utr, .~. + CFQ)
CAARS_H<-update(CFQ, .~. + CAARS_H)

##Testing for interactions
SetSize_by_TargetLoc<-update(dat1utr, .~. + SetSize*TargetLoc)
SetSize_by_DisplayType<-update(SetSize_by_TargetLoc, .~. + SetSize*DisplayType)
TargetLoc_by_DisplayType<-update(SetSize_by_DisplayType, .~. + TargetLoc*DisplayType)

dat1utr_by_TargetLoc<-update(TargetLoc_by_DisplayType, .~. + dat1utr*TargetLoc)
dat1utr_by_DisplayType<-update(dat1utr_by_TargetLoc, .~. + dat1utr*DisplayType)
dat1utr_by_SetSize<-update(dat1utr_by_DisplayType, .~. + dat1utr*SetSize)

SetSize_by_TargetLoc_by_dat1utr<-update(dat1utr_by_SetSize, .~. + SetSize*TargetLoc*dat1utr)
SetSize_by_DisplayType_by_dat1utr<-update(SetSize_by_TargetLoc_by_dat1utr, .~. + SetSize*DisplayType*dat1utr)
TargetLoc_by_DisplayType_by_dat1utr<-update(SetSize_by_DisplayType_by_dat1utr, .~. + TargetLoc*DisplayType*dat1utr)


SetSize_by_TargetLoc_by_DisplayType<-update(TargetLoc_by_DisplayType_by_dat1utr, .~. + SetSize*TargetLoc*DisplayType)
SetSize_by_TargetLoc_by_DisplayType_by_dat1utr<-update(SetSize_by_TargetLoc_by_DisplayType, .~. + SetSize*TargetLoc*DisplayType*dat1utr)


anova(random_intercepts_only, SetSize, TargetLoc, DisplayType, dat1utr,
      SetSize_by_TargetLoc, SetSize_by_DisplayType, TargetLoc_by_DisplayType,
      dat1utr_by_TargetLoc,dat1utr_by_DisplayType,dat1utr_by_SetSize,
      SetSize_by_TargetLoc_by_dat1utr,SetSize_by_DisplayType_by_dat1utr,TargetLoc_by_DisplayType_by_dat1utr,
      SetSize_by_TargetLoc_by_DisplayType,SetSize_by_TargetLoc_by_DisplayType_by_dat1utr)

#######when we ran the single trial model with TargetLoc above is says
####### "fixed-effect model matrix is rank deficient". This means there is not 
####### enough data to fit the model. I think this is because not every 
###### Position has observations accross all of the other vatiables, lets explore:


#Plot RT TargetLoc x dat1utr: dat1utr_by_TargetLoc
source("summarySE.R") 
source("summarySEwithin.R") #function to calculate Std.Er of mean
source("normDataWithin.R")
plotdata <- summarySEwithin(data2, measurevar="RT", withinvars=c("TargetLoc"), betweenvars="dat1utr", idvar="ID")
ggplot(plotdata, aes(x=TargetLoc, y=RT, fill=dat1utr)) +
    geom_bar(position=position_dodge(.9), colour="Black", stat="identity") + 
    geom_errorbar(position=position_dodge(.9), width=.3, aes(ymin=RT-ci, ymax=RT+ci)) + #can change "se" to "ci" if I want to use 95%ci instead
    geom_hline(yintercept=0) +  coord_cartesian(ylim = c(450, 900)) +
    xlab("TargetLoc") + ylab("Reaction-time (ms)") +
    theme(axis.title.x = element_text(face="bold", size=12),
          axis.text.x  = element_text(face="bold", angle=0,  size=12)) +
    theme(axis.title.y = element_text(face="bold", size=12),
          axis.text.y  = element_text(angle=0, vjust=0.5, size=12)) +
    theme(legend.title = element_text(size=11, face="bold")) +
    theme(legend.text = element_text(size = 11, face = "bold")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_blank(), axis.line = element_line(colour = "black"))  







