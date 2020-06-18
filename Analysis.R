# Script created by Frank D'Agostino 2020
# Analysis conducted for the COVID and Climate
# Project as part of the Coronavirus Visualization Team (CVT)

rm(list = ls())
# Clear console
library(ggplot2); library(plyr)
library(foreign); library(dplyr)
library(MASS); library(lme4)
library(glmmTMB); library(gamm4)
library(sandwich); library(lmtest)
library(rdrobust); library(ggpubr)
library(tidyverse); library(broom)
library(AICcmodavg)

#-------------------------------------------------------------------------------
# Data Cleanup
#-------------------------------------------------------------------------------

# Import the partially pre-processed data
air <- read.csv("countyAirPoll.csv"); head(air)

hosp <- read.csv("countyHosp.csv"); head(hosp)

# For air pollution, we have data spanning back to 2010, but we want
# to focus on data regarding 2019 and 2020, so we take a subset
recent <- subset(air, Year >= 2010 & Pollutant == "PM2.5"); head(recent)

# Clean the data by creating a dataframe of average AQI values for
# graphing and easy visualization
berk <- subset(recent, County == "Berkshire"); brist <- subset(recent, County == "Bristol")
# dukes <- subset(recent, County == "Dukes"); 
essex <- subset(recent, County == "Essex")
franklin <- subset(recent, County == "Franklin"); hampden <- subset(recent, County == "Hampden")
hamp <- subset(recent, County == "Hampshire"); middle <- subset(recent, County == "Middlesex")
norf <- subset(recent, County == "Norfolk"); plym <- subset(recent, County == "Plymouth")
suff <- subset(recent, County == "Suffolk"); worc <- subset(recent, County == "Worcester")

# Middlesex and Norfolk have sparse data
dataMean <- data.frame(Countys = c("Berkshire", "Bristol", # "Dukes", 
                                   "Essex", "Franklin", "Hampden", "Hampshire", "Middlesex",
                                    "Norfolk", "Plymouth", "Suffolk", "Worcester"),
                   AQI = c(mean(berk$AQI), mean(brist$AQI), # mean(dukes$AQI),
                           mean(essex$AQI), mean(franklin$AQI), mean(hampden$AQI),
                           mean(hamp$AQI), mean(middle$AQI), mean(norf$AQI),
                           mean(plym$AQI), mean(suff$AQI), mean(worc$AQI))); head(data)

#-------------------------------------------------------------------------------
# Basic Visualization
#-------------------------------------------------------------------------------

# Plotting AQI over time
ggplot(suff, aes(x=Year,y=AQI)) + stat_summary(fun.y = "mean",geom="point") +
  stat_summary(fun.y = "mean",geom="line")  +
  labs(x = "Year", y = "Air Quality Index", shape = "") +
  ggtitle("Progression of Air Quality (PM2.5) over the Last Decade in Suffolk County, MA")
  # + theme(legend.position="bottom") # + geom_vline(xintercept=2020)
  # shape= factor(Pollutant, labels = c("Rest of U.S.", "Arizona")))) +

# Use ggplot to visually demonstrate the relatively small differences
# in AQI (stands for air quality index)
ggplot(dataMean, aes(x=Countys, y=AQI, fill=Countys))+
  geom_bar(stat = "identity", colour="black", width=0.5)+labs(y= "Average AQI", x = "County")+
  ggtitle("Barplot of Average AQI Value at Each Massachusetts County in 2016-2020")+
  theme_classic()+
  scale_fill_brewer(palette="Set3")

# Produce a barplot demonstrating the number of COVID-19 hospitalizations
# per county using the builtin R barplot function
barplot(table(hosp$HospTotal, hosp$County), 
        col = "lightblue", 
        xlab = "Massachusetts County", 
        ylab = "Total Hospitalizations", 
        main="Total Number of COVID-19 Hospitalizations in Each Massachusetts County")

#-------------------------------------------------------------------------------
# Consolidating the Measurements into a Dataframe
#-------------------------------------------------------------------------------

# Using regex, we want to select air pollution data from only the months
# of January-May from 2020, since the hospitalization data spans from
# 4/30/20 to 5/22/20, which is the most recent date
count <- recent$County[grep("^([5][\\/]|[4][\\/]|[3][\\/]|[2][\\/]|[1][\\/]|)\\d+(\\/20)", recent$Date)]; head(count)
poll <- recent$AQI[grep("^([5][\\/]|[4][\\/]|[3][\\/]|[2][\\/]|[1][\\/]|)\\d+(\\/20)", recent$Date)]; head(poll)
mayAir <- data.frame(count, poll); head(mayAir)

# Use the aggregate function to group everything by counties
temp <- aggregate(x=hosp, by=list(hosp$County), FUN=mean); temp
aggHosp <- subset(temp, Group.1 != "Nantucket" & Group.1 != "Barnstable" & Group.1 != "Dukes"); aggHosp

allAir <- data.frame(dataMean$County, dataMean$AQI); head(allAir)
aggAir <- aggregate(x=allAir, by=list(allAir$dataMean.County), FUN=mean); aggAir

# Nantucket, Dukes, and Barnstable do not have enough Air data in April/May 2020,
# and so we remove them from the hospital data to have even comparisons

# We also add the confounding variables to them
c1 <- read.csv("confoundVars.csv"); c1
c2 <- read.csv("confoundPercent.csv"); c2
conf1 <- subset(c1, County != "Nantucket" & County != "Barnstable" & County != "Dukes"); conf1
conf2 <- subset(c2, County != "Nantucket" & County != "Barnstable" & County != "Dukes"); conf2

dataComb <- data.frame(County = aggAir$Group.1,
                       AQI = aggAir$dataMean.AQI, HospTotal = aggHosp$HospTotal,
                       ICU = aggHosp$ICU, Pop = conf1$Total.Pop,
                       Young = conf1$Under.18, Old = conf1$X65...Over,
                       Cancer = conf1$Lung.Cancer, Nonwhite = conf1$Non.White,
                       Poverty = conf1$Poverty.Estimate); dataComb

#-------------------------------------------------------------------------------
# Check for Preliminary Relationships
#-------------------------------------------------------------------------------

# Regress HospTotal with AQI and check for signifance
hosp.lm <- lm(HospTotal~AQI, data=dataComb); hosp.lm;
plot(dataComb$AQI, dataComb$HospTotal)
curve(-2.4*x + 129.7,col = "cyan", lwd=3, add = T)
# P-value is greater than 0.05, suggesting no signifance
# Additionally, more general hospitalizations are associated
# with a decreased AQI
coeftest(hosp.lm, vcov=vcovHC(hosp.lm,type="HC1"))

# Regress ICU with AQI and check for significance
icu.lm <- lm(ICU~AQI, data=dataComb); icu.lm;
plot(dataComb$AQI, dataComb$ICU)
curve(-0.7115*x + 36.3471,col = "cyan", lwd=3, add = T)
# P-value is greater than 0.05, suggesting no signifance
# Additionally, more ICU hospitalizations are associated
# with a decreased AQI
coeftest(icu.lm, vcov=vcovHC(icu.lm,type="HC1"))

# Let's look at nonwhite population
nwAvg <- mean(dataComb$Nonwhite); nwAvg
highNW <- subset(dataComb, Pop > nwAvg); highNW
lowNW <- subset(dataComb, Pop <= nwAvg); lowNW
t.test(highNW$HospTotal, lowNW$HospTotal, equal.values = FALSE)
# Significant value of 0.004807 - extremely interesting!

(ggplot(dataComb, aes(x = HospTotal, y = AQI)) + stat_summary_bin(fun.y='mean', bins=25, color='deeppink4', size=1, alpha=0.7, geom='point', na.rm = TRUE)+labs(y= "Average AQI", x = "COVID-19 Hospitalizations")+ggtitle("Binned Scatter Plot of Average AQI vs. COVID Hospitalizations per County"))
(ggplot(dataComb, aes(x = HospTotal, y = Pop)) + stat_summary_bin(fun.y='mean', bins=25, color='deeppink4', size=1, alpha=0.7, geom='point', na.rm = TRUE)+labs(y= "Population", x = "COVID-19 Hospitalizations")+ggtitle("Binned Scatter Plot of Population vs. COVID Hospitalizations per County"))
(ggplot(dataComb, aes(x = HospTotal, y = Nonwhite)) + stat_summary_bin(fun.y='mean', bins=25, color='deeppink4', size=1, alpha=0.7, geom='point', na.rm = TRUE)+labs(y= "Nonwhite Population", x = "COVID-19 Hospitalizations")+ggtitle("Binned Scatter Plot of Nonwhite Population vs. COVID Hospitalizations per County"))
(ggplot(dataComb, aes(x = HospTotal, y = Poverty)) + stat_summary_bin(fun.y='mean', bins=25, color='deeppink4', size=1, alpha=0.7, geom='point', na.rm = TRUE)+labs(y= "Population Living in Poverty", x = "COVID-19 Hospitalizations")+ggtitle("Binned Scatter Plot of Population in Poverty vs. COVID Hospitalizations per County"))
(ggplot(dataComb, aes(x = HospTotal, y = Old)) + stat_summary_bin(fun.y='mean', bins=25, color='deeppink4', size=1, alpha=0.7, geom='point', na.rm = TRUE)+labs(y= "Population Aged 65+", x = "COVID-19 Hospitalizations")+ggtitle("Binned Scatter Plot of Older Population vs. COVID Hospitalizations per County"))

#-------------------------------------------------------------------------------
# Perform a Permutation Test with Air Quality Data from 2016 to 2020
#-------------------------------------------------------------------------------

# Let's run a permutation test where we only have a sample and know nothing
# about the population, looking at a cutoff of counties above the COVID
# hospitalization average and below it (calculated to be 37.872)
# and testing for AQI
covidAvg <- mean(dataComb$HospTotal); covidAvg
below_avg <- sum(dataComb$AQI*(dataComb$HospTotal <= covidAvg))/sum(dataComb$HospTotal <= covidAvg); below_avg
above_avg <- sum(dataComb$AQI*(dataComb$HospTotal > covidAvg))/sum(dataComb$HospTotal > covidAvg); above_avg
observed_air <- below_avg - above_avg; observed_air

above <- subset(dataComb, HospTotal > covidAvg)
below <- subset(dataComb, HospTotal <= covidAvg)
# Now replace above average counties with a random sample
higher <- sample(dataComb$HospTotal); head(higher)

below_avg <- sum(dataComb$AQI*(higher <= covidAvg))/sum(higher <= covidAvg); below_avg
above_avg <- sum(dataComb$AQI*(higher > covidAvg))/sum(higher > covidAvg); above_avg
below_avg - above_avg

# Repeat this 10000 times
N <- 10000
diffs_air <- numeric(N)
for (i in 1:N){
  # Permuted HospTotal column
  higher <- sample(dataComb$HospTotal); higher 
  below_avg <- sum(dataComb$AQI*(higher <= covidAvg))/sum(higher <= covidAvg); below_avg
  above_avg <- sum(dataComb$AQI*(higher > covidAvg))/sum(higher > covidAvg); above_avg
  # Same likelihood to be positive or negative
  diffs_air[i] <- below_avg - above_avg
}

mean(diffs_air) # Theoretically should be close to zero
hist(diffs_air, breaks = "FD", col="purple", probability=TRUE, xlab="Simulated Differences in County Air Quality Indices", ylab="Density", main="Histogram of Permuted County Air Quality Differences in Regards to Above or Below Average COVID-19 Hospitalizations")

# We can display the observed difference on the histogram
# It's not even on the histogram since it is so large
abline(v = observed_air, col = "green")

# Calculate the probability that a difference this large
# could have arose from a random sample
pvalue_air <- (sum(diffs_air >= observed_air)+1)/(N+1); pvalue_air
# P-value is 0.2339766. Since the p=value is greater than 0.05,
# we fail to reject the null hypothesis since there is a relatively
# high chance the permuted tests were as extreme as the observed differences. 
# Therefore, the evidence is not sufficient enough against the null hypothesis 
# that counties with above average or below average COVID hosptilizations
# are statistically significant in regards to air quality. 

# Now we can run a t-test on this data, by running it for air quality
# for counties with a above average and below average COVID hospitalizations
above <- subset(dataComb, HospTotal > covidAvg); above
below <- subset(dataComb, HospTotal <= covidAvg); below
t.test(above$AQI, below$AQI, equal.values = FALSE)
# The built-in t-test for unequal variances gives p-value = 0.4171, 
# which means that our classical method agrees with our simulation method 
# of a permutation test. However, the t-test seems to give a larger 
# p-value in relation to the permutation test, meaning that the simulation
# methods may be more forgiving in regards to significance.

#-------------------------------------------------------------------------------
# Perform a Permutation Test with COVID Data
#-------------------------------------------------------------------------------

# Let's run a permutation test where we only have a sample and know nothing
# about the population, looking at a cutoff of counties above the air quality
# index average and below it (calculated to be 40.26599)
# and testing for hospitalizations
airAvg <- mean(dataComb$AQI); airAvg
below_avg <- sum(dataComb$HospTotal*(dataComb$AQI <= airAvg))/sum(dataComb$AQI <= airAvg); below_avg
above_avg <- sum(dataComb$HospTotal*(dataComb$AQI > airAvg))/sum(dataComb$AQI > airAvg); above_avg
observed_cov <- below_avg - above_avg; observed_cov

# Now replace above average counties with a random sample
higher <- sample(dataComb$AQI); head(higher)

below_avg <- sum(dataComb$HospTotal*(higher <= airAvg))/sum(higher <= airAvg); below_avg
above_avg <- sum(dataComb$HospTotal*(higher > airAvg))/sum(higher > airAvg); above_avg
below_avg - above_avg

# Repeat this 10000 times
N <- 10000
diffs_cov <- numeric(N)
for (i in 1:N){
  # Permuted AQI column
  higher <- sample(dataComb$AQI)
  below_avg <- sum(dataComb$HospTotal*(higher <= airAvg))/sum(higher <= airAvg)
  above_avg <- sum(dataComb$HospTotal*(higher > airAvg))/sum(higher > airAvg)
  # Same likelihood to be positive or negative
  diffs_cov[i] <- below_avg - above_avg
}

mean(diffs_cov) # Theoretically should be close to zero
hist(diffs_cov, breaks = "FD", col="red", probability=TRUE, xlab="Simulated Differences in County COVID Hospitalizations", ylab="Density", main="Histogram of Permuted County COVID Hospitalization Differences in Regards to Above or Below Average Air Quality Indices")

# We can display the observed difference on the histogram
# It's not even on the histogram since it is so large
abline(v = observed_cov, col = "blue")

# Calculate the probability that a difference this large
# could have arose from a random sample
pvalue_cov <- (sum(diffs_cov >= observed_cov)+1)/(N+1); pvalue_cov
# P-value is 0.4221578. Since the p=value is greater than 0.05,
# we fail to reject the null hypothesis since there is a relatively
# high chance the permuted tests were as extreme as the observed differences. 
# Therefore, the evidence is not sufficient enough against the null hypothesis 
# that counties with above average or below average AQI
# are statistically significant in regards to COVID hospitalizations. 

# Now we can run a t-test on this data, by running it for hospitalizations
# for counties with a above average and below average AQI
above <- subset(dataComb, AQI > airAvg); above
below <- subset(dataComb, AQI <= airAvg); below
t.test(above$HospTotal, below$HospTotal, equal.values = FALSE)
# The built-in t-test for unequal variances gives p-value = 0.8271, 
# which means that our classical method agrees with our simulation method 
# of a permutation test.

#-------------------------------------------------------------------------------
# Perform a Permutation Test with ICU Data
#-------------------------------------------------------------------------------

# Let's run a permutation test where we only have a sample and know nothing
# about the population, looking at a cutoff of counties above the air quality
# index average and below it (calculated to be 40.26599)
# and testing for hospitalizations
airAvg <- mean(dataComb$AQI); airAvg
below_avg <- sum(dataComb$ICU*(dataComb$AQI <= airAvg))/sum(dataComb$AQI <= airAvg); below_avg
above_avg <- sum(dataComb$ICU*(dataComb$AQI > airAvg))/sum(dataComb$AQI > airAvg); above_avg
observed_cov <- below_avg - above_avg; observed_cov

# Now replace above average counties with a random sample
higher <- sample(dataComb$AQI); head(higher)

below_avg <- sum(dataComb$ICU*(higher <= airAvg))/sum(higher <= airAvg); below_avg
above_avg <- sum(dataComb$ICU*(higher > airAvg))/sum(higher > airAvg); above_avg
below_avg - above_avg

# Repeat this 10000 times
N <- 10000
diffs_cov <- numeric(N)
for (i in 1:N){
  # Permuted AQI column
  higher <- sample(dataComb$AQI); higher 
  below_avg <- sum(dataComb$ICU*(higher <= airAvg))/sum(higher <= airAvg); below_avg
  above_avg <- sum(dataComb$ICU*(higher > airAvg))/sum(higher > airAvg); above_avg
  # Same likelihood to be positive or negative
  diffs_cov[i] <- below_avg - above_avg
}

mean(diffs_cov) # Theoretically should be close to zero

# Very interesting because it exhibits bimodal behavior
hist(diffs_cov, breaks = "FD", col="orange", probability=TRUE, xlab="Simulated Differences in County ICU Admissions", ylab="Density", main="Histogram of Permuted County ICU Admissions Differences in Regards to Above or Below Average Air Quality Indices")

# We can display the observed difference on the histogram
# It's not even on the histogram since it is so large
abline(v = observed_cov, col = "purple")

# Calculate the probability that a difference this large
# could have arose from a random sample
pvalue_cov <- (sum(diffs_cov >= observed_cov)+1)/(N+1); pvalue_cov
# P-value is 0.529647. Since the p=value is greater than 0.05,
# we fail to reject the null hypothesis since there is a relatively
# high chance the permuted tests were as extreme as the observed differences. 
# Therefore, the evidence is not sufficient enough against the null hypothesis 
# that counties with above average or below average AQI
# are statistically significant in regards to ICU admissions. 

# Now we can run a t-test on this data, by running it for ICU
# for counties with a above average and below average AQI
above <- subset(dataComb, AQI > airAvg); above
below <- subset(dataComb, AQI <= airAvg); below
t.test(above$ICU, below$ICU, equal.values = FALSE)
# The built-in t-test for unequal variances gives p-value = 0.6949, 
# which means that our classical method agrees with our simulation method 
# of a permutation test.

##################################################################
