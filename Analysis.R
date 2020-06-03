# Script created by Frank D'Agostino 2020
# Analysis conducted for the COVID and Climate
# Project as part of the Coronavirus Visualization Team (CVT)

rm(list = ls())
# Clear console
library(ggplot2)
library(plyr)
library(sandwich)
library(lmtest)

#-------------------------------------------------------------------------------
# Data Cleanup
#-------------------------------------------------------------------------------

# Import the partially pre-processed data
air <- read.csv("countyAirPoll.csv"); head(air)

hosp <- read.csv("countyHosp.csv"); head(hosp)

# For air pollution, we have data spanning back to 2010, but we want
# to focus on data regarding 2019 and 2020, so we take a subset
recent <- subset(air, Year >= 2019); head(recent)

# Clean the data by creating a dataframe of average AQI values for
# graphing and easy visualization
berk <- subset(recent, County == "Berkshire"); brist <- subset(recent, County == "Bristol")
dukes <- subset(recent, County == "Dukes"); essex <- subset(recent, County == "Essex")
franklin <- subset(recent, County == "Franklin"); hampden <- subset(recent, County == "Hampden")
hamp <- subset(recent, County == "Hampshire"); middle <- subset(recent, County == "Middlesex")
norf <- subset(recent, County == "Norfolk"); plym <- subset(recent, County == "Plymouth")
suff <- subset(recent, County == "Suffolk"); worc <- subset(recent, County == "Worchester")

dataMean <- data.frame(Countys = c("Berkshire", "Bristol", "Dukes", "Essex",
                               "Franklin", "Hampden", "Hampshire", "Middlesex",
                               "Norfolk", "Plymouth", "Suffolk", "Worchester"),
                   AQI = c(mean(berk$AQI), mean(brist$AQI), mean(dukes$AQI),
                           mean(essex$AQI), mean(franklin$AQI), mean(hampden$AQI),
                           mean(hamp$AQI), mean(middle$AQI), mean(norf$AQI),
                           mean(plym$AQI), mean(suff$AQI), mean(worc$AQI))); head(data)

#-------------------------------------------------------------------------------
# Basic Visualization
#-------------------------------------------------------------------------------

# Use ggplot to visually demonstrate the relatively small differences
# in AQI (stands for air quality index)
ggplot(dataMean, aes(x=Countys, y=AQI, fill=Countys))+
  geom_bar(stat = "identity", colour="black", width=0.5)+labs(y= "Average AQI", x = "County")+
  ggtitle("Barplot of Average AQI Value at Each Massachusetts County in 2019-2020")+
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
# of April and May from 2020, since the hospitalization data spans from
# 4/30/20 to 5/22/20, which is the most recent date
count <- recent$County[grep("^([5][\\/]|[4][\\/])\\d+(\\/20)", recent$Date)]; head(count)
poll <- recent$AQI[grep("^([5][\\/]|[4][\\/])\\d+(\\/20)", recent$Date)]; head(poll)
mayAir <- data.frame(count, poll); head(mayAir)

# Use the aggregate function to group everything by counties
temp <- aggregate(x=hosp, by=list(hosp$County), FUN=mean); temp
aggHosp <- subset(temp, Group.1 != "Barnstable" & Group.1 != "Nantucket"); aggHosp

aggAir <- aggregate(x=mayAir, by=list(mayAir$count), FUN=mean); aggAir

# Nantucket and Barnstable do not have enough Air data in April/May 2020,
# and so we remove them from the hospital data to have even comparisons
dataComb <- data.frame(County = aggAir$Group.1,
                       AQI = aggAir$poll, HospTotal = aggHosp$HospTotal,
                       ICU = aggHosp$ICU); dataComb

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

#-------------------------------------------------------------------------------
# Perform a Permutation Test
#-------------------------------------------------------------------------------

# Let's run a permutation test where we only have a sample and know nothing
# about the population, looking at a cutoff of counties above the COVID
# hospitalization average and below it (calculated to be 34.74136)
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
  # Permuted AQI column
  higher <- sample(dataComb$HospTotal); higher 
  below_avg <- sum(dataComb$AQI*(higher <= covidAvg))/sum(higher <= covidAvg); below_avg
  above_avg <- sum(dataComb$AQI*(higher > covidAvg))/sum(higher > covidAvg); above_avg
  # Same likelihood to be positive or negative
  diffs_air[i] <- below_avg - above_avg
}

mean(diffs_air) # Theoretically should be close to zero
hist(diffs_air, breaks = "FD", col="purple", probability=TRUE, xlab="Simulated Differences in County COVID Hospitalizations", ylab="Density", main="Histogram of Permuted County Air Quality Differences in Regards to Above or Below Average COVID-19 Hospitalizations")

# We can display the observed difference on the histogram
# It's not even on the histogram since it is so large
abline(v = observed_air, col = "green")

# Calculate the probability that a difference this large
# could have arose from a random sample
pvalue_air <- (sum(diffs_air >= observed_air)+1)/(N+1); pvalue_air
# P-value is 0.2584742. Since the p=value is greater than 0.05,
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
# The built-in t-test for unequal variances gives p-value = 0.4492, 
# which means that our classical method agrees with our simulation method 
# of a permutation test. However, the t-test seems to give a much larger 
# p-value in relation to the permutation test, meaning that the simulation
# methods may be more forgiving in regards to significance.

##################################################################
