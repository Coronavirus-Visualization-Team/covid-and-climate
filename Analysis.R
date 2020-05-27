# Script created by Frank D'Agostino 2020
# Analysis conducted for the COVID and Climate
# Project as part of the Coronavirus Visualization Team (CVT)

rm(list = ls())
# Clear console
library(ggplot2)
library(plyr)

#-------------------------------------------------------------------------------
# Analysis for COVID-19 and Climate
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

# Use ggplot to visually demonstrate the relatively small differences
# in AQI (stands for air quality index)
ggplot(dataMean, aes(x=Countys, y=AQI, fill=Countys))+
  geom_bar(stat = "identity", colour="black", width=0.5)+labs(y= "Average AQI", x = "County")+
  ggtitle("Barplot of Average AQI Value at Each Massachusetts County in 2019-2020")+
  theme_classic()+
  scale_fill_brewer(palette="Set3")

# Produce a barplot demonstrating the number of COVID-19 hospitalizations
# per county using the builtin R barplot function
barplot(table(hosp$HospTotal, hosp$County), col = "lightblue", xlab = "Massachusetts County", ylab = "Total Hospitalizations", main="Total Number of COVID-19 Hospitalizations in Each Massachusetts County")

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


##################################################################
