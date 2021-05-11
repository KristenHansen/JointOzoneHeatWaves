setwd("~/Library/Mobile Documents/com~apple~CloudDocs/OzoneHeatJoint")
source("ozoneheatwavedefs.R")

dataPrep <- function(){
  print("Loading Datasets")
  all <- read.csv("Whole\ CA\ Respi.csv", header = T)
  summary(all$date)
  #popZCTA <- read.csv("2010Census_DemoProfile.csv", header = TRUE)
  all$date2 <- as.Date(all$date, format = "%d%b%Y")
  all.days <- all[which(all$date2 > as.Date("2004-01-01")),]
  # ozone <- read.csv("ca_pollutants_90t17.csv", header = T)
  # ozone$date2 <- as.Date(ozone$date, format = "%d%b%Y")
  # ozone <- ozone[which(ozone$date2 > as.Date("2004-01-01")),]
  ozone = read.csv("Ozone_2004_13.csv")
  ozone$date2 <- as.Date(ozone$date, format = "%Y-%m-%d")
  #ozone <- ozone[which(ozone$date2 > as.Date("2004-01-01")),]
  ozone$max_o3 = ozone$ozone_idw
  ozone$zipcode = ozone$zcta5
  
  temps <- read.csv("ca_zip_temp_04t13.csv", header = TRUE)
  
  temps$date2 <- as.Date(temps$date, format = "%Y-%m-%d")
  temps2 <- temps[which(temps$date2 > as.Date("2004-01-01")),]
  
  temps2$zipcode <- temps2$zipnum
  temps2$zipnum = NULL
  library(tidyr)
  library(dplyr)
  library(ggplot2)
  library(lme4)
  
  all.days <- inner_join(all.days, temps2, by = c("patzip" = "zipcode", "date2"))
  #summary(ozone$date)
  
  
  #data(zipcode)
  #zipcode$zip = readr::parse_number(zipcode$zip)
  
  testset <- inner_join(all.days, ozone, by = c("patzip" = "zipcode", "date2"))
  
  testset$weekday <- weekdays(testset$date2)
  
  testset$month <- format(testset$date2,"%B")
  newtest <- subset(testset, testset$month == "May" | testset$month == "June" | testset$month == "July" | testset$month == "August" | testset$month == "September")
  newtest$month <- as.factor(newtest$month)
  newtest$weekday <- as.factor(newtest$weekday)
  newtest$year <- format(newtest$date2, "%Y")
  newtest <- within(newtest, month <- relevel(month, ref = "May"))
  newtest <- within(newtest, weekday <- relevel(weekday, ref = "Monday"))
  library(dplyr)
  library(Hmisc)
  newtest <- 
    newtest %>%
    group_by(patzip) %>%
    mutate(l1resp = dplyr::lead(resp_disease, n = 1, default = NA))
  
  newtest <- 
    newtest %>%
    group_by(patzip) %>%
    mutate(l2resp = dplyr::lead(resp_disease, n = 2, default = NA))
  
  
  newtest$resp_Avglead2 = rowMedians(as.matrix(newtest[, c("l1resp", "l2resp")]))
  newtest$resp_Avglead2 = ifelse(is.na(newtest$resp_Avglead2), newtest$resp_disease, newtest$resp_Avglead2)
  
  #data(zipcode)
  #zipcode$zip = readr::parse_number(zipcode$zip)
  
  #newtest <- left_join(newtest, zipcode, by = c("patzip" = "zip"))
  return(newtest)
}
newtest <- dataPrep()
newtest <- ozoneDefs(newtest)
newtest <- heatwaveDefs(newtest)

library(dplyr)
data.frame(table(newtest$patzip)) 
#1530 days for each zip code
test_missingness=newtest %>% 
  group_by(patzip) %>% 
  summarise_all(~sum(is.na(max_o3))) %>% 
  transmute(patzip, sumNA = rowSums(.[-1]))

# jointOHWave <- function(x, ozoneDef, heatwaveDef){
#   x$jointOH <- ifelse(ozoneDef == 1 & heatwaveDef == 1, 1, ifelse(heatwaveDef == 1 | ozoneDef == 1, -1, 0))
#   x$onlyHW <- ifelse(ozoneDef == 0 & heatwaveDef == 1, 1, 0)
#   x$onlyOz <- ifelse(ozoneDef == 1 & heatwaveDef == 0, 1, 0)
#   return(x)
# }
# newtest <- jointOHWave(newtest, newtest$o75wave, newtest$hw99_1)
#newtest$max_o3 <- 1000*newtest$max_o3
missingness_byzip = numeric()
j = 1
for(i in unique(newtest$patzip)){
  temp = newtest[which(newtest$patzip == i),]
  missingness_byzip[j]=sum(is.na(temp$max_o3))/nrow(temp)
  j = j + 1
}
#h = missingness_byzip[which(missingness_byzip != 1)]
# Inverse distance in time metric for creating matching 

# Take all the data from the zip code and same year with exposed days

# Use date subtraction in R to compute "distance" between two dates (to each case day) 
#1/distance is your weight


# Take the hospitalizations on those dates and create a weighted average, ie you multiply the hospitalizations
#by the weight and then divide by the sum of the weights 

KernelDistanceMatching <- function(case, control_dataset){
  weightedAvgDiff = numeric(0)
  print(nrow(case))
  for(i in 1:nrow(case)){
    if(i %% 1000 == 0) print(i)
    controls <- control_dataset[which(control_dataset$patzip == case$patzip[i] & control_dataset$year == case$year[i]),]
    controls$distance = as.numeric(abs(case$date2[i] - controls$date2), units = "days")
    controls$weight = 1/controls$distance
    weightedAvgDiff[i] = case$resp_Avglead2[i] - sum(controls$resp_Avglead2*controls$weight)/sum(controls$weight)
  }
  return(weightedAvgDiff)
}



MatchProcedure <- function(x){
  dataMatch <- x
  dataMatch$month <- as.numeric(dataMatch$month)
  
  Xformatch <- dataMatch[,c("patzip","year","date2","hw95_1", "o75wave", "resp_disease", "resp_Avglead2", "latitude", "longitude")]
 
  #Get Contrast for Joint Days
  XJoint <- Xformatch[which(Xformatch$hw95_1 == 1 & Xformatch$o75wave == 1),]
  print(nrow(XJoint))

  zipcodesJoint <- unique(XJoint$patzip)
  XNon <- Xformatch[which(Xformatch$hw95_1 == 0 & Xformatch$o75wave == 0),]
  
  
  zipcodesTotal <- unique(XNon$patzip)
  
 
  
  XJoint$diffMean <- KernelDistanceMatching(XJoint, XNon)
  
  #Get contrast for heat wave only days
  XHeat <- Xformatch[which(Xformatch$hw95_1 == 1 & Xformatch$o75wave == 0),]
  XHeat$diffMean <- KernelDistanceMatching(XHeat, XNon)

  
  #Get contrast for ozone wave only days
  XOzone <- Xformatch[which(Xformatch$hw95_1 == 0 & Xformatch$o75wave == 1),]
  XOzone$diffMean <- KernelDistanceMatching(XOzone,XNon)
  
  return(list(Joint = XJoint, HW = XHeat, OZ = XOzone))
}
  

matched = MatchProcedure(newtest)
#write.csv(matched$Joint, "joint975_2.csv")
#write.csv(matched$HW, "hwonly975_2.csv")
#write.csv(matched$OZ, "ozonly975_2.csv")

matched$Joint$RR <- matched$Joint$resp_disease/(matched$Joint$resp_disease - matched$Joint$diffMean)
matched$HW$RR <- matched$HW$resp_disease/(matched$HW$resp_disease - matched$HW$diffMean) 
matched$OZ$RR <- matched$OZ$resp_disease/(matched$OZ$resp_disease - matched$OZ$diffMean)


# matched$Joint$RR <- matched$Joint$resp_Avglead2/(matched$Joint$resp_Avglead2 - matched$Joint$diffMean)
# matched$HW$RR <- matched$HW$resp_Avglead2/(matched$HW$resp_Avglead2 - matched$HW$diffMean) 
# matched$OZ$RR <- matched$OZ$resp_Avglead2/(matched$OZ$resp_Avglead2 - matched$OZ$diffMean)

matched$Joint <- matched$Joint[which(is.finite(matched$Joint$RR)),]

matched$HW <- matched$HW[which(is.finite(matched$HW$RR)),]

matched$OZ <- matched$OZ[which(is.finite(matched$OZ$RR)),]



RRjoint = aggregate(matched$Joint$RR, list(matched$Joint$patzip), mean)
RRjoint$jointRR = RRjoint$x 
RRjoint$x = NULL

RRhw =aggregate(matched$HW$RR, list(matched$HW$patzip), mean)
RRhw$hwRR = RRhw$x 
RRhw$x = NULL


RRoz = aggregate(matched$OZ$RR, list(matched$OZ$patzip), mean)
RRoz$ozRR = RRoz$x 
RRoz$x = NULL

allRR = inner_join(RRjoint, RRhw, by = c("Group.1"))
allRR = inner_join(allRR, RRoz, by = c("Group.1"))

allRR$RERI = allRR$jointRR - allRR$hwRR - allRR$ozRR + 1
allRR$RERI <- if_else(allRR$RERI < -3,-3,allRR$RERI)
allRR$RERI <- if_else(allRR$RERI > 3, 3, allRR$RERI)

allRR$ZCTA5CE10 <- allRR$Group.1

# write.csv(allRR, "allRR_975_2.csv")
# library(sf)
missing =data.frame(cbind(unique(newtest$patzip), missingness_byzip))
names(missing) <- c("ZCTA5CE10", "missing")
shape <- st_read(dsn = "~/Downloads/ShapefileZCTA/cb_2016_us_zcta510_500k.shp")
m <- merge(shape, missing, by = ("ZCTA5CE10"))
# 
# allRR$hwRR = if_else(allRR$hwRR > 6, 6, allRR$hwRR)
# allRR$ozRR = if_else(allRR$ozRR>6, 6, allRR$ozRR)
# nonEnviro = as.data.frame(matrix(nrow = length(setdiff(unique(newtest$patzip), RRoz$Group.1)),ncol = ncol(allRR)))
# names(nonEnviro) = names(allRR)
# nonEnviro$ZCTA5CE10 = setdiff(unique(newtest$patzip), RRoz$Group.1)
# 
# allRR = rbind(nonEnviro,allRR)
# m <- merge(shape, allRR, by = ("ZCTA5CE10"))
# par(mfrow = c(1,2))
MainStates <- map_data("state")
NewStates <- filter(MainStates,region ==  "california" )
# #hw RR
plot1 <- ggplot()  +
  geom_sf(data = m, aes(fill = missing), colour = NA)  + geom_polygon( data=NewStates, aes(x=long, y=lat, group=group),
                                                                    colour="black", fill = NA)
#
plot1 + scale_fill_gradient2(low = "white",high="red", limits = c(0,0.99)) + theme(legend.key.size = unit(1, "cm")) +labs(fill = "missingness") + ggtitle("Missingness in O3 estimates") + theme(panel.grid = element_blank())
# #oz RR
# plot2 <- ggplot()  +
#   geom_sf(data = m, aes(fill = ozRR), colour = NA) + geom_polygon( data=NewStates, aes(x=long, y=lat, group=group),
#                                                                    colour="black", fill = NA) 
# 
# plot2 <- plot2 + scale_fill_gradient2(low = "blue",mid = "white",
#                          high="red", midpoint = 1, limits = c(0,6)) + theme(legend.key.size = unit(1, "cm")) +
#   labs(fill = "RR") + ggtitle("Ozone relative risk") + theme(panel.grid = element_blank())
# 
# require(gridExtra)
# pdf("RRplotgrid.pdf")
# grid.arrange(plot1, plot2, ncol = 2)
# dev.off()

  