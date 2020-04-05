
#Matching Procedure function
runMatchProc <- function(x){
  
  
  
  dataMatch <- x
  dataMatch$month <- as.numeric(dataMatch$month)
  
  Xformatch <- dataMatch[,c("patzip", "month", "weekend", "year","hw95_1", "o75wave", "resp_disease", "jointOH", "latitude", "longitude", "onlyHW", "onlyOz")]
  
  # line 87 must be changed for the three runs such that we are using cases that are 1/1,1/0 and 0/1
  #Three runs are in 117-119
  XHeat <- Xformatch[which(Xformatch$hw95_1 == 0 & Xformatch$o75wave == 1),]
  print(nrow(XHeat))
  #XHeat <- XHeat[,c("patzip", "month", "weekday")]
  zipcodesJoint <- unique(XHeat$patzip)
  XNon <- Xformatch[which(Xformatch$hw95_1 == 0 & Xformatch$o75wave == 0),]
  #XNon <- XNon[,c("patzip", "month", "weekday")]
  zipcodesTotal <- unique(XNon$patzip)
  
  XHeat$meanNon = XHeat$medianNon = numeric(nrow(XHeat))
  #m.out <- matchit(heatwave ~ patzip + month + weekday, data = Xformatch, method = "exact", ratio = 4)
  #j = 1
  #matchedData <- data.frame(matrix(ncol=11 , nrow= 3*nrow(XHeat)))
  #names(matchedData) = c("patzip", "month", "weekday","hw95_1", "o75wave","resp_disease", "jointOH", "latitude", "longitude", "onlyHW", "onlyOz")
  for(i in 1:nrow(XHeat)){
    if (i %% 1000 == 0) print(i)
    potentialMatch <- XNon[which(XNon$patzip == XHeat$patzip[i]),]
    indices <- ifelse(potentialMatch$month >= XHeat$month[i] - 1 &potentialMatch$month <= XHeat$month[i] + 1 & potentialMatch$weekend == XHeat$weekend[i] & potentialMatch$year == XHeat$year[i], 1,0)
    #print(sum(indices))
    xMatched <- potentialMatch[which(indices == 1),]
    #xMatched <- xMatched[sample(nrow(xMatched), 3, replace = TRUE),]
    XHeat$meanNon[i] <- sum(xMatched$resp_disease)/nrow(xMatched)
    XHeat$medianNon[i] <- median(xMatched$resp_disease)
  }
  
  
  XHeat$diffMed <- XHeat$resp_disease - XHeat$medianNon
  XHeat$diffMean <- XHeat$resp_disease - XHeat$meanNon
  return(list(XHeat, zipcodesJoint, zipcodesTotal))
}

matchedDataset <- runMatchProc(newtest) #CChange 87 to 1/1
matchedDatasetHW <- runMatchProc(newtest) # Change to 1/0
matchedDatasetOZ <- runMatchProc(newtest) #Change to 0/1

matchedDataset[[1]]$RR <- matchedDataset[[1]]$resp_disease/matchedDataset[[1]]$meanNon
matchedDatasetHW[[1]]$RR <- matchedDatasetHW[[1]]$resp_disease/matchedDatasetHW[[1]]$meanNon
matchedDatasetOZ[[1]]$RR <- matchedDatasetOZ[[1]]$resp_disease/matchedDatasetOZ[[1]]$meanNon

matchedDataset <- matchedDataset[[1]][complete.cases(matchedDataset[[1]]),]
matchedDatasetHW <- matchedDatasetHW[[1]][complete.cases(matchedDatasetHW[[1]]),]
matchedDatasetOZ <- matchedDatasetOZ[[1]][complete.cases(matchedDatasetOZ[[1]]),]

matchedDataset <- matchedDataset[which(is.finite(matchedDataset$RR)),]

matchedDatasetHW <- matchedDatasetHW[which(is.finite(matchedDatasetHW$RR)),]
#Not in HW but in other 2 #90822 92055 93429 93453 94020 95064 95075 95234 95680 95717 95728 95934 95970
matchedDatasetOZ <- matchedDatasetOZ[which(is.finite(matchedDatasetOZ$RR)),]

#Not in Joint but in OZ #91371 93546 93603 94575 95536 95537 95836 95983 96120 96146 96148
#90822 92055 93429 93453 94020 95064 95075 95234 95680 95717 95728 95934 95970 91371 93546 93603 94575 95536 95537 95836 95983 96120 96146 96148

geoNeed <- unique(matchedDataset[,c("latitude", "longitude", "patzip")])





RRjoint = aggregate(matchedDataset[,17], list(matchedDataset$patzip), mean)
RRjoint$jointRR = RRjoint$x 
RRjoint$x = NULL

RRhw =aggregate(matchedDatasetHW[,17], list(matchedDatasetHW$patzip), mean)
RRhw$hwRR = RRhw$x 
RRhw$x = NULL


RRoz = aggregate(matchedDatasetOZ[,17], list(matchedDatasetOZ$patzip), mean)
RRoz$ozRR = RRoz$x 
RRoz$x = NULL

allRR = inner_join(RRjoint, RRhw, by = c("Group.1"))
allRR = inner_join(allRR, RRoz, by = c("Group.1"))

allRR$RERI = allRR$jointRR- allRR$hwRR - allRR$ozRR + 1
allRR$RERI <- if_else(allRR$RERI < -3,-3,allRR$RERI)
allRR$RERI <- if_else(allRR$RERI > 3, 3, allRR$RERI)

allRR$signRERI <- if_else(allRR$RERI > 0, "positive", if_else(allRR$RERI == 0, "zero", "negative"))


popZCTA <- read.csv("2010Census_DemoProfile.csv", header = TRUE)

library(tigris)
library("rgeos")
caloutline <- states() %>% filter_state("California")
caloutline <- st_as_sf(caloutline)
allRR$ZCTA5CE10 <- allRR$Group.1
library(sf)
#shape <- st_read(dsn = "~/Downloads/ShapefileZCTA/cb_2016_us_zcta510_500k.shp")
diffjoint = setdiff(unique(matchedDatasetOZ$patzip),unique(matchedDataset$patzip))
diffjoint = as.data.frame(diffjoint)
names(diffjoint) = "ZCTA5CE10"


diffhw = setdiff(unique(matchedDataset$patzip),unique(matchedDatasetHW$patzip))
diffhw = as.data.frame(diffhw)
names(diffhw) = "ZCTA5CE10"
temp <- forPlot[which(is.na(forPlot$max_o3)),]
withoutOZ <- as.data.frame(temp$ZCTA5CE10)
names(withoutOZ) <- "ZCTA5CE10"  
withoutOZ$Group.1 <- withoutOZ$ZCTA5CE10
withoutOZ$jointRR <- rep(NA, nrow(withoutOZ))
withoutOZ$hwRR <- rep(NA, nrow(withoutOZ))
withoutOZ$ozRR <- rep(NA, nrow(withoutOZ))
withoutOZ$RERI <- rep(NA, nrow(withoutOZ))
#withoutOZ$RERI <- rep(NA, nrow(withoutOZ))

library(plyr)
fulldat <- join(allRR,withoutOZ, by = ("ZCTA5CE10"), type = "full")


shape <- st_read(dsn = "~/Downloads/ShapefileZCTA/cb_2016_us_zcta510_500k.shp")
m <- merge(shape, fulldat, by = ("ZCTA5CE10"))

nojoint <- merge(shape,diffjoint, by = ("ZCTA5CE10"))
nohw <- merge(shape, diffhw, by = ("ZCTA5CE10"))
p <- ggplot()  +
  geom_sf(data = m, aes(fill = RERI), colour = NA) +geom_sf(data = nojoint, fill = "grey5", colour = NA) + geom_sf(data = nohw, fill = "grey25", colour = NA)

p + scale_fill_gradient2(low = "purple",mid = "white",
                         high="orange", midpoint = 0, limits = c(-3,3)) + theme(legend.key.size = unit(1, "cm")) +
  labs(fill = "RERI")



library(gridExtra)
#HW RR plot
p <- ggplot(data = m)  + 
  geom_sf(aes(fill = m$hwRR), colour = NA) 

plot1<- p + scale_fill_gradient2(low = "blue",mid = "white",
                                 high="red", midpoint = 1, limits = c(0,3)) + theme(legend.key.size = unit(1, "cm"))+
  labs(fill = "RR (heat 95 %ile)") 

#OZ RR plot
p <- ggplot(data = m)  + 
  geom_sf(aes(fill = m$ozRR), colour = NA) 

plot2<-p + scale_fill_gradient2(low = "blue",mid = "white",
                                high="red", midpoint = 1, limits = c(0,3)) + theme(legend.key.size = unit(1, "cm"))+
  labs(fill = "RR (ozone 75 %ile)")

grid.arrange(plot1, plot2, ncol = 2)

#Joint RR plot
p <- ggplot(data = m) + geom_sf(data = caloutline, fill = "grey", colour = "black") + 
  geom_sf(aes(fill = m$jointRR), colour = NA) 

p + scale_fill_gradient2(low = "blue",mid = "white",
                         high="red", midpoint = 1, limits = c(0,3)) + theme(legend.key.size = unit(1, "cm"))+
  labs(fill = "RR") +
  ggtitle("RR of Resp hosp on joint HW and ozone days")


m <- merge(shape, allRR, by = ("ZCTA5CE10"))
p <- ggplot(data = m) + geom_sf(aes(fill = m$signRERI), colour = NA)
p + theme(legend.key.size = unit(2, "cm"))+
  labs(fill = "RERI") +
  ggtitle("RERI sign for joint Ozone (75%) and Heatwave (95%-1day) days")


diffByZipHeatWave<- aggregate(matchedDataset[,14], list(matchedDataset$patzip), mean)
diffByZipHeatWave[,2] <- ifelse(diffByZipHeatWave[,2] < 0, 0 , diffByZipHeatWave[,2])
diffByZipHeatWave$ZCTA5CE10 = diffByZipHeatWave[,1]
diffByZipHeatWave$latitude = geoNeed$latitude
diffByZipHeatWave$longitude = geoNeed$longitude
#diffByZipHeatWave$population = geoNeed$Total.population
summary(diffByZipHeatWave[,2])
#Medians with mean for zipcode
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-1.98718 -0.72727  0.02564  0.14099  0.64935  5.37662 

#Means
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#-2.1226 -1.0294 -0.7212 -0.3550  0.1714  4.6326 

# library(sf)
# shape <- st_read(dsn = "~/Downloads/ShapefileZCTA/cb_2016_us_zcta510_500k.shp")
# m <- merge(shape, diffByZipHeatWave, by = ("ZCTA5CE10"))
# p <- ggplot(data = m) + geom_sf(aes(fill = m$x), colour = NA)
# p + scale_fill_gradient2(low = "blue",mid = "white",
#                          high="red") + theme(legend.key.size = unit(2, "cm"))
# 
# 
# forPlot <- newtest[!duplicated(newtest[,c('patzip')]),]
# forPlot$ZCTA5CE10 <- forPlot$patzip
# forPlot$pct75_maxo3 <- 1000*forPlot$pct75_maxo3
# forPlot$pct95_maxo3 <- 1000*forPlot$pct95_maxo3
# forPlot$patzip = NULL
# library(sf)
# shape <- st_read(dsn = "~/Downloads/ShapefileZCTA\ 2/cb_2016_us_zcta510_500k.shp")
# m <- merge(shape, forPlot, by = ("ZCTA5CE10"))
# p <- ggplot(data = m) + geom_sf(aes(fill = m$pct95_maxo3), colour = NA)
# p + scale_fill_gradient2(low = "blue",mid = "white",
#                          high="red", midpoint = 60) + theme(legend.key.size = unit(1, "cm"))

# edf. <- as.data.frame(edf.)
# edf.$ZCTA5CE10 <- row.names(edf.)
# edf.$x <- edf.[,1]
# library(sf)
# shape <- st_read(dsn = "~/Downloads/ShapefileZCTA\ 2/cb_2016_us_zcta510_500k.shp")
# m <- merge(shape, edf., by = ("ZCTA5CE10"))
# p <- ggplot(data = m) + geom_sf(aes(fill = m$x), colour = NA)
# p + scale_fill_gradient2(low = "blue", mid = "white",
#                          high="red", midpoint = 0.7) + theme(legend.key.size = unit(1, "cm"))



j = 1
stderrs = numeric(0)
patzip = numeric(0)
for(i in unique(matchedDataset$patzip)){
  print(i)
  j = j + 1
  library(boot)
  
  x <- matchedDataset[which(matchedDataset$patzip == i),]
  if(is.nan(x$RR)){
    stderrs[j] = 0
    patzip[j] = i
  }else{
    meanFunc <- function(x,j){mean(x[j])}
    bootMean <- boot(x$RR[is.finite(x$RR)],meanFunc,1000)
    stderrs[j] <- sd(bootMean$t)
    patzip[j] <- i
  }
}

j = 1
stderrshw = numeric(0)
patziphw = numeric(0)
for(i in unique(matchedDatasetHW$patzip)){
  print(i)
  j = j + 1
  library(boot)
  
  x <- matchedDatasetHW[which(matchedDatasetHW$patzip == i),]
  if(is.nan(x$RR)){
    stderrshw[j] = 0
    patziphw[j] = i
  }else{
    meanFunc <- function(x,j){mean(x[j])}
    bootMean <- boot(x$RR[is.finite(x$RR)],meanFunc,1000)
    stderrshw[j] <- sd(bootMean$t)
    patziphw[j] <- i
  }
}


j = 1
stderrsoz = numeric(0)
patzipoz = numeric(0)
for(i in unique(matchedDatasetOZ$patzip)){
  print(i)
  j = j + 1
  library(boot)
  
  x <- matchedDatasetOZ[which(matchedDatasetOZ$patzip == i),]
  if(is.nan(x$RR)){
    stderrsoz[j] = 0
    patzipoz[j] = i
  }else{
    meanFunc <- function(x,j){mean(x[j])}
    bootMean <- boot(x$RR[is.finite(x$RR)],meanFunc,1000)
    stderrsoz[j] <- sd(bootMean$t)
    patzipoz[j] <- i
  }
}

jointSD <- as.data.frame(cbind(stderrs, patzip))
hwSD <- as.data.frame(cbind(stderrshw, patziphw))
ozSD <- as.data.frame(cbind(stderrsoz, patzipoz))

joinedSD <- inner_join(jointSD, hwSD, by = c("patzip" = "patziphw"))
joinedSD <- inner_join(joinedSD, ozSD, by = c("patzip" = "patzipoz"))

joinedSD$totSD <- sqrt(joinedSD$stderrs^2 + joinedSD$stderrshw^2 + joinedSD$stderrsoz^2)
#total$ZCTA5CE10 <- total$Group.1

total <- inner_join(allRR, joinedSD, by = c("Group.1" = "patzip"))
total$low.ci = total$RERI - 2*(total$totSD)
total$hi.ci = total$RERI + 2*(total$totSD)

summary(total)