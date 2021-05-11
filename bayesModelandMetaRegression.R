#Bayesian Modeling and Meta Regression

#We could do with the variables in the same spatial Bayesian model, or we could do purely spatial Bayesian and then a metaregression

#First Separate 

#Bayesian

#d <- aggregate(allRR$RERI, by = list(zipcode = allRR$patzip), FUN = mean)
#names(d) <- c("id", "Y")
#Read in Data file, this will be different depending on the Heat wave definition used

allRR = read.csv("allRR_975_1.csv", header = T)

d <- data.frame(cbind(allRR$ZCTA5CE10, allRR$RERI))
names(d) <- c("id", "Y")

cases <- d$Y
d# <- left_join(d, population, by =c("id" = "patzip"))

d$ZCTA5CE10 <- d$id


#spBayes
#The acceptance rate depends largely on the proposal distribution. If it has small variance, 
#the ratio of the probabilities between the current point and the proposal will necessarily always be close to 1, giving a high acceptance chance. 
#This is just because the target probability densities we typically work with are locally Lipschitz (a type of smoothness) at small scales, 
#so the probability of two nearby points is similar (informally).

#If your current sample is close to the MAP value, the proposals will have less than one acceptance probability, but it can still be very close to 1.

#As a side note, standard practice is to tune the proposal distribution to get around a 0.2-0.25 acceptance rate.
library("spBayes")
library(MBA)
library(geoR)
library(fields)
library(sp)
library(maptools)
library(rgdal)
library(classInt)
library(lattice)
library(gstat) 
coords <- read.csv("~/Downloads/zip_gps.csv", header = T)
coords <- coords[which(coords$ZCTA5 >= 90001),]

coords$ZCTA5CE10 <- coords$ZCTA5
names(d)
bayesDF<- merge(d, coords, by = c("ZCTA5CE10"))

spdf <- SpatialPointsDataFrame(coords = bayesDF[,6:7], data = bayesDF)

v1 = variogram(Y ~ 1, data = spdf)
plot(v1)
#Assumes Isotropy 

# tau sq is nugget 0.7 ish
# sigma sill 1.5
#phi range 300
n.samples = 10000
bef.sp<-spLM(Y ~ 1, data = bayesDF, coords = as.matrix(bayesDF[,c("FinalLat", "FinalLon")]), starting = list("phi" = 3, "sigma.sq" = 1.0, "tau.sq" = 0.6), tuning = list("phi" = 1, "sigma.sq" = 0.8, "tau.sq" = 0.4), priors = list("phi.Unif" = c(0.001, 4), "sigma.sq.IG" = c(2, 2)), cov.model = "spherical", n.samples = n.samples)
round(summary(mcmc(bef.sp$p.theta.samples))$quantiles, 3)
burn.in <- floor(0.75*n.samples)
bef.sp <- spRecover(bef.sp, start = burn.in)

beta.samples <- bef.sp$p.beta.recover.samples
w.samples = bef.sp$p.w.recover.samples

w.hat.mu = apply(w.samples, 1, mean)
w.hat.sd = apply(w.samples, 1, sd)

bayesDF$w_hat_mu = w.hat.mu
bayesDF$w_hat_sd = w.hat.sd

bayesDF = read.csv("BayesDF.csv")
bayesDF = bayesDF[,c(-1)]
#Interpolated Surface

y.residuals = residuals(lm(Y ~ 1, data = d))
par(mfrow = c(1,2))
surf <- mba.surf(cbind(bayesDF[,c(6:7)], y.residuals), no.X = 300, no.Y = 300, extend = FALSE)$xyz.est
z.lim = range(surf[[3]], na.rm = T)
image.plot(surf, xaxs = "r", yaxs = "r", zlim = z.lim, main = "LM.residuals")
par(mfrow = c(1,1))
surf <- mba.surf(cbind(bayesDF[,c(6:7)], bayesDF$w_hat_mu), no.X = 800, no.Y = 800, extend = TRUE)$xyz.est

library(tigris)
library("rgeos")
library(FRK)
MainStates <- map_data("state")
NewStates <- filter(MainStates,region ==  "california" )
NewStates <- df_to_SpatialPolygons(NewStates, "region", c("long", "lat"), CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
#pdf("bayesPlot.pdf")
image.plot(surf,xaxs = "r", yaxs = "r", zlim = c(-3,3.5))
states <- states(cb = TRUE)
caloutline <- states %>% filter_state("California")
#caloutline <- st_as_sf(caloutline)
plot(NewStates, add = T)
cal.cities = us.cities[which(us.cities$country.etc == "CA"),]
cal.cities = cal.cities[which(cal.cities$pop > 400000),]
cal.cities$name = str_sub(cal.cities$name, 1, nchar(cal.cities$name) -3)
points(cal.cities$long, cal.cities$lat, cex = 1, col = "black", pch = 16)
text(cal.cities$long,cal.cities$lat + 0.2, labels=cal.cities$name, cex=0.9, font=2)
#dev.off()

library("magick")
test=image_read("Lag2dayAvgHospitalizationBayesPlot.png")
p1<-image_ggplot(test)





#MEta Regression
#Dataset 
bayesDF$w_hat_mu = w.hat.mu
bayesDF$w_hat_sd = w.hat.sd

#write.csv(bayesDF, "BayesDF.csv")
demo <- read.csv("CA_demographics.csv", header = T)
pollutants <- read.csv("pm_no2_o3_byzcta_04to13.csv", header = T)
pollutants$ZIP <- pollutants$zipcode
new = merge(demo, pollutants, by = c("ZIP"))


#Linear model for Bayesian results from these demographics

metaregression = left_join(bayesDF, new, by = c("ZCTA5CE10" = "ZIP"))
ac <- read.csv("AC_by_zip.csv")
ac$ZCTA5CE10 = ac$servzip
metaregression = left_join(metaregression, ac, by = "ZCTA5CE10")
forplot <- bayesDF
metaregression$AC = metaregression$AC*100
meta_ac <- lm(w_hat_mu~ AC, data = metaregression)
meta_bike <- lm(w_hat_mu ~ median_income, data = metaregression)
metalm <- lm(w_hat_mu ~ pop_dens_sq_mile + perc_female + black_pct + unemployed_pct + median_income +
               gini_index + foreign_born_pct + no_health_ins_pct + perc_over_65 + perc_not_white + urban +
               treecanopy + parkaccess + commute + automobile + ndvi + mean_pm25 + mean_pm10 + mean_no2, data = metaregression)
summary(metalm)

library(olsrr)
#ols_step_best_subset(metalm)

metaselectlm <- lm(w_hat_mu ~ pop_dens_sq_mile + perc_female  + median_income +
                no_health_ins_pct  + perc_not_white + urban +
                commute + automobile + ndvi + mean_pm25 + mean_pm10 + mean_no2, data = metaregression)
summary(metaselectlm)

metaregression$SNR <- metaregression$w_hat_mu/metaregression$w_hat_sd
library(sf)
shape <- st_read(dsn = "~/Downloads/ShapefileZCTA/cb_2016_us_zcta510_500k.shp")
metaregression$truncSNR <- ifelse(metaregression$SNR < 2 & metaregression$SNR > -2, NA, metaregression$SNR )
forplot$ZCTA5CE10 <- forplot$patzip
m <- merge(shape, metaregression, by = ("ZCTA5CE10"))
MainStates <- map_data("state")
NewStates <- filter(MainStates,region ==  "california" )
p <- ggplot()  +
  geom_sf(data = m, aes(fill = truncSNR), colour = NA)  + geom_polygon( data=NewStates, aes(x=long, y=lat, group=group),
                                                                       colour="black", fill = NA)
require(devtools)
p2<- p + scale_fill_gradient2(low = "blue",mid = "white",
                         high="red", midpoint = 0) + theme(legend.key.height = unit(1.0, "cm"), legend.key.width = unit(0.3, "cm")) +
  scale_x_continuous(label = I) +
  guides(fill = guide_colourbar(barheight = unit( 2.25 , "in" ),
                                ticks.colour = "black",
                                ticks.linewidth = 1, frame.colour = "black",
                                frame.linewidth = 1))+
  scale_y_continuous(label = abs) +
  labs(fill = "SNR") + theme( panel.grid = element_blank(),text = element_text(size=10), axis.title.x=element_blank(), axis.title.y = element_blank(), aspect.ratio = 1, axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=1))
#GRey means non-significant

require(gridExtra)
pdf("BayesGrid.pdf")
grid.arrange(p1, p2, ncol = 2)
dev.off()
