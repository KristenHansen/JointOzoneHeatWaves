bayesDF = read.csv("BayesDF.csv", header = T)

demo <- read.csv("CA_demographics.csv", header = T)
pollutants <- read.csv("pm_no2_o3_byzcta_04to13.csv", header = T)
pollutants$ZIP <- pollutants$zipcode
new = merge(demo, pollutants, by = c("ZIP"))

AC = read.csv("AC_by_zip.csv")
AC$ZIP = AC$servzip
new = merge(new, AC, by = c("ZIP"))
#Linear model for Bayesian results from these demographics

metaregression = left_join(bayesDF, new, by = c("ZCTA5CE10" = "ZIP"))

metalm <- lm(w_hat_mu ~ pop_dens_sq_mile + perc_female + black_pct + unemployed_pct + median_income +
               gini_index + foreign_born_pct + no_health_ins_pct + perc_over_65 + perc_not_white + urban +
               treecanopy + parkaccess + commute + automobile + ndvi + mean_pm25 + mean_pm10 + mean_no2, data = metaregression)
summary(metalm)

meta_AC <- lm(w_hat_mu ~ mean_pm25 + AC, data = metaregression)
0.031418 -1.96*0.011574
meta_AC_2 <- lm(w_hat_mu ~ mean_pm10 + AC, data = metaregression)
0.0133 -1.96*0.0048
meta_AC_3 <- lm(w_hat_mu ~ mean_no2 + AC, data = metaregression)
0.019 - 1.96*0.006

meta_perc_nonwhite = lm(w_hat_mu ~ perc_not_white + median_income, data = metaregression)
summary(meta_perc_nonwhite)
Est_nonwhite = 2.828e-03
CI_nonwhite = c(2.828e-03 - 1.96*2*1.598e-03, 2.828e-03 + 1.96*2*1.598e-03)

meta_perc_black = lm(w_hat_mu ~ black_pct + median_income, data = metaregression)
summary(meta_perc_black)
Est_black = 2.786e-03
CI_black = c(2.786e-03 - 1.96*3.98e-03, 2.786e-03 + 1.96*3.98e-03)

meta_commute = lm(w_hat_mu ~ commute + median_income, data = metaregression)
summary(meta_commute)
Est_commute = -9.822e-03
CI_commute = c(-9.822e-03 - 1.96*2.504e-03, -9.822e-03 + 1.96*2.504e-03)

correlation.matrix=cor(na.omit(demo))
#install.packages("corrplot")
library(corrplot)
corrplot(correlation.matrix, method="color")

#survdata <- read.csv("Public2009Survdata.csv", header = T) #servzip

PCA = prcomp(~ pop_dens_sq_mile + perc_female + black_pct + unemployed_pct + median_income +
         gini_index + foreign_born_pct + no_health_ins_pct + perc_over_65 + perc_not_white + urban +
         treecanopy + parkaccess + commute + automobile + ndvi + mean_pm25 + mean_pm10 + mean_no2 + AC, data = metaregression)
library(ggfortify)
pca.plot <- autoplot(PCA, data = metaregression, colour = 'median_income')
pca.plot

correlation <- cor(metaregression[,c(4,13:36)],use = "complete.obs")
r.eigen <- eigen(correlation)
for (r in r.eigen$values) {
  print(r / sum(r.eigen$values))
}

PCA = princomp(~ tot_pop+pop_dens_sq_mile + perc_female + black_pct + unemployed_pct + median_income +
               gini_index + foreign_born_pct + no_health_ins_pct + perc_over_65 + perc_not_white + urban +
               treecanopy + parkaccess + commute + automobile + ndvi + mean_pm25 + mean_pm10 + mean_no2 + AC,scores = T,cor = T, data = metaregression)
summary(PCA)
#metaregression$pca_composite=PCA$scores[,1]
meta_subset = metaregression[,c(11, 13:29,31:34,36)]
meta_complete=meta_subset[complete.cases(meta_subset),]
meta_complete$pca_composite = (PCA$scores[,1] - mean(PCA$scores[,1]))/sd(PCA$scores[,1])
meta_complete$pca_2 = (PCA$scores[,2] - mean(PCA$scores[,2]))/sd(PCA$scores[,2]) #also important ?? barely smaller in magnitude
meta_comp =lm(w_hat_mu ~ pca_composite, data = meta_complete)

#Create Composite using PCA of all variables 
#Use LM, it is significant

# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)    0.05098    0.02884   1.768  0.07742 .  
# pca_composite  0.10201    0.02885   3.535  0.00043 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.8217 on 810 degrees of freedom
# Multiple R-squared:  0.0152,	Adjusted R-squared:  0.01398 
# F-statistic:  12.5 on 1 and 810 DF,  p-value: 0.0004303

#Look for most important variables
sort(abs(PCA$loadings[,1])) #ndvi, no health ins, perc not white, mean-no2, foreign born, tree canopy, pop dens, median income, mean_pm25, mean pm10, automobile, unemploted, perc over 65
sort(abs(PCA$loadings[,2])) #commute, pop density, automobile, mean_pm10, AC,gini_index, mean_pm25, park_access
#First PC only account for 28 percent of the variation in the data, second hsa 15 percent do we have to think about that? 

