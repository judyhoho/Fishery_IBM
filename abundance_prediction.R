########################################################
# Squid abundance prediction generate
########################################################
## (1) Process raw catch data into CPUE 
## (2) Find optimal GLM for explaining catchability, final GLM: lmer(log(CPUE + const) ~ HP + Light + Tonnage + (1 | Year), data = catch.process)
## (3) Find optimal GAM for predicting abundance on condition of fixed catchability with k-fold cross validation, final GAM: gam(pre_CPUE ~ s(SSTboat, k=3, by = Month) + s(Bathy, k=3, by = Month) + s(Lon, k=3) + s(Lat, k=3) + Month, family = Gamma(link=log), na.rm = TRUE, data = catch.process)

####################################################
# Set working directory and load required packages 
setwd("/Users/judyhoho/Desktop/Gou Lab/Yun")

library(tidyverse)
library(mgcv)
library(lme4)
library(ggplot2)

####################################################
# GLM + GAM (training workflow)

# CPUE calculation
## data input
input_variables <- c("No", "Lat", "Lon","SSTboat", "Bathy", "Year", "Month", "Day", "Lunar.Year", "Lunar.Month", "Lunar.Day", "鎖管", "Tonnage", "HP", "Light", "CTNO", "作業漁法_代號")
catch <- as_tibble(rbind(read.csv("臺灣北部海域火誘網漁獲資料_2019.csv", 
                                  fileEncoding = 'big5')[,input_variables],
                         read.csv("臺灣北部海域火誘網漁獲資料_2020.csv", 
                                  fileEncoding = 'big5')[,input_variables]))
## clear content with NA value and add a column 'squid' as the same as 鎖管 
catch$SSTboat <- as.numeric(catch$SSTboat)
catch <- catch %>%
  filter(作業漁法_代號 == 1 & Light > 0 & Tonnage > 0 & HP > 0 & Bathy > 0 & Bathy <= 300 & SSTboat > 0 &
           鎖管 >= 0 & Lat >= 24 & Lat <= 30 & Lon >= 120 & Lon <= 126)
catch$squid <- catch$鎖管

## process the effort based on each grid (transfer into grids of 0.1°*0.1°)
catch$Lat <- round(catch$Lat, digits = 1)
catch$Lon <- round(catch$Lon, digits = 1)
catch.process <- catch %>%
  group_by(Year, Month, Lon, Lat) %>%
  summarise(CPUE = mean(squid), Light = mean(Light), 
            Tonnage = mean(Tonnage), HP = mean(HP), 
            SSTboat = mean(SSTboat), Bathy = mean(Bathy))

####################################################
# CPUE standardization
## GLMM with the whole data set
## k-fold cross validation for GAM (k=5)

# GLMM
const <- 0.5* min(catch.process$CPUE[catch.process$CPUE > 0])
mod.glmm <- lmer(log(CPUE + const) ~ HP + Light + Tonnage + (1 | Year), data = catch.process)
mod_null <- lmer(log(CPUE + const) ~ 1 + (1 | Year), data = catch.process)
mod1_1 <- lmer(log(CPUE + const) ~ HP + (1 | Year), data = catch.process)
mod1_2 <- lmer(log(CPUE + const) ~ Light + (1 | Year), data = catch.process)
mod1_3 <- lmer(log(CPUE + const) ~ Tonnage + (1 | Year), data = catch.process)
mod2_1 <- lmer(log(CPUE + const) ~ HP + Light + (1 | Year), data = catch.process)
mod2_2 <- lmer(log(CPUE + const) ~ HP + Tonnage + (1 | Year), data = catch.process)
mod2_3 <- lmer(log(CPUE + const) ~ Tonnage + Light + (1 | Year), data = catch.process)
anova(mod_null, mod1_1, mod1_2, mod1_3, mod2_1, mod2_2, mod2_3, mod.glmm) #mod.glmm is the best

res <- residuals(mod.glmm, type = 'response')
fixed.intercept <- fixef(mod.glmm)[['(Intercept)']]
total.intercept <- data.frame(Year = rownames(ranef(mod.glmm)$Year),
                              Intercept = ranef(mod.glmm)$Year + fixed.intercept)
colnames(total.intercept) <- c('Year', 'Intercept')
total.intercept$Year <- as.numeric(total.intercept$Year)
mod.intercept <- lm(Intercept ~ Year, data = total.intercept)
inter.month <- predict(mod.intercept, newdata = data.frame(Year = 2021))
fixed.q <- predict(mod.glmm, newdata = data.frame(HP = mean(catch.process$HP),
                                                  Light = mean(catch.process$Light),
                                                  Tonnage = mean(catch.process$Tonnage),
                                                  Year = as.factor(2020)))
fixed.q <- fixed.q - total.intercept$Intercept[nrow(total.intercept)] + inter.month
catch.process$pre_CPUE <- exp(res + fixed.q) - const

# GAM (k-fold cross validation)
## define RMSE, MAE function
rmse <- function(actual, pred) sqrt(mean((actual - pred)^2))
mae <- function(actual, pred) mean(abs(actual - pred))

set.seed(123)
k <- 5  #number of folds
catch.process$fold <- sample(rep(1:k, length.out = nrow(catch.process))) #create fold IDs

k.df <- data.frame()
aic.df <- data.frame()
catch.test.combined <- data.frame()
for (i in 1:k){
  catch.train <- catch.process %>%
    filter(fold != i)
  catch.test <- catch.process %>%
    filter(fold == i)
  
  #GAM
  catch.train$Month <- as.factor(catch.train$Month)
  catch.train <- catch.train %>%
    filter(SSTboat > 0, Bathy > 0, pre_CPUE >= 0)
  catch.test$Month <- as.factor(catch.test$Month)
  catch.test <- catch.test %>%
    filter(SSTboat > 0, Bathy > 0, pre_CPUE >= 0)
  
  gam1 <- gam(pre_CPUE ~ s(SSTboat, k=3) + s(Bathy, k=3) + Month, 
              family = Gamma(link=log), na.rm = TRUE, data = catch.train)
  gam4 <- gam(pre_CPUE ~ s(SSTboat, k=3, by = Month) + s(Bathy, k=3) + Month, 
              family = Gamma(link=log), na.rm = TRUE, data = catch.train)
  gam5 <- gam(pre_CPUE ~ s(SSTboat, k=3) + s(Bathy, k=3, by = Month) + Month, 
              family = Gamma(link=log), na.rm = TRUE, data = catch.train)
  gam6 <- gam(pre_CPUE ~ s(SSTboat, k=3, by = Month) + s(Bathy, k=3, by = Month) + Month, 
              family = Gamma(link=log), na.rm = TRUE, data = catch.train)
  gam1_2 <- gam(pre_CPUE ~ s(SSTboat, k=3) + s(Bathy, k=3) + s(Lon, k=3) + s(Lat, k=3) + Month, 
                family = Gamma(link=log), na.rm = TRUE, data = catch.train)
  gam4_2 <- gam(pre_CPUE ~ s(SSTboat, k=3, by = Month) + s(Bathy, k=3) + s(Lon, k=3) + s(Lat, k=3) + Month, 
                family = Gamma(link=log), na.rm = TRUE, data = catch.train)
  gam5_2 <- gam(pre_CPUE ~ s(SSTboat, k=3) + s(Bathy, k=3, by = Month) + s(Lon, k=3) + s(Lat, k=3) + Month, 
                family = Gamma(link=log), na.rm = TRUE, data = catch.train)
  gam6_2 <- gam(pre_CPUE ~ s(SSTboat, k=3, by = Month) + s(Bathy, k=3, by = Month) + s(Lon, k=3) + s(Lat, k=3) + Month, 
                family = Gamma(link=log), na.rm = TRUE, data = catch.train)
  
  aic <- AIC(gam1, gam1_2, gam4, gam4_2, gam5, gam5_2, gam6, gam6_2)
  aic.df <- rbind(aic.df, data.frame(k = i,
                                     df = aic$df,
                                     AIC = aic$AIC,
                                     model = c('1', '1_2', '4', '4_2', '5', '5_2', '6', '6_2')))
  
  cwb.combined <- data.frame()
  for (month in 4:9){
    day <- 15
    m <- month
    month <- as.character(month)
    if (nchar(month) == 1) month <- paste0("0", month)
    day <- as.character(day)
    if (nchar(day) == 1) day <- paste0("0", day)
    
    if (file.exists(paste0("氣象局水溫預測資料/2021/", "OCM_FRI_2021", month, day, "10.txt"))){
      cwb <- read.table(paste0("氣象局水溫預測資料/2021/", "OCM_FRI_2021", month, day, "10.txt"))
      colnames(cwb) <- c("Lon", "Lat", "SSTboat", "Bathy", "BT", "MLD")
      cwb$Lon <- round(cwb$Lon, 1)
      cwb$Lat <- round(cwb$Lat, 1)
      cwb <- cwb %>%
        group_by(Lon, Lat) %>%
        summarise(SSTboat = mean(SSTboat), Bathy = mean(Bathy)) %>%
        filter(SSTboat >= 0, Bathy >= 0)
      cwb$Month <- as.factor(m)
      cwb$Pred_gam1 <- predict(gam1, cwb, type = 'response')
      cwb$Pred_gam1_2 <- predict(gam1_2, cwb, type = 'response')
      cwb$Pred_gam4 <- predict(gam4, cwb, type = 'response')
      cwb$Pred_gam4_2 <- predict(gam4_2, cwb, type = 'response')
      cwb$Pred_gam5 <- predict(gam5, cwb, type = 'response')
      cwb$Pred_gam5_2 <- predict(gam5_2, cwb, type = 'response')
      cwb$Pred_gam6 <- predict(gam6, cwb, type = 'response')
      cwb$Pred_gam6_2 <- predict(gam6_2, cwb, type = 'response')
      cwb.combined <- rbind(cwb.combined, cwb)
      
    }
  }
  catch.test <- catch.test %>%
    left_join(cwb.combined, by = c('Month' = 'Month', 'Lon' = 'Lon', 'Lat' = 'Lat'))
  catch.test <- catch.test[!is.na(catch.test$Pred_gam5),]
  catch.test <- catch.test[catch.test$Bathy.y <= 300,] # exclude those with Bathy > 300
  catch.test.combined <- rbind(catch.test.combined, catch.test)
  rmse_gam1 <- rmse(catch.test$pre_CPUE, catch.test$Pred_gam1)
  rmse_gam1_2 <- rmse(catch.test$pre_CPUE, catch.test$Pred_gam1_2)
  rmse_gam4 <- rmse(catch.test$pre_CPUE, catch.test$Pred_gam4)
  rmse_gam4_2 <- rmse(catch.test$pre_CPUE, catch.test$Pred_gam4_2)
  rmse_gam5 <- rmse(catch.test$pre_CPUE, catch.test$Pred_gam5)
  rmse_gam5_2 <- rmse(catch.test$pre_CPUE, catch.test$Pred_gam5_2)
  rmse_gam6 <- rmse(catch.test$pre_CPUE, catch.test$Pred_gam6)
  rmse_gam6_2 <- rmse(catch.test$pre_CPUE, catch.test$Pred_gam6_2)
  mae_gam1 <- mae(catch.test$pre_CPUE, catch.test$Pred_gam1)
  mae_gam1_2 <- mae(catch.test$pre_CPUE, catch.test$Pred_gam1_2)
  mae_gam4 <- mae(catch.test$pre_CPUE, catch.test$Pred_gam4)
  mae_gam4_2 <- mae(catch.test$pre_CPUE, catch.test$Pred_gam4_2)
  mae_gam5 <- mae(catch.test$pre_CPUE, catch.test$Pred_gam5)
  mae_gam5_2 <- mae(catch.test$pre_CPUE, catch.test$Pred_gam5_2)
  mae_gam6 <- mae(catch.test$pre_CPUE, catch.test$Pred_gam6)
  mae_gam6_2 <- mae(catch.test$pre_CPUE, catch.test$Pred_gam6_2)
  k.df <- rbind(k.df, data.frame(K = i,
                                 rmse_gam1 = rmse_gam1,
                                 rmse_gam1_2 = rmse_gam1_2,
                                 rmse_gam4 = rmse_gam4,
                                 rmse_gam4_2 = rmse_gam4_2,
                                 rmse_gam5 = rmse_gam5,
                                 rmse_gam5_2 = rmse_gam5_2,
                                 rmse_gam6 = rmse_gam6,
                                 rmse_gam6_2 = rmse_gam6_2,
                                 mae_gam1 = mae_gam1, 
                                 mae_gam1_2 = mae_gam1_2,
                                 mae_gam4 = mae_gam4, 
                                 mae_gam4_2 = mae_gam4_2,
                                 mae_gam5 = mae_gam5, 
                                 mae_gam5_2 = mae_gam5_2, 
                                 mae_gam6 = mae_gam6,
                                 mae_gam6_2 = mae_gam6_2))
}

compare.df <- aic.df %>%
  group_by(model) %>%
  summarise(AIC = mean(AIC))

compare.df$RMSE <- (colSums(k.df)/5)[2:9]
compare.df$MAE <- (colSums(k.df)/5)[10:17]


##################################################
# GLM + GAM (generate prediction workflow)
## make daily predictions

## data input
input_variables <- c("No", "Lat", "Lon","SSTboat", "Bathy", "Year", "Month", "Day", "Lunar.Year", "Lunar.Month", "Lunar.Day", "鎖管", "Tonnage", "HP", "Light", "CTNO", "作業漁法_代號")
catch <- as_tibble(rbind(read.csv("臺灣北部海域火誘網漁獲資料_2019.csv", 
                                  fileEncoding = 'big5')[,input_variables],
                         read.csv("臺灣北部海域火誘網漁獲資料_2020.csv", 
                                  fileEncoding = 'big5')[,input_variables]))

## clear content with NA value and add a column 'squid' as the same as 鎖管 
catch$SSTboat <- as.numeric(catch$SSTboat)
catch <- catch %>%
  filter(作業漁法_代號 == 1 & Light > 0 & Tonnage > 0 & HP > 0 & Bathy > 0 & Bathy <= 300 & SSTboat > 0 &
           鎖管 >= 0 & Lat >= 24 & Lat <= 30 & Lon >= 120 & Lon <= 126)
catch$squid <- catch$鎖管

## process the effort based on each grid (transfer into grids of 0.1°*0.1°)
catch$Lat <- round(catch$Lat, digits = 1)
catch$Lon <- round(catch$Lon, digits = 1)
catch.process <- catch %>%
  group_by(Year, Month, Lon, Lat) %>%
  summarise(CPUE = mean(squid), Light = mean(Light), 
            Tonnage = mean(Tonnage), HP = mean(HP), 
            SSTboat = mean(SSTboat), Bathy = mean(Bathy))

## GLMM with the whole data set
const <- 0.5 * min(catch.process$CPUE[catch.process$CPUE > 0])
mod.glmm <- lmer(log(CPUE + const) ~ HP + Light + Tonnage + (1 | Year), data = catch.process)
mod_null <- lmer(log(CPUE + const) ~ 1 + (1 | Year), data = catch.process)
mod1_1 <- lmer(log(CPUE + const) ~ HP + (1 | Year), data = catch.process)
mod1_2 <- lmer(log(CPUE + const) ~ Light + (1 | Year), data = catch.process)
mod1_3 <- lmer(log(CPUE + const) ~ Tonnage + (1 | Year), data = catch.process)
mod2_1 <- lmer(log(CPUE + const) ~ HP + Light + (1 | Year), data = catch.process)
mod2_2 <- lmer(log(CPUE + const) ~ HP + Tonnage + (1 | Year), data = catch.process)
mod2_3 <- lmer(log(CPUE + const) ~ Tonnage + Light + (1 | Year), data = catch.process)
anova(mod_null, mod1_1, mod1_2, mod1_3, mod2_1, mod2_2, mod2_3, mod.glmm) #mod.glmm is the best

res <- residuals(mod.glmm, type = 'response')
fixed.intercept <- fixef(mod.glmm)[['(Intercept)']]
total.intercept <- data.frame(Year = rownames(ranef(mod.glmm)$Year),
                              Intercept = ranef(mod.glmm)$Year + fixed.intercept)
colnames(total.intercept) <- c('Year', 'Intercept')
total.intercept$Year <- as.numeric(total.intercept$Year)
mod.intercept <- lm(Intercept ~ Year, data = total.intercept)
inter.month <- predict(mod.intercept, newdata = data.frame(Year = 2021))
fixed.q <- predict(mod.glmm, newdata = data.frame(HP = mean(catch.process$HP),
                                                  Light = mean(catch.process$Light),
                                                  Tonnage = mean(catch.process$Tonnage),
                                                  Year = as.factor(2020)))
fixed.q <- fixed.q - total.intercept$Intercept[nrow(total.intercept)] + inter.month
catch.process$pre_CPUE <- exp(res + fixed.q) - const

## GAM with the whole data set
catch.process$Month <- as.factor(catch.process$Month)
catch.process <- catch.process %>%
  filter(SSTboat > 0, Bathy > 0, pre_CPUE >= 0)
gam <- gam(pre_CPUE ~ s(SSTboat, k=3, by = Month) + s(Bathy, k=3, by = Month) + s(Lon, k=3) + s(Lat, k=3) + Month, 
           family = Gamma(link=log), na.rm = TRUE, data = catch.process)

## generate daily prediction
for (m in 4:9){
  for (day in 1:30){
    month <- as.character(m)
    if (nchar(month) == 1) month <- paste0("0", month)
    day <- as.character(day)
    if (nchar(day) == 1) day <- paste0("0", day)
    
    if (file.exists(paste0("氣象局水溫預測資料/2021/", "OCM_FRI_2021", month, day, "10.txt"))){
      cwb <- read.table(paste0("氣象局水溫預測資料/2021/", "OCM_FRI_2021", month, day, "10.txt"))
      colnames(cwb) <- c("Lon", "Lat", "SSTboat", "Bathy", "BT", "MLD")
      cwb$Month <- as.factor(m)
      Pred_gam <- predict(gam, cwb, type = 'response')
      cwb <- cbind(Pred_gam, cwb)
      write.table(cwb, paste0("map data/predict/predict_2021_revised_final/Pred_2021", month, day, ".txt"))
    }
  }
}



