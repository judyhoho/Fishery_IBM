########################################################
# IBM evaluation : 
## (1) All features used for model evaluation are calculated as below, including turning angle, step length, maximum distance from harbor, Frechet distance, time allocation and fishing location
## (2) Two filter criteria were used to select parameter sets that are further used to compare performance metrics among scenarios
## (3) PCA for parameter values
########################################################


########################################################
# Load packages needed 
########################################################
library(tidyverse)
library(ggplot2)


########################################################
# TA distribution 
## read in raw simulation results and calculate the turning angle between each step
########################################################
Model.name <- 'Optimal'
for (i in 1:27){
  Raw <- readRDS(paste0('/Users/judyhoho/Desktop/IBM/IBM_hexagon/simulation results/', Model.name, '_chunk_', i, '.rds'))
  All.Ship.TA <- data.frame()
  for (r in 1:1000){
    for (s in unique(Raw[[r]]$fishers$Ship)){
      Ship <- Raw[[r]]$fishers %>%
        filter(Ship == s)
      Ship.TA <- Ship[1,]
      for (j in 2:nrow(Ship)){
        if(Ship$Pos.X[j] != Ship$Pos.X[j-1] | Ship$Pos.Y[j] != Ship$Pos.Y[j-1]){
          Ship.TA <- rbind(Ship.TA, Ship[j,])
        }
      }
      Ship.TA$Vector.X <- NA
      Ship.TA$Vector.Y <- NA
      Ship.TA$Angle <- NA
      for (j in 2:nrow(Ship.TA)){
        Ship.TA$Vector.X[j] <- round(Ship.TA$Pos.X[j]-Ship.TA$Pos.X[j-1], 6)
        Ship.TA$Vector.Y[j] <- round(Ship.TA$Pos.Y[j]-Ship.TA$Pos.Y[j-1], 6)
      }
      
      Ship.TA$a1 <- atan2(Ship.TA$Vector.Y, Ship.TA$Vector.X) *180/pi
      Ship.TA$a2 <- c(NA, Ship.TA$a1[1:nrow(Ship.TA)-1])
      Ship.TA$Angle <- ifelse(abs(Ship.TA$a1 - Ship.TA$a2) <= 180, Ship.TA$a1 - Ship.TA$a2, 
                              (360 - abs(Ship.TA$a1 - Ship.TA$a2)) * (-1) * sign(Ship.TA$a1 - Ship.TA$a2))
      Ship.TA$Month <- Raw[[r]]$Params$Month.m
      All.Ship.TA <- rbind(All.Ship.TA, Ship.TA)
      rm(Ship, Ship.TA)
      cat(sprintf("Ship %d of %s is done\n", s, r))
    }
  }
  saveRDS(All.Ship.TA, file = paste0('/Users/judyhoho/Desktop/IBM/IBM_hexagon/simulation results/Organized results/TA/', Model.name, "_", i, ".rds"))
  rm(Raw, All.Ship.TA)
}


########################################################
# TA overlap
## read in VDR to calculate the overlap between real vessels and simulated vessels
## results are stored as .RData called Overlap_Model
########################################################

load('/Users/judyhoho/Desktop/IBM/data/Angle_2021.RData') #VDR

## load in simulated data (processed in previous section) and combine them
Model.name <- "Null"
All.TA <- data.frame() #turning angle data
for (i in c(1:27)){
  Raw.TA <- read_rds(paste0('/Users/judyhoho/Desktop/IBM/IBM_hexagon/simulation results/Organized results/TA/', Model.name, "_", i, ".rds"))
  All.TA <- rbind(All.TA, Raw.TA)
  print(i)
  rm(Raw.TA)
}

Param.grid <- expand.grid(HR = c(0.5, 1, 1.5),
                          Lag = c(4, 12, 24),
                          Pun = c(0.2, 0.4, 0.6, 0.8, 1),
                          Pay = c(0.2, 0.4, 0.6, 0.8, 1),
                          S = 1:20)
Combine.Ship.TA.Overlap <- data.frame()
for (m in 4:9){
  VDR.TA <- All.VDR.TA %>%
    filter(Month == m)
  VDR.D <- density(na.omit(VDR.TA$Angle), from = min(na.omit(VDR.TA$Angle)), to = max(na.omit(VDR.TA$Angle)))
  for (n in 1:nrow(Param.grid)){
    IBM.TA <- All.TA %>%
      filter(Month == m, Pay == Param.grid$Pay[n], 
             Lag == Param.grid$Lag[n], Pun == Param.grid$Pun[n], 
             Fish.Hr == Param.grid$HR[n], S == Param.grid$S[n])
    IBM.D <- density(na.omit(IBM.TA$Angle), from = min(na.omit(IBM.TA$Angle)), to = max(na.omit(IBM.TA$Angle)))
    joint <- pmin(VDR.D$y, IBM.D$y)
    dx <- VDR.D$x[2] - VDR.D$x[1]
    Overlap.area <- sum(joint) * dx
    Combine.Ship.TA.Overlap <- rbind(Combine.Ship.TA.Overlap, data.frame(Overlap = Overlap.area,
                                                                         Pay = Param.grid$Pay[n],
                                                                         Lag = Param.grid$Lag[n],
                                                                         Pun = Param.grid$Pun[n],
                                                                         Fish.Hr = Param.grid$HR[n],
                                                                         S = Param.grid$S[n],
                                                                         Month = m,
                                                                         Model = Model.name))
    print(paste0('TA overlap', n))
    rm(IBM.TA, IBM.D, joint, dx, Overlap.area)
  }
  rm(VDR.TA)
}
save(Combine.Ship.TA.Overlap, file = "/Users/judyhoho/Desktop/IBM/IBM_hexagon/simulation results/Organized results/TA/Overlap_Null.RData") #save results for each scenario

All.TA.Overlap <- data.frame() #combine results of all scenarios
load("/Users/judyhoho/Desktop/IBM/IBM_hexagon/simulation results/Organized results/TA/Overlap_Null.RData")
All.TA.Overlap <- bind_rows(All.TA.Overlap, Combine.Ship.TA.Overlap)

ggplot(All.TA.Overlap) +
  geom_boxplot(aes(x = Model, y = Overlap))


########################################################
# Max distance from harbor
## do z-score transformation according to real vessel data
########################################################

## record maximum distance for each simulated vessel in each iteration
Harbour.Position <- c(121.6333, 25.33270)
Combine.Ship.Max <- data.frame()
Model.name <- 'Optimal'
for (i in 1:27){
  Raw <- readRDS(paste0('/Users/judyhoho/Desktop/IBM/IBM_hexagon/simulation results/', Model.name, '_chunk_retry_', i, '.rds'))
  Max.List <- c()
  for (n in 1:length(Raw)){
    Max.Ship <- Raw[[n]]$fishers
    Max.Ship$Max.Harbor <- sqrt((abs(Max.Ship$Pos.X-Harbour.Position[1])*100)^2 + 
                                  (abs(Max.Ship$Pos.Y-Harbour.Position[2])*111)^2)
    Max.List <- Max.Ship %>%
      group_by(Ship) %>%
      summarise(Max = max(Max.Harbor))
    Combine.Ship.Max <- rbind(Combine.Ship.Max, data.frame(Pay = Max.Ship$Pay[1],
                                                           Lag = Max.Ship$Lag[1],
                                                           Pun = Max.Ship$Pun[1],
                                                           Fish.Hr = Max.Ship$Fish.Hr[1],
                                                           Max = Max.List$Max,
                                                           Ship = Max.List$Ship,
                                                           S = Max.Ship$S[1],
                                                           Month = Raw[[n]]$Params$Month.m,
                                                           Model = Model.name))
    print(paste0('Chunk: ', i, ' Max', n, ' done.'))
    rm(Max.Ship, Max.List)
  }
  rm(Raw)
}
save(Combine.Ship.Max, file = paste0("/Users/judyhoho/Desktop/IBM/IBM_hexagon/simulation results/Organized results/MAX/Max_", Model.name, ".RData"))

## z-score transformation using real VDR
VDR.MAX <- read.csv('/Users/judyhoho/Desktop/Ho Lab 638/R code/IBM/Validation/VDR_Max.csv')[,2:4]
VDR.MAX <- VDR.MAX %>%
  filter(Month > 3 & Month < 10)
VDR.MAX <- VDR.MAX %>%
  group_by(Month) %>%
  summarise(Mean = mean(Max), SD = sd(Max))

load('/Users/judyhoho/Desktop/IBM/IBM_hexagon/simulation results/Organized results/MAX/Max_Optimal.RData')
Optimal <- data.frame()
for (m in 4:9){
  Max.month <- Combine.Ship.Max %>%
    filter(Month == m) %>%
    mutate(Z.score = (Max-VDR.MAX$Mean[m-3])/VDR.MAX$SD[m-3])
  Optimal <- rbind(Optimal, Max.month)
  rm(Max.month)
  print(m)
}

Max.Z.score <- bind_rows(Optimal, Timeindiff, NullResource, Null) #combine results of all scenarios
save(Max.Z.score, file = '/Users/judyhoho/Desktop/IBM/IBM_hexagon/simulation results/Organized results/MAX/Max_All.RData') #save combined results

ggplot(Max.Z.score) +
  geom_boxplot(aes(x = Model, y = Z.score))


########################################################
# Travel distances overlap
## calculate step length for each step and compare that to real VDR
########################################################

## load in simulated angle data
Model.name <- "Optimal"
All.TA <- data.frame() #turning angle data
for (i in c(1:27)){
  Raw.TA <- read_rds(paste0('/Users/judyhoho/Desktop/IBM/IBM_hexagon/simulation results/Organized results/TA/', Model.name, "_retry_", i, ".rds"))
  All.TA <- rbind(All.TA, Raw.TA) 
  print(i)
}

Param.grid <- expand.grid(HR = c(0.5, 1, 1.5),
                          Lag = c(4, 12, 24),
                          Pun = c(0.2, 0.4, 0.6, 0.8, 1),
                          Pay = c(0.2, 0.4, 0.6, 0.8, 1),
                          S = 1:20)

## load in real VDR
load('/Users/judyhoho/Desktop/IBM/data/Angle_2021.RData')
All.VDR.TA$Distance <- NA
All.VDR.TA$Distance[2:nrow(All.VDR.TA)] <- sqrt((All.VDR.TA$Longitude[2:nrow(All.VDR.TA)]-All.VDR.TA$Longitude[1:nrow(All.VDR.TA)-1])^2 + (All.VDR.TA$Latitude[2:nrow(All.VDR.TA)]-All.VDR.TA$Latitude[1:nrow(All.VDR.TA)-1])^2)

## calculate overlap area between simulated vessels and real vessels
Combine.Ship.Dist.Overlap <- data.frame()
for (m in 4:9){
  VDR.m <- All.VDR.TA %>%
    filter(Month == m)
  VDR.m <- VDR.m[!is.na(VDR.m$Vector),]
  VDR.Dist <- density(na.omit(VDR.m$Distance), from = min(na.omit(VDR.m$Distance)), to = max(na.omit(VDR.m$Distance)))
  for (n in 1:nrow(Param.grid)){
    IBM.TA <- All.TA %>%
      filter(Month == m, Pay == Param.grid$Pay[n], 
             Lag == Param.grid$Lag[n], Pun == Param.grid$Pun[n], 
             Fish.Hr == Param.grid$HR[n], S == Param.grid$S[n])
    IBM.TA$Distance <- NA
    IBM.TA$Distance[2:nrow(IBM.TA)] <- sqrt((IBM.TA$Pos.X[2:nrow(IBM.TA)]-IBM.TA$Pos.X[1:nrow(IBM.TA)-1])^2 + (IBM.TA$Pos.Y[2:nrow(IBM.TA)]-IBM.TA$Pos.Y[1:nrow(IBM.TA)-1])^2)
    IBM.TA <- IBM.TA[!is.na(IBM.TA$Vector.X),]
    IBM.Dist <- density(na.omit(IBM.TA$Distance), from = min(na.omit(IBM.TA$Distance)), to = max(na.omit(IBM.TA$Distance)))
    joint <- pmin(VDR.Dist$y, IBM.Dist$y)
    dx <- VDR.Dist$x[2] - VDR.Dist$x[1]
    Overlap.area <- sum(joint) * dx
    Combine.Ship.Dist.Overlap <- rbind(Combine.Ship.Dist.Overlap, 
                                       data.frame(Dist.Overlap = Overlap.area,
                                                  Pay = IBM.TA$Pay[1],
                                                  Lag = IBM.TA$Lag[1],
                                                  Pun = IBM.TA$Pun[1],
                                                  Fish.Hr = IBM.TA$Fish.Hr[1],
                                                  S = IBM.TA$S[1],
                                                  Month = m))
    print(paste0('Month: ', m, ' Param: ', n, ' done.'))
    rm(IBM.TA)
  }
}
Combine.Ship.Dist.Overlap$Model <- Model.name
save(Combine.Ship.Dist.Overlap, file = paste0('/Users/judyhoho/Desktop/IBM/IBM_hexagon/simulation results/Organized results/STEP/Step_', Model.name, '.RData')) #save result for each scenario
All.Step <- bind_rows(Optimal, Timeindiff, NullResource, Null) #combine results of all scenarios
save(Max.Z.score, file = '/Users/judyhoho/Desktop/IBM/IBM_hexagon/simulation results/Organized results/MAX/Step_Combine.RData') #save combined results

All.Step %>%
  ggplot() +
  geom_boxplot(aes(x = Model, y = Dist.Overlap))


########################################################
# Frechet distance 
## calculate Frechet distance (recommend to run on server, the scripts below is conducted on server)
########################################################
library(tidyverse)
library(tableHTML)
library(parallel)
setwd("/home/pojuke/YunHo/IBM_hexagon/Frechet")
source("Frechet_source.R")
load('Angle_2021.RData') #VDR

groups <- readRDS('Optimal_groups.rds') #one simulated vessel in one iteration per group
chunk_size <- 5000
group_indices <- split(seq_along(groups), ceiling(seq_along(groups)/chunk_size))

process_one_group <- function(IBM.track, All.VDR.TA) {
  VDR.TA <- All.VDR.TA %>%
    filter(Month == IBM.track$Month[1])
  bind_rows(lapply(unique(VDR.TA$Ship), function(k) {
    
    VDR.track <- VDR.TA %>%
      filter(Ship == k)
    f <- Frechet(as.matrix(IBM.track[, c(4, 5)]), as.matrix(VDR.track[, c(6, 5)]))
    
    tibble(Month = IBM.track$Month[1],
           Pay = IBM.track$Pay[1],
           Lag = IBM.track$Lag[1],
           Pun = IBM.track$Pun[1],
           Fish.Hr = IBM.track$Fish.Hr[1],
           S = IBM.track$S[1],
           Ship = IBM.track$Ship[1],
           VDR = k,
           Frechet.d = f)
  }))
}

for (i in seq_along(group_indices)) {
  idx <- group_indices[[i]]
  cat(sprintf("Processing chunk %d / %d (%d groups)\n", i, length(group_indices), length(idx)))
  chunk_result <- bind_rows(mclapply(groups[idx], process_one_group,
                                     All.VDR.TA = All.VDR.TA, mc.cores = 6))
  write_csv(chunk_result,
            file = sprintf("Frechet_Optimal_chunk_%03d.csv", i))
  rm(chunk_result)
  gc()
}

## organize the Frechet distance calculated at server (combine them into one)
Combine.Frechet <- data.frame()
for (i in 1:27){
  fname <- sprintf("Frechet_Optimal_chunk_%03d.csv", i)
  Fre <- read.csv(paste0('/Users/judyhoho/Desktop/IBM/IBM_hexagon/simulation results/Organized results/Frechet/', fname))
  print(fname)
  Combine.Frechet <- bind_rows(Combine.Frechet, Fre)
  rm(Fre)
}
Combine.Frechet$Model <- "Optimal"
save(Combine.Frechet, file = '/Users/judyhoho/Desktop/IBM/IBM_hexagon/simulation results/Organized results/Frechet/Frechet_Optimal.RData') #save result for each scenario


########################################################
# Fishing behavior (time allocation)
## calculate the time proportion of fishing to traveling during operation time
########################################################

## calculate that for VDR and save the output
load('/Users/judyhoho/Desktop/Gou Lab/Yun/Lux pattern/VDR_th_2021.RData')
groups <- split(VDR.all, list(VDR.all$Ship, VDR.all$Trip), drop = TRUE)
VDR.Behavior <- data.frame()
for (i in 1:length(groups)){
  VDR.Ship <- groups[[i]]
  df <- VDR.Ship[1,]
  for (j in 2:nrow(VDR.Ship)){
    if (VDR.Ship$Work.Code[j] != VDR.Ship$Work.Code[j-1]){
      df <- bind_rows(df, VDR.Ship[j,])
    }
  }
  VDR.Behavior <- bind_rows(VDR.Behavior, df)
  rm(df)
  print(i)
}
save(VDR.Behavior, file = '/Users/judyhoho/Desktop/IBM/IBM_hexagon/simulation results/Organized results/FISH/VDR_Behavior.RData')
VDR.Behavior <- VDR.Behavior %>%
  filter(Work.Code != 0) %>%
  group_by(Ship, Month, Work.Code) %>%
  summarise(Total = sum(na.omit(Timediff))) %>%
  mutate(Work = ifelse(Work.Code == 2, "Fish", "Travel")) %>%
  group_by(Month, Ship, Work) %>%
  summarise(Total = sum(Total)) %>%
  pivot_wider(names_from = Work, values_from = Total) %>%
  mutate(Fish_Travel_ratio = Fish/Travel)

## calculate that for IBM
IBM.Behavior <- data.frame()
for (k in 1:27){
  IBM <- readRDS(paste0('/Users/judyhoho/Desktop/IBM/IBM_hexagon/simulation results/Optimal_chunk_', k, '.rds'))
  for (i in 1:length(IBM)){
    for (s in 1:5){
      IBM.Ship <- IBM[[i]]$fishers %>% filter(Ship == s)
      df <- IBM.Ship[1,]
      for (j in 2:nrow(IBM.Ship)){
        if (IBM.Ship$Next.Fish[j] != IBM.Ship$Next.Fish[j-1]){
          df <- bind_rows(df, IBM.Ship[j,])
        }
      }
      df <- bind_rows(df, IBM.Ship[nrow(IBM.Ship),])
      df.final <- data.frame()
      for (j in 2:(nrow(df)-1)){
        if ((df$Pos.X[j] != df$Pos.X[j-1] & df$Pos.Y[j] != df$Pos.Y[j-1]) | 
            (df$Pos.X[j] != df$Pos.X[j+1] & df$Pos.Y[j] != df$Pos.Y[j+1])){
          df.final <- bind_rows(df.final, df[j,])
        }
      }
      Fish <- 0
      Travel <- 0
      for (j in 2:nrow(df.final)){
        if (df.final$Cumu.Gain[j] != df.final$Cumu.Gain[j-1]){
          Fish <- Fish + df.final$Cumu.Hour[j] - df.final$Cumu.Hour[j-1]
        } else{Travel <- Travel + df.final$Cumu.Hour[j] - df.final$Cumu.Hour[j-1]}
      }
      IBM.Behavior <- bind_rows(IBM.Behavior, data.frame(Pay = IBM.Ship$Pay[1],
                                                         Lag = IBM.Ship$Lag[1],
                                                         Pun = IBM.Ship$Pun[1],
                                                         Fish.Hr = IBM.Ship$Fish.Hr[1],
                                                         Ship = s,
                                                         S = IBM.Ship$S[1],
                                                         Month = IBM[[i]]$Params$Month.m,
                                                         Travel.time = Travel,
                                                         Fish.time = Fish))
      print(paste0("K = ", k, ", Iteration ", i, "_Ship: ", s))
      rm(IBM.Ship, df, df.final)
    }
  }
  rm(IBM)
}

IBM.Behavior$Model <- "Optimal"
saveRDS(IBM.Behavior, '/Users/judyhoho/Desktop/IBM/IBM_hexagon/simulation results/Organized results/FISH/Behavior_Optimal.rds') #save result for each scenario

## calculate the Z.score for fish time allocation
Optimal <- readRDS( '/Users/judyhoho/Desktop/IBM/IBM_hexagon/simulation results/Organized results/FISH/Behavior_Optimal.rds')
Timeindiff <- readRDS( '/Users/judyhoho/Desktop/IBM/IBM_hexagon/simulation results/Organized results/FISH/Behavior_Timeindiff.rds')
NullResource <- readRDS( '/Users/judyhoho/Desktop/IBM/IBM_hexagon/simulation results/Organized results/FISH/Behavior_NullResource.rds')
Null <- readRDS( '/Users/judyhoho/Desktop/IBM/IBM_hexagon/simulation results/Organized results/FISH/Behavior_Null.rds')
load('/Users/judyhoho/Desktop/IBM/IBM_hexagon/simulation results/Organized results/FISH/VDR_Behavior.RData')
VDR <- VDR.Behavior %>%
  filter(Work.Code != 0) %>%
  group_by(Ship, Month, Work.Code) %>%
  summarise(Total = sum(na.omit(Timediff))) %>%
  mutate(Work = ifelse(Work.Code == 2, "Fish", "Travel")) %>%
  group_by(Month, Ship, Work) %>%
  summarise(Total = sum(Total)) %>%
  pivot_wider(names_from = Work, values_from = Total) %>%
  mutate(ratio = Fish/Travel) %>%
  group_by(Month) %>%
  summarise(Mean = mean(ratio), SD = sd(ratio))
Combine <- bind_rows(Optimal, Timeindiff, NullResource, Null) %>%
  mutate(ratio = Fish.time/Travel.time) %>%
  group_by(Month, Model, Pay, Lag, Pun, Fish.Hr) %>%
  mutate(Z.score = case_when(Month == 4 ~ (ratio-VDR$Mean[2])/VDR$SD[2],
                             Month == 5 ~ (ratio-VDR$Mean[3])/VDR$SD[3],
                             Month == 6 ~ (ratio-VDR$Mean[4])/VDR$SD[4],
                             Month == 7 ~ (ratio-VDR$Mean[5])/VDR$SD[5],
                             Month == 8 ~ (ratio-VDR$Mean[6])/VDR$SD[6],
                             Month == 9 ~ (ratio-VDR$Mean[7])/VDR$SD[7])) %>%
  mutate(RMSE = case_when(Month == 4 ~ sqrt((ratio-VDR$Mean[2])^2),
                          Month == 5 ~ sqrt((ratio-VDR$Mean[3])^2),
                          Month == 6 ~ sqrt((ratio-VDR$Mean[4])^2),
                          Month == 7 ~ sqrt((ratio-VDR$Mean[5])^2),
                          Month == 8 ~ sqrt((ratio-VDR$Mean[6])^2),
                          Month == 9 ~ sqrt((ratio-VDR$Mean[7])^2)))
Combine %>%
  ggplot() +
  geom_boxplot(aes(x = Model, y = Z.score))
saveRDS(Combine, '/Users/judyhoho/Desktop/IBM/IBM_hexagon/simulation results/Organized results/FISH/Behavior_Combine.rds')


########################################################
# Fishing behavior (fishing area overlap)
## calculate area overlap (recommend to run in server, the scripts below is conducted in server)
########################################################

library(tidyverse)
library(sf)
library(parallel)
setwd('/home/pojuke/YunHo/IBM_hexagon')

Calculate.Overlap <- function(groups.IBM, VDR){
  All.Fish.Overlap <- data.frame()
  VDR.m <- VDR %>% filter(Month == groups.IBM$Month[1])
  groups.VDR <- split(VDR.m, list(VDR.m$Ship, VDR.m$Trip), drop = TRUE)
  for (k in 1:length(groups.VDR)){
    IBM.sf <- st_as_sf(groups.IBM, coords = c("X", "Y"), crs = 3857)
    VDR.sf <- groups.VDR[[k]] %>%
      filter(Work.Code == 2) %>% #filter fishing locations
      st_as_sf(coords = c("Longitude", "Latitude"), crs = 3857)
    IBM.poly <- IBM.sf %>%
      st_union() %>%
      st_convex_hull()
    VDR.poly <- VDR.sf %>%
      st_union() %>%
      st_convex_hull()
    poly_overlap <- st_intersection(IBM.poly, VDR.poly)
    if (st_intersects(IBM.poly, VDR.poly, sparse = FALSE) == F){
      J.index <- 0
    } else {
      poly_union <- st_union(IBM.poly, VDR.poly)
      J.index <- as.numeric(st_area(poly_overlap))/as.numeric(st_area(poly_union))}
    All.Fish.Overlap <- rbind(All.Fish.Overlap, data.frame(Month = groups.IBM$Month[1],
                                                           Pay = groups.IBM$Pay[1],
                                                           Lag = groups.IBM$Lag[1],
                                                           Pun = groups.IBM$Pun[1],
                                                           Fish.Hr = groups.IBM$Fish.Hr[1],
                                                           Ship = groups.IBM$Ship[1],
                                                           VDR = groups.VDR[[k]]$Ship[1],
                                                           S = groups.IBM$S[1],
                                                           J.index = J.index))
  }
  rm(VDR.m, groups.VDR, IBM.sf, VDR.sf, IBM.poly, VDR.poly, poly_overlap)
  return(All.Fish.Overlap)
}

load('OVERLAP/Angle_2021.RData')
groups.IBM <- readRDS('OVERLAP/Optimal_groups.rds')
n <- length(groups.IBM)
result <- do.call(rbind, mclapply(seq_along(groups.IBM), function(i) {
  if (i %% 10 == 0)
    cat(sprintf("Progress: %d/%d (%.1f%%)\n",i, n, 100*i/n))
  Calculate.Overlap(groups.IBM[[i]], VDR = All.VDR.TA)
}, mc.cores = 7))
result$Model <- "Optimal"
saveRDS(result, file = "OVERLAP/Overlap_Optimal.rds")

## organize the results calculated on server
Fish.Overlap <- data.frame()
aaa <- readRDS('/Users/judyhoho/Desktop/IBM/IBM_hexagon/simulation results/Organized results/OVERLAP/Overlap_Timeindiff.rds')
Fish.Overlap <- bind_rows(Fish.Overlap, aaa)
saveRDS(Fish.Overlap, file = '/Users/judyhoho/Desktop/IBM/IBM_hexagon/simulation results/Organized results/OVERLAP/Overlap_Combine.rds') #save results for all scenarios combined

## find the combination that generates highest Jaccard index for each iteration
library(gtools)
All.Overlap <- readRDS('/Users/judyhoho/Desktop/IBM/IBM_hexagon/simulation results/Organized results/OVERLAP/Overlap_Combine.rds')
Param.grid <- expand.grid(
  HR = c(0.5, 1, 1.5),
  Lag = c(4, 12, 24),
  Pun = c(0.2, 0.4, 0.6, 0.8, 1),
  Pay = c(0.2, 0.4, 0.6, 0.8, 1),
  S = 1:20,
  Month = 4:9)
IBM.ship <- 1:5
model.name <- "Optimal"
Find.Max <- function(Param, Overlap.df, IBM, model.n){
  Overlap <- Overlap.df %>%
    filter(Month == Param$Month & Pay == Param$Pay & Lag == Param$Lag & 
             Pun == Param$Pun & Fish.Hr == Param$HR & S == Param$S &
             Model == model.n)
  perm <- permutations(n = length(unique(Overlap$VDR)), r = length(IBM), 
                       v = unique(Overlap$VDR))
  all.comb <- do.call(rbind, lapply(seq_len(nrow(perm)), function(i) {
    data.frame(combination = i,
               IBM = IBM, VDR = perm[i, ])})) %>%
    left_join(Overlap, by = c("IBM" = "Ship", "VDR" = "VDR")) %>%
    group_by(combination) %>%
    mutate(J.mean = mean(J.index)) 
  return(max(all.comb$J.mean))
}

Overlap.comb <- data.frame()
for (i in 1:27){
  aaa <- Param.grid[((1000*(i-1)+1):(1000*i)),] %>%
    mutate(J.index = pmap_dbl(list(HR = HR,
                                   Lag = Lag,
                                   Pun = Pun,
                                   Pay = Pay,
                                   S = S,
                                   Month = Month),
                              function(HR, Lag, Pun, Pay, S, Month) {Param <- list(HR = HR,
                                                                                   Lag = Lag,
                                                                                   Pun = Pun,
                                                                                   Pay = Pay,
                                                                                   S = S,
                                                                                   Month = Month)
                              Find.Max(Param = Param,
                                       Overlap.df = All.Overlap,
                                       IBM = IBM.ship,
                                       model.n = model.name)}
    ))
  Overlap.comb <- bind_rows(Overlap.comb, aaa)
  print(i)
  rm(aaa)
}
saveRDS(Overlap.comb, file = paste0('/Users/judyhoho/Desktop/IBM/IBM_hexagon/simulation results/Organized results/OVERLAP/Overlap_comb_', model.name, '.rds')) #save results for each scenario


########################################################
# Compare the performance patterns based on shared parameter sets
## filter each metric by eliminating the worst 15th percentile of parameter values, and take the intersection of remaining parameter combinations after filtering all six metrics  
########################################################


## load in organized results for each evaluation pattern
## Frechet distance
All.Frechet <- readRDS('/Users/judyhoho/Desktop/IBM/IBM_hexagon/simulation results/Organized results/Frechet/Frechet_Combine.rds')
All.Frechet <- All.Frechet %>%
  group_by(Model, Pay, Lag, Pun, Fish.Hr) %>%
  summarise(Mean = mean(Frechet.d), SD = sd(Frechet.d))
## Fishing behaviors (time allocation)
All.Behavior <- readRDS('/Users/judyhoho/Desktop/IBM/IBM_hexagon/simulation results/Organized results/FISH/Behavior_Combine.rds')
All.Behavior <- All.Behavior %>%
  group_by(Model, Pay, Lag, Pun, Fish.Hr) %>%
  summarise(Mean = mean(Z.score), SD = sd(Z.score))
# Fishing behaviors (fishing location)
All.Overlap <- readRDS('/Users/judyhoho/Desktop/IBM/IBM_hexagon/simulation results/Organized results/OVERLAP/Overlap_comb.rds')
All.Overlap <- All.Overlap %>%
  group_by(Model, Pay, Lag, Pun, HR) %>%
  summarise(Mean = mean(J.index), SD = sd(J.index))
colnames(All.Overlap) <- c("Model", "Pay", "Lag", "Pun", "Fish.Hr", "Mean", "SD")
## Maximum distance from harbor
load('/Users/judyhoho/Desktop/IBM/IBM_hexagon/simulation results/Organized results/MAX/Max_All.RData')
All.Max <- Max.Z.score
rm(Max.Z.score)
All.Max <- All.Max %>%
  group_by(Model, Pay, Lag, Pun, Fish.Hr) %>%
  summarise(Mean = mean(Z.score), SD = sd(Z.score))
## Turning angle overlap
All.TA <- readRDS('/Users/judyhoho/Desktop/IBM/IBM_hexagon/simulation results/Organized results/TA/Overlap_Combine.rds')
All.TA <- All.TA %>%
  group_by(Model, Pay, Lag, Pun, Fish.Hr) %>%
  summarise(Mean = mean(Overlap), SD = sd(Overlap))
## Step length overlap
All.Step <- readRDS('/Users/judyhoho/Desktop/IBM/IBM_hexagon/simulation results/Organized results/STEP/Step_Combine.RData')
All.Step <- All.Step %>%
  group_by(Model, Pay, Lag, Pun, Fish.Hr) %>%
  summarise(Mean = mean(Dist.Overlap), SD = sd(Dist.Overlap))

##########################
# Filter to get shared parameter sets
filter.q <- 0.85 #filter quantile
shared <- All.Frechet %>%
  group_by(Model) %>%
  filter(Mean < quantile(Mean, probs = filter.q))

shared <- shared[, 1:5] %>%
  inner_join(All.Overlap, by = c('Model' = 'Model', 'Pay' = 'Pay', 
                                 'Lag' = 'Lag', 'Pun' = 'Pun', 'Fish.Hr' = 'Fish.Hr'))
shared <- shared %>%
  group_by(Model) %>%
  filter(Mean > quantile(Mean, probs = 1- filter.q))

shared <- shared[, 1:5] %>%
  inner_join(All.Behavior, by = c('Model' = 'Model', 'Pay' = 'Pay', 
                              'Lag' = 'Lag', 'Pun' = 'Pun', 'Fish.Hr' = 'Fish.Hr'))
shared <- shared %>%
  group_by(Model) %>%
  filter(abs(Mean) < quantile(abs(Mean), probs = filter.q))

shared <- shared[, 1:5] %>%
  inner_join(All.Max, by = c('Model' = 'Model', 'Pay' = 'Pay', 
                             'Lag' = 'Lag', 'Pun' = 'Pun', 'Fish.Hr' = 'Fish.Hr'))
shared <- shared %>%
  group_by(Model) %>%
  filter(abs(Mean) < quantile(abs(Mean), probs = filter.q))

shared <- shared[, 1:5] %>%
  inner_join(All.TA, by = c('Model' = 'Model', 'Pay' = 'Pay', 
                                   'Lag' = 'Lag', 'Pun' = 'Pun', 'Fish.Hr' = 'Fish.Hr'))
shared <- shared %>%
  group_by(Model) %>%
  filter(Mean > quantile(Mean, probs = 1- filter.q))

shared <- shared[,1:5] %>%
  inner_join(All.Step, by = c('Model' = 'Model', 'Pay' = 'Pay', 
                                 'Lag' = 'Lag', 'Pun' = 'Pun', 'Fish.Hr' = 'Fish.Hr'))
shared <- shared %>%
  group_by(Model) %>%
  filter(Mean > quantile(Mean, probs = 1- filter.q))

Optimal <- shared %>%
  filter(Model == 'Optimal')
NullR <- shared %>%
  filter(Model == 'NullResource')
NullB <- shared %>%
  filter(Model == 'Timeindiff')
Null <- shared %>%
  filter(Model == 'Null')

Final.shared <- inner_join(Optimal, NullR, Null, by = c('Pay' = 'Pay', 'Lag' = 'Lag', 
                                                 'Pun' = 'Pun', 'Fish.Hr' = 'Fish.Hr'))[,2:5]

Final.shared <- inner_join(Final.shared, NullB, by = c('Pay' = 'Pay', 'Lag' = 'Lag', 
                                         'Pun' = 'Pun', 'Fish.Hr' = 'Fish.Hr'))[,1:4]

# Plot all six performance metrics based on shared parameters
FR <- All.Frechet %>%
  inner_join(Final.shared, by = c('Pay' = 'Pay', 'Lag' = 'Lag', 
                           'Pun' = 'Pun', 'Fish.Hr' = 'Fish.Hr')) 
FR$Model <- factor(FR$Model, levels = c('Optimal', 'Timeindiff', 'NullResource', 'Null'))
FR.plot <- ggplot(FR) +
  geom_boxplot(aes(x = Model, y = Mean)) +
  geom_jitter(aes(x = Model, y = Mean), col = 'darkgrey', alpha = 1, size = 0.8) +
  scale_x_discrete(labels = str_wrap(c('Optimal', 'Time-indifferent', 'Random resource', 'Null'), width = 1)) +
  theme(legend.position = 'none') +
  scale_y_reverse("Mean(Frechet distance)")

FJ <- All.Overlap %>%
  inner_join(Final.shared, by = c('Pay' = 'Pay', 'Lag' = 'Lag', 
                                  'Pun' = 'Pun', 'Fish.Hr' = 'Fish.Hr')) 
FJ$Model <- factor(FJ$Model, levels = c('Optimal', 'Timeindiff', 'NullResource', 'Null'))
FJ.plot <- ggplot(FJ) +
  geom_boxplot(aes(x = Model, y = Mean)) +
  geom_jitter(aes(x = Model, y = Mean), col = 'darkgrey', alpha = 1, size = 0.8) +
  scale_x_discrete(labels = str_wrap(c('Optimal', 'Time- indifferent', 'Random resource', 'Null'), width = 5)) +
  theme(legend.position = 'none') +
  labs(y = "Mean(Jaccard index)")

FB <- All.Behavior %>%
  inner_join(Final.shared, by = c('Pay' = 'Pay', 'Lag' = 'Lag', 
                                  'Pun' = 'Pun', 'Fish.Hr' = 'Fish.Hr')) 
FB$Model <- factor(FB$Model, levels = c('Optimal', 'Timeindiff', 'NullResource', 'Null'))
FB.plot <- ggplot(FB) +
  geom_boxplot(aes(x = Model, y = Mean)) +
  geom_jitter(aes(x = Model, y = Mean), col = 'darkgrey', alpha = 1, size = 0.8) +
  scale_x_discrete(labels = str_wrap(c('Optimal', 'Time- indifferent', 'Random resource', 'Null'), width = 5)) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  theme(legend.position = 'none') +
  labs(y = "z(Time proportion)")

MAX <- All.Max %>%
  inner_join(Final.shared, by = c('Pay' = 'Pay', 'Lag' = 'Lag', 
                           'Pun' = 'Pun', 'Fish.Hr' = 'Fish.Hr')) 
MAX$Model <- factor(MAX$Model, levels = c('Optimal', 'Timeindiff', 'NullResource', 'Null'))
MAX.plot <- ggplot(MAX) +
  geom_boxplot(aes(x = Model, y = Mean)) +
  geom_jitter(aes(x = Model, y = Mean), col = 'darkgrey', alpha = 1, size = 0.8) +
  scale_x_discrete(labels = str_wrap(c('Optimal', 'Time- indifferent', 'Random resource', 'Null'), width = 5)) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  theme(legend.position = 'none') +
  labs(y = "z(Maximum distance from harbor)")

TA <- All.TA %>%
  inner_join(Final.shared, by = c('Pay' = 'Pay', 'Lag' = 'Lag', 
                           'Pun' = 'Pun', 'Fish.Hr' = 'Fish.Hr')) 
TA$Model <- factor(TA$Model, levels = c('Optimal', 'Timeindiff', 'NullResource', 'Null'))
TA.plot <- ggplot(TA) +
  geom_boxplot(aes(x = Model, y = Mean)) +
  geom_jitter(aes(x = Model, y = Mean), col = 'darkgrey', alpha = 1, size = 0.8) +
  scale_x_discrete(labels = str_wrap(c('Optimal', 'Time- indifferent', 'Random resource', 'Null'), width = 5)) +
  theme(legend.position = 'none') +
  labs(y = "Mean(Overlap of turning angle)")

STEP <- All.Step %>%
  inner_join(Final.shared, by = c('Pay' = 'Pay', 'Lag' = 'Lag', 
                           'Pun' = 'Pun', 'Fish.Hr' = 'Fish.Hr')) 
STEP$Model <- factor(STEP$Model, levels = c('Optimal', 'Timeindiff', 'NullResource', 'Null'))
STEP.plot <- ggplot(STEP) +
  geom_boxplot(aes(x = Model, y = Mean)) +
  geom_jitter(aes(x = Model, y = Mean), col = 'darkgrey', alpha = 1, size = 0.8) +
  scale_x_discrete(labels = str_wrap(c('Optimal', 'Time- indifferent', 'Random resource', 'Null'), width = 5)) +
  theme(legend.position = 'none') +
  labs(y = "Mean(Overlap of step length)")

library(patchwork)
patchwork <- FR.plot + FJ.plot + FB.plot + MAX.plot  + TA.plot + STEP.plot + plot_layout(nrow = 2) # + plot_spacer()
patchwork + plot_annotation(tag_levels = ('a'), tag_prefix = '(', tag_suffix = ')') & theme_bw() + theme(plot.tag = element_text(size = 12, face = "bold"), panel.grid = element_blank(), axis.title.y = element_text(size = 15, face = "bold"), axis.text.y = element_text(size = 13, face = "bold"), axis.text.x = element_text(size = 13, face = "bold", angle = 30, hjust = 1), axis.title.x = element_blank())


########################################################
# Compare the performance patterns based on respective best-performing parameter sets
## Filter parameter sets by average ranking of all six metrics, without order  
########################################################

All.TA <- All.TA %>%
  group_by(Model) %>%
  mutate(Rank = rank(-Mean))
All.Step <- All.Step %>%
  group_by(Model) %>%
  mutate(Rank = rank(-Mean))
All.Max <- All.Max %>%
  group_by(Model) %>%
  mutate(Rank = rank(abs(Mean)-0))
All.Behavior <- All.Behavior %>%
  group_by(Model) %>%
  mutate(Rank = rank(abs(Mean)-0))
All.Overlap <- All.Overlap %>%
  group_by(Model) %>%
  mutate(Rank = rank(Mean))
All.Frechet <- All.Frechet %>%
  group_by(Model) %>%
  mutate(Rank = rank(Mean))

Ranking <- inner_join(All.TA[,-7], All.Step[,-7], by = c('Model' = 'Model', 'Pay' = 'Pay', 'Lag' = 'Lag', 'Pun' = 'Pun', 'Fish.Hr' = 'Fish.Hr'), keep = FALSE)
Ranking <- inner_join(Ranking, All.Max[,-7], by = c('Model' = 'Model', 'Pay' = 'Pay', 'Lag' = 'Lag', 'Pun' = 'Pun', 'Fish.Hr' = 'Fish.Hr'), keep = FALSE)
Ranking <- inner_join(Ranking, All.Behavior[,-7], by = c('Model' = 'Model', 'Pay' = 'Pay', 'Lag' = 'Lag', 'Pun' = 'Pun', 'Fish.Hr' = 'Fish.Hr'), keep = FALSE)  
Ranking <- inner_join(Ranking, All.Overlap[,-7], by = c('Model' = 'Model', 'Pay' = 'Pay', 'Lag' = 'Lag', 'Pun' = 'Pun', 'Fish.Hr' = 'Fish.Hr'), keep = FALSE) 
Ranking <- inner_join(Ranking, All.Frechet[,-7], by = c('Model' = 'Model', 'Pay' = 'Pay', 'Lag' = 'Lag', 'Pun' = 'Pun', 'Fish.Hr' = 'Fish.Hr'), keep = FALSE)  
Ranking$Final <- rowMeans(Ranking[, c(7,9,11,13,15,17)])
Ranking <- Ranking %>%
  group_by(Model) %>%
  slice_min(Final, n = 23)
colnames(Ranking) <- c('Model', 'Pay', 'Lag', 'Pun', 'Fish.Hr', 'Mean.TA', 'Rank.TA', 'Mean.Step', 'Rank.Step', 'Mean.Max', 'Rank.Max', 'Mean.Behavior', 'Rank.Behavior', 'Mean.Overlap', 'Rank.Overlap', 'Mean.Frechet', 'Rank.Frechet', 'Rank.Final')

# Plot all six performance metrics based on respective best-performing parameters
FR <- All.Frechet %>%
  inner_join(Ranking, by = c('Pay' = 'Pay', 'Lag' = 'Lag', 
                           'Pun' = 'Pun', 'Fish.Hr' = 'Fish.Hr', 'Model' = 'Model')) 
FR$Model <- factor(FR$Model, levels = c('Optimal', 'Timeindiff', 'NullResource', 'Null'))
FR.plot <- ggplot(FR) +
  geom_boxplot(aes(x = Model, y = Mean)) +
  geom_jitter(aes(x = Model, y = Mean), col = 'darkgrey', alpha = 1, size = 0.8) +
  scale_x_discrete(labels = str_wrap(c('Optimal', 'Time- indifferent', 'Random resource', 'Null'), width = 1)) +
  scale_y_reverse("Mean(Frechet distance)") +
  theme(legend.position = 'none') 

FJ <- All.Overlap %>%
  inner_join(Ranking, by = c('Pay' = 'Pay', 'Lag' = 'Lag', 
                                  'Pun' = 'Pun', 'Fish.Hr' = 'Fish.Hr', 'Model' = 'Model')) 
FJ$Model <- factor(FJ$Model, levels = c('Optimal', 'Timeindiff', 'NullResource', 'Null'))
FJ.plot <- ggplot(FJ) +
  geom_boxplot(aes(x = Model, y = Mean)) +
  geom_jitter(aes(x = Model, y = Mean), col = 'darkgrey', alpha = 1, size = 0.8) +
  scale_x_discrete(labels = str_wrap(c('Optimal', 'Time- indifferent', 'Random resource', 'Null'), width = 5)) +
  theme(legend.position = 'none') +
  labs(y = "Mean(Jaccard index)")

FB <- All.Behavior %>%
  inner_join(Ranking, by = c('Pay' = 'Pay', 'Lag' = 'Lag', 
                                  'Pun' = 'Pun', 'Fish.Hr' = 'Fish.Hr', 'Model' = 'Model')) 
FB$Model <- factor(FB$Model, levels = c('Optimal', 'Timeindiff', 'NullResource', 'Null'))
FB.plot <- ggplot(FB) +
  geom_boxplot(aes(x = Model, y = Mean)) +
  geom_jitter(aes(x = Model, y = Mean), col = 'darkgrey', alpha = 1, size = 0.8) +
  scale_x_discrete(labels = str_wrap(c('Optimal', 'Time- indifferent', 'Random resource', 'Null'), width = 5)) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  theme(legend.position = 'none') +
  labs(y = "z(Time proportion)")

MAX <- All.Max %>%
  inner_join(Ranking, by = c('Pay' = 'Pay', 'Lag' = 'Lag', 
                           'Pun' = 'Pun', 'Fish.Hr' = 'Fish.Hr', 'Model' = 'Model')) 
MAX$Model <- factor(MAX$Model, levels = c('Optimal', 'Timeindiff', 'NullResource', 'Null'))
MAX.plot <- ggplot(MAX) +
  geom_boxplot(aes(x = Model, y = Mean)) +
  geom_jitter(aes(x = Model, y = Mean), col = 'darkgrey', alpha = 1, size = 0.8) +
  scale_x_discrete(labels = str_wrap(c('Optimal', 'Time- indifferent', 'Random resource', 'Null'), width = 5)) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  theme(legend.position = 'none') +
  labs(y = "z(Maximum distance from harbor)")

TA <- All.TA %>%
  inner_join(Ranking, by = c('Pay' = 'Pay', 'Lag' = 'Lag', 
                           'Pun' = 'Pun', 'Fish.Hr' = 'Fish.Hr', 'Model' = 'Model')) 
TA$Model <- factor(TA$Model, levels = c('Optimal', 'Timeindiff', 'NullResource', 'Null'))
TA.plot <- ggplot(TA) +
  geom_boxplot(aes(x = Model, y = Mean)) +
  geom_jitter(aes(x = Model, y = Mean), col = 'darkgrey', alpha = 1, size = 0.8) +
  scale_x_discrete(labels = str_wrap(c('Optimal', 'Time- indifferent', 'Random resource', 'Null'), width = 5)) +
  theme(legend.position = 'none') +
  labs(y = "Mean(Overlap of turning angle)")

STEP <- All.Step %>%
  inner_join(Ranking, by = c('Pay' = 'Pay', 'Lag' = 'Lag', 
                           'Pun' = 'Pun', 'Fish.Hr' = 'Fish.Hr', 'Model' = 'Model')) 
STEP$Model <- factor(STEP$Model, levels = c('Optimal', 'Timeindiff', 'NullResource', 'Null'))
STEP.plot <- ggplot(STEP) +
  geom_boxplot(aes(x = Model, y = Mean)) +
  geom_jitter(aes(x = Model, y = Mean), col = 'darkgrey', alpha = 1, size = 0.8) +
  scale_x_discrete(labels = str_wrap(c('Optimal', 'Time- indifferent', 'Random resource', 'Null'), width = 5)) +
  theme(legend.position = 'none') +
  labs(y = "Mean(Ovrelap of step length)")

patchwork <- FR.plot + FJ.plot + FB.plot + MAX.plot + TA.plot + STEP.plot + plot_layout(nrow = 2)
patchwork + plot_annotation(tag_levels = ('a'), tag_prefix = '(', tag_suffix = ')') & theme_bw() + theme(plot.tag = element_text(size = 12, face = "bold"), panel.grid = element_blank(), axis.title.y = element_text(size = 15, face = "bold"), axis.text.y = element_text(size = 13, face = "bold"), axis.text.x = element_text(size = 13, face = "bold", angle = 30, hjust = 1), axis.title.x = element_blank())


########################################################
# PCA of best-performing parameter values among models (respective best-performing parameter sets)
library(plotly)
library(ggfortify)

colnames(Ranking)[1] <- "Scenario"
pca_res <- prcomp(Ranking[,2:5], scale. = TRUE)
autoplot(pca_res, data = Ranking, colour = 'Scenario', loadings = TRUE, loadings.label = TRUE,
         loadings.color = "black", loadings.label.color = "black",
         loadings.label.hjust = 1.2, loadings.label.vjust = 1.2, size = 4, loadings.label.size = 8) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(size = 22, face = 'bold'),
        legend.text = element_text(size = 20, face = "bold"),
        axis.title = element_text(size = 18, face = "bold"),
        axis.text = element_text(size = 16, face = "bold")) 
ggsave(filename = paste0('PCA_ranking.png'),
       path = '/Users/judyhoho/Desktop/IBM/IBM_hexagon/figures', width = 12, height = 10)

pca_res <- prcomp(shared[,2:5], scale. = TRUE)
autoplot(pca_res, data = shared, colour = 'Model', loadings = TRUE, loadings.label = TRUE,
         loadings.color = "black", loadings.label.color = "black",
         loadings.label.hjust = 1.2, loadings.label.vjust = 1.2) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "bottom") 
