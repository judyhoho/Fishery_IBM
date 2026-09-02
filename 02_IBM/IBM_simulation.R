########################################################
# IBM (hexagon version)
########################################################
## (1) The code is divided into 3 parts: load required packages and data, functions definition, and run simulation in parallel 
## (2) Functions included: Decision (determine movement for next step), Renew (resource replenish), func_fishery (simulation for all vessels to finish one trip)
## (3) All model scenarios are included, with two switch one controlling the decision rule and the other controlling the landscape setting

####################################################
# Set working directory and load required packages 
library(tidyverse)
library(sf)
library(terra)
library(parallel)
setwd('/home/pojuke/YunHo/IBM_hexagon')

suppressPackageStartupMessages({ #used to suppress messages in log file
  library(tidyverse)
  library(sf)
  library(terra)
  library(parallel)
})

####################################################
# Load in data needed
## the boundary for Northern Taiwan and China 
Taiwan <- read.delim('data/Taiwan_24.txt', sep = '')
Taiwan <- Taiwan[-c(731:733),] #eliminate data points that would make the polygon not applicable
Taiwan <- rbind(Taiwan, Taiwan[1,]) #make the polygon a closed one
island <- st_sf(geometry = st_sfc(st_polygon(list(as.matrix(Taiwan[, c("lon", "lat")]))), crs = 4326))
island <- island %>% 
  st_transform( , crs = 3826) #polygon
China <- read.delim('data/china_2.txt', sep = '')
China <- China[-c(2753),] 
China <- rbind(China, China[1,]) #make the polygon a closed one
island.china <- st_sf(geometry = st_sfc(st_polygon(list(as.matrix(China[, c("lon", "lat")]))), crs = 4326))
island.china <- island.china %>%
  st_transform( , crs = 3826) #polygon

## 杜氏鎖管分區圖 for the area to give punishment to 
Duvauceli <- read.csv('data/Duv_squid.csv', stringsAsFactors = FALSE)
Duvauceli <- st_as_sf(Duvauceli, wkt = "wkt_geom", crs = 4326) 
Duvauceli <- st_transform(Duvauceli, 3826)

## Load NeighborRef for the data frame of neighboring grids 
load('data/NeighborRef.RData')

## Load in intersect reference
load('data/Ref_Hexagon.RData')

## Set global parameters 
Travel.Velocity <- 0.15
Harbour.Position <- c(121.6333, 25.33270) #the centroid of the nearest grid


####################################################
# Decision function 
####################################################
Decision <- function(Abun.Map = Squid.Abun.Map,
                     Velo = Travel.Velocity,
                     Prop.Capt = 0.5, 
                     Cumulative.Capt = 0,
                     Current.Hr = 0,
                     Current.D = 1,
                     Current.Pos = Harbour.Position,
                     Harbour.Pos = Harbour.Position,
                     Light.cost = 1500, # NT dollar/hr
                     Payoff.pun = 0.3,
                     Target.D = 10,
                     Next.Fish.Time = 0,
                     D.boundary = Duvauceli,
                     Duv.pun = 0.5,
                     Squid.P = 200, # NT dollar/kg
                     Oil.cost = 73.8, # NT dollar/km
                     Fish.hour = 1,
                     Inter.ref = Intersect.ref,
                     MVT = TRUE,
                     Bias.sto = FALSE,
                     Landing.sto = FALSE){
  
  # Calculate payoff map
  Dist.Pay.Map <- Abun.Map
  Current.pos.sf <- st_as_sf(data.frame(Lon = Current.Pos[1], Lat = Current.Pos[2]), coords = c("Lon", "Lat"),
                             crs = 4326) %>%
    st_transform(st_crs(Abun.Map)) %>%
    st_within(Abun.Map)
  Ref <- Inter.ref[,Current.pos.sf[[1]]]
  Dist.Pay.Map <- Dist.Pay.Map[Ref == 0,] #eliminate grids that will cross the land
  Dist.Pay.Map$Distance <- sqrt((Dist.Pay.Map$Lon-Current.Pos[1])^2 + (Dist.Pay.Map$Lat-Current.Pos[2])^2)
  Dist.Pay.Map$Travel.Time <- Dist.Pay.Map$Distance/Velo
  if (MVT){ #TRUE for Optimal/Resource-randomized scenario, FALSE for Time-indifferent/Null scenario
    Dist.Pay.Map$Payoff <- ((Abun.Map$Abun[Ref == 0] * Prop.Capt * Squid.P) - (Oil.cost * Dist.Pay.Map$Distance * 111) - (Fish.hour * Light.cost))/(Fish.hour + Dist.Pay.Map$Travel.Time) 
  } else{Dist.Pay.Map$Payoff <- ((Abun.Map$Abun[Ref == 0] * Prop.Capt * Squid.P) - (Oil.cost * Dist.Pay.Map$Distance * 111) - (Fish.hour * Light.cost))}
  
  Travel.Hour <- Current.Hr + Dist.Pay.Map$Travel.Time
  Dist.Pay.Map$Payoff <- ifelse(Travel.Hour < 2 | Travel.Hour > 14,
                                Dist.Pay.Map$Payoff * Payoff.pun,
                                Dist.Pay.Map$Payoff) #lower payoff outside operation time 
  inside_region <- st_intersects(Dist.Pay.Map$geometry, D.boundary$wkt_geom)
  inside_region <- lengths(inside_region) > 0
  Dist.Pay.Map$Payoff[which(inside_region == T)] <- Dist.Pay.Map$Payoff[which(inside_region == T)] * Duv.pun #lower payoff in the Duv squid area
  Dist.Pay.Map <- Dist.Pay.Map[!is.na(Dist.Pay.Map$Payoff),]
  
  # Stochasticity for payoff (not currently used)
  if (Bias.sto){
    for (j in 1:nrow(Dist.Pay.Map)){
      if (Dist.Pay.Map$Payoff[j] > 0){
        Dist.Pay.Map$Payoff[j] <- runif(1, min = Dist.Pay.Map$Payoff[j] * 0.9,
                                        max = Dist.Pay.Map$Payoff[j] * 1.1)
      }
    }
  }
  
  # Identify the optimal grid
  ## if all payoff is 0 or current day exceeds target day, go home.
  if(all(Dist.Pay.Map$Payoff == 0 | Current.D >= Target.D)){ 
    Best.Grid <- Harbour.Pos; Day.End <- TRUE
    Travel.T <- 
      sqrt((Best.Grid[1]-Current.Pos[1])^2 + (Best.Grid[2]-Current.Pos[2])^2)/Velo
  } else {
    ## determine the best grid (called rough initially for landing stochasticity)
    Best.Rough.Grid <- Dist.Pay.Map[which(Dist.Pay.Map$Payoff == max(Dist.Pay.Map$Payoff)),]
    if(nrow(Best.Rough.Grid) > 1) { 
      Dist.Harbour <- c()
      for (i in 1:nrow(Best.Rough.Grid)){
        Dist.Harbour[i] <- (Best.Rough.Grid$Lon[i] - Current.Pos[1])^2 + (Best.Rough.Grid$Lat[2] - Current.Pos[2])^2
      }
      ## if more than one grid are with highest payoff, choose the closet one 
      Closer <- Best.Rough.Grid[which(Dist.Harbour == min(Dist.Harbour)),]
      Best.Rough.Grid <- Closer
      if(nrow(Closer) > 1){
        ## if they are of the same distance, randomly choose one
        Chosen <- sample(1:nrow(Closer), 1)
        Best.Rough.Grid <- Closer[Chosen, ]
      }
    }
    
    # Stochasticity for landing (not currently used)
    if (Landing.sto){
      ## randomly choose one precise grid within all neighboring grids
      Order <- which(Best.Rough.Grid$Lon == Abun.Map$Lon & Best.Rough.Grid$Lat == Abun.Map$Lat)
      Nneighbor <- length(Abun.Map$NeighborLon[[Order]])
      Chosen <- sample(1:Nneighbor, 1)
      Best.Grid <- c(Abun.Map$NeighborLon[[Order]][Chosen], 
                     Abun.Map$NeighborLat[[Order]][Chosen])
    } else{Best.Grid <- c(Best.Rough.Grid$Lon, Best.Rough.Grid$Lat)}
    
    Day.End <- FALSE
    Travel.T <- 
      sqrt((Best.Grid[1]-Current.Pos[1])^2 + (Best.Grid[2]-Current.Pos[2])^2)/Velo
  } # next target grid is determined as Gest.Grid
  
  # Decide the action for next step
  if(Next.Fish.Time == 1 & Current.Hr >= 2 & Current.Hr <= 14){ #if stay and fish
    Best.Grid <- Current.Pos; Travel.T <- 0
    Cumulative.Capt <- Cumulative.Capt + Abun.Map$Abun[which(as.character(Abun.Map$Lon) == as.character(Current.Pos[1]) & as.character(Abun.Map$Lat) == as.character(Current.Pos[2]))] * Prop.Capt * Fish.hour
    Fish.Time <- Fish.hour; Next.Fish.Time <- 0 #reset Next.Fish.Time
  } else if(all(!Day.End & (Current.Hr+Travel.T >= 2) & (Current.Hr+Travel.T <= 14))){ #if next step would be fishing
    Next.Fish.Time <- 1; Fish.Time <- 0
  } else if(all(Best.Grid == Current.Pos & (Current.Hr <= 2 | Current.Hr >= 14))){ #if already at the best grid but cannot fish due to operating time limit
    Travel.T <- 1; Next.Fish.Time <- 0; Fish.Time <- 0
  } else{Next.Fish.Time <- 0; Fish.Time <- 0} #if next step is not fishing
  
  return(list(Next.Pos = Best.Grid,
              Travel.Time = Travel.T,
              Fish.Time = Fish.Time,
              Next.Fish.Time = Next.Fish.Time,
              Cumu.Capt = Cumulative.Capt,
              Day.End = Day.End))
} # End of the Decision function. 


####################################################
# Renew function 
####################################################
Renew <- function(Deplete.Day = 1,
                  Abun.Map = Squid.Abun.Map,
                  Deplete.list = Deplete,
                  Neighbor.df = neighbor.df){
  
  Deplete.Day <- Deplete.Day + 1
  Deplete.list <- Deplete.list %>%
    group_by(Lon, Lat) %>%
    distinct()
  Deplete.list <- cbind(Deplete.list, 
                        Order = sample(1:nrow(Deplete.list), nrow(Deplete.list),
                                       replace = FALSE)) #assign a random renew order
  Deplete.list <- Deplete.list[order(Deplete.list$Order),]
  for (j in 1:nrow(Deplete.list)){
    ## set the grid that is to be renewed and its neighboring grids
    D.grid <- st_as_sf(data.frame(Lon = Deplete.list$Lon[j], Lat = Deplete.list$Lat[j]), coords = c("Lon", "Lat"), crs = 4326) %>%
      st_transform(crs = 3826) %>% 
      st_within(Abun.Map) %>%
      unlist()
    D.neighborLon <- Neighbor.df$Lon[Abun.Map$neighbors[D.grid][[1]]]
    D.neighborLat <- Neighbor.df$Lat[Abun.Map$neighbors[D.grid][[1]]]
    for (k in 1:length(D.neighborLat)){
      Prob.renew <- as.logical(rbinom(1, size = 1, prob = 0.5)) #50% chance for renew
      if(Prob.renew){
        Prop.renew <- runif(1, 0, 0.1) #proportion of renew is between 0 and 0.1
        R.grid <- st_as_sf(data.frame(Lon = D.neighborLon[k], Lat = D.neighborLat[k]), coords = c("Lon", "Lat"), crs = 4326) %>%
          st_transform(crs = 3826) %>% 
          st_within(Abun.Map) %>%
          unlist()
        Abun.Map$Abun[D.grid] = Abun.Map$Abun[D.grid] + Abun.Map$Abun[R.grid] * Prop.renew
        Abun.Map$Abun[R.grid] = Abun.Map$Abun[R.grid] - Abun.Map$Abun[R.grid] * Prop.renew
      }
    }
  }
  
  ## after renewal, reset the Deplete.list
  Deplete.list <- data.frame(matrix(ncol = 2,nrow = 0, 
                                    dimnames=list(NULL, c('Lon', 'Lat'))))
  return(list(Update.Deplete.D = Deplete.Day,
              Update.Abun = Abun.Map,
              Update.Deplete = Deplete.list))
} # End of the Renew function


####################################################
# Simulation function 
####################################################
# Function for running simulations
func_fishery <- function(month.selected = 6,
                         MVT = TRUE,
                         Resource = TRUE,
                         Pay = 0.2,
                         HR = 0.5,
                         Pun = 0.2,
                         Lag = 0,
                         S = 1){
  
  # Generate the resource distribution map
  ## average from maps of 5 random days in that month, and check whether abundance prediction exist for that month
  Month_list <- c('Jan', 'Feb', 'Mar', 'Apr', 'May', 'June', 'July', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec')
  Background.list <- c('01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30')
  Background.check <- file.exists(paste0('data/predict_2021/Pred_20210', month.selected, Background.list, '.txt')) #prediction is missing for certain days
  Background.list <- Background.list[Background.check] 
  
  Random.five <- sample(Background.list, 5)
  Squid.Abun.Map <- data.frame()
  for (b in Random.five){
    cwb <- as_tibble(read.delim(paste0('data/predict_2021/Pred_20210', month.selected, b,  '.txt'), sep = '')[,1:5])
    cwb <- cwb[!is.na(cwb$Pred_gam),]
    cwb <- cwb %>%
      filter(Lon >= 120, Lon <= 126, Lat >= 24, Lat <= 30) %>%
      filter(!(Lon >= 123 & Lon <= 126 & Lat >= 24 & Lat <= 25.7)) %>%
      filter(Bathy > 0 & Bathy < 300)
    df_sf_proj <- st_as_sf(cwb, coords = c("Lon", "Lat"), crs = 4326) %>%
      st_transform(, crs = 3826)  
    coords <- st_coordinates(df_sf_proj)
    r_ext <- ext(min(coords[,1]), max(coords[,1]), min(coords[,2]), max(coords[,2]))
    res_value <- 1000 
    r <- rast(extent = r_ext, resolution = res_value)
    r <- rasterize(vect(df_sf_proj),  r, field = "Pred_gam", fun = mean)
    r_sf <- st_as_sf(as.points(r, values = TRUE))
    hex_grid <- st_make_grid(r_sf, cellsize = 10000, square = FALSE)
    hex_sf <- st_sf(geometry = hex_grid)
    hex_values <- st_join(hex_sf, r_sf) %>%
      group_by(geometry) %>%
      summarise(value = mean(mean, na.rm = TRUE)) %>%
      drop_na(value)
    st_crs(hex_values) <- 3826
    idx <- st_intersects(hex_values, island) #this returns intersection between all points of hexagon and the point layer
    hex_values <- hex_values[lengths(idx) == 0,] #keep those hexagons that do not overlap with any points
    idx <- st_intersects(hex_values, island.china) 
    hex_values <- hex_values[lengths(idx) == 0,]
    hex_centroids_ll <- st_transform(st_centroid(hex_values$geometry), crs = 4326)
    hex_values$Lon <- st_coordinates(hex_centroids_ll)[,1]
    hex_values$Lat <- st_coordinates(hex_centroids_ll)[,2]
    Squid.Abun.Map <- rbind(Squid.Abun.Map, hex_values[, c(2,3,4,1)])
    rm(cwb, df_sf_proj, r, r_sf, hex_sf, hex_values)
  }
  colnames(Squid.Abun.Map) <- c('Abun', 'Lon', 'Lat', 'geometry')
  Squid.Abun.Map <- Squid.Abun.Map %>%
    group_by(Lon, Lat) %>%
    summarise(Abun = mean(Abun)) %>%
    ungroup() 
  Squid.Abun.Map <- Squid.Abun.Map %>%
    st_join(neighbor.df, join = st_equals) 
  Squid.Abun.Map <- Squid.Abun.Map[,c(1,2,3,4,8)]
  colnames(Squid.Abun.Map) <- c("Lon", "Lat", "Abun", "geometry", "neighbors")
  
  if (Resource == FALSE){#Smooth random resource, for Resource-randomized and Null scenario
    Squid.Abun.Map$Abun_ini <- Squid.Abun.Map$Abun
    Reorder = sample(1:nrow(Squid.Abun.Map), nrow(Squid.Abun.Map))
    Squid.Abun.Map$shuffle <- Squid.Abun.Map$Abun[Reorder]
    hex_coords <- st_transform(Squid.Abun.Map$geometry, crs = 4326)
    hex_coords <- st_coordinates(st_centroid(hex_coords))
    theta <- runif(1, 0.1, 0.5)  # smoothing distance in same units as coordinates
    Squid.Abun.Map$smoothed <- sapply(seq_len(nrow(Squid.Abun.Map)), function(i) {
      dists <- sqrt((hex_coords[,1] - hex_coords[i,1])^2 + (hex_coords[,2] - hex_coords[i,2])^2)
      w <- exp(-0.5 * (dists/theta)^2)
      sum(w * Squid.Abun.Map$shuffle, na.rm=TRUE)/sum(w, na.rm=TRUE)
    })
    Squid.Abun.Map$Abun <- Squid.Abun.Map$smoothed * (sum(Squid.Abun.Map$Abun_ini)/sum(Squid.Abun.Map$smoothed))
  }
  
  # Initial states set up
  Travel.Velocity = 0.15
  Prop.Capture = 0.2
  Current.Abun = Squid.Abun.Map
  Initial.Abun = Current.Abun
  Current.Hour = 0
  Cumulative.Hour = 0
  Current.Position = Harbour.Position 
  Cumulative.Gain = 0
  Cumulative.Distance = 0
  Day = 1
  Target.Day = 10
  Oil.Cost = 73.8
  Light.Oil.Cost = 1500
  Squid.Price = 200
  Payoff.Punishment = Pay
  Fish.Time = 0
  Next.Fish.Time = 0
  Fish.Hour = HR
  Duv.boundary = Duvauceli
  Duv.Punishment = Pun
  Ref.Inter = Intersect.ref
  Day.End = FALSE
  No.Ship = 5
  Leave.Order = 1:5 
  Time.Lag = Lag
  All.Ship = list()
  
  for (i in 1:No.Ship){
    df = data.frame(Day = Day + (Current.Hour + (which(Leave.Order == i)-1) * Time.Lag) %/% 24,
                    Hour = (Current.Hour + (which(Leave.Order == i)-1) * Time.Lag) %% 24,
                    Cumu.Hour = Cumulative.Hour + (which(Leave.Order == i)-1) * Time.Lag,
                    Pos.X = Current.Position[1],
                    Pos.Y = Current.Position[2],
                    Cumu.Gain = Cumulative.Gain,
                    Cumu.Dist = Cumulative.Distance,
                    Fish = Fish.Time,
                    Next.Fish = Next.Fish.Time,
                    End = Day.End)
    All.Ship[[i]] = df
  }
  
  Deplete.D = 1
  Deplete = data.frame(matrix(ncol = 2,nrow = 0, dimnames=list(NULL, c('Lon', 'Lat')))) 
  # End of initial set up
  
  # Simulation begins
  repeat{
    # Decide who's term it is now
    Ships.End = c()
    for (i in 1:length(All.Ship)){
      Ships.End = c(Ships.End, tail(All.Ship[[i]]$End, 1))
    }
    
    Ships.Time = c()
    for (i in 1:length(All.Ship)){
      Ships.Time = c(Ships.Time, tail(All.Ship[[i]]$Cumu.Hour, 1))
    }
    Ships.Time = Ships.Time[!Ships.End]
    
    Term = which(!Ships.End)[which(Ships.Time == min(Ships.Time))]
    if (length(Term) > 1){
      Ship.Move = sample(Term, 1)
    } else {Ship.Move = Term}
    
    # The selected ship makes next movement  
    if (tail(All.Ship[[Ship.Move]]$End, 1) == FALSE){
      Current.Hour = tail(All.Ship[[Ship.Move]]$Hour, 1)
      Cumulative.Hour = tail(All.Ship[[Ship.Move]]$Cumu.Hour, 1)
      Current.Position = c(tail(All.Ship[[Ship.Move]]$Pos.X, 1), 
                           tail(All.Ship[[Ship.Move]]$Pos.Y, 1))
      Cumulative.Gain = tail(All.Ship[[Ship.Move]]$Cumu.Gain, 1)
      Cumulative.Distance = tail(All.Ship[[Ship.Move]]$Cumu.Dist, 1)
      Day = tail(All.Ship[[Ship.Move]]$Day, 1)
      Fish.Time = tail(All.Ship[[Ship.Move]]$Fish, 1)
      Next.Fish.Time = tail(All.Ship[[Ship.Move]]$Next.Fish, 1)
      
      if (Current.Hour >= 14 & Current.Hour <= 22){
        Current.Hour = (Current.Hour + 1) %% 24
        Cumulative.Hour = Cumulative.Hour + 1
        Fish.Time = 0
        Day.End = FALSE
      } else {
        # Make the decision of next step.
        Next.Step = Decision(Abun.Map = Current.Abun,
                             Velo = Travel.Velocity,
                             Prop.Capt = Prop.Capture,
                             Cumulative.Capt = Cumulative.Gain,
                             Current.Hr = Current.Hour,
                             Current.D = Day,
                             Current.Pos = Current.Position,
                             Harbour.Pos = Harbour.Position,
                             Light.cost = Light.Oil.Cost,
                             Payoff.pun = Payoff.Punishment,
                             Target.D = Target.Day,
                             Next.Fish.Time = Next.Fish.Time,
                             D.boundary = Duv.boundary,
                             Duv.pun = Duv.Punishment,
                             Squid.P = Squid.Price,
                             Oil.cost = Oil.Cost,
                             Fish.hour = Fish.Hour,
                             Inter.ref = Ref.Inter,
                             MVT = MVT,
                             Bias.sto = FALSE,
                             Landing.sto = FALSE)
        
        # Perform the decided next step.
        Cumulative.Distance = Cumulative.Distance + sqrt((Next.Step$Next.Pos[1] - Current.Position[1])^2 + (Next.Step$Next.Pos[2] - Current.Position[2])^2)
        Current.Position = Next.Step$Next.Pos
        Day = Day + (Current.Hour + Next.Step$Travel.Time + Next.Step$Fish.Time) %/% 24
        Current.Hour = (Current.Hour + Next.Step$Travel.Time + Next.Step$Fish.Time) %% 24
        Cumulative.Hour = Cumulative.Hour + Next.Step$Travel.Time + Next.Step$Fish.Time
        Next.Fish.Time = Next.Step$Next.Fish.Time
        Fish.Time = Next.Step$Fish.Time
        Day.End = Next.Step$Day.End
        
        if (Next.Step$Fish.Time == Fish.Hour){
          Cumulative.Gain = Next.Step$Cumu.Capt
          Spot = which(as.character(Current.Abun$Lon) == as.character(Current.Position[1]) & as.character(Current.Abun$Lat) == as.character(Current.Position[2]))
          Current.Abun$Abun[Spot] = (1 - Prop.Capture * Next.Step$Fish.Time) * Current.Abun$Abun[Spot] 
          Deplete = rbind(Deplete, data.frame(Lon = Current.Position[1], 
                                              Lat = Current.Position[2]))
          
        }
      } # end of deciding next step
      
      All.Ship[[Ship.Move]] = rbind(All.Ship[[Ship.Move]], 
                                    data.frame(Day = Day,
                                               Hour = Current.Hour,
                                               Cumu.Hour = Cumulative.Hour,
                                               Pos.X = Current.Position[1],
                                               Pos.Y = Current.Position[2],
                                               Cumu.Gain = Cumulative.Gain,
                                               Cumu.Dist = Cumulative.Distance,
                                               Fish = Fish.Time,
                                               Next.Fish = Next.Fish.Time,
                                               End = Day.End)) 
    } # end of selected ship 
    
    
    # Renew resource
    Ship.Morning = c()
    for (i in 1:length(All.Ship)){
      Ship.Morning = c(Ship.Morning, tail(All.Ship[[i]]$Hour, 1) >= 14 & tail(All.Ship[[i]]$Hour, 1) <= 22)
    }
    Ships.End = c()
    for (i in 1:length(All.Ship)){
      Ships.End = c(Ships.End, tail(All.Ship[[i]]$End, 1))
    }
    
    Out.boat = ((24*(Deplete.D - 1)+ 22) %/% Time.Lag) + 1
    if (Out.boat <= length(Ship.Morning)){
      Ship.Morning = Ship.Morning[1:Out.boat]
      Ships.End = Ships.End[1:Out.boat]
    }
    
    if (all(Day == Deplete.D & Ship.Morning & Ships.End != TRUE & nrow(Deplete) > 0)){ #once a day in the morning
      Renew.squid = Renew(Deplete.Day = Deplete.D, Abun.Map = Current.Abun, 
                          Deplete.list = Deplete, Neighbor.df = neighbor.df)
      Deplete.D = Renew.squid$Update.Deplete.D
      Current.Abun = Renew.squid$Update.Abun
      Deplete = Renew.squid$Update.Deplete}

    if(all(Ships.End)) break
  } # end of simulation
  
  for (j in 1:length(All.Ship)){ #record and combine records
    All.Ship[[j]]$Pay <- rep(Payoff.Punishment, nrow(All.Ship[[j]]))
    All.Ship[[j]]$Lag <- rep(Time.Lag, nrow(All.Ship[[j]]))
    All.Ship[[j]]$Pun <- rep(Duv.Punishment, nrow(All.Ship[[j]]))
    All.Ship[[j]]$Fish.Hr <- rep(Fish.Hour, nrow(All.Ship[[j]]))
    All.Ship[[j]]$Ship <- rep(j, nrow(All.Ship[[j]]))
    All.Ship[[j]]$S <- rep(S, nrow(All.Ship[[j]]))
    All.Ship[[j]]$Theta <- rep(theta, nrow(All.Ship[[j]]))
  }
  Combine = data.frame()
  for (k in 1:5){
    Combine = rbind(Combine, as.data.frame(All.Ship[[k]]))
  }
  return(list(Initial.Abun, Combine))
}

####################################################
# Set parameters for simulations
##run in parallel
####################################################
Param.grid <- expand.grid(
  HR = c(0.5, 1, 1.5),
  Lag = c(4, 12, 24),
  Pun = c(0.2, 0.4, 0.6, 0.8, 1),
  Pay = c(0.2, 0.4, 0.6, 0.8, 1),
  S = 1:20,
  Month.m = 4:9,
  MVT = TRUE, #behavior rule
  Resource = TRUE #resource landscape
)

run_simulation <- function(params) {
  Simu.result <- func_fishery(
    month.selected = params$Month.m,
    Pay = params$Pay,
    HR  = params$HR,
    Pun = params$Pun,
    Lag = params$Lag,
    S   = params$S,
    MVT = MVT,
    Resource = Resource
  )
  return(list(Params = params,
              fishers = Simu.result[[2]],
              initial.map = Simu.result[[1]]))
}

####################################################
# Run parallel on server
####################################################
n.cores   <- 7
chunk_sz  <- 1000
n.jobs    <- nrow(Param.grid)
out_dir   <- "/home/pojuke/YunHo/IBM_hexagon"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

chunk_id <- 1
buffer   <- list()

for (start in seq(1, n.jobs, by = chunk_sz)) {
  
  end <- min(start + chunk_sz - 1, n.jobs)
  idx <- start:end
  cat(sprintf("Running jobs %d–%d at %s\n", start, end, Sys.time()))
  
  results_chunk <- mclapply(
    idx,
    function(i) run_simulation(Param.grid[i, ]),
    mc.cores = n.cores
  )
  buffer <- c(buffer, results_chunk)
  
  ## export the results for every 1000 iterations
  saveRDS(
    buffer,
    file = file.path(out_dir, paste0("Null_chunk_", chunk_id, ".rds"))
  )
  cat(sprintf("Saved chunk %d at %s\n", chunk_id, Sys.time()))
  buffer   <- list()   #clear memory
  chunk_id <- chunk_id + 1
}


