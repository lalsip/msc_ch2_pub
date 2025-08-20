# Libraries --------------------

library(tidyverse)
#install.packages('tidyverse', dependencies = TRUE)

library(lubridate) # working with dates
library(hms) # working with times
library(measurements) #deg min sec to dd

library(patchwork) # easier side by side plots
library(ggpubr) # side by side with one legend
library(lattice) # bwplot for homogeneity
library(GGally) # ggpairs for collinearity
library(gstat) # semivariogram for independence

library(sf) # need for reading spatial data, work with spatial vector data (replaces sp)
library(sp) # old version of sf
#library(rgdal)
# 8.15.25 rgdal was removed from cran in 2023-06
# 8.15.25 dld last version (1.6-7) from archive https://cran.r-project.org/src/contrib/Archive/rgdal/
# remotes::install_local(path = './data/rgdal_1.6-7.tar.gz', dependencies = TRUE)
#library(maptools)

library(raster)
library(spatstat)

library(lme4) # glmer.nb (mixed model)
library(MASS) # glm.nb
library(mgcv) # gam and gamm
library(gstat)

library(vegan)  #install.packages("vegan", dependencies = TRUE) 
library(lattice)
source("./data/HighstatLibV13.R")

# Map --------------------

# all functions and methods in sf that operate on spatial data are prefixed by st_ for spatial type

# GSL outline
sf<-st_read("./data/sf/NWT_Detailed_LakesRivers_GSL.shp")
# Reading layer ... using driver `ESRI Shapefile'
# Simple feature collection with 1 feature and 11 fields
# geometry type:  POLYGON
# dimension:      XY
# bbox:           xmin: -116.9824 ymin: 60.82129 xmax: -108.9008 ymax: 62.95337
# geographic CRS: NAD83

# if need to change what part of lake zoomed in on, do here - only some coord work
#ggplot()+
#  geom_sf(data=sf)+
#  coord_sf(xlim=c(-117, -112.5), ylim=c(60.75, 62.55), # bbox
#           expand=FALSE)

# GSL Management zones
mz<-st_read("./data/sf/GSL_managementzones.shp")
# Reading layer ... using driver `ESRI Shapefile'
# Simple feature collection with 5 features and 1 field
# geometry type:  LINESTRING
# dimension:      XY
# bbox:           xmin: -116.1181 ymin: 60.87688 xmax: -113.4238 ymax: 62.02824
# geographic CRS: WGS 84

# GSL rivers/inflows/outflows (file is for canada)
rv<-st_read("./data/sf/Med_Canada_Rivers.shp")
# Reading layer ... using driver `ESRI Shapefile'
# Simple feature collection with 5625 features and 11 fields
# geometry type:  LINESTRING
# dimension:      XY
# bbox:           xmin: -141 ymin: 42.31844 xmax: -54.0504 ymax: 82.79323
# geographic CRS: WGS 84  

# note different CRS: mz and rv use WGS 84, but sf uses NAD 83
# transform CRS for sf to WGS84
sf_WGS84 <- st_transform(sf, 4326) #EPSG 4326 is WGS84
mz_NAD83 <- st_transform(mz, 4269) # NAD 83

#https://stackoverflow.com/questions/60656445/how-to-fix-degree-symbol-not-showing-correctly-in-r-on-linux-fedora-31/60733863#60733863

## create base maps
GSL <- ggplot() + 
  #geom_sf(data = sf) + 
  geom_sf(data = sf_WGS84, colour = "grey") + # lake outline
  geom_sf(data = mz, colour = "Red",size = 0.5)+ # management zones
 geom_sf(data= rv, colour = "light grey") + # inflows
  coord_sf(xlim=c(-117.5, -109), ylim=c(60.5, 63), # bbox
           expand=FALSE) + #expand=FALSE gets rid of gaps on edges
  scale_x_continuous(labels = function(x) paste0(x, '\u00B0')) +
  scale_y_continuous(labels = function(x) paste0(x, '\u00B0')) +
  theme(
    panel.background = element_blank(), # Remove panel background
    axis.line = element_line(colour = "grey"), # Add axis lines
  )

GSLmb <- ggplot() + 
  geom_sf(data=sf_WGS84, colour="grey") +
  geom_sf(data = mz, colour = "Red",size = 0.5)+ # management zones
  geom_sf(data= rv, colour = "light grey") + # inflows
  coord_sf(xlim=c(-117, -112.5), ylim=c(60.75, 62.55), # bbox is different - no east arm
           expand=FALSE)+ # no gaps on axes
  scale_x_continuous(labels = function(x) paste0(x, '\u00B0')) +
  scale_y_continuous(labels = function(x) paste0(x, '\u00B0')) +
  theme(
    panel.background = element_blank(), # Remove panel background
    axis.line = element_line(colour = "grey"), # Add axis lines
  )

# western basin 
GSLwb <- ggplot() + 
  geom_sf(data=sf_WGS84, colour="grey") +
  geom_sf(data = mz, colour = "Red",size = 0.5)+ # management zones
  #geom_sf(data= rv, colour = "light grey") + # inflows
  coord_sf(xlim=c(-117, -115), ylim=c(60.75, 61.8), # bbox is different - no east arm
           expand=FALSE)+ # no gaps on axes
  scale_x_continuous(labels = function(x) paste0(x, '\u00B0')) +
  scale_y_continuous(labels = function(x) paste0(x, '\u00B0')) +
  theme( #axis.text = element_blank(),
    panel.background = element_blank(), # Remove panel background
    axis.line = element_line(colour = "grey")) # Add axis lines


# Boat catch + lim ----------
boat_raw <- read.csv("./data/boatcatch.csv")

# add in NAs where nothing
# https://stackoverflow.com/questions/51214357/how-to-remove-only-rows-that-have-all-na-in-r
boat_raw[boat_raw == ""] <- NA # define NA pattern
# removes rows with no data in them (all col in row NA) - import issue
boat_raw <- boat_raw[rowSums(is.na(boat_raw)) != ncol(boat_raw),] 
# removes cols with no data in them - import issue
boat_raw <- boat_raw[colSums(is.na(boat_raw)) != nrow(boat_raw)]

# check for NA
colSums(is.na(boat_raw))

# working copy
boat <- boat_raw %>% 
  dplyr::select(-note, -notes_dataentry, -scanned, -netID, -crew, -endlat_decdeg, -endlon_decdeg, -endsitedep_m)

# tidy
# convert var into correct formats
boat <- dplyr::rename(boat, site=grid)

boat$site <- as.factor(boat$site)
boat$spp <- as.factor(boat$spp)
boat$gillnet <- as.factor(boat$gillnet) %>% 
  recode_factor("benthic"="bottom")
boat$mesh_mm <- factor(boat$mesh_mm, 
                       levels=c("13", "19", "25", "32", "38", "51", "64", "89", "114", "133", "140"), 
                       ordered=T)
boat$effset <- as.factor(boat$effset)
boat$layer <- as.factor(boat$layer) %>% 
  recode_factor("hypo"="Hypo")

boat$DO_mgL <- as.numeric(boat$DO_mgL)
boat$startlat_decdeg <- as.numeric(boat$startlat_decdeg)
boat$startlon_decdeg <- -boat$startlon_decdeg
boat$totnum <- as.numeric(boat$totnum) # will need to change to numeric once confirm the 'lots' etc numbers

boat$setdt <- paste(boat$setdate, boat$settime, sep=" ") 
boat$setdt <- mdy_hm(boat$setdt) #convert to POSIXct format 
boat$settime <- hms::as_hms(boat$setdt) #time to hms
boat$setdate <- mdy(boat$setdate)

boat$liftdt <- paste(boat$liftdate, boat$lifttime, sep=" ") 
boat$liftdt <- mdy_hm(boat$liftdt) #convert to POSIXct format 
boat$lifttime <- hms::as_hms(boat$liftdt) #time to hms
boat$liftdate <- mdy(boat$liftdate)
boat$setting <- as.factor(boat$setting)

# calculate soak time
boat<- boat %>% 
  mutate(soaktime= as.numeric(difftime(lubridate::ymd_hms(liftdt), 
                                       lubridate::ymd_hms(setdt), units="hours")))

boat$year <- as.factor(year(boat$setdate))

# GSLwb +
#   geom_point(data=boat %>% filter(year==2023),
#              mapping=aes(x=startlon_decdeg, y=startlat_decdeg, colour=gillnet),
#              position=position_dodge(width=0.1), size=3) +
#   xlab(NULL) +
#   ylab(NULL) +
#   theme(legend.position="bottom")

boat <- boat %>% 
  filter(gillnet != "coney") %>% 
  filter(mesh_mm != "133") %>% 
  droplevels()

# need to add in spp that had 0 caught
GSLspp <- data.frame(
  spp = c("LKWF", "LKT", "INCO", "BBT", "LNSK", "RDWF", "LKH", "SJCK", "LCK", "TP", 
          "AG", "NSSBK", "SPT", "NPK", "GDI", "NPD", "LKDS", "WY", "CHUB", "SPSP", 
          "SAUGER", "ALP", "WSK","SHSP", "SLSP"), 
  totnum=rep(0, 25))

mesh_mm=data.frame(mesh_mm=c("13", "19", "25", "32", "38", "51", "64", "89", "114", "140"))

GSLspp <- merge(GSLspp, mesh_mm, by=NULL)
GSLspp$mesh_mm <- factor(GSLspp$mesh_mm, 
                         levels=c("13", "19", "25", "32", "38", "51", "64", "89", "114",  "140"), 
                         ordered=T)
GSLspp$spp <- as.factor(GSLspp$spp)
# filter for LKWF to go faster
GSLspp <- GSLspp %>% 
  filter(spp == "LKWF")

# get a df with the visited sites
boat1 <- boat %>% 
  group_by(setdate, site, gillnet) %>% 
  #filter(spp=="LKWF") %>% 
  group_by()

boatdates <- boat1 %>% 
  group_by (setdate, site, gillnet) %>% 
  dplyr::select(setdate, site, gillnet, setting_m) %>% 
  unique()

# datespp should be all the dates sampled with all the zeros
datespp <- merge(boatdates, GSLspp, by=NULL)

datespp <- datespp %>% 
  arrange(setdate, site, gillnet, mesh_mm)
datespp <- datespp %>% 
  filter(!is.na(mesh_mm))
datespp <- datespp %>% 
  fill(setting_m, .direction=c("down"))
datespp <- datespp %>% 
  unique()

# actual data with zeros added in
boat1_zeros <- right_join(boat1, datespp, by=c("setdate", "gillnet", 
                                               "setting_m", "spp", "site", "mesh_mm"))

boat1_zeros <- boat1_zeros %>% 
  mutate(totnum.x = ifelse(is.na(totnum.x), totnum.y, totnum.x)) 

# this should be the complete dataset with all zeros added 
boat1_zeros <- boat1_zeros %>% 
  group_by(setdate, site, gillnet) %>% 
  fill(settime, liftdate, lifttime, soaktime, startlon_decdeg, startlat_decdeg, startsitedep_m, setdt, 
       liftdt, year, effset, dep_m, temp_C, layer, DO_mgL, DO_perc, pH, turb_FNU, bot, setting, setting_m,
       .direction=c("downup")) %>% 
  group_by() %>% 
  dplyr::select(-totnum.y)


boat1_zeros <- dplyr::rename(boat1_zeros, totnum.boat=totnum.x)

boat1_zeros <- boat1_zeros %>% 
  filter(setting == "bot") %>% 
  filter(!is.na(temp_C)) %>% # removes those w/o YSI
  droplevels()

boat1_zeros <- boat1_zeros %>% 
  arrange(liftdate, site, mesh_mm)

boat1_zeros <- boat1_zeros %>% 
  filter(!is.na(mesh_mm)) %>% 
  unique()

boat1_zeros <- boat1_zeros %>% 
  dplyr::select(-spp, -panelarea_m2, -bot)

boat1_zeros <- boat1_zeros %>% 
  filter(effset == '1')
boat1_zeros <- boat1_zeros %>% 
  filter(gillnet != "pelagic")

boat1_zeros <- boat1_zeros %>% 
  dplyr::select(-lifttime, -settime, -setdt, -liftdt) %>% 
  unique()

# GSLmb +
#   geom_point(aes(x=boat$startlon_decdeg, y=boat$startlat_decdeg))


# Land catch --------------
# catch that was processed

land_raw <- read.csv("./data/proccatch.csv")

# add in NAs where nothing
# https://stackoverflow.com/questions/51214357/how-to-remove-only-rows-that-have-all-na-in-r
land_raw[land_raw == ""] <- NA # define NA pattern
# removes rows with no data in them (all col in row NA) - import issue
land_raw <- land_raw[rowSums(is.na(land_raw)) != ncol(land_raw),] 
# removes cols with no data in them - import issue
land_raw <- land_raw[colSums(is.na(land_raw)) != nrow(land_raw)]

# check for NA
colSums(is.na(land_raw))

# working copy
land <- land_raw %>% 
  dplyr::select(-note, -notes_dataentry, -X, -sampler)

# tidy
# convert var into correct formats
land <- dplyr::rename(land, site=grid)

land$site <- as.factor(land$site)
land$spp <- as.factor(land$spp)
land$gillnet <- as.factor(land$gillnet) %>% recode_factor("benthic"="bottom")
land <- land %>% 
  filter(gillnet != "coney") %>% 
  filter(mesh_mm != "133") %>% 
  droplevels()

land$mesh_mm <- factor(land$mesh_mm, 
                       levels=c("13", "19", "25", "32", "38", "51", "64", "89", "114", "140"), 
                       ordered=T)
land$totnum <- as.numeric(land$totnum) # 'lots' now NA
land$totwg_g <- as.numeric(land$totwg_g)

land$setdate <- mdy(land$setdate)
land$liftdate <- mdy(land$liftdate)

# add zeros
# get a df with the visited sites
land1 <- land %>% 
  group_by(setdate, site, gillnet) %>% 
  filter(spp=="LKWF") %>% 
  group_by()
landdates <- land1 %>% 
  group_by (setdate, site, gillnet) %>% 
  dplyr::select(setdate, site, gillnet, setting_m) %>% 
  unique()

# datespp should be all the dates sampled with all the zeros
datespp_land <- merge(landdates, GSLspp, by=NULL)

datespp_land <- datespp_land %>% 
  filter(!is.na(mesh_mm))

datespp_land <- datespp_land %>% 
  arrange(setdate, site, gillnet, mesh_mm)

# actual data with zeros added in
land1_zeros <- right_join(land1, datespp_land, by=c("setdate", "gillnet", "setting_m", "spp", "site", "mesh_mm"))

land1_zeros <-land1_zeros %>% 
  mutate(totnum.x= ifelse(is.na(totnum.x), totnum.y, totnum.x)) %>% 
  mutate(totwg_g = ifelse(is.na(totwg_g), 0, totwg_g))

# this should be the complete dataset with all zeros added 
land1_zeros <- land1_zeros %>% 
  group_by(setdate, site, gillnet) %>% 
  fill(liftdate, .direction=c("down")) %>% 
  group_by() %>% 
  dplyr::select(-totnum.y)

land1_zeros <- dplyr::rename(land1_zeros, totnum.land=totnum.x)

land1_zeros <- land1_zeros %>% 
  filter(!is.na(mesh_mm))

land1_zeros <- land1_zeros %>% 
  arrange(setdate, site, gillnet, mesh_mm) %>% 
  dplyr::select(-spp, -setting_m)

land1_zeros <- land1_zeros %>% 
  filter(gillnet != "pelagic")

# Combine Boat + lim and Land -----

catch <- left_join(boat1_zeros, land1_zeros, by=c("site", "setdate", "gillnet", "mesh_mm"))


# bathy data all FIGS grids 
bathyFIGS <- read.csv("./data/allcoord.csv") #xlsx file is GSL2016allgridcoord

# rename to shorter name
bathyFIGS$z<- as.numeric(bathyFIGS$Grid_M)
bathyFIGS$z <- -bathyFIGS$z # depths should be negative

bathyFIGS$y<- as.numeric(bathyFIGS$LatN)
bathyFIGS$x<- as.numeric(bathyFIGS$LonW)

a <- bathyFIGS %>% 
  filter(Area == "IE") 

mean(a$z)

# add FMAs
FIGS <- bathyFIGS %>% 
  dplyr::select(Grid, Area) %>% 
  dplyr::rename("site"="Grid", "FMA"="Area")
FIGS$site <- as.factor(FIGS$site)
FIGS$FMA <- as.factor(FIGS$FMA) 
FIGS$FMA <-  recode_factor(FIGS$FMA, "IW"="W", "IE"="E", "II"="2", "III"="3", "IV"="4", "V"="5")

catch <- left_join(catch, FIGS)

catchwb <- catch %>% 
  filter(FMA =="W" | FMA=="E")

# filtering out NAs so that window for dist to shore works
catchwb <- catchwb %>% 
  filter(!is.na(catchwb$startlat_decdeg))

# Distance from nearest shore -------
# 8.15.25 this used rgdal::readOGR(), updated to use sf pkg
# make a shapefile of lake, just 2 zones for now

# ZONE_IW<-readOGR(dsn='C:/Users/alsipl/Desktop/JackDistFromShore', layer='ZONE IW')
# ZONE_II<-readOGR(dsn='C:/Users/alsipl/Desktop/JackDistFromShore', layer='ZONE II')
# ZONE_III<-readOGR(dsn='C:/Users/alsipl/Desktop/JackDistFromShore', layer='ZONE III')
# crs(ZONE_II)<-'+proj=longlat +datum=WGS84 +no_defs'
# ZONE_II.UTM=spTransform(ZONE_II, geo_proj)
# crs(ZONE_III)<-'+proj=longlat +datum=WGS84 +no_defs'
# ZONE_III.UTM=spTransform(ZONE_III, geo_proj)
geo_proj='+proj=utm +zone=11' # 8.20.25 make sure consistent zone as in catchwb$latlon

ZONE_IE <- sf::st_read("./data/sf/ZONE IE.shp")
sf::st_crs(ZONE_IE) <- '+proj=longlat +datum=WGS84 +no_defs'
ZONE_IE.UTM = sf::st_transform(ZONE_IE, geo_proj)

ZONE_IW <- sf::st_read("./data/sf/ZONE IW.shp")
sf::st_crs(ZONE_IW) <- '+proj=longlat +datum=WGS84 +no_defs'
ZONE_IW.UTM = sf::st_transform(ZONE_IW, geo_proj)

# 8.20.25 added in FMAs II and III so boundary of IE not incl
ZONE_II <- sf::st_read("./data/sf/ZONE II.shp")
sf::st_crs(ZONE_II) <- '+proj=longlat +datum=WGS84 +no_defs'
ZONE_II.UTM = sf::st_transform(ZONE_II, geo_proj)

ZONE_III <- sf::st_read("./data/sf/ZONE III.shp")
sf::st_crs(ZONE_III) <- '+proj=longlat +datum=WGS84 +no_defs'
ZONE_III.UTM = sf::st_transform(ZONE_III, geo_proj)

#merge zones
## helps avoid measuring to area boundary rather than shore
FISH.SHP.UTM=sf::st_union(ZONE_IE.UTM, ZONE_IW.UTM)

# plot(FISH.SHP.UTM)
FISH.SHP.UTM=sf::st_union(FISH.SHP.UTM,ZONE_II.UTM)
FISH.SHP.UTM=sf::st_union(FISH.SHP.UTM,ZONE_III.UTM)

# citation()
# citation('spatstat.geom')

catchwb$latlon <- terra::project(x=cbind(c(catchwb$startlon_decdeg), c(catchwb$startlat_decdeg)), 
                                 from="+proj=longlat +datum=WGS84 +no_defs", #8.15.25 from= was missing
                                 to='+proj=utm +zone=11') #8.15.25 JH had as 14, should be 11
#create window from lake
WIN.UTM<- spatstat.geom::as.owin(FISH.SHP.UTM)
# plot(WIN.UTM)

##create ppp object in lake window
p<-spatstat.geom::ppp(catchwb$latlon[,1], catchwb$latlon[,2], window = WIN.UTM)
# 8.20.25 ignore warning about duplicated points

# calculate distance of point to boundary of window
d<-spatstat.geom::bdist.points(p) 

# add distance to shore to dataset
catchwb$distance_to_shore <- d

#double check that it worked
GSLmb +
  geom_point(aes(x=startlon_decdeg, y=startlat_decdeg, 
                 colour=distance_to_shore/1000), #/1000 for km
             data=catchwb)
# 8.20.25 makes sense now that FMAs II and III incl 

# not working? make sure none of coord are missing or negative
## this should be giving distance to ANY shore

# https://stackoverflow.com/questions/26321614/obtaining-the-subset-of-points-which-are-outside-a-specified-irregular-polygon-u#:~:text=If%20your%20main%20goal%20is%20to%20%22find%20out,outside%20the%20window%20you%20specified%20when%20creating%20test.ppp.
#exterior_points <- attr(p, "rejects")


# Distance from southern shore -----
## JH i made a nonsense one by eye, you'll need to make a line shape object from your actual zones for this to work
south<-sf::st_read(dsn='./data/sf/southsouth.shp')
sf::st_crs(south) <-'+proj=longlat +datum=WGS84 +no_defs'
# south.UTM = sf::st_transform(south, geo_proj)
south <- as_Spatial(south) #convert to spatial counterpart # 8.20.25 needed now that not using readOGR
# south2 <- as(south, "Spatial") # this also works for simple features

pnts<-SpatialPointsDataFrame(catchwb[,c("startlon_decdeg","startlat_decdeg")], catchwb,
                               proj4string = CRS('+proj=longlat +datum=WGS84 +no_defs'))


dist.mat <- geosphere::dist2Line(p = pnts, line = south)

# add on dist to S shore to dataset
catchwb <- cbind(catchwb, dist.mat) 

catchwb <- catchwb %>% 
  dplyr::select(-lon, -lat, -ID, -latlon) %>% 
  dplyr::rename("distance_to_south" = "distance")

#double check that it worked
GSLmb +
  geom_point(aes(x=startlon_decdeg, y=startlat_decdeg, 
                 colour=distance_to_south/1000), #/1000 for km
             data=catchwb)
# 8.20.25 works! or at least numbers make sense

# filter for VALID sets
catchwb <- catchwb %>% 
  filter(soaktime <= 30) %>% 
  filter(soaktime >=18) %>% 
  filter(effset==1)

table(catchwb$site, catchwb$year)

# create factors for general visuals
# catchwb <- catchwb %>%
#   mutate(fturb_FNU = ifelse(turb_FNU == 0, "0",
#                             ifelse(turb_FNU <=5, "<5",
#                                    ifelse(turb_FNU <=10, "<10",
#                                           ifelse(turb_FNU <=25, "<25",
#                                                  ifelse(turb_FNU <=50, "<50",
#                                                         ifelse(turb_FNU <=100, "<100", ">100")))))))
# catchwb$fturb_FNU <- factor(catchwb$fturb_FNU, levels=c("0", "<5", "<10", "<25", "<50", "<100", ">100"))



# drop coney
catchwb <- catchwb %>% 
  filter(year != "2022") %>% 
  filter(year != "2023") %>% 
  filter(mesh_mm != "133") 

# drop mids
catchwb <- catchwb %>% 
  filter(setting != "p10") %>% 
  filter(setting != "p20") %>% 
  dplyr::select(-setting)
  
# drop unused levels
catchwb$gillnet <- droplevels(catchwb$gillnet)
catchwb$mesh_mm <- droplevels(catchwb$mesh_mm)
#catchwb$setting <- droplevels(catchwb$setting)


catchwb$site <- droplevels(catchwb$site)
catchwb$FMA <- droplevels(catchwb$FMA)

catchwb <- catchwb %>%
  filter(turb_FNU <= 80)
# don't think it was calibrated past that

# create object with only relevant covariates and response var
catchwb2 <- catchwb %>% 
  dplyr::select(site, gillnet, soaktime, startlat_decdeg, startlon_decdeg, startsitedep_m, mesh_mm, totnum.boat, 
                temp_C, layer, turb_FNU, DO_perc, pH, year, distance_to_shore, distance_to_south, liftdate.x, 
                totnum.land, totwg_g)

catchwb2 <- catchwb2 %>% 
 dplyr::rename("lat" = "startlat_decdeg", "lon" = "startlon_decdeg") %>% 
  dplyr::rename("sitedep_m" = "startsitedep_m") %>% 
  dplyr::rename("fgillnet" = "gillnet") %>% 
  dplyr::rename("fyear" = "year")

# catchwb2 <- catchwb2 %>% 
#   filter(setting != "p10") # shouldn't have been set there
# catchwb2$setting <- droplevels(catchwb2$setting)

LongLatToUTM <- function(x,y,zone, Hemisphere = "north"){
  xy <- data.frame(ID = 1:length(x), X = x, Y = y)
  coordinates(xy) <- c("X", "Y")
  proj4string(xy) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")  
  if (Hemisphere == "north"){
    res <- spTransform(xy, 
                       CRS(paste("+proj=utm +zone=",zone," ellps=WGS84",sep='')))
  }
  
  if (Hemisphere == "south"){
    res <- spTransform(xy, 
                       CRS(paste("+proj=utm +zone=",zone,
                                 " +ellps=WGS84 +datum=WGS84 +units=m +no_defs +south",sep='')))
  }
  
  return(as.data.frame(res))
}

LL <- LongLatToUTM(x = catchwb2$lon, 
                   y = catchwb2$lat, 
                   zone = 11,
                   Hemisphere = "south")

catchwb2$X.utm <- LL$coords.x1 #8.20.25 had to changed from LL$X to LL$coords.x1
catchwb2$Y.utm <- LL$coords.x2 #same as above
catchwb2$Xkm   <- catchwb2$X.utm / 1000 
catchwb2$Ykm   <- catchwb2$Y.utm / 1000

catchwb2 <- catchwb2 %>% 
  dplyr::select(-X.utm, -Y.utm)



# Do we have any NAs?
colSums(is.na(catchwb2))

table(catchwb2$fyear, catchwb2$site)

# what if I don't want to deal with mesh size
catchwb4 <- catchwb2 %>% 
  dplyr::group_by(fyear, site) %>% 
  dplyr::mutate(nettot = sum(totnum.boat)) %>% 
  dplyr::mutate(netarea_m2 = ifelse(fgillnet=="pelagic", 600, 300))  # benthic 300 m2, pel 600
catchwb4$netarea_m2 <- factor(catchwb4$netarea_m2, levels=c("600", "300"))

catchwb4 <- catchwb4 %>% 
  dplyr::group_by(fyear, site) %>% 
  dplyr::mutate(totnum.land = ifelse(is.na(totnum.land), 0, totnum.land), 
                totwg_g = ifelse(is.na(totwg_g), 0, totwg_g))

catchwb4 <- catchwb4 %>% 
  dplyr::group_by(fyear, site) %>% 
  dplyr::mutate(landtot = sum(totnum.land)) %>% 
  dplyr::mutate(landwt_g = sum(totwg_g)) %>% 
  dplyr::mutate(avgwt_g = landwt_g / landtot)

catchwb4 <- catchwb4 %>% 
  dplyr::select(-mesh_mm, -totnum.boat, -totnum.land, -totwg_g) %>% 
  ungroup() %>% 
  unique()

table(catchwb4$fyear, catchwb4$site)

catchwb5 <- catchwb4

catchwb5$year <- as.numeric(catchwb5$fyear)

# Do we have any NAs?
colSums(is.na(catchwb5))

 
# won't exclude yet since need to pinpoint var
dim(catchwb2)

table(catchwb5$year, catchwb5$netarea_m2) # totally off balance - gotta drop 600
catchwb5 <- catchwb5 %>% 
  filter(netarea_m2 == "300") %>% 
  droplevels
table(catchwb5$year, catchwb5$netarea_m2)

table(catchwb5$site)
table(catchwb5$year)
# 1  2  3  4  5  6  7  8 
# 9  6  9 15 23 15 16 18 


# standardize each continuous covariate. ----
 catchwb5 <- catchwb5 %>%
  dplyr::filter(!is.na(temp_C)) %>% 
  dplyr::filter(!is.na(turb_FNU)) #%>% 
#  dplyr::filter(!is.na(DO_perc))
catchwb5$soaktimec <- MyStd(catchwb5$soaktime)
catchwb5$turb_FNUc <- MyStd(catchwb5$turb_FNU)
catchwb5$temp_Cc <- MyStd(catchwb5$temp_C)
#catchwb5$setting_mc <- MyStd(catchwb5$setting_m)
#catchwb5$sitedep_mc<- MyStd(catchwb5$sitedep_m)
catchwb5$distance_to_southc <- MyStd(catchwb5$distance_to_south)
catchwb5$distance_to_shorec <- MyStd(catchwb5$distance_to_shore)
#catchwb5$DO_percc <- MyStd(catchwb5$DO_perc)

colSums(is.na(catchwb5))

catchwb5$year <- as.numeric(catchwb5$fyear)
catchwb5 <- catchwb5 %>% 
   dplyr::select(-netarea_m2)



# subset --------

# nrow(catchwb5) # 111
# length(unique(catchwb5$site)) # 32 sites
# 0.75*32 # 24 - train
# 0.25*32 # 8 - test
# 
# set.seed(1234)
# catchwb5_train <- catchwb5 %>%  
#   group_by(site) %>% 
#   dplyr::summarise() %>% 
#   dplyr::sample_n(24) 
# 
# catchwb5_train <- catchwb5 %>% 
#   dplyr::right_join(catchwb5_train, by='site')
# 
# catchwb5_test <- catchwb5 %>%  
#   dplyr::anti_join(catchwb5_train)
# 
# table(catchwb5_train$site)
# table(catchwb5_test$site)
# 
# catchwb5_test <- catchwb5_test %>% 
#   dplyr::rename('resp_obs' = 'nettot')  %>% 
#   add_column(nettot=NA)
# 
# catchwb5_pred <- catchwb5_train %>% 
#   full_join(catchwb5_test) # now the full dataset where test data has NA as nettot 

# num sites ------
# 9 sites (of 32) were only sampled once 
# remove those sites? 

catchwb6 <- catchwb5 %>% 
  group_by(site) %>% 
  dplyr::mutate(freq=n()) %>% 
  ungroup() %>% 
  dplyr::filter(freq > 1) %>% 
  dplyr::select(-freq)

catchwb6$site <- droplevels(catchwb6$site)
# now have a dataset with sites that were visited multiple times only
nrow(catchwb6) # 103
length(unique(catchwb6$site)) # 24 sites
0.75*24 # 18 - train
0.25*24 # 6 - test

set.seed(1234)
catchwb6_train <- catchwb6 %>%  
  group_by(site) %>% 
  dplyr::summarise() %>% 
  dplyr::sample_n(18) 

catchwb6_train <- catchwb6 %>% 
  dplyr::right_join(catchwb6_train, by='site')

catchwb6_test <- catchwb6 %>%  
  dplyr::anti_join(catchwb6_train)

table(catchwb6_train$site)
table(catchwb6_test$site, catchwb6_test$year)
table(catchwb6_test$site)

catchwb6_test <- catchwb6_test %>% 
  dplyr::rename('resp_obsN' = 'nettot')  %>% 
  add_column(nettot=NA) %>% 
  dplyr::rename('resp_obsB' = 'avgwt_g') %>% 
  add_column(avgwt_g = NA)

catchwb6_pred <- catchwb6_train %>% 
  full_join(catchwb6_test) # now the full dataset where test data has NA as nettot 


#https://stackoverflow.com/questions/36068963/r-how-to-split-a-data-frame-into-training-validation-and-test-sets


