# Y = # LKWF individuals (netcount) -> count data so Poisson or NB
# X = 
# 
# - soak time (NEED for effort)
# - net area (NEED for effort since mix of P and B nets)
# - year
# - turb
# - Dist from S shore
# - Dist from nearest shore
# - site dep 
# - temp


#### WESTERN BASIN ONLY #######
# no coney net #ten panel multi mesh only

## libraries ---------
source("./data/HighstatLibV13.R") # changed
#load("./data/datacleanoutput.RData") #clean data from dataclean2023.R

library(ggplot2)
library(tidyverse)

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

Mode(catchwb6$nettot)

# Does catch (totnum) differ with year and spatial location?
# Is there a relationship of catch and turbidity or thermal strata?


mean(catchwb6_train$soaktime) # 22.95
sd(catchwb6_train$soaktime) # 1.39


## experimental design ------
catchwb6_pred <- catchwb6_pred %>% 
  mutate(set = ifelse(!is.na(nettot), "Training", "Testing"))


GSLwb +
  geom_point(aes(x=lon, y=lat, colour=set), size=4, alpha=0.4, data=catchwb6_pred) +
  labs(x="Longitude", y="Latitude", colour="Dataset") +
    #facet_wrap(~fyear) +
  theme(axis.text = element_text(size=9), 
       panel.spacing = unit(20, "pt"), 
       legend.position = c(0.92, 0.9) )


citation('INLA')

table(catchwb6_train$fyear) # bottom-set nets (incl b)
# 2012 2013 2014 2015 2016 2017 2018 2019 ## catchwb6_train (only repeated sites)
#  6    5    7    9   14   10   13   12 

# outliers ---------
MyVar <- c("soaktime", "Xkm", "Ykm", "nettot", "temp_C", "turb_FNU", "DO_perc", "pH",
           "distance_to_south", "distance_to_shore", "sitedep_m")
Mydotplot(catchwb6_train[, MyVar])
# one catch above 150, fine
# few odd turbs, fine


## normality ----
ggplot(catchwb6_train, aes(sample=nettot)) +
  stat_qq() +
  stat_qq_line()
# not normally distirbuted

# zeros ----
sum(catchwb6$nettot == 0)  #Number of zeros
100 * sum(catchwb6$nettot == 0) / nrow(catchwb6)  #% of zeros
# 0

# collinearity density ------
MyVar <- c("soaktime", "nettot", "temp_C", "turb_FNU", "DO_perc", "pH",
           "distance_to_south", "sitedep_m", "distance_to_shore", "year")
Mypairs(catchwb6_train[, MyVar])

# high collinearity = over 0.3

MyVar <- c("soaktime", "nettot", "temp_C", "turb_FNU", "DO_perc", "pH",
           "distance_to_south", "distance_to_shore", "year")
Mypairs(catchwb6_train[, MyVar])
# no site dep

MyVar <- c("soaktime", "nettot", "temp_C", "turb_FNU", "DO_perc", "pH",
           "distance_to_south",  "year")
Mypairs(catchwb6_train[, MyVar])
# no distance to shore

MyVar <- c("soaktime", "nettot", "temp_C", "turb_FNU", "DO_perc",
           "distance_to_south",  "year")
Mypairs(catchwb6_train[, MyVar])
# no pH

boxplot(nettot ~ site, data = catchwb6)
boxplot(nettot ~ layer, data=catchwb6)
boxplot(nettot ~ year, data=catchwb6)

# turb
ggplot(data = catchwb6_train,  aes(x = turb_FNU, y = nettot))+
  geom_point() +
  geom_smooth() #(method = "lm")
# more or less increases with turb but few pts at higher turbs
# still stands when filter <= 25

# temp
ggplot(data = catchwb6_train,  aes(x = temp_C, y = nettot))+
  geom_point() +
  geom_smooth() 
# slight decrease with temp

# DO
ggplot(data = catchwb6_train,  aes(x = DO_perc, y = nettot))+
  geom_point() +
  geom_smooth()
# slight increase with DO

# pH
ggplot(data = catchwb6_train,  aes(x = pH, y = nettot))+
  geom_point() +
  geom_smooth() 
# slight increase with DO

# distance to any shore
ggplot(data = catchwb6_train,  aes(x = distance_to_shore, y = nettot))+
  geom_point() +
  geom_smooth() 
# slight increase with dist

# distance to south shore
ggplot(data = catchwb6_train,  aes(x = distance_to_south, y = nettot))+
  geom_point() +
  geom_smooth() 
# slight decrease with dist

# sitedep
ggplot(data = catchwb6_train,  aes(x = sitedep_m, y = nettot))+
  geom_point() +
  geom_smooth() 
# increase with site dep to about 30 m then delcine
# range only 3-55 m

# site effect big
ggplot(data = catchwb6_train,  aes(x = site, y = nettot))+
  geom_point() +
  geom_smooth() +#(method = "lm")
  geom_boxplot() 

ggplot(data=catchwb6_train, aes(x=lon, y=lat, colour=nettot)) +
  geom_point(alpha=0.7)

# stripey due to site grid

# year efffect 
ggplot(data = catchwb6_train,  aes(x = fyear, y = nettot))+
  geom_point() +
    geom_boxplot() 



# Relationships density -----
MyVar <- c("soaktime",	"turb_FNU",	"DO_perc", "pH", "sitedep_m",	"distance_to_south", "distance_to_shore",
           "temp_C", "year")
MyMultipanel.ggp2(Z = catchwb6_train, 
                  varx = MyVar, 
                  vary = "nettot", 
                  ylab = "Total num in net",
                  addSmoother = TRUE,
                  addRegressionLine = FALSE,
                  addHorizontalLine = FALSE)
# dist to south, site dep, soak, turb not lookin linear

# spatial
xyplot(Ykm ~ Xkm | factor(year),
       aspect = "iso",
       col = 1, 
       pch = 16,
       data = catchwb6_train)

# create workable set of var (VIF oriented)
MyVar <- c("soaktime", "sitedep_m", "turb_FNU", 
           "temp_C", "distance_to_shore", "distance_to_south", "year")
Mypairs(catchwb6[, MyVar])
corvif(catchwb6_train[,MyVar])
# remove sitedep at 3.25

MyVar <- c("soaktime",  "turb_FNU",
"temp_C", "distance_to_shore", "distance_to_south", "year")
Mypairs(catchwb6[, MyVar])
corvif(catchwb6[,MyVar])
# remove  distance to shore at 1.45

MyVar <- c("soaktime", "turb_FNU",  
           "temp_C", "distance_to_south", "year")
Mypairs(catchwb6[, MyVar])
corvif(catchwb6[,MyVar])



MyMultipanel.ggp2(Z = catchwb6_train, 
                  varx = MyVar, 
                  vary = "nettot", 
                  ylab = "nettot",
                  addSmoother = TRUE,
                  addRegressionLine = FALSE,
                  addHorizontalLine = FALSE)

# Relationships biomass -----
MyVar <- c("soaktime",	"turb_FNU",	"DO_perc", "pH", "sitedep_m",	"distance_to_south", "distance_to_shore",
           "temp_C", "year")
MyMultipanel.ggp2(Z = catchwb6_train, 
                  varx = MyVar, 
                  vary = "avgwt_g", 
                  ylab = "Avg Mass of LKWF",
                  addSmoother = TRUE,
                  addRegressionLine = FALSE,
                  addHorizontalLine = FALSE)
# dist to south, site dep, soak, turb not lookin linear

# spatial
xyplot(Ykm ~ Xkm | factor(year),
       aspect = "iso",
       col = 1, 
       pch = 16,
       data = catchwb6_train)

# create workable set of var (VIF oriented)
MyVar <- c("soaktime", "sitedep_m", "turb_FNU", 
           "temp_C", "distance_to_shore", "distance_to_south", "year")
Mypairs(catchwb6_train[, MyVar])
corvi(catchwb6_train[,MyVar])
# remove sitedep at 3.25

MyVar <- c("soaktime",  "turb_FNU",
           "temp_C", "distance_to_shore", "distance_to_south", "year")
Mypairs(catchwb6_train[, MyVar])
corvif(catchwb6_train[,MyVar])
# remove  distance to shore at 1.45

MyVar <- c("soaktime", "turb_FNU",  
           "temp_C", "distance_to_south", "year")
Mypairs(catchwb6_train[, MyVar])
corvif(catchwb6_train[,MyVar])



MyMultipanel.ggp2(Z = catchwb6_train, 
                  varx = MyVar, 
                  vary = "nettot", 
                  ylab = "nettot",
                  addSmoother = TRUE,
                  addRegressionLine = FALSE,
                  addHorizontalLine = FALSE)


# dependency -----
# since some site were sampled repeatedly opver time, site is a dependency

# independence of Y ------

## are smaller fish younger--
lkwfagewb <- lkwfage %>% 
  dplyr::filter(area == "IW" | area == "IE") %>% 
  droplevels()

lkwfagewb <- lkwfagewb %>% 
  dplyr::filter(gillnet == "benthic") %>% 
  droplevels()
lkwfagewb$area <- factor(lkwfagewb$area, levels=c("IW", "IE", "II", "III", "IV", "V"))

ggplot(mapping=aes(x=age, y=roundwt_g), data=lkwfagewb) +
  geom_point() +
  geom_smooth() +
  xlab("Age (years)") +
  ylab("Round weight (g)") +
  facet_wrap(~area)
