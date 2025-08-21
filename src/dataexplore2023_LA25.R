### Data exploration 
# Y = # LKWF individuals (netcount) -> count data so Poisson or NB
# OR
# Y = mean wt LKWF caught (avgwt_g) -> continous positive data so Gamma

# X = 
# 
# - soak time (NEED for effort)
# - year
# - turb
# - Dist from S shore
# - Dist from nearest shore
# - site dep 
# - temp


# WESTERN BASIN ONLY (FMAs IW and IE)
# ten panel bottom set gillnet only

#STEP 1. Load fxns and data and tidy things #########

## libraries ---------
library(ggplot2)
library(tidyverse)


## functions ----

# Mixed effects models and extensions in ecology with R. (2009).
# Zuur, AF, Ieno, EN, Walker, N, Saveliev, AA, and Smith, GM. Springer.
# https://www.highstat.com/index.php/books2?view=article&id=17&catid=18
source("./data/HighstatLibV13.R") # 
# used the corvif() and helper fxns and Mystd() fxn
corvif <- function(dataz) {
  dataz <- as.data.frame(dataz)
  
  #vif part
  form    <- formula(paste("fooy ~ ",paste(strsplit(names(dataz)," "),collapse=" + ")))
  dataz   <- data.frame(fooy=1 + rnorm(nrow(dataz)) ,dataz)
  lm_mod  <- lm(form,dataz)
  
  cat("\n\nVariance inflation factors\n\n")
  print(myvif(lm_mod))
}


#Support function for corvif. Will not be called by the user
myvif <- function(mod) {
  v <- vcov(mod)
  assign <- attributes(model.matrix(mod))$assign
  if (names(coefficients(mod)[1]) == "(Intercept)") {
    v <- v[-1, -1]
    assign <- assign[-1]
  } else warning("No intercept: vifs may not be sensible.")
  terms <- labels(terms(mod))
  n.terms <- length(terms)
  if (n.terms < 2) stop("The model contains fewer than 2 terms")
  if (length(assign) > dim(v)[1] ) {
    diag(tmp_cor)<-0
    if (any(tmp_cor==1.0)){
      return("Sample size is too small, 100% collinearity is present")
    } else {
      return("Sample size is too small")
    }
  }
  R <- cov2cor(v)
  detR <- det(R)
  result <- matrix(0, n.terms, 3)
  rownames(result) <- terms
  colnames(result) <- c("GVIF", "Df", "GVIF^(1/2Df)")
  for (term in 1:n.terms) {
    subs <- which(assign == term)
    result[term, 1] <- det(as.matrix(R[subs, subs])) * det(as.matrix(R[-subs, -subs])) / detR
    result[term, 2] <- length(subs)
  }
  if (all(result[, 2] == 1)) {
    result <- data.frame(GVIF=result[, 1])
  } else {
    result[, 3] <- result[, 1]^(1/(2 * result[, 2]))
  }
  invisible(result)
}

#Standardize the continuous covariates
MyStd <- function(x) { (x - mean(x)) / sd(x)}

# used the modification to pairs() in base R from
# Yuan Y, Cantoni E, Treble M, Flemming JM (2021).
# Spatiotemporal modeling of bycatch data: methods and a practical guide through 
# a case study in a Canadian Arctic fishery. Canadian Journal of Fisheries and 
# Aquatic Sciences 79(1) https://doi.org/10.1139/cjfas-2020-0267 
panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
}
panel.cor <- function(x, y, digits = 2, prefix = "", ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(round(r,2), nsmall = digits)
  txt <- paste0(prefix, txt)
  text(0.5, 0.5, txt, cex = 3*log(r+1.5))
}


mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}



## read data -----

# load("./data/datacleanoutput.RData") #clean data from dataclean2023.R
# run all of dataclean2024_LA25.R

## standardise -----
# Stadardize continuous variables for the model

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


mode(catchwb6$nettot)

# Does catch (totnum) differ with year and spatial location?
# Is there a relationship of catch and turbidity or thermal strata?


mean(catchwb6_train$soaktime) # 22.95
sd(catchwb6_train$soaktime) # 1.39


# STEP 2. Understand design ##########
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

# 
# citation('INLA')

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

# MyVar <- c("soaktime", "nettot", "temp_C", "turb_FNU", "DO_perc", "pH",
#            "distance_to_south", "sitedep_m", "distance_to_shore", "year")
# Mypairs(catchwb6_train[, MyVar])
# 
# # high collinearity = over 0.3
# 
# MyVar <- c("soaktime", "nettot", "temp_C", "turb_FNU", "DO_perc", "pH",
#            "distance_to_south", "distance_to_shore", "year")
# Mypairs(catchwb6_train[, MyVar])
# # no site dep
# 
# MyVar <- c("soaktime", "nettot", "temp_C", "turb_FNU", "DO_perc", "pH",
#            "distance_to_south",  "year")
# Mypairs(catchwb6_train[, MyVar])
# # no distance to shore
# 
# MyVar <- c("soaktime", "nettot", "temp_C", "turb_FNU", "DO_perc",
#            "distance_to_south",  "year")
# Mypairs(catchwb6_train[, MyVar])
# # no pH




# Relationships numbers -----
par(mfrow=c(2,2),mar=c(2,4,2,2))
boxplot(nettot ~ site, data = catchwb6, main="Site")
boxplot(nettot ~ layer, data=catchwb6, main="Strata")
boxplot(nettot ~ year, data=catchwb6, main="Year")

par(mfrow=c(1,2),mar=c(2,4,2,2))
boxplot(nettot ~ year, data=catchwb6, main="Numbers")
boxplot(avgwt_g ~ year, data=catchwb6, main="Mean weight (g)")

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


# # collinearity numbers ------
panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
}
panel.cor <- function(x, y, digits = 2, prefix = "", ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(round(r,2), nsmall = digits)
  txt <- paste0(prefix, txt)
  text(0.5, 0.5, txt, cex = 3*log(r+1.5))
}
# create workable set of var (VIF oriented)
MyVar <- c("soaktime", "sitedep_m", "turb_FNU", 
           "temp_C", "distance_to_shore", "distance_to_south", "year")
Mypairs(catchwb6[, MyVar])
cor(catchwb6[,MyVar], method="pearson")
pairs(catchwb6[, MyVar], diag.panel=panel.hist,pch=20, lower.panel = panel.cor)

corvif(catchwb6[,MyVar]) # corvif() from highstatlibv13.R
# remove sitedep at 10.53

MyVar <- c("soaktime",  "turb_FNU",
"temp_C", "distance_to_shore", "distance_to_south", "year")
Mypairs(catchwb6[, MyVar])
cor(catchwb6[,MyVar], method="pearson")
pairs(catchwb6[, MyVar], diag.panel=panel.hist,pch=20, lower.panel = panel.cor)

corvif(catchwb6[,MyVar])
# could remove  distance to shore at 2.33 but not needed since <3
# 
# MyVar <- c("soaktime", "turb_FNU",  
#            "temp_C", "distance_to_south", "year")
# Mypairs(catchwb6[, MyVar])
# corvif(catchwb6[,MyVar])



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
