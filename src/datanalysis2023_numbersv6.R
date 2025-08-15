library(lattice)  
library(ggplot2)
library(MASS)
library(mgcv)
library(plyr)
library(tidyr)
library(INLA)


library(sp)
library(gstat)
library(fields)
library(RColorBrewer)
source('C:/Users/alsipl/Documents/MSc/src/SPDEsupport.R')

inla.setOption(num.threads=8) # sets the number of cores being used by INLA

## covariates
# year: numeric, 8 years (2012-19)
# soak time: numeric, range 18-30 wiht most around 24
# net area: factor with 2 levels 

## control.compute
# DIC: baysian measure analagous to AIC. lower is best model fit with fewest covariates.
# DIC knoewn to be unreliable when posterior distribution is not well summarized by its mean
# WAIC: more fully bayseian
# MLIK: probability of the observed data under a given model.
# CPO: probability density of an observed response based on the model fit to the rest of the data. 
#if looking at cpo for each observation, larger values indicate better fit of model to the data. 
#-sum(logCPO) summarizes, with smaller values pointing to better model fit
# config: allows for prediction and simulation
catchwb6_pred <- catchwb6_pred %>% 
  mutate(set = ifelse(!is.na(nettot), "Training", "Testing"))

ggplot(data=catchwb6_pred) +
  geom_point(mapping=aes(x=soaktime, y=nettot), colour="black") +
  geom_point(mapping=aes(x=soaktime, y=resp_obsN), colour="red") +
  geom_smooth(mapping=aes(x=soaktime, y=resp_obsN), colour="red") +
  xlab ("Soak time (h)") +
  ylab ("Total LKWF per net") +
  xlim(18, 27)

# Base Density models ------
# Effort: soak, net area, year

# *Poisson GLM ------
# not used
I1 <- inla(nettot ~ year + soaktimec +
             f(site, model='iid'),
           control.compute = list(dic = TRUE, cpo=TRUE, waic=TRUE, config=TRUE),
           control.family = list(link='log'),
           family = "poisson",
           data = catchwb6_train)
summary(I1)
-sum(log(I1$cpo$cpo)) # smaller is better fit of model to data

# DIC 955
# WAIC 3158
# MLIK -673
# CPO 714

# Posterior mean values and 95% CI for the regression parameters:
Beta1<- I1$summary.fixed[,c("mean", "0.025quant", "0.975quant")] 
print(Beta1, digits = 2)  

#               mean 0.025quant 0.975quant
# (Intercept)  3.581      3.386     3.7717
# year        -0.011     -0.029     0.0067
# soaktimec   -0.055     -0.101    -0.0086
## suggests soak not important (0 is not in the quantile range)

# frequentist style
mu1  <- I1$summary.fitted.values[,"mean"]     #Fitted values
E1   <- (catchwb6_train$nettot - mu1) / sqrt(mu1)      #Pearson residuals
N    <- nrow(catchwb6_train)                             #Sample size
Npar <- nrow(I1$summary.fixed)                #Number of parameters
Dispersion <- sum(E1^2) / (N - Npar)          #Dispersion statistic (should be 1)
Dispersion
# 10 - overdispersed! 
# try NB


range(catchwb6_train$nettot) # 5-130

# why overdispersed
# fitted values versus the Pearson residuals
par(mfrow=c(2,1))
plot(x = mu1,
     y = E1,
     xlab = "Fitted values",
     ylab = "Pearson residuals",
     cex.lab = 1.5)
abline(h = 0, lty = 2)
# fitted ranges 20 to 70 (whole range should be 5 to 130)
# residuals evenly spread
# wide range of residuals (-5, +10)

boxplot(E1 ~ site, 
        ylab = "Pearson residuals",
        data = catchwb6_train)
abline(h = 0, lty = 2)

#residuals vs each covariate (in and not in the model).
# Check for non-linear patterns in the residuals.
MyVar <- c("turb_FNU", "sitedep_m", "temp_C",  "soaktime",  
           "distance_to_south", "distance_to_shore", "year")
catchwb6_train$E1 <- E1
MyMultipanel.ggp2(Z = catchwb6_train, 
                  varx = MyVar, 
                  vary = "E1", 
                  ylab = "Residuals",
                  addSmoother = TRUE,
                  addRegressionLine = FALSE,
                  addHorizontalLine = TRUE)
# non linearity present in aturb and year 

# spatial dependency?
catchwb6_train$Sign <- ifelse(E1 >= 0, "Positive residual", "Negative residual") 
ggplot() +
  theme(text = element_text(size=13)) +
  geom_point(aes(x = Xkm, y = Ykm, col = Sign), shape = 1, size = 4, data = catchwb6_train) +
  xlab("Xkm") +
  ylab("Ykm")  
# clumpy - suggest spatial dependency

MyData2 <- data.frame(E1 = E1, 
                      Xkm = catchwb6_train$Xkm, 
                      Ykm = catchwb6_train$Ykm)
coordinates(MyData2) <- c("Xkm", "Ykm")
V1 <- variogram(E1 ~ 1, 
                MyData2, 
                cressie = TRUE, 
                cutoff = 70)
ggplot() + 
  geom_point(data = V1, aes(x = dist, y= gamma)) +
  geom_smooth(data = V1, se = FALSE, span = 0.9, aes(x = dist, y = gamma)) +
  xlab("Distance") + ylab("Semi-variogram") +
  theme(text = element_text(size = 15)) +
  theme(legend.position="none") 
# YES - up to about 35 km

# par(mfrow=c(1,1))
# acf(catchwb6_train$nettot, catchwb6_train$year)

### obs vs pred 1a ---------

I1a <- inla(nettot ~ year + soaktimec +
              f(site, model='iid'),
           control.compute = list(dic = TRUE, cpo=TRUE, waic=TRUE, config=TRUE),
           family = "poisson",
           control.family=list(link='log'),
           data = catchwb6_pred) # now pred
summary(I1a) # should be the same as before

#extract the predictor values associated with the appended data
catchwb6_test1a <- cbind(catchwb6_test, 
                        I1a$summary.linear.predictor[(nrow(catchwb6_train)+1):nrow(catchwb6_pred),])
head(catchwb6_test1a)

catchwb6_test1a <- reshape:::rename(catchwb6_test1a, c("0.025quant"="lower", "0.975quant"="upper"))

ggplot(catchwb6_test1a, aes(y=mean, x=soaktimec)) +
  geom_point(data=catchwb6_train, aes(y=nettot)) + # obs
  geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.2) + # predicted SD
  geom_line() # predicted
# MUCH smaller than observed

ggplot(catchwb6_test1a, aes(y=mean, x=soaktimec)) +
  geom_point(aes(y=resp_obsN)) + # obs
  geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.2) + # predicted SD
  geom_line()

### sim 1 overdispersion------ 

s1 <- inla(nettot ~ year + soaktimec +
             f(site, model='iid'),
           control.compute=list(config = TRUE),
           family = "poisson",
           data = catchwb6_train)

set.seed(12345)
SimData <- inla.posterior.sample(n = 1, result = s1)
BetaRows <- 82:85  
Betas    <- SimData[[1]]$latent[BetaRows]

X <- model.matrix(~year + soaktimec, data = catchwb6_train)
mu <- exp(X %*% Betas)
Ysim <- rpois(n = nrow(catchwb6_train), lambda = mu)
Es <- (Ysim - mu) /sqrt(mu)
sum(Es^2) / (N - Npar)
# truly Poisson data (simulated from the model). It should be close to 1.
## 1.33

# sum of squared Pearson residuals (divided by N-p) for the original data and model is:
E1 <- (catchwb6_train$nettot - mu1) /sqrt(mu1)
Dispersion <- sum(E1^2) / (N - Npar)
Dispersion    
## 24.17

NSim <- 1000
SimData <- inla.posterior.sample(n = NSim, result = s1)
N      <- nrow(catchwb6_train)
Ysim   <- matrix(nrow = N, ncol = NSim)
mu.sim <- matrix(nrow = N, ncol = NSim)
for (i in 1:NSim){
  Betas      <- SimData[[i]]$latent[BetaRows]
  mu.sim[,i] <- exp(X %*% Betas)
  Ysim[,i]   <- rpois(n = N, lambda = exp(X %*% Betas))
}


E1 <- (catchwb6_train$nettot - mu1) /sqrt(mu1)
Dispersion <- sum(E1^2) / (N - Npar)
Dispersion

# We calculate the same statistic for 
# each of the simulated data sets.
Dispersion.sim <- NULL
for(i in 1:NSim){
  e2 <- (Ysim[,i] - mu.sim[,i]) /sqrt(mu.sim[,i])
  Dispersion.sim[i] <- sum(e2^2) / (N - Npar)
}


# Let us plot the results
par(mfrow = c(1,1), mar = c(5,5,2,2), cex.lab = 1.5)
hist(Dispersion.sim,
     xlab = "Dispersion",
     ylab = "Frequency",
     xlim = c(0.6, 25),
     breaks = 25)

# And visualize the dispersion for the original Poisson GLMM
points(x = Dispersion, 
       y = 0, 
       pch = 16, 
       cex = 3, 
       col = 2)
# super over dispersed

# *NB GLM ------
# not used
I2 <- inla(nettot ~ year + soaktimec +
             f(site, model='iid'),
           control.compute = list(dic = TRUE, cpo=TRUE, waic=TRUE, config=TRUE),
           family = "nbinomial",
           control.family = list(link='log'),
           data = catchwb6_train)
summary(I2)
-sum(log(I2$cpo$cpo)) # smaller is better fit of model to data

# DIC 685
# WAIC 686
# MLIK -358
# CPO 343

# Posterior mean values and 95% CI for the regression parameters:
Beta2<- I2$summary.fixed[,c("mean", "0.025quant", "0.975quant")] 
print(Beta2, digits = 2)  


#                 mean 0.025quant 0.975quant
# (Intercept)  3.781      3.425       4.14
# year        -0.035     -0.099       0.03
# soaktimec   -0.144     -0.329       0.04
# suggest all var impt

# posterior mean of theta
theta.pd <- I2$marginals.hyperpar$`size for the nbinomial observations (1/overdispersion)`
theta.pm <- inla.emarginal(function(x) x, theta.pd)
theta.pm ## 2.5
# posterior distribution of theta
par(mfrow = c(1,1), mar = c(5,5,2,2), cex.lab = 1.5)
plot(theta.pd, 
     type = "l",
     xlim = c(1,5),
     xlab = expression(paste(theta)),
     ylab = expression(paste("Pr(", theta," | data)")))

# overdispersed?
mu2 <- I2$summary.fitted.values[,"mean"]
E2  <- (catchwb6_train$nettot - mu2) /sqrt(mu2 + mu2^2 / theta.pm)
N <- nrow(catchwb6_train)
p <- nrow(I2$summary.fixed) 

Dispersion.nb <- sum(E2^2) / (N - p)
Dispersion.nb  
## 0.97
## will need to check to see if this is okay with simulatin

# why underdispersed
# fitted values versus the Pearson residuals
par(mfrow=c(2,1))
plot(x = mu2,
     y = E2,
     xlab = "Fitted values",
     ylab = "Pearson residuals",
     cex.lab = 1.5)
abline(h = 0, lty = 2)
# fitted ranges 30-60 (whole range should be 1 to 146)
# residuals mostly on lower hald
# wide range of residuals (-1 to +4)

boxplot(E2 ~ site, 
        ylab = "Pearson residuals",
        data = catchwb6_train)
abline(h = 0, lty = 2)

#residuals vs each covariate (in and not in the model).
# Check for non-linear patterns in the residuals.
MyVar <- c("turb_FNU", "sitedep_m", "temp_C", "DO_perc", "pH", "soaktime",  
           "distance_to_south", "distance_to_shore", "year")
catchwb6_train$E2 <- E2
MyMultipanel.ggp2(Z = catchwb6_train, 
                  varx = MyVar, 
                  vary = "E2", 
                  ylab = "Residuals",
                  addSmoother = TRUE,
                  addRegressionLine = FALSE,
                  addHorizontalLine = TRUE)
# non linearity present in all covariates, aside from temp d2S and soak

# spatial dependency?
catchwb6_train$Sign <- ifelse(E2 >= 0, "Positive residual", "Negative residual") 
ggplot() +
  theme(text = element_text(size=13)) +
  geom_point(aes(x = Xkm, y = Ykm, col = Sign), shape = 1, size = 4, data = catchwb6_train) +
  xlab("Xkm") +
  ylab("Ykm")  
# clumpy - suggest spatial dependency

MyData2 <- data.frame(E2 = E2, 
                      Xkm = catchwb6_train$Xkm, 
                      Ykm = catchwb6_train$Ykm)
coordinates(MyData2) <- c("Xkm", "Ykm")
V2 <- variogram(E2 ~ 1, 
                MyData2, 
                cressie = TRUE, 
                cutoff = 70)
ggplot() + 
  geom_point(data = V2, aes(x = dist, y= gamma)) +
  geom_smooth(data = V2, se = FALSE, span = 0.9, aes(x = dist, y = gamma)) +
  xlab("Distance") + ylab("Semi-variogram") +
  theme(text = element_text(size = 15)) +
  theme(legend.position="none") 
# YES - up to about 35 km



### obs vs pred 2a ---------
I2a <- inla(nettot ~ year + soaktimec +
              f(site, model='iid'),
           control.compute = list(dic = TRUE, waic=TRUE, config=TRUE),
           family = "nbinomial",
           control.family = list(link='log'),
           data = catchwb6_pred) # now pred
summary(I2a) # should be the same as before

#extract the predictor values associated with the appended data
catchwb6_test2a <- cbind(catchwb6_test, 
                        I2a$summary.linear.predictor[(nrow(catchwb6_train)+1):nrow(catchwb6_pred),])
head(catchwb6_test2a)

catchwb6_test2a <- reshape:::rename(catchwb6_test2a, c("0.025quant"="lower", "0.975quant"="upper"))

ggplot(catchwb6_test2a, aes(y=mean, x=soaktimec)) +
  geom_point(data=catchwb6_train, aes(y=nettot)) + # obs
  geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.2) + # predicted SD
  geom_line() # predicted
# MUCH smaller than observed

ggplot(catchwb6_test2a, aes(y=mean, x=soaktimec)) +
  geom_point(aes(y=resp_obsN)) + # obs
  geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.2) + # predicted SD
  geom_line()

p <- ggplot()
p <- p + geom_point(data = catchwb6_train, 
                    aes(y = nettot, x = soaktime),
                    size = 1, colour="black")
p <- p + xlab("Soak time (h)") + ylab("Total LKWF per net")
p <- p+ xlim(18,27)
p <- p + theme(text = element_text(size=15)) 
p <- p + geom_point(data = catchwb6_test2a, 
                    aes(x = soaktime, 
                        y = mean), 
                    colour = "red")
# p <- p + geom_ribbon(data = catchwb6_test2, 
#                      aes(x = soaktime, 
#                          ymax = upper, 
#                          ymin = lower),
#                      alpha = 0.2) 
p <- p + geom_smooth(data = catchwb6_test2a, 
                    aes(x = soaktime, 
                        y = mean), 
                    colour = "red")
p

ggplot(data=catchwb6_test2) +
  geom_point(mapping=aes(x=resp_obsN, y=mu, colour=as.factor(catchwb6_test2$year)), size=3) +
  xlab("Observed") +
  ylab("Predicted") +
  labs(colour="Year") +
  geom_abline() +
  xlim(c(0,120)) +
  ylim(c(0,120)) 

### sim 2 underdispersion------ 
s2 <- inla(nettot ~ soaktimec + year + 
             f(site, model='iid'),
                 control.compute = list(config = TRUE),  ##compute TRUE
                 family = "nbinomial", 
                 data = catchwb6_train)

NSim <- 1000
SimData <- inla.posterior.sample(n = NSim, result = s2)
nd <- length(rownames(SimData[[1]]$latent))
LastBeta <- rownames(SimData[[1]]$latent)[nd]
substrRight <- function(x, n){ substr(x, nchar(x)-n+1, nchar(x)) }
Last2Character <- substrRight(LastBeta, 2)
BetasInModel <- rownames(s2$summary.fixed) #<--Change the name 'I1.sim' for your model
MyID         <- function(x){ which(rownames(SimData[[1]]$latent) == x) }
BetasInINLA <- BetasInModel
if (Last2Character == ":1") { BetasInINLA  <- paste(BetasInModel, ":1", sep ="") }
BetaRows <- lapply(BetasInINLA, MyID)
BetaRows <- as.numeric(BetaRows)
BetaRows

X    <- model.matrix(~ soaktimec  + year , data = catchwb6_train)
Xmat <- as.matrix(X)

N    <- nrow(catchwb6_train)                     #Sample size
Ysim <- matrix(nrow = N, ncol = NSim) #Create space for simulated data
mu.i <- matrix(nrow = N, ncol = NSim) #Create space for expected values
theta <- vector(length = NSim)        #Create space for simulated thetas
for (i in 1:NSim){
  Betas    <- SimData[[i]]$latent[BetaRows]       #Get betas  
  eta      <- Xmat %*% Betas #X * beta + a[VesselID]
  theta[i] <- SimData[[i]]$hyperpar["size for the nbinomial observations (1/overdispersion)"]
  mu.i[,i] <- exp(eta)                            #Fitted values
  Ysim[,i] <- rnegbin(n = N,                      #Simulated NB BallanWrasse data
                      mu = mu.i[,i] ,
                      theta = theta[i]) 
}

# Calculate the dispersion statistic for each of these 1,000
# data sets.
DispNB <- vector(length = NSim)             #Create space
N      <- nrow(catchwb6_train)                         #Sample size
Npar  <- nrow(I2$summary.fixed) +1 #Number of parameters
for(i in 1:NSim){
  ei <- (Ysim[,i] - mu.i[,i]) / sqrt(mu.i[,i] + mu.i[,i]^2 / theta[i])
  DispNB[i] <- sum(ei^2) / (N - Npar)
}


# Let's plot this as a table
par(mfrow = c(1,1), mar = c(5,5,2,2), cex.lab = 1.5)
hist(DispNB, 
     xlim = c(0, 3),
     xlab = "Dispersion statistic",
     ylab = "Frequency",
     main = "Simulation results")
points(x = Dispersion.nb, 
       y = 0, 
       pch = 16, 
       cex = 5, 
       col = 2)
# so the dispersion is fine

# NB GAM ------
I3 <- inla(nettot ~ soaktimec +
             f(site, model='iid') +
             f(year, model='rw1'),
           control.compute = list(dic = TRUE, cpo=TRUE, waic=TRUE, config=TRUE),
           family = "nbinomial",
           control.family = list(link='log'),
           data = catchwb6_train)
summary(I3)
-sum(log(I3$cpo$cpo)) # smaller is better fit of model to data

HyperPar.I3 <- bri.hyperpar.summary(I3)
round(HyperPar.I3, digits = 3)
# year and site impt

# Posterior mean values and 95% CI for the regression parameters:
Beta3<- I3$summary.fixed[,c("mean", "0.025quant", "0.975quant")] 
print(Beta3, digits = 2)  

#                   mean 0.025quant 0.975quant
# (Intercept)  3.605       3.47      3.743
# soaktimec   -0.086      -0.25      0.079

# posterior mean of theta
theta.pd <- I3$marginals.hyperpar$`size for the nbinomial observations (1/overdispersion)`
theta.pm <- inla.emarginal(function(x) x, theta.pd)
theta.pm ## 3.35
# posterior distribution of theta
par(mfrow = c(1,1), mar = c(5,5,2,2), cex.lab = 1.5)
plot(theta.pd, 
     type = "l",
     xlim = c(1,6),
     xlab = expression(paste(theta)),
     ylab = expression(paste("Pr(", theta," | data)")))

# overdispersed?
mu3 <- I3$summary.fitted.values[,"mean"]
E3  <- (catchwb6_train$nettot - mu3) /sqrt(mu3 + mu3^2 / theta.pm)
N <- nrow(catchwb6_train)
p <- nrow(I3$summary.fixed) 
Dispersion.nb <- sum(E3^2) / (N - p)
Dispersion.nb  
## 0.96

# why underdispersed
# fitted values versus the Pearson residuals
par(mfrow=c(2,1))
plot(x = mu3,
     y = E3,
     xlab = "Fitted values",
     ylab = "Pearson residuals",
     cex.lab = 1.5)
abline(h = 0, lty = 2)

boxplot(E3 ~ site, 
        ylab = "Pearson residuals",
        data = catchwb6_train)
abline(h = 0, lty = 2)

#residuals vs each covariate (in and not in the model).
# Check for non-linear patterns in the residuals.
MyVar <- c("turb_FNU", "sitedep_m", "temp_C", "DO_perc", "pH", "soaktime",  
           "distance_to_south", "distance_to_shore", "year")
catchwb6_train$E3 <- E3
MyMultipanel.ggp2(Z = catchwb6_train, 
                  varx = MyVar, 
                  vary = "E3", 
                  ylab = "Residuals",
                  addSmoother = TRUE,
                  addRegressionLine = FALSE,
                  addHorizontalLine = TRUE)
# non linearity present in all covariates, aside from temp and soak
# smooth on soaktime didn't work?

# spatial dependency?
catchwb6_train$Sign <- ifelse(E3 >= 0, "Positive residual", "Negative residual") 
ggplot() +
  theme(text = element_text(size=13)) +
  geom_point(aes(x = Xkm, y = Ykm, col = Sign), shape = 1, size = 4, data = catchwb6_train) +
  xlab("Xkm") +
  ylab("Ykm")  
# clumpy - suggest spatial dependency

MyData3 <- data.frame(E3 = E3, 
                      Xkm = catchwb6_train$Xkm, 
                      Ykm = catchwb6_train$Ykm)
coordinates(MyData3) <- c("Xkm", "Ykm")
V3 <- variogram(E3 ~ 1, 
                MyData3, 
                cressie = TRUE, 
                cutoff = 50)
ggplot() + 
  geom_point(data = V3, aes(x = dist, y= gamma)) +
  geom_smooth(data = V3, se = FALSE, span = 0.9, aes(x = dist, y = gamma)) +
  xlab("Distance") + ylab("Semi-variogram") +
  theme(text = element_text(size = 15)) +
  theme(legend.position="none") 
# not really


### obs vs pred 3a ---------

I3a <- inla(nettot ~ soaktimec +
              f(site, model='iid') +
              f(year, model='rw1'),
            control.compute = list(dic = TRUE, cpo=TRUE, waic=TRUE, config=TRUE),
            family = "nbinomial",
            control.family = list(link='log'),
            data = catchwb6_pred) # now pred
summary(I3a) # should be the same as before

#extract the predictor values associated with the appended data
catchwb6_test3a <- cbind(catchwb6_test, 
                         I3a$summary.linear.predictor[(nrow(catchwb6_train)+1):nrow(catchwb6_pred),])
head(catchwb6_test3a)

catchwb6_test3a <- reshape:::rename(catchwb6_test3a, c("0.025quant"="lower", "0.975quant"="upper"))



# And the rest is a matter of plotting it all. 
# First we back-standardise LogMining. It is also an option to back-transform it.
catchwb6_test3a$soaktime2 <- catchwb6_test3a$soaktimec * sd(catchwb6_train$soaktime) + mean(catchwb6_train$soaktime)

catchwb6_pred <- catchwb6_train %>% 
  full_join(catchwb6_test3a)
catchwb6_pred <- catchwb6_pred %>% 
  mutate(set = ifelse(!is.na(nettot), "Training", "Testing"))
catchwb6_pred$set <- as.factor(catchwb6_pred$set)


catchwb6_pred <- catchwb6_pred %>% 
  mutate(nettot=ifelse(!is.na(nettot), nettot, mean))

head(catchwb6_pred)

p <- ggplot()
p <- p + geom_smooth(data = catchwb6_test3a, 
                     aes(x = soaktime, 
                         y = mean), alpha=0.5)
p <- p + geom_point(data = catchwb6_pred, 
                    aes(y = nettot, x = soaktime, col=set),
                    size = 3)
p <- p + xlab("Soak time (h)") + ylab("Total LKWF") + labs ( col = "Dataset")
p <- p + theme(text = element_text(size=15), legend.position = "bottom") 
# p <- p + geom_ribbon(data = catchwb6_test3a,
#                      aes(x = soaktime2,
#                          ymax = upper,
#                          ymin = lower),
#                      alpha = 0.2)

p

p <- ggplot()
p <- p + geom_smooth(data = catchwb6_pred %>% dplyr::filter(set=="Training"), 
                     aes(x = soaktime, 
                         y = nettot), alpha=0.5)
p <- p + geom_point(data = catchwb6_pred, 
                    aes(y = nettot, x = soaktime, col=set),
                    size = 3)
p <- p + xlab("Soak time (h)") + ylab("Total LKWF") + labs ( col = "Dataset")
p <- p + theme(text = element_text(size=15), legend.position = "bottom") 
p <- p+ xlim(18,27)
# p <- p + geom_ribbon(data = catchwb6_test3a,
#                      aes(x = soaktime2,
#                          ymax = upper,
#                          ymin = lower),
#                      alpha = 0.2)

p

ggplot(data=catchwb6_test3a) +
  geom_point(mapping=aes(x=resp_obsN, y=mean), size=3) +
  xlab("Observed") +
  ylab("Predicted") +
  ylim(0,130) + xlim(0,130) +
  theme(text=element_text(size=15)) +
  geom_abline()

rsq <- function (x, y) cor(x, y) ^ 2
rsq(catchwb6_test3a$resp_obsN, catchwb6_test3a$mean)

summary(lm(resp_obsN ~ mean, data=catchwb6_test3a))$r.squared

# mean and SD of predicted  vs obs in train
ggplot(catchwb6_test3a, aes(y=mean, x=soaktime)) +
  geom_point(data=catchwb6_train, aes(y=nettot)) + # obs
  geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.2) + # predicted SD of pred mean num
  geom_point(colour="blue") # predicted
# MUCH smaller than observed

# mean and SD of predicted vs obs in test
ggplot(catchwb6_test3a, aes(y=mean, x=soaktimec)) +
  geom_point(aes(y=resp_obsN)) + # obs
  geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.2) + # predicted SD
  geom_point(colour="blue")

## sim 3 underdispersion -----

s3 <- inla(nettot ~ soaktimec + 
             f(year, model='rw1')+
             f(site, model='iid'),
           control.compute = list(config = TRUE),  ##compute TRUE
           family = "nbinomial", 
           data = catchwb6_train)

set.seed(12345)

NSim <- 1000
SimData <- inla.posterior.sample(n = NSim, result = s3) #CHANGE
nd <- length(rownames(SimData[[1]]$latent))
LastBeta <- rownames(SimData[[1]]$latent)[nd]
substrRight <- function(x, n){ substr(x, nchar(x)-n+1, nchar(x)) }
Last2Character <- substrRight(LastBeta, 2)
BetasInModel <- rownames(s3$summary.fixed) #CHANGE
MyID         <- function(x){ which(rownames(SimData[[1]]$latent) == x) }
BetasInINLA <- BetasInModel
if (Last2Character == ":1") { BetasInINLA  <- paste(BetasInModel, ":1", sep ="") }
BetaRows <- lapply(BetasInINLA, MyID)
BetaRows <- as.numeric(BetaRows)
BetaRows

# Now the random effects.
MyGrep <- function(x, SimulatedData){ 
  # Small function to determine the row numbers of random effects 
  # in the output of inla.posterior.sample
  # This is partly code from Aaron Adamack (Aaron.Adamack@dfo-mpo.gc.ca)
  # SimulatedData is the name of the object containing the simulation results
  # x is an element of BetasInINLA
  names.SimData <- attributes(SimulatedData[[1]]$latent)$dimnames[[1]]
  names.SimData[grep(x, names.SimData)] 
}

# Determine the names of the random effects
names.SimData <- attributes(SimData[[1]]$latent)$dimnames[[1]] 
YearNames     <- names.SimData[grep("year", names.SimData)]   #<-Change 'Site' for your model
YearNames

YearRows <- lapply(YearNames, MyID)
YearRows <- as.numeric(YearRows)
YearRows


SiteNames     <-  names.SimData[grep("site", names.SimData)]
SiteNames

SiteRows <- lapply(SiteNames, MyID)
SiteRows <- as.numeric(SiteRows)
SiteRows

# Start a loop to extract betas and random  effects, calculate
# the fitted values and simulate count data from 
# the model.
N    <- nrow(catchwb6_train)
Ysim <- matrix(nrow = N, ncol = NSim)
mu.i <- matrix(nrow = N, ncol = NSim)
X    <- model.matrix(~ soaktimec, data = catchwb6_train)
X   <- as.matrix(X)

YearID <- as.numeric(as.factor(catchwb6_train$year))
SiteID <- as.numeric(as.factor(catchwb6_train$site))

theta <- vector(length = NSim)        #Create space for simulated thetas



for (i in 1: NSim){
  Betas    <- SimData[[i]]$latent[BetaRows]       #Get simulated betas
  Year_i <- SimData[[i]]$latent[YearRows] #Get simulated random intercepts for year
  Site_ij  <- SimData[[i]]$latent[SiteRows]   #Get simulated random intercepts for site
  
  #eta = X * beta + year_i + Site_ij +
  eta <- X %*% Betas + Year_i[YearID] + Site_ij[SiteID] 
  
  #Get the simulated theta  
  theta[i] <- SimData[[i]]$hyperpar["size for the nbinomial observations (1/overdispersion)"]
  
  # Calculate mu = exp(eta)
  mu.i[,i] <- exp(eta)                        #Fitted values
  Ysim[,i] <- rnegbin(n = N,                  #Simulated count data
                      mu = mu.i[,i] ,
                      theta = theta[i]) 
}


# Now we have 1000 simulated data sets from the model.
# What shall we do with these simulated data sets?
# Calculate the dispersion for each simulated data set
DispNB <- vector(length = NSim) #Create space
N      <- nrow(catchwb6_train)              #Sample size 
Npar   <- length(Betas) + 2+ 1     #Number of regression parameters + theta
for(i in 1:NSim){
  ei <- (Ysim[,i] - mu.i[,i]) / sqrt(mu.i[,i] + mu.i[,i]^2 / theta[i]) #NB Pearson residuals
  DispNB[i] <- sum(ei^2) / (N - Npar)                                   #Dispersion statistic
}


# Let's plot this as a table
par(mfrow = c(1,1), mar = c(5,5,2,2), cex.lab = 1.5)
hist(DispNB, 
     xlim = c(0, 3),
     xlab = "Dispersion statistic",
     ylab = "Frequency",
     main = "Simulation results")
points(x = Dispersion.nb, 
       y = 0, 
       pch = 16, 
       cex = 5, 
       col = 2)


# *NB GAM no site------
I3 <- inla(nettot ~ soaktimec +
             f(year, model='rw1'),
           control.compute = list(dic = TRUE, cpo=TRUE, waic=TRUE, config=TRUE),
           family = "nbinomial",
           control.family = list(link='log'),
           data = catchwb6_train)
summary(I3)
-sum(log(I3$cpo$cpo)) # smaller is better fit of model to data

HyperPar.I3 <- bri.hyperpar.summary(I3)
round(HyperPar.I3, digits = 3)


# Posterior mean values and 95% CI for the regression parameters:
Beta3<- I3$summary.fixed[,c("mean", "0.025quant", "0.975quant")] 
print(Beta3, digits = 2)  

#                   mean 0.025quant 0.975quant
# (Intercept)  3.605       3.47      3.743
# soaktimec   -0.086      -0.25      0.079

# posterior mean of theta
theta.pd <- I3$marginals.hyperpar$`size for the nbinomial observations (1/overdispersion)`
theta.pm <- inla.emarginal(function(x) x, theta.pd)
theta.pm ## 3.35
# posterior distribution of theta
par(mfrow = c(1,1), mar = c(5,5,2,2), cex.lab = 1.5)
plot(theta.pd, 
     type = "l",
     xlim = c(1,6),
     xlab = expression(paste(theta)),
     ylab = expression(paste("Pr(", theta," | data)")))

# overdispersed?
mu3 <- I3$summary.fitted.values[,"mean"]
E3  <- (catchwb6_train$nettot - mu3) /sqrt(mu3 + mu3^2 / theta.pm)
N <- nrow(catchwb6_train)
p <- nrow(I3$summary.fixed) 
Dispersion.nb <- sum(E3^2) / (N - p)
Dispersion.nb  
## 0.96

# why underdispersed
# fitted values versus the Pearson residuals
par(mfrow=c(2,1))
plot(x = mu3,
     y = E3,
     xlab = "Fitted values",
     ylab = "Pearson residuals",
     cex.lab = 1.5)
abline(h = 0, lty = 2)

boxplot(E3 ~ site, 
        ylab = "Pearson residuals",
        data = catchwb6_train)
abline(h = 0, lty = 2)

#residuals vs each covariate (in and not in the model).
# Check for non-linear patterns in the residuals.
MyVar <- c("turb_FNU", "sitedep_m", "temp_C", "DO_perc", "pH", "soaktime",  
           "distance_to_south", "distance_to_shore", "year")
catchwb6_train$E3 <- E3
MyMultipanel.ggp2(Z = catchwb6_train, 
                  varx = MyVar, 
                  vary = "E3", 
                  ylab = "Residuals",
                  addSmoother = TRUE,
                  addRegressionLine = FALSE,
                  addHorizontalLine = TRUE)
# non linearity present in all covariates, aside from temp and soak
# smooth on soaktime didn't work?

# spatial dependency?
catchwb6_train$Sign <- ifelse(E3 >= 0, "Positive residual", "Negative residual") 
ggplot() +
  theme(text = element_text(size=13)) +
  geom_point(aes(x = Xkm, y = Ykm, col = Sign), shape = 1, size = 4, data = catchwb6_train) +
  xlab("Xkm") +
  ylab("Ykm")  
# clumpy - suggest spatial dependency

MyData3 <- data.frame(E3 = E3, 
                      Xkm = catchwb6_train$Xkm, 
                      Ykm = catchwb6_train$Ykm)
coordinates(MyData3) <- c("Xkm", "Ykm")
V3 <- variogram(E3 ~ 1, 
                MyData3, 
                cressie = TRUE, 
                cutoff = 70)
ggplot() + 
  geom_point(data = V3, aes(x = dist, y= gamma)) +
  geom_smooth(data = V3, se = FALSE, span = 0.9, aes(x = dist, y = gamma)) +
  xlab("Distance") + ylab("Semi-variogram") +
  theme(text = element_text(size = 15)) +
  theme(legend.position="none") 
# not really


### obs vs pred 3a ---------

I3a <- inla(nettot ~ soaktimec +
              f(year, model='rw1'),
           control.compute = list(dic = TRUE, cpo=TRUE, waic=TRUE, config=TRUE),
           family = "nbinomial",
           control.family = list(link='log'),
           data = catchwb6_pred) # now pred
summary(I3a) # should be the same as before

#extract the predictor values associated with the appended data
catchwb6_test3a <- cbind(catchwb6_test, 
                        I3a$summary.linear.predictor[(nrow(catchwb6_train)+1):nrow(catchwb6_pred),])
head(catchwb6_test3a)

catchwb6_test3a <- reshape:::rename(catchwb6_test3a, c("0.025quant"="lower", "0.975quant"="upper"))

ggplot(catchwb6_test3a, aes(y=mean, x=soaktimec)) +
  geom_point(data=catchwb6_train, aes(y=nettot)) + # obs
  geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.2) + # predicted SD
  geom_line() # predicted
# MUCH smaller than observed

ggplot(catchwb6_test3a, aes(y=mean, x=soaktimec)) +
  geom_point(aes(y=resp_obsN)) + # obs
  geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.2) + # predicted SD
  geom_line()
catchwb6_test3a$soaktime2 <- catchwb6_test3a$soaktimec * sd(catchwb6_train$soaktime) + mean(catchwb6_train$soaktime)
p <- ggplot()
p <- p + geom_point(data = catchwb6_train, 
                    aes(y = nettot, x = soaktime),
                    size = 1)
p <- p + xlab("soak time (h)") + ylab("total number LKWF")
p <- p + theme(text = element_text(size=15)) 
p <- p + geom_point(data = catchwb6_test3a, 
                    aes(x = soaktime2, 
                        y = mean), 
                    col = "red")
p <- p + geom_ribbon(data = catchwb6_test3a, 
                     aes(x = soaktime2, 
                         ymax = upper, 
                         ymin = lower),
                     alpha = 0.2)
p

ggplot(data=catchwb6_test3a) +
  geom_point(mapping=aes(x=resp_obsN, y=mean)) +
  xlab("Observed") +
  ylab("Predicted") 
## sim 3 underdispersion -----

s3 <- inla(nettot ~ soaktimec + 
             f(year, model='rw1'),
           control.compute = list(config = TRUE),  ##compute TRUE
           family = "nbinomial", 
           data = catchwb6_train)

set.seed(12345)

NSim <- 1000
SimData <- inla.posterior.sample(n = NSim, result = s3) #CHANGE
nd <- length(rownames(SimData[[1]]$latent))
LastBeta <- rownames(SimData[[1]]$latent)[nd]
substrRight <- function(x, n){ substr(x, nchar(x)-n+1, nchar(x)) }
Last2Character <- substrRight(LastBeta, 2)
BetasInModel <- rownames(s3$summary.fixed) #CHANGE
MyID         <- function(x){ which(rownames(SimData[[1]]$latent) == x) }
BetasInINLA <- BetasInModel
if (Last2Character == ":1") { BetasInINLA  <- paste(BetasInModel, ":1", sep ="") }
BetaRows <- lapply(BetasInINLA, MyID)
BetaRows <- as.numeric(BetaRows)
BetaRows

# Now the random effects.
MyGrep <- function(x, SimulatedData){ 
  # Small function to determine the row numbers of random effects 
  # in the output of inla.posterior.sample
  # This is partly code from Aaron Adamack (Aaron.Adamack@dfo-mpo.gc.ca)
  # SimulatedData is the name of the object containing the simulation results
  # x is an element of BetasInINLA
  names.SimData <- attributes(SimulatedData[[1]]$latent)$dimnames[[1]]
  names.SimData[grep(x, names.SimData)] 
}

# Determine the names of the random effects
names.SimData <- attributes(SimData[[1]]$latent)$dimnames[[1]] 
YearNames     <- names.SimData[grep("year", names.SimData)]   #<-Change 'Site' for your model
YearNames

YearRows <- lapply(YearNames, MyID)
YearRows <- as.numeric(YearRows)
YearRows


# Start a loop to extract betas and random  effects, calculate
# the fitted values and simulate count data from 
# the model.
N    <- nrow(catchwb6_train)
Ysim <- matrix(nrow = N, ncol = NSim)
mu.i <- matrix(nrow = N, ncol = NSim)
X    <- model.matrix(~ soaktimec, data = catchwb6_train)
X   <- as.matrix(X)

YearID <- as.numeric(as.factor(catchwb6_train$year))


theta <- vector(length = NSim)        #Create space for simulated thetas



for (i in 1: NSim){
  Betas    <- SimData[[i]]$latent[BetaRows]       #Get simulated betas
  Year_i <- SimData[[i]]$latent[YearRows] #Get simulated random intercepts for year

  
  #eta = X * beta + year_i + Site_ij +
  eta <- X %*% Betas + Year_i[YearID] 
  
  #Get the simulated theta  
  theta[i] <- SimData[[i]]$hyperpar["size for the nbinomial observations (1/overdispersion)"]
  
  # Calculate mu = exp(eta)
  mu.i[,i] <- exp(eta)                        #Fitted values
  Ysim[,i] <- rnegbin(n = N,                  #Simulated count data
                      mu = mu.i[,i] ,
                      theta = theta[i]) 
}


# Now we have 1000 simulated data sets from the model.
# What shall we do with these simulated data sets?
# Calculate the dispersion for each simulated data set
DispNB <- vector(length = NSim) #Create space
N      <- nrow(catchwb6_train)              #Sample size 
Npar   <- length(Betas) + 2+ 1     #Number of regression parameters + theta
for(i in 1:NSim){
  ei <- (Ysim[,i] - mu.i[,i]) / sqrt(mu.i[,i] + mu.i[,i]^2 / theta[i]) #NB Pearson residuals
  DispNB[i] <- sum(ei^2) / (N - Npar)                                   #Dispersion statistic
}


# Let's plot this as a table
par(mfrow = c(1,1), mar = c(5,5,2,2), cex.lab = 1.5)
hist(DispNB, 
     xlim = c(0, 3),
     xlab = "Dispersion statistic",
     ylab = "Frequency",
     main = "Simulation results")
points(x = Dispersion.nb, 
       y = 0, 
       pch = 16, 
       cex = 5, 
       col = 2)

# NB GAM enviro full ------
I2 <- inla(nettot ~  soaktimec +  distance_to_southc + distance_to_shorec + temp_Cc + turb_FNUc+
             f(site, model='iid') +
             f(year, model='rw1'), 
           control.compute = list(dic = TRUE, cpo=TRUE, waic=TRUE, config=TRUE),
           family = "nbinomial",
           control.family = list(link='log'),
           data = catchwb6_train)
summary(I2)
-sum(log(I2$cpo$cpo)) # smaller is better fit of model to data

# Posterior mean values and 95% CI for the regression parameters:
Beta2<- I2$summary.fixed[,c("mean", "0.025quant", "0.975quant")] 
print(Beta2, digits = 2)  
# quantiles that DO NOT contain 0 indicated impt var
# all except intercept contain zero
HyperPar.I2 <- bri.hyperpar.summary(I2)
round(HyperPar.I2, digits = 3)

# posterior mean of theta
theta.pd <- I2$marginals.hyperpar$`size for the nbinomial observations (1/overdispersion)`
theta.pm <- inla.emarginal(function(x) x, theta.pd)
theta.pm ## 2.5
# posterior distribution of theta
par(mfrow = c(1,1), mar = c(5,5,2,2), cex.lab = 1.5)
plot(theta.pd, 
     type = "l",
     xlim = c(1,6),
     xlab = expression(paste(theta)),
     ylab = expression(paste("Pr(", theta," | data)")))

# overdispersed?
mu2 <- I2$summary.fitted.values[,"mean"]
E2  <- (catchwb6_train$nettot - mu2) /sqrt(mu2 + mu2^2 / theta.pm)
N <- nrow(catchwb6_train)
p <- nrow(I2$summary.fixed) 
Dispersion.nb <- sum(E2^2) / (N - p)
Dispersion.nb  


# why underdispersed
# fitted values versus the Pearson residuals
par(mfrow=c(2,1))
plot(x = mu2,
     y = E2,
     xlab = "Fitted values",
     ylab = "Pearson residuals",
     cex.lab = 1.5)
abline(h = 0, lty = 2)

boxplot(E2 ~ site, 
        ylab = "Pearson residuals",
        data = catchwb6_train)
abline(h = 0, lty = 2)

#residuals vs each covariate (in and not in the model).
# Check for non-linear patterns in the residuals.
MyVar <- c("turb_FNU", "sitedep_m", "temp_C", "DO_perc", "pH", "soaktime",  
           "distance_to_south", "distance_to_shore", "year")
catchwb6_train$E2 <- E2
MyMultipanel.ggp2(Z = catchwb6_train, 
                  varx = MyVar, 
                  vary = "E2", 
                  ylab = "Residuals",
                  addSmoother = TRUE,
                  addRegressionLine = FALSE,
                  addHorizontalLine = TRUE)
# non linearity present in all covariates, aside from temp d2S and soak

# spatial dependency?
catchwb6_train$Sign <- ifelse(E2 >= 0, "Positive residual", "Negative residual") 
ggplot() +
  theme(text = element_text(size=13)) +
  geom_point(aes(x = Xkm, y = Ykm, col = Sign), shape = 1, size = 4, data = catchwb6_train) +
  xlab("Xkm") +
  ylab("Ykm")  
# clumpy - suggest spatial dependency

MyData2 <- data.frame(E2 = E2, 
                      Xkm = catchwb6_train$Xkm, 
                      Ykm = catchwb6_train$Ykm)
coordinates(MyData2) <- c("Xkm", "Ykm")
V2 <- variogram(E2 ~ 1, 
                MyData2, 
                cressie = TRUE, 
                cutoff = 70)
ggplot() + 
  geom_point(data = V2, aes(x = dist, y= gamma)) +
  geom_smooth(data = V2, se = FALSE, span = 0.9, aes(x = dist, y = gamma)) +
  xlab("Distance") + ylab("Semi-variogram") +
  theme(text = element_text(size = 15)) +
  theme(legend.position="none") 
# no



### obs vs pred 2a ---------
I2a <- inla(nettot ~  soaktimec +  temp_Cc + turb_FNU + distance_to_southc + distance_to_shorec +
              f(year, model='rw1') +
              f(site, model='iid'),
            control.compute = list(dic = TRUE, cpo=TRUE, waic=TRUE, config=TRUE),
            family = "nbinomial",
            control.family = list(link='log'),
            data = catchwb6_pred) # now pred
summary(I2a) # should be the same as before

#extract the predictor values associated with the appended data
catchwb6_test2a <- cbind(catchwb6_test, 
                         I2a$summary.linear.predictor[(nrow(catchwb6_train)+1):nrow(catchwb6_pred),])
head(catchwb6_test2a)

catchwb6_test2a <- reshape:::rename(catchwb6_test2a, c("0.025quant"="lower", "0.975quant"="upper"))

ggplot(catchwb6_test2a, aes(y=mean, x=temp_Cc)) +
  geom_point(data=catchwb6_train, aes(y=nettot)) + # obs
  geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.2) + # predicted SD
  geom_line() # predicted
# MUCH smaller than observed

ggplot(catchwb6_test2a, aes(y=mean, x=temp_Cc)) +
  geom_point(aes(y=resp_obsN)) + # obs
  geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.2) + # predicted SD
  geom_line()


catchwb6_test2a$soaktime2 <- catchwb6_test2a$soaktimec * sd(catchwb6_train$soaktime) + mean(catchwb6_train$soaktime)
p <- ggplot()
p <- p + geom_point(data = catchwb6_train, 
                    aes(y = nettot, x = soaktime),
                    size = 1)
p <- p + xlab("soak time (h)") + ylab("total number LKWF")
p <- p + theme(text = element_text(size=15)) 
p <- p + geom_point(data = catchwb6_test2a, 
                    aes(x = soaktime2, 
                        y = mean), 
                    col = "red")
p <- p + geom_ribbon(data = catchwb6_test2a, 
                     aes(x = soaktime2, 
                         ymax = upper, 
                         ymin = lower),
                     alpha = 0.2)
p

ggplot(data=catchwb6_test2a) +
  geom_point(mapping=aes(x=resp_obsN, y=mean)) +
  xlab("Observed") +
  ylab("Predicted") 

# *NB GAM temp + turb ------
# not used
I2 <- inla(nettot ~  soaktimec + temp_Cc + turb_FNUc +
             f(site, model='iid') +
             f(year, model='rw1'), 
           control.compute = list(dic = TRUE, cpo=TRUE, waic=TRUE, config=TRUE),
           family = "nbinomial",
           control.family = list(link='log'),
           data = catchwb6_train)
summary(I2)
-sum(log(I2$cpo$cpo)) # smaller is better fit of model to data

# Posterior mean values and 95% CI for the regression parameters:
Beta2<- I2$summary.fixed[,c("mean", "0.025quant", "0.975quant")] 
print(Beta2, digits = 2)  
# quantiles that DO NOT contain 0 indicated impt var
# all except intercept contain zero


# posterior mean of theta
theta.pd <- I2$marginals.hyperpar$`size for the nbinomial observations (1/overdispersion)`
theta.pm <- inla.emarginal(function(x) x, theta.pd)
theta.pm ## 2.5
# posterior distribution of theta
par(mfrow = c(1,1), mar = c(5,5,2,2), cex.lab = 1.5)
plot(theta.pd, 
     type = "l",
     xlim = c(1,6),
     xlab = expression(paste(theta)),
     ylab = expression(paste("Pr(", theta," | data)")))

# overdispersed?
mu2 <- I2$summary.fitted.values[,"mean"]
E2  <- (catchwb6_train$nettot - mu2) /sqrt(mu2 + mu2^2 / theta.pm)
N <- nrow(catchwb6_train)
p <- nrow(I2$summary.fixed) 
Dispersion.nb <- sum(E2^2) / (N - p)
Dispersion.nb  


# why underdispersed
# fitted values versus the Pearson residuals
par(mfrow=c(2,1))
plot(x = mu2,
     y = E2,
     xlab = "Fitted values",
     ylab = "Pearson residuals",
     cex.lab = 1.5)
abline(h = 0, lty = 2)

boxplot(E2 ~ site, 
        ylab = "Pearson residuals",
        data = catchwb6_train)
abline(h = 0, lty = 2)

#residuals vs each covariate (in and not in the model).
# Check for non-linear patterns in the residuals.
MyVar <- c("turb_FNU", "sitedep_m", "temp_C", "DO_perc", "pH", "soaktime",  
           "distance_to_south", "distance_to_shore", "year")
catchwb6_train$E2 <- E2
MyMultipanel.ggp2(Z = catchwb6_train, 
                  varx = MyVar, 
                  vary = "E2", 
                  ylab = "Residuals",
                  addSmoother = TRUE,
                  addRegressionLine = FALSE,
                  addHorizontalLine = TRUE)
# non linearity present in all covariates, aside from temp d2S and soak

# spatial dependency?
catchwb6_train$Sign <- ifelse(E2 >= 0, "Positive residual", "Negative residual") 
ggplot() +
  theme(text = element_text(size=13)) +
  geom_point(aes(x = Xkm, y = Ykm, col = Sign), shape = 1, size = 4, data = catchwb6_train) +
  xlab("Xkm") +
  ylab("Ykm")  
# clumpy - suggest spatial dependency

MyData2 <- data.frame(E2 = E2, 
                      Xkm = catchwb6_train$Xkm, 
                      Ykm = catchwb6_train$Ykm)
coordinates(MyData2) <- c("Xkm", "Ykm")
V2 <- variogram(E2 ~ 1, 
                MyData2, 
                cressie = TRUE, 
                cutoff = 70)
ggplot() + 
  geom_point(data = V2, aes(x = dist, y= gamma)) +
  geom_smooth(data = V2, se = FALSE, span = 0.9, aes(x = dist, y = gamma)) +
  xlab("Distance") + ylab("Semi-variogram") +
  theme(text = element_text(size = 15)) +
  theme(legend.position="none") 
# no



### obs vs pred 2a ---------
I2a <- inla(nettot ~  soaktimec +  temp_Cc + turb_FNU + 
              f(year, model='rw1') +
              f(site, model='iid'),
            control.compute = list(dic = TRUE, cpo=TRUE, waic=TRUE, config=TRUE),
            family = "nbinomial",
            control.family = list(link='log'),
            data = catchwb6_pred) # now pred
summary(I2a) # should be the same as before

#extract the predictor values associated with the appended data
catchwb6_test2a <- cbind(catchwb6_test, 
                         I2a$summary.linear.predictor[(nrow(catchwb6_train)+1):nrow(catchwb6_pred),])
head(catchwb6_test2a)

catchwb6_test2a <- reshape:::rename(catchwb6_test2a, c("0.025quant"="lower", "0.975quant"="upper"))

ggplot(catchwb6_test2a, aes(y=mean, x=temp_Cc)) +
  geom_point(data=catchwb6_train, aes(y=nettot)) + # obs
  geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.2) + # predicted SD
  geom_line() # predicted
# MUCH smaller than observed

ggplot(catchwb6_test2a, aes(y=mean, x=temp_Cc)) +
  geom_point(aes(y=resp_obsN)) + # obs
  geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.2) + # predicted SD
  geom_line()


catchwb6_test2a$soaktime2 <- catchwb6_test2a$soaktimec * sd(catchwb6_train$soaktime) + mean(catchwb6_train$soaktime)
p <- ggplot()
p <- p + geom_point(data = catchwb6_train, 
                    aes(y = nettot, x = soaktime),
                    size = 1)
p <- p + xlab("soak time (h)") + ylab("total number LKWF")
p <- p + theme(text = element_text(size=15)) 
p <- p + geom_point(data = catchwb6_test2a, 
                    aes(x = soaktime2, 
                        y = mean), 
                    col = "red")
p <- p + geom_ribbon(data = catchwb6_test2a, 
                     aes(x = soaktime2, 
                         ymax = upper, 
                         ymin = lower),
                     alpha = 0.2)
p

ggplot(data=catchwb6_test2a) +
  geom_point(mapping=aes(x=resp_obsN, y=mean)) +
  xlab("Observed") +
  ylab("Predicted") 

# *NB GAM temp ------
# not used
I2 <- inla(nettot ~  soaktimec +   temp_Cc + 
             f(site, model='iid') +
             f(year, model='rw1'), 
           control.compute = list(dic = TRUE, cpo=TRUE, waic=TRUE, config=TRUE),
           family = "nbinomial",
           control.family = list(link='log'),
           data = catchwb6_train)
summary(I2)
-sum(log(I2$cpo$cpo)) # smaller is better fit of model to data

# Posterior mean values and 95% CI for the regression parameters:
Beta2<- I2$summary.fixed[,c("mean", "0.025quant", "0.975quant")] 
print(Beta2, digits = 2)  
# quantiles that DO NOT contain 0 indicated impt var
# all except intercept contain zero


# posterior mean of theta
theta.pd <- I2$marginals.hyperpar$`size for the nbinomial observations (1/overdispersion)`
theta.pm <- inla.emarginal(function(x) x, theta.pd)
theta.pm ## 2.5
# posterior distribution of theta
par(mfrow = c(1,1), mar = c(5,5,2,2), cex.lab = 1.5)
plot(theta.pd, 
     type = "l",
     xlim = c(1,6),
     xlab = expression(paste(theta)),
     ylab = expression(paste("Pr(", theta," | data)")))

# overdispersed?
mu2 <- I2$summary.fitted.values[,"mean"]
E2  <- (catchwb6_train$nettot - mu2) /sqrt(mu2 + mu2^2 / theta.pm)
N <- nrow(catchwb6_train)
p <- nrow(I2$summary.fixed) 
Dispersion.nb <- sum(E2^2) / (N - p)
Dispersion.nb  


# why underdispersed
# fitted values versus the Pearson residuals
par(mfrow=c(2,1))
plot(x = mu2,
     y = E2,
     xlab = "Fitted values",
     ylab = "Pearson residuals",
     cex.lab = 1.5)
abline(h = 0, lty = 2)

boxplot(E2 ~ site, 
        ylab = "Pearson residuals",
        data = catchwb6_train)
abline(h = 0, lty = 2)

#residuals vs each covariate (in and not in the model).
# Check for non-linear patterns in the residuals.
MyVar <- c("turb_FNU", "sitedep_m", "temp_C", "DO_perc", "pH", "soaktime",  
           "distance_to_south", "distance_to_shore", "year")
catchwb6_train$E2 <- E2
MyMultipanel.ggp2(Z = catchwb6_train, 
                  varx = MyVar, 
                  vary = "E2", 
                  ylab = "Residuals",
                  addSmoother = TRUE,
                  addRegressionLine = FALSE,
                  addHorizontalLine = TRUE)
# non linearity present in all covariates, aside from temp d2S and soak

# spatial dependency?
catchwb6_train$Sign <- ifelse(E2 >= 0, "Positive residual", "Negative residual") 
ggplot() +
  theme(text = element_text(size=13)) +
  geom_point(aes(x = Xkm, y = Ykm, col = Sign), shape = 1, size = 4, data = catchwb6_train) +
  xlab("Xkm") +
  ylab("Ykm")  
# clumpy - suggest spatial dependency

MyData2 <- data.frame(E2 = E2, 
                      Xkm = catchwb6_train$Xkm, 
                      Ykm = catchwb6_train$Ykm)
coordinates(MyData2) <- c("Xkm", "Ykm")
V2 <- variogram(E2 ~ 1, 
                MyData2, 
                cressie = TRUE, 
                cutoff = 70)
ggplot() + 
  geom_point(data = V2, aes(x = dist, y= gamma)) +
  geom_smooth(data = V2, se = FALSE, span = 0.9, aes(x = dist, y = gamma)) +
  xlab("Distance") + ylab("Semi-variogram") +
  theme(text = element_text(size = 15)) +
  theme(legend.position="none") 
# no



### obs vs pred 2a ---------
I2a <- inla(nettot ~  soaktimec +  temp_Cc +
              f(year, model='rw1') +
              f(site, model='iid'),
            control.compute = list(dic = TRUE, cpo=TRUE, waic=TRUE, config=TRUE),
            family = "nbinomial",
            control.family = list(link='log'),
            data = catchwb6_pred) # now pred
summary(I2a) # should be the same as before

#extract the predictor values associated with the appended data
catchwb6_test2a <- cbind(catchwb6_test, 
                         I2a$summary.linear.predictor[(nrow(catchwb6_train)+1):nrow(catchwb6_pred),])
head(catchwb6_test2a)

catchwb6_test2a <- reshape:::rename(catchwb6_test2a, c("0.025quant"="lower", "0.975quant"="upper"))

ggplot(catchwb6_test2a, aes(y=mean, x=temp_Cc)) +
  geom_point(data=catchwb6_train, aes(y=nettot)) + # obs
  geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.2) + # predicted SD
  geom_line() # predicted
# MUCH smaller than observed

ggplot(catchwb6_test2a, aes(y=mean, x=temp_Cc)) +
  geom_point(aes(y=resp_obsN)) + # obs
  geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.2) + # predicted SD
  geom_line()


catchwb6_test2a$soaktime2 <- catchwb6_test2a$soaktimec * sd(catchwb6_train$soaktime) + mean(catchwb6_train$soaktime)
p <- ggplot()
p <- p + geom_point(data = catchwb6_train, 
                    aes(y = nettot, x = soaktime),
                    size = 1)
p <- p + xlab("soak time (h)") + ylab("total number LKWF")
p <- p + theme(text = element_text(size=15)) 
p <- p + geom_point(data = catchwb6_test2a, 
                    aes(x = soaktime2, 
                        y = mean), 
                    col = "red")
p <- p + geom_ribbon(data = catchwb6_test2a, 
                     aes(x = soaktime2, 
                         ymax = upper, 
                         ymin = lower),
                     alpha = 0.2)
p

ggplot(data=catchwb6_test2a) +
  geom_point(mapping=aes(x=resp_obsN, y=mean)) +
  xlab("Observed") +
  ylab("Predicted") 

# *NB GAM temp d2S ------
# not used
I2 <- inla(nettot ~  soaktimec +  distance_to_southc + temp_Cc +
             f(site, model='iid') +
             f(year, model='rw1'), 
           control.compute = list(dic = TRUE, cpo=TRUE, waic=TRUE, config=TRUE),
           family = "nbinomial",
           control.family = list(link='log'),
           data = catchwb6_train)
summary(I2)
-sum(log(I2$cpo$cpo)) # smaller is better fit of model to data

# Posterior mean values and 95% CI for the regression parameters:
Beta2<- I2$summary.fixed[,c("mean", "0.025quant", "0.975quant")] 
print(Beta2, digits = 2)  
# quantiles that DO NOT contain 0 indicated impt var
# all except intercept contain zero


# posterior mean of theta
theta.pd <- I2$marginals.hyperpar$`size for the nbinomial observations (1/overdispersion)`
theta.pm <- inla.emarginal(function(x) x, theta.pd)
theta.pm ## 2.5
# posterior distribution of theta
par(mfrow = c(1,1), mar = c(5,5,2,2), cex.lab = 1.5)
plot(theta.pd, 
     type = "l",
     xlim = c(1,6),
     xlab = expression(paste(theta)),
     ylab = expression(paste("Pr(", theta," | data)")))

# overdispersed?
mu2 <- I2$summary.fitted.values[,"mean"]
E2  <- (catchwb6_train$nettot - mu2) /sqrt(mu2 + mu2^2 / theta.pm)
N <- nrow(catchwb6_train)
p <- nrow(I2$summary.fixed) 
Dispersion.nb <- sum(E2^2) / (N - p)
Dispersion.nb  


# why underdispersed
# fitted values versus the Pearson residuals
par(mfrow=c(2,1))
plot(x = mu2,
     y = E2,
     xlab = "Fitted values",
     ylab = "Pearson residuals",
     cex.lab = 1.5)
abline(h = 0, lty = 2)

boxplot(E2 ~ site, 
        ylab = "Pearson residuals",
        data = catchwb6_train)
abline(h = 0, lty = 2)

#residuals vs each covariate (in and not in the model).
# Check for non-linear patterns in the residuals.
MyVar <- c("turb_FNU", "sitedep_m", "temp_C", "DO_perc", "pH", "soaktime",  
           "distance_to_south", "distance_to_shore", "year")
catchwb6_train$E2 <- E2
MyMultipanel.ggp2(Z = catchwb6_train, 
                  varx = MyVar, 
                  vary = "E2", 
                  ylab = "Residuals",
                  addSmoother = TRUE,
                  addRegressionLine = FALSE,
                  addHorizontalLine = TRUE)
# non linearity present in all covariates, aside from temp d2S and soak

# spatial dependency?
catchwb6_train$Sign <- ifelse(E2 >= 0, "Positive residual", "Negative residual") 
ggplot() +
  theme(text = element_text(size=13)) +
  geom_point(aes(x = Xkm, y = Ykm, col = Sign), shape = 1, size = 4, data = catchwb6_train) +
  xlab("Xkm") +
  ylab("Ykm")  
# clumpy - suggest spatial dependency

MyData2 <- data.frame(E2 = E2, 
                      Xkm = catchwb6_train$Xkm, 
                      Ykm = catchwb6_train$Ykm)
coordinates(MyData2) <- c("Xkm", "Ykm")
V2 <- variogram(E2 ~ 1, 
                MyData2, 
                cressie = TRUE, 
                cutoff = 70)
ggplot() + 
  geom_point(data = V2, aes(x = dist, y= gamma)) +
  geom_smooth(data = V2, se = FALSE, span = 0.9, aes(x = dist, y = gamma)) +
  xlab("Distance") + ylab("Semi-variogram") +
  theme(text = element_text(size = 15)) +
  theme(legend.position="none") 
# no



### obs vs pred 2a ---------
I2a <- inla(nettot ~  soaktimec +  temp_Cc +  distance_to_southc +
              f(year, model='rw1') +
              f(site, model='iid'),
            control.compute = list(dic = TRUE, cpo=TRUE, waic=TRUE, config=TRUE),
            family = "nbinomial",
            control.family = list(link='log'),
            data = catchwb6_pred) # now pred
summary(I2a) # should be the same as before

#extract the predictor values associated with the appended data
catchwb6_test2a <- cbind(catchwb6_test, 
                         I2a$summary.linear.predictor[(nrow(catchwb6_train)+1):nrow(catchwb6_pred),])
head(catchwb6_test2a)

catchwb6_test2a <- reshape:::rename(catchwb6_test2a, c("0.025quant"="lower", "0.975quant"="upper"))

ggplot(catchwb6_test2a, aes(y=mean, x=temp_Cc)) +
  geom_point(data=catchwb6_train, aes(y=nettot)) + # obs
  geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.2) + # predicted SD
  geom_line() # predicted
# MUCH smaller than observed

ggplot(catchwb6_test2a, aes(y=mean, x=temp_Cc)) +
  geom_point(aes(y=resp_obsN)) + # obs
  geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.2) + # predicted SD
  geom_line()


catchwb6_test2a$soaktime2 <- catchwb6_test2a$soaktimec * sd(catchwb6_train$soaktime) + mean(catchwb6_train$soaktime)
p <- ggplot()
p <- p + geom_point(data = catchwb6_train, 
                    aes(y = nettot, x = soaktime),
                    size = 1)
p <- p + xlab("soak time (h)") + ylab("total number LKWF")
p <- p + theme(text = element_text(size=15)) 
p <- p + geom_point(data = catchwb6_test2a, 
                    aes(x = soaktime2, 
                        y = mean), 
                    col = "red")
p <- p + geom_ribbon(data = catchwb6_test2a, 
                     aes(x = soaktime2, 
                         ymax = upper, 
                         ymin = lower),
                     alpha = 0.2)
p

ggplot(data=catchwb6_test2a) +
  geom_point(mapping=aes(x=resp_obsN, y=mean)) +
  xlab("Observed") +
  ylab("Predicted") 

# *NB GAM tur d2S ------
# not used
I2 <- inla(nettot ~  soaktimec +  distance_to_southc +turb_FNUc+
             f(site, model='iid') +
             f(year, model='rw1'), 
           control.compute = list(dic = TRUE, cpo=TRUE, waic=TRUE, config=TRUE),
           family = "nbinomial",
           control.family = list(link='log'),
           data = catchwb6_train)
summary(I2)
-sum(log(I2$cpo$cpo)) # smaller is better fit of model to data

# Posterior mean values and 95% CI for the regression parameters:
Beta2<- I2$summary.fixed[,c("mean", "0.025quant", "0.975quant")] 
print(Beta2, digits = 2)  
# quantiles that DO NOT contain 0 indicated impt var
# all except intercept contain zero


# posterior mean of theta
theta.pd <- I2$marginals.hyperpar$`size for the nbinomial observations (1/overdispersion)`
theta.pm <- inla.emarginal(function(x) x, theta.pd)
theta.pm ## 2.5
# posterior distribution of theta
par(mfrow = c(1,1), mar = c(5,5,2,2), cex.lab = 1.5)
plot(theta.pd, 
     type = "l",
     xlim = c(1,6),
     xlab = expression(paste(theta)),
     ylab = expression(paste("Pr(", theta," | data)")))

# overdispersed?
mu2 <- I2$summary.fitted.values[,"mean"]
E2  <- (catchwb6_train$nettot - mu2) /sqrt(mu2 + mu2^2 / theta.pm)
N <- nrow(catchwb6_train)
p <- nrow(I2$summary.fixed) 
Dispersion.nb <- sum(E2^2) / (N - p)
Dispersion.nb  


# why underdispersed
# fitted values versus the Pearson residuals
par(mfrow=c(2,1))
plot(x = mu2,
     y = E2,
     xlab = "Fitted values",
     ylab = "Pearson residuals",
     cex.lab = 1.5)
abline(h = 0, lty = 2)

boxplot(E2 ~ site, 
        ylab = "Pearson residuals",
        data = catchwb6_train)
abline(h = 0, lty = 2)

#residuals vs each covariate (in and not in the model).
# Check for non-linear patterns in the residuals.
MyVar <- c("turb_FNU", "sitedep_m", "temp_C", "DO_perc", "pH", "soaktime",  
           "distance_to_south", "distance_to_shore", "year")
catchwb6_train$E2 <- E2
MyMultipanel.ggp2(Z = catchwb6_train, 
                  varx = MyVar, 
                  vary = "E2", 
                  ylab = "Residuals",
                  addSmoother = TRUE,
                  addRegressionLine = FALSE,
                  addHorizontalLine = TRUE)
# non linearity present in all covariates, aside from temp d2S and soak

# spatial dependency?
catchwb6_train$Sign <- ifelse(E2 >= 0, "Positive residual", "Negative residual") 
ggplot() +
  theme(text = element_text(size=13)) +
  geom_point(aes(x = Xkm, y = Ykm, col = Sign), shape = 1, size = 4, data = catchwb6_train) +
  xlab("Xkm") +
  ylab("Ykm")  
# clumpy - suggest spatial dependency

MyData2 <- data.frame(E2 = E2, 
                      Xkm = catchwb6_train$Xkm, 
                      Ykm = catchwb6_train$Ykm)
coordinates(MyData2) <- c("Xkm", "Ykm")
V2 <- variogram(E2 ~ 1, 
                MyData2, 
                cressie = TRUE, 
                cutoff = 70)
ggplot() + 
  geom_point(data = V2, aes(x = dist, y= gamma)) +
  geom_smooth(data = V2, se = FALSE, span = 0.9, aes(x = dist, y = gamma)) +
  xlab("Distance") + ylab("Semi-variogram") +
  theme(text = element_text(size = 15)) +
  theme(legend.position="none") 
# no



### obs vs pred 2a ---------
I2a <- inla(nettot ~  soaktimec +   turb_FNU + distance_to_southc +
              f(year, model='rw1') +
              f(site, model='iid'),
            control.compute = list(dic = TRUE, cpo=TRUE, waic=TRUE, config=TRUE),
            family = "nbinomial",
            control.family = list(link='log'),
            data = catchwb6_pred) # now pred
summary(I2a) # should be the same as before

#extract the predictor values associated with the appended data
catchwb6_test2a <- cbind(catchwb6_test, 
                         I2a$summary.linear.predictor[(nrow(catchwb6_train)+1):nrow(catchwb6_pred),])
head(catchwb6_test2a)

catchwb6_test2a <- reshape:::rename(catchwb6_test2a, c("0.025quant"="lower", "0.975quant"="upper"))

ggplot(catchwb6_test2a, aes(y=mean, x=temp_Cc)) +
  geom_point(data=catchwb6_train, aes(y=nettot)) + # obs
  geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.2) + # predicted SD
  geom_line() # predicted
# MUCH smaller than observed

ggplot(catchwb6_test2a, aes(y=mean, x=temp_Cc)) +
  geom_point(aes(y=resp_obsN)) + # obs
  geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.2) + # predicted SD
  geom_line()


catchwb6_test2a$soaktime2 <- catchwb6_test2a$soaktimec * sd(catchwb6_train$soaktime) + mean(catchwb6_train$soaktime)
p <- ggplot()
p <- p + geom_point(data = catchwb6_train, 
                    aes(y = nettot, x = soaktime),
                    size = 1)
p <- p + xlab("soak time (h)") + ylab("total number LKWF")
p <- p + theme(text = element_text(size=15)) 
p <- p + geom_point(data = catchwb6_test2a, 
                    aes(x = soaktime2, 
                        y = mean), 
                    col = "red")
p <- p + geom_ribbon(data = catchwb6_test2a, 
                     aes(x = soaktime2, 
                         ymax = upper, 
                         ymin = lower),
                     alpha = 0.2)
p

ggplot(data=catchwb6_test2a) +
  geom_point(mapping=aes(x=resp_obsN, y=mean)) +
  xlab("Observed") +
  ylab("Predicted") 

# *NB GAM turb ------
# not used
I2 <- inla(nettot ~  soaktimec +   turb_FNUc+
             f(site, model='iid') +
             f(year, model='rw1'), 
           control.compute = list(dic = TRUE, cpo=TRUE, waic=TRUE, config=TRUE),
           family = "nbinomial",
           control.family = list(link='log'),
           data = catchwb6_train)
summary(I2)
-sum(log(I2$cpo$cpo)) # smaller is better fit of model to data

# Posterior mean values and 95% CI for the regression parameters:
Beta2<- I2$summary.fixed[,c("mean", "0.025quant", "0.975quant")] 
print(Beta2, digits = 2)  
# quantiles that DO NOT contain 0 indicated impt var
# all except intercept contain zero


# posterior mean of theta
theta.pd <- I2$marginals.hyperpar$`size for the nbinomial observations (1/overdispersion)`
theta.pm <- inla.emarginal(function(x) x, theta.pd)
theta.pm ## 2.5
# posterior distribution of theta
par(mfrow = c(1,1), mar = c(5,5,2,2), cex.lab = 1.5)
plot(theta.pd, 
     type = "l",
     xlim = c(1,6),
     xlab = expression(paste(theta)),
     ylab = expression(paste("Pr(", theta," | data)")))

# overdispersed?
mu2 <- I2$summary.fitted.values[,"mean"]
E2  <- (catchwb6_train$nettot - mu2) /sqrt(mu2 + mu2^2 / theta.pm)
N <- nrow(catchwb6_train)
p <- nrow(I2$summary.fixed) 
Dispersion.nb <- sum(E2^2) / (N - p)
Dispersion.nb  


# why underdispersed
# fitted values versus the Pearson residuals
par(mfrow=c(2,1))
plot(x = mu2,
     y = E2,
     xlab = "Fitted values",
     ylab = "Pearson residuals",
     cex.lab = 1.5)
abline(h = 0, lty = 2)

boxplot(E2 ~ site, 
        ylab = "Pearson residuals",
        data = catchwb6_train)
abline(h = 0, lty = 2)

#residuals vs each covariate (in and not in the model).
# Check for non-linear patterns in the residuals.
MyVar <- c("turb_FNU", "sitedep_m", "temp_C", "DO_perc", "pH", "soaktime",  
           "distance_to_south", "distance_to_shore", "year")
catchwb6_train$E2 <- E2
MyMultipanel.ggp2(Z = catchwb6_train, 
                  varx = MyVar, 
                  vary = "E2", 
                  ylab = "Residuals",
                  addSmoother = TRUE,
                  addRegressionLine = FALSE,
                  addHorizontalLine = TRUE)
# non linearity present in all covariates, aside from temp d2S and soak

# spatial dependency?
catchwb6_train$Sign <- ifelse(E2 >= 0, "Positive residual", "Negative residual") 
ggplot() +
  theme(text = element_text(size=13)) +
  geom_point(aes(x = Xkm, y = Ykm, col = Sign), shape = 1, size = 4, data = catchwb6_train) +
  xlab("Xkm") +
  ylab("Ykm")  
# clumpy - suggest spatial dependency

MyData2 <- data.frame(E2 = E2, 
                      Xkm = catchwb6_train$Xkm, 
                      Ykm = catchwb6_train$Ykm)
coordinates(MyData2) <- c("Xkm", "Ykm")
V2 <- variogram(E2 ~ 1, 
                MyData2, 
                cressie = TRUE, 
                cutoff = 70)
ggplot() + 
  geom_point(data = V2, aes(x = dist, y= gamma)) +
  geom_smooth(data = V2, se = FALSE, span = 0.9, aes(x = dist, y = gamma)) +
  xlab("Distance") + ylab("Semi-variogram") +
  theme(text = element_text(size = 15)) +
  theme(legend.position="none") 
# no



### obs vs pred 2a ---------
I2a <- inla(nettot ~  soaktimec +  turb_FNU + 
              f(year, model='rw1') +
              f(site, model='iid'),
            control.compute = list(dic = TRUE, cpo=TRUE, waic=TRUE, config=TRUE),
            family = "nbinomial",
            control.family = list(link='log'),
            data = catchwb6_pred) # now pred
summary(I2a) # should be the same as before

#extract the predictor values associated with the appended data
catchwb6_test2a <- cbind(catchwb6_test, 
                         I2a$summary.linear.predictor[(nrow(catchwb6_train)+1):nrow(catchwb6_pred),])
head(catchwb6_test2a)

catchwb6_test2a <- reshape:::rename(catchwb6_test2a, c("0.025quant"="lower", "0.975quant"="upper"))

ggplot(catchwb6_test2a, aes(y=mean, x=temp_Cc)) +
  geom_point(data=catchwb6_train, aes(y=nettot)) + # obs
  geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.2) + # predicted SD
  geom_line() # predicted
# MUCH smaller than observed

ggplot(catchwb6_test2a, aes(y=mean, x=temp_Cc)) +
  geom_point(aes(y=resp_obsN)) + # obs
  geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.2) + # predicted SD
  geom_line()


catchwb6_test2a$soaktime2 <- catchwb6_test2a$soaktimec * sd(catchwb6_train$soaktime) + mean(catchwb6_train$soaktime)
p <- ggplot()
p <- p + geom_point(data = catchwb6_train, 
                    aes(y = nettot, x = soaktime),
                    size = 1)
p <- p + xlab("soak time (h)") + ylab("total number LKWF")
p <- p + theme(text = element_text(size=15)) 
p <- p + geom_point(data = catchwb6_test2a, 
                    aes(x = soaktime2, 
                        y = mean), 
                    col = "red")
p <- p + geom_ribbon(data = catchwb6_test2a, 
                     aes(x = soaktime2, 
                         ymax = upper, 
                         ymin = lower),
                     alpha = 0.2)
p

ggplot(data=catchwb6_test2a) +
  geom_point(mapping=aes(x=resp_obsN, y=mean)) +
  xlab("Observed") +
  ylab("Predicted") 

# *NB GAM d2S------
# not used
I2 <- inla(nettot ~  soaktimec +  distance_to_southc +
             f(site, model='iid') +
             f(year, model='rw1'), 
           control.compute = list(dic = TRUE, cpo=TRUE, waic=TRUE, config=TRUE),
           family = "nbinomial",
           control.family = list(link='log'),
           data = catchwb6_train)
summary(I2)
-sum(log(I2$cpo$cpo)) # smaller is better fit of model to data

# Posterior mean values and 95% CI for the regression parameters:
Beta2<- I2$summary.fixed[,c("mean", "0.025quant", "0.975quant")] 
print(Beta2, digits = 2)  
# quantiles that DO NOT contain 0 indicated impt var
# all except intercept contain zero


# posterior mean of theta
theta.pd <- I2$marginals.hyperpar$`size for the nbinomial observations (1/overdispersion)`
theta.pm <- inla.emarginal(function(x) x, theta.pd)
theta.pm ## 2.5
# posterior distribution of theta
par(mfrow = c(1,1), mar = c(5,5,2,2), cex.lab = 1.5)
plot(theta.pd, 
     type = "l",
     xlim = c(1,6),
     xlab = expression(paste(theta)),
     ylab = expression(paste("Pr(", theta," | data)")))

# overdispersed?
mu2 <- I2$summary.fitted.values[,"mean"]
E2  <- (catchwb6_train$nettot - mu2) /sqrt(mu2 + mu2^2 / theta.pm)
N <- nrow(catchwb6_train)
p <- nrow(I2$summary.fixed) 
Dispersion.nb <- sum(E2^2) / (N - p)
Dispersion.nb  


# why underdispersed
# fitted values versus the Pearson residuals
par(mfrow=c(2,1))
plot(x = mu2,
     y = E2,
     xlab = "Fitted values",
     ylab = "Pearson residuals",
     cex.lab = 1.5)
abline(h = 0, lty = 2)

boxplot(E2 ~ site, 
        ylab = "Pearson residuals",
        data = catchwb6_train)
abline(h = 0, lty = 2)

#residuals vs each covariate (in and not in the model).
# Check for non-linear patterns in the residuals.
MyVar <- c("turb_FNU", "sitedep_m", "temp_C", "DO_perc", "pH", "soaktime",  
           "distance_to_south", "distance_to_shore", "year")
catchwb6_train$E2 <- E2
MyMultipanel.ggp2(Z = catchwb6_train, 
                  varx = MyVar, 
                  vary = "E2", 
                  ylab = "Residuals",
                  addSmoother = TRUE,
                  addRegressionLine = FALSE,
                  addHorizontalLine = TRUE)
# non linearity present in all covariates, aside from temp d2S and soak

# spatial dependency?
catchwb6_train$Sign <- ifelse(E2 >= 0, "Positive residual", "Negative residual") 
ggplot() +
  theme(text = element_text(size=13)) +
  geom_point(aes(x = Xkm, y = Ykm, col = Sign), shape = 1, size = 4, data = catchwb6_train) +
  xlab("Xkm") +
  ylab("Ykm")  
# clumpy - suggest spatial dependency

MyData2 <- data.frame(E2 = E2, 
                      Xkm = catchwb6_train$Xkm, 
                      Ykm = catchwb6_train$Ykm)
coordinates(MyData2) <- c("Xkm", "Ykm")
V2 <- variogram(E2 ~ 1, 
                MyData2, 
                cressie = TRUE, 
                cutoff = 70)
ggplot() + 
  geom_point(data = V2, aes(x = dist, y= gamma)) +
  geom_smooth(data = V2, se = FALSE, span = 0.9, aes(x = dist, y = gamma)) +
  xlab("Distance") + ylab("Semi-variogram") +
  theme(text = element_text(size = 15)) +
  theme(legend.position="none") 
# no



### obs vs pred 2a ---------
I2a <- inla(nettot ~  soaktimec +  distance_to_southc +
              f(year, model='rw1') +
              f(site, model='iid'),
            control.compute = list(dic = TRUE, cpo=TRUE, waic=TRUE, config=TRUE),
            family = "nbinomial",
            control.family = list(link='log'),
            data = catchwb6_pred) # now pred
summary(I2a) # should be the same as before

#extract the predictor values associated with the appended data
catchwb6_test2a <- cbind(catchwb6_test, 
                         I2a$summary.linear.predictor[(nrow(catchwb6_train)+1):nrow(catchwb6_pred),])
head(catchwb6_test2a)

catchwb6_test2a <- reshape:::rename(catchwb6_test2a, c("0.025quant"="lower", "0.975quant"="upper"))

ggplot(catchwb6_test2a, aes(y=mean, x=temp_Cc)) +
  geom_point(data=catchwb6_train, aes(y=nettot)) + # obs
  geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.2) + # predicted SD
  geom_line() # predicted
# MUCH smaller than observed

ggplot(catchwb6_test2a, aes(y=mean, x=temp_Cc)) +
  geom_point(aes(y=resp_obsN)) + # obs
  geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.2) + # predicted SD
  geom_line()


catchwb6_test2a$soaktime2 <- catchwb6_test2a$soaktimec * sd(catchwb6_train$soaktime) + mean(catchwb6_train$soaktime)
p <- ggplot()
p <- p + geom_point(data = catchwb6_train, 
                    aes(y = nettot, x = soaktime),
                    size = 1)
p <- p + xlab("soak time (h)") + ylab("total number LKWF")
p <- p + theme(text = element_text(size=15)) 
p <- p + geom_point(data = catchwb6_test2a, 
                    aes(x = soaktime2, 
                        y = mean), 
                    col = "red")
p <- p + geom_ribbon(data = catchwb6_test2a, 
                     aes(x = soaktime2, 
                         ymax = upper, 
                         ymin = lower),
                     alpha = 0.2)
p

ggplot(data=catchwb6_test2a) +
  geom_point(mapping=aes(x=resp_obsN, y=mean)) +
  xlab("Observed") +
  ylab("Predicted") 




# NB GAM ind SPDE------
# location exploration
Loc <- cbind(catchwb6_train$Xkm, catchwb6_train$Ykm)
D <- dist(Loc)
par(mfrow = c(1,2), mar = c(5,5,2,2), cex.lab = 1.5)
hist(D, 
     freq = TRUE,
     main = "", 
     xlab = "Distance between sites (km)",
     ylab = "Frequency")

plot(x = sort(D), 
     y = (1:length(D))/length(D), 
     type = "l",
     xlab = "Distance between sites (km)",
     ylab = "Cumulative proportion")


# Create a mesh
RangeGuess <- 15
MaxEdge  <- RangeGuess / 3.5
ConvHull <- inla.nonconvex.hull(Loc)		
mesh    <- inla.mesh.2d(loc = Loc,
                        boundary = ConvHull,
                        max.edge = c(1, 5) * MaxEdge, 
                        cutoff  = MaxEdge / 5)


par(mfrow = c(1, 1), mar=c(0, 0, 2, 0))
plot(mesh)
points(Loc, col = 1, pch = 16, cex = 1.5)

# 2. Define the weighting factors a_ik (also called 
#    the projector matrix).
A <- inla.spde.make.A(mesh, loc = Loc)

# 3. Define the SPDE.
range(catchwb6_train$nettot)
# P(Range < 5 km ) = 0.05
# P(sigma > 2) = 0.05
# SB = exp(u_i)
# some u_i have to be as large as 5.24 to cover 190
# If u_i ~ N(0, sigma_u^2) then it is unlikley that sigma_u > 1.5
#P(sigma > 1.5) = 0.05
M1 <- glm(nettot ~ 1, data = catchwb6_train)
summary(M1)


spde <- inla.spde2.pcmatern(mesh, 
                            prior.range = c(15 , 0.05), 
                            prior.sigma = c(1,  0.05))



# 4. Define the spatial field.
w.index <- inla.spde.make.index(
  name    = 'w', 
  n.spde  = spde$n.spde)

# The size of the mesh is: 
mesh$n




# 5. Make a stack. In this process we tell INLA 
#    at which points on the mesh we sampled the 
#    response variable and the covariates. We 
#    also need to inform INLA (via the stack) at 
#    which points of the mesh we have any other 
#    terms, e.g. the random effects as we know 
#    them from Chapter 8.


# Doing this in INLA means that we get the 95% CI for the smoother.

N    <- nrow(catchwb6_train)
Xcov <- model.matrix(~ soaktimec + year,
                     data = catchwb6_train)
Xcov2 <- Xcov[,-1]

Stack <- inla.stack(
  tag = "Fit",
  data = list(nettot = catchwb6_train$nettot),  
  A = list(A, 1, 1, 1),                 
  effects = list(                 
    w  = w.index,  
    soaktimec = Xcov2[,1],
    year = Xcov2[,2],               
    site = catchwb6_train$site))  


f14 <- formula(catchwb6_train$nettot ~ 1 +
                 soaktimec  +
                 f (year, model='rw1') + 
                 f(site, model = "iid") +
                 f(w, model = spde))


# 7. Run the spatial model in INLA.
I14 <- inla(f14,
            family = "nbinomial", 
            data=inla.stack.data(Stack),
            control.compute = list(dic = TRUE, waic = TRUE, cpo= TRUE, config = TRUE),
            control.family=list(link='log'),
            control.predictor = list(A = inla.stack.A(Stack)))

summary(I14)
-sum(log(I14$cpo$cpo))

# Posterior mean values and 95% CI for the regression parameters:
Beta14 <- I14$summary.fixed[,c("mean", "sd", "0.025quant", "0.975quant")] 
print(Beta14, digits = 2)  

# The hyperparameters:
summary(I14)
HyperPar.I14 <- bri.hyperpar.summary(I14) # brinla function converts precision to sd
round(HyperPar.I14, digits = 3)

Year.pm <- I14$summary.random$year[,c("ID", "mean", "0.025quant", "0.975quant")]
Year.pm

site.pm <- I14$summary.random$site[,c("ID", "mean", "0.025quant", "0.975quant")]




#Here is a visualisation of the spatial correlation.
# Access the posterior mean values of the w.
w.pm <- I14$summary.random$w$mean  
w.sd <- I14$summary.random$w$sd  

# This is a vector of length mesh$n by 1.We can use INLA functions to 
# plot these w.pm.
PlotField2(field = w.pm, 
           mesh = mesh, 
           xlim = range(mesh$loc[,1]), 
           ylim = range(mesh$loc[,2]),
           MyTitle = "Spatial random field")


# Add the sampling locations (in UTM)
points(x = Loc[,1],
       y = Loc[,2], 
       cex = 0.5, 
       col = "black", 
       pch = 16)
# conclusion: 1 area of higher (dark red) expected nettot values and 12 area of lower (dark blue)
# these effects are labeled as spatial correlation as multple neighboting sites show this effect (blue)

# gproj <- inla.mesh.projector(mesh,  dims = c(300, 300))
# g.mean <- inla.mesh.project(gproj, I14$summary.random$w$mean)
# g.sd <- inla.mesh.project(gproj, I14$summary.random$w$sd)
# library(gridExtra)
# grid.arrange(levelplot(g.mean, scales=list(draw=F), xlab='', ylab='', main='mean',col.regions = topo.colors(16)),
#              levelplot(g.sd, scal=list(draw=F), xla='', yla='', main='sd' ,col.regions = topo.colors(16)), nrow=1)


               # levelplot, 
               # at=c(-2, -1.75,-1.5, -1.25, -1,-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75,2), # this was added to put on same scale
               # xlab='', 
               # ylab='',
               # col.regions=topo.colors(16), 
               # scale=list(draw=FALSE)))
## either this is real correlation or I missed a covariate

w.pm        <- I14$summary.random$w$mean  #mu.srf # exp converts from the log link (value of 1 is = exp(1), but now all exp())
w.proj      <- inla.mesh.projector(mesh) 
w.pm100_100 <- inla.mesh.project(w.proj, w.pm)
long        <- as.data.frame(w.proj$lattice$loc)
colnames(long) <- c("EastingKM", "NorthingKM")
long <- long %>%
  mutate(w.pm = as.vector(w.pm100_100))
long


mean <- ggplot() +
  geom_tile(long, mapping=aes(x = EastingKM, y = NorthingKM, 
                fill = w.pm)) +
  geom_point(aes(x=Xkm, y=Ykm), data=catchwb6) +
  xlab("Easting (km)") + 
  ylab("Northing (km)") +
  geom_contour(long, mapping=aes(x = EastingKM, y = NorthingKM, z=w.pm), 
                                  binwidth = 0.025, colour = "black") +
  scale_fill_gradient2(low='darkblue', mid='white', high='darkred') +
  labs(fill="ln(mean SRF)") +
  theme( #axis.text = element_blank(),
    panel.background = element_blank(), # Remove panel background
    axis.line = element_line(colour = "grey")) # Add axis lines
mean

mean + 
  geom_point(aes(x=Xkm, y=Ykm), data=catchwb6_train) 
  
table(catchwb6_train$site)
# hot spots don't seem to be result of more  data

  
w.sd <- I14$summary.random$w$sd 
w.sd100_100 <- inla.mesh.project(w.proj, w.sd)
long        <- as.data.frame(w.proj$lattice$loc)
colnames(long) <- c("EastingKM", "NorthingKM")
long <- long %>%
  mutate(w.sd = as.vector(w.sd100_100))

sd <- ggplot(long, 
            aes(x = EastingKM, 
                y = NorthingKM, 
                fill = w.sd, 
                z = w.sd)) +
  geom_tile() +
  xlab("Easting (km)") + 
  ylab("Northing (km)") +
  geom_contour(binwidth = 0.025, colour = "black") +
  scale_fill_gradient2(low='darkblue', mid='white', high='darkred') +
  labs(fill="ln(sd SRF)")

sd 

grid.arrange(mean, sd, nrow=2)

# how dar for the correlation go?
SpFi.w <- inla.spde2.result(inla=I14, name="w", spde=spde, do.transfer=T) # CHANGE
Kappa <- inla.emarginal(function (x) x, SpFi.w$marginals.kappa[[1]])
Sigma_u <- inla.emarginal(function(x) sqrt(x), SpFi.w$marginals.variance.nominal[[1]])
Range <- inla.emarginal(function(x) x, SpFi.w$marginals.range.nominal[[1]]) 

Sigma_u  #Sigma of the spatial random field.
Range    #Distance at which the correlation diminishes.
LocMesh <- mesh$loc[,1:2]

# Calculate the distance between each vertex.
D <- as.matrix(dist(LocMesh))

# Using the estimated parameters from the model (see above)
# we can calculate the imposed Matern correlation values.
d.vec <- seq(0, max(D), length = 100)      
Cor.M <- (Kappa * d.vec) * besselK(Kappa * d.vec, 1) 
Cor.M[1] <- 1

# Plot the results
par(mfrow = c(1,1), mar = c(5,5,2,2))
plot(x = d.vec, 
     y = Cor.M, 
     pch = 16, 
     type = "l", 
     cex.lab = 1.5,
     xlab = "Distance", 
     ylab = "Correlation",
     xlim = c(0, 80))
abline(h=c(0.8, 0.6, 0.2), lty = 2)
abline(v=c(6, 10.5, 26), col='red')

abline(h=0.6, v=10.5, lty=2)


# How many sites are separated by 0 - 10.5km?
Loc <- cbind(catchwb6_train$Xkm, catchwb6_train$Ykm)
D   <- dist(Loc)
par(mfrow = c(1,1), mar = c(5,5,2,2), cex.lab = 1.5)
hist(D, 
     freq = TRUE,
     main = "", 
     xlab = "Distance between sites (km)",
     ylab = "Frequency")
abline(v = c(10.5), col = 2, lty = 2, lwd = 2)
# no many sites are affected by this correlation!

# posterior mean of theta
theta.pd <- I14$marginals.hyperpar$`size for the nbinomial observations (1/overdispersion)`
theta.pm <- inla.emarginal(function(x) x, theta.pd)
theta.pm 

# posterior distribution of theta
par(mfrow = c(1,1), mar = c(5,5,2,2), cex.lab = 1.5)
plot(theta.pd, 
     type = "l",
     xlim = c(1,7),
     xlab = expression(paste(theta)),
     ylab = expression(paste("Pr(", theta," | data)")))

# Here are the fitted values and Pearson residuals
N <- nrow(catchwb6_train)
mu14 <- I14$summary.fitted.values[1:N,"mean"]           #Fitted values
E14  <- (catchwb6_train$nettot - mu14) /sqrt(mu14 + mu14^2 / theta.pm) #Pearson residuals

p <- nrow(I14$summary.fixed[, c("mean", "0.025quant", "0.975quant")])
Dispersion14 <- sum(E14^2) / (N - p)
Dispersion14
# 0.87


# why underdispersed
# fitted values versus the Pearson residuals
par(mfrow=c(2,1))
plot(x = mu14,
     y = E14,
     xlab = "Fitted values",
     ylab = "Pearson residuals",
     cex.lab = 1.5)
abline(h = 0, lty = 2)

par(mfrow=c(2,1))
plot(x = mu14,
     y = catchwb6_train$nettot,
     xlab = "Fitted values",
     ylab = "observed",
     cex.lab = 1.5)
abline(h = 0, lty = 2)

# Check the Pearson residuals of I2 for any remaining spatial dependency.



MyData14 <- data.frame(E14 = E14, 
                       Xkm = catchwb6_train$Xkm, 
                       Ykm = catchwb6_train$Ykm)
coordinates(MyData14) <- c("Xkm", "Ykm")
V14 <- variogram(E14 ~ 1, 
                 MyData14, 
                 cressie = TRUE, 
                 cutoff = 70)

# Plot the variogram
p <- ggplot()
p <- p + geom_point(data = V14,
                    aes(x = dist, 
                        y = gamma))
p <- p + geom_smooth(data = V13, se = FALSE, span = 0.9, aes(x = dist, y = gamma)) 
p <- p + xlab("Distance") + ylab("Semi-variogram")
p <- p + theme(text = element_text(size = 15)) 
p <- p + theme(legend.position="none") 
p

#residuals vs each covariate (in and not in the model).
# Check for non-linear patterns in the residuals.
MyVar <- c("turb_FNU", "sitedep_m", "temp_C", "DO_perc", "pH", "soaktime",  
           "distance_to_south", "distance_to_shore", "year")
catchwb6_train$E14 <- E14
MyMultipanel.ggp2(Z = catchwb6_train, 
                  varx = MyVar, 
                  vary = "E14", 
                  ylab = "Residuals",
                  addSmoother = TRUE,
                  addRegressionLine = FALSE,
                  addHorizontalLine = TRUE)
# non linearity present in all covariates, aside from temp d2S and soak

# spatial dependency?
catchwb6_train$Sign <- ifelse(E14 >= 0, "Positive residual", "Negative residual") 
ggplot() +
  theme(text = element_text(size=13)) +
  geom_point(aes(x = Xkm, y = Ykm, col = Sign), shape = 1, size = 4, data = catchwb6_train) +
  xlab("Xkm") +
  ylab("Ykm")  
# clumpy - suggest spatial dependency

boxplot(E14 ~ site, 
        ylab = "Pearson residuals",
        data = catchwb6_train)
abline(h = 0, lty = 2)


### obs vs pred 2a ---------
Loc <- cbind(catchwb6_train$Xkm, catchwb6_train$Ykm)
D <- dist(Loc)

# Create a mesh
RangeGuess <- 15
MaxEdge  <- RangeGuess / 3.5
ConvHull <- inla.nonconvex.hull(Loc)		
mesh    <- inla.mesh.2d(loc = Loc,
                        boundary = ConvHull,
                        max.edge = c(1, 5) * MaxEdge, 
                        cutoff  = MaxEdge / 5)


par(mfrow = c(1, 1), mar=c(0, 0, 2, 0))
plot(mesh)
points(Loc, col = 1, pch = 16, cex = 1.5)

# 2. Define the weighting factors a_ik (also called 
#    the projector matrix).
A <- inla.spde.make.A(mesh, loc = Loc)
dim(A)
# 3. Define the SPDE.



spde <- inla.spde2.pcmatern(mesh, 
                            prior.range = c(15 , 0.05), 
                            prior.sigma = c(1,  0.05))



# 4. Define the spatial field.
w.index <- inla.spde.make.index(
  name    = 'w', 
  n.spde  = spde$n.spde)


Xcov <- model.matrix(~ soaktimec + year,
                     data = catchwb6_train)
Xcov2 <- as.data.frame(Xcov)
colnames(Xcov2)

N    <- nrow(catchwb6_train)

Stack <- inla.stack(
  tag = "Fit",
  data = list(y = catchwb6_train$nettot),
  A = list(1,1,1, A),
  effects = list(
    Intercept = rep(1,N),
    Xcov2 = Xcov2[,-1],
    site = catchwb6_train$site, 
    w  = w.index))


f2 <- y ~ -1 + Intercept +
                 soaktimec  +
                 f (year, model='rw1') +
                 f(site, model = "iid") +
                 f(w, model = spde)


I2 <- inla(f2,
            family = "nbinomial",
            data=inla.stack.data(Stack),
            control.compute = list(dic = TRUE, waic = TRUE, config = TRUE),
            control.family=list(link='log'),
            control.predictor = list(A = inla.stack.A(Stack)))

summary(I2)
# 5. Make a stack. In this process we tell INLA 
#    at which points on the mesh we sampled the 
#    response variable and the covariates. We 
#    also need to inform INLA (via the stack) at 
#    which points of the mesh we have any other 
#    terms, e.g. the random effects as we know 
#    them from Chapter 8.

# 1. Create a grid of covariate values for which we want to predict richness.
# MyData <- ddply(catchwb6_test, 
#                 .(fYear), 
#                 summarize,
#                 LogMining.std = seq(min(LogMining.std), 
#                                     max(LogMining.std), 
#                                     length = 50))
# head(MyData)  #This object has 10*50 rows and 2 columns

# equivalent is using my test data




# 2. Make a design matrix Xp

# 2. Make a design matrix Xp
Xp <- model.matrix(~ soaktimec + year, data = catchwb6_test)

# Drop the first column (intercept) from this object, and make sure
# it is a data frame.
Xp2 <- as.data.frame(Xp[,-1])
head(Xp2)



# From here onward, things diverge a little bit from the models
# without spatial correlation. This is due to the 'stack'.

#  1. Make a stack for the data for which predictions are needed.
#  2. Combine this new stack with the stack for the observed data.
#  3. Run INLA with the combined stack.
#  4. Extract the relevant pieces and plot them.


# Step 1. Make a stack for the data in MyData. Recall that MyData was 
#         converted into Xp.
StackPred <- inla.stack(
  tag = "Predictions",
  data = list(y = NA),  
  A = list(1, 1, 1), # a 1 for each effect                 
  effects = list(
    Intercept = rep(1, nrow(Xp2)),
    Xp2 = Xp2,
    site = catchwb6_test$site))   #Covariate values without the intercept

# This stack does not have a spatial field, because we
# only want to have the predictions for the covariates and the 
# intercept. Note that the response variable equals NA. As a result INLA
# will predict these values.


# Step 2: Combine the two stacks.         
Combined.stacks <- inla.stack(Stack, StackPred)	              


# Step 3: Rerun INLA
# This is still the same:
f5 <- y ~ -1 + Intercept +
                 soaktimec  +
                 f (year, model='rw1') + 
                 f(site, model = "iid") +
                 f(w, model = spde)

# Note that we use a different stack! The rest is the same
I5.pred <- inla(f5,
                family = "nbinomial", 
                data = inla.stack.data(Combined.stacks),
                control.compute = list(dic = TRUE, waic=TRUE, config=TRUE),
                control.family=list(link='log'),
                control.predictor = list(
                  link = 1,
                  A = inla.stack.A(Combined.stacks)))


# link = 1: Ensures that for the NAs the correct 
#           link function is used.
# Posterior mean values and 95% CI for the regression parameters:
Beta5.pred<- I5.pred$summary.fixed[,c("mean", "0.025quant", "0.975quant")] 
print(Beta5.pred, digits = 2)

summary(I5.pred)
HyperPar.I5.pred <- bri.hyperpar.summary(I5.pred)
round(HyperPar.I5.pred, digits = 3)

# Step 4: Extract the relevant pieces and plot the results.
# This is a crucial piece of crucial for extracting the correct rows.
# It provides an index for which rows in the combined stack
# belong to the observed data and which rows to the 
# artificial covariate data. 
index.Fit <- inla.stack.index(Combined.stacks,
                              tag = "Fit")$data

index.Pred <- inla.stack.index(Combined.stacks,
                               tag = "Predictions")$data

#And we can extact the correct rows     
F5.fit  <- I5.pred$summary.fitted.values[index.Fit, c(1,3,5)]  
F5.pred <- I5.pred$summary.fitted.values[index.Pred, c(1,3,5)]  


# It is the second set that we need as these are for the
# artificial covariate values. Add them to the MyData object.
catchwb6_test2 <- cbind(catchwb6_test, F5.pred) 
dim(catchwb6_test2)
colnames(catchwb6_test2)
head(catchwb6_test2)

# Rename the column names so that we can use the same ggplot2 code as
# in the Poisson GLM exercise.
catchwb6_test2 <- reshape:::rename(catchwb6_test2, c("0.025quant"="lower", 
                                                     "0.975quant"="upper"))


# And the rest is a matter of plotting it all. 
# First we back-standardise LogMining. It is also an option to back-transform it.
catchwb6_test2$soaktime <- catchwb6_test2$soaktimec * sd(catchwb6_train$soaktime) + mean(catchwb6_train$soaktime)

catchwb6_pred <- catchwb6_train %>% 
  full_join(catchwb6_test2)
catchwb6_pred <- catchwb6_pred %>% 
  mutate(set = ifelse(!is.na(nettot), "Training", "Testing"))
catchwb6_pred$set <- as.factor(catchwb6_pred$set)


catchwb6_pred <- catchwb6_pred %>% 
  mutate(nettot=ifelse(!is.na(nettot), nettot, mean))

head(catchwb6_pred)

p <- ggplot()
p <- p + geom_smooth(data = catchwb6_pred %>% dplyr::filter(set=="Training"), 
                     aes(x = soaktime, 
                         y = nettot), alpha=0.5)
p <- p + geom_point(data = catchwb6_pred, 
                    aes(y = nettot, x = soaktime, col=set),
                    size = 3)
p <- p + xlab("Soak time (h)") + ylab("Total LKWF") + labs ( col = "Dataset")
p <- p + theme(text = element_text(size=15), legend.position = "bottom") 
p <- p+ xlim(18,27)
# p <- p + geom_ribbon(data = catchwb6_test3a,
#                      aes(x = soaktime2,
#                          ymax = upper,
#                          ymin = lower),
#                      alpha = 0.2)

p

p <- ggplot()
p <- p + geom_smooth(data = catchwb6_train, 
                     aes(x = soaktime, 
                         y = nettot), alpha=0.5)
p <- p + geom_point(data = catchwb6_pred, 
                    aes(y = nettot, x = soaktime, col=set),
                    size = 3)
p <- p + xlab("Soak time (h)") + ylab("Total LKWF") + labs ( col = "Dataset")
p <- p + theme(text = element_text(size=15), legend.position = "bottom") 
# p <- p + geom_ribbon(data = catchwb6_test3a,
#                      aes(x = soaktime2,
#                          ymax = upper,
#                          ymin = lower),
#                      alpha = 0.2)

p
## be care


ggplot(data=catchwb6_test2) +
  geom_point(mapping=aes(x=resp_obsN, y=mean), size=3) +
  xlab("Observed") +
  ylab("Predicted") +
  ylim(0,130) + xlim(0,130) +
  theme(text=element_text(size=15)) +
  geom_abline()

rsq <- function (x, y) cor(x, y) ^ 2
rsq(catchwb6_test2$resp_obsN, catchwb6_test2$mean)


summary(lm(resp_obsN ~ mean, data=catchwb6_test2))$r.squared
# temporal issue?
  

## sim 3 underdispersion -----

s3 <- inla(nettot ~ soaktimec +  
             f(year, model='rw1')+
             f(site, model='iid'),
           control.compute = list(config = TRUE),  ##compute TRUE
           family = "nbinomial", 
           data = catchwb6_train)

set.seed(12345)

NSim <- 1000
SimData <- inla.posterior.sample(n = NSim, result = s3) #CHANGE
nd <- length(rownames(SimData[[1]]$latent))
LastBeta <- rownames(SimData[[1]]$latent)[nd]
substrRight <- function(x, n){ substr(x, nchar(x)-n+1, nchar(x)) }
Last2Character <- substrRight(LastBeta, 2)
BetasInModel <- rownames(s3$summary.fixed) #CHANGE
MyID         <- function(x){ which(rownames(SimData[[1]]$latent) == x) }
BetasInINLA <- BetasInModel
if (Last2Character == ":1") { BetasInINLA  <- paste(BetasInModel, ":1", sep ="") }
BetaRows <- lapply(BetasInINLA, MyID)
BetaRows <- as.numeric(BetaRows)
BetaRows

# Now the random effects.
MyGrep <- function(x, SimulatedData){ 
  # Small function to determine the row numbers of random effects 
  # in the output of inla.posterior.sample
  # This is partly code from Aaron Adamack (Aaron.Adamack@dfo-mpo.gc.ca)
  # SimulatedData is the name of the object containing the simulation results
  # x is an element of BetasInINLA
  names.SimData <- attributes(SimulatedData[[1]]$latent)$dimnames[[1]]
  names.SimData[grep(x, names.SimData)] 
}

# Determine the names of the random effects
names.SimData <- attributes(SimData[[1]]$latent)$dimnames[[1]] 
YearNames     <- names.SimData[grep("year", names.SimData)]   #<-Change 'Site' for your model
YearNames

YearRows <- lapply(YearNames, MyID)
YearRows <- as.numeric(YearRows)
YearRows


SiteNames     <-  names.SimData[grep("site", names.SimData)]
SiteNames

SiteRows <- lapply(SiteNames, MyID)
SiteRows <- as.numeric(SiteRows)
SiteRows

# Start a loop to extract betas and random  effects, calculate
# the fitted values and simulate count data from 
# the model.
N    <- nrow(catchwb6_train)
Ysim <- matrix(nrow = N, ncol = NSim)
mu.i <- matrix(nrow = N, ncol = NSim)
X    <- model.matrix(~ soaktimec , data = catchwb6_train)
X   <- as.matrix(X)

YearID <- as.numeric(as.factor(catchwb6_train$year))
SiteID <- as.numeric(as.factor(catchwb6_train$site))

theta <- vector(length = NSim)        #Create space for simulated thetas



for (i in 1: NSim){
  Betas    <- SimData[[i]]$latent[BetaRows]       #Get simulated betas
  Year_i <- SimData[[i]]$latent[YearRows] #Get simulated random intercepts for year
  Site_ij  <- SimData[[i]]$latent[SiteRows]   #Get simulated random intercepts for site
  
  #eta = X * beta + year_i + Site_ij +
  eta <- X %*% Betas + Year_i[YearID] + Site_ij[SiteID] 
  
  #Get the simulated theta  
  theta[i] <- SimData[[i]]$hyperpar["size for the nbinomial observations (1/overdispersion)"]
  
  # Calculate mu = exp(eta)
  mu.i[,i] <- exp(eta)                        #Fitted values
  Ysim[,i] <- rnegbin(n = N,                  #Simulated count data
                      mu = mu.i[,i] ,
                      theta = theta[i]) 
}


# Now we have 1000 simulated data sets from the model.
# What shall we do with these simulated data sets?
# Calculate the dispersion for each simulated data set
DispNB <- vector(length = NSim) #Create space
N      <- nrow(catchwb6_train)              #Sample size 
Npar   <- length(Betas) + 2+ 1     #Number of regression parameters + theta
for(i in 1:NSim){
  ei <- (Ysim[,i] - mu.i[,i]) / sqrt(mu.i[,i] + mu.i[,i]^2 / theta[i]) #NB Pearson residuals
  DispNB[i] <- sum(ei^2) / (N - Npar)                                   #Dispersion statistic
}


# Let's plot this as a table
par(mfrow = c(1,1), mar = c(5,5,2,2), cex.lab = 1.5)
hist(DispNB, 
     xlim = c(0, 3),
     xlab = "Dispersion statistic",
     ylab = "Frequency",
     main = "Simulation results")
points(x = Dispersion14, 
       y = 0, 
       pch = 16, 
       cex = 4, 
       col = 2)


# *NB GAM ind SPDE no site------
# location exploration
Loc <- cbind(catchwb6_train$Xkm, catchwb6_train$Ykm)
D <- dist(Loc)
par(mfrow = c(1,2), mar = c(5,5,2,2), cex.lab = 1.5)
hist(D, 
     freq = TRUE,
     main = "", 
     xlab = "Distance between sites (km)",
     ylab = "Frequency")

plot(x = sort(D), 
     y = (1:length(D))/length(D), 
     type = "l",
     xlab = "Distance between sites (km)",
     ylab = "Cumulative proportion")


# Create a mesh
RangeGuess <- 15
MaxEdge  <- RangeGuess / 3.5
ConvHull <- inla.nonconvex.hull(Loc)		
mesh    <- inla.mesh.2d(loc = Loc,
                        boundary = ConvHull,
                        max.edge = c(1, 5) * MaxEdge, 
                        cutoff  = MaxEdge / 5)


par(mfrow = c(1, 1), mar=c(0, 0, 2, 0))
plot(mesh)
points(Loc, col = 1, pch = 16, cex = 1.5)

# 2. Define the weighting factors a_ik (also called 
#    the projector matrix).
A <- inla.spde.make.A(mesh, loc = Loc)

# 3. Define the SPDE.
range(catchwb6_train$nettot)
# P(Range < 5 km ) = 0.05
# P(sigma > 2) = 0.05
# SB = exp(u_i)
# some u_i have to be as large as 5.24 to cover 190
# If u_i ~ N(0, sigma_u^2) then it is unlikley that sigma_u > 1.5
#P(sigma > 1.5) = 0.05
M1 <- glm(nettot ~ 1, data = catchwb6_train)
summary(M1)


spde <- inla.spde2.pcmatern(mesh, 
                            prior.range = c(15 , 0.05), 
                            prior.sigma = c(1,  0.05))



# 4. Define the spatial field.
w.index <- inla.spde.make.index(
  name    = 'w', 
  n.spde  = spde$n.spde)



# 5. Make a stack. In this process we tell INLA 
#    at which points on the mesh we sampled the 
#    response variable and the covariates. We 
#    also need to inform INLA (via the stack) at 
#    which points of the mesh we have any other 
#    terms, e.g. the random effects as we know 
#    them from Chapter 8.


# Doing this in INLA means that we get the 95% CI for the smoother.

N    <- nrow(catchwb6_train)
Xcov <- model.matrix(~ soaktimec + year,
                     data = catchwb6_train)
Xcov2 <- Xcov[,-1]

Stack <- inla.stack(
  tag = "Fit",
  data = list(nettot = catchwb6_train$nettot),  
  A = list(A, 1, 1),                 
  effects = list(                 
    w  = w.index,  
    soaktimec = Xcov2[,1],
    year = Xcov2[,2]) ) 


f14 <- formula(catchwb6_train$nettot ~ 1 +
                 soaktimec  +
                 f (year, model='rw1') +
                 f(w, model = spde))


# 7. Run the spatial model in INLA.
I14 <- inla(f14,
            family = "nbinomial", 
            data=inla.stack.data(Stack),
            control.compute = list(dic = TRUE, waic = TRUE, config = TRUE),
            control.family=list(link='log'),
            control.predictor = list(A = inla.stack.A(Stack)))

summary(I14)
-sum(log(I14$cpo$cpo))

# Posterior mean values and 95% CI for the regression parameters:
Beta14 <- I14$summary.fixed[,c("mean", "0.025quant", "0.975quant")] 
print(Beta14, digits = 2)  

# The hyperparameters:
summary(I14)
HyperPar.I14 <- bri.hyperpar.summary(I14)
round(HyperPar.I14, digits = 3)


#Here is a visualisation of the spatial correlation.
# Access the posterior mean values of the w.
w.pm <- I14$summary.random$w$mean  
w.sd <- I14$summary.random$w$sd  

# This is a vector of length mesh$n by 1.We can use INLA functions to 
# plot these w.pm.
PlotField2(field = w.pm, 
           mesh = mesh, 
           xlim = range(mesh$loc[,1]), 
           ylim = range(mesh$loc[,2]),
           MyTitle = "Spatial random field")

# Add the sampling locations (in UTM)
points(x = Loc[,1],
       y = Loc[,2], 
       cex = 0.5, 
       col = "black", 
       pch = 16)
# conclusion: 1 area of higher (dark red) expected nettot values and 12 area of lower (dark blue)
# these effects are labeled as spatial correlation as multple neighboting sites show this effect (blue)

## either this is real correlation or I missed a covariate


# how dar for the correlation go?
SpFi.w <- inla.spde2.result(inla=I14, name="w", spde=spde, do.transfer=T) # CHANGE
Kappa <- inla.emarginal(function (x) x, SpFi.w$marginals.kappa[[1]])
Sigma_u <- inla.emarginal(function(x) sqrt(x), SpFi.w$marginals.variance.nominal[[1]])
Range <- inla.emarginal(function(x) x, SpFi.w$marginals.range.nominal[[1]]) 

Sigma_u  #Sigma of the spatial random field.
Range    #Distance at which the correlation diminishes.
LocMesh <- mesh$loc[,1:2]

# Calculate the distance between each vertex.
D <- as.matrix(dist(LocMesh))

# Using the estimated parameters from the model (see above)
# we can calculate the imposed Matern correlation values.
d.vec <- seq(0, max(D), length = 100)      
Cor.M <- (Kappa * d.vec) * besselK(Kappa * d.vec, 1) 
Cor.M[1] <- 1

# Plot the results
par(mfrow = c(1,1), mar = c(5,5,2,2))
plot(x = d.vec, 
     y = Cor.M, 
     pch = 16, 
     type = "l", 
     cex.lab = 1.5,
     xlab = "Distance", 
     ylab = "Correlation",
     xlim = c(0, 80))
abline(h=c(0.8, 0.6, 0.2), lty = 2)
abline(v=c(6, 10.75, 26), col='red')

abline(h=0.6, v=10.75)


# How many sites are separated by 0 - 10.75km?
Loc <- cbind(catchwb6_train$Xkm, catchwb6_train$Ykm)
D   <- dist(Loc)
par(mfrow = c(1,1), mar = c(5,5,2,2), cex.lab = 1.5)
hist(D, 
     freq = TRUE,
     main = "", 
     xlab = "Distance between sites (km)",
     ylab = "Frequency")
abline(v = c(10.75), col = 2, lty = 2, lwd = 2)
# no many sites are affected by this correlation!

# posterior mean of theta
theta.pd <- I14$marginals.hyperpar$`size for the nbinomial observations (1/overdispersion)`
theta.pm <- inla.emarginal(function(x) x, theta.pd)
theta.pm 

# posterior distribution of theta
par(mfrow = c(1,1), mar = c(5,5,2,2), cex.lab = 1.5)
plot(theta.pd, 
     type = "l",
     xlim = c(1,7),
     xlab = expression(paste(theta)),
     ylab = expression(paste("Pr(", theta," | data)")))

# Here are the fitted values and Pearson residuals
N <- nrow(catchwb6_train)
mu14 <- I14$summary.fitted.values[1:N,"mean"]           #Fitted values
E14  <- (catchwb6_train$nettot - mu14) /sqrt(mu14 + mu14^2 / theta.pm) #Pearson residuals

p <- nrow(I14$summary.fixed[, c("mean", "0.025quant", "0.975quant")])
Dispersion14 <- sum(E14^2) / (N - p)
Dispersion14
# 0.87


# why underdispersed
# fitted values versus the Pearson residuals
par(mfrow=c(2,1))
plot(x = mu14,
     y = E14,
     xlab = "Fitted values",
     ylab = "Pearson residuals",
     cex.lab = 1.5)
abline(h = 0, lty = 2)

par(mfrow=c(2,1))
plot(x = mu14,
     y = catchwb6_train$nettot,
     xlab = "Fitted values",
     ylab = "observed",
     cex.lab = 1.5)
abline(h = 0, lty = 2)

# Check the Pearson residuals of I2 for any remaining spatial dependency.
# Make a sample-variogram of the Pearson residuals with distances up to 0.5 km.


MyData14 <- data.frame(E14 = E14, 
                       Xkm = catchwb6_train$Xkm, 
                       Ykm = catchwb6_train$Ykm)
coordinates(MyData14) <- c("Xkm", "Ykm")
V14 <- variogram(E14 ~ 1, 
                 MyData14, 
                 cressie = TRUE, 
                 cutoff = 70)

# Plot the variogram
p <- ggplot()
p <- p + geom_point(data = V14,
                    aes(x = dist, 
                        y = gamma))
p <- p + geom_smooth(data = V13, se = FALSE, span = 0.9, aes(x = dist, y = gamma)) 
p <- p + xlab("Distance") + ylab("Semi-variogram")
p <- p + theme(text = element_text(size = 15)) 
p <- p + theme(legend.position="none") 
p

#residuals vs each covariate (in and not in the model).
# Check for non-linear patterns in the residuals.
MyVar <- c("turb_FNU", "sitedep_m", "temp_C", "DO_perc", "pH", "soaktime",  
           "distance_to_south", "distance_to_shore", "year")
catchwb6_train$E14 <- E14
MyMultipanel.ggp2(Z = catchwb6_train, 
                  varx = MyVar, 
                  vary = "E14", 
                  ylab = "Residuals",
                  addSmoother = TRUE,
                  addRegressionLine = FALSE,
                  addHorizontalLine = TRUE)
# non linearity present in all covariates, aside from temp d2S and soak

# spatial dependency?
catchwb6_train$Sign <- ifelse(E14 >= 0, "Positive residual", "Negative residual") 
ggplot() +
  theme(text = element_text(size=13)) +
  geom_point(aes(x = Xkm, y = Ykm, col = Sign), shape = 1, size = 4, data = catchwb6_train) +
  xlab("Xkm") +
  ylab("Ykm")  
# clumpy - suggest spatial dependency

boxplot(E14 ~ site, 
        ylab = "Pearson residuals",
        data = catchwb6_train)
abline(h = 0, lty = 2)


### obs vs pred 2a ---------
Loc <- cbind(catchwb6_train$Xkm, catchwb6_train$Ykm)
D <- dist(Loc)

# Create a mesh
RangeGuess <- 15
MaxEdge  <- RangeGuess / 3.5
ConvHull <- inla.nonconvex.hull(Loc)		
mesh    <- inla.mesh.2d(loc = Loc,
                        boundary = ConvHull,
                        max.edge = c(1, 5) * MaxEdge, 
                        cutoff  = MaxEdge / 5)


par(mfrow = c(1, 1), mar=c(0, 0, 2, 0))
plot(mesh)
points(Loc, col = 1, pch = 16, cex = 1.5)

# 2. Define the weighting factors a_ik (also called 
#    the projector matrix).
A <- inla.spde.make.A(mesh, loc = Loc)
dim(A)
# 3. Define the SPDE.



spde <- inla.spde2.pcmatern(mesh, 
                            prior.range = c(15 , 0.05), 
                            prior.sigma = c(1,  0.05))



# 4. Define the spatial field.
w.index <- inla.spde.make.index(
  name    = 'w', 
  n.spde  = spde$n.spde)


Xcov <- model.matrix(~ soaktimec + year,
                     data = catchwb6_train)
Xcov2 <- as.data.frame(Xcov)
colnames(Xcov2)

N    <- nrow(catchwb6_train)

Stack <- inla.stack(
  tag = "Fit",
  data = list(y = catchwb6_train$nettot),
  A = list(1,1,A),
  effects = list(
    Intercept = rep(1,N),
    Xcov2 = Xcov2[,-1], 
    w  = w.index))


f2 <- y ~ -1 + Intercept +
  soaktimec  +
  f (year, model='rw1') +
  f(w, model = spde)


I2 <- inla(f2,
           family = "nbinomial",
           data=inla.stack.data(Stack),
           control.compute = list(dic = TRUE, waic = TRUE, config = TRUE),
           control.family=list(link='log'),
           control.predictor = list(A = inla.stack.A(Stack)))

summary(I2)
# 5. Make a stack. In this process we tell INLA 
#    at which points on the mesh we sampled the 
#    response variable and the covariates. We 
#    also need to inform INLA (via the stack) at 
#    which points of the mesh we have any other 
#    terms, e.g. the random effects as we know 
#    them from Chapter 8.

# 1. Create a grid of covariate values for which we want to predict richness.
# MyData <- ddply(catchwb6_test, 
#                 .(fYear), 
#                 summarize,
#                 LogMining.std = seq(min(LogMining.std), 
#                                     max(LogMining.std), 
#                                     length = 50))
# head(MyData)  #This object has 10*50 rows and 2 columns

# equivalent is using my test data




# 2. Make a design matrix Xp

# 2. Make a design matrix Xp
Xp <- model.matrix(~ soaktimec + year, data = catchwb6_test)

# Drop the first column (intercept) from this object, and make sure
# it is a data frame.
Xp2 <- as.data.frame(Xp[,-1])
head(Xp2)



# From here onward, things diverge a little bit from the models
# without spatial correlation. This is due to the 'stack'.

#  1. Make a stack for the data for which predictions are needed.
#  2. Combine this new stack with the stack for the observed data.
#  3. Run INLA with the combined stack.
#  4. Extract the relevant pieces and plot them.


# Step 1. Make a stack for the data in MyData. Recall that MyData was 
#         converted into Xp.
StackPred <- inla.stack(
  tag = "Predictions",
  data = list(y = NA),  
  A = list(1, 1), # a 1 for each effect                 
  effects = list(
    Intercept = rep(1, nrow(Xp2)),
    Xp2 = Xp2))   #Covariate values without the intercept

# This stack does not have a spatial field, because we
# only want to have the predictions for the covariates and the 
# intercept. Note that the response variable equals NA. As a result INLA
# will predict these values.


# Step 2: Combine the two stacks.         
Combined.stacks <- inla.stack(Stack, StackPred)	              


# Step 3: Rerun INLA
# This is still the same:
f5 <- y ~ -1 + Intercept +
  soaktimec  +
  f (year, model='rw1') + 
  f(w, model = spde)

# Note that we use a different stack! The rest is the same
I5.pred <- inla(f5,
                family = "nbinomial", 
                data = inla.stack.data(Combined.stacks),
                control.compute = list(dic = TRUE, waic=TRUE, config=TRUE),
                control.family=list(link='log'),
                control.predictor = list(
                  link = 1,
                  A = inla.stack.A(Combined.stacks)))


# link = 1: Ensures that for the NAs the correct 
#           link function is used.
# Posterior mean values and 95% CI for the regression parameters:
Beta5.pred<- I5.pred$summary.fixed[,c("mean", "0.025quant", "0.975quant")] 
print(Beta5.pred, digits = 2)

summary(I5.pred)
HyperPar.I5.pred <- bri.hyperpar.summary(I5.pred)
round(HyperPar.I5.pred, digits = 3)

# Step 4: Extract the relevant pieces and plot the results.
# This is a crucial piece of crucial for extracting the correct rows.
# It provides an index for which rows in the combined stack
# belong to the observed data and which rows to the 
# artificial covariate data. 
index.Fit <- inla.stack.index(Combined.stacks,
                              tag = "Fit")$data

index.Pred <- inla.stack.index(Combined.stacks,
                               tag = "Predictions")$data

#And we can extact the correct rows     
F5.fit  <- I5.pred$summary.fitted.values[index.Fit, c(1,3,5)]  
F5.pred <- I5.pred$summary.fitted.values[index.Pred, c(1,3,5)]  


# It is the second set that we need as these are for the
# artificial covariate values. Add them to the MyData object.
catchwb6_test2 <- cbind(catchwb6_test, F5.pred) 
dim(catchwb6_test2)
colnames(catchwb6_test2)
head(catchwb6_test2)

# Rename the column names so that we can use the same ggplot2 code as
# in the Poisson GLM exercise.
catchwb6_test2 <- reshape:::rename(catchwb6_test2, c("0.025quant"="lower", 
                                                     "0.975quant"="upper", 
                                                     "mean" = "mu"))


# And the rest is a matter of plotting it all. 
# First we back-standardise LogMining. It is also an option to back-transform it.
catchwb6_test2$soaktime <- catchwb6_test2$soaktimec * sd(catchwb6_train$soaktime) + mean(catchwb6_train$soaktime)

p <- ggplot()
p <- p + geom_point(data = catchwb6_train, 
                    aes(y = nettot, x = soaktime),
                    size = 1, colour="black")
p <- p + xlab("soak time (h)") + ylab("Total LKWF in net")
p <- p + theme(text = element_text(size=15)) 
p <- p + geom_point(data = catchwb6_test2, 
                    aes(x = soaktime, 
                        y = mu), 
                    colour = "red")
p <- p + geom_ribbon(data = catchwb6_test2, 
                     aes(x = soaktime, 
                         ymax = upper, 
                         ymin = lower),
                     alpha = 0.2)
p

ggplot(data=catchwb6_test2) +
  geom_point(mapping=aes(x=resp_obsN, y=mu, colour=as.factor(catchwb6_test2$year)), size=3) +
  xlab("Observed") +
  ylab("Predicted") +
  geom_abline() +
  xlim(c(0,120)) +
  ylim(c(0,120)) 

## sim 3 underdispersion -----

s3 <- inla(nettot ~ soaktimec +  
             f(year, model='rw1')+
             f(site, model='iid'),
           control.compute = list(config = TRUE),  ##compute TRUE
           family = "nbinomial", 
           data = catchwb6_train)

set.seed(12345)

NSim <- 1000
SimData <- inla.posterior.sample(n = NSim, result = s3) #CHANGE
nd <- length(rownames(SimData[[1]]$latent))
LastBeta <- rownames(SimData[[1]]$latent)[nd]
substrRight <- function(x, n){ substr(x, nchar(x)-n+1, nchar(x)) }
Last2Character <- substrRight(LastBeta, 2)
BetasInModel <- rownames(s3$summary.fixed) #CHANGE
MyID         <- function(x){ which(rownames(SimData[[1]]$latent) == x) }
BetasInINLA <- BetasInModel
if (Last2Character == ":1") { BetasInINLA  <- paste(BetasInModel, ":1", sep ="") }
BetaRows <- lapply(BetasInINLA, MyID)
BetaRows <- as.numeric(BetaRows)
BetaRows

# Now the random effects.
MyGrep <- function(x, SimulatedData){ 
  # Small function to determine the row numbers of random effects 
  # in the output of inla.posterior.sample
  # This is partly code from Aaron Adamack (Aaron.Adamack@dfo-mpo.gc.ca)
  # SimulatedData is the name of the object containing the simulation results
  # x is an element of BetasInINLA
  names.SimData <- attributes(SimulatedData[[1]]$latent)$dimnames[[1]]
  names.SimData[grep(x, names.SimData)] 
}

# Determine the names of the random effects
names.SimData <- attributes(SimData[[1]]$latent)$dimnames[[1]] 
YearNames     <- names.SimData[grep("year", names.SimData)]   #<-Change 'Site' for your model
YearNames

YearRows <- lapply(YearNames, MyID)
YearRows <- as.numeric(YearRows)
YearRows


SiteNames     <-  names.SimData[grep("site", names.SimData)]
SiteNames

SiteRows <- lapply(SiteNames, MyID)
SiteRows <- as.numeric(SiteRows)
SiteRows

# Start a loop to extract betas and random  effects, calculate
# the fitted values and simulate count data from 
# the model.
N    <- nrow(catchwb6_train)
Ysim <- matrix(nrow = N, ncol = NSim)
mu.i <- matrix(nrow = N, ncol = NSim)
X    <- model.matrix(~ soaktimec , data = catchwb6_train)
X   <- as.matrix(X)

YearID <- as.numeric(as.factor(catchwb6_train$year))
SiteID <- as.numeric(as.factor(catchwb6_train$site))

theta <- vector(length = NSim)        #Create space for simulated thetas



for (i in 1: NSim){
  Betas    <- SimData[[i]]$latent[BetaRows]       #Get simulated betas
  Year_i <- SimData[[i]]$latent[YearRows] #Get simulated random intercepts for year
  Site_ij  <- SimData[[i]]$latent[SiteRows]   #Get simulated random intercepts for site
  
  #eta = X * beta + year_i + Site_ij +
  eta <- X %*% Betas + Year_i[YearID] + Site_ij[SiteID] 
  
  #Get the simulated theta  
  theta[i] <- SimData[[i]]$hyperpar["size for the nbinomial observations (1/overdispersion)"]
  
  # Calculate mu = exp(eta)
  mu.i[,i] <- exp(eta)                        #Fitted values
  Ysim[,i] <- rnegbin(n = N,                  #Simulated count data
                      mu = mu.i[,i] ,
                      theta = theta[i]) 
}


# Now we have 1000 simulated data sets from the model.
# What shall we do with these simulated data sets?
# Calculate the dispersion for each simulated data set
DispNB <- vector(length = NSim) #Create space
N      <- nrow(catchwb6_train)              #Sample size 
Npar   <- length(Betas) + 2+ 1     #Number of regression parameters + theta
for(i in 1:NSim){
  ei <- (Ysim[,i] - mu.i[,i]) / sqrt(mu.i[,i] + mu.i[,i]^2 / theta[i]) #NB Pearson residuals
  DispNB[i] <- sum(ei^2) / (N - Npar)                                   #Dispersion statistic
}


# Let's plot this as a table
par(mfrow = c(1,1), mar = c(5,5,2,2), cex.lab = 1.5)
hist(DispNB, 
     xlim = c(0, 3),
     xlab = "Dispersion statistic",
     ylab = "Frequency",
     main = "Simulation results")
points(x = Dispersion14, 
       y = 0, 
       pch = 16, 
       cex = 4, 
       col = 2)

# NB GAM envr all ind SPDE------
# location exploration
Loc <- cbind(catchwb6_train$Xkm, catchwb6_train$Ykm)
D <- dist(Loc)
par(mfrow = c(1,2), mar = c(5,5,2,2), cex.lab = 1.5)
hist(D, 
     freq = TRUE,
     main = "", 
     xlab = "Distance between sites (km)",
     ylab = "Frequency")

plot(x = sort(D), 
     y = (1:length(D))/length(D), 
     type = "l",
     xlab = "Distance between sites (km)",
     ylab = "Cumulative proportion")


# Create a mesh
RangeGuess <- 15
MaxEdge  <- RangeGuess / 3.5
ConvHull <- inla.nonconvex.hull(Loc)		
mesh    <- inla.mesh.2d(loc = Loc,
                        boundary = ConvHull,
                        max.edge = c(1, 5) * MaxEdge, 
                        cutoff  = MaxEdge / 5)


par(mfrow = c(1, 1), mar=c(0, 0, 2, 0))
plot(mesh)
points(Loc, col = 1, pch = 16, cex = 1.5)

# 2. Define the weighting factors a_ik (also called 
#    the projector matrix).
A <- inla.spde.make.A(mesh, loc = Loc)

# 3. Define the SPDE.
range(catchwb6_train$nettot)
# P(Range < 5 km ) = 0.05
# P(sigma > 2) = 0.05
# SB = exp(u_i)
# some u_i have to be as large as 5.24 to cover 190
# If u_i ~ N(0, sigma_u^2) then it is unlikley that sigma_u > 1.5
#P(sigma > 1.5) = 0.05
M1 <- glm(nettot ~ 1, data = catchwb6_train)
summary(M1)


spde <- inla.spde2.pcmatern(mesh, 
                            prior.range = c(15 , 0.05), 
                            prior.sigma = c(1,  0.05))



# 4. Define the spatial field.
w.index <- inla.spde.make.index(
  name    = 'w', 
  n.spde  = spde$n.spde)



# 5. Make a stack. In this process we tell INLA 
#    at which points on the mesh we sampled the 
#    response variable and the covariates. We 
#    also need to inform INLA (via the stack) at 
#    which points of the mesh we have any other 
#    terms, e.g. the random effects as we know 
#    them from Chapter 8.


# Doing this in INLA means that we get the 95% CI for the smoother.

Xcov <- model.matrix(~ soaktimec +  temp_Cc + turb_FNUc + distance_to_southc + distance_to_shorec + year ,
                     data = catchwb6_train)
Xcov2 <- as.data.frame(Xcov)
colnames(Xcov2)

N    <- nrow(catchwb6_train)

Stack <- inla.stack(
  tag = "Fit",
  data = list(y = catchwb6_train$nettot),
  A = list(1,1,A),
  effects = list(
    Intercept = rep(1,N),
    Xcov2 = Xcov2[,-1], 
    w  = w.index))


f14 <- y ~ -1 + Intercept +
  soaktimec  + temp_Cc + turb_FNUc + distance_to_southc + distance_to_shorec +
  f (year, model='rw1') +
  f(w, model = spde)


I14 <- inla(f14,
           family = "nbinomial",
           data=inla.stack.data(Stack),
           control.compute = list(dic = TRUE, waic = TRUE, config = TRUE),
           control.family=list(link='log'),
           control.predictor = list(A = inla.stack.A(Stack)))

summary(I14)


# Posterior mean values and 95% CI for the regression parameters:
Beta14 <- I14$summary.fixed[,c("mean", "0.025quant", "0.975quant")] 
print(Beta14, digits = 2)  
## 0 is in all intervals so not important

# The hyperparameters:
summary(I14)
HyperPar.I14 <- bri.hyperpar.summary(I14)
round(HyperPar.I14, digits = 3)


#Here is a visualisation of the spatial correlation.
# Access the posterior mean values of the w.
w.pm <- I14$summary.random$w$mean  
w.sd <- I14$summary.random$w$sd  

# This is a vector of length mesh$n by 1.We can use INLA functions to 
# plot these w.pm.
PlotField2(field = w.pm, 
           mesh = mesh, 
           xlim = range(mesh$loc[,1]), 
           ylim = range(mesh$loc[,2]),
           MyTitle = "Spatial random field")

# Add the sampling locations (in UTM)
points(x = Loc[,1],
       y = Loc[,2], 
       cex = 0.5, 
       col = "black", 
       pch = 16)
# conclusion: 2 area of higher (dark red) expected nettot values and 1 area of lower (dark blue)
# these effects are labeled as spatial correlation as multple neighboting sites show this effect (blue)

## either this is real correlation or I missed a covariate


# how dar for the correlation go?
SpFi.w <- inla.spde2.result(inla=I14, name="w", spde=spde, do.transfer=T) # CHANGE
Kappa <- inla.emarginal(function (x) x, SpFi.w$marginals.kappa[[1]])
sigma_u <- inla.emarginal(function(x) sqrt(x), SpFi.w$marginals.variance.nominal[[1]])
Range <- inla.emarginal(function(x) x, SpFi.w$marginals.range.nominal[[1]]) 

sigma_u  #Sigma of the spatial random field.
Range    #Distance at which the correlation diminishes.
LocMesh <- mesh$loc[,1:2]

# Calculate the distance between each vertex.
D <- as.matrix(dist(LocMesh))

# Using the estimated parameters from the model (see above)
# we can calculate the imposed Matern correlation values.
d.vec <- seq(0, max(D), length = 100)      
Cor.M <- (Kappa * d.vec) * besselK(Kappa * d.vec, 1) 
Cor.M[1] <- 1

# Plot the results
par(mfrow = c(1,1), mar = c(5,5,2,2))
plot(x = d.vec, 
     y = Cor.M, 
     pch = 16, 
     type = "l", 
     cex.lab = 1.5,
     xlab = "Distance", 
     ylab = "Correlation",
     xlim = c(0, 80))
abline(h=c(0.8, 0.6, 0.2), lty = 2)
abline(v=c(6.5, 12, 28.5), col='red')

abline(h=0.6, v=12)


# How many sites are separated by 0 - 10.25km?
Loc <- cbind(catchwb6_train$Xkm, catchwb6_train$Ykm)
D   <- dist(Loc)
par(mfrow = c(1,1), mar = c(5,5,2,2), cex.lab = 1.5)
hist(D, 
     freq = TRUE,
     main = "", 
     xlab = "Distance between sites (km)",
     ylab = "Frequency")
abline(v = c(12), col = 2, lty = 2, lwd = 2)
# no many sites are affected by this correlation!

# posterior mean of theta
theta.pd <- I14$marginals.hyperpar$`size for the nbinomial observations (1/overdispersion)`
theta.pm <- inla.emarginal(function(x) x, theta.pd)
theta.pm 

# posterior distribution of theta
par(mfrow = c(1,1), mar = c(5,5,2,2), cex.lab = 1.5)
plot(theta.pd, 
     type = "l",
     xlim = c(1,7),
     xlab = expression(paste(theta)),
     ylab = expression(paste("Pr(", theta," | data)")))

# Here are the fitted values and Pearson residuals
N <- nrow(catchwb6_train)
mu14 <- I14$summary.fitted.values[1:N,"mean"]           #Fitted values
E14  <- (catchwb6_train$nettot - mu14) /sqrt(mu14 + mu14^2 / theta.pm) #Pearson residuals

p <- nrow(I14$summary.fixed[, c("mean", "0.025quant", "0.975quant")])
Dispersion14 <- sum(E14^2) / (N - p)
Dispersion14
# 0.87


# why underdispersed
# fitted values versus the Pearson residuals
par(mfrow=c(2,1))
plot(x = mu14,
     y = E14,
     xlab = "Fitted values",
     ylab = "Pearson residuals",
     cex.lab = 1.5)
abline(h = 0, lty = 2)

par(mfrow=c(2,1))
plot(x = mu14,
     y = catchwb6_train$nettot,
     xlab = "Fitted values",
     ylab = "observed",
     cex.lab = 1.5)
abline(h = 0, lty = 2)

MyData14 <- data.frame(E14 = E14, 
                       Xkm = catchwb6_train$Xkm, 
                       Ykm = catchwb6_train$Ykm)
coordinates(MyData14) <- c("Xkm", "Ykm")
V14 <- variogram(E14 ~ 1, 
                 MyData14, 
                 cressie = TRUE, 
                 cutoff = 70)

# Plot the variogram
p <- ggplot()
p <- p + geom_point(data = V14,
                    aes(x = dist, 
                        y = gamma))
p <- p + geom_smooth(data = V13, se = FALSE, span = 0.9, aes(x = dist, y = gamma)) 
p <- p + xlab("Distance") + ylab("Semi-variogram")
p <- p + theme(text = element_text(size = 15)) 
p <- p + theme(legend.position="none") 
p



# Check the Pearson residuals of I2 for any remaining spatial dependency.
# Make a sample-variogram of the Pearson residuals with distances up to 0.5 km.
MyData  <- data.frame(E14 = E14, 
                      X  = catchwb6_train$Xkm, 
                      Y  = catchwb6_train$Ykm)
coordinates(MyData)    <- c("X", "Y")
V14 <- variogram(E14 ~ 1, MyData, cressie = TRUE, cutoff = 0.45)

# Plot the variogram
p <- ggplot()
p <- p + geom_point(data = V14,
                    aes(x = dist, 
                        y = gamma))
p <- p + geom_smooth(data = V14,
                     se = FALSE,
                     span = 0.9,
                     aes(x = dist, 
                         y = gamma))
p <- p + xlab("Distance") + ylab("Semi-variogram")
p <- p + theme(text = element_text(size = 15)) 
p <- p + theme(legend.position="none") 
p
boxplot(E14 ~ site, 
        ylab = "Pearson residuals",
        data = catchwb6_train)
abline(h = 0, lty = 2)

#residuals vs each covariate (in and not in the model).
# Check for non-linear patterns in the residuals.
MyVar <- c("turb_FNU", "sitedep_m", "temp_C", "DO_perc", "pH", "soaktime",  
           "distance_to_south", "distance_to_shore", "year")
catchwb6_train$E14 <- E14
MyMultipanel.ggp2(Z = catchwb6_train, 
                  varx = MyVar, 
                  vary = "E14", 
                  ylab = "Residuals",
                  addSmoother = TRUE,
                  addRegressionLine = FALSE,
                  addHorizontalLine = TRUE)
# non linearity present in all covariates, aside from temp d2S and soak

# spatial dependency?
catchwb6_train$Sign <- ifelse(E14 >= 0, "Positive residual", "Negative residual") 
ggplot() +
  theme(text = element_text(size=13)) +
  geom_point(aes(x = Xkm, y = Ykm, col = Sign), shape = 1, size = 4, data = catchwb6_train) +
  xlab("Xkm") +
  ylab("Ykm")  
# clumpy - suggest spatial dependency




### obs vs pred 2a ---------
Loc <- cbind(catchwb6_train$Xkm, catchwb6_train$Ykm)
D <- dist(Loc)

# Create a mesh
RangeGuess <- 15
MaxEdge  <- RangeGuess / 3.5
ConvHull <- inla.nonconvex.hull(Loc)		
mesh    <- inla.mesh.2d(loc = Loc,
                        boundary = ConvHull,
                        max.edge = c(1, 5) * MaxEdge, 
                        cutoff  = MaxEdge / 5)


par(mfrow = c(1, 1), mar=c(0, 0, 2, 0))
plot(mesh)
points(Loc, col = 1, pch = 16, cex = 1.5)

# 2. Define the weighting factors a_ik (also called 
#    the projector matrix).
A <- inla.spde.make.A(mesh, loc = Loc)
dim(A)
# 3. Define the SPDE.



spde <- inla.spde2.pcmatern(mesh, 
                            prior.range = c(15 , 0.05), 
                            prior.sigma = c(1,  0.05))



# 4. Define the spatial field.
w.index <- inla.spde.make.index(
  name    = 'w', 
  n.spde  = spde$n.spde)


Xcov <- model.matrix(~ soaktimec +  temp_Cc + turb_FNUc + distance_to_southc + distance_to_shorec + year ,
                     data = catchwb6_train)
Xcov2 <- as.data.frame(Xcov)
colnames(Xcov2)

N    <- nrow(catchwb6_train)

Stack <- inla.stack(
  tag = "Fit",
  data = list(y = catchwb6_train$nettot),
  A = list(1,1,A),
  effects = list(
    Intercept = rep(1,N),
    Xcov2 = Xcov2[,-1], 
    w  = w.index))


f2 <- y ~ -1 + Intercept +
  soaktimec  + temp_Cc + turb_FNUc + distance_to_southc + distance_to_shorec +
  f (year, model='rw1') +
  f(w, model = spde)


I2 <- inla(f2,
           family = "nbinomial",
           data=inla.stack.data(Stack),
           control.compute = list(dic = TRUE, waic = TRUE, config = TRUE),
           control.family=list(link='log'),
           control.predictor = list(A = inla.stack.A(Stack)))

summary(I2)
# 5. Make a stack. In this process we tell INLA 
#    at which points on the mesh we sampled the 
#    response variable and the covariates. We 
#    also need to inform INLA (via the stack) at 
#    which points of the mesh we have any other 
#    terms, e.g. the random effects as we know 
#    them from Chapter 8.

# 1. Create a grid of covariate values for which we want to predict richness.
# MyData <- ddply(catchwb6_test, 
#                 .(fYear), 
#                 summarize,
#                 LogMining.std = seq(min(LogMining.std), 
#                                     max(LogMining.std), 
#                                     length = 50))
# head(MyData)  #This object has 10*50 rows and 2 columns

# equivalent is using my test data




# 2. Make a design matrix Xp

# 2. Make a design matrix Xp
Xp <- model.matrix(~ soaktimec +  temp_Cc + turb_FNUc + distance_to_southc + distance_to_shorec + year,
                   data = catchwb6_test)

# Drop the first column (intercept) from this object, and make sure
# it is a data frame.
Xp2 <- as.data.frame(Xp[,-1])
head(Xp2)



# From here onward, things diverge a little bit from the models
# without spatial correlation. This is due to the 'stack'.

#  1. Make a stack for the data for which predictions are needed.
#  2. Combine this new stack with the stack for the observed data.
#  3. Run INLA with the combined stack.
#  4. Extract the relevant pieces and plot them.


# Step 1. Make a stack for the data in MyData. Recall that MyData was 
#         converted into Xp.
StackPred <- inla.stack(
  tag = "Predictions",
  data = list(y = NA),  
  A = list(1, 1), # a 1 for each effect                 
  effects = list(
    Intercept = rep(1, nrow(Xp2)),
    Xp2 = Xp2))   #Covariate values without the intercept

# This stack does not have a spatial field, because we
# only want to have the predictions for the covariates and the 
# intercept. Note that the response variable equals NA. As a result INLA
# will predict these values.


# Step 2: Combine the two stacks.         
Combined.stacks <- inla.stack(Stack, StackPred)	              


# Step 3: Rerun INLA
# This is still the same:
f5 <- y ~ -1 + Intercept +
  +   soaktimec  + temp_Cc + turb_FNUc + distance_to_southc + distance_to_shorec +
  +   f (year, model='rw1') +
  +   f(w, model = spde)

# Note that we use a different stack! The rest is the same
I5.pred <- inla(f5,
                family = "nbinomial", 
                data = inla.stack.data(Combined.stacks),
                control.compute = list(dic = TRUE, waic=TRUE, config=TRUE),
                control.family=list(link='log'),
                control.predictor = list(
                  link = 1,
                  A = inla.stack.A(Combined.stacks)))


# link = 1: Ensures that for the NAs the correct 
#           link function is used.

# Posterior mean values and 95% CI for the regression parameters:
Beta5.pred<- I5.pred$summary.fixed[,c("mean", "0.025quant", "0.975quant")] 
print(Beta5.pred, digits = 2)

summary(I5.pred)
HyperPar.I5.pred <- bri.hyperpar.summary(I5.pred)
round(HyperPar.I5.pred, digits = 3)

# Step 4: Extract the relevant pieces and plot the results.
# This is a crucial piece of crucial for extracting the correct rows.
# It provides an index for which rows in the combined stack
# belong to the observed data and which rows to the 
# artificial covariate data. 
index.Fit <- inla.stack.index(Combined.stacks,
                              tag = "Fit")$data

index.Pred <- inla.stack.index(Combined.stacks,
                               tag = "Predictions")$data

#And we can extact the correct rows     
F5.fit  <- I5.pred$summary.fitted.values[index.Fit, c(1,3,5)]  
F5.pred <- I5.pred$summary.fitted.values[index.Pred, c(1,3,5)]  


# It is the second set that we need as these are for the
# artificial covariate values. Add them to the MyData object.
catchwb6_test2 <- cbind(catchwb6_test, F5.pred) 
dim(catchwb6_test2)
colnames(catchwb6_test2)
head(catchwb6_test2)

# Rename the column names so that we can use the same ggplot2 code as
# in the Poisson GLM exercise.
catchwb6_test2 <- reshape:::rename(catchwb6_test2, c("0.025quant"="lower", 
                                                     "0.975quant"="upper", 
                                                     "mean" = "mu"))


# And the rest is a matter of plotting it all. 
# First we back-standardise LogMining. It is also an option to back-transform it.
catchwb6_test2$soaktime <- catchwb6_test2$soaktimec * sd(catchwb6_train$soaktime) + mean(catchwb6_train$soaktime)

p <- ggplot()
p <- p + geom_point(data = catchwb6_train, 
                    aes(y = nettot, x = soaktime),
                    size = 1, colour="black")
p <- p + xlab("soak time (h)") + ylab("Total LKWF in net")
p <- p + theme(text = element_text(size=15)) 
p <- p + geom_point(data = catchwb6_test2, 
                    aes(x = soaktime, 
                        y = mu), 
                    colour = "red")
p <- p + geom_ribbon(data = catchwb6_test2, 
                     aes(x = soaktime, 
                         ymax = upper, 
                         ymin = lower),
                     alpha = 0.2)
p

ggplot(data=catchwb6_test2) +
  geom_point(mapping=aes(x=resp_obsN, y=mu, colour=as.factor(catchwb6_test2$year)), size=3) +
  xlab("Observed") +
  ylab("Predicted") +
  geom_abline() +
  xlim(c(0,120)) +
  ylim(c(0,120)) 

## sim 3 underdispersion -----

s3 <- inla(nettot ~ soaktimec +  temp_Cc +
             f(year, model='rw1')+
             f(site, model='iid'),
           control.compute = list(config = TRUE),  ##compute TRUE
           family = "nbinomial", 
           data = catchwb6_train)

set.seed(12345)

NSim <- 1000
SimData <- inla.posterior.sample(n = NSim, result = s3) #CHANGE
nd <- length(rownames(SimData[[1]]$latent))
LastBeta <- rownames(SimData[[1]]$latent)[nd]
substrRight <- function(x, n){ substr(x, nchar(x)-n+1, nchar(x)) }
Last2Character <- substrRight(LastBeta, 2)
BetasInModel <- rownames(s3$summary.fixed) #CHANGE
MyID         <- function(x){ which(rownames(SimData[[1]]$latent) == x) }
BetasInINLA <- BetasInModel
if (Last2Character == ":1") { BetasInINLA  <- paste(BetasInModel, ":1", sep ="") }
BetaRows <- lapply(BetasInINLA, MyID)
BetaRows <- as.numeric(BetaRows)
BetaRows

# Now the random effects.
MyGrep <- function(x, SimulatedData){ 
  # Small function to determine the row numbers of random effects 
  # in the output of inla.posterior.sample
  # This is partly code from Aaron Adamack (Aaron.Adamack@dfo-mpo.gc.ca)
  # SimulatedData is the name of the object containing the simulation results
  # x is an element of BetasInINLA
  names.SimData <- attributes(SimulatedData[[1]]$latent)$dimnames[[1]]
  names.SimData[grep(x, names.SimData)] 
}

# Determine the names of the random effects
names.SimData <- attributes(SimData[[1]]$latent)$dimnames[[1]] 
YearNames     <- names.SimData[grep("year", names.SimData)]   #<-Change 'Site' for your model
YearNames

YearRows <- lapply(YearNames, MyID)
YearRows <- as.numeric(YearRows)
YearRows


SiteNames     <-  names.SimData[grep("site", names.SimData)]
SiteNames

SiteRows <- lapply(SiteNames, MyID)
SiteRows <- as.numeric(SiteRows)
SiteRows

# Start a loop to extract betas and random  effects, calculate
# the fitted values and simulate count data from 
# the model.
N    <- nrow(catchwb6_train)
Ysim <- matrix(nrow = N, ncol = NSim)
mu.i <- matrix(nrow = N, ncol = NSim)
X    <- model.matrix(~ soaktimec + temp_Cc, data = catchwb6_train)
X   <- as.matrix(X)

YearID <- as.numeric(as.factor(catchwb6_train$year))
SiteID <- as.numeric(as.factor(catchwb6_train$site))

theta <- vector(length = NSim)        #Create space for simulated thetas



for (i in 1: NSim){
  Betas    <- SimData[[i]]$latent[BetaRows]       #Get simulated betas
  Year_i <- SimData[[i]]$latent[YearRows] #Get simulated random intercepts for year
  Site_ij  <- SimData[[i]]$latent[SiteRows]   #Get simulated random intercepts for site
  
  #eta = X * beta + year_i + Site_ij +
  eta <- X %*% Betas + Year_i[YearID] + Site_ij[SiteID] 
  
  #Get the simulated theta  
  theta[i] <- SimData[[i]]$hyperpar["size for the nbinomial observations (1/overdispersion)"]
  
  # Calculate mu = exp(eta)
  mu.i[,i] <- exp(eta)                        #Fitted values
  Ysim[,i] <- rnegbin(n = N,                  #Simulated count data
                      mu = mu.i[,i] ,
                      theta = theta[i]) 
}


# Now we have 1000 simulated data sets from the model.
# What shall we do with these simulated data sets?
# Calculate the dispersion for each simulated data set
DispNB <- vector(length = NSim) #Create space
N      <- nrow(catchwb6_train)              #Sample size 
Npar   <- length(Betas) + 2+ 1     #Number of regression parameters + theta
for(i in 1:NSim){
  ei <- (Ysim[,i] - mu.i[,i]) / sqrt(mu.i[,i] + mu.i[,i]^2 / theta[i]) #NB Pearson residuals
  DispNB[i] <- sum(ei^2) / (N - Npar)                                   #Dispersion statistic
}


# Let's plot this as a table
par(mfrow = c(1,1), mar = c(5,5,2,2), cex.lab = 1.5)
hist(DispNB, 
     xlim = c(0, 3),
     xlab = "Dispersion statistic",
     ylab = "Frequency",
     main = "Simulation results")
points(x = Dispersion.nb, 
       y = 0, 
       pch = 16, 
       cex = 4, 
       col = 2)



# Catch rate ------
## basic CPUE calc 
a <- ggplot(mapping=aes(x=year, y=nettot/soaktime), data=catchwb6_pred) +
  geom_point() +
  geom_smooth() +
  xlab("Year") +
  ylab("CPUE")
  
## CPUE from the selected model (SRF )
b <- ggplot(mapping=aes(x=year, y=mu14/soaktime), data=catchwb6_train) +
  geom_point() +
  geom_smooth() +
  ylim(0,6) +
  xlab("Year") +
  ylab("CPUE")

a/b
