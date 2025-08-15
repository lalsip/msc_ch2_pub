library(lattice)  
library(ggplot2)
library(MASS)
library(mgcv)
library(plyr)
library(tidyr)
library(INLA)

library(glmmTMB)
library(GGally)
# library(brinla)
# install.packages('brinla')

library(ggmap)
library(sp)
library(gstat)
library(fields)
library(RColorBrewer)
source('C:/Users/alsipl/Documents/MSc/src/SPDEsupport.R')

inla.setOption(num.threads=8) # sets the number of cores being used by INLA

catchwb6_pred <- catchwb6_pred %>% 
  mutate(set = ifelse(!is.na(nettot), "Training", "Testing"))

ggplot(data=catchwb6_pred) +
  geom_point(mapping=aes(x=soaktime, y=avgwt_g), colour="black") +
  geom_point(mapping=aes(x=soaktime, y=resp_obsB), colour="red") +
  geom_smooth(mapping=aes(x=soaktime, y=resp_obsB), colour="red") +
  xlab ("Soak time (h)") +
  ylab ("Mean weight (g) of LKWF") +
  xlim(18, 27)

# Base Density models ------
# Effort: soak, net area, year
# *Gamma GLM no site------
I1 <- inla(avgwt_g ~ year + soaktimec ,
           control.compute = list(dic = TRUE, waic=TRUE, config=TRUE),
           control.family = list(link='log'),
           family = "gamma",
           data = catchwb6_train)
summary(I1)


# Posterior mean values and 95% CI for the regression parameters:
Beta1 <- I1$summary.fixed[,c("mean", "sd", "0.025quant", "0.975quant")] 
print(Beta1, digits = 2)  

# Hyperparameters
summary(I1)
HyperPar.I1 <- bri.hyperpar.summary(I1)
round(HyperPar.I1, digits = 3)

# This is the parameter r from the Gamma variance
# dispersion stat - denominator in varaiance of gamma distribution
r  <- I1$summary.hyperpar["Precision parameter for the Gamma observations", "mean"]
r



# Get the fitted values and Pearson residuals
N     <- nrow(catchwb6_train)
mu1   <- I1$summary.fitted.values[1:N,"mean"] 
VarY1 <- mu1^2 / r
E1    <- (catchwb6_train$avgwt_g - mu1) / sqrt(VarY1)

# A Gamma GLM(M) cannot be overdispersed


range(catchwb6_train$avgwt_g) # 65.5-1256.923

# why overdispersed
# fitted values versus the Pearson residuals
par(mfrow=c(1,1))
plot(x = mu1,
     y = E1,
     xlab = "Fitted values",
     ylab = "Pearson residuals",
     cex.lab = 1.5)
abline(h = 0, lty = 2)
# fitted ranges 300 to 700 (whole range should be 65 to 1260)
# residuals evenly spread
# good range of residuals (-2 ro 2)

# Plot residuals versus site
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

I1a <- inla(avgwt_g ~ year + soaktimec +
              f(site, model='iid'),
            control.compute = list(dic = TRUE, cpo=TRUE, waic=TRUE, config=TRUE),
            family = "gamma",
            control.family=list(link='log'),
            data = catchwb6_pred) # now pred
summary(I1a) # should be the same as before

#extract the predictor values associated with the appended data
catchwb6_test1a <- cbind(catchwb6_test, 
                         I1a$summary.linear.predictor[(nrow(catchwb6_train)+1):nrow(catchwb6_pred),])
head(catchwb6_test1a)

catchwb6_test1a <- reshape:::rename(catchwb6_test1a, c("0.025quant"="lower", "0.975quant"="upper"))

ggplot(catchwb6_test1a, aes(y=mean, x=soaktimec)) +
  geom_point(data=catchwb6_train, aes(y=avgwt_g)) + # obs
  geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.2) + # predicted SD
  geom_line() # predicted
# MUCH smaller than observed

ggplot(catchwb6_test1a, aes(y=mean, x=soaktimec)) +
  geom_point(aes(y=resp_obsB)) + # obs
  geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.2) + # predicted SD
  geom_line()

# And the rest is a matter of plotting it all. 
# First we back-standardise LogMining. It is also an option to back-transform it.
catchwb6_test1a$soaktime2 <- catchwb6_test1a$soaktimec * sd(catchwb6_train$soaktime) + mean(catchwb6_train$soaktime)
p <- ggplot()
p <- p + geom_point(data = catchwb6_train, 
                    aes(y = avgwt_g, x = soaktime),
                    size = 1)
p <- p + xlab("soak time (h)") + ylab("total number LKWF")
p <- p + theme(text = element_text(size=15)) 
p <- p + geom_point(data = catchwb6_test1a, 
                    aes(x = soaktime2, 
                        y = mean), 
                    col = "red")
p <- p + geom_ribbon(data = catchwb6_test1a, 
                     aes(x = soaktime2, 
                         ymax = upper, 
                         ymin = lower),
                     alpha = 0.2)
p

ggplot(data=catchwb6_test1a) +
  geom_point(mapping=aes(x=resp_obsB, y=mean)) +
  xlab("Observed") +
  ylab("Predicted") 


# Gamma GLM ------
I1 <- inla(avgwt_g ~ year + soaktimec +
             f(site, model='iid'),
           control.compute = list(dic = TRUE, cpo=TRUE, waic=TRUE, config=TRUE),
           control.family = list(link='log'),
           family = "gamma",
           data = catchwb6_train)
summary(I1)
-sum(log(I1$cpo$cpo)) # smaller is better fit of model to data

# Posterior mean values and 95% CI for the regression parameters:
Beta1 <- I1$summary.fixed[,c("mean", "sd", "0.025quant", "0.975quant")] 
print(Beta1, digits = 2)  

# Hyperparameters
summary(I1)
HyperPar.I1 <- bri.hyperpar.summary(I1)
round(HyperPar.I1, digits = 3)

# This is the parameter r from the Gamma variance
# dispersion stat - denominator in varaiance of gamma distribution
r  <- I1$summary.hyperpar["Precision parameter for the Gamma observations", "mean"]
r



# Get the fitted values and Pearson residuals
N     <- nrow(catchwb6_train)
mu1   <- I1$summary.fitted.values[1:N,"mean"] 
VarY1 <- mu1^2 / r
E1    <- (catchwb6_train$avgwt_g - mu1) / sqrt(VarY1)

# A Gamma GLM(M) cannot be overdispersed


range(catchwb6_train$avgwt_g) # 65.5-1256.923

# why overdispersed
# fitted values versus the Pearson residuals
par(mfrow=c(1,1))
plot(x = mu1,
     y = E1,
     xlab = "Fitted values",
     ylab = "Pearson residuals",
     cex.lab = 1.5)
abline(h = 0, lty = 2)
# fitted ranges 300 to 700 (whole range should be 65 to 1260)
# residuals evenly spread
# good range of residuals (-2 ro 2)

# Plot residuals versus site
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

# no

# par(mfrow=c(1,1))
# acf(catchwb6_train$nettot, catchwb6_train$year)

### obs vs pred 1a ---------

I1a <- inla(avgwt_g ~ year + soaktimec +
              f(site, model='iid'),
            control.compute = list(dic = TRUE, cpo=TRUE, waic=TRUE, config=TRUE),
            family = "gamma",
            control.family=list(link='log'),
            data = catchwb6_pred) # now pred
summary(I1a) # should be the same as before

#extract the predictor values associated with the appended data
catchwb6_test1a <- cbind(catchwb6_test, 
                         I1a$summary.linear.predictor[(nrow(catchwb6_train)+1):nrow(catchwb6_pred),])
head(catchwb6_test1a)

catchwb6_test1a <- reshape:::rename(catchwb6_test1a, c("0.025quant"="lower", "0.975quant"="upper"))

catchwb6_test1a$soaktime2 <- catchwb6_test1a$soaktimec * sd(catchwb6_train$soaktime) + mean(catchwb6_train$soaktime)
ggplot(catchwb6_test1a, aes(y=mean, x=soaktime)) +
  geom_point(data=catchwb6_train, aes(y=avgwt_g)) + # obs
  geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.2) + # predicted SD
  geom_line() # predicted
# MUCH smaller than observed

ggplot(catchwb6_test1a, aes(y=mean, x=soaktime2)) +
  geom_point(aes(y=resp_obsB)) + # obs
  geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.2) + # predicted SD
  geom_line()

p <- ggplot()
p <- p + geom_point(data = catchwb6_train, 
                    aes(y = avgwt_g, x = soaktime),
                    size = 1)
p <- p + xlab("soak time (h)") + ylab("avg wt LKWF")
p <- p + theme(text = element_text(size=15)) 
p <- p + geom_point(data = catchwb6_test1a, 
                    aes(x = soaktime, 
                        y = mean), 
                    col = "red")
p <- p + geom_ribbon(data = catchwb6_test1a, 
                     aes(x = soaktime, 
                         ymax = upper, 
                         ymin = lower),
                     alpha = 0.2)
p

# *Gamma GAM ------
# not used
I1 <- inla(avgwt_g ~ soaktimec +
             f(year, model='rw1') +
             f(site, model='iid'),
           control.compute = list(dic = TRUE, cpo=TRUE, waic=TRUE, config=TRUE),
           control.family = list(link='log'),
           family = "gamma",
           data = catchwb6_train)
summary(I1)
-sum(log(I1$cpo$cpo)) # smaller is better fit of model to data

# Posterior mean values and 95% CI for the regression parameters:
Beta1 <- I1$summary.fixed[,c("mean", "sd", "0.025quant", "0.975quant")] 
print(Beta1, digits = 2)  

# Hyperparameters
summary(I1)
HyperPar.I1 <- bri.hyperpar.summary(I1)
round(HyperPar.I1, digits = 3)

# This is the parameter r from the Gamma variance
# dispersion stat - denominator in varaiance of gamma distribution
r  <- I1$summary.hyperpar["Precision parameter for the Gamma observations", "mean"]
r



# Get the fitted values and Pearson residuals
N     <- nrow(catchwb6_train)
mu1   <- I1$summary.fitted.values[1:N,"mean"] 
VarY1 <- mu1^2 / r
E1    <- (catchwb6_train$avgwt_g - mu1) / sqrt(VarY1)

# A Gamma GLM(M) cannot be overdispersed


range(catchwb6_train$avgwt_g) # 65.5-1256.923

# why overdispersed
# fitted values versus the Pearson residuals
par(mfrow=c(1,1))
plot(x = mu1,
     y = E1,
     xlab = "Fitted values",
     ylab = "Pearson residuals",
     cex.lab = 1.5)
abline(h = 0, lty = 2)
# fitted ranges 300 to 700 (whole range should be 65 to 1260)
# residuals evenly spread
# good range of residuals (-2 ro 2)

# Plot residuals versus site
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

I1a <- inla(avgwt_g ~ year + soaktimec +
              f(site, model='iid'),
            control.compute = list(dic = TRUE, cpo=TRUE, waic=TRUE, config=TRUE),
            family = "gamma",
            control.family=list(link='log'),
            data = catchwb6_pred) # now pred
summary(I1a) # should be the same as before

#extract the predictor values associated with the appended data
catchwb6_test1a <- cbind(catchwb6_test, 
                         I1a$summary.linear.predictor[(nrow(catchwb6_train)+1):nrow(catchwb6_pred),])
head(catchwb6_test1a)

catchwb6_test1a <- reshape:::rename(catchwb6_test1a, c("0.025quant"="lower", "0.975quant"="upper"))

ggplot(catchwb6_test1a, aes(y=mean, x=soaktimec)) +
  geom_point(data=catchwb6_train, aes(y=avgwt_g)) + # obs
  geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.2) + # predicted SD
  geom_line() # predicted
# MUCH smaller than observed

ggplot(catchwb6_test1a, aes(y=mean, x=soaktimec)) +
  geom_point(aes(y=resp_obsB)) + # obs
  geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.2) + # predicted SD
  geom_line()

# *Gamma GLM envr------
# not used
I1 <- inla(avgwt_g ~ year + soaktimec + temp_Cc + turb_FNUc + distance_to_southc + distance_to_shorec +
             f(site, model='iid'),
           control.compute = list(dic = TRUE, cpo=TRUE, waic=TRUE, config=TRUE),
           control.family = list(link='log'),
           family = "gamma",
           data = catchwb6_train)
summary(I1)
-sum(log(I1$cpo$cpo)) # smaller is better fit of model to data

# Posterior mean values and 95% CI for the regression parameters:
Beta1 <- I1$summary.fixed[,c("mean", "sd", "0.025quant", "0.975quant")] 
print(Beta1, digits = 2)  
# d2S and year does not contain zero

# Hyperparameters
summary(I1)
HyperPar.I1 <- bri.hyperpar.summary(I1)
round(HyperPar.I1, digits = 3)

# This is the parameter r from the Gamma variance
# dispersion stat - denominator in varaiance of gamma distribution
r  <- I1$summary.hyperpar["Precision parameter for the Gamma observations", "mean"]
r



# Get the fitted values and Pearson residuals
N     <- nrow(catchwb6_train)
mu1   <- I1$summary.fitted.values[1:N,"mean"] 
VarY1 <- mu1^2 / r
E1    <- (catchwb6_train$avgwt_g - mu1) / sqrt(VarY1)

# A Gamma GLM(M) cannot be overdispersed


range(catchwb6_train$avgwt_g) # 65.5-1256.923

# why overdispersed
# fitted values versus the Pearson residuals
par(mfrow=c(1,1))
plot(x = mu1,
     y = E1,
     xlab = "Fitted values",
     ylab = "Pearson residuals",
     cex.lab = 1.5)
abline(h = 0, lty = 2)
# fitted ranges 300 to 700 (whole range should be 65 to 1260)
# residuals evenly spread
# good range of residuals (-2 ro 2)

# Plot residuals versus site
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

I1a <- inla(avgwt_g ~ year + soaktimec + temp_Cc  + turb_FNUc + distance_to_southc +
              f(site, model='iid'),
            control.compute = list(dic = TRUE, cpo=TRUE, waic=TRUE, config=TRUE),
            family = "gamma",
            control.family=list(link='log'),
            data = catchwb6_pred) # now pred
summary(I1a) # should be the same as before

#extract the predictor values associated with the appended data
catchwb6_test1a <- cbind(catchwb6_test, 
                         I1a$summary.linear.predictor[(nrow(catchwb6_train)+1):nrow(catchwb6_pred),])
head(catchwb6_test1a)

catchwb6_test1a <- reshape:::rename(catchwb6_test1a, c("0.025quant"="lower", "0.975quant"="upper"))

ggplot(catchwb6_test1a, aes(y=mean, x=soaktimec)) +
  geom_point(data=catchwb6_train, aes(y=avgwt_g)) + # obs
  geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.2) + # predicted SD
  geom_line() # predicted
# MUCH smaller than observed

ggplot(catchwb6_test1a, aes(y=mean, x=temp_Cc)) +
  geom_point(data=catchwb6_train, aes(y=avgwt_g)) + # obs
  geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.2) + # predicted SD
  geom_line()

ggplot(catchwb6_test1a, aes(y=mean, x=soaktimec)) +
  geom_point(aes(y=resp_obsB)) + # obs
  geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.2) + # predicted SD
  geom_line()

catchwb6_test1a$soaktime2 <- catchwb6_test1a$soaktimec * sd(catchwb6_train$soaktime) + mean(catchwb6_train$soaktime)
p <- ggplot()
p <- p + geom_point(data = catchwb6_train, 
                    aes(y = avgwt_g, x = soaktime),
                    size = 1)
p <- p + xlab("soak time (h)") + ylab("total number LKWF")
p <- p + theme(text = element_text(size=15)) 
p <- p + geom_point(data = catchwb6_test1a, 
                    aes(x = soaktime2, 
                        y = mean), 
                    col = "red")
p <- p + geom_ribbon(data = catchwb6_test1a, 
                     aes(x = soaktime2, 
                         ymax = upper, 
                         ymin = lower),
                     alpha = 0.2)
p

ggplot(data=catchwb6_test1a) +
  geom_point(mapping=aes(x=resp_obsB, y=mean)) +
  xlab("Observed") +
  ylab("Predicted") 

# Gamma GLM d2S------
# not used
I1 <- inla(avgwt_g ~ year + soaktimec + distance_to_southc +
             f(site, model='iid'),
           control.compute = list(dic = TRUE, cpo=TRUE, waic=TRUE, config=TRUE),
           control.family = list(link='log'),
           family = "gamma",
           data = catchwb6_train)
summary(I1)
-sum(log(I1$cpo$cpo)) # smaller is better fit of model to data

# Posterior mean values and 95% CI for the regression parameters:
Beta1 <- I1$summary.fixed[,c("mean", "sd", "0.025quant", "0.975quant")] 
print(Beta1, digits = 2)  
# d2S does not contain zero

# Hyperparameters
summary(I1)
HyperPar.I1 <- bri.hyperpar.summary(I1)
round(HyperPar.I1, digits = 3)

# This is the parameter r from the Gamma variance
# dispersion stat - denominator in varaiance of gamma distribution
r  <- I1$summary.hyperpar["Precision parameter for the Gamma observations", "mean"]
r



# Get the fitted values and Pearson residuals
N     <- nrow(catchwb6_train)
mu1   <- I1$summary.fitted.values[1:N,"mean"] 
VarY1 <- mu1^2 / r
E1    <- (catchwb6_train$avgwt_g - mu1) / sqrt(VarY1)

# A Gamma GLM(M) cannot be overdispersed


range(catchwb6_train$avgwt_g) # 65.5-1256.923

# why overdispersed
# fitted values versus the Pearson residuals
par(mfrow=c(1,1))
plot(x = mu1,
     y = E1,
     xlab = "Fitted values",
     ylab = "Pearson residuals",
     cex.lab = 1.5)
abline(h = 0, lty = 2)
# fitted ranges 300 to 700 (whole range should be 65 to 1260)
# residuals evenly spread
# good range of residuals (-2 ro 2)

# Plot residuals versus site
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

I1a <- inla(avgwt_g ~ year + soaktimec +  distance_to_southc +
              f(site, model='iid'),
            control.compute = list(dic = TRUE, cpo=TRUE, waic=TRUE, config=TRUE),
            family = "gamma",
            control.family=list(link='log'),
            data = catchwb6_pred) # now pred
summary(I1a) # should be the same as before

#extract the predictor values associated with the appended data
catchwb6_test1a <- cbind(catchwb6_test, 
                         I1a$summary.linear.predictor[(nrow(catchwb6_train)+1):nrow(catchwb6_pred),])
head(catchwb6_test1a)

catchwb6_test1a <- reshape:::rename(catchwb6_test1a, c("0.025quant"="lower", "0.975quant"="upper"))

ggplot(catchwb6_test1a, aes(y=mean, x=soaktimec)) +
  geom_point(data=catchwb6_train, aes(y=avgwt_g)) + # obs
  geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.2) + # predicted SD
  geom_line() # predicted
# MUCH smaller than observed

ggplot(catchwb6_test1a, aes(y=mean, x=temp_Cc)) +
  geom_point(data=catchwb6_train, aes(y=avgwt_g)) + # obs
  geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.2) + # predicted SD
  geom_line()

ggplot(catchwb6_test1a, aes(y=mean, x=soaktimec)) +
  geom_point(aes(y=resp_obsB)) + # obs
  geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.2) + # predicted SD
  geom_line()

catchwb6_test1a$soaktime2 <- catchwb6_test1a$soaktimec * sd(catchwb6_train$soaktime) + mean(catchwb6_train$soaktime)
p <- ggplot()
p <- p + geom_point(data = catchwb6_train, 
                    aes(y = avgwt_g, x = soaktime),
                    size = 1)
p <- p + xlab("soak time (h)") + ylab("total number LKWF")
p <- p + theme(text = element_text(size=15)) 
p <- p + geom_point(data = catchwb6_test1a, 
                    aes(x = soaktime2, 
                        y = mean), 
                    col = "red")
p <- p + geom_ribbon(data = catchwb6_test1a, 
                     aes(x = soaktime2, 
                         ymax = upper, 
                         ymin = lower),
                     alpha = 0.2)
p

ggplot(data=catchwb6_test1a) +
  geom_point(mapping=aes(x=resp_obsB, y=mean)) +
  xlab("Observed") +
  ylab("Predicted") 

# *Gamma GLM turb ------
# not used
I1 <- inla(avgwt_g ~ year + soaktimec +  turb_FNUc +
             f(site, model='iid'),
           control.compute = list(dic = TRUE, cpo=TRUE, waic=TRUE, config=TRUE),
           control.family = list(link='log'),
           family = "gamma",
           data = catchwb6_train)
summary(I1)
-sum(log(I1$cpo$cpo)) # smaller is better fit of model to data

# Posterior mean values and 95% CI for the regression parameters:
Beta1 <- I1$summary.fixed[,c("mean", "sd", "0.025quant", "0.975quant")] 
print(Beta1, digits = 2)  


# Hyperparameters
summary(I1)
HyperPar.I1 <- bri.hyperpar.summary(I1)
round(HyperPar.I1, digits = 3)

# This is the parameter r from the Gamma variance
# dispersion stat - denominator in varaiance of gamma distribution
r  <- I1$summary.hyperpar["Precision parameter for the Gamma observations", "mean"]
r



# Get the fitted values and Pearson residuals
N     <- nrow(catchwb6_train)
mu1   <- I1$summary.fitted.values[1:N,"mean"] 
VarY1 <- mu1^2 / r
E1    <- (catchwb6_train$avgwt_g - mu1) / sqrt(VarY1)

# A Gamma GLM(M) cannot be overdispersed


range(catchwb6_train$avgwt_g) # 65.5-1256.923

# why overdispersed
# fitted values versus the Pearson residuals
par(mfrow=c(1,1))
plot(x = mu1,
     y = E1,
     xlab = "Fitted values",
     ylab = "Pearson residuals",
     cex.lab = 1.5)
abline(h = 0, lty = 2)
# fitted ranges 300 to 700 (whole range should be 65 to 1260)
# residuals evenly spread
# good range of residuals (-2 ro 2)

# Plot residuals versus site
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

I1a <- inla(avgwt_g ~ year + soaktimec +  turb_FNUc + 
              f(site, model='iid'),
            control.compute = list(dic = TRUE, cpo=TRUE, waic=TRUE, config=TRUE),
            family = "gamma",
            control.family=list(link='log'),
            data = catchwb6_pred) # now pred
summary(I1a) # should be the same as before

#extract the predictor values associated with the appended data
catchwb6_test1a <- cbind(catchwb6_test, 
                         I1a$summary.linear.predictor[(nrow(catchwb6_train)+1):nrow(catchwb6_pred),])
head(catchwb6_test1a)

catchwb6_test1a <- reshape:::rename(catchwb6_test1a, c("0.025quant"="lower", "0.975quant"="upper"))

ggplot(catchwb6_test1a, aes(y=mean, x=soaktimec)) +
  geom_point(data=catchwb6_train, aes(y=avgwt_g)) + # obs
  geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.2) + # predicted SD
  geom_line() # predicted
# MUCH smaller than observed

ggplot(catchwb6_test1a, aes(y=mean, x=temp_Cc)) +
  geom_point(data=catchwb6_train, aes(y=avgwt_g)) + # obs
  geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.2) + # predicted SD
  geom_line()

ggplot(catchwb6_test1a, aes(y=mean, x=soaktimec)) +
  geom_point(aes(y=resp_obsB)) + # obs
  geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.2) + # predicted SD
  geom_line()

catchwb6_test1a$soaktime2 <- catchwb6_test1a$soaktimec * sd(catchwb6_train$soaktime) + mean(catchwb6_train$soaktime)
p <- ggplot()
p <- p + geom_point(data = catchwb6_train, 
                    aes(y = avgwt_g, x = soaktime),
                    size = 1)
p <- p + xlab("soak time (h)") + ylab("total number LKWF")
p <- p + theme(text = element_text(size=15)) 
p <- p + geom_point(data = catchwb6_test1a, 
                    aes(x = soaktime2, 
                        y = mean), 
                    col = "red")
p <- p + geom_ribbon(data = catchwb6_test1a, 
                     aes(x = soaktime2, 
                         ymax = upper, 
                         ymin = lower),
                     alpha = 0.2)
p

ggplot(data=catchwb6_test1a) +
  geom_point(mapping=aes(x=resp_obsB, y=mean)) +
  xlab("Observed") +
  ylab("Predicted") 

# *Gamma GLM turb d2S------
# not used
I1 <- inla(avgwt_g ~ year + soaktimec +  turb_FNUc + distance_to_southc +
             f(site, model='iid'),
           control.compute = list(dic = TRUE, cpo=TRUE, waic=TRUE, config=TRUE),
           control.family = list(link='log'),
           family = "gamma",
           data = catchwb6_train)
summary(I1)
-sum(log(I1$cpo$cpo)) # smaller is better fit of model to data

# Posterior mean values and 95% CI for the regression parameters:
Beta1 <- I1$summary.fixed[,c("mean", "sd", "0.025quant", "0.975quant")] 
print(Beta1, digits = 2)  


# Hyperparameters
summary(I1)
HyperPar.I1 <- bri.hyperpar.summary(I1)
round(HyperPar.I1, digits = 3)

# This is the parameter r from the Gamma variance
# dispersion stat - denominator in varaiance of gamma distribution
r  <- I1$summary.hyperpar["Precision parameter for the Gamma observations", "mean"]
r



# Get the fitted values and Pearson residuals
N     <- nrow(catchwb6_train)
mu1   <- I1$summary.fitted.values[1:N,"mean"] 
VarY1 <- mu1^2 / r
E1    <- (catchwb6_train$avgwt_g - mu1) / sqrt(VarY1)

# A Gamma GLM(M) cannot be overdispersed


range(catchwb6_train$avgwt_g) # 65.5-1256.923

# why overdispersed
# fitted values versus the Pearson residuals
par(mfrow=c(1,1))
plot(x = mu1,
     y = E1,
     xlab = "Fitted values",
     ylab = "Pearson residuals",
     cex.lab = 1.5)
abline(h = 0, lty = 2)
# fitted ranges 300 to 700 (whole range should be 65 to 1260)
# residuals evenly spread
# good range of residuals (-2 ro 2)

# Plot residuals versus site
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

I1a <- inla(avgwt_g ~ year + soaktimec +  turb_FNUc + distance_to_southc +
              f(site, model='iid'),
            control.compute = list(dic = TRUE, cpo=TRUE, waic=TRUE, config=TRUE),
            family = "gamma",
            control.family=list(link='log'),
            data = catchwb6_pred) # now pred
summary(I1a) # should be the same as before

#extract the predictor values associated with the appended data
catchwb6_test1a <- cbind(catchwb6_test, 
                         I1a$summary.linear.predictor[(nrow(catchwb6_train)+1):nrow(catchwb6_pred),])
head(catchwb6_test1a)

catchwb6_test1a <- reshape:::rename(catchwb6_test1a, c("0.025quant"="lower", "0.975quant"="upper"))

ggplot(catchwb6_test1a, aes(y=mean, x=soaktimec)) +
  geom_point(data=catchwb6_train, aes(y=avgwt_g)) + # obs
  geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.2) + # predicted SD
  geom_line() # predicted
# MUCH smaller than observed

ggplot(catchwb6_test1a, aes(y=mean, x=temp_Cc)) +
  geom_point(data=catchwb6_train, aes(y=avgwt_g)) + # obs
  geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.2) + # predicted SD
  geom_line()

ggplot(catchwb6_test1a, aes(y=mean, x=soaktimec)) +
  geom_point(aes(y=resp_obsB)) + # obs
  geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.2) + # predicted SD
  geom_line()

catchwb6_test1a$soaktime2 <- catchwb6_test1a$soaktimec * sd(catchwb6_train$soaktime) + mean(catchwb6_train$soaktime)
p <- ggplot()
p <- p + geom_point(data = catchwb6_train, 
                    aes(y = avgwt_g, x = soaktime),
                    size = 1)
p <- p + xlab("soak time (h)") + ylab("total number LKWF")
p <- p + theme(text = element_text(size=15)) 
p <- p + geom_point(data = catchwb6_test1a, 
                    aes(x = soaktime2, 
                        y = mean), 
                    col = "red")
p <- p + geom_ribbon(data = catchwb6_test1a, 
                     aes(x = soaktime2, 
                         ymax = upper, 
                         ymin = lower),
                     alpha = 0.2)
p

ggplot(data=catchwb6_test1a) +
  geom_point(mapping=aes(x=resp_obsB, y=mean)) +
  xlab("Observed") +
  ylab("Predicted") 


# *Gamma GLM turb temp ------
# not used
I1 <- inla(avgwt_g ~ year + soaktimec +  turb_FNUc + temp_Cc +
             f(site, model='iid'),
           control.compute = list(dic = TRUE, cpo=TRUE, waic=TRUE, config=TRUE),
           control.family = list(link='log'),
           family = "gamma",
           data = catchwb6_train)
summary(I1)
-sum(log(I1$cpo$cpo)) # smaller is better fit of model to data

# Posterior mean values and 95% CI for the regression parameters:
Beta1 <- I1$summary.fixed[,c("mean", "sd", "0.025quant", "0.975quant")] 
print(Beta1, digits = 2)  


# Hyperparameters
summary(I1)
HyperPar.I1 <- bri.hyperpar.summary(I1)
round(HyperPar.I1, digits = 3)

# This is the parameter r from the Gamma variance
# dispersion stat - denominator in varaiance of gamma distribution
r  <- I1$summary.hyperpar["Precision parameter for the Gamma observations", "mean"]
r



# Get the fitted values and Pearson residuals
N     <- nrow(catchwb6_train)
mu1   <- I1$summary.fitted.values[1:N,"mean"] 
VarY1 <- mu1^2 / r
E1    <- (catchwb6_train$avgwt_g - mu1) / sqrt(VarY1)

# A Gamma GLM(M) cannot be overdispersed


range(catchwb6_train$avgwt_g) # 65.5-1256.923

# why overdispersed
# fitted values versus the Pearson residuals
par(mfrow=c(1,1))
plot(x = mu1,
     y = E1,
     xlab = "Fitted values",
     ylab = "Pearson residuals",
     cex.lab = 1.5)
abline(h = 0, lty = 2)
# fitted ranges 300 to 700 (whole range should be 65 to 1260)
# residuals evenly spread
# good range of residuals (-2 ro 2)

# Plot residuals versus site
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

I1a <- inla(avgwt_g ~ year + soaktimec +  turb_FNUc + temp_Cc +
              f(site, model='iid'),
            control.compute = list(dic = TRUE, cpo=TRUE, waic=TRUE, config=TRUE),
            family = "gamma",
            control.family=list(link='log'),
            data = catchwb6_pred) # now pred
summary(I1a) # should be the same as before

#extract the predictor values associated with the appended data
catchwb6_test1a <- cbind(catchwb6_test, 
                         I1a$summary.linear.predictor[(nrow(catchwb6_train)+1):nrow(catchwb6_pred),])
head(catchwb6_test1a)

catchwb6_test1a <- reshape:::rename(catchwb6_test1a, c("0.025quant"="lower", "0.975quant"="upper"))

ggplot(catchwb6_test1a, aes(y=mean, x=soaktimec)) +
  geom_point(data=catchwb6_train, aes(y=avgwt_g)) + # obs
  geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.2) + # predicted SD
  geom_line() # predicted
# MUCH smaller than observed

ggplot(catchwb6_test1a, aes(y=mean, x=temp_Cc)) +
  geom_point(data=catchwb6_train, aes(y=avgwt_g)) + # obs
  geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.2) + # predicted SD
  geom_line()

ggplot(catchwb6_test1a, aes(y=mean, x=soaktimec)) +
  geom_point(aes(y=resp_obsB)) + # obs
  geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.2) + # predicted SD
  geom_line()

catchwb6_test1a$soaktime2 <- catchwb6_test1a$soaktimec * sd(catchwb6_train$soaktime) + mean(catchwb6_train$soaktime)
p <- ggplot()
p <- p + geom_point(data = catchwb6_train, 
                    aes(y = avgwt_g, x = soaktime),
                    size = 1)
p <- p + xlab("soak time (h)") + ylab("total number LKWF")
p <- p + theme(text = element_text(size=15)) 
p <- p + geom_point(data = catchwb6_test1a, 
                    aes(x = soaktime2, 
                        y = mean), 
                    col = "red")
p <- p + geom_ribbon(data = catchwb6_test1a, 
                     aes(x = soaktime2, 
                         ymax = upper, 
                         ymin = lower),
                     alpha = 0.2)
p

ggplot(data=catchwb6_test1a) +
  geom_point(mapping=aes(x=resp_obsB, y=mean)) +
  xlab("Observed") +
  ylab("Predicted") 
# *Gamma GLM temp d2S------
# not used
I1 <- inla(avgwt_g ~ year + soaktimec +  distance_to_southc + temp_Cc +
             f(site, model='iid'),
           control.compute = list(dic = TRUE, cpo=TRUE, waic=TRUE, config=TRUE),
           control.family = list(link='log'),
           family = "gamma",
           data = catchwb6_train)
summary(I1)
-sum(log(I1$cpo$cpo)) # smaller is better fit of model to data

# Posterior mean values and 95% CI for the regression parameters:
Beta1 <- I1$summary.fixed[,c("mean", "sd", "0.025quant", "0.975quant")] 
print(Beta1, digits = 2)  


# Hyperparameters
summary(I1)
HyperPar.I1 <- bri.hyperpar.summary(I1)
round(HyperPar.I1, digits = 3)

# This is the parameter r from the Gamma variance
# dispersion stat - denominator in varaiance of gamma distribution
r  <- I1$summary.hyperpar["Precision parameter for the Gamma observations", "mean"]
r



# Get the fitted values and Pearson residuals
N     <- nrow(catchwb6_train)
mu1   <- I1$summary.fitted.values[1:N,"mean"] 
VarY1 <- mu1^2 / r
E1    <- (catchwb6_train$avgwt_g - mu1) / sqrt(VarY1)

# A Gamma GLM(M) cannot be overdispersed


range(catchwb6_train$avgwt_g) # 65.5-1256.923

# why overdispersed
# fitted values versus the Pearson residuals
par(mfrow=c(1,1))
plot(x = mu1,
     y = E1,
     xlab = "Fitted values",
     ylab = "Pearson residuals",
     cex.lab = 1.5)
abline(h = 0, lty = 2)
# fitted ranges 300 to 700 (whole range should be 65 to 1260)
# residuals evenly spread
# good range of residuals (-2 ro 2)

# Plot residuals versus site
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

I1a <- inla(avgwt_g ~ year + soaktimec +  distance_to_southc + temp_Cc +
              f(site, model='iid'),
            control.compute = list(dic = TRUE, cpo=TRUE, waic=TRUE, config=TRUE),
            family = "gamma",
            control.family=list(link='log'),
            data = catchwb6_pred) # now pred
summary(I1a) # should be the same as before

#extract the predictor values associated with the appended data
catchwb6_test1a <- cbind(catchwb6_test, 
                         I1a$summary.linear.predictor[(nrow(catchwb6_train)+1):nrow(catchwb6_pred),])
head(catchwb6_test1a)

catchwb6_test1a <- reshape:::rename(catchwb6_test1a, c("0.025quant"="lower", "0.975quant"="upper"))

ggplot(catchwb6_test1a, aes(y=mean, x=soaktimec)) +
  geom_point(data=catchwb6_train, aes(y=avgwt_g)) + # obs
  geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.2) + # predicted SD
  geom_line() # predicted
# MUCH smaller than observed

ggplot(catchwb6_test1a, aes(y=mean, x=temp_Cc)) +
  geom_point(data=catchwb6_train, aes(y=avgwt_g)) + # obs
  geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.2) + # predicted SD
  geom_line()

ggplot(catchwb6_test1a, aes(y=mean, x=soaktimec)) +
  geom_point(aes(y=resp_obsB)) + # obs
  geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.2) + # predicted SD
  geom_line()

catchwb6_test1a$soaktime2 <- catchwb6_test1a$soaktimec * sd(catchwb6_train$soaktime) + mean(catchwb6_train$soaktime)
p <- ggplot()
p <- p + geom_point(data = catchwb6_train, 
                    aes(y = avgwt_g, x = soaktime),
                    size = 1)
p <- p + xlab("soak time (h)") + ylab("total number LKWF")
p <- p + theme(text = element_text(size=15)) 
p <- p + geom_point(data = catchwb6_test1a, 
                    aes(x = soaktime2, 
                        y = mean), 
                    col = "red")
p <- p + geom_ribbon(data = catchwb6_test1a, 
                     aes(x = soaktime2, 
                         ymax = upper, 
                         ymin = lower),
                     alpha = 0.2)
p

ggplot(data=catchwb6_test1a) +
  geom_point(mapping=aes(x=resp_obsB, y=mean)) +
  xlab("Observed") +
  ylab("Predicted") 


# *Gamma GLM temp------
I1 <- inla(avgwt_g ~ year + soaktimec + temp_Cc +
             f(site, model='iid'),
           control.compute = list(dic = TRUE, cpo=TRUE, waic=TRUE, config=TRUE),
           control.family = list(link='log'),
           family = "gamma",
           data = catchwb6_train)
summary(I1)
-sum(log(I1$cpo$cpo)) # smaller is better fit of model to data

# Posterior mean values and 95% CI for the regression parameters:
Beta1 <- I1$summary.fixed[,c("mean", "sd", "0.025quant", "0.975quant")] 
print(Beta1, digits = 2)  

# Hyperparameters
summary(I1)
HyperPar.I1 <- bri.hyperpar.summary(I1)
round(HyperPar.I1, digits = 3)

# This is the parameter r from the Gamma variance
# dispersion stat - denominator in varaiance of gamma distribution
r  <- I1$summary.hyperpar["Precision parameter for the Gamma observations", "mean"]
r



# Get the fitted values and Pearson residuals
N     <- nrow(catchwb6_train)
mu1   <- I1$summary.fitted.values[1:N,"mean"] 
VarY1 <- mu1^2 / r
E1    <- (catchwb6_train$avgwt_g - mu1) / sqrt(VarY1)

# A Gamma GLM(M) cannot be overdispersed


range(catchwb6_train$avgwt_g) # 65.5-1256.923

# why overdispersed
# fitted values versus the Pearson residuals
par(mfrow=c(1,1))
plot(x = mu1,
     y = E1,
     xlab = "Fitted values",
     ylab = "Pearson residuals",
     cex.lab = 1.5)
abline(h = 0, lty = 2)
# fitted ranges 300 to 700 (whole range should be 65 to 1260)
# residuals evenly spread
# good range of residuals (-2 ro 2)

# Plot residuals versus site
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

I1a <- inla(avgwt_g ~ year + soaktimec + temp_Cc  +
              f(site, model='iid'),
            control.compute = list(dic = TRUE, cpo=TRUE, waic=TRUE, config=TRUE),
            family = "gamma",
            control.family=list(link='log'),
            data = catchwb6_pred) # now pred
summary(I1a) # should be the same as before

#extract the predictor values associated with the appended data
catchwb6_test1a <- cbind(catchwb6_test, 
                         I1a$summary.linear.predictor[(nrow(catchwb6_train)+1):nrow(catchwb6_pred),])
head(catchwb6_test1a)

catchwb6_test1a <- reshape:::rename(catchwb6_test1a, c("0.025quant"="lower", "0.975quant"="upper"))

ggplot(catchwb6_test1a, aes(y=mean, x=soaktimec)) +
  geom_point(data=catchwb6_train, aes(y=avgwt_g)) + # obs
  geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.2) + # predicted SD
  geom_line() # predicted
# MUCH smaller than observed

ggplot(catchwb6_test1a, aes(y=mean, x=temp_Cc)) +
  geom_point(data=catchwb6_train, aes(y=avgwt_g)) + # obs
  geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.2) + # predicted SD
  geom_line()

ggplot(catchwb6_test1a, aes(y=mean, x=soaktimec)) +
  geom_point(aes(y=resp_obsB)) + # obs
  geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.2) + # predicted SD
  geom_line()

catchwb6_test1a$soaktime2 <- catchwb6_test1a$soaktimec * sd(catchwb6_train$soaktime) + mean(catchwb6_train$soaktime)
p <- ggplot()
p <- p + geom_point(data = catchwb6_train, 
                    aes(y = avgwt_g, x = soaktime),
                    size = 1)
p <- p + xlab("soak time (h)") + ylab("avg wt LKWF")
p <- p + theme(text = element_text(size=15)) 
p <- p + geom_point(data = catchwb6_test1a, 
                    aes(x = soaktime2, 
                        y = mean), 
                    col = "red")
p <- p + geom_ribbon(data = catchwb6_test1a, 
                     aes(x = soaktime2, 
                         ymax = upper, 
                         ymin = lower),
                     alpha = 0.2)
p

ggplot(data=catchwb6_test1a) +
  geom_point(mapping=aes(x=resp_obsB, y=mean)) +
  xlab("Observed") +
  ylab("Predicted") 

# Gamma GLM ind srf --------
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
# Step 2 of making a mesh: The primary tool to control the shape of the 
# triangles is max.edge. Other useful arguments are the cutoff and offset. 
# Use an outer area to avoid a boundary effect. The outer area can be less fine 
# than the inner area to reduce computing time.

# Come up with an initial guess for the range:
#  At what distance is the correlation smaller than 0.1?
#  20 km? See also the histogram if the distances again.
#  If we choose 20 km then the correlation will affect around 10%
#  of the combinations. The value of 20 km also complies (vaguely) with
#  the sample variogram of the random intercepts (maybe 30 would have been better).

# So..let's use 20 km for the range?
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
range(catchwb6_train$avgwt_g)
# P(Range < 5 km ) = 0.05
# P(sigma > 2) = 0.05
# SB = exp(u_i)
# some u_i have to be as large as 5.24 to cover 190
# If u_i ~ N(0, sigma_u^2) then it is unlikley that sigma_u > 1.5
#P(sigma > 1.5) = 0.05
M1 <- glm(avgwt_g ~ 1, data = catchwb6_train)
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
Xcov <- as.data.frame(Xcov)

Stack <- inla.stack(
  tag = "Fit",
  data = list(avgwt_g = catchwb6_train$avgwt_g),  
  A = list(1, 1, 1, A),                 
  effects = list(                 
    Intercept = rep(1,N),
    Xcov  = Xcov[,-1],               
    site = catchwb6_train$site, # random effect
    w  = w.index))  # SRF


f14 <- avgwt_g ~ -1 + Intercept + soaktimec  + year +
                f(site, model="iid") +
                 f(w, model = spde)


# 7. Run the spatial model in INLA.
I14 <- inla(f14,
            family = "gamma", 
            data=inla.stack.data(Stack),
            control.compute = list(dic = TRUE,  waic=TRUE, config = TRUE),
            control.family=list(link='log'),
            control.predictor = list(A = inla.stack.A(Stack)))

summary(I14)


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

w.pm        <- I14$summary.random$w$mean  #mu.srf # exp converts from the log link (value of 1 is = exp(1), but now all exp())
w.proj      <- inla.mesh.projector(mesh) 
w.pm100_100 <- inla.mesh.project(w.proj, w.pm)
long        <- as.data.frame(w.proj$lattice$loc)
colnames(long) <- c("EastingKM", "NorthingKM")
long <- long %>%
  mutate(w.pm = as.vector(w.pm100_100))

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


# Here are the fitted values and Pearson residuals
N     <- nrow(catchwb6_train)
mu14   <- I14$summary.fitted.values[1:N,"mean"] 
r     <- I14$summary.hyperpar["Precision parameter for the Gamma observations", "mean"]
VarY14 <- mu14^2 / r
E14    <- (catchwb6_train$avgwt_g - mu14) / sqrt(VarY14)



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
     y = catchwb6_train$avgwt_g,
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
p <- p + geom_smooth(data = V14, se = FALSE, span = 0.9, aes(x = dist, y = gamma)) 
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
  data = list(y = catchwb6_train$avgwt_g),
  A = list(1,1,1,A),
  effects = list(
    Intercept = rep(1,N),
    Xcov2 = Xcov2[,-1], 
    site = catchwb6_train$site,
    w  = w.index))


f2 <- y ~ -1 + Intercept +
  soaktimec  + year+
  f(site, model='iid') +
  f(w, model = spde)


I2 <- inla(f2,
           family = "gamma",
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
  soaktimec  + year + 
  f(site, model="iid") +
  f(w, model = spde)

# Note that we use a different stack! The rest is the same
I5.pred <- inla(f5,
                family = "gamma", 
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
  mutate(set = ifelse(!is.na(avgwt_g), "Training", "Testing"))
catchwb6_pred$set <- as.factor(catchwb6_pred$set)


catchwb6_pred <- catchwb6_pred %>% 
  mutate(avgwt_g=ifelse(!is.na(avgwt_g), avgwt_g, mean))

head(catchwb6_pred)

p <- ggplot()
p <- p + geom_smooth(data = catchwb6_test2, 
                     aes(x = soaktime, 
                         y = mean), alpha=0.5)
p <- p + geom_point(data = catchwb6_pred, 
                    aes(y = avgwt_g, x = soaktime, col=set),
                    size = 3)
p <- p + xlab("Soak time (h)") + ylab("Mean weight (g)") + labs ( col = "Dataset")
p <- p + theme(text = element_text(size=15), legend.position = "bottom") 
# p <- p + geom_ribbon(data = catchwb6_test3a,
#                      aes(x = soaktime2,
#                          ymax = upper,
#                          ymin = lower),
#                      alpha = 0.2)

p


ggplot(data=catchwb6_test2) +
  geom_point(mapping=aes(x=resp_obsB, y=mu, colour=as.factor(catchwb6_test2$year)), size=3) +
  xlab("Observed") +
  ylab("Predicted") +
  geom_abline() +
  xlim(c(0,800)) +
  ylim(c(0,800)) 




# *Gamma GLM envrind srf --------
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
# Step 2 of making a mesh: The primary tool to control the shape of the 
# triangles is max.edge. Other useful arguments are the cutoff and offset. 
# Use an outer area to avoid a boundary effect. The outer area can be less fine 
# than the inner area to reduce computing time.

# Come up with an initial guess for the range:
#  At what distance is the correlation smaller than 0.1?
#  20 km? See also the histogram if the distances again.
#  If we choose 20 km then the correlation will affect around 10%
#  of the combinations. The value of 20 km also complies (vaguely) with
#  the sample variogram of the random intercepts (maybe 30 would have been better).

# So..let's use 20 km for the range?
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
range(catchwb6_train$avgwt_g)
# P(Range < 5 km ) = 0.05
# P(sigma > 2) = 0.05
# SB = exp(u_i)
# some u_i have to be as large as 5.24 to cover 190
# If u_i ~ N(0, sigma_u^2) then it is unlikley that sigma_u > 1.5
#P(sigma > 1.5) = 0.05
M1 <- glm(avgwt_g ~ 1, data = catchwb6_train)
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
Xcov <- model.matrix(~ soaktimec + year + temp_Cc + turb_FNUc + distance_to_southc,
                     data = catchwb6_train)
Xcov <- as.data.frame(Xcov)

Stack <- inla.stack(
  tag = "Fit",
  data = list(avgwt_g = catchwb6_train$avgwt_g),  
  A = list(1, 1, 1, A),                 
  effects = list(                 
    Intercept = rep(1,N),
    Xcov  = Xcov[,-1],               
    site = catchwb6_train$site, # random effect
    w  = w.index))  # SRF


f14 <- avgwt_g ~ -1 + Intercept + soaktimec  + year + temp_Cc + turb_FNUc + distance_to_southc +
   f(site, model="iid") +
  f(w, model = spde)


# 7. Run the spatial model in INLA.
I14 <- inla(f14,
            family = "gamma", 
            data=inla.stack.data(Stack),
            control.compute = list(dic = TRUE,  waic=TRUE, config = TRUE),
            control.family=list(link='log'),
            control.predictor = list(A = inla.stack.A(Stack)))

summary(I14)


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


# Here are the fitted values and Pearson residuals
N     <- nrow(catchwb6_train)
mu14   <- I14$summary.fitted.values[1:N,"mean"] 
r     <- I14$summary.hyperpar["Precision parameter for the Gamma observations", "mean"]
VarY14 <- mu14^2 / r
E14    <- (catchwb6_train$avgwt_g - mu14) / sqrt(VarY14)



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
     y = catchwb6_train$avgwt_g,
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


Xcov <- model.matrix(~ soaktimec + year + temp_Cc + turb_FNUc + distance_to_southc,
                     data = catchwb6_train)
Xcov2 <- as.data.frame(Xcov)
colnames(Xcov2)

N    <- nrow(catchwb6_train)

Stack <- inla.stack(
  tag = "Fit",
  data = list(y = catchwb6_train$avgwt_g),
  A = list(1,1,1,A),
  effects = list(
    Intercept = rep(1,N),
    Xcov2 = Xcov2[,-1], 
    site = catchwb6_train$site,
    w  = w.index))


f2 <- y ~ -1 + Intercept +
  soaktimec  +  year + temp_Cc + turb_FNUc + distance_to_southc +
  f(site, model='iid') +
  f(w, model = spde)


I2 <- inla(f2,
           family = "gamma",
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
Xp <- model.matrix(~ soaktimec + year + temp_Cc + turb_FNUc + distance_to_southc, data = catchwb6_test)

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
  soaktimec  + year + temp_Cc + turb_FNUc + distance_to_southc +
    f(site, model="iid") +
  f(w, model = spde)

# Note that we use a different stack! The rest is the same
I5.pred <- inla(f5,
                family = "gamma", 
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
                    aes(y = avgwt_g, x = soaktime),
                    size = 1, colour="black")
p <- p + xlab("soak time (h)") + ylab("avg wt LKWF")
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
  geom_point(mapping=aes(x=resp_obsB, y=mu, colour=as.factor(catchwb6_test2$year)), size=3) +
  xlab("Observed") +
  ylab("Predicted") +
  geom_abline() +
  xlim(c(0,800)) +
  ylim(c(0,800)) 



# *Gamma GLM all env ind srf --------
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
# Step 2 of making a mesh: The primary tool to control the shape of the 
# triangles is max.edge. Other useful arguments are the cutoff and offset. 
# Use an outer area to avoid a boundary effect. The outer area can be less fine 
# than the inner area to reduce computing time.

# Come up with an initial guess for the range:
#  At what distance is the correlation smaller than 0.1?
#  20 km? See also the histogram if the distances again.
#  If we choose 20 km then the correlation will affect around 10%
#  of the combinations. The value of 20 km also complies (vaguely) with
#  the sample variogram of the random intercepts (maybe 30 would have been better).

# So..let's use 20 km for the range?
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
range(catchwb6_train$avgwt_g)
# P(Range < 5 km ) = 0.05
# P(sigma > 2) = 0.05
# SB = exp(u_i)
# some u_i have to be as large as 5.24 to cover 190
# If u_i ~ N(0, sigma_u^2) then it is unlikley that sigma_u > 1.5
#P(sigma > 1.5) = 0.05
M1 <- glm(avgwt_g ~ 1, data = catchwb6_train)
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
Xcov <- model.matrix(~ soaktimec + year + temp_Cc + turb_FNUc + distance_to_southc + distance_to_shorec,
                     data = catchwb6_train)
Xcov <- as.data.frame(Xcov)

Stack <- inla.stack(
  tag = "Fit",
  data = list(avgwt_g = catchwb6_train$avgwt_g),  
  A = list(1, 1, 1, A),                 
  effects = list(                 
    Intercept = rep(1,N),
    Xcov  = Xcov[,-1],               
    site = catchwb6_train$site, # random effect
    w  = w.index))  # SRF


f14 <- avgwt_g ~ -1 + Intercept + soaktimec  + year + temp_Cc + turb_FNUc + distance_to_southc + 
  distance_to_shorec +
  f(site, model="iid") +
  f(w, model = spde)


# 7. Run the spatial model in INLA.
I14 <- inla(f14,
            family = "gamma", 
            data=inla.stack.data(Stack),
            control.compute = list(dic = TRUE,  waic=TRUE, config = TRUE),
            control.family=list(link='log'),
            control.predictor = list(A = inla.stack.A(Stack)))

summary(I14)


# Posterior mean values and 95% CI for the regression parameters:
Beta14 <- I14$summary.fixed[,c("mean", "0.025quant", "0.975quant")] 
print(Beta14, digits = 2)  
# year, temp, distance to south do not contain 0 so are important

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

w.pm        <- I14$summary.random$w$mean  #mu.srf # exp converts from the log link (value of 1 is = exp(1), but now all exp())
w.proj      <- inla.mesh.projector(mesh) 
w.pm100_100 <- inla.mesh.project(w.proj, w.pm)
long        <- as.data.frame(w.proj$lattice$loc)
colnames(long) <- c("EastingKM", "NorthingKM")
long <- long %>%
  mutate(w.pm = as.vector(w.pm100_100))

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
abline(v=c(5, 9.25, 22), col='red')

abline(h=0.6, v=9.25)



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



# Here are the fitted values and Pearson residuals
N     <- nrow(catchwb6_train)
mu14   <- I14$summary.fitted.values[1:N,"mean"] 
r     <- I14$summary.hyperpar["Precision parameter for the Gamma observations", "mean"]
VarY14 <- mu14^2 / r
E14    <- (catchwb6_train$avgwt_g - mu14) / sqrt(VarY14)



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
     y = catchwb6_train$avgwt_g,
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


Xcov <- model.matrix(~ soaktimec + year + distance_to_southc,
                     data = catchwb6_train)
Xcov2 <- as.data.frame(Xcov)
colnames(Xcov2)

N    <- nrow(catchwb6_train)

Stack <- inla.stack(
  tag = "Fit",
  data = list(y = catchwb6_train$avgwt_g),
  A = list(1,1,1,A),
  effects = list(
    Intercept = rep(1,N),
    Xcov2 = Xcov2[,-1], 
    site = catchwb6_train$site,
    w  = w.index))


f2 <- y ~ -1 + Intercept +
  soaktimec  + year + distance_to_southc +
  f(site, model='iid') +
  f(w, model = spde)


I2 <- inla(f2,
           family = "gamma",
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
Xp <- model.matrix(~ soaktimec + year + distance_to_southc, data = catchwb6_test)

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
  soaktimec  + year + distance_to_southc +
  f(site, model="iid") +
  f(w, model = spde)

# Note that we use a different stack! The rest is the same
I5.pred <- inla(f5,
                family = "gamma", 
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
  mutate(set = ifelse(!is.na(avgwt_g), "Training", "Testing"))
catchwb6_pred$set <- as.factor(catchwb6_pred$set)


catchwb6_pred <- catchwb6_pred %>% 
  mutate(avgwt_g=ifelse(!is.na(avgwt_g), avgwt_g, mean))

head(catchwb6_pred)

p <- ggplot()
p <- p + geom_smooth(data = catchwb6_test2, 
                     aes(x = soaktime, 
                         y = mean), alpha=0.5)
p <- p + geom_point(data = catchwb6_pred, 
                    aes(y = avgwt_g, x = soaktime, col=set),
                    size = 3)
p <- p + xlab("Soak time (h)") + ylab("Mean weight (g)") + labs ( col = "Dataset")
p <- p + theme(text = element_text(size=15), legend.position = "bottom") 
# p <- p + geom_ribbon(data = catchwb6_test3a,
#                      aes(x = soaktime2,
#                          ymax = upper,
#                          ymin = lower),
#                      alpha = 0.2)

p

ggplot(data=catchwb6_test2) +
  geom_point(mapping=aes(x=resp_obsB, y=mu, colour=fyear), size=3) +
  xlab("Observed") +
  ylab("Predicted") +
  geom_abline() +
  xlim(c(0,1200)) +
  ylim(c(0,1200)) +
  labs(colour="Year")





# *Gamma GLM d2S temp ind srf --------
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
# Step 2 of making a mesh: The primary tool to control the shape of the 
# triangles is max.edge. Other useful arguments are the cutoff and offset. 
# Use an outer area to avoid a boundary effect. The outer area can be less fine 
# than the inner area to reduce computing time.

# Come up with an initial guess for the range:
#  At what distance is the correlation smaller than 0.1?
#  20 km? See also the histogram if the distances again.
#  If we choose 20 km then the correlation will affect around 10%
#  of the combinations. The value of 20 km also complies (vaguely) with
#  the sample variogram of the random intercepts (maybe 30 would have been better).

# So..let's use 20 km for the range?
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
range(catchwb6_train$avgwt_g)
# P(Range < 5 km ) = 0.05
# P(sigma > 2) = 0.05
# SB = exp(u_i)
# some u_i have to be as large as 5.24 to cover 190
# If u_i ~ N(0, sigma_u^2) then it is unlikley that sigma_u > 1.5
#P(sigma > 1.5) = 0.05
M1 <- glm(avgwt_g ~ 1, data = catchwb6_train)
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
Xcov <- model.matrix(~ soaktimec + year + temp_Cc + distance_to_southc,
                     data = catchwb6_train)
Xcov <- as.data.frame(Xcov)

Stack <- inla.stack(
  tag = "Fit",
  data = list(avgwt_g = catchwb6_train$avgwt_g),  
  A = list(1, 1, 1, A),                 
  effects = list(                 
    Intercept = rep(1,N),
    Xcov  = Xcov[,-1],               
    site = catchwb6_train$site, # random effect
    w  = w.index))  # SRF


f14 <- avgwt_g ~ -1 + Intercept + soaktimec  + year + temp_Cc + distance_to_southc +
  f(site, model="iid") +
  f(w, model = spde)


# 7. Run the spatial model in INLA.
I14 <- inla(f14,
            family = "gamma", 
            data=inla.stack.data(Stack),
            control.compute = list(dic = TRUE,  waic=TRUE, config = TRUE),
            control.family=list(link='log'),
            control.predictor = list(A = inla.stack.A(Stack)))

summary(I14)


# Posterior mean values and 95% CI for the regression parameters:
Beta14 <- I14$summary.fixed[,c("mean", "0.025quant", "0.975quant")] 
print(Beta14, digits = 2)  
## temp no longer impt (has 0 in inv)
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

w.pm        <- I14$summary.random$w$mean  #mu.srf # exp converts from the log link (value of 1 is = exp(1), but now all exp())
w.proj      <- inla.mesh.projector(mesh) 
w.pm100_100 <- inla.mesh.project(w.proj, w.pm)
long        <- as.data.frame(w.proj$lattice$loc)
colnames(long) <- c("EastingKM", "NorthingKM")
long <- long %>%
  mutate(w.pm = as.vector(w.pm100_100))

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
abline(v=c(5, 9.25, 22), col='red')

abline(h=0.6, v=9.25)



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



# Here are the fitted values and Pearson residuals
N     <- nrow(catchwb6_train)
mu14   <- I14$summary.fitted.values[1:N,"mean"] 
r     <- I14$summary.hyperpar["Precision parameter for the Gamma observations", "mean"]
VarY14 <- mu14^2 / r
E14    <- (catchwb6_train$avgwt_g - mu14) / sqrt(VarY14)



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
     y = catchwb6_train$avgwt_g,
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


Xcov <- model.matrix(~ soaktimec + year + distance_to_southc,
                     data = catchwb6_train)
Xcov2 <- as.data.frame(Xcov)
colnames(Xcov2)

N    <- nrow(catchwb6_train)

Stack <- inla.stack(
  tag = "Fit",
  data = list(y = catchwb6_train$avgwt_g),
  A = list(1,1,1,A),
  effects = list(
    Intercept = rep(1,N),
    Xcov2 = Xcov2[,-1], 
    site = catchwb6_train$site,
    w  = w.index))


f2 <- y ~ -1 + Intercept +
  soaktimec  + year + distance_to_southc +
  f(site, model='iid') +
  f(w, model = spde)


I2 <- inla(f2,
           family = "gamma",
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
Xp <- model.matrix(~ soaktimec + year + distance_to_southc, data = catchwb6_test)

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
  soaktimec  + year + distance_to_southc +
  f(site, model="iid") +
  f(w, model = spde)

# Note that we use a different stack! The rest is the same
I5.pred <- inla(f5,
                family = "gamma", 
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
  mutate(set = ifelse(!is.na(avgwt_g), "Training", "Testing"))
catchwb6_pred$set <- as.factor(catchwb6_pred$set)


catchwb6_pred <- catchwb6_pred %>% 
  mutate(avgwt_g=ifelse(!is.na(avgwt_g), avgwt_g, mean))

head(catchwb6_pred)

p <- ggplot()
p <- p + geom_smooth(data = catchwb6_test2, 
                     aes(x = soaktime, 
                         y = mean), alpha=0.5)
p <- p + geom_point(data = catchwb6_pred, 
                    aes(y = avgwt_g, x = soaktime, col=set),
                    size = 3)
p <- p + xlab("Soak time (h)") + ylab("Mean weight (g)") + labs ( col = "Dataset")
p <- p + theme(text = element_text(size=15), legend.position = "bottom") 
# p <- p + geom_ribbon(data = catchwb6_test3a,
#                      aes(x = soaktime2,
#                          ymax = upper,
#                          ymin = lower),
#                      alpha = 0.2)

p

ggplot(data=catchwb6_test2) +
  geom_point(mapping=aes(x=resp_obsB, y=mu, colour=fyear), size=3) +
  xlab("Observed") +
  ylab("Predicted") +
  geom_abline() +
  xlim(c(0,1200)) +
  ylim(c(0,1200)) +
  labs(colour="Year")

# Gamma GLM d2S ind srf --------
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
# Step 2 of making a mesh: The primary tool to control the shape of the 
# triangles is max.edge. Other useful arguments are the cutoff and offset. 
# Use an outer area to avoid a boundary effect. The outer area can be less fine 
# than the inner area to reduce computing time.

# Come up with an initial guess for the range:
#  At what distance is the correlation smaller than 0.1?
#  20 km? See also the histogram if the distances again.
#  If we choose 20 km then the correlation will affect around 10%
#  of the combinations. The value of 20 km also complies (vaguely) with
#  the sample variogram of the random intercepts (maybe 30 would have been better).

# So..let's use 20 km for the range?
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
range(catchwb6_train$avgwt_g)
# P(Range < 5 km ) = 0.05
# P(sigma > 2) = 0.05
# SB = exp(u_i)
# some u_i have to be as large as 5.24 to cover 190
# If u_i ~ N(0, sigma_u^2) then it is unlikley that sigma_u > 1.5
#P(sigma > 1.5) = 0.05
M1 <- glm(avgwt_g ~ 1, data = catchwb6_train)
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
Xcov <- model.matrix(~ soaktimec + year +  distance_to_southc,
                     data = catchwb6_train)
Xcov <- as.data.frame(Xcov)

Stack <- inla.stack(
  tag = "Fit",
  data = list(avgwt_g = catchwb6_train$avgwt_g),  
  A = list(1, 1, 1, A),                 
  effects = list(                 
    Intercept = rep(1,N),
    Xcov  = Xcov[,-1],               
    site = catchwb6_train$site, # random effect
    w  = w.index))  # SRF


f14 <- avgwt_g ~ -1 + Intercept + soaktimec  + year +  distance_to_southc +
  f(site, model="iid") +
  f(w, model = spde)


# 7. Run the spatial model in INLA.
I14 <- inla(f14,
            family = "gamma", 
            data=inla.stack.data(Stack),
            control.compute = list(dic = TRUE,  waic=TRUE, config = TRUE),
            control.family=list(link='log'),
            control.predictor = list(A = inla.stack.A(Stack)))

summary(I14)


# Posterior mean values and 95% CI for the regression parameters:
Beta14 <- I14$summary.fixed[,c("mean", "0.025quant", "0.975quant")] 
print(Beta14, digits = 2)  
## temp no longer impt (has 0 in inv)
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

w.pm        <- I14$summary.random$w$mean  #mu.srf # exp converts from the log link (value of 1 is = exp(1), but now all exp())
w.proj      <- inla.mesh.projector(mesh) 
w.pm100_100 <- inla.mesh.project(w.proj, w.pm)
long        <- as.data.frame(w.proj$lattice$loc)
colnames(long) <- c("EastingKM", "NorthingKM")
long <- long %>%
  mutate(w.pm = as.vector(w.pm100_100))

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
abline(v=c(5, 9.25, 22), col='red')

abline(h=0.6, v=9.25)



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



# Here are the fitted values and Pearson residuals
N     <- nrow(catchwb6_train)
mu14   <- I14$summary.fitted.values[1:N,"mean"] 
r     <- I14$summary.hyperpar["Precision parameter for the Gamma observations", "mean"]
VarY14 <- mu14^2 / r
E14    <- (catchwb6_train$avgwt_g - mu14) / sqrt(VarY14)



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
     y = catchwb6_train$avgwt_g,
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


Xcov <- model.matrix(~ soaktimec + year + distance_to_southc,
                     data = catchwb6_train)
Xcov2 <- as.data.frame(Xcov)
colnames(Xcov2)

N    <- nrow(catchwb6_train)

Stack <- inla.stack(
  tag = "Fit",
  data = list(y = catchwb6_train$avgwt_g),
  A = list(1,1,1,A),
  effects = list(
    Intercept = rep(1,N),
    Xcov2 = Xcov2[,-1], 
    site = catchwb6_train$site,
    w  = w.index))


f2 <- y ~ -1 + Intercept +
  soaktimec  + year + distance_to_southc +
  f(site, model='iid') +
  f(w, model = spde)


I2 <- inla(f2,
           family = "gamma",
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
Xp <- model.matrix(~ soaktimec + year + distance_to_southc, data = catchwb6_test)

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
  soaktimec  + year + distance_to_southc +
  f(site, model="iid") +
  f(w, model = spde)

# Note that we use a different stack! The rest is the same
I5.pred <- inla(f5,
                family = "gamma", 
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
  mutate(set = ifelse(!is.na(avgwt_g), "Training", "Testing"))
catchwb6_pred$set <- as.factor(catchwb6_pred$set)


catchwb6_pred <- catchwb6_pred %>% 
  mutate(avgwt_g=ifelse(!is.na(avgwt_g), avgwt_g, mean))

head(catchwb6_pred)

p <- ggplot()
p <- p + geom_smooth(data = catchwb6_test2, 
                     aes(x = soaktime, 
                         y = mean), alpha=0.5)
p <- p + geom_point(data = catchwb6_pred, 
                    aes(y = avgwt_g, x = soaktime, col=set),
                    size = 3)
p <- p + xlab("Soak time (h)") + ylab("Mean weight (g)") + labs ( col = "Dataset")
p <- p + theme(text = element_text(size=15), legend.position = "bottom") 
# p <- p + geom_ribbon(data = catchwb6_test3a,
#                      aes(x = soaktime2,
#                          ymax = upper,
#                          ymin = lower),
#                      alpha = 0.2)

p

ggplot(data=catchwb6_test2) +
  geom_point(mapping=aes(x=resp_obsB, y=mu, colour=fyear), size=3) +
  xlab("Observed") +
  ylab("Predicted") +
  geom_abline() +
  xlim(c(0,1200)) +
  ylim(c(0,1200)) +
  labs(colour="Year")


# *Gamma GLM no site ind srf --------
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
# Step 2 of making a mesh: The primary tool to control the shape of the 
# triangles is max.edge. Other useful arguments are the cutoff and offset. 
# Use an outer area to avoid a boundary effect. The outer area can be less fine 
# than the inner area to reduce computing time.

# Come up with an initial guess for the range:
#  At what distance is the correlation smaller than 0.1?
#  20 km? See also the histogram if the distances again.
#  If we choose 20 km then the correlation will affect around 10%
#  of the combinations. The value of 20 km also complies (vaguely) with
#  the sample variogram of the random intercepts (maybe 30 would have been better).

# So..let's use 20 km for the range?
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
range(catchwb6_train$avgwt_g)
# P(Range < 5 km ) = 0.05
# P(sigma > 2) = 0.05
# SB = exp(u_i)
# some u_i have to be as large as 5.24 to cover 190
# If u_i ~ N(0, sigma_u^2) then it is unlikley that sigma_u > 1.5
#P(sigma > 1.5) = 0.05
M1 <- glm(avgwt_g ~ 1, data = catchwb6_train)
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
Xcov <- model.matrix(~ soaktimec + year ,
                     data = catchwb6_train)
Xcov <- as.data.frame(Xcov)

Stack <- inla.stack(
  tag = "Fit",
  data = list(avgwt_g = catchwb6_train$avgwt_g),  
  A = list(1,  1, A),                 
  effects = list(                 
    Intercept = rep(1,N),
    Xcov  = Xcov[,-1], # random effect
    w  = w.index))  # SRF


f14 <- avgwt_g ~ -1 + Intercept + soaktimec  + year +
  f(w, model = spde)


# 7. Run the spatial model in INLA.
I14 <- inla(f14,
            family = "gamma", 
            data=inla.stack.data(Stack),
            control.compute = list(dic = TRUE,  waic=TRUE, config = TRUE),
            control.family=list(link='log'),
            control.predictor = list(A = inla.stack.A(Stack)))

summary(I14)


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
abline(v=c(5, 9.25, 22), col='red')

abline(h=0.6, v=9.25)


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



# Here are the fitted values and Pearson residuals
N     <- nrow(catchwb6_train)
mu14   <- I14$summary.fitted.values[1:N,"mean"] 
r     <- I14$summary.hyperpar["Precision parameter for the Gamma observations", "mean"]
VarY14 <- mu14^2 / r
E14    <- (catchwb6_train$avgwt_g - mu14) / sqrt(VarY14)



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
     y = catchwb6_train$avgwt_g,
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
  data = list(y = catchwb6_train$avgwt_g),
  A = list(1,1,A),
  effects = list(
    Intercept = rep(1,N),
    Xcov2 = Xcov2[,-1],
    w  = w.index))


f2 <- y ~ -1 + Intercept +
  soaktimec  + year +
  f(w, model = spde)


I2 <- inla(f2,
           family = "gamma",
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
  A = list(1,  1), # a 1 for each effect                 
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
  soaktimec  + year +
  f(w, model = spde)

# Note that we use a different stack! The rest is the same
I5.pred <- inla(f5,
                family = "gamma", 
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
                    aes(y = avgwt_g, x = soaktime),
                    size = 1, colour="black")
p <- p + xlab("soak time (h)") + ylab("avg wt LKWF")
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
  geom_point(mapping=aes(x=resp_obsB, y=mu), size=3) +
  xlab("Observed") +
  ylab("Predicted") +
  geom_abline() +
  xlim(c(0,800)) +
  ylim(c(0,800)) 


# *Gamma GLM no site ind srf --------
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
# Step 2 of making a mesh: The primary tool to control the shape of the 
# triangles is max.edge. Other useful arguments are the cutoff and offset. 
# Use an outer area to avoid a boundary effect. The outer area can be less fine 
# than the inner area to reduce computing time.

# Come up with an initial guess for the range:
#  At what distance is the correlation smaller than 0.1?
#  20 km? See also the histogram if the distances again.
#  If we choose 20 km then the correlation will affect around 10%
#  of the combinations. The value of 20 km also complies (vaguely) with
#  the sample variogram of the random intercepts (maybe 30 would have been better).

# So..let's use 20 km for the range?
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
range(catchwb6_train$avgwt_g)
# P(Range < 5 km ) = 0.05
# P(sigma > 2) = 0.05
# SB = exp(u_i)
# some u_i have to be as large as 5.24 to cover 190
# If u_i ~ N(0, sigma_u^2) then it is unlikley that sigma_u > 1.5
#P(sigma > 1.5) = 0.05
M1 <- glm(avgwt_g ~ 1, data = catchwb6_train)
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
Xcov <- model.matrix(~ soaktimec + year ,
                     data = catchwb6_train)
Xcov <- as.data.frame(Xcov)

Stack <- inla.stack(
  tag = "Fit",
  data = list(avgwt_g = catchwb6_train$avgwt_g),  
  A = list(1,  1, A),                 
  effects = list(                 
    Intercept = rep(1,N),
    Xcov  = Xcov[,-1], # random effect
    w  = w.index))  # SRF


f14 <- avgwt_g ~ -1 + Intercept + soaktimec  + year +
  f(w, model = spde)


# 7. Run the spatial model in INLA.
I14 <- inla(f14,
            family = "gamma", 
            data=inla.stack.data(Stack),
            control.compute = list(dic = TRUE,  waic=TRUE, config = TRUE),
            control.family=list(link='log'),
            control.predictor = list(A = inla.stack.A(Stack)))

summary(I14)


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
abline(v=c(5, 9.25, 22), col='red')

abline(h=0.6, v=9.25)


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



# Here are the fitted values and Pearson residuals
N     <- nrow(catchwb6_train)
mu14   <- I14$summary.fitted.values[1:N,"mean"] 
r     <- I14$summary.hyperpar["Precision parameter for the Gamma observations", "mean"]
VarY14 <- mu14^2 / r
E14    <- (catchwb6_train$avgwt_g - mu14) / sqrt(VarY14)



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
     y = catchwb6_train$avgwt_g,
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
  data = list(y = catchwb6_train$avgwt_g),
  A = list(1,1,A),
  effects = list(
    Intercept = rep(1,N),
    Xcov2 = Xcov2[,-1],
    w  = w.index))


f2 <- y ~ -1 + Intercept +
  soaktimec  + year +
  f(w, model = spde)


I2 <- inla(f2,
           family = "gamma",
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
  A = list(1,  1), # a 1 for each effect                 
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
  soaktimec  + year +
  f(w, model = spde)

# Note that we use a different stack! The rest is the same
I5.pred <- inla(f5,
                family = "gamma", 
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
                    aes(y = avgwt_g, x = soaktime),
                    size = 1, colour="black")
p <- p + xlab("soak time (h)") + ylab("avg wt LKWF")
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
  geom_point(mapping=aes(x=resp_obsB, y=mu), size=3) +
  xlab("Observed") +
  ylab("Predicted") +
  geom_abline() +
  xlim(c(0,800)) +
  ylim(c(0,800)) 


#