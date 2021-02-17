# GLM Modelling
install.packages("DHARMa")
library(tidyverse)
library(magrittr)
library(readxl)
library(mgcv)
library(MASS)
library(SuppDists)
library(DHARMa)
library(boot)

dataset1 <- read.csv("Data/train1.csv") 
dataset2_raw <- read.csv("Data/train2.csv")
dataset2 <- dataset2_raw[-c(3,4)] # Get rid of lat and long

# Merge datasets together
merged <- rbind(dataset1, dataset2)

# Remove rows with missing values
model_data <- na.omit(merged)
attach(model_data)

# Create a running average vector of all variables - takes average of every 5 values
cphl_avg5 <- rep(NA, 308835)
depth_avg5 <- rep(NA, 308835)
temp_avg5 <- rep(NA, 308835)
cndc_avg5 <- rep(NA, 308835)
psal_avg5 <- rep(NA, 308835)
dox1_avg5 <- rep(NA, 308835)
cdom_avg5 <- rep(NA, 308835)
vbsc_avg5 <- rep(NA, 308835)
irrad443_avg5 <- rep(NA, 308835)
irrad490_avg5 <- rep(NA, 308835)
irrad555_avg5 <- rep(NA, 308835)
irrad670_avg5 <- rep(NA, 308835)
j = 1

# Take average of every 5 data points in train set
for (i in seq(from=1, to=1544171, by=5)) {
  cphl_avg5[j] = mean(cphl[i:(i+4)])
  depth_avg5[j] = mean(depth[i:(i+4)])
  temp_avg5[j] = mean(temp[i:(i+4)])
  cndc_avg5[j] = mean(cndc[i:(i+4)])
  psal_avg5[j] = mean(psal[i:(i+4)])
  dox1_avg5[j] = mean(dox1[i:(i+4)])
  cdom_avg5[j] = mean(cdom[i:(i+4)])
  vbsc_avg5[j] = mean(vbsc[i:(i+4)])
  irrad443_avg5[j] = mean(irrad443[i:(i+4)])
  irrad490_avg5[j] = mean(irrad490[i:(i+4)])
  irrad555_avg5[j] = mean(irrad555[i:(i+4)])
  irrad670_avg5[j] = mean(irrad670[i:(i+4)])
  j = j + 1
}

# Combine all averaged columns into 1 dataframe
model_data_avg5 <- data.frame(cbind(cphl_avg5,depth_avg5,temp_avg5,cndc_avg5,psal_avg5,dox1_avg5,
                                    cdom_avg5,vbsc_avg5,irrad443_avg5,irrad490_avg5,
                                    irrad555_avg5,irrad670_avg5))

# Only take rows with positive cphl_avg values
model_data_avg5_pos <- model_data_avg5[model_data_avg5$cphl_avg5 > 0.01,]

#Rename columns
names(model_data_avg5_pos)[1] <- "cphl"
names(model_data_avg5_pos)[2] <- "depth"
names(model_data_avg5_pos)[3] <- "temp"
names(model_data_avg5_pos)[4] <- "cndc"
names(model_data_avg5_pos)[5] <- "psal"
names(model_data_avg5_pos)[6] <- "dox1"
names(model_data_avg5_pos)[7] <- "cdom"
names(model_data_avg5_pos)[8] <- "vbsc"
names(model_data_avg5_pos)[9] <- "irrad443"
names(model_data_avg5_pos)[10] <- "irrad490"
names(model_data_avg5_pos)[11] <- "irrad555"
names(model_data_avg5_pos)[12] <- "irrad670"

# Method of Moments Estimation of k = shape and theta = scale parameters for Gamma Distribution
sample_mean <- mean(model_data_avg5_pos$cphl)
sample_var <- var(model_data_avg5_pos$cphl)

k <- (sample_mean)^2/sample_var
theta <- sample_var/sample_mean

# Histogram of cphl-average5 combined with density curve for inverse gaussian and gamma
par(mfrow=c(1,1))
hist(cphl, probability = TRUE, xlab = "Chlorophyll-a Concentration", main = "Histogram of Response Variable Chlorophyll-a")
# Method of Moments estimate of gamma distribution
curve(dgamma(x,shape = k, scale = theta), add = TRUE, col="red")

# Manual AIC (Manual process was required by the UNSW course, to show every step)
attach(model_data_avg5_pos)

fit1 <- glm(cphl~., family = Gamma(link = "inverse"), data = model_data_avg5_pos,
            start = rep(1,12))
fit2 <- glm(cphl~1, start = rep(1,1), family = Gamma(link = "inverse"), data = model_data_avg5_pos)

fit3 <- glm(cphl~depth, start = rep(1,2), family = Gamma(link = "inverse"), data = model_data_avg5_pos)
fit3.1 <- glm(cphl~depth+temp, start = rep(1,3), family = Gamma(link = "inverse"), data = model_data_avg5_pos)
fit3.2 <- glm(cphl~depth+psal, start = rep(1,3), family = Gamma(link = "inverse"), data = model_data_avg5_pos)
fit3.3 <- glm(cphl~depth+dox1, start = rep(1,3), family = Gamma(link = "inverse"), data = model_data_avg5_pos)
fit3.4 <- glm(cphl~depth+cdom, start = rep(1,3), family = Gamma(link = "inverse"), data = model_data_avg5_pos)
fit3.5 <- glm(cphl~depth+vbsc, start = rep(1,3), family = Gamma(link = "inverse"), data = model_data_avg5_pos)
fit3.6 <- glm(cphl~depth+irrad443, start = rep(1,3), family = Gamma(link = "inverse"), data = model_data_avg5_pos)
fit3.7 <- glm(cphl~depth+irrad490, start = rep(1,3), family = Gamma(link = "inverse"), data = model_data_avg5_pos)
fit3.8 <- glm(cphl~depth+irrad555, start = rep(1,3), family = Gamma(link = "inverse"), data = model_data_avg5_pos)
fit3.9 <- glm(cphl~depth+irrad670, start = rep(1,3), family = Gamma(link = "inverse"), data = model_data_avg5_pos)

aic3.1 <- AIC(fit3.1)
aic3.2 <- AIC(fit3.2)
aic3.3 <- AIC(fit3.3)
aic3.4 <- AIC(fit3.4)
aic3.5 <- AIC(fit3.5)
aic3.6 <- AIC(fit3.6) ### WINNER
aic3.7 <- AIC(fit3.7)
aic3.8 <- AIC(fit3.8)
aic3.9 <- AIC(fit3.9)

min3 <- min(c(aic3.1,aic3.2,aic3.3,aic3.4,aic3.5,aic3.6,aic3.7,aic3.8,aic3.9))
(min3 == c(aic3.1,aic3.2,aic3.3,aic3.4,aic3.5,aic3.6,aic3.7,aic3.8,aic3.9))

fit4 <- glm(cphl~depth+irrad443, start = rep(1,3), family = Gamma(link = "inverse"), data = model_data_avg5_pos)
fit4.1 <- glm(cphl~depth+irrad443+temp, start = rep(1,4), family = Gamma(link = "inverse"), data = model_data_avg5_pos)
fit4.2 <- glm(cphl~depth+irrad443+psal, start = rep(1,4), family = Gamma(link = "inverse"), data = model_data_avg5_pos)
fit4.3 <- glm(cphl~depth+irrad443+dox1, start = rep(1,4), family = Gamma(link = "inverse"), data = model_data_avg5_pos)
fit4.4 <- glm(cphl~depth+irrad443+cdom, start = rep(1,4), family = Gamma(link = "inverse"), data = model_data_avg5_pos)
fit4.5 <- glm(cphl~depth+irrad443+vbsc, start = rep(1,4), family = Gamma(link = "inverse"), data = model_data_avg5_pos)
fit4.6 <- glm(cphl~depth+irrad443+irrad490, start = rep(1,4), family = Gamma(link = "inverse"), data = model_data_avg5_pos)
fit4.7 <- glm(cphl~depth+irrad443+irrad555, start = rep(1,4), family = Gamma(link = "inverse"), data = model_data_avg5_pos)
fit4.8 <- glm(cphl~depth+irrad443+irrad670, start = rep(1,4), family = Gamma(link = "inverse"), data = model_data_avg5_pos)

aic4.1 <- AIC(fit4.1)
aic4.2 <- AIC(fit4.2)
aic4.3 <- AIC(fit4.3)
aic4.4 <- AIC(fit4.4) ### WINNER
aic4.5 <- AIC(fit4.5)
aic4.6 <- AIC(fit4.6)
aic4.7 <- AIC(fit4.7)
aic4.8 <- AIC(fit4.8)

min4 <- min(c(aic4.1,aic4.2,aic4.3,aic4.4,aic4.5,aic4.6,aic4.7,aic4.8))
(min4 == c(aic4.1,aic4.2,aic4.3,aic4.4,aic4.5,aic4.6,aic4.7,aic4.8))

fit5 <- glm(cphl~depth+irrad443+cdom, start = rep(1,4), family = Gamma(link = "inverse"), data = model_data_avg5_pos)
fit5.1 <- glm(cphl~depth+irrad443+cdom+temp, start = rep(1,5), family = Gamma(link = "inverse"), data = model_data_avg5_pos)
fit5.2 <- glm(cphl~depth+irrad443+cdom+psal, start = rep(1,5), family = Gamma(link = "inverse"), data = model_data_avg5_pos)
fit5.3 <- glm(cphl~depth+irrad443+cdom+dox1, start = rep(1,5), family = Gamma(link = "inverse"), data = model_data_avg5_pos)
fit5.4 <- glm(cphl~depth+irrad443+cdom+vbsc, start = rep(1,5), family = Gamma(link = "inverse"), data = model_data_avg5_pos)
fit5.5 <- glm(cphl~depth+irrad443+cdom+irrad490, start = rep(1,5), family = Gamma(link = "inverse"), data = model_data_avg5_pos)
fit5.6 <- glm(cphl~depth+irrad443+cdom+irrad555, start = rep(1,5), family = Gamma(link = "inverse"), data = model_data_avg5_pos)
fit5.7 <- glm(cphl~depth+irrad443+cdom+irrad670, start = rep(1,5), family = Gamma(link = "inverse"), data = model_data_avg5_pos)

aic5.1 <- AIC(fit5.1) ### WINNER
aic5.2 <- AIC(fit5.2)
aic5.3 <- AIC(fit5.3)
aic5.4 <- AIC(fit5.4)
aic5.5 <- AIC(fit5.5)
aic5.6 <- AIC(fit5.6)
aic5.7 <- AIC(fit5.7)

min5 <- min(c(aic5.1,aic5.2,aic5.3,aic5.4,aic5.5,aic5.6,aic5.7))
(min5 == c(aic5.1,aic5.2,aic5.3,aic5.4,aic5.5,aic5.6,aic5.7))

fit6 <- glm(cphl~depth+irrad443+cdom+temp, start = rep(1,5), family = Gamma(link = "inverse"), data = model_data_avg5_pos)
fit6.1 <- glm(cphl~depth+irrad443+cdom+temp+psal, start = rep(1,6), family = Gamma(link = "inverse"), data = model_data_avg5_pos)
fit6.2 <- glm(cphl~depth+irrad443+cdom+temp+dox1, start = rep(1,6), family = Gamma(link = "inverse"), data = model_data_avg5_pos)
fit6.3 <- glm(cphl~depth+irrad443+cdom+temp+vbsc, start = rep(1,6), family = Gamma(link = "inverse"), data = model_data_avg5_pos)
fit6.4 <- glm(cphl~depth+irrad443+cdom+temp+irrad490, start = rep(1,6), family = Gamma(link = "inverse"), data = model_data_avg5_pos)
fit6.5 <- glm(cphl~depth+irrad443+cdom+temp+irrad555, start = rep(1,6), family = Gamma(link = "inverse"), data = model_data_avg5_pos)
fit6.6 <- glm(cphl~depth+irrad443+cdom+temp+irrad670, start = rep(1,6), family = Gamma(link = "inverse"), data = model_data_avg5_pos)

aic6.1 <- AIC(fit6.1) ### WINNER
aic6.2 <- AIC(fit6.2)
aic6.3 <- AIC(fit6.3)
aic6.4 <- AIC(fit6.4)
aic6.5 <- AIC(fit6.5)
aic6.6 <- AIC(fit6.6)

min6 <- min(c(aic6.1,aic6.2,aic6.3,aic6.4,aic6.5,aic6.6))
(min6 == c(aic6.1,aic6.2,aic6.3,aic6.4,aic6.5,aic6.6))

fit7 <- glm(cphl~depth+irrad443+cdom+temp+psal, start = rep(1,6), family = Gamma(link = "inverse"), data = model_data_avg5_pos)
fit7.1 <- glm(cphl~depth+irrad443+cdom+temp+psal+dox1, start = rep(1,7), family = Gamma(link = "inverse"), data = model_data_avg5_pos)
fit7.2 <- glm(cphl~depth+irrad443+cdom+temp+psal+vbsc, start = rep(1,7), family = Gamma(link = "inverse"), data = model_data_avg5_pos)
fit7.3 <- glm(cphl~depth+irrad443+cdom+temp+psal+irrad490, start = rep(1,7), family = Gamma(link = "inverse"), data = model_data_avg5_pos)
fit7.4 <- glm(cphl~depth+irrad443+cdom+temp+psal+irrad555, start = rep(1,7), family = Gamma(link = "inverse"), data = model_data_avg5_pos)
fit7.5 <- glm(cphl~depth+irrad443+cdom+temp+psal+irrad670, start = rep(1,7), family = Gamma(link = "inverse"), data = model_data_avg5_pos)

aic7.1 <- AIC(fit7.1) 
aic7.2 <- AIC(fit7.2)
aic7.3 <- AIC(fit7.3) ### WINNER
aic7.4 <- AIC(fit7.4)
aic7.5 <- AIC(fit7.5)

min7 <- min(c(aic7.1,aic7.2,aic7.3,aic7.4,aic7.5))
(min7 == c(aic7.1,aic7.2,aic7.3,aic7.4,aic7.5))

fit8 <- glm(cphl~depth+irrad443+cdom+temp+psal+irrad490, start = rep(1,7), family = Gamma(link = "inverse"), data = model_data_avg5_pos)
fit8.1 <- glm(cphl~depth+irrad443+cdom+temp+psal+irrad490+dox1, start = rep(1,8), family = Gamma(link = "inverse"), data = model_data_avg5_pos)
fit8.2 <- glm(cphl~depth+irrad443+cdom+temp+psal+irrad490+vbsc, start = rep(1,8), family = Gamma(link = "inverse"), data = model_data_avg5_pos)
fit8.3 <- glm(cphl~depth+irrad443+cdom+temp+psal+irrad490+irrad555, start = rep(1,8), family = Gamma(link = "inverse"), data = model_data_avg5_pos)
fit8.4 <- glm(cphl~depth+irrad443+cdom+temp+psal+irrad490+irrad670, start = rep(1,8), family = Gamma(link = "inverse"), data = model_data_avg5_pos)

aic8.1 <- AIC(fit8.1) 
aic8.2 <- AIC(fit8.2)
aic8.3 <- AIC(fit8.3) ### WINNER
aic8.4 <- AIC(fit8.4)

min8 <- min(c(aic8.1,aic8.2,aic8.3,aic8.4))
(min8 == c(aic8.1,aic8.2,aic8.3,aic8.4))

fit9 <- glm(cphl~depth+irrad443+cdom+temp+psal+irrad490+irrad555, start = rep(1,8), family = Gamma(link = "inverse"), data = model_data_avg5_pos)
fit9.1 <- glm(cphl~depth+irrad443+cdom+temp+psal+irrad555+irrad490+dox1, start = rep(1,9), family = Gamma(link = "inverse"), data = model_data_avg5_pos)
fit9.2 <- glm(cphl~depth+irrad443+cdom+temp+psal+irrad490+irrad555+vbsc, start = rep(1,9), family = Gamma(link = "inverse"), data = model_data_avg5_pos)
fit9.3 <- glm(cphl~depth+irrad443+cdom+temp+psal+irrad490+irrad555+irrad670, start = rep(1,9), family = Gamma(link = "inverse"), data = model_data_avg5_pos)

aic9.1 <- AIC(fit9.1) ### WINNER
aic9.2 <- AIC(fit9.2)
aic9.3 <- AIC(fit9.3) 

min9 <- min(c(aic9.1,aic9.2,aic9.3))
(min9 == c(aic9.1,aic9.2,aic9.3))

fit10 <- glm(cphl~depth+irrad443+cdom+temp+psal+irrad555+irrad490+dox1, start = rep(1,9), family = Gamma(link = "inverse"), data = model_data_avg5_pos)
fit10.1 <- glm(cphl~depth+irrad443+cdom+temp+psal+irrad490+irrad555+dox1+vbsc, start = rep(1,10), family = Gamma(link = "inverse"), data = model_data_avg5_pos)
fit10.2 <- glm(cphl~depth+irrad443+cdom+temp+psal+irrad490+irrad555+dox1+irrad670, start = rep(1,10), family = Gamma(link = "inverse"), data = model_data_avg5_pos)

aic10.1 <- AIC(fit10.1) ### WINNER
aic10.2 <- AIC(fit10.2)

min10 <- min(c(aic10.1,aic10.2))
(min10 == c(aic10.1,aic10.2))

fit11 <- glm(cphl~depth+irrad443+cdom+temp+psal+irrad490+irrad555+dox1+vbsc, start = rep(1,10), family = Gamma(link = "inverse"), data = model_data_avg5_pos)
fit11.1 <- glm(cphl~depth+irrad443+cdom+temp+psal+irrad490+irrad555+dox1+vbsc+irrad670, start = rep(1,11), family = Gamma(link = "inverse"), data = model_data_avg5_pos)

min11 <- AIC(fit11.1)

ultimate_min_AIC <- min(c(min3,min4,min5,min6,min7,min8,min9,min10,min11))
(ultimate_min_AIC == c(min3,min4,min5,min6,min7,min8,min9,min10,min11))


####
# Final Models that we compare for AIC and RMSE
####
gamma_invsmall <- fit11.1
gamma_log <- glm(cphl~., family = Gamma(link = "log"), data = model_data_avg5_pos)
gamma_invfull <- glm(cphl~., family=Gamma(link = "inverse"),
                     data = model_data_avg5_pos, start = rep(1,12))
gamma_invno670 <- glm(cphl~.-irrad670, family=Gamma(link = "inverse"),
                      data = model_data_avg5_pos, start = rep(1,11))

final_model <- gamma_invfull


#################################################
###Load Test Sets
test1 <- read.csv("Data/test_autumn.csv")
test2 <- read.csv("Data/test_spring.csv")
test3 <- read.csv("Data/test_summer.csv")
test4 <- read.csv("Data/test_winter.csv")

### RMSE on Gamma Inverse (no cndc) - gamma_invno760
a1 <- test1$cphl
a1.pred <- predict(gamma_invno670, newdata = test1, type = "response")
sum(is.na(a1.pred))
rmse1_no670 <- sqrt(mean((a1-a1.pred)^2))

a2 <- test2$cphl
a2.pred <- predict(gamma_invno670, newdata = test2, type = "response")
sum(is.na(a2.pred))
rmse2_no670 <- sqrt(mean((a2-a2.pred)^2))

a3 <- test3$cphl
a3.pred <- predict(gamma_invno670, newdata = test3, type = "response")
sum(is.na(a3.pred))
rmse3_no670 <- sqrt(mean((a3-a3.pred)^2))

a4 <- test4$cphl
a4.pred <- predict(gamma_invno670, newdata = test4, type = "response")
sum(is.na(a4.pred))
rmse4_no670 <- sqrt(mean((a4-a4.pred)^2))
gamma_invno670.rmse <- c(rmse1_no670,rmse2_no670,rmse3_no670,rmse4_no670)

### RMSE on Gamma Inverse (no cndc) - gamma_invsmall
y1 <- test1$cphl
y1.pred <- predict(gamma_invsmall, newdata = test1, type = "response")
sum(is.na(y1.pred))
rmse1 <- sqrt(mean((y1-y1.pred)^2))

y2 <- test2$cphl
y2.pred <- predict(gamma_invsmall, newdata = test2, type = "response")
sum(is.na(y2.pred))
rmse2 <- sqrt(mean((y2-y2.pred)^2))

y3 <- test3$cphl
y3.pred <- predict(gamma_invsmall, newdata = test3, type = "response")
sum(is.na(y3.pred))
rmse3 <- sqrt(mean((y3-y3.pred)^2))

y4 <- test4$cphl
y4.pred <- predict(gamma_invsmall, newdata = test4, type = "response")
sum(is.na(y4.pred))
rmse4 <- sqrt(mean((y4-y4.pred)^2))
gamma_invsmall.rmse <- c(rmse1,rmse2,rmse3,rmse4)

### RMSE on Gamma Inverse - gamma_invfull
z1 <- test1$cphl
z1.pred <- predict(gamma_invfull, newdata = test1, type = "response")
sum(is.na(z1.pred))
rmse1_invF <- sqrt(mean((z1-z1.pred)^2))

z2 <- test2$cphl
z2.pred <- predict(gamma_invfull, newdata = test2, type = "response")
sum(is.na(z2.pred))
rmse2_invF <- sqrt(mean((z2-z2.pred)^2))

z3 <- test3$cphl
z3.pred <- predict(gamma_invfull, newdata = test3, type = "response")
sum(is.na(z3.pred))
rmse3_invF <- sqrt(mean((z3-z3.pred)^2))

z4 <- test4$cphl
z4.pred <- predict(gamma_invfull, newdata = test4, type = "response")
sum(is.na(z4.pred))
rmse4_invF <- sqrt(mean((z4-z4.pred)^2))
gamma_invfull.rmse <- c(rmse1_invF,rmse2_invF,rmse3_invF,rmse4_invF)

### RMSE on Gamma Log - gamma_log
x1 <- test1$cphl
x1.pred <- predict(gamma_log, newdata = test1, type = "response")
sum(is.na(x1.pred))
rmse1_log <- sqrt(mean((x1-x1.pred)^2))

x2 <- test2$cphl
x2.pred <- predict(gamma_log, newdata = test2, type = "response")
sum(is.na(x2.pred))
rmse2_log <- sqrt(mean((x2-x2.pred)^2))

x3 <- test3$cphl
x3.pred <- predict(gamma_log, newdata = test3, type = "response")
sum(is.na(x3.pred))
rmse3_log <- sqrt(mean((x3-x3.pred)^2))

x4 <- test4$cphl
x4.pred <- predict(gamma_log, newdata = test4, type = "response")
sum(is.na(x4.pred))
rmse4_log <- sqrt(mean((x4-x4.pred)^2))
gamma_log.rmse <- c(rmse1_log,rmse2_log,rmse3_log,rmse4_log)

# Compare all
gamma_invfull.rmse ###### Winner
gamma_invno670.rmse
gamma_invsmall.rmse
gamma_log.rmse

# Thus gamma_invfull is the best glm model and used as the final