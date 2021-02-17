###RANDOM FOREST MODEL--------------------------------------------------------------------------------------
install.packages('rsample')
install.packages('randomForest')
install.packages('ranger')
install.packages('caret')
install.packages('h2o')
library(rsample)      # data splitting 
library(randomForest) # basic implementation
library(ranger)       # a faster implementation of randomForest
library(caret)        # an aggregator package for performing many machine learning models
library(h2o) 
library(dplyr)

# removing rows with missing data, dropping variables that don't contribute
model_data <- na.omit(model_data)
model_data1 <- model_data[-c(2,3,4,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34)]

# set seed for reproducibility, splitting into training and testing sets
ss <- model_data1[sample(nrow(model_data1), 8523), ]
sample_size = floor(0.75*nrow(ss))
sample_size
set.seed(123)
train_ind = sample(seq_len(nrow(model_data1)),size = sample_size)
train = model_data1[train_ind,]
test = model_data1[-train_ind,]

m1 <- randomForest(formula = cphl~., data = train)
predict(m1,data=test)
plot(m1) # can see error rate stabilises around 100 trees or so
# to find mine error (number of trees)
which.min(m1$mse) #yields 477 trees (too many)

# looking at variable importance - INCNODEPURITY (This is a measure of variable importance based on the Gini impurity index used for the calculating the splits in trees. The higher the value of mean decrease accuracy or mean decrease gini score, the higher the importance of the variable to our model)

varImpPlot(m1) 

# can further tune using ntree (number of trees), mtry(the number of variables to randomly sample as candidates at each split),sampsize (samplesize), nodesize,maxnodes


### RANGER
model_data1 <- filter(model_data,cphl_flag == 1 & pres_flag == 1 & depth_flag == 1 & temp_flag == 1 & cndc_flag == 1 & psal_flag == 1 & dox1_flag == 1 & cdom_flag == 1 & vbsc_flag == 1 & bbp_flag == 1 &  irrad443_flag == 1 & irrad490_flag == 1 & irrad555_flag == 1 & irrad670_flag == 1)
# removing rows with missing data, dropping variables that don't contribute
model_data <- na.omit(model_data1)
model_data <- model_data1[-c(2,3,4,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34)]

m1 <- model_data[1:200000,]
m2 <- model_data[200001:389675,]

m1_fit <- ranger(cphl ~ ., data = m1)
# parameters (num.trees = 500,max.depth = 8,probability = TRUE)
m1_fit$num.trees
m1_fit$mtry
#O verall out of bag prediction error (prediciton.error), for regression uses MSE
m1_fit$prediction.error
m1_fit$r.squared

