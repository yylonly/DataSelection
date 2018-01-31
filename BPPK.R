library(prodlim)
library(dplyr)
library(readr)

rm(list = ls())

source("C:/Users/yylonly/Desktop/PharmaceuticsData/trunk/BPPK/MDFIS.R")

isAPIgroupValidationSet = FALSE
isAPIgroupTestingSet = FALSE
validationSetRate = 0.2
testingSetRate = 0.1
xIndex = 2:10
yIndex = 11

# fileName = "WaterSolubility.csv"

# fileName = "logP.csv"

# fileName = "Caco2Permeability.csv"

# fileName = "Absorption.csv"

# fileName = "Bioavailability.csv"

# fileName = "HalfLife.csv"

# fileName = "PlasmaProteinBindingRate.csv"

# fileName = "VolumeofDistribution.csv"

fileName = "ADME.csv"

# create subfolder
fileNamewithOut = tools::file_path_sans_ext(fileName)
dir.create(file.path(fileNamewithOut))



#load data
alldata <- read_csv(fileName)
y <- data.matrix(alldata)[, yIndex]
X <- data.matrix(alldata)[, xIndex]
alldata <- cbind(X, y)

maxs <- apply(alldata, 2, max)
mins <- apply(alldata, 2, min)
ranges <- maxs - mins
means <- apply(alldata, 2, mean)
scaledall <- scale(alldata, center = mins, scale = ranges)


## Get best inital dataset
numbers = dim(scaledall)[1];
testNum = as.integer(numbers*testingSetRate)
validNum = as.integer(numbers*validationSetRate)
trainNum = numbers - testNum - validNum

allIndexes <- NULL
allsumdiss <- NULL

################################### Validation Set #######################################################

## Generate 10000 intial data set and get best one
for (i in 1:10000) {
  ## A random sample of 5 data points
  set.seed(i)
  initalIndexes <- sample(numbers, 5)
  
  TrainningSet <- scaledall[-initalIndexes, ]
  initalTestSet <- scaledall[initalIndexes, ]
  
  allIndexes <- rbind(allIndexes, initalIndexes)
  
  diss <- proxy::dist(initalTestSet, TrainningSet)
  sumdiss <- sum(diss)
  allsumdiss <- c(allsumdiss, sumdiss)
  
}

bestInitalIndex <- allIndexes[which.min(allsumdiss), ]
bestDistance <- min(allsumdiss)

#Begin compute remaining testset
RemainingSetV <- scaledall[-bestInitalIndex, ]
initalSetV <- scaledall[bestInitalIndex, ]

## MD-FIS
validaIndex <- MDFIS(initalSetV, RemainingSetV, n = (validNum - 5), obj = minDiss, alpha = 3, isAPIgroup = isAPIgroupValidationSet)
SelectedSetV <- RemainingSetV[validaIndex, ]
validationSet <- rbind(initalSetV, SelectedSetV)
RemainSetForTestingSelection <- RemainingSetV[-validaIndex, ]

numbers <- dim(RemainSetForTestingSelection)[1];

################################### Testing Set ##############################################
## Generate 10000 intial data set and get best one
for (i in 1:10000) {
  ## A random sample of 5 data points
  set.seed(i)
  initalIndexes <- sample(numbers, 5)
  
  TrainningSet <- RemainSetForTestingSelection[-initalIndexes, ]
  initalTestSet <- RemainSetForTestingSelection[initalIndexes, ]
  
  allIndexes <- rbind(allIndexes, initalIndexes)
  
  diss <- proxy::dist(initalTestSet, TrainningSet)
  sumdiss <- sum(diss)
  allsumdiss <- c(allsumdiss, sumdiss)
  #initalIndexes <- c(14, 30, 46, 54, 91)
  #initalIndexes <- c(5,50,78,99,117)
  #initalIndexes <- c(18,64,65,66,84)
  #initalIndexes <- c(18,64,65,66,74,83,84)
  
}

bestInitalIndex <- allIndexes[which.min(allsumdiss), ]
bestDistance <- min(allsumdiss)

#Begin compute remaining testset
RemainingSetV <- RemainSetForTestingSelection[-bestInitalIndex, ]
initalSetV <- RemainSetForTestingSelection[bestInitalIndex, ]


## MD-FIS
validaIndex <- MDFIS(initalSetV, RemainingSetV, n = (testNum - 5), obj = minDiss, alpha = 3, isAPIgroup = isAPIgroupTestingSet)
SelectedSetV <- RemainingSetV[validaIndex, ]

TestingSet <- rbind(initalSetV, SelectedSetV)
TrainingSet <- RemainingSetV[-validaIndex, ]

## write to file
write.csv(TrainingSet, paste(fileNamewithOut, "trainingset.csv", sep = "/"), row.names = FALSE)
write.csv(TestingSet, paste(fileNamewithOut, "testingset.csv", sep = "/"), row.names = FALSE)
write.csv(validationSet, paste(fileNamewithOut, "validationset.csv", sep = "/"), row.names = FALSE)



