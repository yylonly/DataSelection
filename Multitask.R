library(prodlim)
library(dplyr)
library(readr)

rm(list = ls())

source("MDFIS.R")
source("OverrideFunc.R")

isAPIgroupValidationSet = FALSE
isAPIgroupTestingSet = FALSE
validationSetRate = 0.2
testingSetRate = 0.2
pIndex = 2:9
xIndex = 11:1034
yIndex = 1035:1038

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

dataM <- data.matrix(alldata)

P <- dataM[, pIndex]
y <- dataM[, yIndex]
X <- dataM[, xIndex]
#alldata <- cbind(X, y)
#alldata <- data.matrix(alldata)

# normalization P
maxsP <- apply(P, 2, maxN)
minsP <- apply(P, 2, minN)
rangesP <- maxsP - minsP
meansP <- apply(P, 2, meanN)
scaledP <- scale(P, center = minsP, scale = rangesP)

# normalization y

maxsy <- apply(y, 2, maxN)
minsy <- apply(y, 2, minN)
rangesy <- maxsy - minsy
meansy <- apply(y, 2, meanN)
scaledy <- scale(y, center = minsy, scale = rangesy)

scalePFy <- cbind(scaledP, X, scaledy)
scaledall <- scaledP

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

###for PFy
RemainingPFy <- scalePFy[-bestInitalIndex, ]
initalPFy <- scalePFy[bestInitalIndex, ]

## MD-FIS
validaIndex <- MDFIS(initalSetV, RemainingSetV, n = (validNum - 5), obj = minDiss, alpha = 3, isAPIgroup = isAPIgroupValidationSet)
SelectedSetV <- RemainingSetV[validaIndex, ]
validationSet <- rbind(initalSetV, SelectedSetV)
RemainSetForTestingSelection <- RemainingSetV[-validaIndex, ]

###for PFy
SelectedPFy <- RemainingPFy[validaIndex, ]
validationPFy <- rbind(initalPFy, SelectedPFy)
RemainSetForTestingSelectionPFy <- RemainingPFy[-validaIndex, ]



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

###for PFy
RemainingPFy <- RemainSetForTestingSelectionPFy[-bestInitalIndex, ]
initalPFy <- RemainSetForTestingSelectionPFy[bestInitalIndex, ]

## MD-FIS
validaIndex <- MDFIS(initalSetV, RemainingSetV, n = (testNum - 5), obj = minDiss, alpha = 3, isAPIgroup = isAPIgroupTestingSet)
SelectedSetV <- RemainingSetV[validaIndex, ]
TestingSet <- rbind(initalSetV, SelectedSetV)
TrainingSet <- RemainingSetV[-validaIndex, ]

###for PFy
SelectedPFy <- RemainingPFy[validaIndex, ]
TestingSetPFy <- rbind(initalPFy, SelectedPFy)
TrainingSetPFy <- RemainingPFy[-validaIndex, ]

## write to file
write.csv(TrainingSetPFy, na = "NaN", paste(fileNamewithOut, "trainingset.csv", sep = "/"), row.names = FALSE)
write.csv(validationPFy, na = "NaN", paste(fileNamewithOut, "testingset.csv", sep = "/"), row.names = FALSE)
write.csv(TestingSetPFy, na = "NaN", paste(fileNamewithOut, "validationset.csv", sep = "/"), row.names = FALSE)

