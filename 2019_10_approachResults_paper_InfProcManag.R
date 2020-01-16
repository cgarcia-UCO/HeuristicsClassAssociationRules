#############AUTHORSHIP INFORMATION########
# author: Jose Antonio Delgado Osuna and Carlos Garcia-Martinez
# Affiliation: University Hospital Reina Sofia, University of Córdoba
# email: jose.delgado.osuna@gmail.com, cgarcia@uco.es
# init date: Jul 2019

##############OBJECTIVE####################
# This file provides the source of the experiments carried out for the proposal published in
# José A. Delgado, Carlos García-Martínez, José Gómez Barbadillo, Sebastian Ventura.
# Heuristics for interesting class association rule mining a colorectal cancer database
# Information Processing & Management, 2020


##############LICENSE######################
# MIT License
#
# Copyright (c) 2019 cgarcia-UCO
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

## Required packages
# Packages to read arff files
if (!requireNamespace("foreign", quietly = TRUE))
  install.packages("foreign")
if (!requireNamespace("plyr", quietly = TRUE))
  install.packages("plyr")
if (!requireNamespace("stringr", quietly = TRUE))
  install.packages("stringr")

# Others
library(foreign)
library(plyr) # for rbind.fill
library(stringr) # for str_replace_all

debugSource("src/accValuesRules.R")
debugSource("src/rulesExtOperations.R")
debugSource("src/utils.R")



maxLengthRules <- 5
resultsDir <- "results"

files <- c(
  'Datasets/Ntto_qt_paliativo_no/CA_COLORRECTAL_COMPLICACIONES.arff',
  'Datasets/Ntto_qt_paliativo_no/CCR_SEGUIMIENTO_RECIDIVA.arff'
)
j <- files[1]


for (j in files) {
  filename <- j
  
  #Read the dataset
  suppressWarnings(dataset <- read.arff(filename))

  #Remove uninteresting features from the dataset before extracting the rules  
  uninterestingFeatures <- getUninterestingFeatures(filename)
  uninterestingFeatures2 <-
    getSecondUninterestingFeatures(filename)
  dataset[uninterestingFeatures] <- NULL
  dataset[uninterestingFeatures2] <- NULL
  
  #Get the name and value of the goal feature
  features <- colnames(dataset)
  numFeatures <- length(features)
  goalFeature <- features[numFeatures]
  goalValueIndex <- which.min(table(dataset[goalFeature]))
  goalValue <- names(goalValueIndex)
  
  #This software/proposal works just on datasets with two goal values
  if (length(unique(dataset[goalFeature])[[1]]) == 2) {
    print(filename)
    print(summary(dataset[goalFeature]))
    time1 <- Sys.time()
    
    #Aplication of the proposal to obtain the rules
    rules <-
      getRules(
        dataset,
        maxLength = maxLengthRules,
        justLeaves = TRUE,
        simplifyRules = 1,
        minConfidence = 0.7,
        recursion = maxLengthRules,
        maxTimeRecursion = -1
      )
    time2 <- Sys.time()
    usedTime <- difftime(time2, time1, units = c("secs"))
    cat(paste("Required time: ", as.integer(usedTime), "s", sep = ""))
    
    
    #Slightly format the output into another dataframe
    outputRules <- rules[FALSE, ]
    
    if (nrow(rules) > 0) {
      for (j in 1:nrow(rules)) {
        rule <- rules[j, ]
        outputRules <-
          addRuleToSet(
            rule,
            rule$iteration,
            outputRules,
            simplification = rule$simplification,
            oneCatAltIndex = rule$oneCatAltIndex,
            dataset = dataset,
            goalFeature = goalFeature,
            goalValue = goalValue
          )
      }
    }
    
    
    #Generate columns with the name of the features in the conditions
    outputRules["features"] <- NULL
    for (i in 1:nrow(outputRules)) {
      outputRules[i, "features"] <-
        paste(lapply(rules[i, "rule"][[1]], getFeatureFromCondition), collapse = " & ")
    }
    
    #Write the rules in a csv
    outputRules[[1]] <- apply(outputRules, 1, FUN = branchToString, dataset)
    outputRules <- outputRules[order(outputRules$globalQuality, decreasing = TRUE), ]
    destFile <-
      paste(resultsDir,
            "/",
            str_replace_all(filename, "/", "_"),
            ".csv",
            sep = "")
    writeCSV(outputRules, destFile)
    cat("CSV written\n")
    
    
    #Output som statistic information of the results
    outputInfo(rules, dataset, goalFeature, goalValue)
  }
}
