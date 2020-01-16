#############AUTHORSHIP INFORMATION########
# author: Jose Antonio Delgado Osuna and Carlos Garcia-Martinez
# Affiliation: University Hospital Reina Sofia, University of Córdoba
# email: jose.delgado.osuna@gmail.com, cgarcia@uco.es
# init date: Jul 2019

##############OBJECTIVE####################
# This file provides some functions to measure the quality of descripting rules according
# to the information of a dataset, such as, confidence, lift...
# The data structure of rules is that in file src/rulesFunctions.R


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


debugSource("src/rulesFunctions.R")

#This function returns the number of positive patterns covered by the branch
computeNumCovPostCases <- function(branch,dataset, goalFeature, goalValue){
  string <- addDatasetVar(branchToString(branch, dataset), "dataset")
  pIntersect <- eval(parse(text=paste('sum(dataset$',goalFeature,' == \"',goalValue,'\" & (', string,')',', na.rm = TRUE)', sep = "")))
  return(c(covPostCases = as.numeric(pIntersect)))
}

#This function computes the extended complexity of a rule as the product of the alternatives in its conditions
getExtendedComplexity <- function(branch, dataset){
  length <- 1
  
  for (i in branch){
    
    featureName <- getFeatureFromCondition(i)
    
    if (isCategorical(featureName, dataset)){
      alternatives <- getAlternativesCatCondition(i)
      numAlternatives <- length(alternatives[[1]])
      length <- length * numAlternatives # I am implementing the product of combinations
    }
  }
  
  return(length)
}

#This function computes the percentage of positive patterns covered and the confidence of a rule
computeCovPostCases_confidence<-function(branch,dataset,originalDataset, goalFeature, goalValue){
  string <- addDatasetVar(branchToString(branch, dataset), "dataset")
  pGoalValue <- sum(originalDataset[goalFeature] == goalValue, na.rm = TRUE)
  pString <- eval(parse(text=paste('sum(',string,', na.rm = TRUE)', sep = "")))
  pIntersect <- eval(parse(text=paste('sum(dataset$',goalFeature,' == \"',goalValue,'\" & (', string,')',', na.rm = TRUE)', sep = "")))
  covPosCases <- (pIntersect/pGoalValue)
  confidence <- (pIntersect/pString)
  return(c(covPostCases = as.numeric(covPosCases), covPostCases2 = pIntersect, confidence = as.numeric(confidence)))
}

#This function computes several accuracy metrics of a rule, say: percentage of positive patterns covered, absolute value of the that value,
#percentage of patterns that satisfy the conditions of the rule, confidence of the rule, lift, the division of the product of covPostCases and confidence by
#the extended complexity of the rule, and the extended complexity of the rule
computeAccuracyValues<-function(branch, dataset, originalDataset, goalFeature, goalValue){
  
  string <- addDatasetVar(branchToString(branch, dataset), "dataset")
  numPatterns <- nrow(originalDataset)#eval(parse(text=paste('nrow(', nombredataset,')')))
  pGoalValue <- sum(originalDataset[goalFeature] == goalValue, na.rm = TRUE)
  pString <- eval(parse(text=paste('sum(',string,', na.rm = TRUE)', sep = "")))
  pIntersect <- eval(parse(text=paste('sum(dataset$',goalFeature,' == \"',goalValue,'\" & (', string,')',', na.rm = TRUE)', sep = "")))
  covPosCases <- (pIntersect/pGoalValue)
  confidence <- (pIntersect/pString)
  lift<-round((pIntersect/numPatterns) /
                (pGoalValue/numPatterns * pString/numPatterns),
              digits=2)

  extLength <- getExtendedComplexity(branch, dataset)
  globalQuality<-confidence * covPosCases / extLength

    return(c(covPostCases = as.numeric(covPosCases), covPostCases2 = pIntersect, support = pString,
           confidence = as.numeric(confidence), lift = as.numeric(lift),
           globalQuality = as.numeric(globalQuality),
           extLength = as.numeric(extLength)))
}

#This function computes several accuracy metrics of a rule provided in a string format. The metrics are: percentage of positive patterns covered, absolute value of the that value,
#percentage of patterns that satisfy the conditions of the rule, confidence of the rule, and lift
computeAccuracyValuesFromString<-function(string, dataset, originalDataset, goalFeature, goalValue){
  
  string <- addDatasetVar(string, "dataset")
  numPatterns <- nrow(originalDataset)#eval(parse(text=paste('nrow(', nombredataset,')')))
  pGoalValue <- sum(originalDataset[goalFeature] == goalValue, na.rm = TRUE)
  pString <- eval(parse(text=paste('sum(',string,', na.rm = TRUE)', sep = "")))
  pIntersect <- eval(parse(text=paste('sum(dataset$',goalFeature,' == \"',goalValue,'\" & (', string,')',', na.rm = TRUE)', sep = "")))
  covPosCases <- (pIntersect/pGoalValue)
  confidence <- (pIntersect/pString)
  lift<-round((pIntersect/numPatterns) /
                (pGoalValue/numPatterns * pString/numPatterns),
              digits=2)
  return(c(covPostCases = as.numeric(covPosCases), covPostCases2 = pIntersect, support = pString,
           confidence = as.numeric(confidence), lift = as.numeric(lift)))
}


#This function returns the percentage of patterns that satisfy each condition in the given rule
getIndvProbs <- function(branch, dataset, originalDataset){
  probsIndv <- c()
  numPatterns <- nrow(originalDataset)
  
  for (i in branch){
    auxString <- addDatasetVar(branchToString(i, dataset), "dataset")
    aProb <- eval(parse(text=paste('sum(', auxString, ', na.rm = TRUE)', sep = "")))
    probsIndv <- c(probsIndv, aProb / numPatterns)
  }
  
  return(probsIndv)
}


#This function returns, per each condition in the rule, the confidence of the rule without that condition (in negative)
getConfReductionConditions<-function(branch, dataset, originalDataset, goalFeature, goalValue){
  confidenceBranch <- computeCovPostCases_confidence(branch, dataset, dataset, goalFeature, goalValue)
  confidences <- c()
  numConditions <- length(branch)
  
  if (numConditions > 1){
    for (it_Conditions in 1:numConditions){
      confidenceAux <- computeAccuracyValues(branch[-it_Conditions], dataset, dataset, goalFeature, goalValue)
      confidences <- c(confidences, -1 * confidenceAux[["confidence"]])
    }
  } else {
    confidences <- c(-1 * sum(dataset[goalFeature] == goalValue) / nrow(dataset))
  }
  
  return(confidences)
}

#This function returns, per each condition in the rule, the confidence of the rule with just that condition
getConfidencesEachCondition<-function(branch, dataset, originalDataset, goalFeature, goalValue){
  confidenceBranch <- computeCovPostCases_confidence(branch, dataset, dataset, goalFeature, goalValue)
  confidences <- c()
  numConditions <- length(branch)
  
  if (numConditions > 1){
    for (it_Conditions in 1:numConditions){
      confidenceAux <- computeAccuracyValues(branch[it_Conditions], dataset, dataset, goalFeature, goalValue)
      confidences <- c(confidences, confidenceAux[["confidence"]])
    }
  } else {
    confidences <- c(confidenceBranch["confidence"])
  }
  
  return(confidences)
}

#This function indicates if a ruleset covers or does not cover a pattern
ruleSetCovers <- function(ruleset, case, dataset){
  
  for (aRule in ruleset){
    
    if (ruleCovers(aRule, case, dataset))
      return(TRUE)
  }
  
  return(FALSE)
}

#This function indicates if a rule covers or does not cover a pattern
ruleCovers <- function(rule, case, dataset){
  
  for (condition in rule){
    
    if (!conditionCovers(condition, case, dataset)){
      return(FALSE)
    }
  }
  
  return(TRUE)
}

#This function indicates if a condition is satisfied by a pattern
conditionCovers <- function(condition, case, dataset){
  
  feature <- getFeatureFromCondition(condition)
  levelsF <- levels(dataset[,feature])
  
  if (is.na(case[feature])){
    return(FALSE)
  } else if (isCategorical(feature, dataset)){
    alternatives <- getAlternativesCatCondition(condition)
    caseValue <- paste("\"", levelsF[case[,feature]], "\"", sep="")
    result <- isInList(case[feature], as.list(alternatives[[1]]))
    
    if (!result){
      return(FALSE)
    }
  } else {
    result <- eval(parse(text=paste('case$',condition, sep = "")))
    
    if (!result){
      return(FALSE)
    }
  }
  return(TRUE)
}

#This function returns the percentage of positive patterns that are not covered by any rule in rules
getUncoveredPositives <- function(rules, dataset, goalFeature, goalValue){
  
  posCases <- dataset[dataset[goalFeature] == goalValue,]
  numPosCases <- nrow(posCases)
  numCovered <- 0

  for (iCase in 1:numPosCases){
    case <- posCases[iCase,]
    isCovered <- FALSE
    
    if(ruleSetCovers(rules, case, dataset)){
      numCovered <- numCovered + 1
    }
  }
  
  return (1. - (numCovered / numPosCases))
}
