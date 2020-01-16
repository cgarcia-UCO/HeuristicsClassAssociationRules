#############AUTHORSHIP INFORMATION########
# author: Jose Antonio Delgado Osuna and Carlos Garcia-Martinez
# Affiliation: University Hospital Reina Sofia, University of Córdoba
# email: jose.delgado.osuna@gmail.com, cgarcia@uco.es
# init date: Jul 2019

##############OBJECTIVE####################
#This file provides some extended functions to work with the conditions provided by a CART tree.
#This functions try to optimise the rule somehow, either removing conditions, alterantives, or looking for one categorical alterantive rules

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

library(stringr) #For str_pad


source("src/rulesFunctions.R")
source("src/accValuesRules.R")
source("src/utils.R")

#This function removes posible duplicate conditions from a rule. In particular it removes less restrictive conditions than
#others in the rule. In case of categorical conditions, the alternatives in the most restrictive condition are in those of
#the other, which contains additional alternatives. It assumes that the branch from the CART tree have more restrictive
#conditions near the leaves (as is expected and was observed). In case of numerical conditions, it removes the conditions
#for which there is another nearer the leaves with the same name of the feature and same numerical operator.
removeDuplicateConditions <- function(branch, dataset){
  library(gtools)
  numConditions <- length(branch)
  i <- 1
  
  while (i < numConditions){
    iCondition <- getFeatureFromCondition(branch[i])
    j<-i+1
    
    for (j in (i+1):numConditions){
      jCondition <- getFeatureFromCondition(branch[j])
      sameFeature <- (iCondition == jCondition)

      if (sameFeature && isCategorical(iCondition, dataset)){
        iAlternatives <- getAlternativesCatCondition(branch[i])
        jAlternatives <- getAlternativesCatCondition(branch[j])
        
        if (all (jAlternatives[[1]] %in% iAlternatives[[1]])){
          branch <- branch[-i]
          i <- max(0, i - 1)
          numConditions <- length(branch)
          break()
        }
      } else if (sameFeature){ # The else is for numerical conditions
        if (getOperatorNumCondition(branch[i]) == getOperatorNumCondition(branch[j])){
          branch <- branch[-i]
          i <- max(0, i - 1)
          numConditions <- length(branch)
          break()
        }
      }
    }
    
    i <- i + 1
  }
  
  return(branch)
}

#Next function just wraps common instructions to insert a new state in the search stack of the simplifyBranch function.
#This alloows simplifyBranch function look simpler
simplifyBranch_addState <- function(queueState, queueAcValues, queueCoding, queueDepth, newBranch, newValues, newCoding, currentDepth){
  queueState_x <- queueState
  queueAcValues_x <- queueAcValues
  queueCoding_x <- queueCoding
  queueDepth_x <- queueDepth
  
  queueState_x <- c(list(newBranch),queueState_x)
  queueAcValues_x <- c(list(newValues),queueAcValues_x)
  queueCoding_x <- c(list(newCoding),queueCoding_x)
  queueDepth_x <- c(list(currentDepth+1), queueDepth_x)
  
  eval.parent(substitute(queueState<-queueState_x))
  eval.parent(substitute(queueAcValues<-queueAcValues_x))
  eval.parent(substitute(queueCoding<-queueCoding_x))
  eval.parent(substitute(queueDepth<-queueDepth_x))
}


#This function return rules that are constructed from branch by removing each of its conditions, in case these new rules
#satisfy the minConfidence threshold
getReducedRules <- function(branch, dataset, originalDataset, goalFeature, goalValue,
                            minConfidence = 0.7){
  
  outputRules <- list()
  frontier <- list()
  frontier <- c(frontier, list(branch))
  
  while (length(frontier) >= 1){
    
    currentBranch <- frontier[[1]]
    frontier[[1]] <- NULL
    numConditions <- length(currentBranch)
    
    if (numConditions > 1){
      for (it_Conditions in 1:numConditions){
        auxBranch <- currentBranch[-it_Conditions]
        confidenceAux <- computeCovPostCases_confidence(auxBranch, dataset, dataset, goalFeature, goalValue)
        
        if (confidenceAux[["confidence"]] >= minConfidence){
          opRules <- simplifyBranch(auxBranch, dataset, originalDataset, #In some cases, removing one condition produces a rule that might by optimised with simplifyBranch (i.e., other conditions may be removed too, and the rule only improves)
                                    goalFeature, goalValue, maxSubRulePerRule = 1)
          
          if (length(opRules) > 0){
            for (iOpRule in opRules){
              outputRules <- c(outputRules, list(iOpRule))
              frontier <- c(list(iOpRule), frontier)
            }
          } else {
            outputRules <- c(outputRules, list(auxBranch))
            frontier <- c(list(auxBranch), frontier)
          }
        }
      }
    }
  }
  
  return(outputRules)
}

#This function return rules that are constructed from branch by removing all the alternatives but one from each of
#its categorical conditions, in case these new rules satisfy the minConfidence threshold. The function continues
#constructing new rules from those generated, performing a complete search process
getOneCatValueRules <- function(branch, dataset, originalDataset, goalFeature, goalValue,
                                minConfidence = 0.8, minPositiveCases = 7){
  
  outputRules <- list()
  
  frontier <- list()
  frontier <- c(frontier, list(branch))
  
  while (length(frontier) >= 1){
    
    currentBranch <- frontier[[1]]
    frontier[[1]] <- NULL
    numConditions <- length(currentBranch)
    
    for (it_Conditions in 1:numConditions){
      condition <- currentBranch[it_Conditions]
      featureName <- getFeatureFromCondition(condition)
      
      if (isCategorical(featureName, dataset)){
        alternatives <- getAlternativesCatCondition(condition)
        numAlternatives <- length(alternatives[[1]])
        
        if (numAlternatives > 1){
          i <- alternatives[1]
          
          for (i in 1:numAlternatives){
            auxCondition<-paste(featureName, "=", alternatives[[1]][i], sep = "")
            auxBranch <- currentBranch
            auxBranch[it_Conditions] <- auxCondition
            newValues<-computeCovPostCases_confidence(auxBranch, dataset,
                                                      originalDataset, goalFeature, goalValue)
            
            if (newValues["confidence"] >= minConfidence &&
                newValues["covPostCases2"] >= minPositiveCases){
              opRules <- simplifyBranch(auxBranch, dataset, originalDataset, #In some cases, modifying one condition produces a rule that might by optimised with simplifyBranch (i.e., other conditions may be removed too, and the rule only improves)
                                        goalFeature, goalValue, maxSubRulePerRule = 1)
              
              if (length(opRules) > 0){
                for (iOpRule in opRules){
                  outputRules <- c(outputRules, list(iOpRule))
                  frontier <- c(list(iOpRule), frontier)
                }
              } else {
                outputRules <- c(outputRules, list(auxBranch))
                frontier <- c(list(auxBranch), frontier)
              }
            }
          }
        }
      }
    }
  }
  
  return(outputRules)
}

#This function implements a depth-first search trying to simplify the rule/branch provided
#This depth-first search looks for improvements in covPosCases and confidence when removing conditions and/or alternatives
#At most, it evaluates maxSubRulePerRule - currentDepth alternatives per each state. Evaluating all the posibilities is too expensive,
#and evaluating just the first one improving (maxSubRulePerRule = 1) the current state may probably avoid interesting results.
simplifyBranch <- function(branch1, dataset, originalDataset, goalFeature, goalValue, maxSubRulePerRule = 1){
  library(gtools)
  explored <- list()
  exploredAndBetter <- list()
  queueState <- list()
  queueAcValues <- list()
  queueCoding <- list()
  queueDepth <- list()
  outputRules <- list()
  hits <- 0
  allHits <- 0
  countStates <- 0
  
  branch <- removeDuplicateConditions(branch1, dataset)
  
  guidingMetric <- "covPostCases" # This search process must maintain or improve this metric

  improvedAtLeastOnce <- FALSE
  anyImprovement <- FALSE #This controls that the original rule was, at least once, improved. This variable is different from the previous one becuase the previous one is used for the current optimization branch, instead of the original rule
  
  continueTrying <- TRUE
  currentBranch <- branch
  currentValues <- computeCovPostCases_confidence(currentBranch, dataset, originalDataset, goalFeature, goalValue)

  queueState <- c(queueState, list(currentBranch))
  queueAcValues <- c(queueAcValues, list(currentValues))
  queueCoding <- c(queueCoding, list("-"))
  queueDepth <- c(queueDepth, list(0))
  
  while (length(queueState) >= 1){
    currentBranch <- queueState[[1]]
    currentValues <- queueAcValues[[1]]
    currentCoding <- queueCoding[[1]]
    currentDepth <- queueDepth[[1]]
    queueState[[1]] <- NULL
    queueAcValues[[1]] <- NULL
    queueCoding[[1]] <- NULL
    queueDepth[[1]] <- NULL
    improvedAtLeastOnce <- FALSE #This controls that we are generating further children states. Leaves have this equal to FALSE
    numSubRules <- 0
    
    numConditions <- length(currentBranch)
    conditionsIndexes <- permute(1:numConditions)
    it_Conditions <- conditionsIndexes[1]
    
    for (it_Conditions in conditionsIndexes){
      
      condition <- currentBranch[it_Conditions]
      featureName <- getFeatureFromCondition(condition)
      
      if (isCategorical(featureName, dataset)){
        alternatives <- getAlternativesCatCondition(condition)
        numAlternatives <- length(alternatives[[1]])
        alternativeIndexes <- permute(1:numAlternatives)
        
        #First, try the whole condition removal. Then, try their individual removals
        if (numConditions > 1){
          auxBranch <- currentBranch[-it_Conditions]
          auxCoding <- c(currentCoding, paste(featureName,sep=""))
          auxCoding <- sort(auxCoding)
          countStates <- countStates + 1
          
          if (!isInList(auxCoding, explored)){
            explored<-c(explored,list(auxCoding))
            newValues<-computeCovPostCases_confidence(auxBranch, dataset, originalDataset, goalFeature, goalValue)
            
            improved <- !is.na(newValues[guidingMetric]) & (newValues[guidingMetric] >= currentValues[guidingMetric]) &
              !is.na(newValues["confidence"]) & (newValues["confidence"] >= currentValues["confidence"]) # When removing the whole condition, the confidence might deteriorate. That is why this comprobation
            if (improved){
              exploredAndBetter<-c(exploredAndBetter,list(auxCoding))
              anyImprovement <- TRUE
              simplifyBranch_addState(queueState, queueAcValues, queueCoding, queueDepth, auxBranch, newValues, auxCoding, currentDepth)
              numSubRules <- numSubRules + 1
              improvedAtLeastOnce <- TRUE
              
              if (numSubRules >= max(1,(maxSubRulePerRule - currentDepth))){
                break()
              }
            }
          } else{
            allHits <- allHits + 1
              if (isInList(auxCoding, exploredAndBetter)){
                improvedAtLeastOnce <- TRUE
                hits <- hits + 1
              }
          }
        }
        
        #No improvement with whole condition removal, so try each alternative
        if (improvedAtLeastOnce == FALSE){
          #If there are more than one categorical alternative, try their individual removals.
          if (numAlternatives > 1){
            
            i <- alternativeIndexes[1]
            
            #For each categorical alternative, try its removal
            for (i in alternativeIndexes){
              countStates <- countStates + 1
              auxCondition<-paste(featureName, "=", paste(alternatives[[1]][-i], collapse= ","), sep = "")
              auxBranch <- currentBranch
              auxBranch[it_Conditions] <- auxCondition
              auxCoding <- c(currentCoding, paste(featureName,"!=",alternatives[[1]][i],sep=""))
              auxCoding <- sort(auxCoding)
              
              if(!isInList(auxCoding, explored)){
                explored<-c(explored,list(auxCoding))
                newValues<-computeCovPostCases_confidence(auxBranch, dataset, originalDataset, goalFeature, goalValue)
                
                improved <- newValues[guidingMetric] >= currentValues[guidingMetric]
                if (improved){
                  exploredAndBetter <- c(exploredAndBetter, list(auxCoding))
                  anyImprovement <- TRUE
                  simplifyBranch_addState(queueState, queueAcValues, queueCoding, queueDepth, auxBranch, newValues, auxCoding, currentDepth)
                  numSubRules <- numSubRules + 1
                  improvedAtLeastOnce <- TRUE
                  
                  if (numSubRules >= max(1,(maxSubRulePerRule - currentDepth))){
                    break()
                  }
                }
              } else {
                allHits <- allHits + 1
                if (isInList(auxCoding, exploredAndBetter)){
                  improvedAtLeastOnce <- TRUE
                  hits <- hits + 1
                }
              }
            }
            
            if (numSubRules >= max(1,(maxSubRulePerRule - currentDepth))){
              break()
            }
          } else if (numConditions > 1) { ## I think that this else is now never executed, due to I have moved the whole condition removal intent at the beginning
            auxBranch <- currentBranch[-it_Conditions]
            auxCoding <- c(currentCoding, paste(featureName,sep=""))
            auxCoding <- sort(auxCoding)
            countStates <- countStates + 1
            
            if (!isInList(auxCoding, explored)){
              explored<-c(explored,list(auxCoding))
              newValues<-computeCovPostCases_confidence(auxBranch, dataset, originalDataset, goalFeature, goalValue)
              
              improved <- !is.na(newValues[guidingMetric]) & (newValues[guidingMetric] >= currentValues[guidingMetric]) & 
                !is.na(newValues["confidence"]) & (newValues["confidence"] >= currentValues["confidence"]) # When removing the whole condition, the confidence might deteriorate. That is why this comprobation
              if (improved){
                exploredAndBetter<-c(exploredAndBetter,list(auxCoding))
                anyImprovement <- TRUE
                simplifyBranch_addState(queueState, queueAcValues, queueCoding, queueDepth, auxBranch, newValues, auxCoding, currentDepth)
                numSubRules <- numSubRules + 1
                improvedAtLeastOnce <- TRUE
                
                if (numSubRules >= max(1,(maxSubRulePerRule - currentDepth))){
                  break()
                }
              }
            } else {
              allHits <- allHits + 1
              if (isInList(auxCoding, exploredAndBetter)){
                improvedAtLeastOnce <- TRUE
                hits <- hits + 1
              }
            }
          }
        }
      } else if (numConditions > 1){
        auxBranch <- currentBranch[-it_Conditions]
        auxCoding <- c(currentCoding, paste(featureName,sep=""))
        auxCoding <- sort(auxCoding)
        countStates <- countStates + 1
        
        if (!isInList(auxCoding, explored)){
          explored<-c(explored,list(auxCoding))
          newValues<-computeCovPostCases_confidence(auxBranch, dataset, originalDataset, goalFeature, goalValue)
          
          improved <- !is.na(newValues[guidingMetric]) & (newValues[guidingMetric] >= currentValues[guidingMetric]) & 
            !is.na(newValues["confidence"]) & (newValues["confidence"] >= currentValues["confidence"]) # When removing the whole condition, the confidence might deteriorate. That is why this comprobation
          if (improved){
            exploredAndBetter<-c(exploredAndBetter,list(auxCoding))
            anyImprovement <- TRUE
            simplifyBranch_addState(queueState, queueAcValues, queueCoding, queueDepth, auxBranch, newValues, auxCoding, currentDepth)
            numSubRules <- numSubRules + 1
            improvedAtLeastOnce <- TRUE
            
            if (numSubRules >= max(1,(maxSubRulePerRule - currentDepth))){
              break()
            }
          }
        } else {
          allHits <- allHits + 1
          if (isInList(auxCoding, exploredAndBetter)){
            improvedAtLeastOnce <- TRUE
            hits <- hits + 1
          }
        }
      }
    }
    
    if (!improvedAtLeastOnce){ #This determines that this branch has not been subject of further improvement
      allHits <- allHits + 1
      if (anyImprovement == TRUE){ #This assures that the branch being included is not the initial one
        outputRules <- c(outputRules, list(currentBranch))
        hits <- hits + 1
      }
    }
  }
  
  #print(paste("Mejoras_visitadas:", hits, "paraReglas:", length(outputRules), "nodosRevisitados:", allHits, "nodosTotales:", countStates))
  
  return(outputRules)
}


#This function returns the rules, with goalFeature = goalValue, from a CART tree. It also appends some accuracy metrics to these rules
getRulesFromFit <- function(fit, dataset, goalFeature, goalValue, justLeaves = TRUE,
                             minConfidence = 0.65, minPostCasesCovered = 0.05){
  library(stringr)
  goalValueIndex <- which(levels(dataset[[goalFeature]]) == goalValue)
  i<-rownames(fit[["frame"]])[2]
  rules <- data.frame(rule = character(),
                      length = numeric(),
                      covPostCases = numeric(),
                      covPostCases2 = numeric(),
                      support = numeric(),
                      confidence = numeric(),
                      lift = numeric(),
                      globalQuality = numeric(),
                      extLength = numeric(),
                      optimised = character(),
                      treeNode = numeric(),
                      stringsAsFactors = FALSE)
  
  for (i in rownames(fit[["frame"]])[-1]){
    
    if (justLeaves == FALSE || fit[["frame"]][i,"var"] == "<leaf>"){
      
      if (fit[["frame"]][i,"yval"] == goalValueIndex){
        branch <- getBranch(fit, i)
        length <- length(branch)
        accValues<-computeAccuracyValues(branch, dataset, dataset, goalFeature, goalValue)
        
        if (accValues["covPostCases"] > 0){ #In some cases, the rule covers no positives, so it should be discarded
          
          if (accValues["confidence"] >= minConfidence &&
              accValues["covPostCases"] >= minPostCasesCovered){
            newRule<-c(list(rule = list(branch),
                            length = as.numeric(length)),
                       as.list(accValues), as.list("From Tree node"),
                       treeNode = as.numeric(i))
            rules[nrow(rules)+1,] <- newRule
          }
        }
      }
    }
  }
  
  return(rules)
}


#This function returns the earlier depth at which a feature is used in any condition in a CART tree
getDepthFeature <- function(fit, feature){
  depth <- trunc(log2(
    as.integer(rownames(fit[["frame"]])[min(which(fit[["frame"]][["var"]] == feature))])))
  return(depth)
}



#This function get rules from CART trees that are recursively constructed by removing the deciding features.
#This deletion produces diverse CART trees (diverse combinations of conditions) and diverse rules therefore.
#This function has the particularity of removing the features per depths (all the features at the same depth). I
#tried another version where the deciding features were removed one by one. This produced slightly different results
#but consumed more computational time, because more recursive calls were generated.
getRulesRecursiveTrees <- function(dataset, maxLength = 4, justLeaves = TRUE, minConfidence = 0.7,
                                   call = "", time1 = Sys.time(), recursion = 3, exploringDepth = 0,
                                   depth = 0, maxTime = 180,
                                   goalFeature, goalValue){
  features <- colnames(dataset)
  numFeatures <- length(features)

  rules <- data.frame(rule = character(), # 1
                      confConditions = character(), # 2
                      confReductions = character(), # 3
                      length = numeric(), # 4
                      covPostCases = numeric(), # 5
                      covPostCases2 = numeric(), # 6
                      support = numeric(), # 7
                      confidence = numeric(), # 8
                      lift = numeric(), # 9
                      globalQuality = numeric(), # 14
                      extLength = numeric(), # 16
                      iteration = character(), # 17
                      treeNode = numeric(), # 18
                      simplification = numeric(), # 19
                      oneCatAltIndex = numeric(), # 20
                      optimised = character(), # 21
                      stringsAsFactors = FALSE)
  
  fit <- getRPartTree(dataset, goalFeature, maxdepth = maxLength)
  
  usedFeatures <- fit[["frame"]][["var"]]
  usedFeatures <- usedFeatures[usedFeatures != "<leaf>" & !is.na(usedFeatures)]
  usedFeatures <- unique(usedFeatures)
  
  
  #Print progress
  consoleWidth <- getConsoleWidth()
  time2 <- Sys.time()
  cat(str_pad(paste(as.integer(difftime(time2, time1, units=c("secs"))), "s / ",
                    maxTime, " : Call: ", call, sep=""),
              consoleWidth-2, "right"))
  cat("\r")
  #print progress
  
  
  
  signMaxTime <- sign(maxTime)
  thisCallMaxTime <- abs(maxTime)
  
  if (length(usedFeatures) > 0 && recursion > 0 && (signMaxTime * as.integer(difftime(time2, time1, units=c("secs")))) < thisCallMaxTime){
    newRules<-getRulesFromFit(fit, dataset, goalFeature, goalValue,
                               justLeaves = justLeaves, minConfidence = (minConfidence - 0.05))
    
    #Include the rules in the set
    if(nrow(newRules) > 0){
      for (j in 1:nrow(newRules)){
        rule <- newRules[j,]
        rules <- addRuleToSet(rule, iteration = call, rules, dataset = dataset, goalFeature = goalFeature,
                              goalValue = goalValue)
      }
    }
    
    
    for (exploredDepth in exploringDepth:recursion){
      removingIndexes <- c()
      
      for (feature in usedFeatures){
        depthFeature <- getDepthFeature(fit, feature)
        
        if (depthFeature == exploredDepth){
          subIndex <- which(names(dataset) == feature)
          removingIndexes <- c(removingIndexes, subIndex)
        }
      }
      
      if (length(removingIndexes) > 0){
        newRules <- getRulesRecursiveTrees(dataset[,-removingIndexes], maxLength = maxLengthRules, justLeaves = TRUE,
                                           call = paste(call, exploredDepth, sep=""), time1 = time1,
                                           recursion = recursion, exploringDepth = exploredDepth, depth = depth + 1,
                                           maxTime = maxTime, goalFeature = goalFeature, goalValue = goalValue,
                                           minConfidence = minConfidence)

        if (depth == 0){
          if(nrow(newRules) > 0){
            for (j in 1:nrow(newRules)){
              rule <- newRules[j,]
              rules <- addRuleToSet(rule, iteration = rule$iteration, rules, dataset = dataset,
                                    goalFeature = goalFeature,
                                    goalValue = goalValue)
            }
          }
        } else {
          rules <- rbind(rules,newRules)
        }
      }
    }
  }
  
  if (depth == 0){
    time2 <- Sys.time()
    cat(str_pad(paste("Tree nodes: ", as.integer(difftime(time2, time1, units=c("secs"))), "s", sep=""),
                consoleWidth-2, "right"))
    cat("\n")
  }
  return(rules)
}

#This function implements the complete workflow to generate potentially interesting rules from the dataset:
#RISTRE -> SIMPLIFICATION -> OneCatRules -> ShorterRules
getRules <- function(dataset, maxLength = 4, justLeaves = TRUE, simplifyRules = 1,
                      minConfidence = 0.7, recursion = 2,
                      maxTimeRecursion = 180){
  
  features <- colnames(dataset)
  numFeatures <- length(features)
  goalFeature <- features[numFeatures]
  goalValueIndex <- which.min(table(dataset[goalFeature]))
  goalValue <- names(goalValueIndex)

  newRules <- getRulesRecursiveTrees(dataset, maxLength = maxLength, justLeaves = justLeaves, minConfidence = minConfidence,
                                     recursion = recursion, maxTime = maxTimeRecursion, goalFeature = goalFeature,
                                     goalValue = goalValue)
  
  rules <- data.frame(rule = character(), # 1
                      confConditions = character(), # 2
                      confReductions = character(), # 3
                      length = numeric(), # 4
                      covPostCases = numeric(), # 5
                      covPostCases2 = numeric(), # 6
                      support = numeric(), # 7
                      confidence = numeric(), # 8
                      lift = numeric(), # 9
                      globalQuality = numeric(), # 14
                      extLength = numeric(), # 16
                      iteration = character(), # 17
                      treeNode = numeric(), # 18
                      simplification = numeric(), # 19
                      oneCatAltIndex = numeric(), # 20
                      optimised = character(), # 21
                      stringsAsFactors = FALSE)

  
  #Include the rules in the set
  if(nrow(newRules) > 0){
    for (j in 1:nrow(newRules)){
      rule <- newRules[j,]
      rules <- addRuleToSet(rule, iteration=rule$iteration, rules,
                            dataset = dataset, goalFeature = goalFeature,
                            goalValue = goalValue)
    }
  }
  
  
  numOriginaRules <- nrow(rules)
  
  ##Optimise the rules
  if (simplifyRules > 0 && numOriginaRules > 0){
    cat("Simplifying rules\n")
    time2 <- Sys.time()
    for (iRule in 1:numOriginaRules){
      cat(paste("\r", sprintf("%3d",iRule), "/", numOriginaRules, "original rules"))
      rule <- rules[iRule,]
      branch <- rule$rule[[1]]
      opBranchs <- simplifyBranch(branch, dataset, dataset, goalFeature = goalFeature,
                                  goalValue = goalValue, maxSubRulePerRule = simplifyRules)
      
      if (length(opBranchs) > 0){
        indexSimplification <- 1
        for (iOpRule in opBranchs){
          alreadyInRules <- isInList(iOpRule, rules[[1]])
          
          if (!alreadyInRules){
            opAccValues<-computeAccuracyValues(iOpRule, dataset, dataset, goalFeature, goalValue)
            pack <- c(list(rule = list(iOpRule),
                           length = as.numeric(length(iOpRule))),
                      as.list(opAccValues), optimised = as.character("Simplified"),
                      treeNode = as.numeric(rule$treeNode))
            rules <- addRuleToSet(pack, rule$iteration, rules, simplification = indexSimplification
                                  , dataset = dataset, goalFeature = goalFeature,
                                  goalValue = goalValue)
          }
          
          indexSimplification <- indexSimplification + 1
        }
      }
    }
    time1 <- Sys.time()
    cat(paste("Simplification ( ", simplifyRules, "): ", (nrow(rules) - nrow(newRules)), " en ", as.integer(difftime(time1, time2, units=c("secs"))), "s", sep=""))
    cat("\n")
  }
  
  #Get OneCatAlternativeRules from just optimised rules
  numRules = nrow(rules)
  if (simplifyRules > 0 && numRules > 0){
    cat("Obtaining oneCatAlt rules\n")
    for (iRule in (numOriginaRules+1):numRules){
      cat(paste("\r", sprintf("%3d",(iRule - numOriginaRules)), "/", (numRules - numOriginaRules), "simplified rules"))
      rule <- rules[iRule,]
      branch <- rule$rule[[1]]
      singleCatValueRules <- getOneCatValueRules(branch, dataset, dataset,
                                                 goalFeature, goalValue)
      
      if (length(singleCatValueRules) > 0){
        indexOneCatAlt <- 1
        for (iSingleCatVRule in singleCatValueRules){
          iSingleCatVRuleValues<-computeAccuracyValues(iSingleCatVRule, dataset, dataset, goalFeature, goalValue)
          pack <- c(list(rule = list(iSingleCatVRule),
                         length = as.numeric(length(iSingleCatVRule))),
                    as.list(iSingleCatVRuleValues), optimised = as.character("SingleCatAlt"),
                    treeNode = as.numeric(rule$treeNode))
          rules <- addRuleToSet(pack, rule$iteration, rules, simplification = rule$simplification,
                                oneCatAltIndex = indexOneCatAlt, dataset = dataset, goalFeature = goalFeature,
                                goalValue = goalValue)
          indexOneCatAlt <- indexOneCatAlt + 1
        }
      }
    }
    cat("\n")
  }
  
  #Get Reduced rules from optimised rules onward
  numRules = nrow(rules)
  if (simplifyRules > 0 && numRules > 0){
    cat("Obtaining even shorter rules\n")
    for (iRule in (numOriginaRules+1):numRules){
      cat(paste("\r", sprintf("%3d",(iRule - numOriginaRules)), "/", (numRules - numOriginaRules), "optimised rules"))
      rule <- rules[iRule,]
      branch <- rule$rule[[1]]
      shorterRules <- getReducedRules(branch, dataset, dataset, goalFeature, goalValue)
      
      if (length(shorterRules) > 0){
        indexOneCatAlt <- 1
        for (iShorterRule in shorterRules){
          alreadyInRules <- isInList(iShorterRule, rules[[1]])
          
          if (!alreadyInRules){
            iShorterRuleValues<-computeAccuracyValues(iShorterRule, dataset, dataset, goalFeature, goalValue)
            pack <- c(list(rule = list(iShorterRule),
                           length = as.numeric(length(iShorterRule))),
                      as.list(iShorterRuleValues), optimised = as.character("ShorterRule"),
                      treeNode = as.numeric(rule$treeNode))
            rules <- addRuleToSet(pack, rule$iteration, rules, simplification = rule$simplification,
                                  oneCatAltIndex = indexOneCatAlt, dataset = dataset, goalFeature = goalFeature,
                                  goalValue = goalValue)
            indexOneCatAlt <- indexOneCatAlt + 1
          }
        }
      }
    }
    cat("\n")
  }

  rules <- rules[rules$confidence >= minConfidence,]
  
  return(rules)
}
