#############AUTHORSHIP INFORMATION########
# author: Jose Antonio Delgado Osuna and Carlos Garcia-Martinez
# Affiliation: University Hospital Reina Sofia, University of Córdoba
# email: jose.delgado.osuna@gmail.com, cgarcia@uco.es
# init date: Jul 2019

##############OBJECTIVE####################
# This file provides some basic functions to work with the conditions provided by
# a CART tree. In fact, most of them works with a vector structure with a condition in each position as a string,
# in the format [<index>] "<nameFeature> <numerical_operator> <Value>" or [<index>] "<nameFeature> <categorical_operator> <value>(,<value>)*"
# For instance:
# [1] "P_N_IQ_MES< 1.5"
# [2] "PROCEDIMIENTO=Amputacion Abdominoperineal,Colectomia Subtotal,Colectomia Total"
# [3] "CIRUJANO=13,14,15"  

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

#This function generates a CART tree from a dataset.
getRPartTree <- function(dataset, goalFeature, maxdepth = 30, minsplit = 6){
  if (!requireNamespace("rpart", quietly = TRUE))
    install.packages("rpart")
  library(rpart)
  
  if (!requireNamespace("rpart.plot", quietly = TRUE))
    install.packages("rpart.plot")
  library(rpart.plot)
  fit <- rpart(as.formula(paste(goalFeature, "~.")),data=dataset,method="class",maxdepth = maxdepth,
               xval = 1, minsplit = minsplit)
  return(fit)
}

#This function just insert the name of the dataset in the string generated from the concatenation of conditions of a rule
addDatasetVar <- function(string, nameDatasetVar){
  strings <- strsplit(string, "&")
  output <- paste(nameDatasetVar, "$", strings[[1]][1], sep="")
  
  for (i in strings[[1]][-1]){
    output <- paste(output, " & ", nameDatasetVar, "$", trimws(i), sep="")
  }
  
  return (output)
}

#This function returns a branch of a CART tree, from root to leaf.
getBranch <- function(fit, node){
  x_aux<-path.rpart(fit, node = node, pretty = 0, print.it = FALSE)
  branch <- eval(parse(text=paste("x_aux$\"",toString(node), "\"", sep = "")))
  return(branch[-1])
}

#This function returns the name of the feature used in a condition
getFeatureFromCondition <- function(condition){
  subsubstring <- strsplit(condition, "=|<|>")
  featureName <- subsubstring[[1]][1]
  return(featureName)
}

#This function returns the set of alternatives in a categorical condition as an array
getAlternativesCatCondition <- function(condition){
  pos<-gregexpr("=", condition)
  subsubstring <- strsplit(
    substring(condition,pos[[1]]+1),
    ",")
  return(subsubstring)
}

#This function returns the operator of a numerical condition
getOperatorNumCondition <- function(condition){
  pos<-gregexpr("<|>", condition)
  return(gsub(" ","",substring(condition,pos[[1]], pos[[1]]+1), fixed= TRUE))
}

#This function returns the posible values of a categorical feature according to the dataset
getLevelsFeature <- function(feature, dataset){
  
  return(levels(dataset[[feature]]))
}

#This function indicates if a feature is categorical or it is not
isCategorical <- function(featureName, dataset){
  return (length(getLevelsFeature(featureName, dataset)) > 1)
}

#This function generates a string, almost ready to be evaluated according to different statistics metrics, from a rule represented as an array of conditions
branchToString <- function(branch, dataset){
  if (typeof(branch) == "list"){
    branch <- branch[[1]]
  }
  
  featureName <- getFeatureFromCondition(branch[1])
  
  if (length(levels(dataset[[featureName]])) < 1){
    substring <- paste(branch[1], sep = "")
  } else {
    pos<-gregexpr("=", branch[1])
    pos<-pos[[1]][[1]] #In some cases, there are more than one "=" character. In particular, I found some values for some features with the '=' character
    substring <- paste(substring(branch[1],1,pos-1), " %in% c(", sep = "")
    subsubstring <- strsplit(
      substring(branch[1],pos+1),
      ",")
    substring <- paste(substring, "\"", subsubstring[[1]][1], "\"", sep = "")
    
    for (k in subsubstring[[1]][-1]){
      substring <- paste(substring, ", \"", k, "\"", sep = "")
    }
    
    substring <- paste(substring, ")", sep = "")
  }
  
  for (j in branch[-1]){
    subsubstring <- strsplit(j, "=|<|>")
    substring <- paste(substring, " & ", sep = "")
    featureName <- subsubstring[[1]][1]
    
    if (length(levels(dataset[[featureName]])) < 1){
      substring <- substring <- paste(substring, j, sep = "")
    } else {
      pos<-gregexpr("=", j)
      pos <- pos[[1]][[1]]
      substring <- paste(substring, 
                         substring(j,1,pos-1),
                         " %in% c(", sep = "")
      subsubstring <- strsplit(
        substring(j,pos+1),
        ",")
      substring <- paste(substring, "\"", subsubstring[[1]][1], "\"", sep = "")
      
      for (k in subsubstring[[1]][-1]){
        substring <- paste(substring, ", \"", k, "\"", sep = "")
      }
      
      substring <- paste(substring, ")", sep = "")
    }
  }
  
  return(substring)  
}

#This function returns the number of distinct features used in a ruleset
numDistinctFeatures <- function(rules){
  set <- list()

  for (aRule in rules){
    for (condition in aRule){
      feature <- getFeatureFromCondition(condition)

      if (!isInList(feature, set)){
        set <- c(set, feature)
      }
    }
  }
  
  return(length(set))
}

#This function returns the distribution of the used features over the rules in rules
getFeaturesDistribution <- function(rules){
  everyFeature <- c()
  
  for (aRule in rules){
    for (condition in aRule){
      feature <- getFeatureFromCondition(condition)
      everyFeature <- c(everyFeature, feature)
    }
  }
  
  result <- quantile(summary(as.factor(everyFeature)))
  result <- round(result / length(rules) * 100, 1)
  return(rev(result))
}

#This function returns the rules that use a particular condition. This is useful in case you know that a condition makes the rule trivial
whichRulesUseCondition <- function(ruleSet, condition){
  
  output <- c()
  numRules <- length(ruleSet)
  
  for (iRule in 1:numRules){
    aRule <- ruleSet[iRule]
    
    if (isInList(condition, aRule[[1]]))
      output <- c(output, iRule)
  }
  
  return(output)
}

#This function returns the rules that use any of the provided conditions. This is useful in case you know that a condition makes the rules trivial
whichRulesUseConditions <- function(ruleSet, conditions){
  output <- c()
  
  for (condition in conditions){
    aux <- whichRulesUseCondition(ruleSet, condition)
    output <- c(output, aux)
  }
  
  return(output)
}
