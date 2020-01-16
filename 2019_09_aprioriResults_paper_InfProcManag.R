#############AUTHORSHIP INFORMATION########
# author: Jose Antonio Delgado Osuna and Carlos Garcia-Martinez
# Affiliation: University Hospital Reina Sofia, University of Córdoba
# email: jose.delgado.osuna@gmail.com, cgarcia@uco.es
# init date: Jul 2019

##############OBJECTIVE####################
# This file provides the source of the apriori experiments carried out for the work published in
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

#Read dataset
library(methods)
if (!requireNamespace("foreign", quietly = TRUE))
  install.packages("foreign")
library(foreign)
if (!requireNamespace("arules", quietly = TRUE))
  install.packages("arules")
suppressWarnings(library(arules))
source('src/utils.R')

#This function reads the dataset provided its filename, and discretizes it for the application of apriori.
#It is specific of the Colorectal Cancer dataset, so that is the reason of specific settings on particular features
getDataset <- function(filename){
  suppressWarnings(dataset <- read.arff(filename))
  uninterestingFeatures <- getUninterestingFeatures(filename)
  uninterestingFeatures2 <- getSecondUninterestingFeatures(filename)
  dataset[uninterestingFeatures] <- NULL
  dataset[uninterestingFeatures2] <- NULL
  dataDiscretized <- discretizeDF(dataset, met = 
                                    list(CEA=list(method="interval"),
                                         P_GLASGOW=list(method="interval"),
                                         P_N_IQ_MES=list(method="interval"),
                                         GANGLIOS_AFECTADOS=list(method="interval"),
                                         INVAS_MESORR_MM=list(method="interval")))
  return(dataDiscretized)
}

#This function print some information of the shortest rules from the given set
printShortestRules <- function(rules){
  
  if (length(rules) > 0){
    rules_info <-
      data.frame(
        lengthLHS = size(lhs(rules)),
        LHS = labels(lhs(rules)), 
        # RHS = labels(rhs(rules)),          
        quality(rules)
      )
    print(length(rules))
    print(head(rules_info[order(rules_info$lengthLHS), ], n=43))
    print(table(rules_info$lengthLHS))
  }
}

#This function computes and print the usage distribution of the features appearing in the rules
attrsDistribution <- function(rules){
  nonZeros<-which(rowSums(lhs(rules)@data) > 0)
  print(length(nonZeros))
  # summary(sort(rowSums(lhs(rules)@data)[nonZeros]))
  print(summary(sort(rowSums(lhs(rules)@data)[nonZeros])/length(rules)))
  orden<-order(rowSums(lhs(rules)@data)[nonZeros], decreasing = TRUE)[1:20]
  print(colnames(lhs(rules))[nonZeros][orden])
  print(round(sort(rowSums(lhs(rules)@data)[nonZeros], decreasing = TRUE)[1:5]/length(rules)*100))
}

#This function computes the percentage of patterns that are not covered by any rule
uncoveredPatterns<-function(rules,dataset){
  if (length(rules) > 0){
    x<-is.subset(rules, as(dataset,"transactions"))
    return((nrow(dataset) - sum(apply(x,2,any)))/nrow(dataset))
  } else {
    return(0)
  }
}

#This function print some informatio of the rules obtained
printResults<-function(rules,dataset){
  printShortestRules(rules[!is.redundant(rules)])
  attrsDistribution(rules[!is.redundant(rules)])
  uncoveredPatterns(rules[!is.redundant(rules)], dataset)
}

#This function returns the name of the feature used in a condition of a rule
getAttName<-function(condition){
  pos <- as.integer(gregexpr(pattern ='=',condition)[[1]][[1]])
  if (pos > 0)
    name <- substr(condition,1,pos-1)
  else
    name <- condition
  
  return(name)
}

#This function runs the apriori experiment provided the values of support, confidence and maxlength of the rules
testCase<-function(dataset, support, confidence, maxlen, verbose, goalFeature, goalValue){
  rules <- apriori(dataset,
                   parameter=list(supp=support, originalSupport = FALSE, ext = TRUE, conf=confidence, target="rules",
                                  maxlen=maxlen, maxtime=3600),
                   appearance = list(rhs=c(paste(goalFeature,"=",goalValue,sep=""))), control = list(verbose=verbose))
  numRules <- length(rules)
  
  if (numRules > 0){
    nonRedundant <- rules[!is.redundant(rules)]
    numNonRedundant <- length(nonRedundant)
    tableNumAttrs <- table(size(lhs(nonRedundant)))
    nonZeros<-which(rowSums(lhs(nonRedundant)@data) > 0)
    numTotalAttributes <- length(rowSums(lhs(nonRedundant)@data))
    usedConds<-length(nonZeros)
    nameUsedAttrs <-colnames(lhs(nonRedundant))[nonZeros]
    nameUsedAttrs <- unlist(lapply(nameUsedAttrs, getAttName))
    usedAttrs<-length(unique(nameUsedAttrs))
    prueba<-data.frame(cond=colnames(lhs(nonRedundant))[nonZeros], attr=nameUsedAttrs, value=rowSums(lhs(nonRedundant)@data)[nonZeros])
    
    attrsDist<-summary(sort(rowSums(lhs(nonRedundant)@data)[nonZeros], decreasing = TRUE)/length(nonRedundant))#No se usa
    attrsDist2<-sort(rowSums(lhs(nonRedundant)@data)[nonZeros], decreasing = TRUE)/length(nonRedundant)#No se usa
    attrsDist3<-summary(aggregate(prueba$value, by=list(Category=prueba$attr), FUN=sum)$x)/length(nonRedundant)
    uncoveredPatterns_value<-round(uncoveredPatterns(nonRedundant, dataset[dataset[goalFeature]==goalValue,])*100)
    minSupport<-min(quality(nonRedundant)$lhs.support)
    meanSupport<-mean(quality(nonRedundant)$lhs.support)
    minConfidence<-min(quality(nonRedundant)$confidence)
    meanConfidence<-mean(quality(nonRedundant)$confidence)
    
    auxString <- paste(
      percent(minSupport, digits = 1), ":", percent(meanSupport, digits = 1),
      " & ",
      percent(minConfidence, digits = 1), ":", percent(meanConfidence, digits = 1),
                       " & ", numNonRedundant, ' & \\', '{', sep="")
    
    if ("1" %in% names(tableNumAttrs))
      auxString <- paste(auxString, tableNumAttrs["1"],sep="")
    else
      auxString <- paste(auxString, "0", sep="")
    
    for (aux in c("2","3","4","5")){
      if (aux %in% names(tableNumAttrs))
        auxString<-paste(auxString, ",", tableNumAttrs[aux], sep="")
      else
        auxString<-paste(auxString, ",0", sep="")
    }
    
    auxString<-paste(auxString, "\\} & ", usedAttrs, "/", usedConds, "/", numTotalAttributes, " & \\{", percent(attrsDist3["Max."], digits = 0), sep="")
    
    for (aux in rev(attrsDist3[c("Min.", "1st Qu.", "Median", "3rd Qu.")])){
      auxString <- paste(auxString, ",", percent(aux, digits = 0), sep="")
    }
    
    auxString<-paste(auxString, "\\}\\% & ", uncoveredPatterns_value, "\\% \\\\\n", sep="")
    
    cat(auxString)
  } else {
    cat(paste(percent(support, digits = 1), " & ", percent(confidence, digits = 0), " & ",
              "0 & & & & \\\\\n", sep=""))
  }
  
  if (verbose)
    printResults(rules,dataset[dataset[goalFeature]==goalValue,])
  
  return(rules)
}

#This function formats a given number as a percentage
#copied from https://stackoverflow.com/questions/7145826/how-to-format-a-number-as-percentage-in-r
percent <- function(x, digits = 2, format = "f", ...) {
  paste0(formatC(100 * x, format = format, digits = digits, ...), "\\%")
}


#Experiments on complications
filename <- "Datasets/Ntto_qt_paliativo_no/CA_COLORRECTAL_COMPLICACIONES.arff"
dataset <- getDataset(filename)
features <- colnames(dataset)
numFeatures <- length(features)
goalFeature <- features[numFeatures]
goalValueIndex <- which.min(table(dataset[goalFeature]))
goalValue <- names(goalValueIndex)

cat(filename)
cat(paste0('LHS Support',
           ' & Confidence',
           # ' & Max len.',
           # ' & N.rules',
           ' & N. non-redundant',
           ' & |antecedent|',
           ' & Distinct. attrs.',
           ' & Attr. freq.',
           ' & %Uncover. pos \n'))

maxlen<-6

testCase(dataset, 0.05, 0.7, maxlen, FALSE, goalFeature, goalValue)

for (support in c(0.025))
#for (support in c(0.05,0.025)) #No rules were obtained for 0.05 and 0.9 and 0.8
  for (confidence in c(0.9,0.8,0.7))
    testCase(dataset, support, confidence, maxlen, FALSE, goalFeature, goalValue)

maxlen<-5 #Executions with maxlen<-6 and support <- 0.01 consumed more than 18GB and were interrupted
for (support in c(0.01, 0.005))
  for (confidence in c(0.9,0.8,0.7))
    testCase(dataset, support, confidence, maxlen, FALSE, goalFeature, goalValue)


#Experiments for recurrences
filename <- "Datasets/Ntto_qt_paliativo_no/CCR_SEGUIMIENTO_RECIDIVA.arff"
dataset <- getDataset(filename)
features <- colnames(dataset)
numFeatures <- length(features)
goalFeature <- features[numFeatures]
goalValueIndex <- which.min(table(dataset[goalFeature]))
goalValue <- names(goalValueIndex)

cat(filename)
cat(paste0('LHS Support',
           ' & Confidence',
           # ' & Max len.',
           # ' & N.rules',
           ' & N. non-redundant',
           ' & |antecedent|',
           ' & Distinct. attrs.',
           ' & Attr. freq.',
           ' & %Uncover. pos \n'))

maxlen<-6
for (support in c(0.05))
  for (confidence in c(0.9,0.8,0.7))
    testCase(dataset, support, confidence, maxlen, FALSE, goalFeature, goalValue)

maxlen<-5 
for (support in c(0.025,0.01, 0.005))
  for (confidence in c(0.9,0.8,0.7))
    testCase(dataset, support, confidence, maxlen, FALSE, goalFeature, goalValue)
