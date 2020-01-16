#############AUTHORSHIP INFORMATION########
# author: Jose Antonio Delgado Osuna and Carlos Garcia-Martinez
# Affiliation: University Hospital Reina Sofia, University of Córdoba
# email: jose.delgado.osuna@gmail.com, cgarcia@uco.es
# init date: Jul 2019

##############OBJECTIVE####################
#This file provides some general and useful functions

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

#This function replaces tabs per ampersands in a string
replaceTabsPerAmpersands <- function (string){
  gsub("\t"," & ",string, fixed= TRUE)
}

#This function checks if one object is or is not in a list (it implements an ordered check, so (3,2) is not the same as (2,3))
isInList <- function(element, list){
  
  for (possible in list){
    if (length(possible) != length(element))
      next
    
    if (prod(element == possible) == 1)
      return(TRUE)
  }
  return(FALSE)
}

#This function implements a copy to clipboard function
#Use example: copyToClipboard(apply(<data.frame>,1,FUN = branchToString, dataset))
copyToClipboard <- function(data){
  con<-pipe("xclip -i -selection clipboard","w")
  write.table(data, con, sep="\t")
  close(con)
}

#This function write a ruleset into a csv file.
#This function is very specific for the information we want to analyse. Its aim is just to alleviate the code of another function
writeCSV <- function(rules, filename){
  
  ##Dataframe-format for the csv output
  output <- data.frame(trivials = character(), rule = character(),
                       features = character(),
                       covPostCases = numeric(), confidence = numeric(), lift = numeric(),
                       globalQuality = numeric(),
                       withoutTrivial = character(),
                       optimised = character(),
                       iteration = character(),
                       treeNode = numeric(),
                       simplification = numeric(),
                       oneCatAltIndex = numeric(),
                       extLength = numeric(),
                       covPostCases2 = numeric(),
                       support = numeric(),
                       stringsAsFactors = FALSE)
  
  indexRow<-2
  indexM_eColumn <- letters[2 * maxLengthRules+5+1]
  
  emptyRow <- " "
  for (indexConditions in 2:maxLengthRules){
    emptyRow <- paste(emptyRow, "\t \t ", sep="");
  }
  
  for (indexRule in 1:nrow(rules)){
    rule <- rules[indexRule,]
    stringRule<-as.character(rule$rule)
    stringFeatures <- as.character(rule$features)
    numConditions<-sum(gregexpr("&", stringRule)[[1]] >= 0) + 1
    
    #Add some empty cells, to...
    {
      #skip the features
      for (i in 2:maxLengthRules){
        rule["confReductions"] <- paste(rule$confReductions, " & ", sep = "")
        rule["confConditions"] <- paste(rule$confConditions, " & ", sep = "")
      }
      
      #skip the conditions up to the maximum number of conditions
      if (maxLengthRules > numConditions){
        for (iNumConditions in 1:(maxLengthRules - numConditions)){
          stringRule <- paste(stringRule," & ", '" "', sep = "")
          stringFeatures <- paste(stringFeatures," & ", '" "', sep = "")
          #rule["rule"] <- paste(rule["rule"]," & ", '" "', sep = "")
          rule["confReductions"] <- paste(rule$confReductions, " & ", sep = "")
          rule["confConditions"] <- paste(rule$confConditions, " & ", sep = "")
        }
      }
    }
    
    
    rule["rule"]<-gsub(" & ","\t",stringRule, fixed= TRUE)
    rule["features"]<-gsub(" & ", "\t",stringFeatures, fixed= TRUE)
    # rule["rule"]<-gsub(" & ","\t",rule["rule"], fixed= TRUE)
    rule["confReductions"] <- gsub(" & ","\t",rule$confReductions, fixed= TRUE)
    rule["confConditions"] <- gsub(" & ","\t",rule$confConditions, fixed= TRUE)
    
    indexConditionLetter<-3
    indexFeatureLetter <- indexConditionLetter + maxLengthRules
    formula<-paste("=if(or(A", indexRow, " = \"R\"", sep="")
    
    for (indexConditions in 1:maxLengthRules){
      formula<-paste(formula,
                     ";isnumber(match(", letters[indexConditionLetter], indexRow, ";a$1:a$1000;0))",
                     ";isnumber(match(", letters[indexFeatureLetter], indexRow, ";a$1:a$1000;0))",
                     sep="")
      indexConditionLetter <- indexConditionLetter + 1
      indexFeatureLetter <- indexFeatureLetter + 1
    }
    formula<-paste(formula,");0;",indexM_eColumn, indexRow, ")", sep="")
    
    output <- rbind.fill(output, data.frame(trivials = NA,
                                            rule[c("rule","features","covPostCases","confidence","lift","globalQuality")],
                                            withoutTrivial = formula, rule["optimised"], iteration = rule["iteration"],
                                            treeNode = rule["treeNode"], simplification = rule["simplification"],
                                            oneCatAltIndex = rule["oneCatAltIndex"],
                                            extLength = rule["extLength"],
                                            covPostCases2 = rule["covPostCases2"], support = rule["support"]))
    output <- rbind.fill(output, data.frame(trivials = NA, rule = rule["confConditions"][[1]],
                                            withoutTrivial = formula,
                                            optimised = rule["optimised"], iteration = rule["iteration"],
                                            treeNode = rule["treeNode"], simplification = rule["simplification"],
                                            oneCatAltIndex = rule["oneCatAltIndex"]))
    output <- rbind.fill(output, data.frame(trivials = NA, rule = rule[c("confReductions")][[1]],
                                            withoutTrivial = formula,
                                            optimised = rule["optimised"], iteration = rule["iteration"],
                                            treeNode = rule["treeNode"], simplification = rule["simplification"],
                                            oneCatAltIndex = rule["oneCatAltIndex"]))
    output <- rbind.fill(output, data.frame(rule = emptyRow,
                                            withoutTrivial = formula, optimised = rule["optimised"], iteration = rule["iteration"],
                                            treeNode = rule["treeNode"], simplification = rule["simplification"],
                                            oneCatAltIndex = rule["oneCatAltIndex"]))
    
    indexRow <- indexRow + 4
  }
  
  write("trivials	empty	rule	c2\tc3\tc4\tc5\tf1\tf2\tf3\tf4\tf5\tcovPostCases	confidence	lift	globalQuality	withoutTrivial	optimised	iteration	treeNode	simplification	oneCatAltIndex	extLength	covPostCases2	support", filename)
  write.table(output, filename, sep="\t", quote = FALSE, na='""', append=TRUE, col.names = FALSE)
}

#This function adds a rule to a ruleset. It checks that the rule was not previously in the ruleset and append some accuracy values.
addRuleToSet <- function(rule, iteration, set, simplification = NA, oneCatAltIndex = NA,
                         dataset = NULL, goalFeature = NULL, goalValue = NULL){
  branch <- rule$rule[[1]]
  
  if (!isInList(branch, set[[1]])){
    confReductions<-getConfReductionConditions(branch, dataset, dataset, goalFeature, goalValue)
    confidencesConditions<-getConfidencesEachCondition(branch, dataset, dataset, goalFeature, goalValue)
    newRule <- c(list(rule = as.list(rule[["rule"]]),
                      confConditions = as.character(paste(confidencesConditions, collapse=" & ")),
                      confReductions = as.character(paste(confReductions, collapse=" & ")),
                      length = as.numeric(rule$length),
                      covPostCases = as.numeric(rule$covPostCases),
                      covPostCases2 = as.numeric(rule$covPostCases2),
                      support = as.numeric(rule$support),
                      confidence = as.numeric(rule$confidence),
                      lift = as.numeric(rule$lift),
                      globalQuality = as.numeric(rule$globalQuality),
                      extLength = as.numeric(rule$extLength),
                      iteration = as.character(iteration),
                      treeNode = as.numeric(rule$treeNode),
                      simplification = as.numeric(simplification),
                      oneCatAltIndex = as.numeric(oneCatAltIndex),
                      optimised = as.character(rule$optimised)
    ))
    
    set[nrow(set)+1,] <- newRule
  }
  
  return(set)
}

#This function outputs, to the standard output, some statistics of the rules obtained
outputInfo <- function(rules, dataset, goalFeature, goalValue){
  cat("\\hline\n")
  cat(' & LHS Supp. & Conf. & Number& $\\vert$ antecedent $\\vert$ & Distinct & Attr & \\%Uncovered \\\\\n')
  cat("Id & Min : Avg & Min : Avg & rules & \\{1,2,3,4,5\\} & used attrs & use freq. & positives\\\\\n")
  
  #Tree nodes
  cat("\\hline\n")
  cat(paste("Tree & ",
            #"-:",
            format(round(
              min(rules[rules$optimised == "From Tree node","support"]) / nrow(dataset) * 100, 1),
              # min(rules[rules$optimised == "From Tree node","covPostCases2"]) / nrow(dataset) * 100, 1),
              nsmall = 1),":",
            format(round(
              mean(rules[rules$optimised == "From Tree node","support"]) / nrow(dataset) * 100, 1),
              # mean(rules[rules$optimised == "From Tree node","covPostCases2"]) / nrow(dataset) * 100, 1),
              nsmall = 1),
            "\\% & ", sep = ""))
  cat(paste("70\\%:", format(round(mean(rules[rules$optimised == "From Tree node","confidence"]) * 100, 0), nsmall = 0), "\\% & ", sep=""))
  cat(paste(sum(rules$optimised == "From Tree node"), " & ", sep=""))
  cat(paste("\\{", sum(rules$optimised == "From Tree node" & rules$length == 1),
            ",",sum(rules$optimised == "From Tree node" & rules$length == 2),
            ",",sum(rules$optimised == "From Tree node" & rules$length == 3),
            ",",sum(rules$optimised == "From Tree node" & rules$length == 4),
            ",",sum(rules$optimised == "From Tree node" & rules$length == 5),"\\} & ",sep=""))
  cat(paste(numDistinctFeatures(rules[rules$optimised == "From Tree node", "rule"]), " & ", sep=""))
  cat(paste("\\{", paste(getFeaturesDistribution(rules[rules$optimised == "From Tree node", "rule"]), collapse=","), "\\}\\% & ", sep=""))
  cat(paste(round(
    getUncoveredPositives(rules[rules$optimised == "From Tree node", "rule"], dataset, goalFeature, goalValue) * 100),
    "\\%", sep=""))
  cat("\\\\\n")
  
  #Simplified rules
  cat("\\hline\n")
  cat(paste("Simplified & ",
            #"-:",
            format(
              min(mean(rules[rules$optimised == "Simplified","support"]) / nrow(dataset) * 100, 1),
              # min(mean(rules[rules$optimised == "Simplified","covPostCases2"]) / nrow(dataset) * 100, 1),
              nsmall = 1),":",
            format(round(
              mean(rules[rules$optimised == "Simplified","support"]) / nrow(dataset) * 100, 1),
              # mean(rules[rules$optimised == "Simplified","covPostCases2"]) / nrow(dataset) * 100, 1),
              nsmall = 1),
            "\\% & ", sep = "")
      )
  cat(paste("70\\%:", format(round(mean(rules[rules$optimised == "Simplified","confidence"]) * 100, 0), nsmall = 0), "\\% & ", sep=""))
  cat(paste(sum(rules$optimised == "Simplified"), " & ", sep=""))
  cat(paste("\\{", sum(rules$optimised == "Simplified" & rules$length == 1),
            ",",sum(rules$optimised == "Simplified" & rules$length == 2),
            ",",sum(rules$optimised == "Simplified" & rules$length == 3),
            ",",sum(rules$optimised == "Simplified" & rules$length == 4),
            ",",sum(rules$optimised == "Simplified" & rules$length == 5),"\\} & ",sep=""))
  cat(paste(numDistinctFeatures(rules[rules$optimised == "Simplified", "rule"]), " & ", sep=""))
  cat(paste("\\{", paste(getFeaturesDistribution(rules[rules$optimised == "Simplified", "rule"]), collapse=","), "\\}\\% & ", sep=""))
  cat(paste(
    round(
      getUncoveredPositives(
        rules[rules$optimised == "Simplified", "rule"],
        dataset, goalFeature, goalValue) * 100),
    ":",
    round(
      getUncoveredPositives(
        rules[rules$optimised == "From Tree node" | rules$optimised == "Simplified", "rule"],
        dataset, goalFeature, goalValue) * 100),
    "\\%",
    sep=""))
  cat("\\\\\n")
  
  #SingleCatAlt rules
  cat("\\hline\n")
  cat(paste("SingleCatAlt & ",
            #0.5:,
            format(round(
              min(rules[rules$optimised == "SingleCatAlt","support"]) / nrow(dataset) * 100, 1),
              # min(rules[rules$optimised == "SingleCatAlt","covPostCases2"]) / nrow(dataset) * 100, 1),
              nsmall = 1),":",
            format(round(
              mean(rules[rules$optimised == "SingleCatAlt","support"]) / nrow(dataset) * 100, 1),
              # mean(rules[rules$optimised == "SingleCatAlt","covPostCases2"]) / nrow(dataset) * 100, 1),
              nsmall = 1),
            "\\% & ", sep = ""))
  cat(paste("70\\%:", format(round(mean(rules[rules$optimised == "SingleCatAlt","confidence"]) * 100, 0), nsmall = 0), "\\% & ", sep=""))
  cat(paste(sum(rules$optimised == "SingleCatAlt"), " & ", sep=""))
  cat(paste("\\{", sum(rules$optimised == "SingleCatAlt" & rules$length == 1),
            ",",sum(rules$optimised == "SingleCatAlt" & rules$length == 2),
            ",",sum(rules$optimised == "SingleCatAlt" & rules$length == 3),
            ",",sum(rules$optimised == "SingleCatAlt" & rules$length == 4),
            ",",sum(rules$optimised == "SingleCatAlt" & rules$length == 5),"\\} & ",sep=""))
  cat(paste(numDistinctFeatures(rules[rules$optimised == "SingleCatAlt", "rule"]), " & ", sep=""))
  cat(paste("\\{", paste(getFeaturesDistribution(rules[rules$optimised == "SingleCatAlt", "rule"]), collapse=","), "\\}\\% & ", sep=""))
  cat(paste(
    round(
      getUncoveredPositives(
        rules[rules$optimised == "SingleCatAlt", "rule"],
        dataset, goalFeature, goalValue) * 100), ":",
    round(
      getUncoveredPositives(
        rules[rules$optimised == "From Tree node" | rules$optimised == "Simplified" | rules$optimised == "SingleCatAlt", "rule"],
        dataset, goalFeature, goalValue) * 100), "\\%",
    sep=""))
  cat("\\\\\n")

  
  #Shorter rules
  cat("\\hline\n")
  cat(paste("ShorterRules & ",
            #0.5:",
            format(round(
              min(rules[rules$optimised == "ShorterRule","support"]) / nrow(dataset) * 100, 1),
              # min(rules[rules$optimised == "ShorterRule","covPostCases2"]) / nrow(dataset) * 100, 1),
              nsmall = 1),":",
            format(round(
              mean(rules[rules$optimised == "ShorterRule","support"]) / nrow(dataset) * 100, 1),
              # mean(rules[rules$optimised == "ShorterRule","covPostCases2"]) / nrow(dataset) * 100, 1),
              nsmall = 1),
            "\\% & ", sep = ""))
  cat(paste("70\\%:", format(round(mean(rules[rules$optimised == "ShorterRule","confidence"]) * 100, 0), nsmall = 0), "\\% & ", sep=""))
  cat(paste(sum(rules$optimised == "ShorterRule"), " & ", sep=""))
  cat(paste("\\{", sum(rules$optimised == "ShorterRule" & rules$length == 1),
            ",",sum(rules$optimised == "ShorterRule" & rules$length == 2),
            ",",sum(rules$optimised == "ShorterRule" & rules$length == 3),
            ",",sum(rules$optimised == "ShorterRule" & rules$length == 4),
            ",",sum(rules$optimised == "ShorterRule" & rules$length == 5),"\\} & ",sep=""))
  cat(paste(numDistinctFeatures(rules[rules$optimised == "ShorterRule", "rule"]), " & ", sep=""))
  cat(paste("\\{", paste(getFeaturesDistribution(rules[rules$optimised == "ShorterRule", "rule"]), collapse=","), "\\}\\% & ", sep=""))
  cat(paste(
    round(
      getUncoveredPositives(
        rules[rules$optimised == "ShorterRule", "rule"],
        dataset, goalFeature, goalValue) * 100), ":",
    round(
      getUncoveredPositives(
        rules[rules$optimised == "From Tree node" | rules$optimised == "Simplified" | rules$optimised == "SingleCatAlt" | rules$optimised == "ShorterRule", "rule"],
        dataset, goalFeature, goalValue) * 100), "\\%",
    sep=""))
  cat("\\\\\n")
  
  
  cat("\\hline\n")
  cat(paste("Total & ",
            #"-:",
            format(round(
              min(rules$support) / nrow(dataset) * 100, 1),
              # min(rules$covPostCases2) / nrow(dataset) * 100, 1),
            nsmall = 1),":",
            format(round(
              mean(rules$support) / nrow(dataset) * 100, 1),
              # mean(rules$covPostCases2) / nrow(dataset) * 100, 1),
              nsmall = 1),
            "\\% & ", sep = ""))
  cat(paste("70\\%:", format(round(mean(rules$confidence) * 100, 0), nsmall = 0), "\\% & ", sep=""))
  cat(paste(nrow(rules), " & ", sep=""))
  cat(paste("\\{", sum(rules$length == 1),
            ",",sum(rules$length == 2),
            ",",sum(rules$length == 3),
            ",",sum(rules$length == 4),
            ",",sum(rules$length == 5),"\\} & ",sep=""))
  cat(paste(numDistinctFeatures(rules$rule), " & ", sep=""))
  cat(paste("\\{", paste(getFeaturesDistribution(rules$rule), collapse=","), "\\}\\% & ", sep=""))
  cat(paste(round(
    getUncoveredPositives(rules$rule, dataset, goalFeature, goalValue) * 100),
    "\\%", sep=""))
  cat("\\\\\n")
  
  cat("\\hline\n")
}

#This function returns a set of uninteresting conditions according to the used dataset. This information was provided by the expert
getUninterestingConditions <- function(nameDataset){
  if (nameDataset == "Datasets/Ntto_qt_paliativo_no/CA_COLORRECTAL_COMPLICACIONES.arff"){
    return(c("TIEMPO_HASTA_TM< 75.5", "REINGRESO=No Aplicable,Si", "TIEMPO_HASTA_RL< 78.5",
             "ESTANCIA>=12.5", "P_N_IQ_MES>=1.5", "TIEMPO_HASTA_EX_ACT< 82.5", "TIEMPO_HASTA_EX_ACT< 75.5",
             "EXITUS=TRUE", "SEG_EXITUS=TRUE"))
  } else if (nameDataset == "Datasets/Ntto_qt_paliativo_no/CCR_SEGUIMIENTO_RECIDIVA.arff"){
    return(c())
  }
  
  return(NULL)
}

#This function returns a set of uninteresting features according to the used dataset. This information was provided by the expert
getUninterestingFeatures <- function(nameDataset){
  if (nameDataset == "Datasets/Ntto_qt_paliativo_no/CA_COLORRECTAL_COMPLICACIONES.arff"){
    return(c(
      "CARCINOMATOSIS",
      "CAUSA_EXITUS",
      "CIRUGIA_CARCINOMATO",
      "CIRUGIA_COMPL",
      "CIRUGIA_MXH",
      "CIRUGIA_MXP",
      "CIRUGIA_RL",
      "CIRUGIA_TM",
      "COMPLICACION_ESTOMA",
      "ESTANCIA",
      "EXITUS",
      "FIN_SEGUIMIENTO",
      "LOCALIZACION_RL",
      "LOCALIZACION_RD",
      "LOCALIZACION_TM",
      "MX_CEREBRO",
      "MX_GANGLIOS_RETRO",
      "MX_HIGADO",
      "MX_HUESO",
      "MX_PULMON",
      "ONCOLOGO",
      "OTRAS_MX",
      "PERDIDO_SEGUIMIENTO",
      "P_N_IQ_MES",
      "REINGRESO",
      "RT_POST",
      "SEG_EXITUS",
      "SEG_LOCALIZACION",
      "TIEMPO_HASTA_EX_ACT",
      "TIEMPO_HASTA_RD",
      "TIEMPO_HASTA_RL",
      "TIEMPO_HASTA_TM",
      "TM",
      "TTO_ADYUVANTE",
      "TTO_NO_QUIRURGICO",
      "TTO_SUSPENDIDO",
      "VISTO_ONCOLOGIA",
      "ZZZ" #artificial and non-existent feature
             ))
  } else if (nameDataset == "Datasets/Ntto_qt_paliativo_no/CCR_SEGUIMIENTO_RECIDIVA.arff"){
    return(c(
      "LOCALIZACION_RD"
      ,"LOCALIZACION_RL"
      , "TIEMPO_HASTA_RD"
      , "TIEMPO_HASTA_RL"
      , "CIRUGIA_RL" #También?
      
      ,"CIRUGIA_MXH"
      , "CIRUGIA_MXP"
      , "MXH_RESECADAS"
      , "MX_HIGADO"
      , "MX_PULMON"
      , "MX_HUESO"
      , "MX_CEREBRO"
      , "MX_GANGLIOS_RETRO"
      , "MESTASTASIS"
      , "OTRAS_MX"
      , "TIEMPO_HASTA_TM"

      , "TIEMPO_HASTA_EX_ACT"
      , "SEG_EXITUS"
      , "FIN_SEGUIMIENTO"
      , "CAUSA_EXITUS"
      , "EXITUS"
      ))
  } else if (nameDataset == "Datasets/Ntto_qt_paliativo_no/CA_COLORRECTAL_MORTALIDAD.arff"){
    return(c('EXITUS'
             , 'TIEMPO_HASTA_EX_ACT'
             , 'CAUSA_EXITUS'))
  } else if (nameDataset == "Datasets/Ntto_qt_paliativo_no_AND_reseccion_si/CCR_SEGUIMIENTO_RECIDIVA.arff"){
    return(c('TIEMPO_HASTA_TM',
             'TIEMPO_HASTA_RL'))
  }
  
  return(NULL)
}

#This function returns a set of other uninteresting features according to the used dataset. This information was provided by the expert
getSecondUninterestingFeatures <- function(nameDataset){
  if (nameDataset == "Datasets/Ntto_qt_paliativo_no/CA_COLORRECTAL_COMPLICACIONES.arff"){
    return(c("ORIGEN",
             "ZZZ"))
  } else if (nameDataset == "Datasets/Ntto_qt_paliativo_no/CCR_SEGUIMIENTO_RECIDIVA.arff"){
    return(c("ZZZ"
    ))
  }
}

#This function returns the width of the console where the program is being executed    
getConsoleWidth <- function(){
  consoleWidth <- as.integer(Sys.getenv('COLUMNS'))
  if(is.na(consoleWidth))
    consoleWidth <- as.integer(options('width'))
  
  return(consoleWidth)
}

#This function write rules into a file
writeRules <- function(rules, filename){
  
  numRules <- length(rules)

  for (iRule in 1:numRules){
    aRule <- rules[[iRule]]
    numConditions <- length(aRule)
    string <- aRule[[1]]
    
    if (numConditions > 1){
      for (iCondition in 2:numConditions){
        string <- paste(string, aRule[[iCondition]], sep="\t")
      }
    }
    
    write(string, filename, append=TRUE)
  }
}

#This function read rules from a file
readRules <- function(filename){
  rules <- list()
  data <- read.delim(filename, header = FALSE, sep = "\n", stringsAsFactors = FALSE)
  numRules <- nrow(data)
  
  for (i in 1:numRules){
    conditions <- strsplit(data[i,], split="\t")
    rules <- c(rules, conditions)
  }
  
  rules
}

#This function format a number in percentage
#Function from https://stackoverflow.com/questions/7145826/how-to-format-a-number-as-percentage-in-r
percent <- function(x, digits = 2, format = "f", ...) {
  paste0(formatC(100 * x, format = format, digits = digits, ...), "\\%")
}
