library(mice)
library(lubridate)
library(tibble)
library (survRM2)
library(readxl)
library ("ggsci")
library (glmnet)
library(optmatch)
library(dplyr)
library("gsDesign")

## Formula's based on Rubin's 1981: http://fmwww.bc.edu/RePEc/bocode/c/carlin.pdf
Rubin <- function (coef, var, alpha = 0.05){
  
  # Variances
  m <- length (coef)
  Bmi <- mean (coef) # pooled estimate
  W <- mean (var) # variance
  B <- (1/(m - 1)) * sum ((coef - Bmi)^2) # between imputation variance
  Vmi <- W + ((1+ (1/m)) * B) # total variance
  
  # Statistic:
  t <- Vmi^(-0.5) * Bmi
  df <- (m - 1)*(1 + (W/((1+(1/m))* B)))^2
  
  # Put results in dataframe & return:
  RESULT <- data.frame (EST = Bmi,
                        SE = sqrt (Vmi),
                        EST.lo = Bmi + (sqrt (Vmi) * qt (1-alpha/2, df = df)),
                        EST.up = Bmi - (sqrt (Vmi) * qt (1-alpha/2, df = df)),
                        t = t,
                        df = df,
                        p = 2*pt (t, df = df, lower.tail = if (t < 0) {T} else {F}))
  
  return (RESULT)
}

## pooled Kaplan Meier curve
poolKM <- function (IMP, maxFU, byFU = 0.1,
                    varSTIME = "STIME", varSTATUS = "STATUS"){
  
  ## Define time horizon
  km <- seq (0, maxFU, by = byFU)
  
  ## Number of imputations
  impN <- length (IMP)
  
  ## Estimate survival per imputated dataset
  Surv <- lapply (1:impN, function (ii){
    
    ## Select data
    dtemp <- IMP[[ii]]
    
    ## fit survfit:
    form <- as.formula (paste0 ("Surv (", varSTIME, ", ", varSTATUS, ") ~ 1")) 
    x <- survfit (form, data = dtemp, se = T)
    
    ## Estimates: complementary log - log #The errrors have a standard extreme value-distribution or double-exponential distribution
    Si <- summary (x, time = km)$surv #along time sequences, 0, 1,..., maxFU, potential survival rate
    Qi <- log (- log (1 - Si))
    Ui <- summary (x, time = km)$std.err^2 # along time sequences, 0, 1,..., maxFU, se^2
    Ui <- Ui / (log (1 - Si) * (1 - Si))^2
    
    ## Adjust NaN
    Qi[Qi == Inf] <- NA
    Ui[is.na (Ui)] <- 0
    
    data.frame (Qi = Qi,
                Ui = Ui)
    
  })
  
  ## Pool by Rubins
  # Reference: https://sci-hub.se/10.1186/s12874-015-0048-4
  Umu <- rowMeans (sapply (Surv, function (d){d$Ui})) # mean variance
  Qmu.between <- rowMeans (sapply (Surv, function (d){
    ifelse (is.na (d$Qi), d$Qi[!is.na (d$Qi)][1], d$Qi)
  })) # mean estimate adjusted for between imputation variance (should be low and not depend on missing)
  Bmu <- (1/(impN - 1)) * rowSums (sweep (sapply (Surv, function (d){
    ifelse (is.na (d$Qi), d$Qi[!is.na (d$Qi)][1], d$Qi)
  }), 1, Qmu.between)^2) # between variance
  Vmu <- Umu + (1 + (1/impN))*Bmu
  Qmu <- rowMeans (sapply (Surv, function (d){
    ifelse (is.na (d$Qi), log (- log (1 - 0.999999)), d$Qi)
  }))
  lower <- Qmu + qnorm (0.025)*sqrt (Vmu)
  upper <- Qmu + qnorm (0.975)*sqrt (Vmu)
  
  ## Back-transform:
  lower <- 1 - exp (-exp (lower))
  upper <- 1 - exp (-exp (upper))
  Qmu <- 1 - exp (-exp (Qmu))
  
  ## Clean & return data:
  result <- data.frame (time = km,
                        surv = Qmu,
                        lower = lower,
                        upper = upper)

  ## Small adjustments:
  result$surv[result$surv >= 0.999999] <- 1
  result$lower[result$lower >= 0.999999] <- 1
  result$upper[result$upper >= 0.999999] <- 1

  result[result$surv > 0.9999, c ("surv", "lower", "upper")] <- 1

  return (result)

}

## interim data for randomized clinical trial data 
IntrimData <- function (RCT, Date){
  
  RCTtemp <- RCT
  RCTtemp <- RCTtemp[RCTtemp$DATEINCL <= Date, ]
    
  ## Restate Survival data
  RCTtemp$STATUS_Death <- ifelse (RCTtemp$STATUS_Death == 0, 0,
                                  ifelse (RCTtemp$DATEDEATH <= Date, 1, 0))
  RCTtemp$STIME_Death <- ifelse (RCTtemp$STATUS_Death == 1, 
                                ifelse (RCTtemp$DATEDEATH >= Date, 
                                 as.numeric (ymd(Date) - ymd(RCTtemp$DATEDX))/(365.25/12),
                                 RCTtemp$STIME_Death), 
                                as.numeric (ymd(Date) - ymd(RCTtemp$DATEDX))/(365.25/12))
  
  return (RCTtemp)
  
}

## get matches for each imputated data
getMatches <- function (formula, 
                        Registry, 
                        RCT, 
                        method, 
                        distance,
                        plot = F,
                        GSD = F,
                        ratio = 1,
                        original = RCT){
  
  Registry$TRT <- 0
  Registry$GRP <- 0
  RCT$GRP <- 1
  
  ## extract formular related variables
  vars <- c(all.vars(formula), "ID","STIME_Death", "STATUS_Death")
  data <- rbind (RCT[, vars], Registry[, vars])
  rownames (data) <- 1:nrow (data)
  match <- matchit (formula = formula, 
                    data = data, 
                    method = method, 
                    distance = distance,
                    ratio = ratio)
  
  ## Define matches
  match.ii <- names (match$match.matrix[, 1])
  
  R <- data.frame (TrialID = data[rownames (data) %in% match.ii,]$ID,
                   TrialPS = match$distance[match.ii])
  
  for (i in 1:ratio){
    
    match.jj <- match$match.matrix[, i]
    
    R <- cbind(R,data.frame (MatchID = data[rownames (data) %in% match.jj,]$ID,
                             MatchPS = match$distance[match.jj]))
    
    names(R)[ncol(R)-1] = paste("MatchID", i)
    names(R)[ncol(R)] = paste("MatchPS", i)
    
  }
   
  return(R)
  
}

## censor rules for historical data
censorImp <- function (d, censDate, maxFU){
  
  ## Cap events after trial end
  d$STATUS <- ifelse (is.na (d$DATEEVENT), 0,
                      ifelse (d$DATEEVENT <= as.Date (censDate), 1, 0))
  d$STIME <- ifelse (d$DATECENS <= as.Date (censDate),
                     d$DATECENS - d$DATEDX,
                     as.Date (censDate) - d$DATEDX)
  d$STIME <- as.numeric (d$STIME)/(365.25/12)
  d$STATUS <- ifelse (d$STIME > maxFU, 0, d$STATUS)
  d$STIME <- ifelse (d$STIME > maxFU, maxFU, d$STIME)
  
  ## Similar for Death
  d$STATUS_Death <- ifelse (is.na (d$DATEDEATH), 0,
                            ifelse (d$DATEDEATH <= as.Date (censDate), 1, 0))
  d$STIME_Death <- ifelse (d$DATECENS_Death <= as.Date (censDate),
                           d$DATECENS_Death - d$DATEDX,
                           as.Date (censDate) - d$DATEDX)
  d$STIME_Death <- as.numeric (d$STIME_Death)/(365.25/12)
  d$STATUS_Death <- ifelse (d$STIME_Death > maxFU, 0, d$STATUS_Death)
  d$STIME_Death <- ifelse (d$STIME_Death > maxFU, maxFU, d$STIME_Death)
  
  ## Return data
  d
  
}

## Summary table
sumTable <- function(Registry, EliRegistry) {
  df <- data.frame(varName = c(), Registry = c(), EliRegistry = c(), p = c())
  for (col in colnames(Registry)){
    
    if (class(Registry[,col][[1]]) == "character"){
      dat <- as.factor(Registry[,col][!is.na(Registry[,col])])
      r = 0
      category = 0
      for (l in levels(dat)){
        if (sum(dat==l)/length(dat) > r) {
          category = l
          r = format(round(sum(dat==l)/length(dat), 2), nsmall = 2) }
      }
      
      dat2 <- as.factor(EliRegistry[,col][!is.na(EliRegistry[,col])])
      r2 = 0
      category2 = 0
      for (l in levels(dat2)){
        if (sum(dat2==l)/length(dat2) > r2) {
          category2 = l
          r2 = format(round(sum(dat2==l)/length(dat2), 2), nsmall = 2) }
        }
        
      val <- paste(category,r,sep='-')
      val2 <- paste(category2,r2,sep='-')
      pval <- NA

      df <- rbind(df,data.frame(varName = col, Registry = val, EliRegistry = val2, p = pval))
    }
    
    else if (class(Registry[,col][[1]]) == "numeric") {
      val <- format(round(mean(Registry[,col][!is.na(Registry[,col])]), 2), nsmall = 2) 
      val2 <- format(round(mean(EliRegistry[,col][!is.na(EliRegistry[,col])]), 2), nsmall = 2) 
      pval <- ifelse(t.test(Registry[,col][!is.na(Registry[,col])], EliRegistry[,col][!is.na(EliRegistry[,col])])$p.value < 0.05, 
                     paste(format(round(t.test(Registry[,col][!is.na(Registry[,col])], EliRegistry[,col][!is.na(EliRegistry[,col])])$p.value, 3), nsmall = 3)  ,"*"), "")
      df <- rbind(df,data.frame(varName = col, Registry = val, EliRegistry = val2, p = pval))
    }
    
    else if (class(Registry[,col][[1]]) == "date") {
      val <- mean.Date(as.Date(Registry[,col][!is.na(Registry[,col])], format=c("%m-%d-%Y")))
      val2 <- mean.Date(as.Date(EliRegistry[,col][!is.na(EliRegistry[,col])], format=c("%m-%d-%Y")))
      pval <- NA
      df <- rbind(df,data.frame(varName = col, Registry = val, EliRegistry = val2, p = pval))
    }
  }
  return(df)
}