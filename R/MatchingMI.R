library(lubridate)
library(tibble)
library (survRM2)
library("MatchIt")
library (glmnet)
library(optmatch)
library(dplyr)
library(ggplot2)


#######################################################################
######                                                           ######
######                         Matching                          ######
######                                                           ######
#######################################################################

####### Propensity score matches ######
#' Get Matches and test equivalence
#'
#' Get Matches based on formula and propensity score matching
#' @param formula formula for matches
#' @param data.registry eligible censoring historical data for matching
#' @param data.rct Randomized Clinical Trial data
#' @param method method for propensity score matching; nearest by default
#' @param ratio match ratio; 1 by default
#' @param distance propensity score estimation method; glm by default; eg, lasso, randomforest ect.
#' @param id.var ID variable 
#' @param stime.var survival time variable 
#' @param status.var status variable 
#' @param alpha alpha ratio, default as 0.2
#' @param plot F as no plot required; F by default; options: dist, caliper, spagetti, forest
#' @return equivalence test rmst and hr 
#' @return matched pairs and propensity score
#' @return matched dataset
#' @return plot-propensity score distribution, caliper distance, spagetti plot and forest plot
#' @examples 
#' get.Matched (formula = as.formula (GRP ~ log (DISDUR) + as.factor(SEX) + BMI), data.registry, data.rct, method = "nearest", ratio = 1, distance = "glm", id.var = "ID", stime.var = "STIME_Death", status.var = "STATUS_Death", alpha = 0.2, plot = F)
#' get.Matched (formula = as.formula (GRP ~ log (DISDUR) + BMI), Registry = Registry, data.registry, data.rct, plot = "dist")
#' @export

get.Matched <- function (formula, 
                         data.registry, 
                         data.rct,
                         method = "nearest",
                         ratio = 1,
                         distance = "glm",
                         id.var = "ALSNR",
                         stime.var = "STIME_Death",
                         status.var = "STATUS_Death",
                         alpha = 0.2,
                         plot = F,
                         ...){
  
  ## Prepare data
  data.registry$TRT <- 0
  data.registry$GRP <- 0
  data.rct$GRP <- 1
  
  ## Get unique variables in formula
  vars <- c (all.vars (as.formula (formula)), "TRT", id.var, stime.var, status.var)
  
  if ((all (vars %in% colnames (data.registry)) == F) | (all (vars %in% colnames (data.rct) == F))){
    stop ("Not all variables are present in data registry or RCT")
  }
  
  ## Combine datasets for matching
  data <- rbind (data.rct[, vars],
                 data.registry[, vars])
  rownames (data) <- 1:nrow (data)
  
  ## MatchIt package function:
  match <- matchit (formula = formula, 
                    data = data, 
                    method = method, 
                    distance = distance,
                    ratio = ratio)
  
  ## Define matches
  match.ii <- names (match$match.matrix[, 1])
  
  ## Create dataset
  R <- data.frame (TrialID = data[match.ii, id.var],
                   TrialPS = match$distance[match.ii])
  
  for (i in 1:ratio){
    
    match.jj <- match$match.matrix[, i]
    
    R <- cbind (R,
                data.frame (MatchID = data[match.jj, id.var],
                            MatchPS = match$distance[match.jj]))

    if (i > 1){
      names(R)[ncol(R)-1] <- paste0 ("MatchID", i)
      names(R)[ncol(R)] <- paste0 ("MatchPS", i) 
    } else {
      names(R)[ncol(R)-1] <- "MatchID"
      names(R)[ncol(R)] <- "MatchPS"
    }
  }
  
  ## Plot functions: 
  if(plot == "dist"){
    
    # plot of distribution PS score
    PS.trl <- R$TrialPS
    PS.con <- R$MatchPS
    hist(PS.trl, col = rgb(1,0,0,0.5), xlab="Probability of patients", 
         ylab="Percentage of patients", main="Distribution of propensity scores" )
    
    hist(PS.con, col = rgb(0,0,1,0.5), add = TRUE)
    legend("topright", legend=c("Randomized","Matched external"), col=c(rgb(1,0,0,0.5), 
                                                          rgb(0,0,1,0.5)), pt.cex=2, pch=15 )
    if(ratio > 1) {
      
      greys <- grep("gr[ea]y", colours(), value = TRUE)
      
      for (i in 2:ratio){
        
        hist (paste("PS.con", i, sep = ""), col = greys[i], add = TRUE)
        legend("topright", legend=c("Randomized","Matched external", "Other matched external"), col=c(rgb(1,0,0,0.5), 
                                                                            rgb(0,0,1,0.5), greys[1]), pt.cex=2, pch=15 )
      }
    }  
   }else if (plot == "caliper"){
    # plot of ordered caliper distances
    hist(abs(PS.trl-PS.con), xlab="Caliper distance", 
         ylab="Percentage of patients", main="Caliper distance between pairs" )
  
    }else if (plot == "spagetti"){
    # plot of matches & segments [spagetti]
    PS.trl <- R$TrialPS
    PS.con <- R$MatchPS
    seg.col <- ifelse (abs (PS.trl - PS.con) > 0.05, "darkred",
                       ifelse (abs (PS.trl - PS.con) > 0.01, "darkorange", "darkgreen"))
  
    plot (NULL, xlim = c (0.75, 2.25), ylim = c (0,1),  xaxt = "n", xlab = "", ylab = "Propensity Score")
    mtext(c("Randomized",  "Matched external"), side=1, at=c(1, 2))
    points (rep (1, length (PS.trl)), PS.trl, pch = 16, col = adjustcolor (seg.col, .5))
    points (rep (2, length (PS.con)), PS.con, pch = 16, col = adjustcolor (seg.col, .5))
    legend("topright", legend = c("Moderate", "Good", "Excellent"),
           pch = 16,col = c("darkred", "darkorange", "darkgreen"))
    
    segments (x0 = rep (1, length (PS.trl)), 
              x1 = rep (2, length (PS.con)),
              y0 = PS.trl, y1 = PS.con,
              col = adjustcolor (seg.col, .5)) # 1 trial, 2 matched  
  
    if(ratio > 1) {
    
    for (i in 2:ratio){
      
      points (rep (2, length (PS.trl)), paste("PS.con", i, sep = ""), 
              pch = 16, col = adjustcolor (seg.col, .5))
      
      segments (x0 = rep (1, length (PS.trl)), 
                x1 = rep (2, length (PS.trl)),
                y0 = PS.trl, y1 = paste("PS.con", i, sep = ""),
                col = adjustcolor (seg.col, .5)) # 1 trial, 2 matched  
      }
    }  
  }else if(plot == "forest") {
    # Forest plot of standardized differences after matching
    d.mean = c() 
    d.lower = c()
    d.upper = c()
    
    for (col in all.vars (as.formula (formula))[-1]){

      if (class(data[,col][[1]]) == "numeric") {
        r <- data[match.jj, col]/data[match.ii, col]
        d.mean = c(d.mean, mean(r[,1]))
        d.lower = c(d.lower, mean(r[,1]) + qnorm (0.025)*sd(r[,1]))
        d.upper = c(d.upper, mean(r[,1]) + qnorm (0.975)*sd(r[,1]))
          
      } else{
        d.mean = c(d.mean, mean(data[match.jj, col] == data[match.ii, col]))
        d.lower = c(d.lower, mean(data[match.jj, col] == data[match.ii, col]) + qnorm (0.025)*sd(data[match.jj, col] == data[match.ii, col]))
        d.upper = c(d.upper, mean(data[match.jj, col] == data[match.ii, col]) + qnorm (0.975)*sd(data[match.jj, col] == data[match.ii, col]))
      }
    }
    
    d <- data.frame( mean = d.mean,
                     lower = d.lower,
                     upper = d.upper,
                     variables = all.vars (as.formula (formula))[-1]
    )
      forest <- d|> ggplot(aes(y = variables)) + theme_classic() +
      geom_point(aes(x = mean), shape=15, size=3) +
      geom_linerange(aes(xmin=lower, xmax=upper)) +
      geom_vline(xintercept = 1, linetype="dashed") +
      labs(x="Matched external/Randomized", y="")
      print(forest)
  }
  
  ## Testing both RMST & HR
  vars.test <- c (stime.var, status.var, "GRP")
  test.data <- rbind (data[match.jj, ],
                      data[data$GRP == 1 & data$TRT == 0, ]) # randomised placebo patients
  
  ## Models:
  m <- coxph (as.formula (paste0 ("Surv (", stime.var, ", ", status.var, ") ~ GRP")), data = test.data)
  rm <- rmst2 (time = test.data[, stime.var][[1]],
               status = test.data[, status.var][[1]],
               arm = test.data$GRP, alpha = 0.2)
  
  ## Test results:
  Test.Results <- data.frame (hr = as.numeric (coef (m)),
                              hr.se = as.numeric (sqrt (vcov (m))),
                              md = as.numeric (rm$unadjusted.result[1,][1]),
                              md.se = (as.numeric (rm$unadjusted.result[1,][3]) - as.numeric (rm$unadjusted.result[1,][1]))/qnorm (0.975),
                              
                              rmst.trial = rm$RMST.arm1$result["RMST", "Est."],
                              rmst.se.trial = rm$RMST.arm1$result["RMST", "se"],
                              
                              rmst.registry = rm$RMST.arm0$result["RMST", "Est."],
                              rmst.se.registry = rm$RMST.arm0$result["RMST", "se"],
                              
                              n_external = table (test.data$GRP)[1],
                              n_randomised = table (test.data$GRP)[2],
                              n_allexternal = nrow (data.registry))
  
  ## Equivalence test:
  Equivalence <- data.frame (rmst.trial = Test.Results$rmst.trial,
                             rmst.trial.lb = Test.Results$rmst.trial - qnorm (1 - (alpha/2))*Test.Results$rmst.se.trial,
                             rmst.trial.ub = Test.Results$rmst.trial + qnorm (1 - (alpha/2))*Test.Results$rmst.se.trial,
                             
                             rmst.registry = Test.Results$rmst.registry,
                             rmst.registry.lb = Test.Results$rmst.registry - qnorm (1 - (alpha/2))*Test.Results$rmst.se.registry,
                             rmst.registry.ub = Test.Results$rmst.registry + qnorm (1 - (alpha/2))*Test.Results$rmst.se.registry)
  Equivalence$equal <- dplyr::between (x = Equivalence$rmst.trial, left = Equivalence$rmst.registry.lb, right = Equivalence$rmst.registry.ub) &
    dplyr::between (x = Equivalence$rmst.registry, left = Equivalence$rmst.trial.lb, right = Equivalence$rmst.trial.ub) 
  
  ## Kaplan meier curves
  km <- survfit (as.formula (paste0 ("Surv (", stime.var, ", ", status.var, ") ~ GRP")), data = test.data)
  
  ## Treatment comparisons
  
  # Primary models
  m.prim <- coxph (as.formula (paste0 ("Surv (", stime.var, ", ", status.var, ") ~ TRT")), data = data.rct)
  m.augment <- coxph (as.formula (paste0 ("Surv (", stime.var, ", ", status.var, ") ~ TRT")),
                      data = rbind (data[data[, id.var] %in% R$MatchID, ],
                                    data[data$GRP == 1, ]))
  m.nomatch <- coxph (as.formula (paste0 ("Surv (", stime.var, ", ", status.var, ") ~ TRT")),
                      data = data)
  
  ## Results
  TRT.Results <- data.frame (log.hr.primary = as.numeric (coef (m.prim)),
                             log.hr.se = as.numeric (sqrt (vcov (m.prim))),
                             hr.primary = exp (coef (m.prim)),
                             hr.primary.lb = exp (confint (m.prim))[1],
                             hr.primary.ub = exp (confint (m.prim))[2],
                             pval.primary = drop1 (m.prim, test = "Chisq")["TRT", "Pr(>Chi)"], 
                             
                             log.hr.augmented = as.numeric (coef (m.augment)),
                             log.hr.se.augmented = as.numeric (sqrt (vcov (m.augment))),
                             hr.augmented = exp (coef (m.augment)),
                             hr.augmented.lb = exp (confint (m.augment))[1],
                             hr.augmented.ub = exp (confint (m.augment))[2],
                             pval.augmented = drop1 (m.augment, test = "Chisq")["TRT", "Pr(>Chi)"],
                             
                             log.hr.nomatch = as.numeric (coef (m.nomatch)),
                             log.hr.se.nomatch = as.numeric (sqrt (vcov (m.nomatch))),
                             hr.nomatch = exp (coef (m.nomatch)),
                             hr.nomatch.lb = exp (confint (m.nomatch))[1],
                             hr.nomatch.ub = exp (confint (m.nomatch))[2],
                             pval.nomatch = drop1 (m.nomatch, test = "Chisq")["TRT", "Pr(>Chi)"])
  
  ## Create output files:
  L <- list (Test.Results,
             Equivalence,
             TRT.Results,
             km,
             R,
             rbind (data[data[, id.var] %in% R$MatchID, ],
                    data[data$GRP == 1, ]))
  
  names (L) <- c ("Comparison", "Equivalence", "Treatment", "Curves", "MatchSet", "Dataset")
  return (L)
  
}
