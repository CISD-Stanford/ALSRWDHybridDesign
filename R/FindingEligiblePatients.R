library(lubridate)
library(tibble)

#######################################################################
######                                                           ######
######              Finding Eligible Historical Patients         ######
######                                                           ######
#######################################################################

####### EligiblePatients ######
#' finding eligible patients in historical dataset
#'
#' Data preparation for eligible historical dataset with inclusion and exclusion criteria
#' The function can select eligible patients with none, one or multiple inclusion or exclusion criteria. 
#' If neither inclusion nor exclusion criteria is provided, the function returns the original registry dataset.
#' 
#' @param Registry historical dataset
#' @param inclusion inclusion criteria, eg, "AGE >= 18"; multiple inclusion criteria, eg, list("AGE >= 18", "AGE < 64", "TEST1 + TEST2 > 50")
#' @param exclusion exclusion criteria, eg, "DISDUR > 36"; multiple inclusion criteria, eg, list("DISDUR > 36", "AGE >= 64")
#' @param rm.na remove na in the inclusion or exclusion variables, default as F
#' @param summary statistics of submitted full registry and eligible patients, default as F
#' @return EligiblePatients 
#' @examples 
#' EligiblePatients (Registry = Registry, inclusion = list("AGE >= 18", "AGE < 64"), exclusion = "DISDUR > 36")
#' EligiblePatients (Registry = Registry, exclusion = "DISDUR > 36", rm.na = T)
#' EligiblePatients (Registry = Registry, summary = T)
#' @export

EligiblePatients <- function(Registry, 
                             inclusion=NULL, 
                             exclusion=NULL,
                             rm.na = F,
                             summary = F){
  
  ## Eligibility criteria
  
  Eligible <- function(Registry, inclusion, exclusion){
    
    selection <- function (data, criteria){
      
      return (lapply(criteria, function(x){
        
        ## Not remove NA
        if(rm.na == T) {
          
          ## Conditions on two variables
          if(strsplit(x," ")[[1]][2] == "+" | strsplit(x," ")[[1]][2] == "-" ){
            
            eval(parse(text=paste("Registry$",strsplit(x," ")[[1]][1], strsplit(x," ")[[1]][2],
                                  "Registry$",strsplit(x," ")[[1]][3],strsplit(x," ")[[1]][4], 
                                  strsplit(x," ")[[1]][5], sep = '')))
          
            }else{
            
              eval(parse(text=paste("Registry$",x,sep = ''))) 
            }
        
          }else{
          
            ## Conditions on two variables
            if(strsplit(x," ")[[1]][2] == "+" | strsplit(x," ")[[1]][2] == "-" ){
              
              eval(parse(text=paste("Registry[!is.na(Registry$",strsplit(x," ")[[1]][1], "),]$", 
                                    strsplit(x," ")[[1]][1],strsplit(x," ")[[1]][2], 
                                    "Registry[!is.na(Registry$",strsplit(x," ")[[1]][3], "),]$",
                                    strsplit(x," ")[[1]][3], strsplit(x," ")[[1]][4], strsplit(x," ")[[1]][5], sep = '')))
              
            }else{
              
              eval(parse(text=paste("Registry[!is.na(Registry$",strsplit(x," ")[[1]][1], "),]$", x,sep = '')))
              
            }
          
            }
        
      }))
    }
    
    # intersect
    inc <- selection(Registry, inclusion)
    exc <- selection(Registry, exclusion)
    
    if(!is.null(inclusion)){
      Eligibility <- inc[[1]]
      
      }else if(!is.null(exclusion)){
        Eligibility <- 1-exc[[1]]
      }
    
    if(!is.null(inclusion)){
      for (i in 1:length(inc)){
        Eligibility <- Eligibility*inc[[i]]
      }
    }
    
    if(!is.null(exclusion)){
      for (i in 1:length(exc)){
        Eligibility <- Eligibility * (1-exc[[i]])
      }
    }
    
    Eligibility[is.na(Eligibility)]=0
    
    return (as.logical(Eligibility))
    
  }

  ## Select eligible patients
    
  if((is.null(inclusion)) && (is.null(exclusion))){
      
    if(summary == F){
      
      return(Registry)
      
      } else {
         
         df <- sumTable(Registry, Registry)
         return(df)
         
    }
      
  }else{
          
    EliIndex <- Eligible (Registry, inclusion, exclusion)
    EliRegistry <- Registry[EliIndex, ]
          
    if(summary == F){
          
          return(EliRegistry)
      
          } else {
       
          df <- sumTable(Registry, EliRegistry)
          return(df)
          
          }
  }
}


