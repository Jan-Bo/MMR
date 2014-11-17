#' Contribution Calculations
#' 
#' Calculates the contribution
#' 
#' @param object An MMR.Optim object for which prediction is desired
#' @param combine A logical to if all predictions to be combined
#' @param combineMethod A character of {"mean", "MAPE", "MSE", "RMSE"}
#' 
#' @author Gabriel Mohanna
#' 
#' @export
#' 
MMR.Contributions <- function(MMR.object,
                              combine=TRUE,
                              combineMethod="MSE") {
    # Compute each trial contributions by variable
    trialsContrib <- lapply(MMR.object,
                            function(nlsObject){
                                trialContrib(nlsObject)
                            })
    
    # Get a list of unique coefficient names
    uniqueCoef <- unique(unlist(sapply(trialsContrib, 
                                       function(x) {
                                           if(sum(is.na(x))) {
                                               character()
                                           } else {
                                               names(x)
                                           }
                                       })))
    names(uniqueCoef) <- uniqueCoef
    
    # Get each variable contributions for each trial
    coefContrib <- lapply(uniqueCoef,
                          function(coef) {
                              sapply(trialsContrib,
                                     function(trial) {
                                         if(coef %in% names(trial)) {
                                             sum(trial[[coef]])
                                         } else {
                                             NA
                                         }
                                     })
                          })
    
    if(combine) {
        # Get weights
        if(combineMethod=="mean") {
            weights <- sapply(MMR.object,
                              function(nlsObject) {
                                  ifelse(class(nlsObject)=="nls", 1, NA)
                              })
        } else if(combineMethod=="MAPE") {
            weights <- 1 / getModelErrors(MMR.object, "MAPE", training = F)
        } else if(combineMethod=="MSE") {
            weights <- 10 / getModelErrors(MMR.object, "MSE", training = F)
        } else if(combineMethod=="RMSE") {
            weights <- 10 / sqrt(getModelErrors(MMR.object, "MSE", training = F))
        }
        
        # Get OOB indicators
        #inTrainSamples <- 
        #    mclapply(MMR.object,
        #             function(nlsObject) {
        #                 if(class(nlsObject)=="nls") {
        #                     nlsObject$model$inTrain
        #                 } else {
        #                     NA
        #                 }
        #             })
        
        #inTrainDF <- 
        #    do.call(cbind.data.frame, 
        #            lapply(inTrainSamples, 
        #                   function(x) {
        #                       if(sum(is.na(x))) {
        #                           rep(NA, max.items)
        #                       } else {
        #                           ifelse(x, 0, 1)
        #                       }
        #                   }))
                        
        sapply(coefContrib, 
               function(coef) {
                   sum(coef * weights, na.rm = T) / sum(weights, na.rm=T)
               })
    } else {
        coefContrib
    }
}


#' Contribution Calculations for One Trial
#' 
#' Calculates the contribution for one trial
#' 
#' @param object An MMR.Optim object for which prediction is desired
#' 
#' @author Gabriel Mohanna
#' 
trialContrib <- function(trial) {
    if(class(trial)!="nls") {
        NA
    } else {
        # Set up new.data along with intercept
        inTrain  <- trial$model$inTrain
        new.data <- c(list(Intercept=rep(1, length(inTrain))), 
                      getModeledData(trial))
        
        # Get model cofficients & adstocks
        coef.names     <- c("Intercept", paste0(names(new.data[-1]), ".coef"   ))
        adstocks.names <- c("Intercept", paste0(names(new.data[-1]), ".adstock"))
        
        coef.est <- sapply(coef.names,
                           function(x) coef(trial)[x], 
                           USE.NAMES = F)
        adstock.est <- sapply(adstocks.names, 
                              function(x){
                                  if(x %in% names(coef(trial))){
                                      coef(trial)[x]
                                  } else {
                                      zero <- 0
                                      names(zero) <- x
                                      zero
                                  }
                              } , 
                              USE.NAMES = F)
        adstock.est["Intercept"] <- 0
        
        # Caculate weekly contributions per variable
        mapply(function(var, coefficient, adstocks){
            coefficient*adstock(var, adstocks)
        }, new.data, coef.est, adstock.est,
        SIMPLIFY = F)
    }
}
