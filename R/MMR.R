#' Nonlinear Least Squares Parameter Estimation
#' 
#' Determines nonlinear least-squares estimates marketing-mix model
#' 
#' @param dep.var A character string for dependent variable
#' @param indp.vars A character vector of variable names to be tried
#' @param data A dataset from in which to evaluate the \code{indp.vars}
#' @param algorithm A character string specifying the algorithm to use. Only "LM" based on Levenberg-Marquardt and "port" algorithms are implemented.
#' @param cvTrials A list of logicals to indicate which variable are part of training data. MMR.Optim will repeat for every elements of \code{cvTrials}.
#' @param keep.vars An optional character vector of variable names that will be kept in each model. \code{keep.vars} must be part of \code{indep.vars}.
#' @param adstock.vars An optional character vector of variable names to be adstocked. Default min max is between 0 & 0.999. \code{adstock.vars} must be part of \code{indp.vars}.
#' @param adstock.min An optional named numeric vector of minimum adstock for each \code{adstock.vars}. Names of \code{adstock.min} must be part of \code{adstock.vars}.
#' @param adstock.max An optional named numeric vector of maximum adstock for each \code{adstock.vars}. Names of \code{adstock.max} must be part of \code{adstock.vars}.
#' @param transform.vars An optional named character vector of functional transformations for variables in \code{indp.vars}. Values can only be {'Power', 'PowerInv', 'NegExp', 'Log'}. Names of \code{transform.vars} must be part of \code{indp.vars}.
#' @param transform.min An optional named numeric vector of minimum transformation for each \code{transform.vars}.
#' @param transform.max An optional named numeric vector of maximum transformation for each \code{transform.vars}.
#' @param lags Currently not used.
#' @param indp.try A positive integer indicating the variables to be randomly selected as independent variables in each \code{cvTrials} trials. It can also be a function.
#' 
#' @details
#' MMR.Optim determines the nonlinear least-squares estimates of marketing-mix model using optimization techniques such as nls & nlsLM.
#' This approach estimates model coefficients as well as the adstock transformations and any dimensioning returns curves for each variable. This has the advantage of objectively finding the coefficients & transformations. Transformation can be modified or restricted if the user specifies optional min/max values.
#' In addition, many models can be built by using the \code{cvTrials}, which is a re-sampled data points identifier. Models are then averaged according to their out-of-sample values.
#' Finally, not all variables need to be tried at the same time. If the user specifies \code{indp.try} then the program will randomly select variables from the \code{indp.vars} for each \code{cvTrials}.
#' 
#' @author Gabriel Mohanna
#' 
#' @export
#' 
MMR.Optim <- function(dep.var,
                      indp.vars=character(),
                      data, 
                      algorithm="LM", 
                      cvTrials,
                      keep.vars=character(), 
                      adstock.vars=character(),
                      adstock.min=numeric(),
                      adstock.max=numeric(),
                      transform.vars=character(),
                      transform.min=numeric(),
                      transform.max=numeric(),
                      lags=numeric(),
                      coef.min=numeric(),
                      coef.max=numeric(),
                      indp.try =length(indp.vars)-length(keep.vars)) {
    # Error checking
    if(missing(dep.var))
        stop("Argument 'dep.var' is missing.")
    
    if(length(indp.vars)==0)
        stop("No indpendent variable(s) to try.")
    
    if(missing(data))
        stop("Argument 'data' is missing.")
    
    if(!exists(deparse(substitute(data))))
        stop("The dataset '", deparse(substitute(data)), "' doesn't exists.")
    
    if(sum(!(c(dep.var, indp.vars) %in% names(data))))
        stop("The following variable(s) are not part of '", deparse(substitute(data)), "' dataset:\n",
             paste(c(dep.var, indp.vars)[!(c(dep.var, indp.vars) %in% names(data))], collapse=", "))
    
    if(length(transform.vars)) {
        if(sum(!(names(transform.vars) %in% indp.vars))) {
            stop("The following 'transform.vars' variable(s) are not part of 'indp.vars':\n",
                 paste(names(transform.vars)[!(names(transform.vars) %in% indp.vars)], collapse=", "))
        }
        
        transform.vars <- tolower(transform.vars)
        if(sum(!(transform.vars %in% c("power", "powerinv", "log", "negexp")))) {
            stop("The following transform.varsation(s) are not recognized:\n",
                 paste(as.character(transform.vars[!(transform.vars %in% c("power", "powerinv", "log", "negexp"))]), collapse = ", "),
                 "\n\nOnly {'Power', 'PowerInv', 'Log', 'NegExp'} are allowed.", call.=F)
        }
    }
    
    # Create complete list of vars, transformations and min/max values
    # Create adstock min/max
    adstock.min  <- completeParameters(indp.vars, completeParameters(adstock.vars, adstock.min, "adstock.min", 0  ), "adstock.min", 0)
    adstock.max  <- completeParameters(indp.vars, completeParameters(adstock.vars, adstock.max, "adstock.max", .99), "adstock.max", 0)
    
    # Create transform min/max
    transform.vars <- completeParameters(indp.vars, transform.vars, "transform.vars", "Identity")
    transform.min  <- completeParameters(names(transform.vars), transform.min , "transform.min" , -Inf)
    transform.max  <- completeParameters(names(transform.vars), transform.max , "transform.max" ,  Inf)
    
    # Create lags
    lags <- completeParameters(indp.vars, lags, "lags", 0)
    
    # Craete coefficient min/max
    coef.min <- completeParameters(indp.vars, coef.min, "coef.min", -Inf)
    coef.max <- completeParameters(indp.vars, coef.max, "coef.max",  Inf)
    
    # Cross-validation trials
    if(missing(cvTrials)){
        cvTrials <- lapply(createInTrainFolds(data, k=1), function(x) !x)
    }
    
    # Run optimizations over trials & extract parameters coefficients
    optm.LM <-
        mclapply(cvTrials, 
                 FUN= function(inTrain){
                     indp.in <- indp.try
                     if(class(indp.try)=="function")
                         indp.in <- indp.try()
                     indp.in <- min(max(indp.in, 0), length(indp.vars)-length(keep.vars))
                     
                     # Run optimization
                     optmFit <- MMR.Optim.InTrain(dep.var,
                                                  indp.vars,
                                                  data, 
                                                  algorithm,
                                                  inTrain,
                                                  keep.vars,
                                                  adstock.vars,
                                                  adstock.min,
                                                  adstock.max,
                                                  transform.vars,
                                                  transform.min,
                                                  transform.max,
                                                  lags,
                                                  coef.min,
                                                  coef.max,
                                                  indp.in)})
}


#' Nonlinear Least Squares Parameter Estimation
#' 
#' Determines nonlinear least-squares estimates marketing-mix model for a single trial.
#' 
#' @param dep.var A character string for dependent variable
#' @param indp.vars A character vector of variable names to be tried
#' @param data A dataset from in which to evaluate the \code{indp.vars}
#' @param algorithm A character string specifying the algorithm to use. Only "LM" based on Levenberg-Marquardt and "port" algorithms are implemented.
#' @param inTrain A logical list indicating which observations are used for training.
#' @param keep.vars An optional character vector of variable names that will be kept in each model. \code{keep.vars} must be part of \code{indep.vars}.
#' @param adstock.vars An optional character vector of variable names to be adstocked. Default min max is between 0 & 0.999. \code{adstock.vars} must be part of \code{indp.vars}.
#' @param adstock.min An optional named numeric vector of minimum adstock for each \code{adstock.vars}. Names of \code{adstock.min} must be part of \code{adstock.vars}.
#' @param adstock.max An optional named numeric vector of maximum adstock for each \code{adstock.vars}. Names of \code{adstock.max} must be part of \code{adstock.vars}.
#' @param transform.vars An optional named character vector of functional transformations for variables in \code{indp.vars}. Values can only be {'Power', 'PowerInv', 'NegExp', 'Log'}. Names of \code{transform.vars} must be part of \code{indp.vars}.
#' @param transform.min An optional named numeric vector of minimum transformation for each \code{transform.vars}.
#' @param transform.max An optional named numeric vector of maximum transformation for each \code{transform.vars}.
#' @param lags Currently not used.
#' @param indp.try A positive integer indicating the variables to be randomly selected as independent variables in each \code{cvTrials} trials. It can also be a function.
#' 
#' @examples
#' \donttest{MMR.Optim.InTrain("sales", c("Trend", "Season"), c("ad1", "ad2", "ad3"), ad.data, sample(c(T, F), nrow(ad.data)))}
#' 
MMR.Optim.InTrain <- function(dep.var, 
                              indp.vars,
                              data, 
                              algorithm, 
                              inTrain,
                              keep.vars,
                              adstock.vars,
                              adstock.min,
                              adstock.max,
                              transform.vars,
                              transform.min,
                              transform.max,
                              lags,
                              coef.min,
                              coef.max,
                              indp.try) {
    ## QA
    #cat("indp.vars: "); cat(indp.vars)
    #cat("\nkeep.vars: "); cat(keep.vars)
    #cat("\n\n")
    
    # Determine which variable is present in data for the given inTrain sample
    inTrainSum <- rowsum(data[, indp.vars[indp.vars %in% names(data)], drop=F], inTrain)
    
    if(!("TRUE" %in% rownames(inTrainSum))) {
        present.vars <- character()
        
    } else {
        inTrainSum <- t(inTrainSum)[,"TRUE", drop=F]
        ## QA
        #print(inTrainSum)
        #print("")
        present.vars <- rownames(inTrainSum)#[inTrainSum > 0]
    }
    
    keep.vars <- present.vars[ (present.vars %in% keep.vars)]
    temp.vars <- present.vars[!(present.vars %in% keep.vars)]
    
    ## QA
    #cat("present.vars: "); cat(present.vars)
    #cat("\nkeep.vars: "); cat(keep.vars)
    #cat("\ntemp.vars: "); cat(temp.vars)
    #cat("\n\n")
    
    # Randomly select variables from present; keep the keep.vars
    select.vars <- sample(temp.vars, min(indp.try, length(temp.vars)), F)
    model.vars  <- c(keep.vars, select.vars)
    model.vars  <- model.vars[order(match(model.vars, indp.vars))]
    
    ## QA
    #cat("model.vars: "); cat(model.vars) 
    #cat("\n\n")
    
    # Get parameters starting, min & max values
    if(sum(length(model.vars))==0) {
        formula0 <- paste(dep.var, "~", "1")
    } else {
        formula0 <- paste(dep.var, "~", 
                          paste(model.vars, collapse=" + "))
    }
    
    modFit0 <- lm(formula0, data=data[inTrain, ])
    
    ## QA
    #cat("formula0:  "); cat(formula0)
    #print(modFit0)
    #cat("\n\n")
    
    coef.start <- coef(modFit0)
    coef.NA    <- names(coef.start[is.na(coef.start)])
    coef.start <- coef.start[!(names(coef.start) %in% coef.NA)]    
    coef.start.NoInt <- coef.start[names(coef.start) != "(Intercept)"]
    
    # Subset variable transformations, min & max to only include names of coef.start wo the intercept
    adstock.min  <- adstock.min [names(adstock.min ) %in% names(coef.start)]
    adstock.max  <- adstock.max [names(adstock.max ) %in% names(coef.start)]
    
    transform.vars <- transform.vars[names(transform.vars) %in% names(coef.start)]
    transform.min  <- transform.min [names(transform.min ) %in% names(coef.start)]
    transform.max  <- transform.max [names(transform.max ) %in% names(coef.start)]
    
    lags <- lags[names(lags) %in% names(coef.start)]
    
    coef.min <- coef.min[names(coef.min) %in% names(coef.start)]
    coef.max <- coef.max[names(coef.max) %in% names(coef.start)]
    
    ## QA
    #cat("coef.start:\n"); print(coef.start.NoInt)
    #cat("adstock.min:\n"); print(adstock.min)
    #cat("adstock.max:\n"); print(adstock.max)
    #cat("transform.vars:\n"); print(transform.vars)
    #cat("transform.min:\n"); print(transform.min)
    #cat("transform.max:\n"); print(transform.max)
    #cat("lags:\n"); print(lags)
    #cat("coef.min:\n"); print(coef.min)
    #cat("coef.max:\n"); print(coef.max)
    
    ## QA
    #cat("NA coef.start: \n")
    #print(coef.start[is.na(coef.start)])
    #cat("\n")
    
    # Build formula
    formula2 <- MMR.Optim.Formula(dep.var, 
                                  names(coef.start.NoInt),
                                  adstock.min,
                                  adstock.max,
                                  transform.vars,
                                  transform.min,
                                  transform.max)
    
    ## QA
    #cat("Formula:\n"); cat(formula2, "\n")
    
    coef.min <- c(Intercept=-Inf, coef.min)
    coef.max <- c(Intercept= Inf, coef.max)
    
    coef.start <- pmin(pmax(coef.start, coef.min), coef.max)
    
    names(coef.start)[1]  <- "Intercept"
    names(coef.start)     <- paste0(names(coef.start ), ".coef")
    names(coef.min  )     <- paste0(names(coef.min   ), ".coef")
    names(coef.max  )     <- paste0(names(coef.max   ), ".coef")
    
    min.max.selector <- adstock.min != adstock.max
    adstock.start    <- adstock.min[min.max.selector]
    adstock.min      <- adstock.min[min.max.selector]
    adstock.max      <- adstock.max[min.max.selector]
    
    if(length(adstock.start)) {
        names(adstock.start)  <- paste0(names(adstock.min), ".adstock")
        names(adstock.min  )  <- paste0(names(adstock.min), ".adstock")
        names(adstock.max  )  <- paste0(names(adstock.max), ".adstock")
    }
    
    transform.start <- mapply(
        function(tr.var, tr.min, tr.max) {
            tr.var <- tolower(tr.var)
            if(tr.var %in% c("identity", "log") | tr.min==tr.max) {
                NA
            } else {
                if(tr.var=="power") {
                    max(min(1, tr.max), tr.min)
                } else if(tr.var=="powerinv") {
                    min(max(0, tr.min), tr.max)
                } else if(tr.var=="negexp") {
                    min(max(.005, tr.min), tr.max)
                }
            }
        },
        transform.vars, transform.min, transform.max)
    transform.start <- transform.start[!is.na(transform.start)]
    
    #print(rbind(transform.min, transform.max))
    #print(transform.vars)
    #print(!(tolower(transform.vars) %in% c("log", "identity")));
    #print((transform.min != transform.max))
    
    min.max.selector <- !(tolower(transform.vars) %in% c("log", "identity")) & (transform.min != transform.max)
    transform.min    <- transform.min[min.max.selector]
    transform.max    <- transform.max[min.max.selector]
    
    if(length(transform.start)) {
        names(transform.start) <- paste0(names(transform.start), ".power")
        names(transform.min  ) <- paste0(names(transform.min  ), ".power")
        names(transform.max  ) <- paste0(names(transform.max  ), ".power")
    }
    
    
    ## QA
    #cat("\ncoef:\n"); print(rbind(coef.start, coef.min, coef.max))
    #cat("adstock:\n"); print(rbind(adstock.start, adstock.min, adstock.max))
    #cat("transform:\n"); print(rbind(transform.start, transform.min, transform.max))
    
    nls.start <- c(coef.start, adstock.start, transform.start)
    nls.min   <- c(coef.min  , adstock.min  , transform.min)
    nls.max   <- c(coef.max  , adstock.max  , transform.max)
    
    ## QA
    #cat("\nnls.start: \n"); print(nls.start)
    #cat("\nnls.min:   \n"); print(nls.min  )
    #cat("\nnls.max:   \n"); print(nls.max  )
    #cat("\n")
    
    data$inTrain <- inTrain
    
    ## QA
    #cat("\n-----------------------------------------------------------------------")
    #cat("\n-----------------------------------------------------------------------\n")
    
    # Run optimization
    if(algorithm=="port") {
        try(nls(formula2,
                data      = data, 
                algorithm = "port", 
                start     = nls.start, 
                lower     = nls.min, 
                upper     = nls.max,
                model     = T,
                control   = nls.control(maxiter=500)),
            silent = T)
    } else {
        try(nlsLM(formula2,
                  data      = data, 
                  algorithm = "LM", 
                  start     = nls.start, 
                  lower     = nls.min, 
                  upper     = nls.max,
                  model     = T,
                  control   = nls.lm.control(maxiter=500)),
            silent = F)
    }
}


#' Write Optim Formula
#' 
#' Writes the inTrain formula such as "sales ~ inTrain * (Trend + ad1)
#' 
#' @param dep.var A numeric vector
#' @param coef.vars A character vector
#' @param adstock.min A numeric vector
#' @param adstock.max A numeric vector
#' @param transform.vars A character vector
#' @param transform.min A numeric vector
#' @param transform.max A numeric vector
#' 
#' @author Gabriel Mohanna
#' 
MMR.Optim.Formula <- function(dep.var, 
                              coef.vars,
                              adstock.min,
                              adstock.max,
                              transform.vars,
                              transform.min,
                              transform.max) {
    # Apply adstock transformations
    adstock.list <-
        mapply(
            function(coef.name, ads.min, ads.max) {
                if(ads.min==ads.max) {
                    if(ads.max==0) {
                        coef.name
                    } else {
                        paste0("adstock(", coef.name, ", ", ads.max, ")")
                    }
                } else {
                    paste0("adstock(", coef.name, ", ", coef.name, ".adstock)")
                }
            },
            coef.vars, adstock.min, adstock.max,
            SIMPLIFY = F)
    
    # Build other transformations
    transform.list <-
        mapply(
            function(coef.name, ads.var, tr.var, tr.min, tr.max) {
                tr.var <- tolower(tr.var)
                if(tr.var=="identity") {
                    paste0(coef.name, ".coef*", ads.var)
                } else if(tr.var=="power") {
                    if(tr.min==tr.max) {
                        if(tr.max==1) {
                            paste0(coef.name, ".coef*", ads.var)
                        } else {
                            paste0(coef.name, ".coef*(", ads.var, "^", tr.max, ")")
                        }
                    } else {
                        paste0(coef.name, ".coef*(", ads.var, "^", coef.name, ".power)")
                    }
                } else if(tr.var=="powerinv") {
                    if(tr.min==tr.max) {
                        if(tr.max==0) {
                            paste0(coef.name, ".coef*", ads.var)
                        } else {
                            paste0(coef.name, ".coef/", "(", tr.max, " + 1/", ads.var, ")")
                        }
                    } else {
                        paste0(coef.name, ".coef/", "(", coef.name, ".power + 1/", ads.var, ")")
                    }
                } else if(tr.var=="log") {
                    paste0(coef.name, ".coef*", "log(1+", ads.var, ")")
                } else if(tr.var=="negexp") {
                    if(tr.min==tr.max) {
                        paste0(coef.name, ".coef*", "(1 - exp(-", tr.max, " * ", ads.var, "))")
                    } else {
                        paste0(coef.name, ".coef*", "(1 - exp(-", coef.name, ".power * ", ads.var, "))")
                    }
                }
            },
            coef.vars, adstock.list, transform.vars, transform.min, transform.max,
            SIMPLIFY = F)
    
    paste0(dep.var, " ~ ",
           "inTrain * (Intercept.coef", 
           paste(rep(" +", length(transform.list)), transform.list, collapse = ""),
           ")")
}


#' Model Prediction
#' 
#' Predicts 
#' 
#' @param object An MMR.Optim object for which prediction is desired
#' @param newdata something
#' @param combine A logical to indicate if all predictions to be combined
#' @param combineMethod A character of {"mean", "MAPE", "MSE", "RMSE"}
#' @param combineTest A logical to indicate to use inTest observations
#' 
#' @author Gabriel Mohanna
#' 
#' @export
#' 
MMR.Predict <- function(MMR.object, 
                        newdata,
                        combine=TRUE,
                        combineMethod="MSE",
                        combineTest=TRUE){
    if(!(combineMethod %in% c("mean", "MAPE", "MSE", "RMSE"))){
        stop("'combineMethod' can only be {'mean', 'MAPE', 'MSE', or 'RMSE'}")
    }
    
    useExistingData <- missing(newdata)
    
    optmFit.predict <- 
        mclapply(MMR.object,
                 FUN= function(optmFit){
                     if(class(optmFit) == "nls") {
                         if(useExistingData) {
                             new.data <- getModeledData(optmFit)
                         } else {
                             new.data <- newdata
                         }
                         
                         new.data$inTrain <- rep(T, length(new.data[[1]]))
                         
                         predict(optmFit, newdata = new.data)
                     } else {
                         NA
                     }
                 })
    
    if(combine) {
        # Convert predictions to data frame by trial
        max.items <- max(sapply(optmFit.predict, length))
        predictDF <- do.call(cbind.data.frame, 
                             lapply(optmFit.predict, 
                                    function(x) {
                                        if(sum(is.na(x))) {
                                            rep(NA, max.items)
                                        } else {
                                            x
                                        }
                                    }))
        
        # Get weights
        if(combineMethod=="mean") {
            weights <- sapply(optmFit.predict,
                              function(x) {
                                  ifelse(sum(is.na(x)), NA, 1)
                              })
        } else if(combineMethod=="MAPE") {
            weights <- 1 / getModelErrors(MMR.object, "MAPE", training = !combineTest)
        } else if(combineMethod=="MSE") {
            weights <- 10 / getModelErrors(MMR.object, "MSE", training = !combineTest)
        } else if(combineMethod=="RMSE") {
            weights <- 10 / sqrt(getModelErrors(MMR.object, "MSE", training = !combineTest))
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
        
        # Calculate weighted mean
        rowSums(predictDF * matrix(weights, nrow(predictDF), 
                                   ncol = ncol(predictDF), 
                                   byrow = T), 
                na.rm=T) / sum(weights, na.rm=T)
    } else {
        optmFit.predict
    }
}


#' Returns data used in optimization
#' 
#' Returns data used in optimization without inTrain or intercept
#' 
#' @param nlsObject An nls object.
#' 
#' @author Gabriel Mohanna
#' 
#' @export
#' 
getModeledData <- function(nlsObject) {
    new.data  <- nlsObject$model
    notDataNames <- c(names(new.data[1]), "inTrain", names(coef(nlsObject)))
    dataNames <- names(new.data[!(names(new.data) %in% notDataNames)])
    new.data[dataNames]
}


#' Prints Coefficient Estimates for MMR.Optim
#' 
#' Summarize the coefficient values for MMR.Optim
#' 
#' @param object An MMR.Optim object for which extraction of model coeeficients is meaningful
#' 
#' @author Gabriel Mohanna
#' 
#' @export
#' 
MMR.Coef <- function(object) {
    #Extract optimized parameters estimates
    optmFit.DF <- 
        mclapply(object,
                 FUN= function(optmFit) {
                     # Extract optimizated parameters cofficients
                     if(class(optmFit) == "nls") {
                         inTrain <- optmFit$model$inTrain
                         
                         new.data <- getModeledData(optmFit)
                         
                         # optmFit$model[[1]] is the first element in the model data frame and contains the dependant variable
                         new.data$inTrain <- inTrain
                         MSE.train  <- NA
                         MAPE.train <- NA
                         if(sum(inTrain)){
                             predicted  <- predict(optmFit, newdata= new.data)
                             MSE.train  <- sum(((optmFit$model[[1]] - predicted)* inTrain)^2)/sum(inTrain)
                             MAPE.train <- sum(abs(predicted/optmFit$model[[1]]-1)*inTrain)/sum(inTrain)
                         }
                         
                         new.data$inTrain <- !inTrain
                         MSE.test  <- NA
                         MAPE.test <- NA
                         if(sum(!inTrain)) {
                             predicted <- predict(optmFit, newdata= new.data)
                             MSE.test  <- sum(((optmFit$model[[1]] - predicted)*!inTrain)^2)/sum(!inTrain)
                             MAPE.test <- sum(abs(predicted/optmFit$model[[1]]-1)*!(inTrain))/sum(!inTrain)
                         }
                         
                         optmFit.return <- c(coef(optmFit), MSE.train, MSE.test, MAPE.train, MAPE.test)
                         optmFit.DF <- data.frame(Param=c(names(coef(optmFit)), "MSE.train", "MSE.test", "MAPE.train", "MAPE.test"),
                                                  optmFit=optmFit.return)
                         names(optmFit.DF)[2] <- paste0("optmFit", sample(1:1000000, 1))
                     } else {
                         optmFit.DF <- data.frame(Param=c("MSE.train", "MSE.test", "MAPE.train", "MAPE.test"), optmFit=NA)
                         names(optmFit.DF)[2] <- paste0("optmFit", sample(1:1000000, 1))
                     }
                     
                     optmFit.DF
                 })
    
    # Convert list of optimizations into data frame
    if(length(optmFit.DF)) {
        optmFit.DF <- Reduce(function(x, y) merge(x, y, by="Param", all=T), optmFit.DF)
    } else {
        optmFit.DF <- optmFit.DF$fold1
    }
    names(optmFit.DF) <- c("Param", names(object))
    
    ## Sort parameters data frame
    #coef.order <- c("Intercept", 
    #                paste0(c(fixed.vars, random.vars), ".coef"), 
    #                paste0(random.vars, ".adstock"),
    #                c("rss.train", "rss.test"))
    
    #optm.LM <- optm.LM[order(match(optm.LM$Param, coef.order)), ]
    
    # Return file optimized parameters data frame
    optmFit.DF
}


#' MMR.Coef.Plot
#' 
#' Some description
#' 
#' @param object An MMR.Optim object for which extraction of model coeeficients is meaningful
#' @param Param A character vector or one or more variables on which adstock is to be plotted
#' @param y A character of {"mean", "MAPE", "MSE", "RMSE"}
#' 
#' @author Gabriel Mohanna
#' 
#' @export
#' 
MMR.Coef.Plot <- function(object, 
                          Param, 
                          y) {
    if(missing(object))
        stop("Argument \"object\" is missing.")
    if(missing(Param))
        stop("Argument \"Param\" is missing.")
    if(!missing(y) && !(y %in% c("MAPE", "MSE", "RMSE")))
        stop("Argument \"y\" can only be 'MAPE', 'MSE', 'or RMSE'.")
    
    MMR.Coef <- MMR.Coef(object)
    
    if(sum(grepl("\\.coef$", MMR.Coef$Param))) {
        MMR.Coef <- MMR.Coef(object)
        Adstock <- melt(MMR.Coef[grepl("\\.coef$", MMR.Coef$Param), ], id.vars="Param", variable.name = "Trial", value.name = "Adstock")
        
        if(missing(y)) {
            Adstock$Param <- sub("\\.coef$", "", Adstock$Param)
            Adstock <- Adstock[order(Adstock$Param),]
            Adstock <- Adstock[complete.cases(Adstock),]
            Adstock <- Adstock[Adstock$Param %in% Param, ]
            
            unique.Param <- unique(Adstock$Param)
            Adstock$Param <- factor(Adstock$Param,levels = unique.Param[order(unique.Param, decreasing = T)])
            ggplot(Adstock, aes(Param, Adstock)) + 
                geom_boxplot() + coord_flip() + 
                xlab("Param") + ylab("Coefficient") + ggtitle("Coefficients") +
                theme_bw() + theme(plot.title = element_text(size=16, face="bold", vjust=2),
                                   axis.title.x = element_text(vjust=-0.35))
        } else {
            if(y=="MAPE") {
                MAPE.test <- melt(MMR.Coef[grepl("MAPE\\.test", MMR.Coef$Param), ], id.vars="Param", variable.name = "Trial", value.name = "MAPE.test")
            } else {
                MAPE.test <- melt(MMR.Coef[grepl("MSE\\.test" , MMR.Coef$Param), ], id.vars="Param", variable.name = "Trial", value.name = "MAPE.test")
            }
            
            if(y=="RMSE") 
                MAPE.test$MAPE.test <- sqrt(MAPE.test$MAPE.test)
            MAPE.test$MAPA <- 1/MAPE.test$MAPE.test
            
            Adstock.MAPE <- merge(Adstock, MAPE.test[,-1], by="Trial", all=T)
            Adstock.MAPE$Param <- sub("\\.coef$", "", Adstock.MAPE$Param)
            
            Adstock.MAPE <- Adstock.MAPE[order(Adstock.MAPE$Param),]
            Adstock.MAPE <- Adstock.MAPE[complete.cases(Adstock.MAPE),]
            Adstock.MAPE <- Adstock.MAPE[Adstock.MAPE$Param %in% Param, ]
            
            ggplot(Adstock.MAPE, aes(Adstock, MAPA)) + 
                geom_point(alpha=.5) + facet_grid(Param ~ .) + 
                xlab("Coefficient") + ylab(paste0("1 / ", y)) + ggtitle(paste0("Coefficients vs. 1/", y)) +
                theme_bw() + theme(plot.title = element_text(size=16, face="bold", vjust=2),
                                   axis.title.x = element_text(vjust=-0.35))
            #+ geom_vline(xintercept = Adstock.MAPE$Adstock[which.max(Adstock.MAPE$MAPA)])
        }  
        #ggplot(Adstock.MAPE, aes(MAPA, Adstock)) + geom_boxplot(aes(Adstock)) + geom_point() + coord_flip() + facet_grid(Param ~ .)
    } else {
        stop(paste("No adstock to show."))
    }
}


#' Adstock Rates Plot
#' 
#' Some description
#' 
#' @param object An MMR.Optim object for which extraction of model coeeficients is meaningful
#' @param Param A character vector or one or more variables on which adstock is to be plotted
#' @param y A character of {"mean", "MAPE", "MSE", "RMSE"}
#' 
#' @author Gabriel Mohanna
#' 
#' @export
#' 
MMR.Adstock.Plot <- function(object, 
                             Param, 
                             y) {
    if(missing(object))
        stop("Argument \"object\" is missing.")
    if(missing(Param))
        stop("Argument \"Param\" is missing.")
    if(!missing(y) && !(y %in% c("MAPE", "MSE", "RMSE")))
        stop("Argument \"y\" can only be 'MAPE', 'MSE', 'or RMSE'.")
    
    MMR.Coef <- MMR.Coef(object)
    
    if(sum(grepl("\\.adstock$", MMR.Coef$Param))) {
        MMR.Coef <- MMR.Coef(object)
        Adstock <- melt(MMR.Coef[grepl("\\.adstock$", MMR.Coef$Param), ], id.vars="Param", variable.name = "Trial", value.name = "Adstock")
        
        if(missing(y)) {
            Adstock$Param <- sub("\\.adstock$", "", Adstock$Param)
            Adstock <- Adstock[order(Adstock$Param),]
            Adstock <- Adstock[complete.cases(Adstock),]
            Adstock <- Adstock[Adstock$Param %in% Param, ]
            
            unique.Param <- unique(Adstock$Param)
            Adstock$Param <- factor(Adstock$Param,levels = unique.Param[order(unique.Param, decreasing = T)])
            ggplot(Adstock, aes(Param, Adstock)) + 
                geom_boxplot() + ylim(0, 1) + coord_flip() + 
                xlab("Param") + ylab("Adstock") + ggtitle("Adstock Rates Boxplot") +
                theme_bw() + theme(plot.title = element_text(size=16, face="bold", vjust=2),
                                   axis.title.x = element_text(vjust=-0.35))
        } else {
            if(y=="MAPE") {
                MAPE.test <- melt(MMR.Coef[grepl("MAPE\\.test", MMR.Coef$Param), ], id.vars="Param", variable.name = "Trial", value.name = "MAPE.test")
            } else {
                MAPE.test <- melt(MMR.Coef[grepl("MSE\\.test" , MMR.Coef$Param), ], id.vars="Param", variable.name = "Trial", value.name = "MAPE.test")
            }
            
            if(y=="RMSE") 
                MAPE.test$MAPE.test <- sqrt(MAPE.test$MAPE.test)
            MAPE.test$MAPA <- 1/MAPE.test$MAPE.test
            
            Adstock.MAPE <- merge(Adstock, MAPE.test[,-1], by="Trial", all=T)
            Adstock.MAPE$Param <- sub("\\.adstock$", "", Adstock.MAPE$Param)
            
            Adstock.MAPE <- Adstock.MAPE[order(Adstock.MAPE$Param),]
            Adstock.MAPE <- Adstock.MAPE[complete.cases(Adstock.MAPE),]
            Adstock.MAPE <- Adstock.MAPE[Adstock.MAPE$Param %in% Param, ]
            
            ggplot(Adstock.MAPE, aes(Adstock, MAPA)) + 
                geom_point(alpha=.5) + facet_grid(Param ~ .) + 
                ylab(paste0("1 / ", y)) + ggtitle(paste0("Adstock Rates vs. 1/", y)) + xlim(0, 1) +
                theme_bw() + theme(plot.title = element_text(size=16, face="bold", vjust=2),
                                   axis.title.x = element_text(vjust=-0.35))
            #+ geom_vline(xintercept = Adstock.MAPE$Adstock[which.max(Adstock.MAPE$MAPA)])
        }  
        #ggplot(Adstock.MAPE, aes(MAPA, Adstock)) + geom_boxplot(aes(Adstock)) + geom_point() + coord_flip() + facet_grid(Param ~ .)
    } else {
        stop(paste("No adstock to show."))
    }
}


#' Calculates Model Error
#' 
#' Some description
#' 
#' @param MMR.object An MMR.Optim object for which extraction of models error is meaningful
#' @param type A character of {"MSE", "MAPE"}
#' @param training A logical to indicate inTrain vs. inTest
#' 
#' @author Gabriel Mohanna
#' 
#' @export
#' 
getModelErrors <- function(MMR.object, type="MSE", training=T) {
    if(!(type %in% c("MSE", "MAPE")))
        stop(paste0("Argument type = '", type, "' is not recognized."))
    
    sapply(MMR.object,
           function(nlsObject) {
               errorVal <- NA
               
               if(class(nlsObject)=="nls") {
                   if(training) {
                       inTrain <- nlsObject$model$inTrain
                   } else {
                       inTrain <- !nlsObject$model$inTrain
                   }
                   new.data <- getModeledData(nlsObject)
                   
                   if(sum(inTrain)) {
                       new.data$inTrain <- inTrain
                       actual <- nlsObject$model[[1]]
                       predicted <- predict(nlsObject, newdata=new.data)
                       
                       if(type=="MSE") {
                           errorVal <- sum(((actual - predicted)*inTrain)^2)/sum(inTrain)
                       } else if(type=="MAPE") {
                           errorVal <- sum(abs(predicted / actual - 1)*inTrain)/sum(inTrain)
                       }
                   }
               }
               
               errorVal
           },
           simplify=T)
}


#' Returns a Data Frame for a List
#' 
#' Returns a Data Frame for a List of Single Vectors
#' 
#' @param mylist A list of single vecors
#' 
#' @author Gabriel Mohanna
#' 
#' @export
#' 
list2df <- function(mylist) {
    max.items <- max(sapply(mylist, length))
    
    do.call(cbind.data.frame, 
            lapply(mylist, 
                   function(x) {
                       if(sum(is.na(x))==1) {
                           rep(NA, max.items)
                       } else {
                           x
                       }
                   }))
}

#' Create a complete set of parameters
#' 
#' Some docs
#' 
#' @param indp.vars A character list of variable names
#' @param some.vars A vector (numeric or character)
#' @param some.vars.name A string indicating the name of 'some.vars'
#' @param var.init Intitial values to completed for all 'indp.vars' not in 'some.vars'
#' 
#' @author Gabriel Mohanna
#' 
completeParameters <- function(indp.vars, some.vars, some.vars.name, var.init) {
    if(missing(some.vars)) {
        some.vars <- rep(var.init, length(indp.vars))
        names(some.vars) <- indp.vars
    } else {
        if(sum(!(names(some.vars) %in% indp.vars))) {
            stop("The following '",some.vars.name, "' variable(s) are not part of 'indp.vars':\n",
                 paste(names(some.vars)[!(names(some.vars) %in% indp.vars)], collapse=", "), call.=F)
        } else if(sum(!(indp.vars %in% names(some.vars)))) {
            some.vars.add <- rep(var.init, length(indp.vars) - length(some.vars))
            #print(some.vars.name)
            #print(indp.vars); print(names(some.vars)); print(indp.vars[!(indp.vars %in% names(some.vars))]); print(!(indp.vars %in% names(some.vars))); print(some.vars.add)
            names(some.vars.add) <- indp.vars[!(indp.vars %in% names(some.vars))]
            some.vars <- c(some.vars, some.vars.add)
            some.vars <- some.vars[order(match(names(some.vars), indp.vars))]
        }
    }
    
    some.vars
}
