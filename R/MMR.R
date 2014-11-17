#' Nonlinear Least Squares Parameter Estimation
#' 
#' Determines nonlinear least-squares estimates marketing-mix model
#' 
#' @param dep.var Dependent variable
#' @param indp.vars List of variable names to be tried
#' @param min.adstock An optional list as a numeric vector of minimum adstock for each \code{indep.vars}
#' @param max.adstock An optional list as a numeric vector of maximum adstock for each \code{indep.vars}
#' @param keep.vars An optional list of variable names that will be kept in each model. \code{keep.vars} must be part of \code{indep.vars}
#' @param transform TBD
#' @param data A dataset from in which to evaluate the \code{indp.vars}
#' @param algorithm Character string specifying the algorithm to use. Only "LM" based on Levenberg-Marquardt and "port" algorithms are implented.
#' @param cvTrials A list of logicals to indicate which variable are part of training data. MMR.Optim will repeat for every elements of \code{cvTrials}
#' @param indp.vars.try Number of variables randomly sample as independant variables in each \code{cvTrials} trials
#' 
#' @details
#' MMR.Optim determines the nonlinear least-squares estimates of marketing-mix model using optimization techniques such as nls & nlsLM. 
#' The current functionality, for now, focuses on estimating the adstock effect with future plans to also determine dimensioning returns using power transformation or something alike.
#' 
#' @author Gabriel Mohanna
#' 
#' @export
#' 
MMR.Optim <- function(dep.var,
                      indp.vars=character(),
                      min.adstock,
                      max.adstock,
                      keep.vars=character(), 
                      transform=character(),
                      data, 
                      algorithm="LM", 
                      cvTrials,
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

    if(!missing(min.adstock)) {
        if(sum(!(names(min.adstock) %in% indp.vars))) {
            stop("The following 'min.adstock' variable(s) are not part of 'indp.vars':\n",
                 paste(names(min.adstock)[!(names(min.adstock) %in% indp.vars)], collapse=", "))
        } else if(sum(!(indp.vars %in% names(min.adstock)))) {
            min.adstock.add <- rep(0, length(indp.vars) - length(min.adstock))
            names(min.adstock.add) <- indp.vars[!(indp.vars %in% names(min.adstock))]
            min.adstock <- c(min.adstock, min.adstock.add)
        }
    }
    if(!missing(max.adstock)) {
        if(sum(!(names(max.adstock) %in% indp.vars))) {
            stop("The following 'max.adstock' variable(s) are not part of 'indp.vars':\n",
                 paste(names(max.adstock)[!(names(max.adstock) %in% indp.vars)], collapse=", "))
        } else if(sum(!(indp.vars %in% names(max.adstock)))) {
            max.adstock.add <- rep(0.99, length(indp.vars) - length(max.adstock))
            names(max.adstock.add) <- indp.vars[!(indp.vars %in% names(max.adstock))]
            max.adstock <- c(max.adstock, max.adstock.add)
        }
    }
    
    if(missing(min.adstock)) {
        min.adstock <- rep(0, length(indp.vars))
        names(min.adstock) <- indp.vars
    }
    if(missing(max.adstock)){
        max.adstock <- rep(.99, length(indp.vars))
        names(max.adstock) <- indp.vars
    }
    
    if(length(transform)) {
        if(sum(!(names(transform) %in% indp.vars))) {
            stop("The following 'transform' variable(s) are not part of 'indp.vars':\n",
                 paste(names(transform)[!(names(transform) %in% indp.vars)], collapse=", "))
        }
        
        transform <- tolower(transform)
        if(sum(!(transform %in% c("power", "powerinv", "log", "negexp")))) {
            stop("The following transformation(s) are not recognized:\n",
                 paste(as.character(transform[!(transform %in% c("power", "powerinv", "log", "negexp"))]), collapse = ", "),
                 "\n\nOnly {'Power', 'PowerInv', 'Log', 'NegExp'} are allowed.")
        }
    }
    
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
                                                  min.adstock,
                                                  max.adstock,
                                                  keep.vars, 
                                                  transform,
                                                  data, 
                                                  algorithm,
                                                  inTrain,
                                                  indp.in)})
}


#' Nonlinear Least Squares Parameter Estimation
#' 
#' Determines nonlinear least-squares estimates marketing-mix model for a single trial.
#' 
#' @param dep.var Dependent variable
#' @param indp.vars List of variable names to be tried
#' @param min.adstock An optional list as a numeric vector of minimum adstock for each \code{indp.vars}
#' @param max.adstock An optional list as a numeric vector of maximum adstock for each \code{indp.vars}
#' @param keep.vars An optional list of variable names that will be kept in each model. \code{keep.vars} must be part of \code{indep.vars}
#' @param transform TBD
#' @param data A dataset from in which to evaluate the \code{indp.vars}
#' @param algorithm The algorithm to use {"LM", "port"}
#' @param inTrain A logical list indicating which observations are used for training
#' @param indp.try A positive decimal to indicate the probability of including all random.vars
#' 
#' @examples
#' \donttest{MMR.Optim.InTrain("sales", c("Trend", "Season"), c("ad1", "ad2", "ad3"), ad.data, sample(c(T, F), nrow(ad.data)))}
#' 
MMR.Optim.InTrain <- function(dep.var, 
                              indp.vars,
                              min.adstock,
                              max.adstock,
                              keep.vars, 
                              transform,
                              data, 
                              algorithm, 
                              inTrain,
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
    coef.min   <- rep(-Inf, length(coef.start))
    coef.max   <- rep( Inf, length(coef.start))
    
    ## QA
    #cat("NA coef.start: \n")
    #print(coef.start[is.na(coef.start)])
    #cat("\n")
    
    if(length(coef.start)==1) {
        names(coef.start) <- c("Intercept")
    } else {
        names(coef.start) <- c("Intercept", 
                               paste0(names(coef.start)[2:(length(coef.start))], 
                                      ".coef"))        
    }
    names(coef.min  ) <- names(coef.start)
    names(coef.max  ) <- names(coef.start)

    ## QA
    #cat("min.adstock: "); print(min.adstock)
    #cat("max.adstock: "); print(max.adstock)
    #cat("\n\n")
    
    # Split model vars into fixed.vars & random.vars
    model.vars  <- model.vars[!(model.vars %in% coef.NA)]
    fixed.vars  <- model.vars[min.adstock[model.vars]==max.adstock[model.vars]]
    random.vars <- model.vars[min.adstock[model.vars]!=max.adstock[model.vars]]
    
    # Order variables
    fixed.vars  <- fixed.vars [order(match(fixed.vars , indp.vars))]
    random.vars <- random.vars[order(match(random.vars, indp.vars))]
    
    ## QA
    #cat("fixed.vars: "); cat(fixed.vars) 
    #cat("\nrandom.vars: "); cat(random.vars) 
    #cat("\n\n")

    adstocks.start <- min.adstock[random.vars]
    adstocks.min   <- min.adstock[random.vars]
    adstocks.max   <- max.adstock[random.vars]
    
    ## QA
    #print(length(adstocks.start))
    
    if(length(adstocks.start)) {
        names(adstocks.start) <- paste0(names(adstocks.start), ".adstock")
        names(adstocks.min  ) <- paste0(names(adstocks.min  ), ".adstock")
        names(adstocks.max  ) <- paste0(names(adstocks.max  ), ".adstock")
    }
    
    transform.start <- transform[(names(transform) %in% c(fixed.vars, random.vars)) & !(transform %in% "log")]
    transform.min <- rep(-Inf, length(transform.start))
    transform.max <- rep( Inf, length(transform.start))
    if(length(transform.start)) {
        transform.start <- sapply(transform.start, 
                                  function(x) {
                                      if(tolower(x)=="power") {
                                          1
                                      } else if(tolower(x)=="powerinv") {
                                          0
                                      } else if(tolower(x)=="negexp") {
                                          1
                                      }
                                  })

        transform.min <- rep(-Inf, length(transform.start))
        transform.max <- rep( Inf, length(transform.start))
        
        names(transform.start) <- paste0(names(transform.start), ".power")
        names(transform.min  ) <- paste0(names(transform.start), ".power")
        names(transform.max  ) <- paste0(names(transform.start), ".power")
    } else {
        transform.start <- numeric()
    }
    
    nls.start <- c(coef.start, adstocks.start, transform.start)
    nls.min   <- c(coef.min  , adstocks.min  , transform.min)
    nls.max   <- c(coef.max  , adstocks.max  , transform.max)
      
    ## QA
    #cat("\nnls.start: \n"); print(nls.start)
    #cat("\nnls.min:   \n"); print(nls.min  )
    #cat("\nnls.max:   \n"); print(nls.max  )
    #cat("\n")

    # Build formula
    formula2 <- MMR.Optim.Formula(dep.var, fixed.vars, min.adstock[fixed.vars], random.vars, transform)

    ## QA
    #cat("formula2: \n"); cat(formula2)
    
    data$inTrain <- inTrain
    
    ## QA
    #cat("\n-----------------------------------------------------------------------")
    #cat("\n-----------------------------------------------------------------------\n\n")

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
#' @param fixed A character vector of variable names that won't be adstocked
#' @param fixed.adstocked A numeric vector same length as \code{fixed}
#' @param random A character vector of variable names to be adstocked
#' 
#' @author Gabriel Mohanna
#' 
MMR.Optim.Formula <- function(dep.var, 
                              fixed=character(), 
                              fixed.adstock=rep(0, length(fixed)),
                              random=character(),
                              transform=character()) {
    # Apply adstock transformations
    fixed <- mapply(
        function(var.name, var.value) {
            if(var.value) {
                paste0("adstock(", var.name, ", ", var.value, ")")
            } else {
                paste0(var.name)
            }
        },
        fixed, fixed.adstock,
        SIMPLIFY = F)
    
    random <- mapply(
        function(var.name) {
            paste0("adstock(", var.name, ", ", var.name, ".adstock)")
        },
        random,
        SIMPLIFY = F)
    
    adstock.transfrom <- append(fixed, random)
    
    # Apply other functional transformations were they apply
    other.transform <- mapply(
        function(var.name, var.adstock) {
            if(var.name %in% names(transform)) {
                var.transform <- transform[var.name]
                if(var.transform=="power") {
                    paste0(var.name, ".coef*(", var.adstock, "^", var.name, ".power)")
                } else if(var.transform=="negexp") {
                    paste0(var.name, ".coef*", "(1 - exp(-", var.name, ".power * ", var.adstock, "))")
                } else if(var.transform=="log") {
                    paste0(var.name, ".coef*", "log(1+", var.adstock, ")")
                } else if(var.transform=="powerinv") {
                    paste0(var.name, ".coef/", "(", var.name, ".power + 1/", var.adstock, ")")
                }
            } else {
                paste0(var.name, ".coef*", var.adstock)
            }
        },
        names(adstock.transfrom), adstock.transfrom,
        SIMPLIFY = F)
    
    paste0(dep.var, " ~ ",
           "inTrain * (Intercept", 
           paste(rep(" +", length(other.transform)), other.transform, collapse = ""),
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
