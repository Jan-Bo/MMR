#' Training observations selection.
#' 
#' Randomly assigns each observation into training set without replacement
#' 
#' @param data A data frame
#' @param p The percentage of data that goes to training. This can be a constact or a function.
#' @param times The number of partitions to create
#' @param type A character indicating if indicator of inTrain is boolean. Anything else is will produce 1 or 0.
#' 
#' @author Gabriel Mohanna
#' 
#' @details
#' Creates a series of training observations that are randomly selected from the 'data' with probability 'p' without replacement.
#' 
#' @examples
#' createInTrainSamples(data=iris)
#' createInTrainSamples(data=iris, p=.7)
#' createInTrainSamples(data=iris, times=20)
#' createInTrainSamples(data=iris, times=20, p=.7)
#' 
#' @export
#' 
createInTrainSamples <- function(data, times=1, p=0.5, type="bool"){
    if(missing(data))
        stop("Argument \"data\" is missing.")
    
    inTrainObs <- TRUE; outTrainObs <- FALSE
    if(type!="bool") {
        inTrainObs <- 1; outTrainObs <- 0
    } 
    
    samples <- rep(list(rep(outTrainObs, nrow(data))), times)
    samples <- lapply(samples, 
                      function(trial){
                          inTrainProb = p
                          if(class(p)=="function")
                              inTrainProb = p()
                          
                          sapply(trial, 
                                 function(obs){
                                     ifelse(runif(1)<=inTrainProb, inTrainObs, outTrainObs)
                                 })
                      })
    names(samples) <- paste0("sample", 1:times)
    
    samples
}

#' Training observations selection.
#' 
#' Sequentially divides the set of observations into k groups, or folds
#' 
#' @param data A data frame
#' @param k An integer for the number of folds
#' @param type A character indicating if indicator of inTrain is boolean. Anything else is will produce 1 or 0.
#' 
#' @author Gabriel Mohanna
#' 
#' @examples
#' createInTrainFolds(data=iris)
#' createInTrainFolds(data=iris, 5)
#' 
#' @details
#' Observations are divided sequentially into set of k groups, or folds, of approximately equal size. The first fold is treated as a validation set and the other k-1 folds are used for training. This process is repeated for each of the k folds.
#' 
#' @export
#' 
createInTrainFolds <- function(data, k=10, type="bool"){
    if(missing(data))
        stop("Argument \"data\" is missing.")
    
    inTrainObs <- TRUE; outTrainObs <- FALSE
    if(type!="bool") {
        inTrainObs <- 1; outTrainObs <- 0
    } 
    
    k.size <- ceiling(nrow(data)/k)
    folds <- c(1:k)
    folds <- lapply(folds, 
                     function(fold){
                         obss <- c(rep(inTrainObs, nrow(data)))
                         obss[min((k.size*(fold-1)+1), nrow(data)):min(k.size*fold, nrow(data))] <- outTrainObs
                         obss
                     })
    names(folds) <- paste0("fold", 1:k)
    
    folds
}


#' Trials Generations
#' 
#' Generates a list of inTrain tials
#' 
#' @param data A data frame
#' @param cvControl A list of cross-validation selection parameters.
#' 
#' @author Gabriel Mohanna
#' 
MMR.Get.Trials <- function(data, cvControl=cv.control("none")){
    if(missing(data))
        stop("Argument \"data\" is missing.")
    
    if(cvControl$method=="random") {
        trials <- createInTrainSamples(data, times=cvControl$times, p=cvControl$p)
    } else if(cvControl$method=="k-fold") {
        trials <- createInTrainFolds(data, k=cvControl$k)
    } else if(cvControl$method=="leave-one-out") {
        trials <- createInTrainFolds(data, k=nrow(data))
    } else {
        trials <- lapply(createInTrainFolds(data, k=1), function(x) !x)
    } 
    
    trials
} 


#' Observations Splitting Functions
#' 
#' Not yet implemented
#' 
#' @param data A data frame
#' @param initialWindow Something
#' @param horizon Something
#' @param fixedWindow Something
#' 
#' @author Gabriel Mohanna
#' 
#' @details
#' <empty> for now
createInTrainTimeSlices <- function(data, initialWindow, horizon=1, fixedWindow=T) {
    stop("createInTrainTimeSlices is not yet implemented!")
}
