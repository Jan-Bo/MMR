#' Advertising Data Simulation
#' 
#' Creates a random advertising data
#' 
#' @param weeks Number of weeks, or observations
#' @param trend An indicator of a trend sequance
#' @param season An indicator of seasonal effect
#' @param ads.count Count of advertising medium or campaigns
#' @param ads.cor Additial ads that are correlated with 'ads.count'
#' @param ads.noise A count of non-related ads with no effects
#' @param adstock.max The maximum adstock rate to give to each 'ads.count'
#' 
#' @author Gabriel Mohanna
#' 
#' @details
#' Simulates a randomly generated advertising data
#' 
#' @examples
#' # Simulate 2 years of data with trend, seasonality and 3 different ads.
#' sim.ads(weeks=104, trend=TRUE, season=TRUE, ads.count=3)
#' 
#' # Simulate 2 years of data with 3 different ads and 1 unrelated variable.
#' sim.ads(weeks=104, trend=FALSE, season=FALSE, ads.count=3, ads.noise=1)
#' 
#' @export
sim.ads <- function(weeks, trend=T, season=T, ads.count=0, sequantial=T, adstock.max=0, ads.cor=0, ads.noise=0){
    # Create ad.data dataframe and true fit model
    ad.data        <- data.frame(week=1:weeks)
    modelsFit.coef <- data.frame(Param=numeric(), True=numeric())
    data.summary   <- data.frame(Param=numeric(), Start=numeric(), End=numeric(), Coef=numeric(), Adstock=numeric(), Contribution=numeric(), Beta=numeric(), Form=numeric())
    
    beta.ctr <- 0
    
    # Add intercept
    Intercept.coef   <- rep(runif(1, 0, 200), 1)
    Intercept.impact <- Intercept.coef * rep(1, weeks)
    
    modelsFit.coef[beta.ctr+1, ] <- list("Intercept", Intercept.coef)
    data.summary[beta.ctr+1  , ] <- list("Intercept", 1, weeks, Intercept.coef, 0, sum(Intercept.impact), paste0("b", beta.ctr), paste0("b", beta.ctr))
    beta.ctr = beta.ctr + 1
    
    # Add trend
    if(trend){
        ad.data$Trend <- c(1:weeks)
        Trend.coef    <- runif(1, 10, 50)
        Trend.impact  <- Trend.coef * ad.data$Trend
        
        modelsFit.coef[beta.ctr+1, ] <- list("Trend", Trend.coef)
        data.summary[beta.ctr+1  , ] <- list("Trend", 1, weeks, Trend.coef, 0, sum(Trend.impact), paste0("b", beta.ctr), paste0("b", beta.ctr, "*", "Trend"))
        beta.ctr = beta.ctr + 1
    } else {
        Trend.impact <- rep(0, weeks)
    }
    
    # Add seasonality
    if(season){
        # Introductory Time Series with R Page 102
        ad.data$Season <- sin(2*pi*(1:weeks)/52)
        ad.data$Season <- ad.data$Season+abs(min(ad.data$Season))
        Season.coef    <- runif(1, 800, 1200)
        Season.impact  <- Season.coef * ad.data$Season
        
        modelsFit.coef[beta.ctr+1, ] <- list("Season", Season.coef)
        data.summary[beta.ctr+1  , ] <- list("Season", 1, weeks, Season.coef, 0, sum(Season.impact), paste0("b", beta.ctr), paste0("b", beta.ctr, "*", "Season"))
        beta.ctr = beta.ctr + 1
    } else {
        Season.impact <- rep(0, weeks)
    }
    
    # Create empty vector for ads.impact
    ads.impact <- rep(0, weeks)
    
    # Create ads with different campaigns
    if(ads.count > 0) {
        for(i in 1:ads.count) {
            # Create ad
            ad            <- rep(0, weeks)
            if(sequantial) {
                start <- round(runif(1, weeks/ads.count*(i-1) + 1, weeks/ads.count*i))
            } else {
                start <- round(runif(1, 1, weeks))
            }
            end       <- round(runif(1, start, weeks))
            
            ## QA
            #cat(c("start: ", start, "\n"))
            #cat(c("end: ", end, "\n"))
            #cat(c("series: ", ad, "\n"))
            #cat(c("ads: ", ad[start:end], "\n"))
            #cat("--------------------------------------\n")
            
            ad[start:end] <- pmax(runif(end - start + 1, runif(1, -100, -1), runif(1, 1, 1000)), 0)
            
            # Create parameters & impact
            coef    <- runif(1, 2, 15)
            adstock <- runif(1, 0, adstock.max) * rbinom(1, 1, .7)
            #power   <- runif(1, 0, 1          ) ^ rbinom(1, 1, .25)
            power   <- 1
            impact  <- coef * MMR.Transform(ad, adstock, power)
            
            # Combine to full ads impact
            ads.impact <- ads.impact + impact
            
            # Add to ad.data, modelsFit.coef & data.summary
            ad.data[, paste0("ad", i)] <- ad
            
            data.summary[beta.ctr+1, ] <- 
                list(paste0("ad", i), 
                     start, 
                     end, 
                     coef,
                     adstock,
                     sum(impact),
                     paste0("b", beta.ctr),
                     paste0("b", beta.ctr, "*", "adstock.transform(ad", i, ", a", i,")"))
            
            beta.ctr = beta.ctr + 1
            
            #rm(list=c(paste0("ad", i, ".start"), paste0("ad", i, ".end")))
        }
    }
    
    # Create correlated ads
    if(ads.cor > 0) {
        for(j in 1:ads.cor) {
            
        }        
    }
    
    # Create noise ads that is unrelated to the sales data
    if(ads.noise > 0) {
        for(i in 1:ads.noise) {
            # Create ad
            ad            <- rep(0, weeks)
            start         <- round(runif(1, 1, weeks))
            end           <- round(runif(1, start, weeks))
            ad[start:end] <- pmax(runif(end - start + 1, runif(1, -100, -1), runif(1, 1, 1000)), 0)
            
            # Add to ad.data, modelsFit.coef & data.summary
            ad.data[, paste0("noise", i)] <- ad
            
            data.summary[beta.ctr+1, ] <- 
                list(paste0("noise", i), 
                     start, 
                     end, 
                     0,
                     0,
                     0,
                     paste0("b", beta.ctr),
                     paste0("b", beta.ctr, "*", "adstock.transform(noise", i, ", a", ads.count+ads.cor+i,")"))
            
            beta.ctr = beta.ctr + 1
        }        
    }
    
    # Create error
    error <- rnorm(weeks, 0, 10)
    
    # Create final sales & dataset
    ad.data$sales <- Intercept.impact + Trend.impact + Season.impact + ads.impact + error
    
    # Calculate % contributions
    data.summary$Contrib.Pct <- data.summary$Contribution / sum(data.summary$Contribution)
    data.summary <- data.summary[,c(1:6,9,7:8)]
    
    # Build final return list
    data.summary <<- data.summary
    
    par(mfrow = c(1,1), mar=c(3,3,4,1))
    barplot(t(ad.data[,paste0("ad", c(1:ads.count))]), 
            col=1:ads.count, border=F, main="Sales vs. Advertising", axes=F)
    axis(1, tck=F, pos=0, labels=NA)
    axis(2, tck=F, pos=0, labels=NA, las=1)
    par(new=T)
    plot(x=ad.data$week, y=ad.data$sales, type="l", 
         main="Sales vs. Advertising", xlab="Time (Usually Weeks)", ylab="Sales", frame=F, axes=F)
    
    print(data.summary[, 1:7])
    
    return(ad.data)
}
