
library(progress)
library(MASS)

# BF calculation function
Bf<-function(sd, obtained, dfdata, meanoftheory, sdtheory, dftheory, tail = 2)
{
  area <- 0
  normarea <- 0
  theta <- meanoftheory - 10 * sdtheory
  incr <- sdtheory / 200
  for (A in -2000:2000){
    theta <- theta + incr
    dist_theta <- dnorm(theta, meanoftheory, sdtheory)
    dist_theta <- dt((theta-meanoftheory)/sdtheory, df=dftheory)
    if(identical(tail, 1)){
      if (theta <= 0){
        dist_theta <- 0
      } else {
        dist_theta <- dist_theta * 2
      }
    }
    height <- dist_theta * dt((obtained-theta)/sd, df = dfdata)
    area <- area + height * incr
    normarea <- normarea + dist_theta*incr
  }
  LikelihoodTheory <- area/normarea
  Likelihoodnull <- dt(obtained/sd, df = dfdata)
  BayesFactor <- LikelihoodTheory / Likelihoodnull
  BayesFactor
}





# function simulating studies and sequential analyses
simulate_PLC_study <- function(N,
                               num_conditions,
                               mean_exp,
                               SD_exp,
                               mean_con,
                               SD_con,
                               correlation){

pb$tick() # progress bar tick (useful to track progress of long analyses)

# generate two correlated random variables with mean = 0 and sd = 1
data = as.data.frame(  mvrnorm(N, mu = c(0, 0),
                               Sigma = matrix(c(1, correlation,
                                                correlation, 1),
                                              ncol = 2, byrow = TRUE),
                               empirical = F))
names(data) = c("outcome_exp", "outcome_con")

# transform these variables to match mean and sd from pilot study
data[,"outcome_exp"] = data[,"outcome_exp"] * SD_exp + mean_exp
data[,"outcome_con"] = data[,"outcome_con"] * SD_con + mean_con



# generate simulated data based on pilot study derived mean and sd data
# data[,"outcome_exp"] = rnorm(N, mean = mean_exp, sd = SD_exp)
# data[,"outcome_con"] = rnorm(N, mean = mean_con, sd = SD_con)

# extract Bayes Factor
# the analysis uses sequential analysis and optional stopping if BF < 0.333 or BF > 3.
# the Bf code was imported from the pilot_analysis.Rmd file
# doing analysis after every participant after reaching N = 20
#for(i in 20:N){
 for(i in seq(20,N,10)){ # use this line if you have to test large sample sizes, to reduce procesing time. This simulates analysis after every 10 participants
  data_current = data[1:i,]
  effect <- t.test(data[,"outcome_exp"], data[,"outcome_con"], paired=TRUE)
  BF = Bf(sd = as.numeric(effect$estimate[1]/effect$statistic[1]), obtained = as.numeric(effect$estimate[1]), dfdata = effect$parameter,  meanoftheory = 0, sdtheory = 50, dftheory = 10000, tail = 1)
  if(BF > 3){break}
  if(BF < 0.3333){break}
}

mean_dif = mean(data_current[,"outcome_exp"] - data_current[,"outcome_con"])
sd_dif = sd(data_current[,"outcome_exp"] - data_current[,"outcome_con"])
cor = cor(data[,"outcome_exp"], data[,"outcome_con"])

# the function returns the Bayes Factor value at the stopping point of the study
return(c(BF, mean_dif, sd_dif, cor))
}







iterations = 10000

### original pilot study data
# mean_exp = 50.61047 # mean Stroop interference in the volitional condition
# SD_exp = 54.79365 # sd of Stroop interference in the volitional condition
# mean_con = 25.66169 # mean Stroop interference in the suggestion condition
# SD_con = 49.86554 # sd of Stroop interference in the suggestion condition


pb <- progress_bar$new(
  format = " simulation progress [:bar] :percent eta: :eta",
  total = iterations, clear = FALSE, width= 60)

out = replicate(iterations, 
                simulate_PLC_study(N = 80,
                                   num_conditions = 2,
                                   mean_exp = 50.61047,
                                   SD_exp = 54.79365,
                                   mean_con = 25.66169,
                                   SD_con = 49.86554,
                                   correlation = 0.1729514))

out_table = as.data.frame(t(out))
names(out_table) = c("BF", "mean_dif", "sd_dif", "cor")

# just to check if the simulation produced the desired mean difference and correlation
# the original study had the following parameters:
# mean_dif = 24.94878 # mean difference of Stroop interference between the suggestion and violation condition
# sd_dif = 67.40772 # sd of difference of Stroop interference between the suggestion and violation condition
# cor = 0.1729514 # correlation of the Stroop interference between the suggestion and violation condition

mean(out_table[,"mean_dif"])
mean(out_table[,"sd_dif"])
mean(out_table[,"cor"])

mean(out_table[,"BF"] > 3) # 0.87
mean(out_table[,"BF"] < 0.3333) # 0.0069


### original study data
# assuming SD = 50 in stroop interference in both groups, and means equal in the two conditions


pb <- progress_bar$new(
  format = " simulation progress [:bar] :percent eta: :eta",
  total = iterations, clear = FALSE, width= 60)

out = replicate(iterations, 
                simulate_PLC_study(N = 80,
                                   num_conditions = 2,
                                   mean_exp = 25,
                                   SD_exp = 50,
                                   mean_con = 25,
                                   SD_con = 50,
                                   correlation = 0.1729514))


out_table = as.data.frame(t(out))
names(out_table) = c("BF", "mean_dif", "sd_dif", "cor")

mean(out_table[,"BF"] > 3) # 0.015
mean(out_table[,"BF"] < 0.3333) # 0.8146
