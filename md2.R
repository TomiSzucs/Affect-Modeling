library(RWiener)
library(dplyr)
library(loo)
library(coda)
library(rjags)
library(runjags)
load.module("wiener")
load.module("lecuyer")
#set.seed(1998)

setwd("C:/Users/u0166113/OneDrive - KU Leuven/Documents/Affectometrics/Study 1")
data <- read.table("choice_18-06-2024_cleaned_rt3000ms100ms.csv", header = TRUE, sep = ",")  # Adjust 'sep' based on your file's delimiter
#true.parameters <- read.table("M4 true parameters.txt", header = TRUE, sep = ",")  # Adjust 'sep' based on your file's delimiter

##############################################preprocessing###################################################
##############################################################################################################
#exclude participants who produce consecutive responses
pair = paste0(data$participant,data$affect)
consecutive_response = rle(pair)
summary(consecutive_response$lengths)
consecutive_pair <-consecutive_response$values[consecutive_response$lengths >= 15]
consecutive_pair = sapply(consecutive_pair, function(x) substring(x, 1, nchar(x) - 1))
consecutive_pair = as.numeric(unique(consecutive_pair))

data =  data[ !(data$participant %in% consecutive_pair),]


#assign an index from 1 to n to participants
data$Sub = as.numeric(factor(data$participant, levels = unique(data$participant)))

# Separate the column 'outcome_report' into two columns 'outcome.report.win' and 'outcome.report.loss'
data <- cbind(data, do.call(rbind, strsplit(gsub("\\[|\\]", "", data$outcome_report), ",")))
# Rename the columns
names(data)[names(data) == "1"] <- "outcome.report.win"
names(data)[names(data) == "2"] <- "outcome.report.loss"         

#Separate "probability" into "probability.win" and "probability.loss"
data$probability.win = ifelse(data$probability %in% c(1,2,3), 0 ,1)
data$probability.loss = ifelse(data$probability %in% c(1,2,3), 1,0)
#data$probability.win = data$probability
#data$probability.loss = 16 - data$probability


#create prediction error for winning & losing condition
data$prediction.error.win = ifelse(data$outcome >0, data$probability.loss ,0)
data$prediction.error.loss = ifelse(data$outcome < 0, data$probability.win ,0)
#create prediction error for md5
data$prediction.error = ifelse(data$outcome>0, data$probability.loss, data$probability.win)
#convert true outcome to its absolute value
data$outcome.absolute.win = ifelse(data$outcome>0, abs(data$outcome), 0)
data$outcome.absolute.loss = ifelse(data$outcome<0, abs(data$outcome), 0)

#convert reaction time from ms to s 
data$rt = data$rt/1000
names(data)[names(data) == "affect"] <- "Decision"
data$Decision.lag.1 <- c(NA, data$Decision[-nrow(data)])
data$Decision.lag.1[which(!duplicated(data$Sub))] <- 0

data_model_md2=list(
  #allow for different numbers of trials
  nTrials = length(data$rt),
  nSubject = length(unique(data$Sub)),
  subject = data$Sub,
  o.true = data$outcome,
  o.abs = abs(data$outcome),
  o.win = as.numeric(data$outcome.report.win),
  o.loss = as.numeric(data$outcome.report.loss),
  o.abs.win = abs(as.numeric(data$outcome.report.win)),
  o.abs.loss = abs(as.numeric(data$outcome.report.loss)),
  p.win = data$probability.win,
  p.loss = data$probability.loss,
  pe = data$prediction.error,
  y = ifelse(data$Decision>0, data$rt, -data$rt),
  Dec_previous = data$Decision.lag.1,
  o.dummy = ifelse(data$outcome>0,1,0),
  nVarZ = 3,
  nVarV = 5,
  RscalZ = 5 ,  
  RscalV = 7 ,  
  RmatZ = diag(x=1,nrow = 3),
  RmatV = diag(x=1,nrow = 5))



model_String_md2 ='
model {
  
  #likelihood function
  for (t in 1:nTrials){
                 
      #model the response bias
      p[t] = Z[subject[t],1] + Z[subject[t],2] * (p.win[t]*o.win[t] +
                     p.loss[t]*o.loss[t]) + Z[subject[t],3]* Dec_previous[t]
      beta[t] = 0.8*(1/(1+exp(-p[t])))+0.1
      
      #model the drift rate
      delta[t] = V[subject[t],1] + V[subject[t],2] * o.dummy[t] + V[subject[t],3] * o.abs.win[t] + V[subject[t],4] * o.abs.loss[t] + 
           V[subject[t],5] * Dec_previous[t]
  
  
      y[t] ~ dwiener(alpha[subject[t]], 
                 tau[subject[t]], 
                 beta[t], 
                 delta[t])
  
  }
  
  #looic
  for (t in 1:nTrials){
    
     ld[t] = logdensity.wiener(y[t], alpha[subject[t]], 
                 tau[subject[t]], 
                 beta[t], 
                 delta[t])
    
  }
  
    
  #priors
  for (i in 1:nSubject) {
    
    Z[i,1:nVarZ] ~ dmnorm(MuZ[1:nVarZ], InvCovMatZ[1:nVarZ,1:nVarZ])
    V[i,1:nVarV] ~ dmnorm(MuV[1:nVarV], InvCovMatV[1:nVarV,1:nVarV])
    
    alpha[i] ~ dnorm(mualpha,precalpha)
    
    #beta distribution
    tau[i] ~ dnorm(mutau,prectau)T(0.001, 3)

  }

  #hyperpriors
  
  for (varIdx in 1:nVarZ) { 
  MuZ[varIdx] ~ dnorm( 0 , 1e-3 ) 
  }
  
  for (varIdx in 1:nVarV) { 
  MuV[varIdx] ~ dnorm( 0 , 1e-3 ) 
  }
  
  InvCovMatZ ~ dwish(RmatZ[1:nVarZ,1:nVarZ], RscalZ )
  InvCovMatV ~ dwish(RmatV[1:nVarV,1:nVarV], RscalV )
  CovMatZ <- inverse( InvCovMatZ )
  CovMatV <- inverse( InvCovMatV )
  
  for ( varIdx in 1:nVarZ ) { 
    SigmaZ[varIdx] <- sqrt(CovMatZ[varIdx,varIdx]) 
  }
  for ( varIdx in 1:nVarV ) { 
    SigmaV[varIdx] <- sqrt(CovMatV[varIdx,varIdx]) 
  }
    
  for ( varIdx1 in 1:nVarZ ) { 
    for ( varIdx2 in 1:nVarZ ) {
      RhoZ[varIdx1,varIdx2] <- ( CovMatZ[varIdx1,varIdx2]
                               / (SigmaZ[varIdx1]*SigmaZ[varIdx2]) )
    } }
  
  for ( varIdx1 in 1:nVarV ) { 
    for ( varIdx2 in 1:nVarV ) {
      RhoV[varIdx1,varIdx2] <- ( CovMatV[varIdx1,varIdx2]
                               / (SigmaV[varIdx1]*SigmaV[varIdx2]) )
  } }
   
  mualpha ~ dunif(.1, 10) 
  mutau ~ dunif(.0001, 3)
  
  precalpha  ~ dgamma(.001, .001)
  prectau ~ dgamma(.001, .001)
}


'

writeLines(model_String_md2, con = 'Model D2.txt')

#init function
initfunction_md2 <- function(chain){
  return(list(
    mualpha = 1,
    mutau  = .001,
    precalpha = runif(1, 1, 100),
    prectau   = runif(1, 1, 100),
    .RNG.name = "lecuyer::RngStream"))
}




#Create list of parameters to be monitored
parameters_md2 <- c("alpha", "tau", "Z","V",
                    "mualpha", "mutau", "MuZ", "MuV",
                    "precalpha", "prectau", "SigmaZ","SigmaV", "RhoZ", "RhoV", "ld")





#Run the model in runjags
startTime = proc.time()
jagsModel_md2 <- run.jags(method = "parallel",
                          model = 'Model D2.txt',
                          monitor = parameters_md2,
                          data = data_model_md2,
                          inits = initfunction_md2,
                          n.chains = 3,
                          adapt = 3000, #how long the samplers "tune"
                          burnin = 3000, #how long of a burn in
                          sample = 3000,
                          thin = 3, #thin if high autocorrelation to avoid huge files
                          modules = c("wiener", "lecuyer"),
                          summarise = F,
                          plots = T)
stopTime = proc.time()
elapsedTime = stopTime - startTime
show(elapsedTime/60) #Tells how long it took to run analysis

save(jagsModel_md2, file = "jags_sample_MD2.RData")