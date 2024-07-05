library(RWiener)
library(dplyr)
library(loo)
library(coda)
library(rjags)
library(runjags)
load.module("wiener")
load.module("lecuyer")
set.seed(2024)

setwd("C:/Users/u0166113/OneDrive - KU Leuven/Documents/Affectometrics/Study 1")
data <- read.table("choice_18-06-2024_cleaned_rt3000ms100ms.csv", header = TRUE, sep = ",")  # Adjust 'sep' based on your file's delimiter
#true.parameters <- read.table("M4 true parameters.txt", header = TRUE, sep = ",")  # Adjust 'sep' based on your file's delimiter

#load("jags_sample_md3.RData")
#traceplot(jagsModel_md3$mcmc[,"MuV[1]"])
#traceplot(jagsModel_md3$mcmc[,"MuV[2]"])
#traceplot(jagsModel_md3$mcmc[,"MuV[3]"])
#traceplot(jagsModel_md3$mcmc[,"MuV[4]"])
#traceplot(jagsModel_md3$mcmc[,"MuV[5]"])
#traceplot(jagsModel_md3$mcmc[,"MuV[6]"])
#traceplot(jagsModel_md3$mcmc[,"MuV[7]"])
#traceplot(jagsModel_md3$mcmc[,"MuV[8]"])
#traceplot(jagsModel_md3$mcmc[,"MuV[9]"])
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



data_model_md3=list(
  #allow for different numbers of trials
  nTrials = length(data$rt),
  nSubject = length(unique(data$Sub)),
  subject = data$Sub,
  o.dummy = ifelse(data$outcome>0, 1, 0),
  o.abs.win = data$outcome.absolute.win,
  o.abs.loss = data$outcome.absolute.loss,
  pe.win = data$prediction.error.win,
  pe.loss = data$prediction.error.loss,
  o.pe.win = (data$outcome.absolute.win-mean(data$outcome.absolute.win)) *(data$prediction.error.win-mean(data$prediction.error.win)),
  o.pe.loss = (data$outcome.absolute.loss-mean(data$outcome.absolute.loss)) *(data$prediction.error.loss-mean(data$prediction.error.loss)),
  prob.win = data$probability.win,
  prob.loss = data$probability.loss,
  o.win = as.numeric(data$outcome.report.win),
  o.loss = as.numeric(data$outcome.report.loss),
  y = ifelse(data$Decision==1, data$rt, -data$rt),
  Dec_previous = data$Decision.lag.1,
  nVarZ = 3,
  nVarV = 9,
  RscalZ = 5 ,  
  RscalV = 11 ,  
  RmatZ = diag(x=1,nrow = 3),
  RmatV = diag(x=1,nrow = 9))


model_String_md3 ='
model {
  
  #likelihood function
  for (t in 1:nTrials){
                 
      #model the response bias
      p[t] = Z[subject[t],1] + Z[subject[t],2] * (prob.win[t]*o.win[t] +
                     prob.loss[t]*o.loss[t]) + Z[subject[t],3]* Dec_previous[t]
      beta[t] = 0.8*(1/(1+exp(-p[t])))+0.1
      
      #model the drift rate
      delta[t] = V[subject[t],1] + V[subject[t],2] * o.dummy[t] + V[subject[t],3] * o.abs.win[t] + V[subject[t],6] * o.abs.loss[t] + 
                    V[subject[t],4] * pe.win[t] + V[subject[t],7] * pe.loss[t] +
                    V[subject[t],5] * o.pe.win[t] + V[subject[t],8] * o.pe.loss[t] +
                    V[subject[t],9] * Dec_previous[t]
  
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

writeLines(model_String_md3, con = 'md3.txt')



#init function
initfunction_md3 <- function(chain){
  return(list(
    mualpha = 1,
    mutau  = 0.01,
    prectau = 10,
    .RNG.name = "lecuyer::RngStream"))
}

#Create list of parameters to be monitored
parameters_md3 <- c("alpha", "tau", "Z","V",
                    "mualpha", "mutau", "MuZ", "MuV",
                    "precalpha", "prectau", "SigmaZ","SigmaV", "RhoZ", "RhoV","ld")




#Run the model in runjags
startTime = proc.time()
jagsModel_md3 <- run.jags(method = "parallel",
                          model = 'md3.txt',
                          monitor = parameters_md3,
                          data = data_model_md3,
                          n.chains = 3,
                          adapt = 15000, #how long the samplers "tune"
                          burnin = 3000, #how long of a burn in
                          sample = 3000,
                          inits = initfunction_md3,
                          thin = 3, #thin if high autocorrelation to avoid huge files
                          modules = c("wiener", "lecuyer"),
                          summarise = F,
                          plots = T)
stopTime = proc.time()
elapsedTime = stopTime - startTime
print(elapsedTime/60) #Tells how long it took to run analysis

save(jagsModel_md3, file = "jags_sample_md3_Tomi.RData")

traceplot(jagsModel_md3$mcmc[,"MuV[1]"])
traceplot(jagsModel_md3$mcmc[,"MuV[2]"])
traceplot(jagsModel_md3$mcmc[,"MuV[3]"])
traceplot(jagsModel_md3$mcmc[,"MuV[4]"])
traceplot(jagsModel_md3$mcmc[,"MuV[5]"])
traceplot(jagsModel_md3$mcmc[,"MuV[6]"])
traceplot(jagsModel_md3$mcmc[,"MuV[7]"])
traceplot(jagsModel_md3$mcmc[,"MuV[8]"])
traceplot(jagsModel_md3$mcmc[,"MuV[9]"])

mcmc_md3 = as.matrix(as.mcmc.list(jagsModel_md3$mcmc), chains = F)

#group level parameters plot
par(mfrow = c(3, 3), cex.main = 3, cex.axis=1.5)
options(repr.plot.width = 20, repr.plot.height = 8)  
hist(mcmc_md3[, "MuZ[1]"], breaks = 100, main = expression(mu[bold(z[0])]), xlab = expression(mu[bold(z[0])]))
hist(mcmc_md3[, "MuZ[2]"], breaks = 100, main = expression(mu[bold(z[1])]), xlab = expression(mu[bold(z[1])]))
hist(mcmc_md3[, "MuZ[3]"], breaks = 100, main = expression(mu[bold(z[phi])]), xlab = expression(mu[bold(z[phi])]))
hist(mcmc_md3[, "MuV[1]"], breaks = 100, main = expression(mu[bold(v[L0])]), xlab = expression(mu[bold(v[L0])]))
hist(mcmc_md3[, "MuV[2]"], breaks = 100, main = expression(mu[bold(v[w0])]), xlab = expression(mu[bold(v[w0])]))
hist(mcmc_md3[, "MuV[3]"], breaks = 100, main = expression(mu[bold(v[w1])]), xlab = expression(mu[bold(v[w1])]))
hist(mcmc_md3[, "MuV[4]"], breaks = 100, main = expression(mu[bold(v[w2])]), xlab = expression(mu[bold(v[w2])]))
hist(mcmc_md3[, "MuV[5]"], breaks = 100, main = expression(mu[bold(v[w3])]), xlab = expression(mu[bold(v[w3])]))
hist(mcmc_md3[, "MuV[6]"], breaks = 100, main = expression(mu[bold(v[L1])]), xlab = expression(mu[bold(v[L1])]))
hist(mcmc_md3[, "MuV[7]"], breaks = 100, main = expression(mu[bold(v[L2])]), xlab = expression(mu[bold(v[L2])]))
hist(mcmc_md3[, "MuV[8]"], breaks = 100, main = expression(mu[bold(v[L3])]), xlab = expression(mu[bold(v[L3])]))
hist(mcmc_md3[, "MuV[9]"], breaks = 100, main = expression(mu[bold(v[phi])]), xlab = expression(mu[bold(v[phi])]))
hist(mcmc_md3[, "mualpha"], breaks = 100, main = expression(mu[bold(alpha)]), xlab = expression(mu[bold(alpha)]))
hist(mcmc_md3[, "mutau"], breaks = 100, main = expression(mu[bold(tau)]), xlab = expression(mu[bold(tau)]))
