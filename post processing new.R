library(RWiener)
library(dplyr)
library(loo)
library(coda)
library(rjags)
library(ggplot2)
library(MCMCvis)
library(stringr)
library(tidyr)
#set.seed(1998)

setwd("C:/Users/u0166113/OneDrive - KU Leuven/Documents/Affectometrics/Study 1")
setwd("C:/Users/u0166113/OneDrive - KU Leuven/Documents/Affectometrics/Study 2")

participant_estimate_extractor <- function(mcmc_obj){
  # Filter columns with names matching the format "beta[i,j]"
  filtered_columns <- grep("^beta\\[\\d+,\\d+\\]$", colnames(mcmc_obj), value = TRUE)
  parsed_names <- str_match(filtered_columns, "beta\\[(\\d+),(\\d+)\\]")
  participants <- as.integer(parsed_names[,2])
  parameters <- as.integer(parsed_names[,3])
  
  # Subset the matrix to keep only these columns
  filtered_mcmc <- mcmc_obj[, filtered_columns]
  posterior_mean <- colMeans(filtered_mcmc)
  
  mean_df <- data.frame(
    participant = participants,
    parameter = parameters,
    mean = posterior_mean
  )
  
  mean_df_wide <- pivot_wider(mean_df, names_from = parameter, values_from = mean, names_prefix = "beta_")
  return(mean_df_wide)
}

participant_estimate_extractor_drift <- function(mcmc_obj){
  # Filter columns with names matching the format "beta[i,j]"
  filtered_columns <- grep("^V\\[\\d+,\\d+\\]$", colnames(mcmc_obj), value = TRUE)
  parsed_names <- str_match(filtered_columns, "V\\[(\\d+),(\\d+)\\]")
  participants <- as.integer(parsed_names[,2])
  parameters <- as.integer(parsed_names[,3])
  
  # Subset the matrix to keep only these columns
  filtered_mcmc <- mcmc_obj[, filtered_columns]
  posterior_mean <- colMeans(filtered_mcmc)
  
  mean_df <- data.frame(
    participant = participants,
    parameter = parameters,
    mean = posterior_mean
  )
  
  mean_df_wide <- pivot_wider(mean_df, names_from = parameter, values_from = mean, names_prefix = "beta_")
  return(mean_df_wide)
}

participant_estimate_extractor_PT <- function(mcmc_obj){
  # Filter columns with names matching the format "beta[i,j]"
  filtered_columns <- grep("^beta\\[\\d+,\\d+\\]$", colnames(mcmc_obj), value = TRUE)
  parsed_names <- str_match(filtered_columns, "beta\\[(\\d+),(\\d+)\\]")
  participants <- as.integer(parsed_names[,2])
  parameters <- as.integer(parsed_names[,3])
  
  # Subset the matrix to keep only these columns
  filtered_mcmc <- mcmc_obj[, filtered_columns]
  posterior_mean <- colMeans(filtered_mcmc)
  
  mean_df <- data.frame(
    participant = participants,
    parameter = parameters,
    mean = posterior_mean
  )
  
  mean_df_wide <- pivot_wider(mean_df, names_from = parameter, values_from = mean, names_prefix = "beta_")
  return(mean_df_wide)
}

print(colnames(mcmc_md1))
# Export parameter estimates
dataframes <- list(
  m1 = pp_estimates_m1,
  m2 = pp_estimates_m2,
  m3 = pp_estimates_m3,
  m4 = pp_estimates_m4,
  m4star = pp_estimates_m4star,
  m5 = pp_estimates_m5
)

dataframes <- list(m2 = pp_estimates_m2)

dataframes_drift <- list(
  md1 = pp_estimates_md1,
  md2 = pp_estimates_md2,
  md3 = pp_estimates_md3,
  md4 = pp_estimates_md4,
  md4star = pp_estimates_md4star,
  md5 = pp_estimates_md5
)

# Write each dataframe to a separate CSV file
for (name in names(dataframes)) {
  write.csv(dataframes[[name]], paste0(name, "_parameter_estimates.csv"), row.names = FALSE)
}
for (name in names(dataframes_drift)) {
  write.csv(dataframes_drift[[name]], paste0(name, "_parameter_estimates.csv"), row.names = FALSE)
}
################################################Model 1###########################################
##################################################################################################
load("jags_sample_m1.RData")

#look into the traceplot of one variable
#traceplot(jagsModel_m1$mcmc[,"Mu[3]"])

mcmc_m1 = as.matrix(as.mcmc.list(jagsModel_m1$mcmc), chains = F)
pp_estimates_m1 <- participant_estimate_extractor(mcmc_m1)
quantiles_m1 <- apply(mcmc_m1, 2, function(x) quantile(x, c(0.025, 0.975)))
posterior_mean_m1 =  colMeans(mcmc_m1)
posterior_mean_m1["Mu[1]"] - 50
posterior_mean_m1["Mu[2]"]
posterior_mean_m1["Mu[3]"]
posterior_mean_m1["Sigma[1]"]
posterior_mean_m1["Sigma[2]"]
posterior_mean_m1["Sigma[3]"]
names(posterior_mean_m1)
#group level parameters plot
par(mfrow = c(2, 2), cex.main = 3, cex.axis=1.5)
options(repr.plot.width = 20, repr.plot.height = 8)  
hist(mcmc_m1[, "Mu[1]"], breaks = 100, main = expression(mu[bold(beta[0])]), xlab = expression(mu[bold(beta[0])]),col="skyblue",border ="skyblue")
hist(mcmc_m1[, "Mu[2]"], breaks = 100, main = expression(mu[bold(beta[1])]), xlab = expression(mu[bold(beta[1])]),col="skyblue",border ="skyblue")
hist(mcmc_m1[, "Mu[3]"], breaks = 100, main = expression(mu[bold(beta[phi])]), xlab = expression(mu[bold(beta[phi])]),col="skyblue",border ="skyblue")
mtext("M1", line = 0.5, outer = TRUE)
length(mcmc_m1[, "Mu[1]"])
#check convergence
#exclude ld to keep it fast
exclude_ld =  colnames(jagsModel_m1$mcmc[[1]])[grepl("ld", colnames(jagsModel_m1$mcmc[[1]]))]
ld_m1 = list()

for (i in 1:length(jagsModel_m1$mcmc)){
  ld_m1[[i]] = jagsModel_m1$mcmc[[i]][,colnames(jagsModel_m1$mcmc[[i]]) %in% exclude_ld]
  jagsModel_m1$mcmc[[i]] = jagsModel_m1$mcmc[[i]][,!colnames(jagsModel_m1$mcmc[[i]]) %in% exclude_ld]
}

R_hat_m1 = gelman.diag(jagsModel_m1$mcmc, multivariate = F)
mcmc_list = as.mcmc.list(jagsModel_m1)
length(which(as.numeric(R_hat_m1$psrf[,"Point est."])>1.10))

#calculate elpd
ld_m1_test = array(c(ld_m1[[1]],ld_m1[[2]],ld_m1[[3]]),dim=c(dim(ld_m1[[1]])[1],3,dim(ld_m1[[1]])[2]))
loo_m1 = loo(ld_m1_test)
pareto_k_table(loo_m1)
ld_m1 = apply(as.matrix(ld_m1), c(1, 2), mean)
################################################Model 2###########################################
##################################################################################################
load("jags_sample_m2.RData")
load("jags_sample_m2_nonlinear_simulation.Rdata")

#look into the traceplot of one variable
#traceplot(jagsModel_m2$mcmc[,"Mu[3]"])

mcmc_m2 = as.matrix(as.mcmc.list(jagsModel_m2$mcmc), chains = F)

pp_estimates_m2 <- participant_estimate_extractor(mcmc_m2)
quantiles_m2 <- apply(mcmc_m2, 2, function(x) quantile(x, c(0.025, 0.975)))
posterior_mean_m2 =  colMeans(mcmc_m2)
posterior_mean_m2["Mu[1]"] - 50  # beta_L0
posterior_mean_m2["Mu[1]"] + posterior_mean_m2["Mu[2]"] - 50  # beta_w0
posterior_mean_m2["Mu[3]"]  # beta_w1
posterior_mean_m2["Mu[4]"]  # beta_L1
posterior_mean_m2["Mu[5]"]  # beta_phi

#group level parameters plot
par(mfrow = c(3, 3), cex.main = 3, cex.axis=1.5)
options(repr.plot.width = 20, repr.plot.height = 8)  
hist(mcmc_m2[, "Mu[1]"] - 50, breaks = 100, main = expression(mu[bold(beta[L0])]), xlab = expression(mu[bold(beta[L0])]),col="skyblue",border ="skyblue")
hist(mcmc_m2[, "Mu[1]"]+mcmc_m2[, "Mu[2]"]-50, breaks = 100, main = expression(mu[bold(beta[w0])]), xlab = expression(mu[bold(beta[w0])]),col="skyblue",border ="skyblue")
hist(mcmc_m2[, "Mu[3]"], breaks = 100, main = expression(mu[bold(beta[w1])]), xlab = expression(mu[bold(beta[w1])]),col="skyblue",border ="skyblue")
hist(mcmc_m2[, "Mu[4]"], breaks = 100, main = expression(mu[bold(beta[L1])]), xlab = expression(mu[bold(beta[L1])]),col="skyblue",border ="skyblue")
hist(mcmc_m2[, "Mu[5]"], breaks = 100, main = expression(mu[bold(beta[phi])]), xlab = expression(mu[bold(beta[phi])]),col="skyblue",border ="skyblue")


#hypothesis
hist(abs(mcmc_m2[, "Mu[1]"]+mcmc_m2[, "Mu[2]"]-50)-abs(mcmc_m2[, "Mu[1]"]-50), breaks = 100, main = expression("|" *mu[bold(beta[w0])]*"|"-"|" *mu[bold(beta[L0])]*"|"), xlab = expression("|" *mu[bold(beta[w0])]*"|"-"|" *mu[bold(beta[L0])]*"|"), col="skyblue",border ="skyblue")
hist(abs(mcmc_m2[, "Mu[3]"])-abs(mcmc_m2[, "Mu[4]"]), breaks = 100, main = expression("|" *mu[bold(beta[w1])]*"|"-"|" *mu[bold(beta[L1])]*"|"), xlab = expression("|" *mu[bold(beta[w1])]*"|"-"|" *mu[bold(beta[L1])]*"|"), col="skyblue",border ="skyblue")
#check convergence
#exclude ld to keep it fast
exclude_ld =  colnames(jagsModel_m2$mcmc[[1]])[grepl("ld", colnames(jagsModel_m2$mcmc[[1]]))]
for (i in 1:length(jagsModel_m2$mcmc)){
  jagsModel_m2$mcmc[[i]] = jagsModel_m2$mcmc[[i]][,!colnames(jagsModel_m2$mcmc[[i]]) %in% exclude_ld]
}

R_hat_m2 = gelman.diag(jagsModel_m2$mcmc, multivariate = F)
length(which(as.numeric(R_hat_m2$psrf[,"Point est."])>1.10))


################################################Model 3###########################################
##################################################################################################
load("jags_sample_m3.RData")

#look into the traceplot of one variable
#traceplot(jagsModel_m3$mcmc[,"Mu[3]"])

mcmc_m3 = as.matrix(as.mcmc.list(jagsModel_m3$mcmc), chains = F)
pp_estimates_m3 <- participant_estimate_extractor(mcmc_m3)
#quantiles_m3 <- apply(mcmc_m3, 2, function(x) quantile(x, c(0.025, 0.975)))
posterior_mean_m3 =  colMeans(mcmc_m3)
posterior_mean_m3["Mu[1]"] - 50  # beta_L0
posterior_mean_m3["Mu[1]"] + posterior_mean_m2["Mu[2]"] - 50  # beta_w0
posterior_mean_m3["Mu[3]"]  # beta_w1
posterior_mean_m3["Mu[4]"]  # beta_w2
posterior_mean_m3["Mu[5]"]  # beta_w3
posterior_mean_m3["Mu[6]"]  # beta_L1
posterior_mean_m3["Mu[7]"]  # beta_L2
posterior_mean_m3["Mu[8]"]  # beta_L3
posterior_mean_m3["Mu[9]"]  # beta_phi

#group level parameters plot
par(mfrow = c(4, 3), cex.main = 3, cex.axis=1.5)
options(repr.plot.width = 20, repr.plot.height = 8)  
hist(mcmc_m3[, "Mu[1]"], breaks = 100, main = expression(mu[bold(beta[L0])]), xlab = expression(mu[bold(beta[L0])]),col="skyblue",border ="skyblue")
hist(mcmc_m3[, "Mu[1]"]+mcmc_m3[, "Mu[2]"], breaks = 100, main = expression(mu[bold(beta[w0])]), xlab = expression(mu[bold(beta[w0])]),col="skyblue",border ="skyblue")
hist(mcmc_m3[, "Mu[3]"], breaks = 100, main = expression(mu[bold(beta[w1])]), xlab = expression(mu[bold(beta[w1])]),col="skyblue",border ="skyblue")
hist(mcmc_m3[, "Mu[4]"], breaks = 100, main = expression(mu[bold(beta[w2])]), xlab = expression(mu[bold(beta[w2])]),col="skyblue",border ="skyblue")
hist(mcmc_m3[, "Mu[5]"], breaks = 100, main = expression(mu[bold(beta[w3])]), xlab = expression(mu[bold(beta[w3])]),col="skyblue",border ="skyblue")
hist(mcmc_m3[, "Mu[6]"], breaks = 100, main = expression(mu[bold(beta[L1])]), xlab = expression(mu[bold(beta[L1])]),col="skyblue",border ="skyblue")
hist(mcmc_m3[, "Mu[7]"], breaks = 100, main = expression(mu[bold(beta[L2])]), xlab = expression(mu[bold(beta[L2])]),col="skyblue",border ="skyblue")
hist(mcmc_m3[, "Mu[8]"], breaks = 100, main = expression(mu[bold(beta[L3])]), xlab = expression(mu[bold(beta[L3])]),col="skyblue",border ="skyblue")
hist(mcmc_m3[, "Mu[9]"], breaks = 100, main = expression(mu[bold(beta[phi])]), xlab = expression(mu[bold(beta[phi])]),col="skyblue",border ="skyblue")

#hypothesis
hist(abs(mcmc_m3[, "Mu[4]"])-abs(mcmc_m3[, "Mu[7]"]), breaks = 100, main = expression("|" *mu[bold(beta[w2])]*"|"-"|" *mu[bold(beta[L2])]*"|"), xlab = expression("|" *mu[bold(beta[w2])]*"|"-"|" *mu[bold(beta[L2])]*"|"), col="skyblue",border ="skyblue")
hist(abs(mcmc_m3[, "Mu[5]"])-abs(mcmc_m3[, "Mu[8]"]), breaks = 100, main = expression("|" *mu[bold(beta[w3])]*"|"-"|" *mu[bold(beta[L3])]*"|"), xlab = expression("|" *mu[bold(beta[w3])]*"|"-"|" *mu[bold(beta[L3])]*"|"), col="skyblue",border ="skyblue")
#check convergence
#exclude ld to keep it fast
exclude_ld =  colnames(jagsModel_m3$mcmc[[1]])[grepl("ld", colnames(jagsModel_m3$mcmc[[1]]))]
for (i in 1:length(jagsModel_m3$mcmc)){
  jagsModel_m3$mcmc[[i]] = jagsModel_m3$mcmc[[i]][,!colnames(jagsModel_m3$mcmc[[i]]) %in% exclude_ld]
}

R_hat_m3 = gelman.diag(jagsModel_m3$mcmc, multivariate = F)
length(which(as.numeric(R_hat_m3$psrf[,"Point est."])>1.10))


################################################Model 4###########################################
##################################################################################################
load("jags_sample_m4.RData")

#look into the traceplot of one variable
#traceplot(jagsModel_m4$mcmc[,"Mu[3]"])

mcmc_m4 = as.matrix(as.mcmc.list(jagsModel_m4$mcmc), chains = F)
pp_estimates_m4 <- participant_estimate_extractor(mcmc_m4)
quantiles_m4 <- apply(mcmc_m4, 2, function(x) quantile(x, c(0.025, 0.975)))
posterior_mean_m4 =  colMeans(mcmc_m4)
posterior_mean_m4["Mu[1]"]  # beta_0
posterior_mean_m4["Mu[2]"]  # beta_1
2*(1/(1+exp(-posterior_mean_m4["Mu[3]"])))  # gamma
posterior_mean_m4["Mu[4]"]  # beta_phi

quantiles_m4[, "Mu[1]"] - 50  # beta_0
quantiles_m4[, "Mu[2]"]  # beta_1
2*(1/(1+exp(-quantiles_m4[, "Mu[3]"])))  # gamma
quantiles_m4[, "Mu[4]"]  # beta_phi

#group level parameters plot
par(mfrow = c(2, 2), cex.main = 3, cex.axis=1.5)
options(repr.plot.width = 20, repr.plot.height = 8)  
hist(mcmc_m4[, "Mu[1]"], breaks = 100, main = expression(mu[bold(beta[0])]), xlab = expression(mu[bold(beta[0])]))
hist(mcmc_m4[, "Mu[2]"], breaks = 100, main = expression(mu[bold(beta[1])]), xlab = expression(mu[bold(beta[1])]))
hist(2*(1/(1+exp(-mcmc_m4[, "Mu[3]"]))), breaks = 100, main = expression(mu[bold(gamma)]), xlab = expression(mu[bold(gamma)]))
hist(mcmc_m4[, "Mu[4]"], breaks = 100, main = expression(mu[bold(beta[phi])]), xlab = expression(mu[bold(beta[phi])]))

#check convergence
#exclude ld to keep it fast
exclude_ld =  colnames(jagsModel_m4$mcmc[[1]])[grepl("ld", colnames(jagsModel_m4$mcmc[[1]]))]
ld_m4 = list()
for (i in 1:length(jagsModel_m4$mcmc)){
  ld_m4[[i]] = jagsModel_m4$mcmc[[i]][,colnames(jagsModel_m4$mcmc[[i]]) %in% exclude_ld]
  jagsModel_m4$mcmc[[i]] = jagsModel_m4$mcmc[[i]][,!colnames(jagsModel_m4$mcmc[[i]]) %in% exclude_ld]
}


R_hat_m4 = gelman.diag(jagsModel_m4$mcmc, multivariate = F)
length(which(as.numeric(R_hat_m4$psrf[,"Point est."])>1.10))

#calculate elpd
ld_m4_test = array(c(ld_m4[[1]],ld_m4[[2]],ld_m4[[3]]),dim=c(dim(ld_m4[[1]])[1],3,dim(ld_m4[[1]])[2]))
loo_m4 = loo(ld_m4_test)
pareto_k_table(loo_m1)
ld_m1 = apply(as.matrix(ld_m1), c(1, 2), mean)

################################################Model 4star###########################################
##################################################################################################
load("jags_sample_m4star.RData")

#look into the traceplot of one variable
#traceplot(jagsModel_m4star$mcmc[,"Mu[3]"])

mcmc_m4star = as.matrix(as.mcmc.list(jagsModel_m4star$mcmc), chains = F)
pp_estimates_m4star <- participant_estimate_extractor(mcmc_m4star)
#quantiles_m4star <- apply(mcmc_m4star, 2, function(x) quantile(x, c(0.025, 0.975)))
#posterior_mean_m4star =  colMeans(mcmc_m4star)

#group level parameters plot
par(mfrow = c(2, 2), cex.main = 3, cex.axis=1.5)
options(repr.plot.width = 20, repr.plot.height = 8)  
hist(mcmc_m4star[, "Mu[1]"], breaks = 100, main = expression(mu[bold(beta[0])]), xlab = expression(mu[bold(beta[0])]))
hist(mcmc_m4star[, "Mu[2]"], breaks = 100, main = expression(mu[bold(beta[1])]), xlab = expression(mu[bold(beta[1])]))
hist(mcmc_m4star[, "Mu[3]"], breaks = 100, main = expression(mu[bold(beta[2])]), xlab = expression(mu[bold(beta[2])]))

#check convergence
#exclude ld to keep it fast
exclude_ld =  colnames(jagsModel_m4star$mcmc[[1]])[grepl("ld", colnames(jagsModel_m4star$mcmc[[1]]))]
for (i in 1:length(jagsModel_m4star$mcmc)){
  jagsModel_m4star$mcmc[[i]] = jagsModel_m4star$mcmc[[i]][,!colnames(jagsModel_m4star$mcmc[[i]]) %in% exclude_ld]
}

#calculate elpd
ld_m1_test = array(c(ld_m1[[1]],ld_m1[[2]],ld_m1[[3]]),dim=c(dim(ld_m1[[1]])[1],3,dim(ld_m1[[1]])[2]))
loo_m1 = loo(ld_m1_test)

R_hat_m4star = gelman.diag(jagsModel_m4star$mcmc, multivariate = F)
length(which(as.numeric(R_hat_m4star$psrf[,"Point est."])>1.10))


################################################Model 5###########################################
##################################################################################################
load("jags_sample_m5.RData")

#look into the traceplot of one variable
#traceplot(jagsModel_m5$mcmc[,"Mu[3]"])

mcmc_m5 = as.matrix(as.mcmc.list(jagsModel_m5$mcmc), chains = F)
pp_estimates_m5 <- participant_estimate_extractor(mcmc_m5)
#quantiles_m5 <- apply(mcmc_m5, 2, function(x) quantile(x, c(0.025, 0.975)))
#posterior_mean_m5 =  colMeans(mcmc_m5)

#group level parameters plot
par(mfrow = c(2, 2), cex.main = 3, cex.axis=1.5)
options(repr.plot.width = 20, repr.plot.height = 8)  
hist(mcmc_m5[, "Mu[1]"], breaks = 100, main = expression(mu[bold(beta[0])]), xlab = expression(mu[bold(beta[0])]))
hist(mcmc_m5[, "Mu[2]"], breaks = 100, main = expression(mu[bold(beta[1])]), xlab = expression(mu[bold(beta[1])]))
hist(mcmc_m5[, "Mu[3]"], breaks = 100, main = expression(mu[bold(beta[2])]), xlab = expression(mu[bold(beta[2])]))

#check convergence
#exclude ld to keep it fast
exclude_ld =  colnames(jagsModel_m5$mcmc[[1]])[grepl("ld", colnames(jagsModel_m5$mcmc[[1]]))]

ld_m5 = list()
for (i in 1:length(jagsModel_m5$mcmc)){
  ld_m5[[i]] = jagsModel_m5$mcmc[[i]][,colnames(jagsModel_m5$mcmc[[i]]) %in% exclude_ld]
  jagsModel_m5$mcmc[[i]] = jagsModel_m5$mcmc[[i]][,!colnames(jagsModel_m5$mcmc[[i]]) %in% exclude_ld]
}

R_hat_m5 = gelman.diag(jagsModel_m5$mcmc, multivariate = F)
length(which(as.numeric(R_hat_m5$psrf[,"Point est."])>1.10))


################################################Model D1###########################################
##################################################################################################
load("jags_sample_md1.RData")

#look into the traceplot of one variable
traceplot(jagsModel_md1$mcmc[,"MuZ[3]"])

mcmc_md1 = as.matrix(as.mcmc.list(jagsModel_md1$mcmc), chains = F)
colnames(mcmc_md1)[1500:1550]
pp_estimates_md1 <- participant_estimate_extractor_drift(mcmc_md1)

quantiles_md1 <- apply(mcmc_md1, 2, function(x) quantile(x, c(0.025, 0.975)))
posterior_mean_md1 =  colMeans(mcmc_md1)

posterior_mean_md1["MuV[1]"]  # beta_0
posterior_mean_md1["MuV[2]"]  # beta_1
posterior_mean_md1["MuV[3]"]  # beta_phi
posterior_mean_md1["MuZ[1]"]  # z0
posterior_mean_md1["MuZ[2]"]  # z2
posterior_mean_md1["MuZ[3]"]  # z3

quantiles_md1[, "MuV[1]"]  # beta_0
quantiles_md1[, "MuV[2]"]  # beta_1
quantiles_md1[, "MuV[3]"]  # beta_phi
quantiles_md1[, "MuZ[1]"]  # z0
quantiles_md1[, "MuZ[2]"]  # z1
quantiles_md1[, "MuZ[3]"]  # z2

#group level parameters plot
par(mfrow = c(3, 3), cex.main = 3, cex.axis=1.5)
options(repr.plot.width = 20, repr.plot.height = 8)  
hist(mcmc_md1[, "MuZ[1]"], breaks = 100, main = expression(mu[bold(z[0])]), xlab = expression(mu[bold(z[0])]),col="skyblue",border ="skyblue")
hist(mcmc_md1[, "MuZ[2]"], breaks = 100, main = expression(mu[bold(z[1])]), xlab = expression(mu[bold(z[1])]),col="skyblue",border ="skyblue")
hist(mcmc_md1[, "MuZ[3]"], breaks = 100, main = expression(mu[bold(z[phi])]), xlab = expression(mu[bold(z[phi])]),col="skyblue",border ="skyblue")
hist(mcmc_md1[, "MuV[1]"], breaks = 100, main = expression(mu[bold(v[0])]), xlab = expression(mu[bold(v[0])]),col="skyblue",border ="skyblue")
hist(mcmc_md1[, "MuV[2]"], breaks = 100, main = expression(mu[bold(v[1])]), xlab = expression(mu[bold(v[1])]),col="skyblue",border ="skyblue")
hist(mcmc_md1[, "MuV[3]"], breaks = 100, main = expression(mu[bold(v[phi])]), xlab = expression(mu[bold(v[phi])]),col="skyblue",border ="skyblue")
hist(mcmc_md1[, "mualpha"], breaks = 100, main = expression(mu[bold(alpha)]), xlab = expression(mu[bold(alpha)]),col="skyblue",border ="skyblue")
hist(mcmc_md1[, "mutau"], breaks = 100, main = expression(mu[bold(tau)]), xlab = expression(mu[bold(tau)]),col="skyblue",border ="skyblue")

#check convergence
#exclude ld to keep it fast
exclude_ld =  colnames(jagsModel_md1$mcmc[[1]])[grepl("ld", colnames(jagsModel_md1$mcmc[[1]]))]
ld_md1 = list()

for (i in 1:length(jagsModel_md1$mcmc)){
  ld_md1[[i]] = jagsModel_md1$mcmc[[i]][,colnames(jagsModel_md1$mcmc[[i]]) %in% exclude_ld]
  jagsModel_md1$mcmc[[i]] = jagsModel_md1$mcmc[[i]][,!colnames(jagsModel_md1$mcmc[[i]]) %in% exclude_ld]
}

#R_hat_md1 = gelman.diag(jagsModel_md1$mcmc, multivariate = F)
#mcmc_list = as.mcmc.list(jagsModel_md1)
#length(which(as.numeric(R_hat_md1$psrf[,"Point est."])>1.10))

#calculate elpd
ld_md1_test = array(c(ld_md1[[1]],ld_md1[[2]],ld_md1[[3]]),dim=c(dim(ld_md1[[1]])[1],3,dim(ld_md1[[1]])[2]))
loo_md1 = loo(ld_md1_test)
#pareto_k_table(loo_md1)
#ld_md1 = apply(as.matrix(ld_md1), c(1, 2), mean)

loo_compare(loo_md1,loo_md4)
################################################Model D2###########################################
##################################################################################################
#load("jags_sample_md2.RData")
#load("jags_sample_MD2_scale_dichotomized.RData")
load("jags_sample_MD2.RData")

#look into the traceplot of one variable
traceplot(jagsModel_md2$mcmc[,"MuV[5]"])

mcmc_md2 = as.matrix(as.mcmc.list(jagsModel_md2$mcmc), chains = F)
pp_estimates_md2 <- participant_estimate_extractor_drift(mcmc_md2)
quantiles_md2 <- apply(mcmc_md2, 2, function(x) quantile(x, c(0.025, 0.975)))
posterior_mean_md2 =  colMeans(mcmc_md2)


posterior_mean_md2["MuV[1]"]  # beta_L0
posterior_mean_md2["MuV[1]"] + posterior_mean_md2["MuV[2]"]  # beta_w0
posterior_mean_md2["MuV[3]"]  # beta_w1
posterior_mean_md2["MuV[4]"]  # beta_L1
posterior_mean_md2["MuV[5]"]  # beta_phi
posterior_mean_md2["MuZ[1]"]  # z0
posterior_mean_md2["MuZ[2]"]  # z2
posterior_mean_md2["MuZ[3]"]  # z3

quantiles_md2[, "MuV[1]"]  # beta_L0
quantiles_md2[, "MuV[1]"] + quantiles_md2[, "MuV[2]"]  # beta_w0
quantiles_md2[, "MuV[3]"]  # beta_w1
quantiles_md2[, "MuV[4]"]  # beta_L1
quantiles_md2[, "MuV[5]"]  # beta_phi
quantiles_md2[, "MuZ[1]"]  # z0
quantiles_md2[, "MuZ[2]"]  # z1
quantiles_md2[, "MuZ[3]"]  # z2

#group level parameters plot
par(mfrow = c(4, 3), cex.main = 3, cex.axis=1.5)
options(repr.plot.width = 20, repr.plot.height = 8)  
hist(mcmc_md2[, "MuZ[1]"], breaks = 100, main = expression(mu[bold(z[0])]), xlab = expression(mu[bold(z[0])]),col="skyblue",border ="skyblue")
hist(mcmc_md2[, "MuZ[2]"], breaks = 100, main = expression(mu[bold(z[1])]), xlab = expression(mu[bold(z[1])]),col="skyblue",border ="skyblue")
hist(mcmc_md2[, "MuZ[3]"], breaks = 100, main = expression(mu[bold(z[phi])]), xlab = expression(mu[bold(z[phi])]),col="skyblue",border ="skyblue")
hist(mcmc_md2[, "MuV[1]"], breaks = 100, main = expression(mu[bold(v[L0])]), xlab = expression(mu[bold(v[L0])]),col="skyblue",border ="skyblue")
hist(mcmc_md2[, "MuV[1]"]+mcmc_md2[, "MuV[2]"], breaks = 100, main = expression(mu[bold(v[w0])]), xlab = expression(mu[bold(v[w0])]),col="skyblue",border ="skyblue")
hist(mcmc_md2[, "MuV[3]"], breaks = 100, main = expression(mu[bold(v[w1])]), xlab = expression(mu[bold(v[w1])]),col="skyblue",border ="skyblue")
hist(mcmc_md2[, "MuV[4]"], breaks = 100, main = expression(mu[bold(v[L1])]), xlab = expression(mu[bold(v[L1])]),col="skyblue",border ="skyblue")
hist(mcmc_md2[, "MuV[5]"], breaks = 100, main = expression(mu[bold(v[phi])]), xlab = expression(mu[bold(v[phi])]),col="skyblue",border ="skyblue")
hist(mcmc_md2[, "mualpha"], breaks = 100, main = expression(mu[bold(alpha)]), xlab = expression(mu[bold(alpha)]),col="skyblue",border ="skyblue")
hist(mcmc_md2[, "mutau"], breaks = 100, main = expression(mu[bold(tau)]), xlab = expression(mu[bold(tau)]),col="skyblue",border ="skyblue")


hist(abs(mcmc_md2[, "MuV[1]"]+mcmc_md2[, "MuV[2]"])-abs(mcmc_md2[, "MuV[1]"]), breaks = 100, main = expression("|" *mu[bold(v[w0])]*"|"-"|" *mu[bold(v[L0])]*"|"), xlab = expression("|" *mu[bold(v[w0])]*"|"-"|" *mu[bold(v[L0])]*"|"), col="skyblue",border ="skyblue")
hist(abs(mcmc_md2[, "MuV[3]"])-abs(mcmc_md2[, "MuV[4]"]), breaks = 100, main = expression("|" *mu[bold(v[w1])]*"|"-"|" *mu[bold(v[L1])]*"|"), xlab = expression("|" *mu[bold(v[w1])]*"|"-"|" *mu[bold(v[L1])]*"|"), col="skyblue",border ="skyblue")

quantile(abs(mcmc_md2[, "MuV[1]"]+mcmc_md2[, "MuV[2]"])-abs(mcmc_md2[, "MuV[1]"]), c(0.025, 0.975))
quantile(abs(mcmc_md2[, "MuV[3]"])-abs(mcmc_md2[, "MuV[4]"]), c(0.025, 0.975))


#check convergence
#exclude ld to keep it fast
exclude_ld =  colnames(jagsModel_md2$mcmc[[1]])[grepl("ld", colnames(jagsModel_md2$mcmc[[1]]))]
ld_md2 = list()

for (i in 1:length(jagsModel_md2$mcmc)){
  ld_md2[[i]] = jagsModel_md2$mcmc[[i]][,colnames(jagsModel_md2$mcmc[[i]]) %in% exclude_ld]
  jagsModel_md2$mcmc[[i]] = jagsModel_md2$mcmc[[i]][,!colnames(jagsModel_md2$mcmc[[i]]) %in% exclude_ld]
}

R_hat_md2 = gelman.diag(jagsModel_md2$mcmc, multivariate = F)
mcmc_list = as.mcmc.list(jagsModel_md2)
length(which(as.numeric(R_hat_md2$psrf[,"Point est."])>1.10))

#calculate elpd
#ld_md2_test = array(c(ld_md2[[1]],ld_md2[[2]],ld_md2[[3]]),dim=c(dim(ld_md2[[1]])[1],3,dim(ld_md2[[1]])[2]))
#loo_md2 = loo(ld_md2_test)
#pareto_k_table(loo_md2)
#ld_md2 = apply(as.matrix(ld_md2), c(1, 2), mean)

########################################################################################
#######################Model D3#########################################################
load("jags_sample_md3_Tomi.RData")
load("jags_sample_md3_v2.RData")

#look into the traceplot of one variable
traceplot(jagsModel_md3$mcmc[,"MuV[2]"])

mcmc_md3 = as.matrix(as.mcmc.list(jagsModel_md3$mcmc), chains = F)
pp_estimates_md3 <- participant_estimate_extractor_drift(mcmc_md3)
quantiles_md3 <- apply(mcmc_md3, 2, function(x) quantile(x, c(0.025, 0.975)))
posterior_mean_md3 =  colMeans(mcmc_md3)

posterior_mean_md3["MuV[1]"]  # beta_L0
posterior_mean_md3["MuV[1]"] + posterior_mean_md3["MuV[2]"]  # beta_w0
posterior_mean_md3["MuV[3]"]  # beta_w1
posterior_mean_md3["MuV[4]"]  # beta_w2
posterior_mean_md3["MuV[5]"]  # beta_w3
posterior_mean_md3["MuV[6]"]  # beta_L1
posterior_mean_md3["MuV[7]"]  # beta_L2
posterior_mean_md3["MuV[8]"]  # beta_L3
posterior_mean_md3["MuV[9]"]  # beta_phi
posterior_mean_md3["MuZ[1]"]  # z0
posterior_mean_md3["MuZ[2]"]  # z2
posterior_mean_md3["MuZ[3]"]  # z3

quantiles_md3[, "MuV[1]"]  # beta_L0
quantiles_md3[, "MuV[1]"] + quantiles_md2[, "MuV[2]"]  # beta_w0
quantiles_md3[, "MuV[3]"]  # beta_w1
quantiles_md3[, "MuV[4]"]  # beta_w2
quantiles_md3[, "MuV[5]"]  # beta_w3
quantiles_md3[, "MuV[6]"]  # beta_L1
quantiles_md3[, "MuV[7]"]  # beta_L2
quantiles_md3[, "MuV[8]"]  # beta_L3
quantiles_md3[, "MuV[9]"]  # beta_phi
quantiles_md3[, "MuZ[1]"]  # z0
quantiles_md3[, "MuZ[2]"]  # z1
quantiles_md3[, "MuZ[3]"]  # z2


#group level parameters plot
par(mfrow = c(5, 3), cex.main = 3, cex.axis=1.5)
options(repr.plot.width = 20, repr.plot.height = 8)  
hist(mcmc_md3[, "MuZ[1]"], breaks = 100, main = expression(mu[bold(z[0])]), xlab = expression(mu[bold(z[0])]),col="skyblue",border ="skyblue")
hist(mcmc_md3[, "MuZ[2]"], breaks = 100, main = expression(mu[bold(z[1])]), xlab = expression(mu[bold(z[1])]),col="skyblue",border ="skyblue")
hist(mcmc_md3[, "MuZ[3]"], breaks = 100, main = expression(mu[bold(z[phi])]), xlab = expression(mu[bold(z[phi])]),col="skyblue",border ="skyblue")
hist(mcmc_md3[, "MuV[1]"], breaks = 100, main = expression(mu[bold(v[L0])]), xlab = expression(mu[bold(v[L0])]),col="skyblue",border ="skyblue")
hist(mcmc_md3[, "MuV[1]"]+mcmc_md3[, "MuV[2]"], breaks = 100, main = expression(mu[bold(v[w0])]), xlab = expression(mu[bold(v[w0])]),col="skyblue",border ="skyblue")
hist(mcmc_md3[, "MuV[3]"], breaks = 100, main = expression(mu[bold(v[w1])]), xlab = expression(mu[bold(v[w1])]),col="skyblue",border ="skyblue")
hist(mcmc_md3[, "MuV[4]"], breaks = 100, main = expression(mu[bold(v[w2])]), xlab = expression(mu[bold(v[w2])]),col="skyblue",border ="skyblue")
hist(mcmc_md3[, "MuV[5]"], breaks = 100, main = expression(mu[bold(v[w3])]), xlab = expression(mu[bold(v[w3])]),col="skyblue",border ="skyblue")
hist(mcmc_md3[, "MuV[6]"], breaks = 100, main = expression(mu[bold(v[L1])]), xlab = expression(mu[bold(v[L1])]),col="skyblue",border ="skyblue")
hist(mcmc_md3[, "MuV[7]"], breaks = 100, main = expression(mu[bold(v[L2])]), xlab = expression(mu[bold(v[L2])]),col="skyblue",border ="skyblue")
hist(mcmc_md3[, "MuV[8]"], breaks = 100, main = expression(mu[bold(v[L3])]), xlab = expression(mu[bold(v[L3])]),col="skyblue",border ="skyblue")
hist(mcmc_md3[, "MuV[9]"], breaks = 100, main = expression(mu[bold(v[phi])]), xlab = expression(mu[bold(v[phi])]),col="skyblue",border ="skyblue")
hist(mcmc_md3[, "mualpha"], breaks = 100, main = expression(mu[bold(alpha)]), xlab = expression(mu[bold(alpha)]),col="skyblue",border ="skyblue")
hist(mcmc_md3[, "mutau"], breaks = 100, main = expression(mu[bold(tau)]), xlab = expression(mu[bold(tau)]),col="skyblue",border ="skyblue")

#hypothesis
hist(abs(mcmc_md3[, "MuV[4]"])-abs(mcmc_md3[, "MuV[7]"]), breaks = 100, main = expression("|" *mu[bold(v[w2])]*"|"-"|" *mu[bold(v[L2])]*"|"), xlab = expression("|" *mu[bold(v[w2])]*"|"-"|" *mu[bold(v[L2])]*"|"), col="skyblue",border ="skyblue")
hist(abs(mcmc_md3[, "MuV[5]"])-abs(mcmc_md3[, "MuV[8]"]), breaks = 100, main = expression("|" *mu[bold(v[w3])]*"|"-"|" *mu[bold(v[L3])]*"|"), xlab = expression("|" *mu[bold(v[w3])]*"|"-"|" *mu[bold(v[L3])]*"|"), col="skyblue",border ="skyblue")

exclude_ld =  colnames(jagsModel_md3$mcmc[[1]])[grepl("ld", colnames(jagsModel_md3$mcmc[[1]]))]
ld_md3 = list()

for (i in 1:length(jagsModel_md3$mcmc)){
  ld_md3[[i]] = jagsModel_md3$mcmc[[i]][,colnames(jagsModel_md3$mcmc[[i]]) %in% exclude_ld]
  jagsModel_md3$mcmc[[i]] = jagsModel_md3$mcmc[[i]][,!colnames(jagsModel_md3$mcmc[[i]]) %in% exclude_ld]
}

R_hat_md3 = gelman.diag(jagsModel_md3$mcmc, multivariate = F)
#mcmc_list = as.mcmc.list(jagsModel_md3)
length(which(as.numeric(R_hat_md3$psrf[,"Point est."])>1.10))

################################################Model D4###########################################
##################################################################################################
load("jags_sample_md4.RData")

#look into the traceplot of one variable
traceplot(jagsModel_md4$mcmc[,"MuZ[3]"])

mcmc_md4 = as.matrix(as.mcmc.list(jagsModel_md4$mcmc), chains = F)
pp_estimates_md4 <- participant_estimate_extractor_drift(mcmc_md4)
quantiles_md4 <- apply(mcmc_md4, 2, function(x) quantile(x, c(0.025, 0.975)))
posterior_mean_md4 =  colMeans(mcmc_md4)

#group level parameters plot
par(mfrow = c(3, 3), cex.main = 3, cex.axis=1.5)
options(repr.plot.width = 20, repr.plot.height = 8)  
hist(mcmc_md4[, "MuZ[1]"], breaks = 100, main = expression(mu[bold(z[0])]), xlab = expression(mu[bold(z[0])]),col="skyblue",border ="skyblue")
hist(mcmc_md4[, "MuZ[2]"], breaks = 100, main = expression(mu[bold(z[1])]), xlab = expression(mu[bold(z[1])]),col="skyblue",border ="skyblue")
hist(mcmc_md4[, "MuZ[3]"], breaks = 100, main = expression(mu[bold(z[phi])]), xlab = expression(mu[bold(z[phi])]),col="skyblue",border ="skyblue")
hist(mcmc_md4[, "MuV[1]"], breaks = 100, main = expression(mu[bold(v[0])]), xlab = expression(mu[bold(v[0])]),col="skyblue",border ="skyblue")
hist(mcmc_md4[, "MuV[2]"], breaks = 100, main = expression(mu[bold(v[1])]), xlab = expression(mu[bold(v[1])]),col="skyblue",border ="skyblue")
hist(mcmc_md4[, "MuV[3]"], breaks = 100, main = expression(mu[bold(v[phi])]), xlab = expression(mu[bold(v[phi])]),col="skyblue",border ="skyblue")
hist(mcmc_md4[, "mualpha"], breaks = 100, main = expression(mu[bold(alpha)]), xlab = expression(mu[bold(alpha)]),col="skyblue",border ="skyblue")
hist(mcmc_md4[, "mutau"], breaks = 100, main = expression(mu[bold(tau)]), xlab = expression(mu[bold(tau)]),col="skyblue",border ="skyblue")
hist(2* (1/(1+exp(-mcmc_md4[, "mugamma"]))), breaks = 100, main = expression(mu[bold(gamma)]), xlab = expression(mu[bold(gamma)]),col="skyblue",border ="skyblue")


#check convergence
#exclude ld to keep it fast
exclude_ld =  colnames(jagsModel_md4$mcmc[[1]])[grepl("ld", colnames(jagsModel_md4$mcmc[[1]]))]
ld_md4 = list()

for (i in 1:length(jagsModel_md4$mcmc)){
  ld_md4[[i]] = jagsModel_md4$mcmc[[i]][,colnames(jagsModel_md4$mcmc[[i]]) %in% exclude_ld]
  jagsModel_md4$mcmc[[i]] = jagsModel_md4$mcmc[[i]][,!colnames(jagsModel_md4$mcmc[[i]]) %in% exclude_ld]
}

R_hat_md4 = gelman.diag(jagsModel_md4$mcmc, multivariate = F)
length(which(as.numeric(R_hat_md4$psrf[,"Point est."])>1.10))

#calculate elpd
ld_md4_test = array(c(ld_md4[[1]],ld_md4[[2]],ld_md4[[3]]),dim=c(dim(ld_md4[[1]])[1],3,dim(ld_md4[[1]])[2]))
loo_md4 = loo(ld_md4_test)
#pareto_k_table(loo_md4)
#ld_md4 = apply(as.matrix(ld_md4_test), c(1, 2), mean)


################################################Model D4star###########################################
##################################################################################################
load("jags_sample_md4star.RData")

#look into the traceplot of one variable
traceplot(jagsModel_md4star$mcmc[,"MuZ[3]"])

mcmc_md4star = as.matrix(as.mcmc.list(jagsModel_md4star$mcmc), chains = F)
pp_estimates_md4star <- participant_estimate_extractor_drift(mcmc_md4star)
quantiles_md4star <- apply(mcmc_md4star, 2, function(x) quantile(x, c(0.025, 0.975)))
posterior_mean_md4star =  colMeans(mcmc_md4star)

#group level parameters plot
par(mfrow = c(3, 3), cex.main = 3, cex.axis=1.5)
options(repr.plot.width = 20, repr.plot.height = 8)  
hist(mcmc_md4star[, "MuZ[1]"], breaks = 100, main = expression(mu[bold(z[0])]), xlab = expression(mu[bold(z[0])]))
hist(mcmc_md4star[, "MuZ[2]"], breaks = 100, main = expression(mu[bold(z[1])]), xlab = expression(mu[bold(z[1])]))
hist(mcmc_md4star[, "MuZ[3]"], breaks = 100, main = expression(mu[bold(z[phi])]), xlab = expression(mu[bold(z[phi])]))
hist(mcmc_md4star[, "MuV[1]"], breaks = 100, main = expression(mu[bold(v[0])]), xlab = expression(mu[bold(v[0])]))
hist(mcmc_md4star[, "MuV[2]"], breaks = 100, main = expression(mu[bold(v[1])]), xlab = expression(mu[bold(v[1])]))
hist(mcmc_md4star[, "MuV[3]"], breaks = 100, main = expression(mu[bold(v[phi])]), xlab = expression(mu[bold(v[phi])]))
hist(mcmc_md4star[, "mualpha"], breaks = 100, main = expression(mu[bold(alpha)]), xlab = expression(mu[bold(alpha)]))
hist(mcmc_md4star[, "mutau"], breaks = 100, main = expression(mu[bold(tau)]), xlab = expression(mu[bold(tau)]))


#check convergence
#exclude ld to keep it fast
exclude_ld =  colnames(jagsModel_md4star$mcmc[[1]])[grepl("ld", colnames(jagsModel_md4star$mcmc[[1]]))]
ld_md4star = list()

for (i in 1:length(jagsModel_md4star$mcmc)){
  #ld_md4star[[i]] = jagsModel_md4star$mcmc[[i]][,colnames(jagsModel_md4star$mcmc[[i]]) %in% exclude_ld]
  jagsModel_md4star$mcmc[[i]] = jagsModel_md4star$mcmc[[i]][,!colnames(jagsModel_md4star$mcmc[[i]]) %in% exclude_ld]
}

R_hat_md4star = gelman.diag(jagsModel_md4star$mcmc, multivariate = F)
length(which(as.numeric(R_hat_md4star$psrf[,"Point est."])>1.10))

#calculate elpd
ld_md4star_test = array(c(ld_md4star[[1]],ld_md4star[[2]],ld_md4star[[3]]),dim=c(dim(ld_md4star[[1]])[1],3,dim(ld_md4star[[1]])[2]))
loo_md4star = loo(ld_md4star_test)
#pareto_k_table(loo_md4star)
#ld_md4star = apply(as.matrix(ld_md4star), c(1, 2), mean)

################################################Model D5###########################################
##################################################################################################
load("jags_sample_md5.RData")

#look into the traceplot of one variable
traceplot(jagsModel_md5$mcmc[,"MuZ[3]"])

mcmc_md5 = as.matrix(as.mcmc.list(jagsModel_md5$mcmc), chains = F)
pp_estimates_md5 <- participant_estimate_extractor_drift(mcmc_md5)
quantiles_md5 <- apply(mcmc_md5, 2, function(x) quantile(x, c(0.025, 0.975)))
posterior_mean_md5 =  colMeans(mcmc_md5)

#group level parameters plot
par(mfrow = c(3, 3), cex.main = 3, cex.axis=1.5)
options(repr.plot.width = 20, repr.plot.height = 8)  
hist(mcmc_md5[, "MuZ[1]"], breaks = 100, main = expression(mu[bold(z[0])]), xlab = expression(mu[bold(z[0])]))
hist(mcmc_md5[, "MuZ[2]"], breaks = 100, main = expression(mu[bold(z[1])]), xlab = expression(mu[bold(z[1])]))
hist(mcmc_md5[, "MuZ[3]"], breaks = 100, main = expression(mu[bold(z[phi])]), xlab = expression(mu[bold(z[phi])]))
hist(mcmc_md5[, "MuV[1]"], breaks = 100, main = expression(mu[bold(v[0])]), xlab = expression(mu[bold(v[0])]))
hist(mcmc_md5[, "MuV[2]"], breaks = 100, main = expression(mu[bold(v[1])]), xlab = expression(mu[bold(v[1])]))
hist(mcmc_md5[, "MuV[3]"], breaks = 100, main = expression(mu[bold(v[phi])]), xlab = expression(mu[bold(v[phi])]))
hist(mcmc_md5[, "mualpha"], breaks = 100, main = expression(mu[bold(alpha)]), xlab = expression(mu[bold(alpha)]))
hist(mcmc_md5[, "mutau"], breaks = 100, main = expression(mu[bold(tau)]), xlab = expression(mu[bold(tau)]))


#check convergence
#exclude ld to keep it fast
exclude_ld =  colnames(jagsModel_md5$mcmc[[1]])[grepl("ld", colnames(jagsModel_md5$mcmc[[1]]))]
ld_md5 = list()

for (i in 1:length(jagsModel_md5$mcmc)){
  ld_md5[[i]] = jagsModel_md5$mcmc[[i]][,colnames(jagsModel_md5$mcmc[[i]]) %in% exclude_ld]
  jagsModel_md5$mcmc[[i]] = jagsModel_md5$mcmc[[i]][,!colnames(jagsModel_md5$mcmc[[i]]) %in% exclude_ld]
}

R_hat_md5 = gelman.diag(jagsModel_md5$mcmc, multivariate = F)
length(which(as.numeric(R_hat_md5$psrf[,"Point est."])>1.10))

#calculate elpd
ld_md5_test = array(c(ld_md5[[1]],ld_md5[[2]],ld_md5[[3]]),dim=c(dim(ld_md5[[1]])[1],3,dim(ld_md5[[1]])[2]))
loo_md5 = loo(ld_md5_test)
#pareto_k_table(loo_md5)
#ld_md5 = apply(as.matrix(ld_md5), c(1, 2), mean)



#######################################Model comparison###########################################
##################################################################################################

ld_m4star_test = (ld_m4star[[1]]+ld_m4star[[2]]+ld_m4star[[3]])/3
loo_m4star = loo(as.matrix(ld_m4star_test))

ld_m1_test = (ld_m1[[1]]+ld_m1[[2]]+ld_m1[[3]])/3
ld_m2_test = (ld_m2[[1]]+ld_m2[[2]]+ld_m2[[3]])/3
ld_m3_test = (ld_m3[[1]]+ld_m3[[2]]+ld_m3[[3]])/3
ld_m4_test = (ld_m4[[1]]+ld_m4[[2]]+ld_m4[[3]])/3
ld_m4star_test = (ld_m4star[[1]]+ld_m4star[[2]]+ld_m4star[[3]])/3
ld_m5_test = (ld_m5[[1]]+ld_m5[[2]]+ld_m5[[3]])/3

loo_m1 = loo(as.matrix(ld_m1_test))
loo_m2 = loo(as.matrix(ld_m2_test))
loo_m3 = loo(as.matrix(ld_m3_test))
loo_m4 = loo(as.matrix(ld_m4_test))
loo_m4star = loo(as.matrix(ld_m4star_test))
loo_m5 = loo(as.matrix(ld_m5_test))

# fit(M4) > fit(M1)
loo_compare(loo_m4, loo_m1, loo_m5)

# fit(M4) > fit(M4*)
loo_compare(loo_m4,loo_m4star)

# fit(M4) > fit(M5)
loo_compare(loo_m4, loo_m5)

# fit(M5) > fit(M4*)
loo_compare(loo_m4star, loo_m5)


#######################################Model comparison DDM#######################################
##################################################################################################

ld_md1_test = (ld_md1[[1]]+ld_md1[[2]]+ld_md1[[3]])/3
ld_md2_test = (ld_md2[[1]]+ld_md2[[2]]+ld_md2[[3]])/3
ld_md3_test = (ld_md3[[1]]+ld_md3[[2]]+ld_md3[[3]])/3
ld_md4star_test = (ld_md4star[[1]]+ld_md4star[[2]]+ld_md4star[[3]])/3
ld_md4_test = (ld_md4[[1]]+ld_md4[[2]]+ld_md4[[3]])/3
ld_md5_test = (ld_md5[[1]]+ld_md5[[2]]+ld_md5[[3]])/3

loo_md1 = loo(as.matrix(ld_md1_test))
loo_md2 = loo(as.matrix(ld_md2_test))
loo_md3 = loo(as.matrix(ld_md3_test))
loo_md4 = loo(as.matrix(ld_md4_test))
loo_md4star = loo(as.matrix(ld_md4star_test))
loo_md5 = loo(as.matrix(ld_md5_test))

# fit(M4) > fit(M1)
loo_compare(loo_md4, loo_md1, loo_md5)

# fit(M4) > fit(M4*)
loo_compare(loo_md4,loo_md4star)

# fit(M4) > fit(M5)
loo_compare(loo_md4, loo_md5)

# fit(M5) > fit(M4*)
loo_compare(loo_md4star,loo_md5)

