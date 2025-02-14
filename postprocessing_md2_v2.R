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



################################################Model D2###########################################
##################################################################################################
#load("jags_sample_md2.RData")
#load("jags_sample_MD2_scale_dichotomized.RData")
load("jags_sample_MD2_task_v2.RData")
jagsModel_md2_v2 <- jagsModel_md2
#look into the traceplot of one variable
traceplot(jagsModel_md2_v2$mcmc[,"MuV[5]"])

mcmc_md2_v2 = as.matrix(as.mcmc.list(jagsModel_md2_v2$mcmc), chains = F)
quantiles_md2V2 <- apply(mcmc_md2_v2, 2, function(x) quantile(x, c(0.025, 0.975)))
posterior_mean_md2v2 =  colMeans(mcmc_md2_v2)


posterior_mean_md2v2["MuV[1]"]  # beta_L0
posterior_mean_md2v2["MuV[1]"] + posterior_mean_md2v2["MuV[2]"]  # beta_w0
posterior_mean_md2v2["MuV[3]"]  # beta_w1
posterior_mean_md2v2["MuV[4]"]  # beta_L1
posterior_mean_md2v2["MuV[5]"]  # beta_phi
posterior_mean_md2v2["MuZ[1]"]  # z0
posterior_mean_md2v2["MuZ[2]"]  # z2
posterior_mean_md2v2["MuZ[3]"]  # z3

quantiles_md2V2[, "MuV[1]"]  # beta_L0
quantiles_md2V2[, "MuV[1]"] + quantiles_md2V2[, "MuV[2]"]  # beta_w0
quantiles_md2V2[, "MuV[3]"]  # beta_w1
quantiles_md2V2[, "MuV[4]"]  # beta_L1
quantiles_md2V2[, "MuV[5]"]  # beta_phi
quantiles_md2V2[, "MuZ[1]"]  # z0
quantiles_md2V2[, "MuZ[2]"]  # z1
quantiles_md2V2[, "MuZ[3]"]  # z2

#group level parameters plot
par(mfrow = c(4, 3), cex.main = 3, cex.axis=1.5)
options(repr.plot.width = 20, repr.plot.height = 8)  
hist(mcmc_md2_v2[, "MuZ[1]"], breaks = 100, main = expression(mu[bold(z[0])]), xlab = expression(mu[bold(z[0])]),col="skyblue",border ="skyblue")
hist(mcmc_md2_v2[, "MuZ[2]"], breaks = 100, main = expression(mu[bold(z[1])]), xlab = expression(mu[bold(z[1])]),col="skyblue",border ="skyblue")
hist(mcmc_md2_v2[, "MuZ[3]"], breaks = 100, main = expression(mu[bold(z[phi])]), xlab = expression(mu[bold(z[phi])]),col="skyblue",border ="skyblue")
hist(mcmc_md2_v2[, "MuV[1]"], breaks = 100, main = expression(mu[bold(v[L0])]), xlab = expression(mu[bold(v[L0])]),col="skyblue",border ="skyblue")
hist(mcmc_md2_v2[, "MuV[1]"]+mcmc_md2_v2[, "MuV[2]"], breaks = 100, main = expression(mu[bold(v[w0])]), xlab = expression(mu[bold(v[w0])]),col="skyblue",border ="skyblue")
hist(mcmc_md2_v2[, "MuV[3]"], breaks = 100, main = expression(mu[bold(v[w1])]), xlab = expression(mu[bold(v[w1])]),col="skyblue",border ="skyblue")
hist(mcmc_md2_v2[, "MuV[4]"], breaks = 100, main = expression(mu[bold(v[L1])]), xlab = expression(mu[bold(v[L1])]),col="skyblue",border ="skyblue")
hist(mcmc_md2_v2[, "MuV[5]"], breaks = 100, main = expression(mu[bold(v[phi])]), xlab = expression(mu[bold(v[phi])]),col="skyblue",border ="skyblue")
hist(mcmc_md2_v2[, "mualpha"], breaks = 100, main = expression(mu[bold(alpha)]), xlab = expression(mu[bold(alpha)]),col="skyblue",border ="skyblue")
hist(mcmc_md2_v2[, "mutau"], breaks = 100, main = expression(mu[bold(tau)]), xlab = expression(mu[bold(tau)]),col="skyblue",border ="skyblue")


hist(abs(mcmc_md2_v2[, "MuV[1]"]+mcmc_md2_v2[, "MuV[2]"])-abs(mcmc_md2_v2[, "MuV[1]"]), breaks = 100, main = expression("|" *mu[bold(v[w0])]*"|"-"|" *mu[bold(v[L0])]*"|"), xlab = expression("|" *mu[bold(v[w0])]*"|"-"|" *mu[bold(v[L0])]*"|"), col="skyblue",border ="skyblue")
hist(abs(mcmc_md2_v2[, "MuV[3]"])-abs(mcmc_md2_v2[, "MuV[4]"]), breaks = 100, main = expression("|" *mu[bold(v[w1])]*"|"-"|" *mu[bold(v[L1])]*"|"), xlab = expression("|" *mu[bold(v[w1])]*"|"-"|" *mu[bold(v[L1])]*"|"), col="skyblue",border ="skyblue")

quantile(abs(mcmc_md2_v2[, "MuV[1]"]+mcmc_md2_v2[, "MuV[2]"])-abs(mcmc_md2_v2[, "MuV[1]"]), c(0.025, 0.975))
quantile(abs(mcmc_md2_v2[, "MuV[3]"])-abs(mcmc_md2_v2[, "MuV[4]"]), c(0.025, 0.975))


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
