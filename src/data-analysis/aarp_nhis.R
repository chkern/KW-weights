# Setup
library(survey)
library(ranger)
library(partykit)
library(gbm)

# Set directory
setwd("C:/Users/wangl29/Box/Research/Lingxiao Projects/Machine learning methods/Nonprobability weighting")
# Read aarp data
aarp_syn = read.table("aarp_orig.txt", head=T)# change to "aarp_orig.txt" for the original aarp data
# Read nhis data
nhis_m = read.table("nhis.txt", head=T)
# Load R functions for pseudo weights calculation
source("weighting_functions1119.R")
# Check variable names in the two data sets. Please see "variable dictionary.xlsx" for data dictionary
names(aarp_syn)
names(nhis_m)
# Number of records in AARP and NHIS
n_c=dim(aarp_syn)[1]
n_s=dim(nhis_m)[1]

# NHIS estimate of 9-year all-cause mortality
est_nhis=sum(nhis_m$mtlty*nhis_m$elig_wt)/sum(nhis_m$elig_wt)
age_c.m_nhis = outer(nhis_m$age_c4, c(1:4), function(a, b)as.integer(a==b))
grp_wt.m_nihs = nhis_m$elig_wt*age_c.m_nhis
est_nhis = c(est_nhis, apply(nhis_m$mtlty*grp_wt.m_nihs, 2, sum)/apply(grp_wt.m_nihs, 2, sum))*100
round(est_nhis, 2)

# Combine NHIS and AARP data 
psa_dat = rbind(nhis_m, aarp_syn)
# Name of data source indicator in the combined sample 
rsp_name="trt" # 1 for AARP, 0 for NHIS

#Fitted propensity score model
Formula_fit = as.formula("trt ~ age+sex+as.factor(race)+as.factor(martl)+educ+bmi+as.factor(smk1)+phys+health+
                  sex:as.factor(martl)+educ:health+as.factor(race):educ+age:phys+age:as.factor(martl)+
                 phys:health+as.factor(race):health+age:as.factor(race)+sex:as.factor(smk1)+
                 as.factor(race):as.factor(smk1)+age:health+educ:phys+age:as.factor(smk1)+age:bmi+
                 educ:as.factor(smk1)+as.factor(martl):health+sex:educ+age:educ+as.factor(smk1):phys+
                 as.factor(race):as.factor(martl)+sex:health+sex:as.factor(race)+as.factor(race):phys+
                 bmi:as.factor(smk1)+as.factor(martl):phys+as.factor(smk1):health+bmi:phys+as.factor(martl):educ+
                 sex:bmi+educ:bmi+sex:phys")
# unweighted propenstiy score model
svyds = svydesign(ids =~1, weight = rep(1, n_c+n_s), data = psa_dat)
lgtreg = svyglm(Formula_fit, family = binomial, design = svyds)
p_score = predict.glm(lgtreg, type = "response")
# Propensity scores for the cohort
p_score.c = as.data.frame(p_score[psa_dat[,rsp_name]==1])
# Propensity scores for the survey sample
p_score.s = as.data.frame(p_score[psa_dat[,rsp_name]==0])
#### Fit logistic regression model to cohort and weighted survey sample
ds = svydesign(ids=~1, weight = ~ wt, data = psa_dat)
lgtreg.w = svyglm(Formula_fit, family = binomial, design = ds)
# Predict propensity scores
p_score.w = predict.glm(lgtreg.w, type = "response")
p_score.w.c = as.data.frame(p_score.w[psa_dat[,rsp_name]==1])
###########################################################################
# Add more propensity score calculation method here,                      #
# Save the propensity scores in data frames p_score.c and p_score.w.c     #
# for unweighted (PSAS and KW) and weighted (IPSW) analyses respectively  #
###########################################################################
# Data prep
psa_dat$trt_f <- as.factor(psa_dat$trt)
psa_dat$sex_f <- as.factor(psa_dat$sex)
psa_dat$race_f <- as.factor(psa_dat$race)
psa_dat$martl_f <- as.factor(psa_dat$martl)
psa_dat$smk1_f <- as.factor(psa_dat$smk1)
cols <- ncol(model.matrix(trt_f ~ age+sex_f+race_f+martl_f+educ+bmi+smk1_f+phys+health, data = psa_dat))
# Model-based recursive partitioning (MOB)
mob <- glmtree(trt_f ~ age+sex_f+race_f+martl_f+educ+bmi+smk1_f+phys+health | age+sex_f+race_f+martl_f+educ, 
               data = psa_dat,
               family = binomial,
               alpha = 0.001,
               minsplit = 1000,
               maxdepth = 5,
               verbose = T,
               prune = "AIC")
p_score_mob <- predict(mob, psa_dat, type = "response")
# Propensity scores for the cohort
p_score.c = cbind(p_score.c, p_score_mob[psa_dat[,rsp_name]==1])
p_score.s = cbind(p_score.s, p_score_mob[psa_dat[,rsp_name]==0])
# Random Forest (RF)
rf <- ranger(trt_f ~ age+sex_f+race_f+martl_f+educ+bmi+smk1_f+phys+health,
             data = psa_dat,
             splitrule = "gini",
             num.trees = 500,
             mtry = floor(sqrt(cols)),
             min.node.size = 15,
             write.forest = T,
             probability = T,
             save.memory = F)
p_score_rf <- predict(rf, psa_dat, type = "response")
# Propensity scores for the cohort
p_score.c = cbind(p_score.c, p_score_rf$predictions[, 2][psa_dat[,rsp_name]==1])
p_score.s = cbind(p_score.s, p_score_rf$predictions[, 2][psa_dat[,rsp_name]==0])

# Extremely Randomized Trees (XTREE)
extr <- ranger(trt_f ~ age+sex_f+race_f+martl_f+educ+bmi+smk1_f+phys+health,
               data = psa_dat,
               splitrule = "extratrees",
               num.random.splits = 1,
               num.trees = 500,
               mtry = floor(sqrt(cols)),
               min.node.size = 15,
               write.forest = T,
               probability = T,
               save.memory = F)
p_score_extr <- predict(extr, psa_dat, type = "response")
# Propensity scores for the cohort
p_score.c = cbind(p_score.c, p_score_extr$predictions[, 2][psa_dat[,rsp_name]==1])
p_score.s = cbind(p_score.s, p_score_extr$predictions[, 2][psa_dat[,rsp_name]==0])
# Gradient Boosting (GBM)
boost <- gbm(trt ~ age+sex_f+race_f+martl_f+educ+bmi+smk1_f+phys+health,
             data = psa_dat,
             distribution = "bernoulli",
             n.trees = 250,
             interaction.depth = 3,
             shrinkage = 0.05,
             bag.fraction = 0.5,
             n.cores = 4,
             verbose = TRUE)
p_score_boost <- predict(boost, psa_dat, n.trees = 250, type = "response")
# Propensity scores for the cohort
p_score.c = cbind(p_score.c, p_score_boost[psa_dat[,rsp_name]==1])
p_score.s = cbind(p_score.s, p_score_boost[psa_dat[,rsp_name]==0])

cor(p_score.c)

# Record the number of methods used for propensity score calculation
n_pw=dim(p_score.w.c)[2]  # from weighted model, for IPSW
n_p=dim(p_score.c)[2]   # from unweighted model, for PSAS, and KW

# Setup a matrix storing the weighted estimates
est=matrix(0, (n_pw+2*n_p), 5)
# Convert categorial age group variable to a matrix of dummy variables
age_c.m = outer(aarp_syn$age_c4, c(1:4), function(a, b)as.integer(a==b))

###########################################################################
#### Calcute pseudo weights and weighted estimates 
# Inverse of propensity score method
for (i in 1:n_pw){
  # calculate IPSW weights
  aarp_syn$ipsw=ipsw.wt(p_score.c = p_score.w.c[,i], svy.wt = nhis_m$elig_wt)
  # estimate overall estimate of mortality rate
  est[i,1] = sum(aarp_syn$mtlty*aarp_syn$ipsw)/sum(aarp_syn$ipsw)*100
  # estimate mortality rate by age group 
  grp_wt.m = aarp_syn$ipsw*age_c.m
  est[i,c(2:5)] = apply(aarp_syn$mtlty*grp_wt.m, 2, sum)/apply(grp_wt.m, 2, sum)*100
  # name the IPSW weights in the data frame of AARP
  names(aarp_syn)[dim(aarp_syn)[2]]=paste0("ipsw.", i, collapse = "")
}
# Sub-classification
for(i in 1:n_p){
  # calculate PSAS weights
  # Use propensity score estimated by different methods for cohort and survey sample 
  aarp_syn$psas = psas.wt(p_score.c = p_score.c[,i], p_score.s = p_score.s[,i], svy.wt= nhis_m$elig_wt, nclass = 5)$pswt
  # estimate overall estimate of mortality rate
  est[(i+n_pw),1] = sum(aarp_syn$mtlty*aarp_syn$psas)/sum(aarp_syn$psas)*100
  # estimate mortality rate by age group 
  grp_wt.m = aarp_syn$psas*age_c.m
  est[(i+n_pw),c(2:5)] = apply(aarp_syn$mtlty*grp_wt.m, 2, sum)/apply(grp_wt.m, 2, sum)*100
  # name the PSAS weights in the data frame of AARP
  names(aarp_syn)[dim(aarp_syn)[2]]=paste0("psas.", i, collapse = "")
}
  # Kernel-weighting 
for(i in 1:n_p){
  # calculate KW weights
  # Use propensity score estimated by different methods for cohort and survey sample 
  aarp_syn$kw = kw.wt(p_score.c = p_score.c[,i], p_score.s = p_score.s[,i], svy.wt= nhis_m$elig_wt, Large=T)$pswt
  est[(i+n_pw+n_p),1] = sum(aarp_syn$mtlty*aarp_syn$kw)/sum(aarp_syn$kw)*100
  grp_wt.m = aarp_syn$kw*age_c.m
  est[(i+n_pw+n_p),c(2:5)] = apply(aarp_syn$mtlty*grp_wt.m, 2, sum)/apply(grp_wt.m, 2, sum)*100
  names(aarp_syn)[dim(aarp_syn)[2]]=paste0("kw.", i, collapse = "")
  print(i)
}
  aarp_syn$kw = kw.wt(p_score.c = p_score.c_all[,1], 
                      p_score.s = p_score.s_all[,1], 
                      mtch_v=c(nhis_m$age_c4, aarp_syn$age_c4), 
                      svy.wt= nhis_m$elig_wt, Large=T)$pswt

# Naive cohort estimate 
est.cht=mean(aarp_syn$mtlty)*100
est.cht=c(est.cht, apply(aarp_syn$mtlty*age_c.m, 2, mean)/apply(age_c.m, 2, mean)*100)
est=rbind(est.cht, est)
# Relative difference from weighted NHIS estimate
rel.diff = t((t(est)-est_nhis)/est_nhis*100)
colnames(rel.diff)=c("Overall", "50-54", "55-59", "60-64", "64+")
# Please change the row names (weighting method)
rownames(rel.diff)= c("NIH-AARP", "IPSW", "PSAS-Logit", "PSAS-MOB", "PSAS-RF", "PSAS-XTREE", "PSAS-GBM", "KW-Logit", "KW-MOB", "KW-RF", "KW-XTREE", "KW-GBM")
round(rel.diff, 3)
# bias reduction%
bias.r = t((rel.diff[1,]-t(rel.diff))/rel.diff[1,])*100
colnames(bias.r)=c("Overall", "50-54", "55-59", "60-64", "64+")
# Please change the row names (weighting method)
rownames(bias.r)= c("NIH-AARP", "IPSW", "PSAS-Logit", "PSAS-MOB", "PSAS-RF", "PSAS-XTREE", "PSAS-GBM", "KW-Logit", "KW-MOB", "KW-RF", "KW-XTREE", "KW-GBM")
round(bias.r, 3) 
