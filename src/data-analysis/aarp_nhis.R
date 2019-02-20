# Setup
library(ggplot2)
library(gridExtra)
library(survey)
library(cobalt)
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
source("weighting_functions.R")
# Check variable names in the two data sets. Please see "variable dictionary.xlsx" for data dictionary
names(aarp_syn)
names(nhis_m)
# Number of records in AARP and NHIS
n_c=dim(aarp_syn)[1]
n_s=dim(nhis_m)[1]

# Combine NHIS and AARP data 
psa_dat = rbind(nhis_m, aarp_syn)
psa_dat$wt = c(nhis_m$elig_wt, rep(1, n_c))
psa_dat$trt = c(rep(0, n_s), rep(1, n_c))
# Name of data source indicator in the combined sample 
rsp_name="trt" # 1 for AARP, 0 for NHIS

# Data prep
psa_dat$trt_f <- as.factor(psa_dat$trt)
psa_dat$sex_f <- as.factor(psa_dat$sex)
psa_dat$race_f <- as.factor(psa_dat$race)
psa_dat$martl_f <- as.factor(psa_dat$martl)
psa_dat$smk1_f <- as.factor(psa_dat$smk1)
cols <- ncol(model.matrix(trt_f ~ age+sex_f+race_f+martl_f+educ+bmi+smk1_f+phys+health, data = psa_dat))
covars <- c("age", "sex_f", "race_f", "martl_f", "educ", "bmi", "smk1_f", "phys", "health")

# Covariate balance before adjustment
tab_pre_adjust <- bal.tab(psa_dat[, covars], treat = psa_dat$trt, weights = psa_dat$wt, s.d.denom = "pooled", binary = "std")
summary(abs(tab_pre_adjust$Balance[, "Diff.Adj"]))

###########################################################################
#### Calculate propensity scores and pseudo weights based on logistic regression

#Fitted propensity score model
Formula_fit = as.formula("trt ~ age+sex+as.factor(race)+as.factor(martl)+educ+bmi+as.factor(smk1)+phys+health+
                  sex:as.factor(martl)+educ:health+as.factor(race):educ+age:phys+age:as.factor(martl)+
                 phys:health+as.factor(race):health+age:as.factor(race)+sex:as.factor(smk1)+
                 as.factor(race):as.factor(smk1)+age:health+educ:phys+age:as.factor(smk1)+age:bmi+
                 educ:as.factor(smk1)+as.factor(martl):health+sex:educ+age:educ+as.factor(smk1):phys+
                 as.factor(race):as.factor(martl)+sex:health+sex:as.factor(race)+as.factor(race):phys+
                 bmi:as.factor(smk1)+as.factor(martl):phys+as.factor(smk1):health+bmi:phys+as.factor(martl):educ+
                 sex:bmi+educ:bmi+sex:phys")
# unweighted propensity score model
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

# calculate IPSW weights
aarp_syn$ipsw = ipsw.wt(p_score.c = p_score.w.c[,1], svy.wt = nhis_m$elig_wt)
# calculate PSAS weights
aarp_syn$psas = psas.wt(p_score.c = p_score.c[,1], p_score.s = p_score.s[,1], svy.wt = nhis_m$elig_wt, nclass = 5)$pswt
# calculate KW weights
aarp_syn$kw.1 = kw.wt(p_score.c = p_score.c[,1], p_score.s = p_score.s[,1], svy.wt = nhis_m$elig_wt, Large=T)$pswt
# Save propensity scores
psa_dat$ps.1 = p_score

###########################################################################
#### Calculate propensity scores and KW weights based on ML methods

#### Model-based recursive partitioning (MOB)

# Set try-out values and prepare loop
psa_dat$wt_kw <- psa_dat$wt
tune_maxdepth <- 2:10
p_scores <- data.frame(matrix(ncol = length(tune_maxdepth), nrow = nrow(psa_dat)))
smds <- rep(NA, length(tune_maxdepth)+1)
smds[1] <- mean(abs(tab_pre_adjust$Balance[, "Diff.Adj"]))
i <- 0

# Loop over try-out values
repeat {
  i <- i+1
  print(i)
  # Run model
  maxdepth <- tune_maxdepth[i]
  mob <- glmtree(trt_f ~ age+sex_f+race_f+martl_f+educ+bmi+smk1_f+phys+health | age+sex_f+race_f+martl_f+educ+bmi+smk1_f+phys+health, 
                 data = psa_dat,
                 family = binomial,
                 alpha = 0.05,
                 minsplit = 1000,
                 maxdepth = maxdepth)
  p_scores[, i] <- predict(mob, psa_dat, type = "response")
  p_score_c <- p_scores[psa_dat$trt == 1, i]
  p_score_s <- p_scores[psa_dat$trt == 0, i]
  # Calculate KW weights
  aarp_syn$kw <- kw.wt(p_score.c = p_score_c, p_score.s = p_score_s, svy.wt = nhis_m$elig_wt, Large = T)$pswt
  # Calculate covariate balance
  psa_dat$wt_kw[psa_dat$trt == 1] <- aarp_syn$kw
  smds[i+1] <- mean(abs(bal.tab(psa_dat[, covars], treat = psa_dat$trt, weights = psa_dat$wt_kw, s.d.denom = "pooled", binary = "std")$Balance[, "Diff.Adj"]))
  # Save KW weights of current iteration 
  names(aarp_syn)[dim(aarp_syn)[2]] <- paste0("kw.mob.", i)
  # Check improvement in covariate balance
  if (smds[i] - smds[i+1] < 0.001 | length(tune_maxdepth) == i){
    break
  }
}

# Select KW weights with best average covariate balance
names(aarp_syn)[names(aarp_syn) == paste0("kw.mob.", ifelse(i == 1, i, i-1))] <- "kw.2"
aarp_syn[, grep("kw.mob", names(aarp_syn))] <- NULL
# Save propensity scores
psa_dat$ps.2 <- p_scores[, ifelse(i == 1, i, i-1)]

#### Random Forest (RF)

# Set try-out values and prepare loop
psa_dat$wt_kw <- psa_dat$wt
tune_mtry <- c(floor(sqrt(cols)), floor(log(cols)))
p_scores <- data.frame(matrix(ncol = length(tune_mtry), nrow = nrow(psa_dat)))
smds <- rep(NA, length(tune_mtry))

# Loop over try-out values
for (i in seq_along(tune_mtry)) {
  print(i)
  # Run model
  mtry <- tune_mtry[i]
  rf <- ranger(trt_f ~ age+sex_f+race_f+martl_f+educ+bmi+smk1_f+phys+health,
               data = psa_dat,
               splitrule = "gini",
               num.trees = 500,
               mtry = mtry,
               min.node.size = 15,
               probability = T)
  p_scores[, i] <- predict(rf, psa_dat, type = "response")$predictions[, 2]
  p_score_c <- p_scores[psa_dat$trt == 1, i]
  p_score_s <- p_scores[psa_dat$trt == 0, i]
  # Calculate KW weights
  aarp_syn$kw <- kw.wt(p_score.c = p_score_c, p_score.s = p_score_s, svy.wt = nhis_m$elig_wt, Large = T)$pswt
  # Calculate covariate balance
  psa_dat$wt_kw[psa_dat$trt == 1] <- aarp_syn$kw
  smds[i] <- mean(abs(bal.tab(psa_dat[, covars], treat = psa_dat$trt, weights = psa_dat$wt_kw, s.d.denom = "pooled", binary = "std")$Balance[, "Diff.Adj"]))
  # Save KW weights of current iteration
  names(aarp_syn)[dim(aarp_syn)[2]] <- paste0("kw.rf.", i)
}

# Select KW weights with best average covariate balance
best <- which.min(smds)
names(aarp_syn)[names(aarp_syn) == paste0("kw.rf.", best)] <- "kw.3"
aarp_syn[, grep("kw.rf", names(aarp_syn))] <- NULL
# Save propensity scores
psa_dat$ps.3 <- p_scores[, best]

#### Extremely Randomized Trees (XTREE)

# Set try-out values and prepare loop
psa_dat$wt_kw <- psa_dat$wt
tune_mtry <- c(floor(sqrt(cols)), floor(log(cols)))
p_scores <- data.frame(matrix(ncol = length(tune_mtry), nrow = nrow(psa_dat)))
smds <- rep(NA, length(tune_mtry))

# Loop over try-out values
for (i in seq_along(tune_mtry)) {
  print(i)
  # Run model
  mtry <- tune_mtry[i]
  xtree <- ranger(trt_f ~ age+sex_f+race_f+martl_f+educ+bmi+smk1_f+phys+health,
                  data = psa_dat,
                  splitrule = "extratrees",
                  num.random.splits = 1,
                  num.trees = 500,
                  mtry = mtry,
                  min.node.size = 15,
                  probability = T)
  p_scores[, i] <- predict(xtree, psa_dat, type = "response")$predictions[, 2]
  p_score_c <- p_scores[psa_dat$trt == 1, i]
  p_score_s <- p_scores[psa_dat$trt == 0, i]
  # Calculate KW weights
  aarp_syn$kw <- kw.wt(p_score.c = p_score_c, p_score.s = p_score_s, svy.wt = nhis_m$elig_wt, Large = T)$pswt
  # Calculate covariate balance
  psa_dat$wt_kw[psa_dat$trt == 1] <- aarp_syn$kw
  smds[i] <- mean(abs(bal.tab(psa_dat[, covars], treat = psa_dat$trt, weights = psa_dat$wt_kw, s.d.denom = "pooled", binary = "std")$Balance[, "Diff.Adj"]))
  # Save KW weights of current iteration 
  names(aarp_syn)[dim(aarp_syn)[2]] <- paste0("kw.xtree.", i)
}

# Select KW weights with best average covariate balance
best <- which.min(smds)
names(aarp_syn)[names(aarp_syn) == paste0("kw.xtree.", best)] <- "kw.4"
aarp_syn[, grep("kw.xtree", names(aarp_syn))] <- NULL
# Save propensity scores
psa_dat$ps.4 <- p_scores[, best]

#### Gradient Boosting (GBM)

# Set try-out values and prepare loop
psa_dat$wt_kw <- psa_dat$wt
tune_idepth <- 1:3
tune_ntree <- c(50, 100, 250, 500, 1000, 2000)
p_scores_o <- data.frame(matrix(ncol = length(tune_idepth), nrow = nrow(psa_dat)))
p_scores_i <- data.frame(matrix(ncol = length(tune_ntree), nrow = nrow(psa_dat)))
smds_o <- rep(NA, length(tune_idepth))
smds_i <- rep(NA, length(tune_ntree)+1)
smds_i[1] <- mean(abs(tab_pre_adjust$Balance[, "Diff.Adj"]))

# Outer loop over try-out values
for (i in seq_along(tune_idepth)) {
  print(i)
  idepth <- tune_idepth[i] 
  j <- 0
  # Inner loop over try-out values
  repeat {
    j <- j+1
    # Run model
    ntree <- tune_ntree[j]
    boost <- gbm(trt ~ age+sex_f+race_f+martl_f+educ+bmi+smk1_f+phys+health,
                 data = psa_dat,
                 distribution = "bernoulli",
                 n.trees = ntree,
                 interaction.depth = idepth,
                 shrinkage = 0.05,
                 bag.fraction = 1)
    p_scores_i[, j] <- predict(boost, psa_dat, n.trees = ntree, type = "response")
    p_score_c <- p_scores_i[psa_dat$trt == 1, j]
    p_score_s <- p_scores_i[psa_dat$trt == 0, j]
    # Calculate KW weights
    aarp_syn$kw <- kw.wt(p_score.c = p_score_c, p_score.s = p_score_s, svy.wt = nhis_m$elig_wt, Large = T)$pswt
    # Calculate covariate balance
    psa_dat$wt_kw[psa_dat$trt == 1] <- aarp_syn$kw
    smds_i[j+1] <- mean(abs(bal.tab(psa_dat[, covars], treat = psa_dat$trt, weights = psa_dat$wt_kw, s.d.denom = "pooled", binary = "std")$Balance[, "Diff.Adj"]))
    # Save KW weights of current iteration 
    names(aarp_syn)[dim(aarp_syn)[2]] <- paste0("kw.gbm.i", j)
    # Check improvement in covariate balance
    if (smds_i[j] - smds_i[j+1] < 0.001 | length(tune_ntree) == j){
      break
    }
  }
  names(aarp_syn)[names(aarp_syn) == paste0("kw.gbm.i", ifelse(j == 1, j, j-1))] <- paste0("kw.gbm.o", i)
  aarp_syn[, grep("kw.gbm.i", names(aarp_syn))] <- NULL
  p_scores_o[i] <- p_scores_i[, ifelse(j == 1, j, j-1)]
  smds_o[i] <- smds_i[ifelse(j == 1, j+1, j)]
}

# Select KW weights with best average covariate balance
best <- which.min(smds_o)
names(aarp_syn)[names(aarp_syn) == paste0("kw.gbm.o", best)] <- "kw.5"
aarp_syn[, grep("kw.gbm.o", names(aarp_syn))] <- NULL
# Save propensity scores
psa_dat$ps.5 <- p_scores_o[, best]

###########################################################################
#### Calcute weighted estimates 

# NHIS estimate of 9-year all-cause mortality
est_nhis = sum(nhis_m$mtlty*nhis_m$elig_wt)/sum(nhis_m$elig_wt)
age_c.m_nhis = outer(nhis_m$age_c4, c(1:4), function(a, b)as.integer(a==b))
grp_wt.m_nihs = nhis_m$elig_wt*age_c.m_nhis
est_nhis = c(est_nhis, apply(nhis_m$mtlty*grp_wt.m_nihs, 2, sum)/apply(grp_wt.m_nihs, 2, sum))*100
round(est_nhis, 2)

# Record the number of methods used for propensity score calculation
n_pw = length(grep("ipsw", names(aarp_syn), value = T))  # from weighted model, for IPSW
n_p1 = length(grep("psas", names(aarp_syn), value = T)) # from unweighted models, for PSAS
n_p2 = length(grep("kw", names(aarp_syn), value = T)) # from unweighted models, for KW

# Setup a matrix storing the weighted estimates
est = matrix(0, (n_pw+n_p1+n_p2), 5)
# Convert categorial age group variable to a matrix of dummy variables
age_c.m = outer(aarp_syn$age_c4, c(1:4), function(a, b)as.integer(a==b))

#### Inverse of propensity score method
# estimate overall estimate of mortality rate
est[n_pw,1] = sum(aarp_syn$mtlty*aarp_syn$ipsw)/sum(aarp_syn$ipsw)*100
# estimate mortality rate by age group 
grp_wt.m = aarp_syn$ipsw*age_c.m
est[n_pw,c(2:5)] = apply(aarp_syn$mtlty*grp_wt.m, 2, sum)/apply(grp_wt.m, 2, sum)*100

#### Sub-classification
# estimate overall estimate of mortality rate
est[n_pw+n_p1,1] = sum(aarp_syn$mtlty*aarp_syn$psas)/sum(aarp_syn$psas)*100
# estimate mortality rate by age group 
grp_wt.m = aarp_syn$psas*age_c.m
est[n_pw+n_p1,c(2:5)] = apply(aarp_syn$mtlty*grp_wt.m, 2, sum)/apply(grp_wt.m, 2, sum)*100

#### Kernel-weighting 
for(i in 1:n_p2){
  # estimate overall estimate of mortality rate
  est[(i+n_pw+n_p1),1] = sum(aarp_syn$mtlty*eval(parse(text = paste0("aarp_syn$kw.", i))))/sum(eval(parse(text = paste0("aarp_syn$kw.", i))))*100
  # estimate mortality rate by age group 
  grp_wt.m = eval(parse(text = paste0("aarp_syn$kw.", i)))*age_c.m
  est[(i+n_pw+n_p1),c(2:5)] = apply(aarp_syn$mtlty*grp_wt.m, 2, sum)/apply(grp_wt.m, 2, sum)*100
}

###########################################################################
#### Compare weighted estimates and balance

# Naive cohort estimate 
est.cht=mean(aarp_syn$mtlty)*100
est.cht=c(est.cht, apply(aarp_syn$mtlty*age_c.m, 2, mean)/apply(age_c.m, 2, mean)*100)
est=rbind(est.cht, est)
# Relative difference from weighted NHIS estimate
rel.diff = t((t(est)-est_nhis)/est_nhis*100)
colnames(rel.diff)=c("Overall", "50-54", "55-59", "60-64", "64+")
# Please change the row names (weighting method)
rownames(rel.diff)= c("NIH-AARP", "IPSW", "PSAS", "KW-Logit", "KW-MOB", "KW-RF", "KW-XTREE", "KW-GBM")
round(rel.diff, 3)
# bias reduction%
bias.r = t((rel.diff[1,]-t(rel.diff))/rel.diff[1,])*100
colnames(bias.r)=c("Overall", "50-54", "55-59", "60-64", "64+")
# Please change the row names (weighting method)
rownames(bias.r)= c("NIH-AARP", "IPSW", "PSAS", "KW-Logit", "KW-MOB", "KW-RF", "KW-XTREE", "KW-GBM")
round(bias.r, 3)

# Covariate balance before and after adjustment (KW)
psa_dat$educ_f <- as.factor(psa_dat$educ)
psa_dat$health_f <- as.factor(psa_dat$health)
covars2 <- c("age", "sex_f", "race_f", "martl_f", "educ_f", "bmi", "smk1_f", "phys", "health_f")
tab_post_adjust_smd <- data.frame(matrix(ncol = n_p2+1, nrow = 19))
aarp_syn$kw.0 <- aarp_syn$elig_wt
p_vars <- list()

for(i in 0:n_p2){
  psa_dat$wt_kw <- psa_dat$wt
  psa_dat$wt_kw[psa_dat$trt == 1] <- eval(parse(text = paste0("aarp_syn$kw.", i)))
  tab_post_adjust <- bal.tab(psa_dat[, covars], treat = psa_dat$trt, weights = psa_dat$wt_kw, s.d.denom = "pooled", binary = "std")
  tab_post_adjust_smd[, i+1] <- abs(tab_post_adjust$Balance[, "Diff.Adj"])
  for(var in covars2){
    p_vars[[paste0(i, var)]] <- bal.plot(psa_dat[, covars2], treat = psa_dat$trt, weights = psa_dat$wt_kw, var.name = var) +
      labs(title = "", subtitle = "")
  }
}

p_ps_u <- list()
p_ps_w <- list()

for(i in 1:n_p2){
  psa_dat$wt_kw <- psa_dat$wt
  p_ps_u[[i]] <- bal.plot(psa_dat, treat = psa_dat$trt, weights = psa_dat$wt_kw, var.name = paste0("ps.", i)) +
    labs(title = "", subtitle = "")
  psa_dat$wt_kw[psa_dat$trt == 1] <- eval(parse(text = paste0("aarp_syn$kw.", i)))
  p_ps_w[[i]] <- bal.plot(psa_dat, treat = psa_dat$trt, weights = psa_dat$wt_kw, var.name = paste0("ps.", i)) +
    labs(title = "", subtitle = "")
}

# Covariate balance plots (KW)
p0 <- grid.arrange(grobs = p_vars[c("0age", "0sex_f", "0race_f", "0martl_f", "0educ_f", "0bmi", "0smk1_f", "0phys", "0health_f")], top = "Distributional Balance: Unadjusted")
p1 <- grid.arrange(grobs = p_vars[c("1age", "1sex_f", "1race_f", "1martl_f", "1educ_f", "1bmi", "1smk1_f", "1phys", "1health_f")], top = "Distributional Balance: KW with logit regression")
p2 <- grid.arrange(grobs = p_vars[c("2age", "2sex_f", "2race_f", "2martl_f", "2educ_f", "2bmi", "2smk1_f", "2phys", "2health_f")], top = "Distributional Balance: KW with MOB")
p3 <- grid.arrange(grobs = p_vars[c("3age", "3sex_f", "3race_f", "3martl_f", "3educ_f", "3bmi", "3smk1_f", "3phys", "3health_f")], top = "Distributional Balance: KW with RF")
p4 <- grid.arrange(grobs = p_vars[c("4age", "4sex_f", "4race_f", "4martl_f", "4educ_f", "4bmi", "4smk1_f", "4phys", "4health_f")], top = "Distributional Balance: KW with XTREE")
p5 <- grid.arrange(grobs = p_vars[c("5age", "5sex_f", "5race_f", "5martl_f", "5educ_f", "5bmi", "5smk1_f", "5phys", "5health_f")], top = "Distributional Balance: KW with GBM")

# Propensity score plots (KW)
p6 <- grid.arrange(grobs = p_ps_u, top = "Distribution of propensity scores: Weighted NHIS, unweighted AARP")
p7 <- grid.arrange(grobs = p_ps_w, top = "Distribution of propensity scores: Weighted NHIS, KW-weighted AARP")
