# Setup
library(ggplot2)
library(gridExtra)
library(survey)
library(cobalt)
library(ranger)
library(partykit)
library(gbm)
library(mboost)

# Set directory
#setwd("C:/Users/wangl29/Box/Research/Lingxiao Projects/Machine learning methods/Nonprobability weighting")
#setwd("/Users/yanli/Box/Lingxiao Projects/Machine learning methods/Nonprobability weighting")
setwd("/home/wangl29/kw_ml")

# Read aarp data
aarp_syn = read.table("aarp_orig_all.txt", head=T)
# Read nhis data
nhis_m = read.table("nhis_all.txt", head=T)
# Load R functions for pseudo weights calculation
source("weighting_functions.R")
source("output.R")
# Check variable names in the two data sets. Please see "variable dictionary.xlsx" for data dictionary
names(aarp_syn)
names(nhis_m)
# Number of records in AARP and NHIS
n_c=dim(aarp_syn)[1]
n_s=dim(nhis_m)[1]

# Combine NHIS and AARP data 
psa_dat = rbind(nhis_m, aarp_syn)
#psa_dat$wt = c(nhis_m$elig_wt, rep(1, n_c))
psa_dat$wt = c(nhis_m$wt, rep(1, n_c))
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

#Formula_fit = as.formula("trt ~ age+sex+as.factor(race)+as.factor(martl)+educ+bmi+as.factor(smk1)+phys+health")
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

set.seed(0.3240913255)
###########################################################################
#### Calculate propensity scores and KW weights based on ML methods

####################################################################################
####                  Model-based recursive partitioning (MOB)                  ####
####################################################################################
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
  smds[i+1] <- mean(abs(bal.tab(psa_dat[, covars], treat = psa_dat$trt, weights = psa_dat$wt_kw, s.d.denom = "pooled", 
                                method = "weighting", binary = "std")$Balance[, "Diff.Adj"]))
  # Save KW weights of current iteration 
  names(aarp_syn)[dim(aarp_syn)[2]] <- paste0("kw.mob.", i)
  # Check improvement in covariate balance
  if (abs(smds[i] - smds[i+1]) < 0.001 | length(tune_maxdepth) == i){
    print(paste0("mob", i))
    break
  }
}

# Select KW weights with best average covariate balance
# Select KW weights with best average covariate balance
best <- which.min(smds[2:(length(tune_maxdepth)+1)])
print(best)
names(aarp_syn)[names(aarp_syn) == paste0("kw.mob.", best)] <- "kw.2"
aarp_syn[, grep("kw.mob", names(aarp_syn))] <- NULL
# Save propensity scores
psa_dat$ps.2 <- p_scores[, best]

####################################################################################
####                             Random Forest (RF)                             ####
####################################################################################
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
  smds[i] <- mean(abs(bal.tab(psa_dat[, covars], treat = psa_dat$trt, weights = psa_dat$wt_kw, s.d.denom = "pooled", 
                              method = "weighting", binary = "std")$Balance[, "Diff.Adj"]))
  # Save KW weights of current iteration
  names(aarp_syn)[dim(aarp_syn)[2]] <- paste0("kw.rf.", i)
  print(paste0("rf", i))
}

# Select KW weights with best average covariate balance
best <- which.min(smds)
names(aarp_syn)[names(aarp_syn) == paste0("kw.rf.", best)] <- "kw.3"
aarp_syn[, grep("kw.rf", names(aarp_syn))] <- NULL
# Save propensity scores
psa_dat$ps.3 <- p_scores[, best]



###########################################################################
##                   Extremely Randomized Trees (XTREE)                  ##
###########################################################################

# Set try-out values and prepare loop
psa_dat$wt_kw <- psa_dat$wt
tune_mtry <- c(floor(sqrt(cols)), floor(log(cols)))
p_scores <- data.frame(matrix(ncol = length(tune_mtry), nrow = nrow(psa_dat)))
smds <- rep(NA, length(tune_mtry))

# Loop over try-out values
for (i in seq_along(tune_mtry)) {
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
  smds[i] <- mean(abs(bal.tab(psa_dat[, covars], treat = psa_dat$trt, weights = psa_dat$wt_kw, s.d.denom = "pooled", 
                              method = "weighting", binary = "std")$Balance[, "Diff.Adj"]))
  # Save KW weights of current iteration 
  names(aarp_syn)[dim(aarp_syn)[2]] <- paste0("kw.xtree.", i)
  print(paste0("xtree", i))
}

# Select KW weights with best average covariate balance
best <- which.min(smds)
names(aarp_syn)[names(aarp_syn) == paste0("kw.xtree.", best)] <- "kw.4"
aarp_syn[, grep("kw.xtree", names(aarp_syn))] <- NULL
# Save propensity scores
psa_dat$ps.4 <- p_scores[, best]

####################################################################################
####                          Gradient Boosting (GBM)                           ####
####################################################################################
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
    smds_i[j+1] <- mean(abs(bal.tab(psa_dat[, covars], treat = psa_dat$trt, weights = psa_dat$wt_kw, s.d.denom = "pooled", 
                                    method = "weighting", binary = "std")$Balance[, "Diff.Adj"]))
    # Save KW weights of current iteration 
    names(aarp_syn)[dim(aarp_syn)[2]] <- paste0("kw.gbm.i", j)
    # Check improvement in covariate balance
    if (abs(smds_i[j] - smds_i[j+1]) < 0.001 | length(tune_ntree) == j){
      print(paste0("gbm", j))
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

####################################################################################
####                     Model-based Boosting (mboost)                          ####
####################################################################################
# Set try-out values and prepare loop
c_covars=c("age", "educ", "bmi", "phys", "health")
f_covars=c("sex_f", "race_f", "martl_f", "smk1_f")
c_covars_c = paste0(c_covars, "_c")
psa_dat[, c_covars_c] <- lapply(psa_dat[, c_covars],scale, scale = F)
covars_rnm = c("int", c_covars_c, f_covars)

twoway_int = function(x, y){
  paste_fun = function(x, y) paste0("bols(", x, ", by = ", y, ", intercept = FALSE)")
  two_in = outer(x, y, FUN=paste_fun)
  paste0(t(two_in)[lower.tri(t(two_in))], collapse = "+")
}

psa_dat$wt_kw <- psa_dat$wt
psa_dat$int = rep(1, length(psa_dat$trt))
tune_mstop <- c(50, 100, 250, 500)
p_scores <- data.frame(matrix(ncol = length(tune_mstop), nrow = nrow(psa_dat)))
smds <- rep(NA, length(tune_mstop)+1)
smds[1] <- mean(abs(tab_pre_adjust$Balance[, "Diff.Adj"]))
i <- 0

# Loop over try-out values
repeat {
  i <- i+1
  # Run model
  mstop <- tune_mstop[i]
  mboost <- gamboost(as.formula(paste("trt_f ~", paste0(" bols(", covars_rnm, ", intercept = FALSE)", collapse = "+"), "+",
                                      paste0(" bbs(", c_covars_c, ", center = TRUE, df = 1, knots = 20)", collapse = "+"), "+",
                                      twoway_int(covars_rnm[-1], covars_rnm[-1]))),
                     control = boost_control(mstop = mstop, 
                                             nu = 0.1),
                     family = Binomial(link = "logit"),
                     data = psa_dat)
  p_scores[, i] <- as.numeric(predict(mboost, data=psa_dat, type = "response"))
  p_score_c <- p_scores[psa_dat$trt == 1, i]
  p_score_s <- p_scores[psa_dat$trt == 0, i]
  # Calculate KW weights
  aarp_syn$kw <- kw.wt(p_score.c = p_score_c, p_score.s = p_score_s, svy.wt = nhis_m$elig_wt, Large = T)$pswt
  # Calculate covariate balance
  psa_dat$wt_kw[psa_dat$trt == 1] <- aarp_syn$kw
  smds[i+1] <- mean(abs(bal.tab(psa_dat[, covars], treat = psa_dat$trt, weights = psa_dat$wt_kw, 
                                s.d.denom = "pooled", binary = "std", method = "weighting")$Balance[, "Diff.Adj"]))
  # Save KW weights of current iteration 
  names(aarp_syn)[dim(aarp_syn)[2]] <- paste0("kw.mboost.", i)
  # Check improvement in covariate balance
  if (abs(smds[i] - smds[i+1]) < 0.001 | length(tune_mstop) == i){
    print(paste0("mboost", i))
    break
  }
}

# Select KW weights with best average covariate balance
best <- which.min(smds[2:(length(tune_mstop)+1)])
names(aarp_syn)[names(aarp_syn) == paste0("kw.mboost.", best)] <- "kw.6"
aarp_syn[, grep("kw.mboost", names(aarp_syn))] <- NULL
# Save propensity scores
psa_dat$ps.6 <- p_scores[, best]

###########################################################################
##                     Conditional Random Forests (CRF)                  ##
###########################################################################
# Set try-out values and prepare loop
tune_mincriterion <- c(0.99, 0.95, 0.9)
psa_dat$wt_kw <- psa_dat$elig_wt
p_scores <- NULL
smds <- rep(NA, length(tune_mincriterion))
for(k in 1:200){
  p_nrow = ifelse(k<200,2700,1714)
  p_scores.tmp  <- data.frame(matrix(ncol = length(tune_mincriterion), nrow = p_nrow))
  # Loop over try-out values
  for (i in seq_along(tune_mincriterion)){ 
    psa_dat$trt_f = as.factor(psa_dat$trt_f)
    minc <- tune_mincriterion[i]
    crf <- cforest(trt_f ~ age+sex_f+race_f+martl_f+educ+bmi+smk1_f+phys+health,
                   data = psa_dat,
                   control = ctree_control(mincriterion = minc),
                   ntree = 100)
    is.factor(psa_dat$trt_f)
    p_scores.tmp[, i] <- predict(crf, newdata = psa_dat[((k-1)*2700+1):min(k*2700, nrow(psa_dat)),], type = "prob")[, 2]
  }
  p_scores = rbind(p_scores, p_scores.tmp)
}

for(i in 1:3){
  p_score_c <- p_scores[psa_dat$trt == 1, i]
  p_score_s <- p_scores[psa_dat$trt == 0, i]
  # Calculate KW weights
  aarp_syn$kw <- kw.wt(p_score.c = p_score_c, p_score.s = p_score_s, svy.wt = nhis_m$elig_wt, Large = T)$pswt
  # Calculate covariate balance
  psa_dat$wt_kw[psa_dat$trt == 1] <- aarp_syn$kw
  smds[i] <- mean(abs(bal.tab(psa_dat[, covars], treat = psa_dat$trt, weights = psa_dat$wt_kw,
                              s.d.denom = "pooled", binary = "std", method = "weighting")$Balance[, "Diff.Adj"]))
  # Save KW weights of current iteration 
  names(aarp_syn)[dim(aarp_syn)[2]] <- paste0("kw.crf.", i)
  print(paste0("crf", i))
}
# Select KW weights with best average covariate bala
best <- which.min(smds)
names(aarp_syn)[names(aarp_syn) == paste0("kw.crf.", best)] <- "kw.7"
aarp_syn[, grep("kw.crf", names(aarp_syn))] <- NULL
# Save propensity scores
psa_dat$ps.7 <- p_scores[, best]






#################################################
################################################
# relative difference of Weighted Estimates of 9-year mortality
mtlty_age = output(aarp_syn =aarp_syn, nhis_m,mtlty,byvar="age_c4", d=T)$reldiff
mtlty_age = cbind(mtlty_age, apply(abs(mtlty_age)[,-1],1,mean))
colnames(mtlty_age) = c("Overall", "50-54", "55-59", "60-64", "64+", "Average")
round(mtlty_age, 2)





