# Setup
library(survey)
library(cobalt)
library(ranger)
library(partykit)
library(gbm)
library(mboost)

twoway_int = function(x, y){
  paste_fun = function(x, y) paste0("bols(", x, ", by = ", y, ", intercept = FALSE)")
  two_in = outer(x, y, FUN=paste_fun)
  paste0(t(two_in)[lower.tri(t(two_in))], collapse = "+")
}
#setwd("/Users/lingxiaowang/Google Drive/Machine learning methods for KW weights/Simulations")
#setwd("C:/Users/wangl29/Google Drive/Machine learning methods for KW weights/Simulations")
setwd("/home/wangl29/kw_ml")
# Load R functions for pseudo weights calculation
source("weighting_functions.R")
source("subfunctions.R")
seed1 = seednumber1[,1]
seed2 = seednumber1[,1]
#seeds = read.table("seed.txt", header = T)
#seed1 = seeds[, 1]
#seed2 = seeds[, 2]


# Population generation
N =100000
set.seed(9762342)
v_w = as.data.frame(matrix(rnorm(N*10), N, 10))
names(v_w) = c(paste0(rep("v", 8), c(1:6, 8, 9), sep = ""), "w7", "w10")
v_w$w5 = (v_w$v1*0.16+v_w$v5*0.84)*1.171
v_w$w6 = (v_w$v2*0.67+v_w$v6*0.33)*1.353
v_w$w8 = (v_w$v3*0.18+v_w$v8*0.82)*1.177
v_w$w9 = (v_w$v4*0.67+v_w$v9*0.33)*1.317
pop=v_w[,c(1:4, 11:14, 9, 10)]
pop = pop[,order(as.numeric(substr(names(pop), 2, 3)))]
names(pop) = paste0(rep("w", 10), c(1:10), sep = "")
#round(cor(pop), 1)
beta  = c(0, 1, 1, 1.5, 1.5, -.8, -.5, .7)
n.beta=length(beta)
alpha = c(-2.5, 1, 1, 1, 1, .71, -.19, .26)
n.alpha = length(alpha)
odds_y = exp(as.matrix(cbind(1, pop[,c(1:4, 8:10)]))%*%matrix(alpha, n.alpha, 1))
py = odds_y/(1+odds_y)
pop$y=as.numeric((runif(N)<py)); mean(pop$y)
pop$w2_2 = pop$w2^2
pop$w4_2 = pop$w4^2
pop$w7_2 = pop$w7^2
pop$w1_w3 = .5*pop$w1*pop$w3
pop$w2_w4 = .7*pop$w2*pop$w4
pop$w4_w5 = .5*pop$w4*pop$w5
pop$w5_w6 = .5*pop$w5*pop$w6
pop$w3_w5 = .5*pop$w3*pop$w5
pop$w4_w6 = .7*pop$w4*pop$w6
pop$w5_w7 = .5*pop$w5*pop$w7
pop$w1_w6 = .5*pop$w1*pop$w6
pop$w2_w3 = .7*pop$w2*pop$w3
pop$w3_w4 = .5*pop$w3*pop$w4
pop$w3_2_w5_2 = 0.6*pop$w3^2*pop$w5^2
pop$w1_w2_w3  = 0.6*pop$w1*pop$w2*pop$w3
pop$w4_w5_w7  = 0.6*pop$w4*pop$w5*pop$w7


pop$w1_c = cut(pop$w1, c(-5, -1, 0, 1, 6))
pop$w2_c = cut(pop$w2, c(-5, -1, 0, 1, 6))
pop$w3_c = cut(pop$w3, c(-5, -1, 0, 1, 6))
pop$w4_c = cut(pop$w4, c(-5, -1, 0, 1, 6))


pop$w1_s = 2*sqrt(abs(pop$w1+pop$w2+pop$w3))
pop$w2_s = 2*sqrt(abs(pop$w2+pop$w3+pop$w4))
pop$w3_s = 2*sqrt(abs(pop$w3+pop$w4+pop$w5))
pop$w4_s = 2*sqrt(abs(pop$w4+pop$w5+pop$w7))
#round(cor(pop[,c(1:24,32:35)]),2)

pop$w1_n = pop$w1+rnorm(N)/10+pop$w8/10
pop$w2_n = pop$w2+rnorm(N)+pop$w9/10
pop$w3_n = pop$w3+rnorm(N)+pop$w10/10
pop$w4_n = pop$w4+rnorm(N)+pop$y/2
#round(cor(pop[,c(1:24,36:39)]),2)


ds.m_9 = model.matrix(as.formula(paste("y~",paste(names(pop)[28:31],collapse="+"))), data = pop)
beta_9 = c(0, -0.8, -0.5, 0.5, 0.2, 2, 0.5, 3, 1, 0.5, -1, 2, 0.5)

covars   = names(pop)[1:7]
c_covars = names(pop)[1:7]
f_covars = NULL
if(sum(!(covars%in%c_covars))) f_covars=covars[!(covars%in%c_covars)]

odds=matrix(0, N, 10)
odds[,1]  = exp(as.matrix(cbind(1, pop[,c(1:7)]))%*%matrix(beta, n.beta, 1))
odds[,2]  = exp(as.matrix(cbind(1, pop[,c(1:7, 12)]))%*%matrix(beta[c(1:8, 3)], n.beta+1, 1))
odds[,3]  = exp(as.matrix(cbind(1, pop[,c(1:7, 12:14)]))%*%matrix(beta[c(1:8, 3, 5, 8)], n.beta+3, 1))
odds[,4]  = exp(as.matrix(cbind(1, pop[,c(1:7, 15:18)]))%*%matrix(beta[c(1:8, 2, 3, 5, 6)], n.beta+4, 1))
odds[,5]  = exp(as.matrix(cbind(1, pop[,c(1:7, 12, 15:18)]))%*%matrix(beta[c(1:8, 3, 2, 3, 5, 6)], n.beta+5, 1))
odds[,6]  = exp(as.matrix(cbind(1, pop[,c(1:7, 15, 16, 19:24, 17, 18)]))%*%matrix(beta[c(1:8, 2:6, 2:6)], n.beta+10, 1))
odds[,7]  = exp(as.matrix(cbind(1, pop[,c(1:7, 12:14, 15, 16, 19:24, 17, 18)]))%*%matrix(beta[c(1:8, 3, 5, 8, 2:6, 2:6)], n.beta+13, 1))
odds[,8]  = exp(as.matrix(cbind(1, pop[,c(1:7, 12:14, 15:18, 25:27)]))%*%matrix(beta[c(1:8, 3, 5, 8, 2, 3, 5, 6, 4, 2, 3)], n.beta+10, 1))
#odds[,9]  = exp(as.matrix(cbind(1, pop[,c(1:4, 32:35)]))%*%matrix(beta[c(1:5, 2:5)], 9, 1))
odds[,9] = exp(ds.m_9%*%matrix(beta_9, length(beta_9), 1))
odds[,10] = exp(as.matrix(cbind(1, pop[,c(36:39)]))%*%matrix(beta[1:5], 5, 1))

#q = as.data.frame(odds/(1+odds))
q_c = cbind(odds[,1]^0.3,   odds[,2]^0.25,  odds[,3]^0.25,  odds[,4]^0.27,  odds[,5]^0.25, 
            odds[,6]^0.22,  odds[,7]^0.17,  odds[,8]^0.18,  odds[,9]^0.4,   odds[,10]^0.23)
q_s = cbind(odds[,1]^-0.2,  odds[,2]^-0.2,  odds[,3]^-0.2,  odds[,4]^-0.17, odds[,5]^-0.17, 
            odds[,6]^-0.15, odds[,7]^-0.13, odds[,8]^-0.11, odds[,9]^-0.26, odds[,10]^-0.15)

Formulas = c("trt~w1+w2+w3+w4+w5+w6+w7",
             "trt~w1+w2+w3+w4+w5+w6+w7+w2_2",
             "trt~w1+w2+w3+w4+w5+w6+w7+w2_2+w4_2+w7_2",
             "trt~w1+w2+w3+w4+w5+w6+w7+w1_w3+w2_w4+w4_w5+w5_w6",
             "trt~w1+w2+w3+w4+w5+w6+w7+w2_2+w1_w3+w2_w4+w4_w5+w5_w6",
             "trt~w1+w2+w3+w4+w5+w6+w7+w1_w3+w2_w4+w4_w5+w5_w6+w3_w5+w4_w6+w5_w7+w1_w6+w2_w3+w3_w4",
             "trt~w1+w2+w3+w4+w5+w6+w7+w2_2+w4_2+w7_2+w1_w3+w2_w4+w4_w5+w5_w6+w3_w5+w4_w6+w5_w7+w1_w6+w2_w3+w3_w4",
             "trt~w1+w2+w3+w4+w5+w6+w7+w2_2+w4_2+w7_2+w1_w3+w2_w4+w4_w5+w5_w6+w3_w5+w4_w6+w5_w7+w1_w6+w2_w3+w3_w4+w3_2_w5_2+w1_w2_w3+w4_w5_w7",
             #"trt~w1+w2+w3+w4+w1_s+w2_s+w3_s+w4_s",
             "trt~w1_c+w2_c+w3_c+w4_c",
             "trt~w1_n+w2_n+w3_n+w4_n")

fm.p = update(as.formula(Formulas[1]), ~.^2)



NSIMU=5
n_c = 5000
n_s = 5000
est = array(0, c(NSIMU, 15, length(Formulas)),
            dimnames = list(c(1:NSIMU), c("cht", "cht_w", "svy_w", 
                                          "ipsw_t", "ipsw_m", "ipsw_i", 
                                          "kw_t", "kw_m", "kw_i", "kw_mob", 
                                          "kw_rf", "kw_crf", "kw_xtree", 
                                          "kw_gbm", "kw_mboost"),
                            paste0(rep("model", 10), c(1:10), sep=""))
)
wt_m = est
wt_v = est
smd_all = est

p_score.c=NULL
p_score.s=NULL
for (k in 1:10){
#k=1
  for(simu in 1:NSIMU){
#simu=1
    samp.c = samp.slct(seed = seed1[simu], 
                       fnt.pop = pop, 
                       n = n_c, 
                       dsgn = "pps", 
                       size = q_c[, k])
    wt_m[simu,2, k] = mean(samp.c$wt)
    wt_v[simu,2, k] = var(samp.c$wt)  
    sqrt(wt_v[simu,2, k])/wt_m[simu,2, k]
    est[simu, 1, k] = mean(samp.c$y)
    est[simu, 2, k] = sum(samp.c$y*samp.c$wt)/sum(samp.c$wt)
    samp.s = samp.slct(seed = seed1[simu], 
                       fnt.pop = pop, 
                       n = n_s, 
                       dsgn = "pps", 
                       size = q_s[, k])
    wt_m[simu,3, k] = mean(samp.s$wt)
    wt_v[simu,3, k] = var(samp.s$wt)
    sqrt(wt_v[simu,3, k])/wt_m[simu,3, k]
    est[simu, 3, k] = sum(samp.s$y*samp.s$wt)/sum(samp.s$wt)
     # Combine survey and cohort data 
    psa_dat = rbind(samp.c, samp.s)
    psa_dat$elig_wt = c(rep(1, n_c), samp.s$wt)
    psa_dat$trt_n =c(rep(1, n_c), rep(0, n_s))
    psa_dat$trt   = as.factor(psa_dat$trt_n)
    # Name of data source indicator in the combined sample 
    rsp_name="trt" # 1 for AARP, 0 for NHIS
    # Covariate balance before adjustment
    tab_pre_adjust <- bal.tab(psa_dat[, covars], treat = psa_dat$trt, weights = psa_dat$elig_wt, 
                              s.d.denom = "pooled", binary = "std", method="weighting")

    smd_all[simu,1, k] = mean(abs(tab_pre_adjust$Balance[, "Diff.Adj"]))

    smd_all[simu,2, k] = mean(abs(bal.tab(psa_dat[, covars], treat = psa_dat$trt, weights = psa_dat$wt, 
                                          s.d.denom = "pooled", binary = "std", method="weighting")$Balance[, "Diff.Adj"]))
    
    # IPSW method
    #True propensity score model
    ds = svydesign(ids=~1, weight = ~ elig_wt, data = psa_dat)
    lgtreg.w = svyglm(as.formula(Formulas[k]), family = binomial, design = ds)
    # Predict propensity scores
    p_score.w = lgtreg.w$fitted.values
    p_score.w.c = p_score.w[psa_dat[,rsp_name]==1]
    #Fitted propensity score model (main effects only)
    ds = svydesign(ids=~1, weight = ~ elig_wt, data = psa_dat)
    lgtreg.w = svyglm(as.formula(Formulas[1]), family = binomial, design = ds)
    # Predict propensity scores
    p_score.w = lgtreg.w$fitted.values
    p_score.w.c = cbind(p_score.w.c, p_score.w[psa_dat[,rsp_name]==1])
    
    #Fitted propensity score model (main effects+pairwise interactions)
    ds = svydesign(ids=~1, weight = ~ elig_wt, data = psa_dat)
    lgtreg.w = svyglm(fm.p, family = binomial, design = ds)
    # Predict propensity scores
    p_score.w = lgtreg.w$fitted.values
    p_score.w.c = cbind(p_score.w.c, p_score.w[psa_dat[,rsp_name]==1])
    
    
    n_pw=dim(p_score.w.c)[2]   # from weighted model, for IPSW
    ipsw = as.data.frame(matrix(0, n_c, n_pw))
    names(ipsw)=paste0(rep("ipsw", n_pw), c(1:n_pw), sep = "")
    for(i in 1:n_pw){
      # calculate KW weights
      ipsw[,i] = ipsw.wt(p_score.c = p_score.w.c[,i], svy.wt = samp.s$wt)
      psa_dat$wt_kw = psa_dat$elig_wt
      psa_dat$wt_kw[psa_dat$trt == 1] <- ipsw[,i]
      smd_all[simu,3+i, k] = mean(abs(bal.tab(psa_dat[, covars], treat = psa_dat$trt, weights = psa_dat$wt_kw, 
                                              s.d.denom = "pooled", binary = "std", method = "weighting")$Balance[, "Diff.Adj"]))
      est[simu,3+i, k]  = sum(samp.c$y*ipsw[,i])/sum(ipsw[,i])
      wt_m[simu,3+i, k] = mean(ipsw[,i])
      wt_v[simu,3+i, k] = var(ipsw[,i])
    }
    # Kernel weighting method
    kw = as.data.frame(matrix(0, n_c, 9))
    ###########################################################################
    ##                      True propensity score model                      ##                    
    ###########################################################################
    svyds = svydesign(ids =~1, weight = rep(1, n_c+n_s), data = psa_dat)
    lgtreg = svyglm(as.formula(Formulas[k]), family = binomial, design = svyds)
    p_score = lgtreg$fitted.values
    # Propensity scores for the cohort
    p_score.c = data.frame(lgst = p_score[psa_dat[,rsp_name]==1])
    # Propensity scores for the survey sample
    p_score.s = data.frame(lgst = p_score[psa_dat[,rsp_name]==0])
    kw[,1] = kw.wt(p_score.c = p_score.c[,1], p_score.s = p_score.s[,1], 
                   svy.wt= samp.s$wt, Large=F)$pswt
    psa_dat$wt_kw[psa_dat$trt == 1] <- kw[,1]
    smd_all[simu,7, k] = mean(abs(bal.tab(psa_dat[, covars], treat = psa_dat$trt, weights = psa_dat$wt_kw, 
                                            s.d.denom = "pooled", binary = "std", method = "weighting")$Balance[, "Diff.Adj"]))
    wt_m[simu,7, k] = mean(kw[,1])
    wt_v[simu,7, k] = var(kw[,1])                      
    est[simu, 7, k] = sum(samp.c$y*kw[,1])/sum(kw[,1])
    ###########################################################################
    ##                           Main effects only                           ##
    ###########################################################################
    svyds = svydesign(ids =~1, weight = rep(1, n_c+n_s), data = psa_dat)
    lgtreg = svyglm(as.formula(Formulas[1]), family = binomial, design = svyds)
    p_score = lgtreg$fitted.values
    # Propensity scores for the cohort
    p_score.c = cbind(p_score.c, lgst_m = p_score[psa_dat[,rsp_name]==1])
    # Propensity scores for the survey sample
    p_score.s = cbind(p_score.s, lgst_m = p_score[psa_dat[,rsp_name]==0])
    kw[,2] = kw.wt(p_score.c = p_score.c[,2], p_score.s = p_score.s[,2], 
                   svy.wt= samp.s$wt, Large=F)$pswt
    psa_dat$wt_kw[psa_dat$trt == 1] <- kw[,2]
    smd_all[simu,8, k] = mean(abs(bal.tab(psa_dat[, covars], treat = psa_dat$trt, weights = psa_dat$wt_kw, 
                                          s.d.denom = "pooled", binary = "std", method = "weighting")$Balance[, "Diff.Adj"]))
    wt_m[simu,8, k] = mean(kw[,2])
    wt_v[simu,8, k] = var(kw[,2])                      
    est[simu, 8, k] = sum(samp.c$y*kw[, 2])/sum(kw[, 2])
    ###########################################################################
    ##                 Main effects + Pairwise Interactions                  ##
    ###########################################################################
    svyds = svydesign(ids =~1, weight = rep(1, n_c+n_s), data = psa_dat)
    lgtreg = svyglm(fm.p, family = binomial, design = svyds)
    p_score = lgtreg$fitted.values
    # Propensity scores for the cohort
    p_score.c = cbind(p_score.c, lgst_m = p_score[psa_dat[,rsp_name]==1])
    # Propensity scores for the survey sample
    p_score.s = cbind(p_score.s, lgst_m = p_score[psa_dat[,rsp_name]==0])
    kw[,3] = kw.wt(p_score.c = p_score.c[,3], p_score.s = p_score.s[,3], 
                   svy.wt= samp.s$wt, Large=F)$pswt
    psa_dat$wt_kw[psa_dat$trt == 1] <- kw[,3]
    smd_all[simu,9, k] = mean(abs(bal.tab(psa_dat[, covars], treat = psa_dat$trt, weights = psa_dat$wt_kw, 
                                          s.d.denom = "pooled", binary = "std", method = "weighting")$Balance[, "Diff.Adj"]))
    wt_m[simu,9, k] = mean(kw[,3])
    wt_v[simu,9, k] = var(kw[,3])                      
    est[simu, 9, k] = sum(samp.c$y*kw[, 3])/sum(kw[, 3])

    cols <- ncol(model.matrix(trt~w1+w2+w3+w4+w5+w6+w7, data = psa_dat))
    ###########################################################################
    ##                Model-based recursive partitioning (MOB)               ##
    ###########################################################################
    # Set try-out values and prepare loop
    tune_maxdepth <- 2:10
    psa_dat$wt_kw <- psa_dat$elig_wt
    p_score       <- data.frame(matrix(ncol = length(tune_maxdepth), nrow = nrow(psa_dat)))
    p_score_c.tmp <- data.frame(matrix(ncol = length(tune_maxdepth), nrow = n_c))
    p_score_s.tmp <- data.frame(matrix(ncol = length(tune_maxdepth), nrow = n_s))
    smds <- rep(NA, length(tune_maxdepth)+1)
    smds[1] <- mean(abs(tab_pre_adjust$Balance[, "Diff.Adj"]))
    i <- 0
    # Loop over try-out values
    repeat {
      i <- i+1
      # Run model
      maxdepth <- tune_maxdepth[i]
      mob <- glmtree(trt~w1+w2+w3+w4+w5+w6+w7| w1+w2+w3+w4+w5+w6+w7, 
                     data = psa_dat,
                     family = binomial,
                     alpha = 0.05,
                     minsplit = NULL,
                     maxdepth = maxdepth)
      p_score[, i]       <- predict(mob, psa_dat, type = "response")
      p_score_c.tmp[, i] <- p_score[psa_dat$trt == 1, i]
      p_score_s.tmp[, i] <- p_score[psa_dat$trt == 0, i]
      # Calculate KW weights
      samp.c$kw <- kw.wt(p_score.c = p_score_c.tmp[,i], p_score.s = p_score_s.tmp[,i], 
                         svy.wt = samp.s$wt, Large=F)$pswt
      # Calculate covariate balance
      psa_dat$wt_kw[psa_dat$trt == 1] <- samp.c$kw
      smds[i+1] <- mean(abs(bal.tab(psa_dat[, covars], treat = psa_dat$trt, weights = psa_dat$wt_kw, 
                                    s.d.denom = "pooled", binary = "std", method = "weighting")$Balance[, "Diff.Adj"]))
      # Save KW weights of current iteration 
      names(samp.c)[dim(samp.c)[2]] <- paste0("kw.mob.", i)
      # Check improvement in covariate balance
      if (abs(smds[i] - smds[i+1]) < 0.001 | length(tune_maxdepth) == i){
        print(paste0("mob", i))
        break
      }
    }
    # Select KW weights with best average covariate balance
    best <- which.min(smds[2:(length(tune_maxdepth)+1)])
    p_score.c = cbind(p_score.c, mob = p_score_c.tmp[, best])
    p_score.s = cbind(p_score.s, mob = p_score_s.tmp[, best])
    kw[, 4] = samp.c[,names(samp.c) == paste0("kw.mob.", best)]
    smd_all[simu,10, k] = smds[1+best]
    wt_m[simu,10, k] = mean(kw[,4])
    wt_v[simu,10, k] = var(kw[,4])                      
    est[simu, 10, k] = sum(samp.c$y*kw[, 4])/sum(kw[, 4])
    samp.c[, grep("kw.mob", names(samp.c))] <- NULL
    ###########################################################################
    ##                           Random Forest (RF)                          ##
    ###########################################################################
    # Set try-out values and prepare loop
    tune_mtry <- c(floor(sqrt(cols)), floor(log(cols)))
    psa_dat$wt_kw <- psa_dat$elig_wt
    p_score       <- data.frame(matrix(ncol = length(tune_mtry), nrow = nrow(psa_dat)))
    p_score_c.tmp <- data.frame(matrix(ncol = length(tune_mtry), nrow = n_c))
    p_score_s.tmp <- data.frame(matrix(ncol = length(tune_mtry), nrow = n_s))
    smds <- rep(NA, length(tune_mtry))
    # Loop over try-out values, calculate KW weights and covariate balance
    for (i in seq_along(tune_mtry)) {
      mtry <- tune_mtry[i]
      rf <- ranger(trt~w1+w2+w3+w4+w5+w6+w7,
                   data = psa_dat,
                   splitrule = "gini",
                   num.trees = 500,
                   mtry = mtry,
                   min.node.size = 15,
                   probability = T)
      p_score[, i]       <- predict(rf, psa_dat, type = "response")$predictions[, 2]
      p_score_c.tmp[, i] <- p_score[psa_dat$trt == 1, i]
      p_score_s.tmp[, i] <- p_score[psa_dat$trt == 0, i]
      # Calculate KW weights
      samp.c$kw <- kw.wt(p_score.c = p_score_c.tmp[,i], p_score.s = p_score_s.tmp[,i], 
                         svy.wt = samp.s$wt, Large=F)$pswt
      # Calculate covariate balance
      psa_dat$wt_kw[psa_dat$trt == 1] <- samp.c$kw
      smds[i] <- mean(abs(bal.tab(psa_dat[, covars], treat = psa_dat$trt, weights = psa_dat$wt_kw, 
                                  s.d.denom = "pooled", binary = "std", method = "weighting")$Balance[, "Diff.Adj"]))
      # Save KW weights of current iteration
      names(samp.c)[dim(samp.c)[2]] <- paste0("kw.rf.", i)
      print(paste0("rf", i))
    }
    # Select KW weights with best average covariate balance
    best <- which.min(smds)
    p_score.c = cbind(p_score.c, rf = p_score_c.tmp[, best])
    p_score.s = cbind(p_score.s, rf = p_score_s.tmp[, best])	
    kw[, 5] = samp.c[,names(samp.c) == paste0("kw.rf.", best)]
    smd_all[simu,11, k] = smds[best]
    wt_m[simu,11, k] = mean(kw[,5])
    wt_v[simu,11, k] = var(kw[,5])                      
    est[simu, 11, k] = sum(samp.c$y*kw[, 5])/sum(kw[, 5])
    #names(samp.c)[names(samp.c) == paste0("kw.", "rf.", best, collapse = "")] <- "kw.5"
    samp.c[, grep("kw.rf", names(samp.c))] <- NULL
    ###########################################################################
    ##                     Conditional Random Forests (CRF)                  ##
    ###########################################################################
    # Set try-out values and prepare loop
    tune_mincriterion <- c(0.99, 0.95, 0.9)
    psa_dat$wt_kw <- psa_dat$elig_wt
    p_score  <- data.frame(matrix(ncol = length(tune_mincriterion), nrow = nrow(psa_dat)))
    p_score_c.tmp <- data.frame(matrix(ncol = length(tune_mincriterion), nrow = n_c))
    p_score_s.tmp <- data.frame(matrix(ncol = length(tune_mincriterion), nrow = n_s))
    smds <- rep(NA, length(tune_mincriterion))
    # Loop over try-out values
    for (i in seq_along(tune_mincriterion)){ 
      minc <- tune_mincriterion[i]
      crf <- cforest(trt~w1+w2+w3+w4+w5+w6+w7,
                     data = psa_dat,
                     control = ctree_control(mincriterion = minc),
                     ntree = 100)
      p_score[, i]       <- predict(crf, newdata = psa_dat, type = "prob")[, 2]
      p_score_c.tmp[, i] <- p_score[psa_dat$trt == 1, i]
      p_score_s.tmp[, i] <- p_score[psa_dat$trt == 0, i]
      # Calculate KW weights
      samp.c$kw <- kw.wt(p_score.c = p_score_c.tmp[,i], p_score.s = p_score_s.tmp[,i], 
                         svy.wt = samp.s$wt, Large=F)$pswt
      # Calculate covariate balance
      psa_dat$wt_kw[psa_dat$trt == 1] <- samp.c$kw
      smds[i] <- mean(abs(bal.tab(psa_dat[, covars], treat = psa_dat$trt, weights = psa_dat$wt_kw,
                                  s.d.denom = "pooled", binary = "std", method = "weighting")$Balance[, "Diff.Adj"]))
      # Save KW weights of current iteration 
      names(samp.c)[dim(samp.c)[2]] <- paste0("kw.crf.", i)
      print(paste0("crf", i))
    }
    # Select KW weights with best average covariate balance
    best <- which.min(smds)
    p_score.c = cbind(p_score.c, crf = p_score_c.tmp[, best])
    p_score.s = cbind(p_score.s, crf = p_score_s.tmp[, best])	
    kw[, 6] = samp.c[,names(samp.c) == paste0("kw.crf.", best)]
    smd_all[simu,12, k] = smds[best]
    wt_m[simu,12, k] = mean(kw[,6])
    wt_v[simu,12, k] = var(kw[,6])                      
    est[simu, 12, k] = sum(samp.c$y*kw[, 6])/sum(kw[, 6])
    #names(samp.c)[names(samp.c) == paste0("kw.", "crf.", best, collapse = "")] <- "kw.5"
    samp.c[, grep("kw.crf", names(samp.c))] <- NULL
    ###########################################################################
    ##                   Extremely Randomized Trees (XTREE)                  ##
    ###########################################################################
    # Set try-out values and prepare loop
    tune_mtry <- c(floor(sqrt(cols)), floor(log(cols)))
    psa_dat$wt_kw <- psa_dat$elig_wt
    p_score  <- data.frame(matrix(ncol = length(tune_mtry), nrow = nrow(psa_dat)))
    p_score_c.tmp <- data.frame(matrix(ncol = length(tune_mtry), nrow = n_c))
    p_score_s.tmp <- data.frame(matrix(ncol = length(tune_mtry), nrow = n_s))
    smds <- rep(NA, length(tune_mtry))
    # Loop over try-out values
    for (i in seq_along(tune_mtry)){ 
      mtry <- tune_mtry[i]
      xtree <- ranger(trt_n~w1+w2+w3+w4+w5+w6+w7,
                      data = psa_dat,
                      splitrule = "extratrees",
                      num.random.splits = 1,
                      num.trees = 500,
                      mtry = mtry,
                      min.node.size = 15,
                      probability = T)
      p_score[, i]       <- predict(xtree, psa_dat, type = "response")$predictions[, 2]
      p_score_c.tmp[, i] <- p_score[psa_dat$trt == 1, i]
      p_score_s.tmp[, i] <- p_score[psa_dat$trt == 0, i]
      # Calculate KW weights
      samp.c$kw <- kw.wt(p_score.c = p_score_c.tmp[,i], p_score.s = p_score_s.tmp[,i], 
                         svy.wt = samp.s$wt, Large=F)$pswt
      # Calculate covariate balance
      psa_dat$wt_kw[psa_dat$trt == 1] <- samp.c$kw
      smds[i] <- mean(abs(bal.tab(psa_dat[, covars], treat = psa_dat$trt, weights = psa_dat$wt_kw,
                                  s.d.denom = "pooled", binary = "std", method = "weighting")$Balance[, "Diff.Adj"]))
      # Save KW weights of current iteration 
      names(samp.c)[dim(samp.c)[2]] <- paste0("kw.xtree.", i)
      print(paste0("xtree", i))
    }
    # Select KW weights with best average covariate balance
    best <- which.min(smds)
    p_score.c = cbind(p_score.c, xtree = p_score_c.tmp[, best])
    p_score.s = cbind(p_score.s, xtree = p_score_s.tmp[, best])	
    kw[, 7] = samp.c[,names(samp.c) == paste0("kw.xtree.", best)]
    smd_all[simu,13, k] = smds[best]
    wt_m[simu,13, k] = mean(kw[,7])
    wt_v[simu,13, k] = var(kw[,7])                      
    est[simu, 13, k] = sum(samp.c$y*kw[, 7])/sum(kw[, 7])
    #names(samp.c)[names(samp.c) == paste0("kw.", "xtree.", best, collapse = "")] <- "kw.7"
    samp.c[, grep("kw.xtree", names(samp.c))] <- NULL
    ###########################################################################
    ##                        Gradient Boosting (GBM)                        ##
    ###########################################################################
    # Set try-out values and prepare loop
    #psa_dat$wt_kw <- psa_dat$elig_wt
    tune_idepth <- 1:5
    tune_ntree <- c(100, 250, 500, 1000, 2000)
    p_score_o_c <- data.frame(matrix(ncol = length(tune_idepth), nrow = n_c))
    p_score_o_s <- data.frame(matrix(ncol = length(tune_idepth), nrow = n_s))
    p_scores_i  <- data.frame(matrix(ncol = length(tune_ntree),  nrow = nrow(psa_dat)))
    p_score_i_s <- data.frame(matrix(ncol = length(tune_ntree),  nrow = n_c))
    p_score_i_c <- data.frame(matrix(ncol = length(tune_ntree),  nrow = n_s))
    smds_o <- rep(NA, length(tune_idepth))
    smds_i <- rep(NA, length(tune_ntree)+1)
    smds_i[1] <- mean(abs(tab_pre_adjust$Balance[, "Diff.Adj"]))
    # Outer loop over try-out values
    for (i in seq_along(tune_idepth)) {
      #print(i)
      idepth <- tune_idepth[i] 
      j <- 0
      # Inner loop over try-out values
      repeat {
        j <- j+1
        # Run model
        ntree <- tune_ntree[j]
        boost <- gbm(trt_n~w1+w2+w3+w4+w5+w6+w7, data = psa_dat,
                     distribution = "bernoulli",
                     n.trees = ntree,
                     interaction.depth = idepth,
                     shrinkage = 0.05,
                     bag.fraction = 1)
        p_scores_i[, j] <- predict(boost, psa_dat, n.trees = ntree, type = "response")
        p_score_i_c[, j] <- p_scores_i[psa_dat$trt == 1, j]
        p_score_i_s[, j] <- p_scores_i[psa_dat$trt == 0, j]
        # Calculate KW weights
        samp.c$kw <- kw.wt(p_score.c = p_score_i_c[, j], p_score.s = p_score_i_s[,j], 
                           svy.wt = samp.s$wt, Large=F)$pswt
        # Calculate covariate balance
        psa_dat$wt_kw[psa_dat$trt == 1] <- samp.c$kw
        smds_i[j+1] <- mean(abs(bal.tab(psa_dat[, covars], treat = psa_dat$trt, weights = psa_dat$wt_kw, s.d.denom = "pooled", 
                                        binary = "std", method = "weighting")$Balance[, "Diff.Adj"]))
        # Save KW weights of current iteration 
        names(samp.c)[dim(samp.c)[2]] <- paste0("kw.gbm.i", j)
        # Check improvement in covariate balance
        if (abs(smds_i[j] - smds_i[j+1]) < 0.001 | length(tune_ntree) == j){
          print(paste0("gbm", j))
          break
        }
      } 
      best <- which.min(smds_i[2:(length(tune_ntree)+1)])
      names(samp.c)[names(samp.c) == paste0("kw.gbm.i", best)] <- paste0("kw.gbm.o", i)
      samp.c[, grep("kw.gbm.i", names(samp.c))] <- NULL
      p_score_o_c[,i] <- p_score_i_c[, best]
      p_score_o_s[,i] <- p_score_i_s[, best]
      smds_o[i] <- min(smds_i[2:(length(tune_ntree)+1)], na.rm = T)
    }
    # Select KW weights with best average covariate balance
    best <- which.min(smds_o)
    p_score.c = cbind(p_score.c, gbm = p_score_o_c[, best])
    p_score.s = cbind(p_score.s, gbm = p_score_o_s[, best])	
    kw[, 8] = samp.c[,names(samp.c) == paste0("kw.gbm.o", best)]
    smd_all[simu,14, k] = smds_o[best]
    wt_m[simu,14, k] = mean(kw[,8])
    wt_v[simu,14, k] = var(kw[,8])                      
    est[simu, 14, k] = sum(samp.c$y*kw[, 8])/sum(kw[, 8])
    samp.c[, grep("kw.gbm.o", names(samp.c))] <- NULL
    ##########################################################################
    #                     Model-based Boosting (mboost)                     ##
    ##########################################################################
    #Set try-out values and prepare loop
    c_covars_c = paste0(c_covars, "_c")
    psa_dat[, c_covars_c] <- lapply(psa_dat[, c_covars],scale, scale = F)
    covars_rnm = c("int", c_covars_c, f_covars)
    psa_dat$wt_kw <- psa_dat$elig_wt
    psa_dat$int = rep(1, length(psa_dat$trt))
    tune_mstop <- c(50, 100, 250, 500)
    p_score       <- data.frame(matrix(ncol = length(tune_mstop), nrow = nrow(psa_dat)))
    p_score_c.tmp <- data.frame(matrix(ncol = length(tune_mstop), nrow = n_c))
    p_score_s.tmp <- data.frame(matrix(ncol = length(tune_mstop), nrow = n_s))
    smds <- rep(NA, length(tune_mstop)+1)
    smds[1] <- mean(abs(tab_pre_adjust$Balance[, "Diff.Adj"]))
    i <- 0
    # Loop over try-out values
    repeat {
      i <- i+1
      # Run model
      mstop <- tune_mstop[i]
      mboost <- gamboost(as.formula(paste("trt ~", paste0(" bols(", covars_rnm, ", intercept = FALSE)", collapse = "+"), "+",
                                          paste0(" bbs(", c_covars_c, ", center = TRUE, df = 1, knots = 20)", collapse = "+"), "+",
                                          twoway_int(covars_rnm[-1], covars_rnm[-1]))),
                         control = boost_control(mstop = mstop, 
                                                 nu = 0.1),
                         family = Binomial(link = "logit"),
                         data = psa_dat)
      p_score[,i]       <- as.numeric(predict(mboost, psa_dat, type = "response"))
      p_score_c.tmp[,i] <- p_score[psa_dat$trt == 1, i]
      p_score_s.tmp[,i] <- p_score[psa_dat$trt == 0, i]
      # Calculate KW weights
      samp.c$kw <- kw.wt(p_score.c = p_score_c.tmp[,i], p_score.s = p_score_s.tmp[,i], svy.wt = samp.s$wt, Large = F)$pswt
      # Calculate covariate balance
      psa_dat$wt_kw[psa_dat$trt == 1] <- samp.c$kw
      smds[i+1] <- mean(abs(bal.tab(psa_dat[, covars], treat = psa_dat$trt, weights = psa_dat$wt_kw, 
                                    s.d.denom = "pooled", binary = "std", method = "weighting")$Balance[, "Diff.Adj"]))
      # Save KW weights of current iteration 
      names(samp.c)[dim(samp.c)[2]] <- paste0("kw.mboost.", i)
      # Check improvement in covariate balance
      if (abs(smds[i] - smds[i+1]) < 0.001 | length(tune_mstop) == i){
        print(i)
        break
      }
    }
    # Select KW weights with best average covariate balance
    best <- which.min(smds[2:(length(tune_mstop)+1)])
    p_score.c = cbind(p_score.c, mboost = p_score_c.tmp[, best])
    p_score.s = cbind(p_score.s, mboost = p_score_s.tmp[, best])
    kw[, 9] = samp.c[,names(samp.c) == paste0("kw.mboost.", best)]
    smd_all[simu,15, k] = smds[best+1]
    wt_m[simu,15, k] = mean(kw[,9])
    wt_v[simu,15, k] = var(kw[,9])                      
    est[simu, 15, k] = sum(samp.c$y*kw[, 9])/sum(kw[, 9])
    samp.c[, grep("kw.mboost", names(samp.c))] <- NULL
    print(paste0("simu", simu))
  }
  print(paste0("model", k))
}

seed_k = cbind(node = kk, seed= seed1)

#for(k in 1:7) print(apply((sqrt(wt_v)/wt_m)[,,k], 2, mean))
#for(k in 1:7) print((apply(est[,,k], 2, mean)-mean(pop$y))/mean(pop$y)*100)
#for(k in 1:7) print(apply(est[,,k], 2, var)*1e4)

#save(wt_v, file = "wt_v.rda")
#save(wt_m, file = "wt_m.rda")
#save(est, file = "est.rda")

write.table(as.data.frame.table(wt_v), "/home/wangl29/kw_ml_r/wt_v.txt", append = T, row.names = F, col.names = F)
write.table(as.data.frame.table(wt_m), "/home/wangl29/kw_ml_r/wt_m.txt", append = T, row.names = F, col.names = F)
write.table(as.data.frame.table(est) , "/home/wangl29/kw_ml_r/est.txt" , append = T, row.names = F, col.names = F)
write.table(as.data.frame.table(smd_all) , "/home/wangl29/kw_ml_r/smd_all.txt" , append = T, row.names = F, col.names = F)

