# Setup
library(survey)
library(cobalt)
library(ranger)
library(partykit)
library(gbm)

setwd("/Users/lingxiaowang/Google Drive/Machine learning methods for KW weights/Simulations")
#setwd("C:/Users/wangl29/Google Drive/Machine learning methods for KW weights/Simulations")
#setwd("/home/wangl29/retrosp")
# Load R functions for pseudo weights calculation
source("weighting_functions.R")
source("subfunctions.R")
#seed1 = seednumber1[,1]
#seed2 = seednumber1[,1]
seeds = read.table("seed.txt", header = T)
seed1 = seeds[, 1]
seed2 = seeds[, 2]

# Population generation
N =50000
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
round(cor(pop), 1)
beta  = c(0, 1, 1, 1.5, 1.5, -.8, -.5, .7)
n.beta=length(beta)
alpha = c(-2.5, 1, 1, 1, 1, .71, -.19, .26)
n.alpha = length(alpha)
odds_y = exp(as.matrix(cbind(1, pop[,c(1:4, 8:10)]))%*%matrix(alpha, n.alpha, 1))
py = odds_y/(1+odds_y)
pop$y=as.numeric((runif(N)<py)); pop$y
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
covars = names(pop)[1:7]
odds=matrix(0, N, 7)
odds[,1] = exp(as.matrix(cbind(1, pop[,c(1:7)]))%*%matrix(beta, n.beta, 1))
odds[,2] = exp(as.matrix(cbind(1, pop[,c(1:7, 12)]))%*%matrix(beta[c(1:8, 3)], n.beta+1, 1))
odds[,3] = exp(as.matrix(cbind(1, pop[,c(1:7, 12:14)]))%*%matrix(beta[c(1:8, 3, 5, 8)], n.beta+3, 1))
odds[,4] = exp(as.matrix(cbind(1, pop[,c(1:7, 15:18)]))%*%matrix(beta[c(1:8, 2, 3, 5, 6)], n.beta+4, 1))
odds[,5] = exp(as.matrix(cbind(1, pop[,c(1:7, 12, 15:18)]))%*%matrix(beta[c(1:8, 3, 2, 3, 5, 6)], n.beta+5, 1))
odds[,6] = exp(as.matrix(cbind(1, pop[,c(1:7, 15, 16, 19:24, 17, 18)]))%*%matrix(beta[c(1:8, 2:6, 2:6)], n.beta+10, 1))
odds[,7] = exp(as.matrix(cbind(1, pop[,c(1:7, 12:14, 15, 16, 19:24, 17, 18)]))%*%matrix(beta[c(1:8, 3, 5, 8, 2:6, 2:6)], n.beta+13, 1))
#q = as.data.frame(odds/(1+odds))
q_c = cbind(odds[,1]^0.3, odds[,2]^0.25, odds[,3]^0.2, odds[,4]^0.27, odds[,5]^0.25, odds[,6]^0.22, odds[,7]^0.17)

q_s = cbind(odds[,1]^-0.2, odds[,2]^-0.2, odds[,3]^-0.15, odds[,4]^-0.17, odds[,5]^-0.17, odds[,6]^-0.15, odds[,7]^-0.13)

#names(q) = paste0(rep("q", 7), c(1:7), sep = "")
Formulas = c("trt~w1+w2+w3+w4+w5+w6+w7",
             "trt~w1+w2+w3+w4+w5+w6+w7+w2_2",
             "trt~w1+w2+w3+w4+w5+w6+w7+w2_2+w4_2+w7_2",
             "trt~w1+w2+w3+w4+w5+w6+w7+w1_w3+w2_w4+w4_w5+w5_w6",
             "trt~w1+w2+w3+w4+w5+w6+w7+w2_2+w1_w3+w2_w4+w4_w5+w5_w6",
             "trt~w1+w2+w3+w4+w5+w6+w7+w1_w3+w2_w4+w4_w5+w5_w6+w3_w5+w4_w6+w5_w7+w1_w6+w2_w3+w3_w4",
             "trt~w1+w2+w3+w4+w5+w6+w7+w2_2+w4_2+w7_2+w1_w3+w2_w4+w4_w5+w5_w6+w3_w5+w4_w6+w5_w7+w1_w6+w2_w3+w3_w4")
NSIMU=1000
n_c = 1000
n_s = 1000
est = array(0, c(NSIMU, 11, 7),
            dimnames = list(c(1:NSIMU), c("Naive", "wtd chrt", "wtd svy", "IPSW(true)", "IPSW(main)", "KW (true)", 
                                          "KW(main)", "KW (MOB)", "KW (RF)", "KW (XTree)", "KW (GBM)"),
                            paste0(rep("model", 7), c(1:7), sep=""))
)

wt_m = est
wt_v = est
for (k in 5:7){
  for(simu in 1:NSIMU){
    samp.c = samp.slct(seed = seed1[simu], 
	                     fnt.pop = pop, 
	                     n = n_c, 
	                     dsgn = "pps", 
	                     size = q_c[, k])
	  wt_m[simu,2, k] = mean(samp.c$wt)
	  wt_v[simu,2, k] = var(samp.c$wt)                  
	  est[simu, 1, k] = mean(samp.c$y)
	  est[simu, 2, k] = sum(samp.c$y*samp.c$wt)/sum(samp.c$wt)
	  samp.s = samp.slct(seed = seed1[simu], 
	                     fnt.pop = pop, 
	                     n = n_s, 
	                     dsgn = "pps", 
	                     size = q_s[, k])
	  wt_m[simu,3, k] = mean(samp.s$wt)
	  wt_v[simu,3, k] = var(samp.s$wt)
	  est[simu, 3, k] = sum(samp.s$y*samp.s$wt)/sum(samp.s$wt)
  
	  # Combine survey and cohort data 
	  psa_dat = rbind(samp.c, samp.s)
	  psa_dat$elig_wt = c(rep(1, n_c), samp.s$wt)
	  psa_dat$trt_n =c(rep(1, n_c), rep(0, n_s))
    psa_dat$trt   = as.factor(psa_dat$trt_n)
	  # Name of data source indicator in the combined sample 
	  rsp_name="trt" # 1 for AARP, 0 for NHIS
	  # Covariate balance before adjustment
	  tab_pre_adjust <- bal.tab(psa_dat[, covars], treat = psa_dat$trt, weights = psa_dat$wt, 
	                            s.d.denom = "pooled", binary = "std", method="weighting")
	  
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

    n_pw=dim(p_score.w.c)[2]   # from weighted model, for IPSW
    ipsw = as.data.frame(matrix(0, n_c, n_pw))
    names(ipsw)=paste0(rep("ipsw", n_pw), c(1:n_pw), sep = "")
    for(i in 1:n_pw){
      # calculate KW weights
      ipsw[,i] = ipsw.wt(p_score.c = p_score.w.c[,i], svy.wt = samp.s$wt)
      est[simu,3+i, k] = sum(samp.c$y*ipsw[,i])/sum(ipsw[,i])
      wt_m[simu,3+i, k] = mean(ipsw[,i])
	    wt_v[simu,3+i, k] = var(ipsw[,i])
    }

	  # Kernel weighting method
      kw = as.data.frame(matrix(0, n_c, 6))
	  #True propensity score model
	  svyds = svydesign(ids =~1, weight = rep(1, n_c+n_s), data = psa_dat)
	  lgtreg = svyglm(as.formula(Formulas[k]), family = binomial, design = svyds)
	  p_score = lgtreg$fitted.values
	  # Propensity scores for the cohort
	  p_score.c = p_score[psa_dat[,rsp_name]==1]
	  # Propensity scores for the survey sample
	  p_score.s = p_score[psa_dat[,rsp_name]==0]
      kw[,1] = kw.wt(p_score.c = p_score.c, p_score.s = p_score.s, 
                     svy.wt= samp.s$wt, Large=F)$pswt
	  wt_m[simu,6, k] = mean(kw[,1])
	  wt_v[simu,6, k] = var(kw[,1])                      
      est[simu, 6, k] = sum(samp.c$y*kw[,1])/sum(kw[,1])
	  
	  #Fitted propensity score model (main effects only)
	  svyds = svydesign(ids =~1, weight = rep(1, n_c+n_s), data = psa_dat)
	  lgtreg = svyglm(as.formula(Formulas[1]), family = binomial, design = svyds)
	  p_score = lgtreg$fitted.values
	  # Propensity scores for the cohort
	  p_score.c = cbind(p_score.c, p_score[psa_dat[,rsp_name]==1])
	  # Propensity scores for the survey sample
	  p_score.s = cbind(p_score.s, p_score[psa_dat[,rsp_name]==0])
  
      kw[,2] = kw.wt(p_score.c = p_score.c[,2], p_score.s = p_score.s[,2], 
                    svy.wt= samp.s$wt, Large=F)$pswt
	  wt_m[simu,7, k] = mean(kw[,2])
	  wt_v[simu,7, k] = var(kw[,2])                      
      est[simu, 7, k] = sum(samp.c$y*kw[, 2])/sum(kw[, 2])
	  
	  cols <- ncol(model.matrix(trt~w1+w2+w3+w4+w5+w6+w7, data = psa_dat))
	  
	  # Model-based recursive partitioning (MOB)
	  # Set try-out values and prepare loo
	  tune_maxdepth <- 2:10
	  psa_dat$wt_kw <- psa_dat$wt
	  p_scores.tmp  <- data.frame(matrix(ncol = length(tune_maxdepth), nrow = nrow(psa_dat)))
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
	    p_scores.tmp[, i]  <- predict(mob, psa_dat, type = "response")
	    p_score_c.tmp[, i] <- p_scores.tmp[psa_dat$trt == 1, i]
	    p_score_s.tmp[, i] <- p_scores.tmp[psa_dat$trt == 0, i]
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
	    if (smds[i] - smds[i+1] < 0.001 | length(tune_maxdepth) == i){
	      break
	    }
	  }
	  # Select KW weights with best average covariate balance
	  p_score.c = cbind(p_score.c, p_score_c.tmp[, ifelse(i == 1, i, i-1)])
	  p_score.s = cbind(p_score.s, p_score_s.tmp[, ifelse(i == 1, i, i-1)])
	  kw[, 3] = samp.c[,names(samp.c) == paste0("kw.mob.", ifelse(i == 1, i, i-1))]
	  wt_m[simu,8, k] = mean(kw[,3])
	  wt_v[simu,8, k] = var(kw[,3])                      
	  est[simu, 8, k] = sum(samp.c$y*kw[, 3])/sum(kw[, 3])
	  #names(samp.c)[names(samp.c) == paste0("kw.", "mob.", best, collapse = "")] <- "kw.3"
	  samp.c[, grep("kw.mob", names(samp.c))] <- NULL
	
	  #### Random Forest (RF)
	  # Set try-out values and prepare loop
	  tune_mtry <- c(floor(sqrt(cols)), floor(log(cols)))
	  psa_dat$wt_kw <- psa_dat$wt
	  p_scores.tmp <- data.frame(matrix(ncol = length(tune_mtry), nrow = nrow(psa_dat)))
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
      p_scores.tmp[, i]  <- predict(rf, psa_dat, type = "response")$predictions[, 2]
      p_score_c.tmp[, i] <- p_scores.tmp[psa_dat$trt == 1, i]
      p_score_s.tmp[, i] <- p_scores.tmp[psa_dat$trt == 0, i]
      # Calculate KW weights
      samp.c$kw <- kw.wt(p_score.c = p_score_c.tmp[,i], p_score.s = p_score_s.tmp[,i], 
                         svy.wt = samp.s$wt, Large=F)$pswt
      # Calculate covariate balance
      psa_dat$wt_kw[psa_dat$trt == 1] <- samp.c$kw
      smds[i] <- mean(abs(bal.tab(psa_dat[, covars], treat = psa_dat$trt, weights = psa_dat$wt_kw, 
                                  s.d.denom = "pooled", binary = "std", method = "weighting")$Balance[, "Diff.Adj"]))
      # Save KW weights of current iteration
      names(samp.c)[dim(samp.c)[2]] <- paste0("kw.rf.", i)
    }
    # Select KW weights with best average covariate balance
    best <- which.min(smds)
	  p_score.c = cbind(p_score.c, p_score_c.tmp[, best])
	  p_score.s = cbind(p_score.s, p_score_s.tmp[, best])	
	  kw[, 4] = samp.c[,names(samp.c) == paste0("kw.rf.", best)]
	  wt_m[simu,9, k] = mean(kw[,4])
	  wt_v[simu,9, k] = var(kw[,4])                      
	  est[simu, 9, k] = sum(samp.c$y*kw[, 4])/sum(kw[, 4])
    #names(samp.c)[names(samp.c) == paste0("kw.", "rf.", best, collapse = "")] <- "kw.4"
    samp.c[, grep("kw.rf", names(samp.c))] <- NULL
    
    #### Extremely Randomized Trees (XTREE)
    # Set try-out values and prepare loop
    #tune_mtry <- c(floor(sqrt(cols)), floor(log(cols)))
    #psa_dat$wt_kw <- psa_dat$wt
    #p_scores  <- data.frame(matrix(ncol = length(tune_mtry), nrow = nrow(psa_dat)))
    #p_score_c <- data.frame(matrix(ncol = length(tune_mtry), nrow = n_c))
    #p_score_s <- data.frame(matrix(ncol = length(tune_mtry), nrow = n_s))
    #smds <- rep(NA, length(tune_mtry))
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
      p_scores.tmp[, i]  <- predict(xtree, psa_dat, type = "response")$predictions[, 2]
      p_score_c.tmp[, i] <- p_scores.tmp[psa_dat$trt == 1, i]
      p_score_s.tmp[, i] <- p_scores.tmp[psa_dat$trt == 0, i]
      # Calculate KW weights
      samp.c$kw <- kw.wt(p_score.c = p_score_c.tmp[,i], p_score.s = p_score_s.tmp[,i], 
                         svy.wt = samp.s$wt, Large=F)$pswt
      # Calculate covariate balance
      psa_dat$wt_kw[psa_dat$trt == 1] <- samp.c$kw
      smds[i] <- mean(abs(bal.tab(psa_dat[, covars], treat = psa_dat$trt, weights = psa_dat$wt_kw,
                                  s.d.denom = "pooled", binary = "std", method = "weighting")$Balance[, "Diff.Adj"]))
      # Save KW weights of current iteration 
      names(samp.c)[dim(samp.c)[2]] <- paste0("kw.xtree.", i)
    }
    # Select KW weights with best average covariate balance
    best <- which.min(smds)
	  p_score.c = cbind(p_score.c, p_score_c.tmp[, best])
	  p_score.s = cbind(p_score.s, p_score_s.tmp[, best])	
	  kw[, 5] = samp.c[,names(samp.c) == paste0("kw.xtree.", best)]
	  wt_m[simu,10, k] = mean(kw[,5])
	  wt_v[simu,10, k] = var(kw[,5])                      
	  est[simu, 10, k] = sum(samp.c$y*kw[, 5])/sum(kw[, 5])
    #names(samp.c)[names(samp.c) == paste0("kw.", "xtree.", best, collapse = "")] <- "kw.5"
    samp.c[, grep("kw.xtree", names(samp.c))] <- NULL
    
    #### Gradient Boosting (GBM)
    # Set try-out values and prepare loop
    #psa_dat$wt_kw <- psa_dat$wt
    tune_idepth <- 1:5
    tune_ntree <- c(50, 100, 250, 500, 1000, 2000)
    p_scores_o  <- data.frame(matrix(ncol = length(tune_idepth), nrow = nrow(psa_dat)))
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
        if (smds_i[j] - smds_i[j+1] < 0.001 | length(tune_ntree) == j){
          break
        }
      } 
      names(samp.c)[names(samp.c) == paste0("kw.gbm.i", ifelse(j == 1, j, j-1))] <- paste0("kw.gbm.o", i)
      samp.c[, grep("kw.gbm.i", names(samp.c))] <- NULL
      p_scores_o[,i] <- p_scores_i[, ifelse(j == 1, j, j-1)]
      p_score_o_c[,i] <- p_score_i_c[, ifelse(j == 1, j, j-1)]
      p_score_o_s[,i] <- p_score_i_s[, ifelse(j == 1, j, j-1)]
      smds_o[i] <- smds_i[ifelse(j == 1, j+1, j)]
    }
    # Select KW weights with best average covariate balance
    best <- which.min(smds_o)
	  p_score.c = cbind(p_score.c, p_score_o_c[, best])
	  p_score.s = cbind(p_score.s, p_score_o_s[, best])	
	  kw[, 6] = samp.c[,names(samp.c) == paste0("kw.gbm.o", best)]
	  wt_m[simu,11, k] = mean(kw[,6])
	  wt_v[simu,11, k] = var(kw[,6])                      
	  est[simu, 11, k] = sum(samp.c$y*kw[, 6])/sum(kw[, 6])
    #names(samp.c)[names(samp.c) == paste0("kw.gbm.o", best)] <- "kw.6"
    samp.c[, grep("kw.gbm.o", names(samp.c))] <- NULL    
    print(simu)
  }
  print(paste0("model", k))
}


for(k in 1:7) print(apply((sqrt(wt_v)/wt_m)[,,k], 2, mean))
for(k in 1:7) print((apply(est[,,k], 2, mean)-mean(pop$y))/mean(pop$y)*100)
for(k in 1:7) print(apply(est[,,k], 2, var)*1e4)

save(wt_v, file = "wt_v.rda")
save(wt_m, file = "wt_m.rda")
save(est, file = "est.rda")

write.table(as.data.frame.table(wt_v), "wt_v.txt", append = T, row.names = F, col.names = F)
write.table(as.data.frame.table(wt_m), "wt_m.txt", append = T, row.names = F, col.names = F)
write.table(as.data.frame.table(est) , "est.txt" , append = T, row.names = F, col.names = F)

