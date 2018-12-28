# Setup
library(survey)
library(ranger)
library(partykit)
library(gbm)

setwd("/Users/lingxiaowang/Google Drive/Machine learning methods for KW weights/Simulations")
# Load R functions for pseudo weights calculation
source("weighting_functions1119.R")
source("subfunctions.R")
seeds = read.table("seed.txt", header = T)
seed1 = seeds[, 1]
seed2 = seeds[, 2]

# Population generation
 N =3000
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
beta  = c(0, .8, -.25, .6, -.4, -.8, -.5, .7)
n.beta=length(beta)
alpha = c(-2.5, .3, -.36, -.73, -.2, .71, -.19, .26)
n.alpha = length(alpha)
odds_y = exp(as.matrix(cbind(1, pop[,c(1:4, 8:10)]))%*%matrix(alpha, n.alpha, 1))
py = odds_y/(1+odds_y)
pop$y=as.numeric((runif(N)<py))
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

odds=matrix(0, N, 7)
odds[,1] = exp(as.matrix(cbind(1, pop[,c(1:7)]))%*%matrix(beta, n.beta, 1))
odds[,2] = exp(as.matrix(cbind(1, pop[,c(1:7, 12)]))%*%matrix(beta[c(1:8, 3)], n.beta+1, 1))
odds[,3] = exp(as.matrix(cbind(1, pop[,c(1:7, 12:14)]))%*%matrix(beta[c(1:8, 3, 5, 8)], n.beta+3, 1))
odds[,4] = exp(as.matrix(cbind(1, pop[,c(1:7, 15:18)]))%*%matrix(beta[c(1:8, 2, 3, 5, 6)], n.beta+4, 1))
odds[,5] = exp(as.matrix(cbind(1, pop[,c(1:7, 12, 15:18)]))%*%matrix(beta[c(1:8, 3, 2, 3, 5, 6)], n.beta+5, 1))
odds[,6] = exp(as.matrix(cbind(1, pop[,c(1:7, 15, 16, 19:24, 17, 18)]))%*%matrix(beta[c(1:8, 2:6, 2:6)], n.beta+10, 1))
odds[,7] = exp(as.matrix(cbind(1, pop[,c(1:7, 12:14, 15, 16, 19:24, 17, 18)]))%*%matrix(beta[c(1:8, 3, 5, 8, 2:6, 2:6)], n.beta+13, 1))
q= as.data.frame(odds/(1+odds))
names(q) = paste0(rep("q", 7), c(1:7), sep = "")
Formulas = c("trt~w1+w2+w3+w4+w5+w6+w7",
             "trt~w1+w2+w3+w4+w5+w6+w7+w2_2",
             "trt~w1+w2+w3+w4+w5+w6+w7+w2_2+w4_2+w7_2",
             "trt~w1+w2+w3+w4+w5+w6+w7+w1_w3+w2_w4+w4_w5+w5_w6",
             "trt~w1+w2+w3+w4+w5+w6+w7+w2_2+w1_w3+w2_w4+w4_w5+w5_w6",
             "trt~w1+w2+w3+w4+w5+w6+w7+w1_w3+w2_w4+w4_w5+w5_w6+w3_w5+w4_w6+w5_w7+w1_w6+w2_w3+w3_w4",
             "trt~w1+w2+w3+w4+w5+w6+w7+w2_2+w4_2+w7_2+w1_w3+w2_w4+w4_w5+w5_w6+w3_w5+w4_w6+w5_w7+w1_w6+w2_w3+w3_w4")
NSIMU=100
n_c = 150
n_s = 150
est = array(0, c(NSIMU, 11, 7),
            dimnames = list(c(1:NSIMU),
                            c("Naive", "wtd chrt", "wtd svy", "IPSW (true)", "IPSW(main)", "KW (true)", "KW(main)", "KW (MOB)", "KW (RF)", "KW (XTree)", "KW (GBM)"),
                            paste0(rep("model", 7), c(1:7), sep=""))
)
for (k in 1:7){
for(simu in 1:NSIMU){
	samp.c = samp.slct(seed = seed1[simu], 
	                   fnt.pop = pop, 
	                   n = n_c, 
	                   dsgn = "pps", 
	                   size = q[,paste0("q", k, collaps="")])
	est[simu, 1, k] = mean(samp.c$y)
	est[simu, 2, k] = sum(samp.c$y*samp.c$wt)/sum(samp.c$wt)
	samp.s = samp.slct(seed = seed1[simu], 
	                   fnt.pop = pop, 
	                   n = n_s, 
	                   dsgn = "pps", 
	                   size = 1/q[,paste0("q", k, collaps="")])
	est[simu, 3, k] = sum(samp.s$y*samp.s$wt)/sum(samp.s$wt)
	# Combine NHIS and AARP data 
	psa_dat = rbind(samp.c, samp.s)
	psa_dat$elig_wt = c(rep(1, n_c), samp.s$wt)
	psa_dat$trt_n =c(rep(1, n_c), rep(0, n_s))
    psa_dat$trt   = as.factor(psa_dat$trt_n)
	# Name of data source indicator in the combined sample 
	rsp_name="trt" # 1 for AARP, 0 for NHIS
	
	# IPSW method
	#True propensity score model
	ds = svydesign(ids=~1, weight = ~ elig_wt, data = psa_dat)
	lgtreg.w = svyglm(as.formula(Formulas[k]), family = binomial, design = ds)
	# Predict propensity scores
	p_score.w = predict.glm(lgtreg.w, type = "response")
	p_score.w.c = p_score.w[psa_dat[,rsp_name]==1]
	#Fitted propensity score model (main effects only)
	ds = svydesign(ids=~1, weight = ~ elig_wt, data = psa_dat)
	lgtreg.w = svyglm(as.formula(Formulas[1]), family = binomial, design = ds)
	# Predict propensity scores
	p_score.w = predict.glm(lgtreg.w, type = "response")
	p_score.w.c = cbind(p_score.w.c, p_score.w[psa_dat[,rsp_name]==1])

    n_pw=dim(p_score.w.c)[2]   # from unweighted model, for PSAS, and KW
    ipsw = as.data.frame(matrix(0, n_c, n_pw))
    names(ipsw)=paste0(rep("ipsw", n_pw), c(1:n_pw), sep = "")
    for(i in 1:n_pw){
    	# calculate KW weights
      	ipsw[,i] = ipsw.wt(p_score.c = p_score.w.c[,i], svy.wt = samp.s$wt)
      	est[simu,3+i, k] = sum(samp.c$y*ipsw[,i])/sum(ipsw[,i])
      	}

	# Kernel weighting method
	#True propensity score model
	svyds = svydesign(ids =~1, weight = rep(1, n_c+n_s), data = psa_dat)
	lgtreg = svyglm(as.formula(Formulas[k]), family = binomial, design = svyds)
	p_score = predict.glm(lgtreg, type = "response")
	# Propensity scores for the cohort
	p_score.c = p_score[psa_dat[,rsp_name]==1]
	# Propensity scores for the survey sample
	p_score.s = p_score[psa_dat[,rsp_name]==0]
	
	#Fitted propensity score model (main effects only)
	svyds = svydesign(ids =~1, weight = rep(1, n_c+n_s), data = psa_dat)
	lgtreg = svyglm(as.formula(Formulas[1]), family = binomial, design = svyds)
	p_score = predict.glm(lgtreg, type = "response")
	# Propensity scores for the cohort
	p_score.c = cbind(p_score.c, p_score[psa_dat[,rsp_name]==1])
	# Propensity scores for the survey sample
	p_score.s = cbind(p_score.s, p_score[psa_dat[,rsp_name]==0])
	
	cols <- ncol(model.matrix(trt~w1+w2+w3+w4+w5+w6+w7, data = psa_dat))
	# Model-based recursive partitioning (MOB)
	mob <- glmtree(trt~w1+w2+w3+w4+w5+w6+w7| w1+w2+w3+w4+w5+w6+w7, 
               data = psa_dat,
               family = binomial,
               alpha = 0.001,
               minsplit = 100,
               maxdepth = 5,
               #verbose = T,
               prune = "AIC")
    p_score_mob <- predict(mob, psa_dat, type = "response")
    # Propensity scores for the cohort
    p_score.c = cbind(p_score.c, p_score_mob[psa_dat[,rsp_name]==1])
    p_score.s = cbind(p_score.s, p_score_mob[psa_dat[,rsp_name]==0])
    # Random Forest (RF)
    rf <- ranger(trt~w1+w2+w3+w4+w5+w6+w7,
             data = psa_dat,
             splitrule = "gini",
             num.trees = 50,
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
    extr <- ranger(trt~w1+w2+w3+w4+w5+w6+w7,
               data = psa_dat,
               splitrule = "extratrees",
               num.random.splits = 1,
               num.trees = 50,
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
    boost <- gbm(trt_n~w1+w2+w3+w4+w5+w6+w7,
             data = psa_dat,
             distribution = "bernoulli",
             n.trees = 50,
             interaction.depth = 3,
             shrinkage = 0.05,
             bag.fraction = 0.5,
             n.cores = 4,
             #verbose = TRUE
             )
    p_score_boost <- predict(boost, psa_dat, n.trees = 50, type = "response")
    # Propensity scores for the cohort
    p_score.c = cbind(p_score.c, p_score_boost[psa_dat[,rsp_name]==1])
    p_score.s = cbind(p_score.s, p_score_boost[psa_dat[,rsp_name]==0])
    #cor(p_score.c)
    # Record the number of methods used for propensity score calculation
    n_p=dim(p_score.c)[2]   # from unweighted model, for PSAS, and KW
    kw = as.data.frame(matrix(0, n_c, n_p))
    names(kw)=paste0(rep("kw", n_p), c(1:n_p), sep = "")
    for(i in 1:n_p){
    	# calculate KW weights
      	kw[,i] = kw.wt(p_score.c = p_score.c[,i], p_score.s = p_score.s[,i], svy.wt= samp.s$wt, Large=F)$pswt
      	est[simu,5+i, k] = sum(samp.c$y*kw[,i])/sum(kw[,i])
      	}
    print(simu)
}
print(paste0("model", k))
}

for(k in 1:7) print((apply(est[,,k], 2, mean)-mean(pop$y))/mean(pop$y)*100)
