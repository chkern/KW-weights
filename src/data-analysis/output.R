
output<-function(aarp_syn, nhis.m, var_name, byvar=F, d=F){
  #aarp_syn =aarp_wt; var_name=var.vec[i];byvar="age"; d=dvec[i]
  nhis_wt = (d==T)*nhis_m$elig_wt+(d==F)*nhis_m$wt
  est_nhis.ovl = sum(nhis_m[,grep(var_name,names(nhis_m))]*nhis_wt)/sum(nhis_wt)*100
  
  if(byvar!=F){
    v.byvar=nhis_m[,grep(byvar,names(nhis_m))]
    byvar.m_nhis = outer(v.byvar, c(1:nlevels(as.factor(v.byvar))), function(a, b) as.integer(a==b))
    grp_wt.m_nihs = nhis_wt*byvar.m_nhis
    est_nhis.byvar = apply(nhis_m[,grep(var_name,names(nhis_m))]*grp_wt.m_nihs, 2, sum)/apply(grp_wt.m_nihs, 2, sum)*100
    round(est_nhis.byvar, 2)
  }
  
  # estimate overall estimate of mortality rate
  est.ovl = rep(0,(n_pw+n_p1+n_p2))
  est.ovl[n_pw] = sum(aarp_syn[,grep(var_name,names(aarp_syn))]*aarp_syn$ipsw)/sum(aarp_syn$ipsw)*100
  est.ovl[n_pw+n_p1] = sum(aarp_syn[,grep(var_name,names(aarp_syn))]*aarp_syn$psas)/sum(aarp_syn$psas)*100
  
  if(byvar!=F){  # estimate mortality rate by age group 
    
    v.byvar=aarp_syn[,grep(byvar,names(aarp_syn))]
    ncat.byvar=nlevels(as.factor(v.byvar))
    # Setup a matrix storing the weighted estimates
    est.byvar = matrix(0, (n_pw+n_p1+n_p2), ncat.byvar)
    
    byvar.m = outer(v.byvar, c(1:ncat.byvar), function(a, b) as.integer(a==b))
    #### IPSW
    grp_wt.m = aarp_syn$ipsw*byvar.m
    est.byvar[n_pw,c(1:(ncat.byvar))] = apply(aarp_syn[,grep(var_name,names(aarp_syn))]*grp_wt.m, 2, sum)/apply(grp_wt.m, 2, sum)*100
    #### Sub-classification
    grp_wt.m = aarp_syn$psas*byvar.m
    est.byvar[n_pw+n_p1,c(1:(ncat.byvar))] = apply(aarp_syn[,grep(var_name,names(aarp_syn))]*grp_wt.m, 2, sum)/apply(grp_wt.m, 2, sum)*100
  }
  
  #### Kernel-weighting 
  for(i in 1:n_p2){
    # estimate overall estimate of mortality rate
    est.ovl[(i+n_pw+n_p1)] = sum(aarp_syn[,grep(var_name,names(aarp_syn))]*eval(parse(text = paste0("aarp_syn$kw.", i))))/sum(eval(parse(text = paste0("aarp_syn$kw.", i))))*100
    # estimate mortality rate by age group 
    if(byvar!=F){
      grp_wt.m = eval(parse(text = paste0("aarp_syn$kw.", i)))*byvar.m
      est.byvar[(i+n_pw+n_p1),c(1:(ncat.byvar))] = apply(aarp_syn[,grep(var_name,names(aarp_syn))]*grp_wt.m, 2, sum)/apply(grp_wt.m, 2, sum)*100
    }
  }
  
  # Naive cohort estimate 
  est.cht=mean(aarp_syn[,grep(var_name,names(aarp_syn))])*100
  est.ovl=c(est.cht,est.ovl)
  if(byvar!=F) {
    est.cht.byvar=apply(aarp_syn[,grep(var_name,names(aarp_syn))]*byvar.m, 2, mean)/apply(byvar.m, 2, mean)*100
    est.byvar=rbind(est.cht.byvar, est.byvar[1:(n_pw+n_p1+n_p2),])
  }
  
  
  # Relative difference from weighted NHIS estimate
  rel.diff = t((t(est.ovl)-est_nhis.ovl)/est_nhis.ovl*100)
  colnames(rel.diff)=c("Overall")
  if(byvar!=F) {rel.diff.byvar = t((t(est.byvar)-est_nhis.byvar)/est_nhis.byvar*100)
  rel.diff = cbind(rel.diff,rel.diff.byvar)
  }
  
  # Please change the row names (weighting method)
  # rownames(rel.diff)= c("NIH-AARP", "IPSW", "PSAS", "KW-Logit", "KW-MOB", "KW-RF", "KW-XTREE", "KW-GBM")
  rownames(rel.diff)= c("NIH-AARP", "IPSW", "PSAS", "KW-Logit", "KW-MOB", "KW-RF", "KW-XTREE",  
                      "KW-GBM","KW-mboost", "KW-CRF")

  # bias reduction%
  bias.r = t((rel.diff[1,]-t(rel.diff))/rel.diff[1,])*100
  
  # Please change the row names (weighting method)
  #rownames(bias.r)= c("NIH-AARP", "IPSW", "PSAS", "KW-Logit", "KW-MOB", "KW-RF", "KW-XTREE", "KW-GBM")
  rownames(bias.r)= c("NIH-AARP", "IPSW", "PSAS", "KW-Logit", "KW-MOB", "KW-RF", "KW-XTREE",  
                      "KW-GBM","KW-mboost", "KW-CRF")
  
  list(reldiff=round(rel.diff, 3),biasrdc=round(bias.r,3))
}



############################################################################################
# NHIS estimate of 9-year all-cause mortality
est_nhis = sum(nhis_m$mtlty*nhis_m$elig_wt)/sum(nhis_m$elig_wt)
age_c.m_nhis = outer(nhis_m$age_c4, c(1:4), function(a, b)as.integer(a==b))
grp_wt.m_nihs = nhis_m$elig_wt*age_c.m_nhis
est_nhis = c(est_nhis, apply(nhis_m$mtlty*grp_wt.m_nihs, 2, sum)/apply(grp_wt.m_nihs, 2, sum))*100
round(est_nhis, 2)

###########################################################################

#aarp_wt = read.table("aarp_orig_Oriweights.txt", head=TRUE)# change to "aarp_orig.txt" for the original aarp data
aarp_wt = read.table("aarp_wt.txt",head=TRUE)
# Record the number of methods used for propensity score calculation
n_pw = length(grep("ipsw", names(aarp_wt), value = T))  # from weighted model, for IPSW
n_p1 = length(grep("psas", names(aarp_wt), value = T)) # from unweighted models, for PSAS
n_p2 = length(grep("kw", names(aarp_wt), value = T)) # from unweighted models, for KW

f_covars=c("sex_f", "race_f", "martl_f", "smk1_f")
####compare balance
# Covariate balance after adjustment (KW)
tab_post_adjust_smd <- data.frame(matrix(ncol = n_p2, nrow = length(covars)))
for(i in 1:n_p2){
  psa_dat$wt_kw <- psa_dat$wt
  psa_dat$wt_kw[psa_dat$trt == 1] <- eval(parse(text = paste0("aarp_wt$kw.", i)))
  ds_kw <- svydesign(ids = ~1, weight = ~ wt_kw, data = psa_dat)
  tab_post_adjust <- svyCreateTableOne(vars = covars, factorVars = f_covars, strata = "trt", data = ds_kw, test = FALSE)
  tab_post_adjust_smd[,i] <- rbind(attr(tab_post_adjust$ContTable, "smd"), attr(tab_post_adjust$CatTable, "smd"))
print(i)
}
summary(tab_post_adjust_smd)
############################################################################
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

