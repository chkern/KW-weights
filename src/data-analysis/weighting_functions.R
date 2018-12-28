#################################################################################################
## ipsw.wt is a function calculating pseudo weights using IPSW methods                         ##
## INPUT:  p_score.c - predicted propensity score for cohort                                   ##
##         svy.wt    - a vector of survey weights                                              ##
## OUTPUT: pswt      - IPSW pseudo weights                                                     ##
#################################################################################################

ipsw.wt = function(p_score.c, svy.wt){
	  pswt = as.vector((1-p_score.c)/p_score.c)
	  pswt = pswt/sum(pswt)*sum(svy.wt)
	  pswt
	  }

#################################################################################################
## psas.wt is a function calculating pseudo weights using PSAS methods                         ##
## INPUT:  p_score.c - predicted propensity score for cohort                                   ##
##         p_score.s - predicted propensity score for survey                                   ##
##         svy.wt    - a vector of survey weights                                              ##
##         nclass    - number of subclasses (by percentiles) for sample division               ##
## OUTPUT: pswt      - PSAS pseudo weights                                                     ##
##         nclass    - actual number of subclasses used (empty classes will be combined)       ##
## WARNINGS:                                                                                   ##
##        if fewer than 2 cohort units in one or more subclasses                               ##
##        "Extreme weights may occur due to limited number(<=2) of cohort units in some cells" ##
##        if there are subclasses including no cohort units                                    ##
##        "Empty subclasses were combined with the neighbor subclass."                         ##
#################################################################################################
psas.wt = function(p_score.c, p_score.s, svy.wt, nclass){
	nclass0 = nclass
	m = length(p_score.c)
	n = length(p_score.s)
	p_score = c(p_score.c, p_score.s)
	trt = c(rep(1, m), rep(0, n))
	p_score.q = quantile(p_score, prob = seq(0, 1, length = (nclass+1)))
	p_score.q.u = unique(p_score.q)
	nclass = length(p_score.q.u)-1  
	subclass = cut(p_score, breaks = p_score.q.u, include.lowest = T)  
    levels(subclass) = c(1: nclass)
    nclass.c = length(unique(subclass[trt==1]))
    nclass.s = length(unique(subclass[trt==0]))
	  while (nclass.c!= nclass.s){
    nclass = min(nclass.c, nclass.s)
    p_score.q = quantile(p_score, prob = seq(0, 1, length = (nclass+1)))
    p_score.q.u = unique(p_score.q)
    nclass = length(p_score.q.u)-1  
    subclass = cut(p_score, breaks = p_score.q.u, include.lowest = T)  
    levels(subclass) = c(1: nclass)
    nclass.c = length(unique(subclass[trt==1]))
    nclass.s = length(unique(subclass[trt==0]))
  }
  # Assign pseudo weights to the cohort units
  svy_N = aggregate(svy.wt, by=list(subclass[trt==0]), FUN = sum)[,2]
  cht_n = aggregate(rep(1, m), by=list(subclass[trt==1]), FUN = sum)[,2]
  if(sum(cht_n<=2)>0) warning("Extreme weights may occur due to limited number(<=2) of cohort units in some cells")
  if(nclass<nclass0) warning("Empty cells were combined with the neighbor cells.")
  wt_f  = svy_N/cht_n
  pswt = rep(wt_f, cht_n)
  return(list(pswt = pswt, nclass = nclass))
}

###################################################################################################
## kw.wt is a function calculating pseudo weights using KW methods                               ##
## INPUT:  p_score.c - predicted propensity score for cohort                                     ##
##         p_score.s - predicted propensity score for survey                                     ##
##         svy.wt    - a vector of survey weights                                                ##
##         h         - bandwidth parameter                                                       ##
##                     (will be calculated corresponding to kernel function if not specified).   ##
##         krnfun    - kernel function                                                           ##
##                     "triang": triangular density on (-3, 3)                                     ##
##                     "dnorm":  standard normal density                                           ##
##                     "dnorm_t":truncated standard normal densigy on (-3, 3)                      ##
##         Large     - if the cohort size is so large that it has to be divided into pieces      ##         
##         rm.s      - removing unmatched survey units or not. Default is FALSE                  ##                                  
## OUTPUT: psd.wt    - KW pseudo weights                                                         ##
##         sum_0.s   - number of unmatched survey sample units                                   ##
## WARNINGS:                                                                                     ##
##         If there are unmatched survey sample units, the program gives                         ##
##         "The input bandwidth h is too small. Please choose a larger one!"                     ##
##           If rm.s=T, the program deletes unmatched survey sample units, and gives             ##
##           a warning "records in the prob sample were not used because of a small bandwidth"   ##
##           If rm.s=F, the program evenly distribute weights of unmatched survey sample units   ## 
##           to all cohot units.                                                                 ##
###################################################################################################
kw.wt = function(p_score.c, p_score.s, svy.wt, h=NULL, mtch_v = NULL, krn="triang", Large = F, rm.s = F){
  # get the name of kernel function
  # calculate bandwidth according to the kernel function
  #triangular density
  if(krn=="triang")h = bw.nrd0(p_score.c)/0.9*0.8586768
  if(krn=="dnorm"|krn=="dnorm_t")h = bw.nrd0(p_score.c)
  krnfun = get(krn)
  # create signed distance matrix    
  m = length(p_score.c)
  n = length(p_score.s)
    if (Large == F){
    sgn_dist_mtx = outer(p_score.s, p_score.c, FUN = "-")
    krn_num = krnfun(sgn_dist_mtx/h)
    if(is.null(mtch_v)){
      adj_m = 1
    }else{adj_m=outer(mtch_v[1:n], mtch_v[(n+1):(n+m)], FUN='==')}
    krn_num = krn_num*adj_m
    row.krn = rowSums(krn_num)
    sum_0.s = (row.krn==0)
    delt.svy = sum(sum_0.s)
    if(delt.svy>0){
      warning('The input bandwidth h is too small. Please choose a larger one!')
      if(rm.s == T){
        warning(paste(sum(sum_0.s), "records in the prob sample were not used because of a small bandwidth"))
        row.krn[sum_0.s]=1
      }else{
        krn_num[sum_0.s,]= 1
        row.krn[sum_0.s] = m
      } 
    }    
    row.krn = rowSums(krn_num)
    krn = krn_num/row.krn
    # QC: column sums should be 1 if rm.s=F; otherwise, some column sums could be 0
    #round(sum(rowSums(krn)),0) == round(dim(sgn_dist_mtx)[1],0) # TRUE
    #sum(rowSums(krn))== (dim(sgn_dist_mtx)[1]-sum(sum_0.s)) # TRUE
    # Final pseudo weights
    pswt_mtx = krn*svy.wt
    # QC: row sums should be weights for the survey sample
    # If rm.s = T, some of the survey sample units are not used.
    #svy.wt.u = svy.wt
    #svy.wt.u[sum_0.s]=0
    #sum(round(rowSums(pswt_mtx), 10) == round(svy.wt.u, 10))== dim(sgn_dist_mtx)[1]   #True
    psd.wt = colSums(pswt_mtx)
    # QC: sum of pseudo weights is equal to sum of survey sample weights
    #sum(psd.wt) == sum(svy.wt.u)       #True
  }else{
    psd.wt = rep(0, m)
    grp_size =  floor(n/50)
    up = c(seq(0, n, grp_size)[2:50], n)
    lw = seq(1, n, grp_size)[-51]
    delt.svy = 0
    for(g in 1:50){
      sgn_dist_mtx = outer(p_score.s[lw[g]:up[g]], p_score.c, FUN = "-")
      krn_num = krnfun(sgn_dist_mtx/h)
      if(is.null(mtch_v)){
        adj_m = 1
      }else{adj_m=outer(mtch_v[lw[g]:up[g]], mtch_v[(n+1):(n+m)], FUN='==')}
      krn_num = krn_num*adj_m
      row.krn = rowSums(krn_num)
      sum_0.s = (row.krn==0)
      delt.svy = delt.svy + (sum(sum_0.s)>0)
      if((sum(sum_0.s)>0)){
        warning('The input bandwidth h is too small. Please choose a larger one!')
        if(rm.s == T){
          warning(paste(sum(sum_0.s), "records in the prob sample were not used because of a small bandwidth"))
          row.krn[sum_0.s]=1
        }else{
          krn_num[sum_0.s,]= 1
          row.krn[sum_0.s] = m
        }
      }
      row.krn = rowSums(krn_num)
      krn = krn_num/row.krn
      # QC: column sums should be 1 if rm.s=F; otherwise, some column sums could be 0
      #round(sum(rowSums(krn)),0) == round(dim(sgn_dist_mtx)[1],0) # TRUE
      #sum(rowSums(krn))== (dim(sgn_dist_mtx)[1]-sum(sum_0.s)) # TRUE
      # Final psuedo weights
      pswt_mtx = krn*svy.wt[lw[g]:up[g]]
      # QC: row sums should be weights for the survey sample
      # If rm.s = T, some of the survey sample units are not used.
      #svy.wt.u = svy.wt
      #svy.wt.u[sum_0.s]=0
      #sum(round(rowSums(pswt_mtx), 10) == round(svy.wt.u, 10))== dim(sgn_dist_mtx)[1]   #True
      psd.wt = colSums(pswt_mtx) + psd.wt
    }
  }

  return(list(pswt = psd.wt, delt.svy = delt.svy, h = h))
} # end of kw.wt

#####################################################################################################################
## cmb_dat is a function for data preparation including removing missing values in cohort and survey sample,       ##
## combining the two samples                                                                                       ##
## INPUT:  chtsamp - cohort (a data frame including covariates of the propensity score model)                      ##
##         svysamp - survey sample (a data frame including weight variable and                                     ##
##                                  covariates of the propensity score model)                                      ##
##         svy_wt  - variable name of the survey weight (should be included in svysamp)                            ##
##         Formula - propensity score model, with the format of trt~covariate1+covariate2+...                      ##
## OUTPUT: cmb_dat - combined sample of non-missing cohort and survey sample units including varibles              ##
##                   covariates in propensity score model                                                          ##
##                   "trt" indicator of data source (1 for cohort units, 0 for survey sample units)                ##
##                   "wt" weight varaible (1 for cohort units, survey sample weights)                              ##
##         chtsamp - complete  sample (in terms of covariates)                                                     ##
##         svysamp - complete survey sample (in terms of covariates)                                               ##
## WARNINGS:                                                                                                       ##
##         1, missing values in covariates are not allowed. Records with missing values in the cohort are removed. ##
##         2, missing values in covariates are not allowed. Records with missing values in the survey sample       ## 
##         are removed. The complete cases are reweighted. Missing completely at random is assumed.                ##
#####################################################################################################################

cmb_dat = function(chtsamp,svysamp,svy_wt, Formula){
  # Get names of the response variable and predictors for the propensity score estimation model
  Fml_names = all.vars(Formula)
  # Name of the response variable
  rsp_name = Fml_names[1]
  # Name of the predictors
  mtch_var = Fml_names[-1]
  
  # Remove incomplete records in the cohort, if there are any.   
  chtsamp_sub = as.data.frame(chtsamp[, mtch_var])
  if(sum(is.na(chtsamp_sub))>0){
    cmplt.indx = complete.cases(chtsamp_sub)
    chtsamp_sub = chtsamp_sub[cmplt.indx, ]
    chtsamp = chtsamp[cmplt.indx,]
    warning("Missing values in covariates are not allowed. Records with missing values in the cohort are removed.")
  } 
  # Remove incomplete survey sample, if there are any.
  svysamp_sub = as.data.frame(svysamp[, mtch_var])
  svy_wt.vec = c(svysamp[, svy_wt])
  if(sum(is.na(svysamp_sub))>0){
    cmplt.indx = complete.cases(svysamp_sub)
    svysamp_sub = svysamp_sub[cmplt.indx, ]
    svy_wt.vec = sum(svysamp[, svy_wt])/sum(svy_wt.vec[cmplt.indx])*svy_wt.vec[cmplt.indx]
    warning("Missing values in covariates are not allowed. Records with missing values in the survey sample are removed. 
            The complete cases are reweighted. Missing completely at random is assumed.")
  }# end removing incomplete records
  # size of cohort (complete)
  m = dim(chtsamp_sub)[1]
  # size of survey sample (complete)
  n = dim(svysamp_sub)[1]
  # Combine the two complete samples
  # Set outcome variable for propensity score model 
  chtsamp_sub[,rsp_name] = 1                 # z=1 for cohort sample  
  svysamp_sub[,rsp_name] = 0                 # z=0 for survey sample
  names(chtsamp_sub) = c(mtch_var, rsp_name) # unify the variable names in cohort and survey sample
  names(svysamp_sub) = c(mtch_var, rsp_name)
  cmb_dat = rbind(chtsamp_sub, svysamp_sub)  # combine the two samples
  cmb_dat$wt = c(rep(1, m), svysamp[, svy_wt])
  return(list(cmb_dat = cmb_dat, 
              chtsamp = chtsamp, 
              svysamp = svysamp
              )
         )
}

#####################################################################################################################
## cmb_dat1 is a function combining the cohort and survey sample                                                   ##
## INPUT:  chtsamp - cohort (a data frame including covariates of the propensity score model)                      ##
##         svysamp - survey sample (a data frame including weight variable and                                     ##
##                                  covariates of the propensity score model)                                      ##
##         svy_wt  - variable name of the survey weight (should be included in svysamp)                            ##
##         Formula - propensity score model, with the format of trt~covariate1+covariate2+...                      ##
## OUTPUT: cmb_dat - combined sample of non-missing cohort and survey sample units including varibles              ##
##                   covariates in propensity score model                                                          ##
##                   "trt" indicator of data source (1 for cohort units, 0 for survey sample units)                ##
##                   "wt" weight varaible (1 for cohort units, survey sample weights)                              ##
#####################################################################################################################

cmb_dat1 = function(chtsamp,svysamp,svy_wt, Formula){
  # Get names of the response variable and predictors for the propensity score estimation model
  Fml_names = all.vars(Formula)
  # Name of the response variable
  rsp_name = Fml_names[1]
  # Name of the predictors
  mtch_var = Fml_names[-1]
  
  chtsamp_sub = as.data.frame(chtsamp[, mtch_var])
  svysamp_sub = as.data.frame(svysamp[, mtch_var])
  svy_wt.vec = c(svysamp[, svy_wt])
  m = dim(chtsamp_sub)[1]
  # size of survey sample (complete)
  n = dim(svysamp_sub)[1]
  # Combine the two complete samples
  # Set outcome variable for propensity score model 
  chtsamp_sub[,rsp_name] = 1                 # z=1 for cohort sample  
  svysamp_sub[,rsp_name] = 0                 # z=0 for survey sample
  names(chtsamp_sub) = c(mtch_var, rsp_name) # unify the variable names in cohort and survey sample
  names(svysamp_sub) = c(mtch_var, rsp_name)
  cmb_dat = rbind(chtsamp_sub, svysamp_sub)  # combine the two samples
  cmb_dat$wt = c(rep(1, m), svy_wt.vec)
  cmb_dat
  }

### Kernels ##
# triangular density on (-3, 3)
triang = function(x){
  x[abs(x)>3]=3
  1/3-abs(x)/3^2
}
# triangular densigh on (-2.5, 2.5)
triang_2 = function(x){
  x[abs(x)>2.5]=2.5
  1/2.5-abs(x)/2.5^2
}
# normal density (0, sigma=3)
dnorm_3 = function(x) dnorm(x, sd=3)
# truncated normal on (-3, 3)
dnorm_t = function(x){
  c = integrate(dnorm, -3, 3)$value
  y=dnorm(x)/c
  y[y<=dnorm(3)/c]=0
  y
}
