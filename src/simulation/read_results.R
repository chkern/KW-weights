setwd("C:/Users/wangl29/Google Drive/Machine learning methods for KW weights/Simulations/simu0419")
est     = read.table("est.txt") 
wt_m    = read.table("wt_m.txt")
wt_v    = read.table("wt_v.txt")
smd_all = read.table("smd_all.txt")
m.pop_y = 0.18383
NSIMU=dim(est)[1]/10/15;NSIMU
mthd = c("cht", "cht_w", "svy_w", "ipsw_t", "ipsw_m", "ipsw_i", 
         "kw_t", "kw_m", "kw_i", "kw_mob", "kw_rf", "kw_crf", 
         "kw_xtree","kw_gbm", "kw_mboost")
model = paste0(rep("model", 10), c(1:10), sep="")
est = est[order(est$V3, est$V2, est$V1),]
est.array = array(est$V4, c(NSIMU, 15, 10))
dimnames(est.array) = list(c(1:NSIMU), levels(est$V2), levels(est$V3))
wt_v = wt_v[order(wt_v$V3, wt_v$V2, wt_v$V1),]
wt_v.array = array(wt_v$V4, c(NSIMU, 15, 10))
dimnames(wt_v.array) = list(c(1:NSIMU), levels(est$V2), levels(est$V3))
wt_m = wt_m[order(wt_m$V3, wt_m$V2, wt_m$V1),]
wt_m.array = array(wt_m$V4, c(NSIMU, 15, 10))
dimnames(wt_m.array) = list(c(1:NSIMU), levels(est$V2), levels(est$V3))
smd_all = smd_all[order(smd_all$V3, smd_all$V2, smd_all$V1),]
smd_all.array = array(smd_all$V4, c(NSIMU, 15, 10))
dimnames(smd_all.array) = list(c(1:NSIMU), levels(est$V2), levels(est$V3))

relb = empv = mse = cv = smd = NULL
for(k in 1:10) relb = rbind(relb,round((apply(est.array[,,k], 2, mean)-m.pop_y)/m.pop_y*100, 3))
relb = relb[match(model, levels(est$V3)),match(mthd, levels(est$V2))]
for(k in 1:10) empv = rbind(empv, round(apply(est.array[,,k], 2, var)*1e4, 3))
empv = empv[match(model, levels(est$V3)),match(mthd, levels(est$V2))]
for(k in 1:10) mse = rbind(mse, (round((apply(est.array[,,k], 2, mean)-m.pop_y)^2+apply(est.array[,,k], 2, var), 7)*1e4))
mse = mse[match(model, levels(est$V3)),match(mthd, levels(est$V2))]
for(k in 1:10) cv = rbind(cv, round(apply((sqrt(wt_v.array)/wt_m.array)[,,k], 2, mean), 3))
cv = cv[match(model, levels(est$V3)),match(mthd, levels(est$V2))]
for(k in 1:10) smd = rbind(smd, round(apply(smd_all.array[,,k], 2, mean), 3))
smd = smd[match(model, levels(est$V3)),match(mthd, levels(est$V2))]

relb;empv;mse;cv;smd

