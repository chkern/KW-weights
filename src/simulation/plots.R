# Setup
setwd("C:/Users/wangl29/Google Drive/Machine learning methods for KW weights/Simulations/simu0419")
library(tidyverse)
library(ggplot2)
est     = read.table("est.txt")
wt_m    = read.table("wt_m.txt")
wt_v    = read.table("wt_v.txt")
smd_all = read.table("smd_all.txt")
m.pop_y <- 0.18383

est$Wt_mean = wt_m[,4]
est$Wt_var = wt_v[,4]
est$CV = sqrt(wt_v[,4])/wt_m[,4]
est$Smd = smd_all[,4]
names(est)[1:4]=c("Simu", "Method", "Model", "Est")
est = est[!(est$Method%in%c("cht", "cht_w", "svy_w","ipsw_m", "ipsw_t", "kw_xtree")|est$Model%in%c("model2", "model4", "model5")),]

est$Method = factor(est$Method,
                    levels = c(#"cht", "cht_w", "ipsw_t", 
                               "kw_t", "kw_m", "kw_i",
                               "kw_mob", "kw_rf", "kw_crf", "kw_gbm", "kw_mboost", "ipsw_i"),
                    labels = c(#"Cht (naïve)", "Cht (wtd)", "IPSW (true)", 
                               "KW (true)", "KW (main)", "KW (pairwise)", "KW (MOB)", "KW (RF)", 
                               "KW (CRF)", "KW (GBM)","KW (MBoost)","IPSW (pairwise)")
                  )
est$Model = factor(est$Model,
                   levels = paste0("model", c(1, 3, 6:10)),
                   labels = c("M(0, 0, 0)", 
                              "M(0, 1, 0)",
                              "M(1, 0, 0)", 
                              "M(1, 1, 0)",
                              "M(2, 2, 0)", 
                              "M(0, 0, 1)",
                              "M(0, 0, 2)"))

est_stat = data.frame(bias = aggregate(est$Est, list(est$Model, est$Method), mean)[,3]-m.pop_y,
                      evar = aggregate(est$Est, list(est$Model, est$Method), var)[,3]*1e4)
est_stat$relb = est_stat$bias/m.pop_y*100
est_stat$mse.t = round(est_stat$bias^2*1e4+est_stat$evar,2)
est_stat$mse.p = est_stat$mse.t
trans_m = function(i, dat = est_stat) round(6+1.55/(41-6)*(dat$mse.t[i]-6), 2)
est_stat$mse.p[rev(order(est_stat$mse.t))[1:10]] = trans_m(rev(order(est_stat$mse.t))[1:10])
est_stat$Model = factor(rep(levels(est$Model),  9), levels(est$Model))
est_stat$Method = factor(rep(levels(est$Method), each = 7), levels(est$Method))

est_stat$Wt_mean = aggregate(est$Wt_mean, list(est$Model, est$Method), mean)[,3]
est_stat$Wt_var = aggregate(est$Wt_var, list(est$Model, est$Method), mean)[,3]
est_stat$CV.t = round(aggregate(est$CV, list(est$Model, est$Method), mean)[,3], 2)
est_stat$CV.p = est_stat$CV.t
trans_c = function(i, dat = est_stat) round(3+1/(9-3)*(dat$CV.t[i]-3), 2)
est_stat$CV.p[rev(order(est_stat$CV.t))[1:8]] = trans_c(rev(order(est_stat$CV.t))[1:8])

est_stat$Smd.t = round(aggregate(est$Smd, list(est$Model, est$Method), mean)[,3]*100, 2)
est_stat$Smd.p = est_stat$Smd.t
trans_s = function(i, dat = est_stat) round(4+2/(22-4)*(dat$Smd.t[i]-4), 2)
est_stat$Smd.p[rev(order(est_stat$Smd.t))[1:14]] = trans_s(rev(order(est_stat$Smd.t))[1:14])
est_stat_avg = data.frame(bias    = NA,
                          evar    = round(aggregate(est_stat$evar, list(est_stat$Method), mean)[,2], 2),
                          relb    = c(4.96,	7.36,	5.30,	5.13,	10.33,	4.57,	3.84,	3.69, 10.26),
                          mse.t   = round(aggregate(est_stat$mse.t, list(est_stat$Method), mean)[,2], 2),
                          mse.p   = NA, Model = rep("Average", 9), Method = unique(est_stat$Method), 
                          Wt_mean = round(aggregate(est_stat$Wt_mean, list(est_stat$Method), mean)[,2], 2), 
                          Wt_var  = round(aggregate(est_stat$Wt_var, list(est_stat$Method), mean)[,2], 2),
                          CV.t    = round(aggregate(est_stat$CV.t, list(est_stat$Method), mean)[,2], 2), 
                          CV.p    = NA,
                          Smd.t   = round(aggregate(est_stat$Smd.t, list(est_stat$Method), mean)[,2], 2),
                          Smd.p   = NA
                          )
est_stat_avg$mse.p = est_stat_avg$mse.t
est_stat_avg$mse.p[rev(order(est_stat_avg$mse.t))[1:2]] = trans_m(rev(order(est_stat_avg$mse.t))[1:2], dat = est_stat_avg)
est_stat_avg$CV.p = est_stat_avg$CV.t
est_stat_avg$CV.p[rev(order(est_stat_avg$CV.t))[1]] = trans_c(rev(order(est_stat_avg$CV.t))[1], dat = est_stat_avg)
est_stat_avg$Smd.p = est_stat_avg$Smd.t
est_stat_avg$Smd.p[rev(order(est_stat_avg$Smd.t))[1:2]] = trans_s(rev(order(est_stat_avg$Smd.t))[1:2], dat = est_stat_avg)
est_stat_a = rbind(est_stat, est_stat_avg)
est_stat_a$color = as.factor(ifelse(est_stat_a$Model=="Average", 1, 0))

table(est_stat_a$Model, useNA = "ifany")
# Estimate by method and model (boxplots over simu)

ggplot(aes(x = Model, y = Est), data = est) +
  geom_boxplot(outlier.size = 0.1) +
  stat_summary(fun.y=mean, geom="point", shape=20, size=2, color="blue", fill="blue")+
  facet_wrap(~ Method, ncol = 3) +
  geom_hline(yintercept = m.pop_y, color="red") +
  coord_flip(ylim = c(0.1, 0.25)) +
  labs(x = "", y = "") +
  scale_x_discrete(limits = rev(levels(est$Model))) +
  scale_fill_manual(guide="none")+
  theme_bw() +
  theme(axis.text.x = element_text(vjust=1),
        legend.position = "none")

ggsave("simu_p01.png", width = 11, height = 11)

# Squared error by model and method (mean and std over simu)
ggplot(est_stat_a) +
  geom_bar(stat="identity", aes(x = Model, y = mse.p, fill=color)) +
  geom_text(aes(x = Model, y = mse.p, label = mse.t), size=3, hjust = -0.02)+
  facet_wrap(~ Method, ncol = 3) +
  geom_hline(yintercept = 0) +
  coord_flip() +
  labs(x = "", y = "") +
  scale_x_discrete(limits = rev(levels(est_stat_a$Model))) +
  scale_fill_manual(values=c("grey40", "black"))+
  scale_y_continuous(limit = c(0, 8.5), 
                     breaks = c(seq(0, 4, 2), seq(6, 7.5, 0.5)),
                     labels=c(seq(0, 6, 2), rep("", 2), 40))+  
  theme_bw() +
  theme(axis.text.x = element_text(vjust=1),
        legend.position = "none")

ggsave("simu_p03.png", width = 11, height = 11)

# CV of weights by model and method (mean and std over simu)

ggplot(est_stat_a) +
  geom_bar(stat="identity", aes(x = Model, y = CV.p, fill=color)) +
  geom_text(aes(x = Model, y = CV.p, label = CV.t), size=3, hjust = -0.02)+
  facet_wrap(~ Method, ncol = 3) +
  geom_hline(yintercept = 0) +
  coord_flip() +
  labs(x = "", y = "") +
  scale_x_discrete(limits = rev(levels(est_stat_a$Model))) +
  scale_fill_manual(values=c("grey40", "black"))+
  scale_y_continuous(limit = c(0,4.5),
                     breaks = c(1:2, seq(3, 4, 0.16)),
                     labels=c(1:3, rep("", 5), 9))+  
  theme_bw() +
  theme(axis.text.x = element_text(vjust=1),
        legend.position = "none")


# smd
ggplot(est_stat_a) +
  geom_bar(stat="identity", aes(x = Model, y = Smd.p, fill=color)) +
  geom_text(aes(x = Model, y = Smd.p, label = Smd.t), size=3, hjust = -0.02)+
  facet_wrap(~ Method, ncol = 3) +
  geom_hline(yintercept = 0) +
  coord_flip() +
  labs(x = "", y = "") +
  scale_x_discrete(limits = rev(levels(est_stat_a$Model))) +
  scale_fill_manual(values=c("grey40", "black"))+
  scale_y_continuous(limit = c(0,6.7),
                     breaks = c(0, 2, seq(4, 6, 0.4)),
                     labels=c(seq(0, 4, 2), rep("", 4), 22))+  
  theme_bw() +
  theme(axis.text.x = element_text(vjust=1),
        legend.position = "none")

