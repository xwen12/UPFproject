# Aim 3: baseline levels of UPF-related lipids and incident CHD (Cox frailty models).
# Input: analysis datasets (1.3), 2.1 results. Output: result/ (used by 8.1 and 8.1.1).

rm(list = ls())

proj_dir <- "path/to/proj16_UPF_rework"
raw_dir  <- "path/to/SHS_raw_data"
library(dplyr)
rawData <- raw_dir
readyData <- file.path(proj_dir, "data")
result <-file.path(proj_dir, "result")
tabsFigs <- file.path(proj_dir, "result", "tabsFigs")
setwd(readyData)
load("lipidomics_p4.Rda")
load("lipidomics_p5.Rda")
load("dat_AimOne_p4.Rda")
setwd(result)
key_lipids_risk <- read.table('UPF_lipids_risk.csv',sep=',',header=T)
key_lipids_prot <- read.table('UPF_lipids_prot.csv',sep=',',header=T)

lipidomics_df <- as.data.frame(lipidomics_p4)
temp_1 <- lipidomics_df %>%
  inner_join(dat_AimOne_p4,by = c("ID" = "IDNO" ) )
dt <- temp_1 %>%
  mutate(nova4_kcal_P = scale(nova4_kcal_P,scale=TRUE,center = TRUE),
         S4G0 = scale(S4G0,scale=TRUE,center = TRUE),
         S4DT_KCAL = scale(S4DT_KCAL,scale=TRUE,center = TRUE),
         S4AGE = scale(S4AGE,scale=TRUE,center = TRUE),
         S4EDU = scale(S4EDU,scale=TRUE,center = TRUE),
         S4BMI = scale(S4BMI,scale=TRUE,center = TRUE),
         PDM3 = scale(PDM3,scale=TRUE,center = TRUE))

for(i in 2:1543){
  dt[,i]<-scale((dt[,i]),scale=TRUE,center=TRUE)
}

selected_lipids_names <- key_lipids_risk$lipids_name
n_lipids=length(selected_lipids_names)

lipids_select <- dt[, names(dt) %in% selected_lipids_names]

beta<-rep(0,n_lipids)
se<-rep(0,n_lipids)
HR<-rep(0,n_lipids)
upper_CI<-rep(0,n_lipids)
lower_CI<-rep(0,n_lipids)
pvals<-rep(0,n_lipids)
dt[,'FID']=as.factor(dt[,'FID'])

library(survival)
for(i in 1:n_lipids){
  lipid<-lipids_select[,i]
  fit <-coxph(Surv(TIME_chd,chd) ~ lipid+S4AGE+SEX+CENTER+S4EDU+S4BMI+S4SMOKE+S4ETOH+PDM3+LIPIDLOW+S4ADADM2020+S4USHTN+S4DT_KCAL+frailty(FID), data = dt)
  a=summary(fit)
  beta[i]      <- a$coefficients["lipid", "coef"]
  se[i]        <- a$coefficients["lipid", "se(coef)"]
  HR[i]        <- exp(beta[i])
  lower_CI[i]  <- exp(beta[i] - 1.96 * se[i])
  upper_CI[i]  <- exp(beta[i] + 1.96 * se[i])
  pvals[i]     <- a$coefficients["lipid", "p"]
}
qvals<-p.adjust(pvals,method="BH")
new2=data.frame(colnames(lipids_select),beta,se,HR,lower_CI,upper_CI,pvals,qvals)
colnames(new2)[1] <- "lipids_name"
write.table(new2,'longitu_chd_riskLipids.csv',sep=',',row.names = F, col.names=T)

selected_lipids_names <-key_lipids_prot$lipids_name
n_lipids=length(selected_lipids_names)

lipids_select <- dt[, names(dt) %in% selected_lipids_names]

beta<-rep(0,n_lipids)
se<-rep(0,n_lipids)
HR<-rep(0,n_lipids)
upper_CI<-rep(0,n_lipids)
lower_CI<-rep(0,n_lipids)
pvals<-rep(0,n_lipids)
dt[,'FID']=as.factor(dt[,'FID'])
library(lme4)

for(i in 1:n_lipids){
  lipid<-lipids_select[,i]
  fit <-coxph(Surv(TIME_chd,chd) ~ lipid+S4AGE+SEX+CENTER+S4EDU+S4BMI+S4SMOKE+S4ETOH+PDM3+LIPIDLOW+S4ADADM2020+S4USHTN+S4DT_KCAL+frailty(FID), data = dt)
  a=summary(fit)
  beta[i]      <- a$coefficients["lipid", "coef"]
  se[i]        <- a$coefficients["lipid", "se(coef)"]
  HR[i]        <- exp(beta[i])
  lower_CI[i]  <- exp(beta[i] - 1.96 * se[i])
  upper_CI[i]  <- exp(beta[i] + 1.96 * se[i])
  pvals[i]     <- a$coefficients["lipid", "p"]
}
qvals<-p.adjust(pvals,method="BH")
new2=data.frame(colnames(lipids_select),beta,se,HR,lower_CI,upper_CI,pvals,qvals)
colnames(new2)[1] <- "lipids_name"
write.table(new2,'longitu_chd_protLipids.csv',sep=',',row.names = F, col.names=T)

setwd(result)
res_risk <- read.csv('longitu_chd_riskLipids.csv')
res_prot <- read.csv('longitu_chd_protLipids.csv')
q_all <- p.adjust(c(res_risk$pvals, res_prot$pvals), method = "BH")
res_risk$qvals_withinSet <- res_risk$qvals
res_prot$qvals_withinSet <- res_prot$qvals
res_risk$qvals <- q_all[seq_len(nrow(res_risk))]
res_prot$qvals <- q_all[nrow(res_risk) + seq_len(nrow(res_prot))]
write.table(res_risk,'longitu_chd_riskLipids.csv',sep=',',row.names = F, col.names=T)
write.table(res_prot,'longitu_chd_protLipids.csv',sep=',',row.names = F, col.names=T)
