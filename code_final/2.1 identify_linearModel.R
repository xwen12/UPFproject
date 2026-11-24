# Aim 2: UPF intake and 1,542 lipids, linear mixed models on pooled baseline + follow-up data.
# Output: result/P4P5all.csv (all lipids) and result/UPF_lipids_all.csv (42 named lipids, FDR < 0.05).

rm(list = ls())

proj_dir <- "path/to/proj16_UPF_rework"
raw_dir  <- "path/to/SHS_raw_data"
library(dplyr)
library(tidyr)

rawData <- raw_dir
readyData <- file.path(proj_dir, "data")
result <-file.path(proj_dir, "result")
tabsFigs <- file.path(proj_dir, "result", "tabsFigs")

setwd(readyData)
load("lipidomics_p4.Rda")
load("lipidomics_p5.Rda")
load("dat_AimOne_p4.Rda")
load("dat_AimOne_p5.Rda")

dat_AimOne_p4 <- dat_AimOne_p4 %>%
  select(-chd, -cvd,-TIME_chd,-TIME_cvd,-starts_with("incident"))

lipidomics_df_p4 <- as.data.frame(lipidomics_p4)
temp_1 <- lipidomics_df_p4 %>%
  inner_join(dat_AimOne_p4,by = c("ID" = "IDNO" ) )

lipidomics_df_p5 <- as.data.frame(lipidomics_p5)
temp_2 <- lipidomics_df_p5 %>%
  inner_join(dat_AimOne_p5,by = c("ID" = "IDNO" ) )

temp_1[,"Phase"] <- "Phase4"
temp_2[,"Phase"] <- "Phase5"

names(temp_1) <- sub('^S4', '', names(temp_1))
names(temp_2) <- sub('^S5', '', names(temp_2))
colnames(temp_1)[colnames(temp_1) == "PDM3"] <- "PA"
colnames(temp_2)[colnames(temp_2) == "ped_avg3day"] <- "PA"
temp_2$BOOKNUM <- NULL
temp_3 <- rbind(temp_1, temp_2)
temp_3$Phase <- as.factor(temp_3$Phase)

totaln <- length(unique(temp_3$ID))
dat <- temp_3 %>%
  mutate(across(c(ends_with("_P"), AGE, EDU, BMI, PA,DT_KCAL),
                ~ as.numeric(scale(.x, center = TRUE, scale = TRUE))))

dt <- dat
first_ind <- which(colnames(dt) == "AC(10:0)")
last_ind <- which(colnames(dt) == "Unknown_1230")
for(i in first_ind:last_ind){
  dt[,i]<-scale((dt[,i]),scale=TRUE,center=TRUE)
}

selected_lipids_names <- colnames(dt)[first_ind:last_ind]
(n_lipids=length(selected_lipids_names))
lipids_select<-dt[,first_ind:last_ind]

beta<-rep(0,n_lipids)
se<-rep(0,n_lipids)
upper_CI<-rep(0,n_lipids)
lower_CI<-rep(0,n_lipids)
pvals<-rep(0,n_lipids)
ve_o<-rep(0,n_lipids)
dt[,'FID']=as.factor(dt[,'FID'])
dt[,'ID']=as.factor(dt[,'ID'])
library(lme4)
library(lmerTest)
library(geepack)
library(r2glmm)
for(i in 1:n_lipids){
  lipid<-lipids_select[,i]
  fit<-lmer(lipid~nova4_kcal_P+AGE+SEX+CENTER+EDU+BMI+SMOKE+PA+ETOH+LIPIDLOW+DT_KCAL+Phase+(1|FID)+(1|ID),
            data = dt)
  a=coef(summary(fit))
  beta[i]<-a[2,1]
  se[i]<-a[2,2]
  lower_CI[i]<-beta[i]-1.96*a[2,2]
  upper_CI[i]<-beta[i]+1.96*a[2,2]
  pvals[i]<-a[2,5]
  ve=r2beta(fit,method = "sgv")
  ve_o[i]=ve[ve[,1]=='nova4_kcal_P','Rsq']
}
qvals = p.adjust(pvals, method = "BH")
known = stringr::str_detect(colnames(lipids_select), "Unknown", negate = T)
new=data.frame(lipids_name = colnames(lipids_select),lower_CI,upper_CI,pvals,qvals,beta,se,ve_o,known)
setwd(result)
write.table(new,'P4P5all.csv',sep=',',row.names = F, col.names=T)

setwd(result)
P4P5all <- read.table('P4P5all.csv',sep=',',header=T)

a <- P4P5all %>%
  filter(qvals <0.05)
UPF_lipids_risk <- P4P5all %>%
  filter(qvals <0.05, known,beta>0)
write.table(UPF_lipids_risk,'UPF_lipids_risk.csv',sep=',',row.names = F, col.names=T)

UPF_lipids_prot <- P4P5all %>%
  filter(qvals <0.05, known,beta<0)
write.table(UPF_lipids_prot,'UPF_lipids_prot.csv',sep=',',row.names = F, col.names=T)

UPF_lipids_all <- P4P5all %>%
  filter(qvals <0.05, known)

write.table(UPF_lipids_all,'UPF_lipids_all.csv',sep=',',row.names = F, col.names=T)
