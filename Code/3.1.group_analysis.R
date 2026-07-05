library(readr)
library(dplyr)
library(ggplot2)
library(car)
library(emmeans)
library(effectsize)

setwd("~/BNU_Helab/data/nm_subtyping_adhd/analysis_v2/CT_subtypes_v2/brain")

# Load data
ct = read.csv('./ct_merged_z.csv')
colnames(ct)[1:9]
ct[c('subgroups','sex','site','group','drug')] = lapply(ct[c('subgroups','sex','site','group','drug')], factor)
dem = ct[,c(1:7)]

# Three group ANCOVA
ancova_ct_3g<- data.frame()
for (i in c(8:226)){
  fix.effect <- aov(ct[[i]]~age+sex+site+drug+subgroups,data = ct)
  fix.effect.res <- summary(fix.effect)
  fix.effect.res  <- as.data.frame(fix.effect.res[[1]])
  Diag_effec <- fix.effect.res[c("subgroups"),c(4,5)]
  row.names(Diag_effec)<- names(ct)[i]
  ancova_ct_3g <- rbind.data.frame(ancova_ct_3g,Diag_effec)}
rm(i,fix.effect,fix.effect.res,Diag_effec)
# FDR 
ancova_ct_3g$P_FDR <- p.adjust(ancova_ct_3g$`Pr(>F)`,method = 'fdr',n=length(ancova_ct_3g$`Pr(>F)`))
# save
write.csv(ancova_ct_3g,"./deviation/ancova_ct_3g.csv")

# select significant roi from three group ANCOVA
sig_roi_ct = row.names(subset(ancova_ct_3g, P_FDR < 0.05))
ct_sig = ct[, colnames(ct) %in% sig_roi_ct]
ct_sig = cbind(dem,ct_sig)
# post hoc
ph.t.ct = data.frame(subgroups0_1=numeric(0),subgroups0_2=numeric(0), 
                     subgroups1_2=numeric(0))
ph.p.ct= data.frame(subgroups0_1=numeric(0),subgroups0_2=numeric(0), 
                    subgroups1_2=numeric(0))
for (i in c(8:225)){
  posthoc <- emmeans(aov(ct_sig[[i]]~age+sex+drug+subgroups,data = ct_sig),
                     pairwise~subgroups, adjust='fdr')
  posthoc_contr <- as.data.frame(posthoc[["contrasts"]])
  t = data.frame(subgroups0_1=numeric(0),subgroups0_2=numeric(0),subgroups1_2=numeric(0))
  p = data.frame(subgroups0_1=numeric(0),subgroups0_2=numeric(0),subgroups1_2=numeric(0))
  t[1,'subgroups0_1'] = posthoc_contr$t.ratio[1]
  t[1,'subgroups0_2'] = posthoc_contr$t.ratio[2]
  t[1,'subgroups1_2'] = posthoc_contr$t.ratio[3]
  p[1,'subgroups0_1'] = posthoc_contr$p.value[1]
  p[1,'subgroups0_2'] = posthoc_contr$p.value[2]
  p[1,'subgroups1_2'] = posthoc_contr$p.value[3]
  row.names(t) = names(ct_sig)[i]
  row.names(p) = names(ct_sig)[i]
  ph.t.ct <- rbind.data.frame(ph.t.ct,t)
  ph.p.ct <- rbind.data.frame(ph.p.ct,p)}
rm(i,posthoc,posthoc_contr,t,p)
# FDR posthoc p
ph.p.ct0_1 = p.adjust(ph.p.ct$subgroups0_1,method = 'fdr',n=length(ph.p.ct$subgroups0_1))
ph.p.ct0_2 = p.adjust(ph.p.ct$subgroups0_2,method = 'fdr',n=length(ph.p.ct$subgroups0_2))
ph.p.ct1_2 = p.adjust(ph.p.ct$subgroups1_2,method = 'fdr',n=length(ph.p.ct$subgroups1_2))
ph.p.ct_fdr = data.frame(subgroups0_1 = ph.p.ct0_1,
                         subgroups0_2 = ph.p.ct0_2,
                         subgroups1_2 = ph.p.ct1_2,
                         row.names = row.names(ph.p.ct))
rm(ph.p.ct0_1,ph.p.ct0_2,ph.p.ct1_2)
write.csv(ph.t.ct,"./ph_t_CT.csv")
write.csv(ph.p.ct_fdr,"./ph_p_CT_fdr.csv")


# Pairwise subgroups analysis
t.ct = data.frame(subgroups0_1=numeric(0),subgroups0_2=numeric(0), subgroups1_2=numeric(0))
p.ct= data.frame(subgroups0_1=numeric(0),subgroups0_2=numeric(0), subgroups1_2=numeric(0))

for (i in c(8:226)){
  posthoc <- emmeans(aov(ct[[i]]~age+sex+site+drug+subgroups,data = ct),
                     pairwise~subgroups, adjust='fdr')
  posthoc_contr <- as.data.frame(posthoc[["contrasts"]])
  t = data.frame(subgroups0_1=numeric(0),subgroups0_2=numeric(0),subgroups1_2=numeric(0))
  p = data.frame(subgroups0_1=numeric(0),subgroups0_2=numeric(0),subgroups1_2=numeric(0))
  t[1,'subgroups0_1'] = posthoc_contr$t.ratio[1]
  t[1,'subgroups0_2'] = posthoc_contr$t.ratio[2]
  t[1,'subgroups1_2'] = posthoc_contr$t.ratio[3]
  p[1,'subgroups0_1'] = posthoc_contr$p.value[1]
  p[1,'subgroups0_2'] = posthoc_contr$p.value[2]
  p[1,'subgroups1_2'] = posthoc_contr$p.value[3]
  row.names(t) = names(ct)[i]
  row.names(p) = names(ct)[i]
  t.ct <- rbind.data.frame(t.ct,t)
  p.ct <- rbind.data.frame(p.ct,p)}
rm(i,posthoc,posthoc_contr,t,p)
# FDR p_value
p.ct0_1 = p.adjust(p.ct$subgroups0_1,method = 'fdr',n=length(p.ct$subgroups0_1))
p.ct0_2 = p.adjust(p.ct$subgroups0_2,method = 'fdr',n=length(p.ct$subgroups0_2))
p.ct1_2 = p.adjust(p.ct$subgroups1_2,method = 'fdr',n=length(p.ct$subgroups1_2))
p.ct_fdr = data.frame(subgroups0_1 = p.ct0_1,
                      subgroups0_2 = p.ct0_2,
                      subgroups1_2 = p.ct1_2,
                      row.names = row.names(p.ct))
rm(p.ct0_1,p.ct0_2,p.ct1_2)
write.csv(t.ct,"./deviation/tval_CT_pairwise.csv")
write.csv(p.ct_fdr,"./deviation/pval_fdr_CT_pairwise.csv")


# Case-Control
ct_2g = data.frame()
for (i in c(8:226)){
  fix.effect <- aov(ct[[i]]~age+sex+site+drug+group,data = ct)
  fix.effect.res <- emmeans(fix.effect,pairwise~group)
  dignose_effe  <- as.data.frame(fix.effect.res[["contrasts"]])
  row.names(dignose_effe) <- names(ct)[i]
  ct_2g <- rbind.data.frame(ct_2g,dignose_effe)}
rm(i,fix.effect,fix.effect.res,dignose_effe)
# FDR 
ct_2g$P_FDR <- p.adjust(ct_2g$p.value,method = 'fdr',n=length(ct_2g$p.value))
write.csv(ct_2g,"./deviation/tval_CT_casecontrol.csv")


# Effect size 
# case-control
tval = read.csv('./combat/tstats_ct_casecontrol.csv')
res = aov(ct[[10]]~age+sex+drug+group,data = ct) # check the df.residuals
dval = t_to_d(t=tval$t.ratio,paired = F,ci = 0.95,alternative = 'two.sided',n = 408,df = 403)
dval$roi = tval$X
write.csv(dval,"./combat/ES_ct_casecontrol.csv",row.names = F)

# Three group ANCOVA
fval = read.csv('./deviation/ancova_ct_3g.csv')
summary(aov(ct[[10]]~age+sex+site+drug+subgroups,data = ct)) # check the df_group and df_residuals(error)
# eta_sqr = F * df1/(F*df1 + df.error)
eta_val = data.frame(roi = fval$X, eta.val = fval$F.value * 2/(fval$F.value*2 + 684))
write.csv(eta_val,"./deviation/ES_ct_3g.csv",row.names = F)

# Pairwise t-value
tval = read.csv('./deviation/tval_CT_pairwise.csv')
res = emmeans(aov(ct[[10]]~age+sex+site+drug+subgroups,data = ct),pairwise~subgroups, adjust='fdr')# check the df_residuals
res[["emmeans"]] # check the df
dval_10 = t_to_d(t= -1*(tval$subgroups0_1),paired = F,ci = 0.95,alternative = 'two.sided',df = 684)
dval_20 = t_to_d(t= -1*(tval$subgroups0_2),paired = F,ci = 0.95,alternative = 'two.sided',df = 684)
dval_12 = t_to_d(t= tval$subgroups1_2,paired = F,ci = 0.95,alternative = 'two.sided',df = 684)
dval = data.frame(roi = tval$X,
                  d.subtype1_0 = dval_10$d,
                  d.subtype2_0 = dval_20$d,
                  d.subtype1_2 = dval_12$d)
write.csv(dval,"./deviation/ES_ct_pairwise.csv",row.names = F)
