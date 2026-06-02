#KN-  liquid 2026 05 08, 19
#FINAL PLOTS FOR LIQUID PAPER - NP ONLY - clinical features
Sys.setenv(LANG = "en")
library(tidyverse)
library(car)
library(reshape2) #plot
library(ggprism) #plot
library(FSA) #posthoc
library(corrplot) #correlations
library(survival)
library(survminer)
library(timeROC)
library(survivalROC)
library(pROC)
library(gt)
library(glmnet)
library(brglm2)
library(htmlwidgets)
library(webshot)
library(magick)
library(patchwork)
#read RDS
LIQUID_DF_final <- readRDS("C:/Users/Ieva/rprojects/OTHER DATA/KN_LIQUID/liquid_20260415.RDS")
#leave lavage only
#notch2 the fullest data
LAVAGE_df <- LIQUID_DF_final%>%
  filter(!is.na(NOTCH2_NP)) #103 cases
#create a new groupings
LAVAGE_df <- LAVAGE_df %>%
  mutate(
    TYPE_BENIGN2 = if_else(TYPE %in% c("RSS", "BENIGN"),
                           "BENIGN",
                           TYPE)
  )
table(LAVAGE_df$TYPE_BENIGN2, useNA = "a") #now 31 benign
LAVAGE_df <- LAVAGE_df %>%
  mutate(
    TYPE_BENIGN3 = if_else(TYPE_BENIGN2 %in% c("HGSOC", "OTHER"),
                           "OC",
                           TYPE_BENIGN2)
  )
table(LAVAGE_df$TYPE_BENIGN3, useNA = "a") #now 31 benign
DATA <- c("NOTCH2_NP","CTNNB1_NP","DLL1_NP","HES1_NP" )

#SURVIVAL FULL KN###################################
#make df of only cases with survival
ALL_SURV_LIQUID <- LAVAGE_df %>%
  filter(!is.na(OS), !is.na(STATUS)) 
table(ALL_SURV_LIQUID$STATUS, useNA = "a") #28 dead, 62 NOT
#make factor of gene expresssion
DATA_names_P <- c("NOTCH2_NP","CTNNB1_NP","DLL1_NP","HES1_NP",
                  "NOTCH2_URINE", "CTNNB1_URINE","DLL1_URINE","HES1_URINE",
                  "NOTCH2_TUMOR","CTNNB1_TUMOR","DLL1_TUMOR","HES1_TUMOR" )
ALL_SURV_LIQUID <- ALL_SURV_LIQUID %>%
  mutate(
    across(all_of(DATA_names_P),
           ~ factor(if_else(. > median(., na.rm = TRUE), "High", "Low")),
           .names = "{.col}_f")
  )
#make only oc surv df
OC_SURV_LIQUID <- ALL_SURV_LIQUID %>%
  filter(TYPE %in%c("HGSOC", "OTHER")) 
table(OC_SURV_LIQUID$STATUS, useNA = "a") #22 is dead 38 not

##NP, ALL KN######################################################################
##univariable cox, KN, NP all ##################
cox_model_notch2_np <- coxph(Surv(OS, STATUS) ~ NOTCH2_NP, data = OC_SURV_LIQUID)
summary(cox_model_notch2_np)

cox_model_CTNNB1_np <- coxph(Surv(OS, STATUS) ~ CTNNB1_NP, data = OC_SURV_LIQUID)
summary(cox_model_CTNNB1_np)

cox_model_DLL1_np <- coxph(Surv(OS, STATUS) ~ DLL1_NP, data = OC_SURV_LIQUID)
summary(cox_model_DLL1_np)# 

cox_model_HES1_np <- coxph(Surv(OS, STATUS) ~ HES1_NP, data = OC_SURV_LIQUID)
summary(cox_model_HES1_np) #hr = 1.304 ci (1.022 - 1.664)  p =0.0327 

## Fit survival curves#################
#ctnnb1
OC_SURV_LIQUID$CTNNB1<-  OC_SURV_LIQUID$CTNNB1_NP_f 
p_ctnnb1 <- ggsurvplot(survfit(Surv(OS, STATUS) ~ CTNNB1, data = OC_SURV_LIQUID), 
           data = OC_SURV_LIQUID, pval = F,
           xlab = "Time (months)",palette = c(  "#E64164", "#002060"), 
           title=expression(  "Overall survival by " * italic("CTNNB1") * " expression in UL, all OC cases"))
p_ctnnb1$plot <- p_ctnnb1$plot + labs(subtitle = "HR = 0.989; 95% CI: 0.603 - 1.624; Log-rank  p = 0.966")
print(p_ctnnb1)
#notch2
OC_SURV_LIQUID$NOTCH2 <-  OC_SURV_LIQUID$NOTCH2_NP_f 
p_notch2 <- ggsDLL1p_notch2 <- ggsurvplot(survfit(Surv(OS, STATUS) ~ NOTCH2, data = OC_SURV_LIQUID), 
           data = OC_SURV_LIQUID, pval = F,
           xlab = "Time (months)",palette = c(  "#E64164", "#002060"), 
           title=expression(  "Overall survival by " * italic("NOTCH2") * " expression in UL, all OC cases"))
p_notch2$plot <- p_notch2$plot + labs(subtitle = "HR = 1.142; 95% CI: 0.798 - 1.634; Log-rank  p = 0.468")
print(p_notch2)

#DLL1
OC_SURV_LIQUID$DLL1 <- OC_SURV_LIQUID$DLL1_NP_f
p_dll1 <- ggsurvplot(survfit(Surv(OS, STATUS) ~ DLL1, data = OC_SURV_LIQUID), 
           data = OC_SURV_LIQUID, pval = F,
           xlab = "Time (months)",palette = c(  "#E64164", "#002060"), 
           title=expression(  "Overall survival by " * italic("DLL1") * " expression in UL, all OC cases"))
p_dll1$plot <- p_dll1$plot + labs(subtitle = "HR = 0.681; 95% CI: 0.450 - 1.031; Log-rank  p = 0.070")
print(p_dll1)


#HES1
OC_SURV_LIQUID$HES1 <- OC_SURV_LIQUID$HES1_NP_f
p_hes1 <- ggsurvplot(survfit(Surv(OS, STATUS) ~ HES1, data = OC_SURV_LIQUID), 
           data = OC_SURV_LIQUID, pval = F,
           xlab = "Time (months)",
           palette = c(  "#E64164", "#002060"), 
           title=expression(  "Overall survival by " * italic("HES1") * " expression in UL, all OC cases"))
p_hes1$plot <- p_hes1$plot + labs(subtitle = "HR = 1.304; 95% CI: 1.022 - 1.664; Log-rank  p = 0.033")
print(p_hes1)

combined_plot <-
( p_hes1$plot | p_dll1$plot ) /
( p_notch2$plot | p_ctnnb1$plot )

ggsave(
  filename = "c:/Users/Ieva/rprojects/outputs_all/LIQUID/OC_NP_survival_combined20260519.png",
  plot = combined_plot,
  width = 15,
  height = 10,
  dpi = 300
)
