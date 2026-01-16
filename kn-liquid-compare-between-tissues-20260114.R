#KN-  liquid 2026 01 14
#NP, URINE, PLASMA only the 15 cases data
Sys.setenv(LANG = "en")
library(readxl)
library(tidyverse)
library(car)
library(reshape2) #plot
library(ggprism) #plot
library(survival)
library(survminer) #ggsurvplot
library(pROC)
library(gt)
library(dplyr)
#read RDS
LIQUID_DF <- readRDS("C:/Users/Ieva/rprojects/OTHER DATA/KN_LIQUID/liquid_20260114.RDS") #NEW WITH UPDATED NP kn-97
#make 15 case df
#remove cases that is NA in their type: 
LIQUID_DF_15 <- LIQUID_DF %>%
  filter(!is.na(HES1_URINE))

#COMPARE BETWEEN TISSUES / liquid biopsies########################
#notch2
NOTCH2 <- c("NOTCH2_URINE", "NOTCH2_NP", "NOTCH2_P_norm","NOTCH2_TUMOR")
notch2_long <- LIQUID_DF_15  %>%
  dplyr::select(c("Laboratorinis kodas", NOTCH2)) %>%
  pivot_longer(
    cols = NOTCH2,
    names_to = "gene",
    values_to = "expression"
  )
#set order
notch2_long$gene <- factor(
  notch2_long$gene,
  levels = c("NOTCH2_P_norm", "NOTCH2_URINE", "NOTCH2_NP", "NOTCH2_TUMOR") 
)

ggplot(notch2_long, aes(x = gene, y = expression, group = `Laboratorinis kodas`)) +
  geom_point(size = 2) +
  geom_line(alpha = 0.5) +
  theme_classic() +
  labs(
    x = "Tissue",
    y = "Gene expression",
    title = "Paired gene expression across tissues"
  )

#hes1
hes1 <- c("HES1_URINE", "HES1_NP", "HES1_P_norm","HES1_TUMOR")
hes1_long <- LIQUID_DF_15  %>%
  dplyr::select(c("Laboratorinis kodas", hes1)) %>%
  pivot_longer(
    cols = hes1,
    names_to = "gene",
    values_to = "expression"
  )

#set order
hes1_long$gene <- factor(
  hes1_long$gene,
  levels = c("HES1_P_norm","HES1_URINE", "HES1_NP", "HES1_TUMOR") 
)

ggplot(hes1_long, aes(x = gene, y = expression, group = `Laboratorinis kodas`)) +
  geom_point(size = 2) +
  geom_line(alpha = 0.5) +
  theme_classic() +
  labs(
    x = "Tissue",
    y = "Gene expression",
    title = "Paired gene expression across tissues"
  )

#CTNNB1
CTNNB1 <- c("CTNNB1_URINE", "CTNNB1_NP", "CTNNB1_P_norm","CTNNB1_TUMOR")
CTNNB1_long <- LIQUID_DF_15  %>%
  dplyr::select(c("Laboratorinis kodas", CTNNB1)) %>%
  pivot_longer(
    cols = CTNNB1,
    names_to = "gene",
    values_to = "expression"
  )
#set order
CTNNB1_long$gene <- factor(
  CTNNB1_long$gene,
  levels = c( "CTNNB1_P_norm","CTNNB1_URINE", "CTNNB1_NP","CTNNB1_TUMOR") 
)


ggplot(CTNNB1_long, aes(x = gene, y = expression, group = `Laboratorinis kodas`)) +
  geom_point(size = 2) +
  geom_line(alpha = 0.5) +
  theme_classic() +
  labs(
    x = "Tissue",
    y = "Gene expression",
    title = "Paired gene expression across tissues"
  )

#dll1
DLL1 <- c("DLL1_URINE", "DLL1_NP", "DLL1_P_norm","DLL1_TUMOR")
DLL1_long <- LIQUID_DF_15  %>%
  dplyr::select(c("Laboratorinis kodas", DLL1)) %>%
  pivot_longer(
    cols = DLL1,
    names_to = "gene",
    values_to = "expression"
  )

#set order
DLL1_long$gene <- factor(
  DLL1_long$gene,
  levels = c("DLL1_P_norm", "DLL1_URINE", "DLL1_NP", "DLL1_TUMOR") 
)

ggplot(DLL1_long, aes(x = gene, y = expression, group = `Laboratorinis kodas`)) +
  geom_point(size = 2) +
  geom_line(alpha = 0.5) +
  theme_classic() +
  labs(
    x = "Tissue",
    y = "Gene expression",
    title = "Paired gene expression across tissues"
  )

#NORMALCY, VARIANCE ###################################################
by(LIQUID_DF_15$NOTCH2_NP, LIQUID_DF_15$TYPE, shapiro.test)
by(LIQUID_DF_15$DLL1_NP, LIQUID_DF_15$TYPE, shapiro.test)
by(LIQUID_DF_15$HES1_NP, LIQUID_DF_15$TYPE, shapiro.test)
by(LIQUID_DF_15$CTNNB1_NP, LIQUID_DF_15$TYPE, shapiro.test)

by(LIQUID_DF_15$NOTCH2_URINE, LIQUID_DF_15$TYPE, shapiro.test)
#by(LIQUID_DF_15$DLL1_URINE, LIQUID_DF_15$TYPE, shapiro.test) #Not enough
by(LIQUID_DF_15$HES1_URINE, LIQUID_DF_15$TYPE, shapiro.test)
by(LIQUID_DF_15$CTNNB1_URINE, LIQUID_DF_15$TYPE, shapiro.test)

by(LIQUID_DF_15$NOTCH2_TUMOR, LIQUID_DF_15$TYPE, shapiro.test)#OTHER not normal
by(LIQUID_DF_15$DLL1_TUMOR, LIQUID_DF_15$TYPE, shapiro.test)
by(LIQUID_DF_15$HES1_TUMOR, LIQUID_DF_15$TYPE, shapiro.test)
by(LIQUID_DF_15$CTNNB1_TUMOR, LIQUID_DF_15$TYPE, shapiro.test)

#test variance
car::leveneTest(LIQUID_DF_15$NOTCH2_NP ~ LIQUID_DF_15$TYPE, center = median)
car::leveneTest(LIQUID_DF_15$DLL1_NP ~ LIQUID_DF_15$TYPE, center = median)
car::leveneTest(LIQUID_DF_15$HES1_NP ~ LIQUID_DF_15$TYPE, center = median)
car::leveneTest(LIQUID_DF_15$CTNNB1_NP ~ LIQUID_DF_15$TYPE, center = median)

car::leveneTest(LIQUID_DF_15$NOTCH2_URINE ~ LIQUID_DF_15$TYPE, center = median)
car::leveneTest(LIQUID_DF_15$DLL1_URINE ~ LIQUID_DF_15$TYPE, center = median) #DLL1_URINE 
car::leveneTest(LIQUID_DF_15$HES1_URINE ~ LIQUID_DF_15$TYPE, center = median)
car::leveneTest(LIQUID_DF_15$CTNNB1_URINE ~ LIQUID_DF_15$TYPE, center = median)

car::leveneTest(LIQUID_DF_15$NOTCH2_TUMOR ~ LIQUID_DF_15$TYPE, center = median)
car::leveneTest(LIQUID_DF_15$DLL1_TUMOR ~ LIQUID_DF_15$TYPE, center = median)
car::leveneTest(LIQUID_DF_15$HES1_TUMOR ~ LIQUID_DF_15$TYPE, center = median)
car::leveneTest(LIQUID_DF_15$CTNNB1_TUMOR ~ LIQUID_DF_15$TYPE, center = median)

#t.tests##################################
wilcox.test(LIQUID_DF_15$NOTCH2_TUMOR ~ LIQUID_DF_15$TYPE)
t.test(LIQUID_DF_15$DLL1_TUMOR ~ LIQUID_DF_15$TYPE, var.equal = TRUE) #0.01185
t.test(LIQUID_DF_15$HES1_TUMOR ~ LIQUID_DF_15$TYPE, var.equal = TRUE) 
t.test(LIQUID_DF_15$CTNNB1_TUMOR ~ LIQUID_DF_15$TYPE, var.equal = TRUE)

t.test(LIQUID_DF_15$NOTCH2_URINE ~ LIQUID_DF_15$TYPE, var.equal = TRUE)
#t.test(LIQUID_DF_15$DLL1_URINE  ~ LIQUID_DF_15$TYPE, var.equal = FALSE) #not enough
t.test(LIQUID_DF_15$HES1_URINE  ~ LIQUID_DF_15$TYPE, var.equal = TRUE) #
t.test(LIQUID_DF_15$CTNNB1_URINE  ~ LIQUID_DF_15$TYPE, var.equal = TRUE)

t.test(LIQUID_DF_15$NOTCH2_NP ~ LIQUID_DF_15$TYPE, var.equal = TRUE)
t.test(LIQUID_DF_15$DLL1_NP  ~ LIQUID_DF_15$TYPE, var.equal = FALSE)
t.test(LIQUID_DF_15$HES1_NP  ~ LIQUID_DF_15$TYPE, var.equal = TRUE) 
t.test(LIQUID_DF_15$CTNNB1_NP  ~ LIQUID_DF_15$TYPE, var.equal = TRUE)

#boxplots TUMOR 15 cases#########################################
#melt table for expression
TUMOR_table <- melt(LIQUID_DF_15[, c(12,31:34)], id.vars="TYPE",
                measure.vars=c("NOTCH2_TUMOR","CTNNB1_TUMOR","DLL1_TUMOR","HES1_TUMOR"))
#Tribble 2 groups for grade
each.vs.ref_sig2_g <-  tibble::tribble(
  ~group1, ~group2, ~p.adj,   ~y.position, ~variable,
  "HGSOC",   "OTHER", 0.01185, -2, "DLL1_TUMOR", #stjudent's
)
#make order
TUMOR_table$variable <- factor(
  TUMOR_table$variable,
  levels = c("NOTCH2_TUMOR","CTNNB1_TUMOR","DLL1_TUMOR","HES1_TUMOR")  
)

#TUMOR boxplot
custom_colors_grade <- c("HGSOC" = "deeppink","OTHER" = "lightblue") 
TUMOR_plot <- ggplot(TUMOR_table, aes(x=TYPE , y=value, fill = variable)) +
  geom_boxplot( outlier.shape = NA , alpha=0.3, aes(fill = TYPE )) +
  geom_jitter(aes(color = TYPE ), size=1, alpha=0.5) +
  ylab(label = expression("Gene expression, normalized to " * italic("GAPDH"))) + 
  facet_wrap(.~ variable, nrow = 2, scales = "free") +
  add_pvalue(each.vs.ref_sig2_g, label = "p.adj") + #pvalue
  theme_minimal()+
  theme(
    strip.text.x = element_text(
      size = 12, face = "bold.italic"
    ),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5))+
  labs(x=NULL)+
  stat_boxplot(geom ='errorbar')+
  scale_fill_manual(values = custom_colors_grade) +
  scale_color_manual(values = custom_colors_grade) +
  scale_y_continuous(labels = function(x) 
    gsub("-", "\u2212", as.character(x))) #add long "-" signs

TUMOR_plot

#boxplots URINE ############################
#melt table for expression
URINE_table <- melt(LIQUID_DF_15[, c(12,19:22)], id.vars="TYPE",
                    measure.vars=c("NOTCH2_URINE","CTNNB1_URINE","DLL1_URINE","HES1_URINE"))

#TUMOR boxplot
custom_colors_grade <- c("HGSOC" = "lightpink","OTHER" = "darkblue") 
URINE_plot <- ggplot(URINE_table, aes(x=TYPE , y=value, fill = variable)) +
  geom_boxplot( outlier.shape = NA , alpha=0.3, aes(fill = TYPE )) +
  geom_jitter(aes(color = TYPE ), size=1, alpha=0.5) +
  ylab(label = expression("Gene expression, normalized to  " * italic("GAPDH"))) + 
  facet_wrap(.~ variable, nrow = 2, scales = "free") +
  #add_pvalue(each.vs.ref_sig2_g, label = "p.adj") + #pvalue
  theme_minimal()+
  theme(
    strip.text.x = element_text(
      size = 12, face = "bold.italic"
    ),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5))+
  labs(x=NULL)+
  stat_boxplot(geom ='errorbar')+
  scale_fill_manual(values = custom_colors_grade) +
  scale_color_manual(values = custom_colors_grade) +
  scale_y_continuous(labels = function(x) 
    gsub("-", "\u2212", as.character(x))) #add long "-" signs

URINE_plot

#boxplots NP 15 cases############################################ 
#melt table for expression
NP_table <- melt(LIQUID_DF_15[, c(12,15:18)], id.vars="TYPE",
                 measure.vars=c("NOTCH2_NP","CTNNB1_NP","DLL1_NP","HES1_NP"))
#NP boxplot
custom_colors_grade <- c("HGSOC" = "red","OTHER" = "turquoise") 
NP_plot <- ggplot(NP_table, aes(x=TYPE , y=value, fill = variable)) +
  geom_boxplot( outlier.shape = NA , alpha=0.3, aes(fill = TYPE )) +
  geom_jitter(aes(color = TYPE ), size=1, alpha=0.5) +
  ylab(label = expression("Gene expression, normalized to  " * italic("GAPDH"))) + 
  facet_wrap(.~ variable, nrow = 2, scales = "free") +
  #add_pvalue(each.vs.ref_sig2_g, label = "p.adj") + #pvalue
  theme_minimal()+
  theme(
    strip.text.x = element_text(
      size = 12, face = "bold.italic"
    ),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5))+
  labs(x=NULL)+
  stat_boxplot(geom ='errorbar')+
  scale_fill_manual(values = custom_colors_grade) +
  scale_color_manual(values = custom_colors_grade) +
  scale_y_continuous(labels = function(x) 
    gsub("-", "\u2212", as.character(x))) #add long "-" signs

NP_plot

#SURVIVAL###################################
#make df of only cases with survival
ALL_SURV_LIQUID <- LIQUID_DF_15 %>%
  filter(!is.na(OS), !is.na(STATUS)) #58 observations
table(ALL_SURV_LIQUID$STATUS, useNA = "a") #2 dead, 12 NOT
#make factor of gene expresssion
DATA_names_P <- c("NOTCH2_NP","CTNNB1_NP","DLL1_NP","HES1_NP",
                  "NOTCH2_URINE", "CTNNB1_URINE","DLL1_URINE","HES1_URINE",
                  "NOTCH2_P_norm","CTNNB1_P_norm", "DLL1_P_norm","HES1_P_norm",
                  "NOTCH2_TUMOR","CTNNB1_TUMOR","DLL1_TUMOR","HES1_TUMOR" )
ALL_SURV_LIQUID <- ALL_SURV_LIQUID %>%
  mutate(
    across(all_of(DATA_names_P),
           ~ factor(if_else(. > median(., na.rm = TRUE), "High", "Low")),
           .names = "{.col}_f")
  )
## Fit survival curves np, 15 cases#################
ggsurvplot(survfit(Surv(OS, STATUS) ~ CTNNB1_NP_f, data = ALL_SURV_LIQUID), 
           data = ALL_SURV_LIQUID, pval = TRUE,
           title ="Overall Survival by CTNNB1 NP expression in 15 cases")

ggsurvplot(survfit(Surv(OS, STATUS) ~ NOTCH2_NP_f, data = ALL_SURV_LIQUID), 
           data = ALL_SURV_LIQUID, pval = TRUE, 
           title = "Overall Survival by NOTCH2 NP expression in 15 cases")

ggsurvplot(survfit(Surv(OS, STATUS) ~ DLL1_NP_f, data = ALL_SURV_LIQUID), 
           data = ALL_SURV_LIQUID, pval = TRUE,
           title = "Overall Survival by DLL1 NP expression in 15 cases")

ggsurvplot(survfit(Surv(OS, STATUS) ~ HES1_NP_f, data = ALL_SURV_LIQUID), 
           data = ALL_SURV_LIQUID, pval = TRUE,
           title ="Overall Survival by HES1 NP expression in 15 cases" )

## Fit survival curves URINE#################
ggsurvplot(survfit(Surv(OS, STATUS) ~ CTNNB1_URINE_f, data = ALL_SURV_LIQUID), 
           data = ALL_SURV_LIQUID, pval = TRUE,
           title ="Overall Survival by CTNNB1 URINE expression in 15 cases")

ggsurvplot(survfit(Surv(OS, STATUS) ~ NOTCH2_URINE_f, data = ALL_SURV_LIQUID), 
           data = ALL_SURV_LIQUID, pval = TRUE,
           title ="Overall Survival by NOTCH2 URINE expression in 15 cases")

ggsurvplot(survfit(Surv(OS, STATUS) ~ DLL1_URINE_f, data = ALL_SURV_LIQUID), 
           data = ALL_SURV_LIQUID, pval = TRUE,
           title ="Overall Survival by DLL1 URINE expression in 15 cases")

ggsurvplot(survfit(Surv(OS, STATUS) ~ HES1_URINE_f, data = ALL_SURV_LIQUID), 
           data = ALL_SURV_LIQUID, pval = TRUE,
           title ="Overall Survival by HES1 URINE expression in 15 cases")

#FISHER URINE ####################################################
#make factors in liquid df 15
table(LIQUID_DF_15$TYPE)
LIQUID_DF_15$TYPE <- droplevels(LIQUID_DF_15$TYPE)
LIQUID_DF_15 <- LIQUID_DF_15 %>%
  mutate(
    across(all_of(DATA_names_P),
           ~ factor(if_else(. > median(., na.rm = TRUE), "High", "Low")),
           .names = "{.col}_f")
  )
#fisher tests on urine
fisher.test(table(LIQUID_DF_15$TYPE, LIQUID_DF_15$NOTCH2_URINE_f))
fisher.test(table(LIQUID_DF_15$TYPE, LIQUID_DF_15$CTNNB1_URINE_f))
fisher.test(table(LIQUID_DF_15$TYPE, LIQUID_DF_15$HES1_URINE_f))
fisher.test(table(LIQUID_DF_15$TYPE, LIQUID_DF_15$DLL1_URINE_f))

##PLOT URINE NOTCH2########################
# Fisher's exact test
tab <- table(LIQUID_DF_15$TYPE, LIQUID_DF_15$NOTCH2_URINE_f)
ft <- fisher.test(tab)
pval <- signif(ft$p.value, 3)  # rounded p-value

# Convert to percentages by column
tab_perc <- prop.table(tab, margin = 2) * 100

# Stacked bar plot with italic NOTCH2 in title
bp <- barplot(tab_perc, col = c("skyblue", "salmon"),
              legend = TRUE,
              main = expression("URINE expression of " * italic(NOTCH2)),
              sub = paste0("Fisher's exact test p = ", pval),
              ylab = "Percentage",
              beside = FALSE)

# Add percentages on bars
for(i in 1:nrow(tab_perc)){
  text(x = bp, y = apply(tab_perc[1:i, , drop = FALSE], 2, sum) - tab_perc[i, ]/2,
       labels = paste0(round(tab_perc[i, ], 1), "%"), cex = 0.8, col = "black")
}

##PLOT URINE HES1####################
# Fisher's exact test
tab1 <- table(LIQUID_DF_15$TYPE, LIQUID_DF_15$HES1_URINE_f)
ft1 <- fisher.test(tab1)
pval1 <- signif(ft1$p.value, 3)  # rounded p-value

# Convert to percentages by column
tab_perc1 <- prop.table(tab1, margin = 2) * 100

# Stacked bar plot with italic HES1 in title
bp1 <- barplot(tab_perc1, col = c("skyblue", "salmon"),
               legend = TRUE,
               main = expression("URINE expression of " * italic(HES1)),
               sub = paste0("Fisher's exact test p = ", pval1),
               ylab = "Percentage",
               beside = FALSE)

# Add percentages on bars
for(i in 1:nrow(tab_perc1)){
  text(x = bp1, y = apply(tab_perc1[1:i, , drop = FALSE], 2, sum) - tab_perc1[i, ]/2,
       labels = paste0(round(tab_perc1[i, ], 1), "%"), cex = 0.8, col = "black")
}

##PLOT URINE CTNNB1########################
# Fisher's exact test
tab2 <- table(LIQUID_DF_15$TYPE, LIQUID_DF_15$CTNNB1_URINE_f)
ft2 <- fisher.test(tab2)
pval2 <- signif(ft2$p.value, 3)  # rounded p-value

# Convert to percentages by column
tab_perc2 <- prop.table(tab2, margin = 2) * 100

# Stacked bar plot with italic CTNNB1 in title
bp2 <- barplot(tab_perc2, col = c("skyblue", "salmon"),
               legend = TRUE,
               main = expression("URINE expression of " * italic(CTNNB1)),
               sub = paste0("Fisher's exact test p = ", pval2),
               ylab = "Percentage",
               beside = FALSE)

# Add percentages on bars
for(i in 1:nrow(tab_perc2)){
  text(x = bp2, y = apply(tab_perc2[1:i, , drop = FALSE], 2, sum) - tab_perc2[i, ]/2,
       labels = paste0(round(tab_perc2[i, ], 1), "%"), cex = 0.8, col = "black")
}

##PLOT URINE DLL1########################
# Fisher's exact test
tab3 <- table(LIQUID_DF_15$TYPE, LIQUID_DF_15$DLL1_URINE_f)
ft3 <- fisher.test(tab3)
pval3 <- signif(ft3$p.value, 3)  # rounded p-value

# Convert to percentages by column
tab_perc3 <- prop.table(tab2, margin = 2) * 100

# Stacked bar plot with italic DLL1 in title
bp3 <- barplot(tab_perc3, col = c("skyblue", "salmon"),
               legend = TRUE,
               main = expression("URINE expression of " * italic(DLL1)),
               sub = paste0("Fisher's exact test p = ", pval3),
               ylab = "Percentage",
               beside = FALSE)

# Add percentages on bars
for(i in 1:nrow(tab_perc3)){
  text(x = bp3, y = apply(tab_perc3[1:i, , drop = FALSE], 2, sum) - tab_perc3[i, ]/2,
       labels = paste0(round(tab_perc3[i, ], 1), "%"), cex = 0.8, col = "black")
}

#FISHER NP, 15 CASES#####################
#fisher tests on NP, ONLY THE CASES WITH OTHER LIQUID SAMPLES
fisher.test(table(LIQUID_DF_15$TYPE, LIQUID_DF_15$NOTCH2_NP_f))
fisher.test(table(LIQUID_DF_15$TYPE, LIQUID_DF_15$CTNNB1_NP_f))
fisher.test(table(LIQUID_DF_15$TYPE, LIQUID_DF_15$HES1_NP_f))
fisher.test(table(LIQUID_DF_15$TYPE, LIQUID_DF_15$DLL1_NP_f))

#PLASMA SURVIVAL##########################################
table(ALL_SURV_LIQUID$NOCTH2_P, ALL_SURV_LIQUID$STATUS)
table(ALL_SURV_LIQUID$DLL1_P, ALL_SURV_LIQUID$STATUS)
table(ALL_SURV_LIQUID$HES1_P, ALL_SURV_LIQUID$STATUS)
table(ALL_SURV_LIQUID$CTNNB1_P, ALL_SURV_LIQUID$STATUS)

## Fit survival curves PLASMA
ggsurvplot(survfit(Surv(OS, STATUS) ~ CTNNB1_P, data = ALL_SURV_LIQUID), 
           data = ALL_SURV_LIQUID, pval = TRUE,
           title ="Overall Survival by CTNNB1 PLASMA expression in 15 cases")
ggsurvplot(survfit(Surv(OS, STATUS) ~ HES1_P, data = ALL_SURV_LIQUID), 
           data = ALL_SURV_LIQUID, pval = TRUE,
           title ="Overall Survival by HES1 PLASMA expression in 15 cases")
ggsurvplot(survfit(Surv(OS, STATUS) ~ DLL1_P, data = ALL_SURV_LIQUID), 
           data = ALL_SURV_LIQUID, pval = TRUE,
           title ="Overall Survival by DLL1 PLASMA expression in 15 cases")
ggsurvplot(survfit(Surv(OS, STATUS) ~ NOCTH2_P, data = ALL_SURV_LIQUID), 
           data = ALL_SURV_LIQUID, pval = TRUE,
           title ="Overall Survival by NOTCH2 PLASMA expression in 15 cases")

##univariable cox, KN, plasma categorical data##################
cox_model_notch2_plasma <- coxph(Surv(OS, STATUS) ~ NOCTH2_P, data = ALL_SURV_LIQUID)
summary(cox_model_notch2_plasma)

cox_model_CTNNB1_plasma <- coxph(Surv(OS, STATUS) ~ CTNNB1_P, data = ALL_SURV_LIQUID)
summary(cox_model_CTNNB1_plasma)

#cox_model_DLL1_plasma <- coxph(Surv(OS, STATUS) ~ DLL1_P, data = ALL_SURV_LIQUID)
#summary(cox_model_DLL1_plasma)# not enough data

cox_model_HES1_plasma <- coxph(Surv(OS, STATUS) ~ HES1_P, data = ALL_SURV_LIQUID)
summary(cox_model_HES1_plasma) #mega dideli skaičiai

##univariable cox, KN, plasma numbered data##################
cox_model_HES1_plasma <- coxph(Surv(OS, STATUS) ~ HES1_P_norm, data = ALL_SURV_LIQUID)
summary(cox_model_HES1_plasma) #more managable skaičiai

##median survival plasma, KN categorical data###########################
hes1_fit <- survfit(Surv(OS, STATUS) ~ HES1_P, data = ALL_SURV_LIQUID)
summary(hes1_fit)$table #gives months
summary(hes1_fit, times = 12)$surv #gives survival prob at 1 year
summary(hes1_fit, times = 36)$surv #gives survival prob at 3 yrs
summary(hes1_fit, times = 60)$surv #gives survival prob at 5 yrs
#what is median survival overall?
fit_all <- survfit(Surv(OS, STATUS) ~ 1, data = ALL_SURV_LIQUID)
summary(fit_all)$table #not reached

##median survival plasma, KN numbered data###########################
hes1_fit_norm <- survfit(Surv(OS, STATUS) ~ HES1_P_norm, data = ALL_SURV_LIQUID)
summary(hes1_fit_norm)$table #gives months
summary(hes1_fit_norm, times = 12)$surv #gives survival prob at 1 year
summary(hes1_fit_norm, times = 36)$surv #gives survival prob at 3 yrs
summary(hes1_fit_norm, times = 60)$surv #gives survival prob at 5 yrs
#what is median survival overall?
fit_all <- survfit(Surv(OS, STATUS) ~ 1, data = ALL_SURV_LIQUID)
summary(fit_all)$table #not reached
#ant tiek mažai data kad net neskaičiuoja kaip numbered
ALL_SURV_LIQUID$HES1_P_norm

#PLASMA FIHSER TESTS###########################################
#fisher tests on urine
fisher.test(table(LIQUID_DF_15$TYPE, LIQUID_DF_15$NOCTH2_P))
fisher.test(table(LIQUID_DF_15$TYPE, LIQUID_DF_15$CTNNB1_P))
fisher.test(table(LIQUID_DF_15$TYPE, LIQUID_DF_15$HES1_P))
#fisher.test(table(LIQUID_DF_15$TYPE, LIQUID_DF_15$DLL1_P))#NOT ENOUGH DATA

##PLOT PLASMA NOTCH2########################
# Fisher's exact test
tab6 <- table(LIQUID_DF_15$TYPE, LIQUID_DF_15$NOCTH2_P)
ft6 <- fisher.test(tab6)
pval6 <- signif(ft6$p.value, 3)  # rounded p-value

# Convert to percentages by column
tab_perc6 <- prop.table(tab6, margin = 2) * 100

# Stacked bar plot with italic NOTCH2 in title
bp6 <- barplot(tab_perc6, col = c("skyblue", "salmon"),
               legend = TRUE,
               main = expression("PLASMA expression of " * italic(NOTCH2)),
               sub = paste0("Fisher's exact test p = ", pval6),
               ylab = "Percentage",
               beside = FALSE)

# Add percentages on bars
for(i in 1:nrow(tab_perc6)){
  text(x = bp6, y = apply(tab_perc6[1:i, , drop = FALSE], 2, sum) - tab_perc6[i, ]/2,
       labels = paste0(round(tab_perc6[i, ], 1), "%"), cex = 0.8, col = "black")
}
##PLOT PLASMA CTNNB1########################
# Fisher's exact test
tab7 <- table(LIQUID_DF_15$TYPE, LIQUID_DF_15$CTNNB1_P)
ft7 <- fisher.test(tab7)
pval7 <- signif(ft7$p.value, 3)  # rounded p-value

# Convert to percentages by column
tab_perc7 <- prop.table(tab7, margin = 2) * 100

# Stacked bar plot with italic CTNNB1 in title
bp7 <- barplot(tab_perc7, col = c("skyblue", "salmon"),
               legend = TRUE,
               main = expression("PLASMA expression of " * italic(CTNNB1)),
               sub = paste0("Fisher's exact test p = ", pval7),
               ylab = "Percentage",
               beside = FALSE)

# Add percentages on bars
for(i in 1:nrow(tab_perc7)){
  text(x = bp7, y = apply(tab_perc7[1:i, , drop = FALSE], 2, sum) - tab_perc7[i, ]/2,
       labels = paste0(round(tab_perc7[i, ], 1), "%"), cex = 0.8, col = "black")
}
##PLOT PLASMA HES1########################
# Fisher's exact test
tab8 <- table(LIQUID_DF_15$TYPE, LIQUID_DF_15$HES1_P)
ft8 <- fisher.test(tab8)
pval8 <- signif(ft8$p.value, 3)  # rounded p-value

# Convert to percentages by column
tab_perc8<- prop.table(tab8, margin = 2) * 100

# Stacked bar plot with italic HES1 in title
bp8 <- barplot(tab_perc8, col = c("skyblue", "salmon"),
               legend = TRUE,
               main = expression("PLASMA expression of " * italic(HES1)),
               sub = paste0("Fisher's exact test p = ", pval8),
               ylab = "Percentage",
               beside = FALSE)

# Add percentages on bars
for(i in 1:nrow(tab_perc8)){
  text(x = bp8, y = apply(tab_perc8[1:i, , drop = FALSE], 2, sum) - tab_perc8[i, ]/2,
       labels = paste0(round(tab_perc8[i, ], 1), "%"), cex = 0.8, col = "black")
}

#NP ROC, 15 CASES ##############################################################
LIQUID_DF_15$TYPE <- droplevels(LIQUID_DF_15$TYPE) #drop unused
LIQUID_DF_15$TYPE  <- relevel(LIQUID_DF_15$TYPE , ref = "OTHER") #set control
roc_results_np_15 <- lapply(c("NOTCH2_NP","CTNNB1_NP","DLL1_NP","HES1_NP"), function(col) {
  roc(response = LIQUID_DF_15$TYPE, predictor = LIQUID_DF_15[[col]])})
names(roc_results_np_15) <- c("NOTCH2_NP","CTNNB1_NP","DLL1_NP","HES1_NP")
roc_results_np_15
#extract the aucs
AUC_results_np_15 <- sapply(roc_results_np_15, function(roc_obj) {auc(roc_obj)})
AUC_results_np_15 #extracted aucs

#roc figure 
roc_plot <- function() {
  par(pty = "s") #sets square
  plot.roc(roc_results_np_15[["NOTCH2_NP"]], print.auc = F, col = "#dcbeff",
           cex.main=0.8, 
           main ="HGOSC atskyrimas nuo kitų KV, NP, 15 atvejų",
           xlab = "1 - Specifiškumas", 
           ylab = "Jautrumas", 
           legacy.axes = T) #title
  lines(roc_results_np_15[["CTNNB1_NP"]], col = "#911eb4", lwd =2) 
  lines(roc_results_np_15[["DLL1_NP"]], col ="#ffd8b1", lwd =2) 
  lines(roc_results_np_15[["HES1_NP"]], col = "#42d4f4", lwd =2) 
  legend("bottomright", legend = c( expression(italic("NOTCH2")),
                                    expression(italic("CTNNB1")),
                                    expression(italic("DLL1")), 
                                    expression(italic("HES1"))
  ),
  
  col = c("#dcbeff", "#911eb4", "#ffd8b1", "#42d4f4"), lty = 1, 
  cex = 0.6, lwd =3)
}
#plot
roc_plot()

#roc table
coords_results_np_15<- lapply(roc_results_np_15, function(roc_obj) {
  coords(roc_obj, "best", ret = c("threshold", "accuracy", "sensitivity",
                                  "specificity"),
         transpose = FALSE)
})
coords_results_np_15
#create df
results_np_15<- data.frame(
  Predictor = c("NOTCH2_NP","CTNNB1_NP","DLL1_NP","HES1_NP"),
  AUC = AUC_results_np_15,
  do.call(rbind,coords_results_np_15) 
)
results_np_15
#lithuanize it 
colnames(results_np_15) <- c("Biožymuo", "plotas po kreive", "slenkstinė vertė", 
                             "tikslumas", "jautrumas", "specifiškumas")
rownames(results_np_15) <- c("NOCTH2", "CTNNB1", "DLL1", "HES1")
results_np_15$Biožymuo <- c("NOCTH2", "CTNNB1", "DLL1", "HES1")
#nice formating of the Table metrics for ROC OC
gt_table_np_15 <- results_np_15 %>%
  gt() %>%
  tab_header(
    title = "ROC kriterijai", 
    subtitle = "HGOSC atskyrimas nuo kitų KV, NP, 15 atvejų") %>%
  fmt_number(
    columns = everything(),
    decimals = 3
  ) %>%
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_body(columns = vars(Biožymuo))
  )
#show
gt_table_np_15
#URINE roc, 15 CASES ##############################################################
roc_results_urine_15 <- lapply(c("NOTCH2_URINE","CTNNB1_URINE","DLL1_URINE","HES1_URINE"), function(col) {
  roc(response = LIQUID_DF_15$TYPE, predictor = LIQUID_DF_15[[col]])})
names(roc_results_urine_15) <- c("NOTCH2_URINE","CTNNB1_URINE","DLL1_URINE","HES1_URINE")
roc_results_urine_15
#extract the aucs
AUC_results_urine_15 <- sapply(roc_results_urine_15, function(roc_obj) {auc(roc_obj)})
AUC_results_urine_15 #extracted aucs

#roc figure 
roc_plot2 <- function() {
  par(pty = "s") #sets square
  plot.roc(roc_results_urine_15[["NOTCH2_URINE"]], print.auc = F, col = "#dcbeff",
           cex.main=0.8, 
           main ="HGOSC atskyrimas nuo kitų KV, urine, 15 atvejų",
           xlab = "1 - Specifiškumas", 
           ylab = "Jautrumas", 
           legacy.axes = T) #title
  lines(roc_results_urine_15[["CTNNB1_URINE"]], col = "#911eb4", lwd =2) 
  lines(roc_results_urine_15[["DLL1_URINE"]], col ="#ffd8b1", lwd =2) 
  lines(roc_results_urine_15[["HES1_URINE"]], col = "#42d4f4", lwd =2) 
  legend("bottomright", legend = c( expression(italic("NOTCH2")),
                                    expression(italic("CTNNB1")),
                                    expression(italic("DLL1")), 
                                    expression(italic("HES1"))
  ),
  
  col = c("#dcbeff", "#911eb4", "#ffd8b1", "#42d4f4"), lty = 1, 
  cex = 0.6, lwd =3)
}
#plot
roc_plot2()

#roc table
coords_results_URINE_15<- lapply(roc_results_urine_15, function(roc_obj) {
  coords(roc_obj, "best", ret = c("threshold", "accuracy", "sensitivity",
                                  "specificity"),
         transpose = FALSE)
})
coords_results_URINE_15
#create df
results_URINE_15<- data.frame(
  Predictor = c("NOTCH2_URINE","CTNNB1_URINE","DLL1_URINE","HES1_URINE"),
  AUC = AUC_results_urine_15,
  do.call(rbind,coords_results_URINE_15) 
)
results_URINE_15
#lithuanize it 
colnames(results_URINE_15) <- c("Biožymuo", "plotas po kreive", "slenkstinė vertė", 
                                "tikslumas", "jautrumas", "specifiškumas")
rownames(results_URINE_15) <- c("NOCTH2", "CTNNB1", "DLL1", "HES1")
results_URINE_15$Biožymuo <- c("NOCTH2", "CTNNB1", "DLL1", "HES1")
#nice formating of the Table metrics for ROC OC
gt_table_URINE_15 <- results_URINE_15 %>%
  gt() %>%
  tab_header(
    title = "ROC kriterijai", 
    subtitle = "HGOSC atskyrimas nuo kitų KV, URINE, 15 atvejų") %>%
  fmt_number(
    columns = everything(),
    decimals = 3
  ) %>%
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_body(columns = vars(Biožymuo))
  )
#show
gt_table_URINE_15
