#KN-  liquid 2026 04 17, 2026 05 05, 07, 11
#FINAL PLOTS FOR LIQUID PAPER - NP ONLY
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
#read RDS
LIQUID_DF_final <- readRDS("C:/Users/Ieva/rprojects/OTHER DATA/KN_LIQUID/liquid_20260415.RDS")
#leave lavage only
#notch2 the fullest data
LAVAGE_df <- LIQUID_DF_final%>%
  filter(!is.na(NOTCH2_NP)) #103 cases

#fix RSS to RRS
levels(LAVAGE_df$TYPE)[levels(LAVAGE_df$TYPE) == "RSS"] <- "RRS"

#NP GENE COMPARISONS BETWEEN ALL GROUPS################################
##test normalcy NP ###########################
DATA <- c("NOTCH2_NP","CTNNB1_NP","DLL1_NP","HES1_NP" )
results_normalcy <- lapply(DATA, function(v) {
  by(LAVAGE_df[[v]], LAVAGE_df$TYPE, function(x) {
    if (length(na.omit(x)) >= 3)
      shapiro.test(x)$p.value
    else NA
  })
})
names(results_normalcy) <- DATA
results_normalcy #not normal NOTCH2 EC, HES1 benign

##test variance full NP#############################
variance_results <- lapply(DATA, function(v) {
  # Select variable and TYPE
  df <- LAVAGE_df[, c(v, "TYPE")]
  
  # Remove NA
  df <- df[!is.na(df[[v]]), ]
  
  # Only test if at least 2 groups have data
  if(length(unique(df$TYPE)) > 1) {
    test <- leveneTest(df[[v]] ~ df$TYPE, center = median)
    pval <- test[1, "Pr(>F)"]
  } else {
    pval <- NA
  }
  
  data.frame(
    variable = v,
    levene_p = pval,
    equal_variance = ifelse(!is.na(pval) & pval > 0.05, TRUE, FALSE)
  )
})

# Combine into a single table
variance_results <- bind_rows(variance_results)

# View the results
variance_results #all true for NP

##ANOVA for NP##############################################
#ANOVA (for normal data) - CTNNB1
anova_ctnnb1 <- aov(CTNNB1_NP ~ TYPE, data = LAVAGE_df)
summary(anova_ctnnb1) #significant
TukeyHSD(anova_ctnnb1) #significant: 
#ENDOMETRIAL CANCER-BENIGN p.adj = 0.0012000, 
#HGSOC-ENDOMETRIAL CANCER p.adj = 0.0015951, 
#RRS-ENDOMETRIAL CANCER p.adj = 0.0003289

#ANOVA (for normal data) - DLL1 
anova_dll1 <- aov(DLL1_NP ~ TYPE, data = LAVAGE_df)
summary(anova_dll1) #significant
TukeyHSD(anova_dll1) #significant: RRS-ENDOMETRIAL CANCER p.adj = 0.0015694
#significant:
#ENDOMETRIAL CANCER-BENIGN p = 0.0104,
#HGSOC-ENDOMETRIAL CANCER p = 0.0017,         
#RRS-ENDOMETRIAL CANCER p = 0.0092      

#ANOVA (for not normal data) - HES1
kruskal.test(HES1_NP ~ TYPE, data = LAVAGE_df) #significant
pairwise.wilcox.test(LAVAGE_df$HES1_NP, LAVAGE_df$TYPE, p.adjust.method = "bonferroni") 
#siginifcant: HGSOC-ENDOMETRIAL CANCER p = 0.013
dunnTest(HES1_NP ~ TYPE,
         data = LAVAGE_df,
         method = "bonferroni")
#siginifcant: HGSOC-ENDOMETRIAL CANCER p = 0.02355395 


#ANOVA (for not normal data) - NOTCH2
kruskal.test(NOTCH2_NP ~ TYPE, data = LAVAGE_df) #significant 
pairwise.wilcox.test(LAVAGE_df$NOTCH2_NP, LAVAGE_df$TYPE,
                     p.adjust.method = "bonferroni")
#siginifcant: BENIGN -ENDOMETRIAL CANCER p = 0.0104  
#siginifcant: HGSOC-ENDOMETRIAL CANCER p = 0.0017    
#siginifcant: RRS-ENDOMETRIAL CANCER p = 0.0092   
dunnTest(NOTCH2_NP ~ TYPE,
         data = LAVAGE_df,
         method = "bonferroni")
#significant: BENIGN - ENDOMETRIAL CANCER p = 0.027194259
#siginifcant: HGSOC-ENDOMETRIAL CANCER p = 0.001706016
#siginifcant: RRS-ENDOMETRIAL CANCER p = 00.007731407
##boxplot full np ##########################################
table(LAVAGE_df$TYPE, useNA = "a")
#melt table for expression
GroupNP_table <- melt(LAVAGE_df[, c(12,15:18)], id.vars="TYPE",  measure.vars=c("NOTCH2_NP","CTNNB1_NP","DLL1_NP","HES1_NP"))
#fix names
GroupNP_table <- GroupNP_table %>%
  mutate(TYPE = dplyr::recode(TYPE,
                              "ENDOMETRIAL CANCER" = "EC",
                              "OTHER" = "OTHER OC"))
#make p values
each.vs.ref_sig <- tibble::tribble(
  ~group1, ~group2, ~p.adj,   ~y.position, ~variable,
  "EC",   "BENIGN", 0.002, -2, "CTNNB1_NP",
  "EC",   "HGSOC", 0.003, -1, "CTNNB1_NP",
  "EC",   "RRS", 0.001, -1.5, "CTNNB1_NP",
  "EC",   "RRS", 0.002, -6, "DLL1_NP",
  "EC",   "BENIGN", 0.027, -2, "NOTCH2_NP",#dunn
  "EC",   "HGSOC", 0.002, -1, "NOTCH2_NP",#dunn
  "EC",   "RRS", 0.008, -1.5, "NOTCH2_NP",#dunn
  "EC",   "HGSOC", 0.024, -2, "HES1_NP"#dunn
)

TYPE_FULL_plot <- ggplot(GroupNP_table, aes(x=TYPE , y=value, fill = variable)) +
  geom_boxplot( outlier.shape = NA , alpha=0.3, aes(fill = TYPE )) +
  geom_jitter(aes(color = TYPE ), size=1, alpha=0.5) +
  ylab(label = expression("Gene expression, normalized to  " * italic("GAPDH"))) + 
  facet_wrap(.~ variable, nrow = 2, scales = "free",
             labeller = labeller(
               variable = c(
                 "CTNNB1_NP" = "CTNNB1 expression in uterine lavage",
                 "DLL1_NP" = "DLL1 expression in uterine lavage",
                 "HES1_NP" = "HES1 expression in uterine lavage",
                 "NOTCH2_NP" = "NOTCH2 expression in uterine lavage"
               )) ) +
  add_pvalue(each.vs.ref_sig, label = "p.adj") + #pvalue
  theme_minimal()+
  theme(
    strip.text.x = element_text(
      size = 12, face = "bold.italic"
    ),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5))+
  labs(x=NULL)+
  stat_boxplot(geom ='errorbar')+
  scale_y_continuous(labels = function(x) 
    gsub("-", "\u2212", as.character(x)))+
  scale_fill_manual(values = c(
  "EC" = "#E64B35",
  "OTHER OC" = "#4DBBD5",
  "BENIGN" = "#00A087",
  "HGSOC" = "#3C5488",
  "RRS" = "#B7E4C7"
  )) +
  scale_color_manual(values = c(
    "EC" = "#E64B35",
    "OTHER OC" = "#4DBBD5",
    "BENIGN" = "#00A087",
    "HGSOC" = "#3C5488",
    "RRS" = "#B7E4C7"
  ))

TYPE_FULL_plot
#save
ggsave("C:/Users/Ieva/rprojects/outputs_all/LIQUID/NP_FULL_20260619.png",
       plot = TYPE_FULL_plot,
       width = 18,
       height = 16,
       units = "cm",
       dpi = 400)

#RRS+BENING grouped NP #################################################
table(LAVAGE_df$TYPE, useNA = "a") #assess the situation
#create a new grouping where RRS and BENIGN is one
LAVAGE_df <- LAVAGE_df %>%
  mutate(
    TYPE_BENIGN2 = if_else(TYPE %in% c("RRS", "BENIGN"),
                           "BENIGN",
                           TYPE)
  )
table(LAVAGE_df$TYPE_BENIGN2, useNA = "a") #now 31 benign

##test normalcy of my variables, RRS+BENIGN#####################################
results_normalcy2 <- lapply(DATA, function(v) {
  by(LAVAGE_df[[v]], LAVAGE_df$TYPE_BENIGN2, function(x) {
    if (length(na.omit(x)) >= 3)
      shapiro.test(x)$p.value
    else NA
  })
})
names(results_normalcy2) <- DATA
results_normalcy2 #Not normal NP NOTCH2 endometrial cancer p = 0.01915093

##test variance RRS+BENIGN###############################
# Loop through all variables in DATA_names
variance_results2 <- lapply(DATA, function(v) {
  # Select variable and TYPE
  df <- LAVAGE_df[, c(v, "TYPE_BENIGN2")]
  
  # Remove NA
  df <- df[!is.na(df[[v]]), ]
  
  # Only test if at least 2 groups have data
  if(length(unique(df$TYPE_BENIGN2)) > 1) {
    test <- leveneTest(df[[v]] ~ df$TYPE_BENIGN2, center = median)
    pval <- test[1, "Pr(>F)"]
  } else {
    pval <- NA
  }
  
  data.frame(
    variable = v,
    levene_p = pval,
    equal_variance = ifelse(!is.na(pval) & pval > 0.05, TRUE, FALSE)
  )
})

# Combine into a single table
variance_results2 <- bind_rows(variance_results2)

# View the results
variance_results2 #all equal

##ANOVA, RRS+BENIGN###########################
#ANOVA (for normal data) - CTNNB1, DLL1 and HES1
# CTNNB1
anova_ctnnb1_BENIGN2 <- aov(CTNNB1_NP ~ TYPE_BENIGN2, data = LAVAGE_df)
summary(anova_ctnnb1_BENIGN2) #significant
TukeyHSD(anova_ctnnb1_BENIGN2) 
#significant: ENDOMETRIAL CANCER-BENIGN 0.0000576, HGSOC-ENDOMETRIAL CANCER 0.0009424
#DLL1
anova_dll1_BENIGN2 <- aov(DLL1_NP ~ TYPE_BENIGN2, data = LAVAGE_df)
summary(anova_dll1_BENIGN2) #significant
TukeyHSD(anova_dll1_BENIGN2) 
#significant: BENING-ENDOMETRIAL CANCER 0.0027565, HGSOC-ENDOMETRIAL CANCER 0.0471971
#HES1
anova_hes1_BENIGN2 <- aov(HES1_NP ~ TYPE_BENIGN2, data = LAVAGE_df)
summary(anova_hes1_BENIGN2) #not significant
TukeyHSD(anova_hes1_BENIGN2) 
#significant: HGSOC-BENIGN 0.0340468, HGSOC-ENDOMETRIAL CANCER 0.0142440
#ANOVA (for not normal data) - NOTCH2 
#NOTCH2
kruskal.test(NOTCH2_NP ~ TYPE_BENIGN2, data = LAVAGE_df) #significant
pairwise.wilcox.test(LAVAGE_df$NOTCH2_NP, LAVAGE_df$TYPE_BENIGN2, p.adjust.method = "bonferroni")
#significant: BENING-ENDOMETRIAL CANCER 0.00095 ,
##HGSOC-ENDOMETRIAL CANCER 0.00099,
##OTHER-ENDOMETRIAL CANCER 0.04504   
dunnTest(NOTCH2_NP ~ TYPE_BENIGN2,
         data = LAVAGE_df,
         method = "bonferroni")
#significant: BENING-ENDOMETRIAL CANCER 0.002364269,
##HGSOC-ENDOMETRIAL CANCER 0.001023609,
##OTHER-ENDOMETRIAL CANCER   0.154074655 #not significant anymore
##boxplot full np, RRS+BENIGN ###########################
#melt table for expression
GroupNP_table_BENIGN2 <- melt(LAVAGE_df[, c(38,15:18)], id.vars="TYPE_BENIGN2",
                              measure.vars=c("NOTCH2_NP","CTNNB1_NP","DLL1_NP","HES1_NP"))
#fix names
GroupNP_table_BENIGN2 <- GroupNP_table_BENIGN2 %>%
  mutate(TYPE_BENIGN2 = dplyr::recode(TYPE_BENIGN2,
                              "ENDOMETRIAL CANCER" = "EC",
                              "OTHER" = "OTHER OC",
                              "BENIGN" = "BENIGN + RRS"))
#make p values RRS+benign
each.vs.ref_sig_BENIGN2 <- tibble::tribble(
  ~group1, ~group2, ~p.adj,   ~y.position, ~variable,
  "EC",   "BENIGN + RRS", 0.001, -2, "CTNNB1_NP",#
  "EC",   "HGSOC", 0.001, -1, "CTNNB1_NP",#
  "EC",   "BENIGN + RRS", 0.002 , -2, "NOTCH2_NP",#dunn
  "EC",   "HGSOC", 0.001, -1, "NOTCH2_NP",#dunn
#  "EC",   "OTHER OC", 0.04504, -1.5, "NOTCH2_NP",#wilcox, dunn not significant after adj
  "EC",   "HGSOC", 0.014, -2, "HES1_NP",#
  "BENIGN + RRS",   "HGSOC", 0.034, -1.5, "HES1_NP",#
  "EC",   "HGSOC", 0.047, -6, "DLL1_NP",#
  "EC",   "BENIGN + RRS", 0.003, -7, "DLL1_NP"#
)

#make figure RRS+benign

TYPE_RRS_BENIGN_plot <- ggplot(GroupNP_table_BENIGN2, aes(x=TYPE_BENIGN2 , y=value, fill = variable)) +
  geom_boxplot( outlier.shape = NA , alpha=0.3, aes(fill = TYPE_BENIGN2 )) +
  geom_jitter(aes(color = TYPE_BENIGN2 ), size=1, alpha=0.5) +
  ylab(label = expression("Gene expression, normalized to  " * italic("GAPDH"))) + 
  facet_wrap(.~ variable, nrow = 2, scales = "free",
             labeller = labeller(
               variable = c(
                 "CTNNB1_NP" = "CTNNB1 expression in uterine lavage",
                 "DLL1_NP" = "DLL1 expression in uterine lavage",
                 "HES1_NP" = "HES1 expression in uterine lavage",
                 "NOTCH2_NP" = "NOTCH2 expression in uterine lavage"
               )) ) +
  add_pvalue(each.vs.ref_sig_BENIGN2, label = "p.adj") + #pvalue
  theme_minimal()+
  theme(
    strip.text.x = element_text(
      size = 12, face = "bold.italic"
    ),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5))+
  labs(x=NULL)+
  stat_boxplot(geom ='errorbar')+
  scale_y_continuous(labels = function(x) 
    gsub("-", "\u2212", as.character(x)))+
  scale_fill_manual(values = c(
    "EC" = "#E64B35",
    "OTHER OC" = "#4DBBD5",
    "BENIGN + RRS" = "#00A087",
    "HGSOC" = "#3C5488"
  )) +
  scale_color_manual(values = c(
    "EC" = "#E64B35",
    "OTHER OC" = "#4DBBD5",
    "BENIGN + RRS" = "#00A087",
    "HGSOC" = "#3C5488"
  ))

TYPE_RRS_BENIGN_plot
#save RRS+benign
ggsave("C:/Users/Ieva/rprojects/outputs_all/LIQUID/NP_RRS_BEN_20260619.png",
       plot = TYPE_RRS_BENIGN_plot,
       width = 18,
       height = 16,
       units = "cm",
       dpi = 400)

#HGSOC+OTHER grouped NP #################################################
table(LAVAGE_df$TYPE_BENIGN2, useNA = "a") #assess the situation
#create a new grouping where RRS and BENIGN is one
LAVAGE_df <- LAVAGE_df %>%
  mutate(
    TYPE_BENIGN3 = if_else(TYPE_BENIGN2 %in% c("HGSOC", "OTHER"),
                           "OC",
                           TYPE_BENIGN2)
  )
table(LAVAGE_df$TYPE_BENIGN3, useNA = "a") #now 61 OC
##test normalcy of my variables, OC#####################################
results_normalcy3 <- lapply(DATA, function(v) {
  by(LAVAGE_df[[v]], LAVAGE_df$TYPE_BENIGN3, function(x) {
    if (length(na.omit(x)) >= 3)
      shapiro.test(x)$p.value
    else NA
  })
})
names(results_normalcy3) <- DATA
results_normalcy3 #Not normal NP CTNNB1 OC p = 0.02547908, NOTCH2 EC p = 0.01915093
##test variance OC###############################
# Loop through all variables in DATA_names
variance_results3 <- lapply(DATA, function(v) {
  # Select variable and TYPE
  df <- LAVAGE_df[, c(v, "TYPE_BENIGN3")]
  
  # Remove NA
  df <- df[!is.na(df[[v]]), ]
  
  # Only test if at least 2 groups have data
  if(length(unique(df$TYPE_BENIGN3)) > 1) {
    test <- leveneTest(df[[v]] ~ df$TYPE_BENIGN3, center = median)
    pval <- test[1, "Pr(>F)"]
  } else {
    pval <- NA
  }
  
  data.frame(
    variable = v,
    levene_p = pval,
    equal_variance = ifelse(!is.na(pval) & pval > 0.05, TRUE, FALSE)
  )
})

# Combine into a single table
variance_results3 <- bind_rows(variance_results3)

# View the results
variance_results3 #ctnnb1 not equal (but also not normal)
##ANOVA, HGSOC+OTHERN###########################
#ANOVA (for normal data) - DLL1 and HES1
# DLL1
anova_DLL1_BENIGN3<- aov(DLL1_NP ~ TYPE_BENIGN3, data = LAVAGE_df)
summary(anova_DLL1_BENIGN3) #significant
TukeyHSD(anova_DLL1_BENIGN3) 
#significant: ENDOMETRIAL CANCER-BENIGN 0.0013537;
#OC-ENDOMETRIAL CANCER 0.0283727
#HES1
anova_hes1_BENIGN3 <- aov(HES1_NP ~ TYPE_BENIGN3, data = LAVAGE_df)
summary(anova_hes1_BENIGN3) #signficant
TukeyHSD(anova_hes1_BENIGN3) 
#significant: oc-endometrial 0.0072412, OC - benign 0.0163029
#ANOVA (for not normal data) - NOTCH2 / ctnnb1
#NOTCH2
kruskal.test(NOTCH2_NP ~ TYPE_BENIGN3, data = LAVAGE_df) #significant
pairwise.wilcox.test(LAVAGE_df$NOTCH2_NP, LAVAGE_df$TYPE_BENIGN3, p.adjust.method = "bonferroni")
#significant:ENDOMETRIAL CANCER-BENIGN  0.00047 ;
# ENDOMETRIAL CANCER-OC 0.00114  
dunnTest(NOTCH2_NP ~ TYPE_BENIGN3,
         data = LAVAGE_df,
         method = "bonferroni")
#significant:ENDOMETRIAL CANCER-BENIGN ENDOMETRIAL CANCER 0.0011821343 ;
# ENDOMETRIAL CANCER-OC 0.0009872979  

#CTNNB1
kruskal.test(CTNNB1_NP ~ TYPE_BENIGN3, data = LAVAGE_df) #significant
pairwise.wilcox.test(LAVAGE_df$CTNNB1_NP, LAVAGE_df$TYPE_BENIGN3, p.adjust.method = "bonferroni")
#significant::ENDOMETRIAL CANCER-BENIGN ENDOMETRIAL CANCER 0.0062  ;
# ENDOMETRIAL CANCER-OC 0.0245   
dunnTest(CTNNB1_NP ~ TYPE_BENIGN3,
         data = LAVAGE_df,
         method = "bonferroni")
#significant::ENDOMETRIAL CANCER-BENIGN ENDOMETRIAL CANCER 0.002021517  ;
# ENDOMETRIAL CANCER-OC 0.043300035   

##boxplot full np, HGSOC+OTHERS#######################
#melt table for expression
GroupNP_table_BENIGN3 <- melt(LAVAGE_df[, c(39,15:18)], id.vars="TYPE_BENIGN3",
                              measure.vars=c("NOTCH2_NP","CTNNB1_NP","DLL1_NP","HES1_NP"))
#fix names
GroupNP_table_BENIGN3 <- GroupNP_table_BENIGN3 %>%
  mutate(TYPE_BENIGN3 = dplyr::recode(TYPE_BENIGN3,
                                      "ENDOMETRIAL CANCER" = "EC",
                                      "BENIGN" = "BENIGN + RRS"))
#make p values RRS+benign
each.vs.ref_sig_BENIGN3 <- tibble::tribble(
  ~group1, ~group2, ~p.adj,   ~y.position, ~variable,
  "EC",   "BENIGN + RRS", 0.007, -5, "DLL1_NP",#
  "EC",   "OC", 0.028, -4, "DLL1_NP",#
  "EC",   "OC", 0.007, -2, "HES1_NP",#
  "BENIGN + RRS",   "OC", 0.016, -1, "HES1_NP",#
  "EC",   "OC", 0.001, -2, "NOTCH2_NP",#dunn
  "EC",   "BENIGN + RRS", 0.001 , -3, "NOTCH2_NP",#dunn
  "EC",   "BENIGN + RRS", 0.002 , -2, "CTNNB1_NP",#dunn
  "EC",   "OC", 0.043, -1, "CTNNB1_NP"#dunn
)
#make figure RRS+benign

TYPE_RRS_BENIGN_plot3 <- ggplot(GroupNP_table_BENIGN3, aes(x=TYPE_BENIGN3 , y=value, fill = variable)) +
  geom_boxplot( outlier.shape = NA , alpha=0.3, aes(fill = TYPE_BENIGN3 )) +
  geom_jitter(aes(color = TYPE_BENIGN3 ), size=1, alpha=0.5) +
  ylab(label = expression("Gene expression, normalized to  " * italic("GAPDH"))) + 
  facet_wrap(.~ variable, nrow = 2, scales = "free",
             labeller = labeller(
               variable = c(
                 "CTNNB1_NP" = "CTNNB1 expression in uterine lavage",
                 "DLL1_NP" = "DLL1 expression in uterine lavage",
                 "HES1_NP" = "HES1 expression in uterine lavage",
                 "NOTCH2_NP" = "NOTCH2 expression in uterine lavage"
               )) ) +
  add_pvalue(each.vs.ref_sig_BENIGN3, label = "p.adj") + #pvalue
  theme_minimal()+
  theme(
    strip.text.x = element_text(
      size = 12, face = "bold.italic"
    ),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5))+
  labs(x=NULL)+
  stat_boxplot(geom ='errorbar')+
  scale_y_continuous(labels = function(x) 
    gsub("-", "\u2212", as.character(x)))+
  scale_fill_manual(values = c(
    "EC" = "#E64B35",
    "BENIGN + RRS" = "#00A087",
    "OC" = "#3C5488"
  )) +
  scale_color_manual(values = c(
    "EC" = "#E64B35",
    "BENIGN + RRS" = "#00A087",
    "OC" = "#3C5488"
  ))

TYPE_RRS_BENIGN_plot3
#save RRS+benign
ggsave("C:/Users/Ieva/rprojects/outputs_all/LIQUID/NP_RRS_BEN_HGSOC_OTHER_20260619.png",
       plot = TYPE_RRS_BENIGN_plot3,
       width = 18,
       height = 16,
       units = "cm",
       dpi = 400)
#ROC########################################################
##ROC HGSOC VS BENIGN ########################################
#make hgsoc VS BENIGN DF
HGSOC_BENIGN_DF<- LAVAGE_df %>%
  filter(TYPE%in% c("BENIGN", "HGSOC"))%>% #filter for right samples
  dplyr::select("TYPE", "NOTCH2_NP","CTNNB1_NP","DLL1_NP","HES1_NP"
  )
HGSOC_BENIGN_DF$TYPE <- droplevels(HGSOC_BENIGN_DF$TYPE) #drop unused
HGSOC_BENIGN_DF$TYPE  <- relevel(HGSOC_BENIGN_DF$TYPE , ref = "BENIGN") #set control
table(HGSOC_BENIGN_DF$TYPE, useNA = "a")
roc_results_np_HGOSC_benign<- lapply(c("NOTCH2_NP","CTNNB1_NP","DLL1_NP","HES1_NP"), function(col) {
  roc(response = HGSOC_BENIGN_DF$TYPE, predictor = HGSOC_BENIGN_DF[[col]])})
names(roc_results_np_HGOSC_benign) <- c("NOTCH2_NP","CTNNB1_NP","DLL1_NP","HES1_NP")
roc_results_np_HGOSC_benign
#extract the aucs
auc_values_np_HGSOC_benign <- sapply(roc_results_np_HGOSC_benign, function(roc_obj) {auc(roc_obj)})
auc_values_np_HGSOC_benign #extracted aucs
#roc figure 
roc_plot3 <- function() {
  par(pty = "s") #sets square
  plot.roc(roc_results_np_HGOSC_benign[["NOTCH2_NP"]], print.auc = F, col = "#dcbeff",
           cex.main=0.8, 
           main ="Uterine lavage biomarkers HGSOC vs Benign samples",
           # xlab = "1 - Specifiškumas", #lithuanian version
           # ylab = "Jautrumas", 
           legacy.axes = T) #title
  lines(roc_results_np_HGOSC_benign[["CTNNB1_NP"]], col = "#911eb4", lwd =2) 
  lines(roc_results_np_HGOSC_benign[["DLL1_NP"]], col ="#ffd8b1", lwd =2) 
  lines(roc_results_np_HGOSC_benign[["HES1_NP"]], col = "#42d4f4", lwd =2) 
  legend("bottomright", legend = c( expression(italic("NOTCH2")),
                                    expression(italic("CTNNB1")),
                                    expression(italic("DLL1")), 
                                    expression(italic("HES1"))
  ),
  
  col = c("#dcbeff", "#911eb4", "#ffd8b1", "#42d4f4"), lty = 1, 
  cex = 0.7, lwd =3)
}
#plot
roc_plot3()
#roc table
coords_results_np_HGSOC_benign<- lapply(roc_results_np_HGOSC_benign, function(roc_obj) {
  coords(roc_obj, "best", ret = c("threshold", "accuracy", "sensitivity",
                                  "specificity"),
         transpose = FALSE)
})
coords_results_np_HGSOC_benign
#create df
results_np_HGSOC_benign<- data.frame(
  Predictor = c("NOTCH2_NP","CTNNB1_NP","DLL1_NP","HES1_NP"),
  AUC = auc_values_np_HGSOC_benign,
  do.call(rbind,coords_results_np_HGSOC_benign) 
)
results_np_HGSOC_benign
#change names to remove _np
results_np_HGSOC_benign$Predictor <- c("NOTCH2", #fixed
                                       "CTNNB1", "DLL1", "HES1")
#nice formating of the Table metrics for ROC OC
gt_table_np_HGSOC_ben <- results_np_HGSOC_benign %>%
  gt() %>%
  tab_header(
    title = "ROC measures for uterine lavage biomarkers", 
    subtitle = "Benign vs HGSOC") %>%
  fmt_number(
    columns = everything(),
    decimals = 3
  ) %>%
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_body(columns = vars(Predictor))
  )
#show
gt_table_np_HGSOC_ben

#save plot
# Save the plot as a PNG file
png("C:/Users/Ieva/rprojects/outputs_all/LIQUID/NP_HGSOC_BENIGN_ROC_20260618.png",
    width = 10, height = 10, res = 300, units = "cm")
roc_plot3()
mtext(
  "A",
  side = 3,
  adj = -0.2,
  line = 1,
  cex = 1.2,
  font = 2
)
dev.off()
#there is no other convenient way to save gt outputs
gtsave(gt_table_np_HGSOC_ben,vwidth = 10000,   
       filename = "C:/Users/Ieva/rprojects/outputs_all/LIQUID/NP_HGSOC_BENIGN_ROCtable_20260618.png")
#import images
roc_image2   <- image_read("C:/Users/Ieva/rprojects/outputs_all/LIQUID/NP_HGSOC_BENIGN_ROC_20260618.png")
table_image2 <- image_read("C:/Users/Ieva/rprojects/outputs_all/LIQUID/NP_HGSOC_BENIGN_ROCtable_20260618.png")

# resize table to match ROC image width
table_image2 <- image_resize(table_image2,
                             paste0(image_info(roc_image2)$width, "x"))

# OR resize ROC image to match table width
# roc_image2 <- image_resize(roc_image2,
#                            paste0(image_info(table_image2)$width, "x"))

# combine vertically
combined <- image_append(c(roc_image2, table_image2), stack = TRUE)

# save
image_write(combined,
            "C:/Users/Ieva/rprojects/outputs_all/LIQUID/NP_HGSOC_BENIGN_ROCcombined_20260618.png")
##ROC EC VS BENIGN ########################################
#make EC VS BENIGN DF
EC_BENIGN_DF<- LAVAGE_df %>%
  filter(TYPE%in% c("BENIGN", "ENDOMETRIAL CANCER"))%>% #filter for right samples
  dplyr::select("TYPE", "NOTCH2_NP","CTNNB1_NP","DLL1_NP","HES1_NP"
  )
EC_BENIGN_DF$TYPE <- droplevels(EC_BENIGN_DF$TYPE) #drop unused
EC_BENIGN_DF$TYPE  <- relevel(EC_BENIGN_DF$TYPE , ref = "BENIGN") #set control
table(EC_BENIGN_DF$TYPE, useNA = "a")
roc_results_np_EC_benign<- lapply(c("NOTCH2_NP","CTNNB1_NP","DLL1_NP","HES1_NP"), function(col) {
  roc(response = EC_BENIGN_DF$TYPE, predictor = EC_BENIGN_DF[[col]])})
names(roc_results_np_EC_benign) <- c("NOTCH2_NP","CTNNB1_NP","DLL1_NP","HES1_NP")
roc_results_np_EC_benign
#extract the aucs
auc_values_np_EC_benign <- sapply(roc_results_np_EC_benign, function(roc_obj)
{auc(roc_obj)})
auc_values_np_EC_benign #extracted aucs
#roc figure 
roc_plot2 <- function() {
  par(pty = "s") #sets square
  plot.roc(roc_results_np_EC_benign[["NOTCH2_NP"]], print.auc = F, col = "#dcbeff",
           cex.main=0.8, 
           main ="Uterine lavage biomarkers EC vs Benign samples",
           # xlab = "1 - Specifiškumas", #lithuanian version
           # ylab = "Jautrumas", 
           legacy.axes = T) #title
  lines(roc_results_np_EC_benign[["CTNNB1_NP"]], col = "#911eb4", lwd =2) 
  lines(roc_results_np_EC_benign[["DLL1_NP"]], col ="#ffd8b1", lwd =2) 
  lines(roc_results_np_EC_benign[["HES1_NP"]], col = "#42d4f4", lwd =2) 
  legend("bottomright", legend = c( expression(italic("NOTCH2")),
                                    expression(italic("CTNNB1")),
                                    expression(italic("DLL1")), 
                                    expression(italic("HES1"))
  ),
  
  col = c("#dcbeff", "#911eb4", "#ffd8b1", "#42d4f4"), lty = 1, 
  cex = 0.7, lwd =3)
}
#plot
roc_plot2()
#roc table
coords_results_np_EC_benign<- lapply(roc_results_np_EC_benign, function(roc_obj) {
  coords(roc_obj, "best", ret = c("threshold", "accuracy", "sensitivity",
                                  "specificity"),
         transpose = FALSE)
})
coords_results_np_EC_benign
#create df
results_np_EC_benign<- data.frame(
  Predictor = c("NOTCH2_NP","CTNNB1_NP","DLL1_NP","HES1_NP"),
  AUC = auc_values_np_EC_benign,
  do.call(rbind,coords_results_np_EC_benign) 
)
results_np_EC_benign
#change names to remove _np
results_np_EC_benign$Predictor <- c("NOTCH2", #fixed
                                    "CTNNB1", "DLL1", "HES1")
#nice formating of the Table metrics for ROC OC
gt_table_np_EC_ben <- results_np_EC_benign %>%
  gt() %>%
  tab_header(
    title = "ROC measures for uterine lavage biomarkers", 
    subtitle = "Benign vs endometrial cancer") %>%
  fmt_number(
    columns = everything(),
    decimals = 3
  ) %>%
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_body(columns = vars(Predictor))
  )
#show
gt_table_np_EC_ben

#save plot
# Save the plot as a PNG file
png("C:/Users/Ieva/rprojects/outputs_all/LIQUID/NP_EC_BENIGN_ROC_20260618.png",
    width = 10, height = 10, res = 300, units = "cm")
roc_plot2()
mtext(
  "B",
  side = 3,
  adj = -0.2,
  line = 1,
  cex = 1.2,
  font = 2
)
dev.off()
#there is no other convieneat way to save gt outputs
gtsave(gt_table_np_EC_ben,vwidth = 10000,   
       filename = "C:/Users/Ieva/rprojects/outputs_all/LIQUID/NP_EC_BENIGN_ROCtable_20260618.png")
#import images
roc_image3  <- image_read("C:/Users/Ieva/rprojects/outputs_all/LIQUID/NP_EC_BENIGN_ROC_20260618.png")
table_image3 <- image_read("C:/Users/Ieva/rprojects/outputs_all/LIQUID/NP_EC_BENIGN_ROCtable_20260618.png")

# resize table to match ROC image width
table_image3 <- image_resize(table_image3,
                             paste0(image_info(roc_image3)$width, "x"))
# combine vertically
combined3 <- image_append(c(roc_image3, table_image3), stack = TRUE)

# save
image_write(combined3,
            "C:/Users/Ieva/rprojects/outputs_all/LIQUID/NP_EC_BENIGN_ROCcombined_20260618.png")
##ROC HGSOC VS BENIGN RRS ##########################
#make hgsoc VS BENIGN and RRS DF
HGSOC_BENIGN_RRS_DF<- LAVAGE_df %>%
  filter(TYPE_BENIGN2%in% c("BENIGN", "HGSOC"))%>% #filter for right samples
  dplyr::select("TYPE_BENIGN2", "NOTCH2_NP","CTNNB1_NP","DLL1_NP","HES1_NP"
  )
HGSOC_BENIGN_RRS_DF$TYPE_BENIGN2 <- factor(HGSOC_BENIGN_RRS_DF$TYPE_BENIGN2)
HGSOC_BENIGN_RRS_DF$TYPE_BENIGN2 <- droplevels(HGSOC_BENIGN_RRS_DF$TYPE_BENIGN2) #drop unused
HGSOC_BENIGN_RRS_DF$TYPE_BENIGN2  <- relevel(HGSOC_BENIGN_RRS_DF$TYPE_BENIGN2 , ref = "BENIGN") #set control
table(HGSOC_BENIGN_RRS_DF$TYPE_BENIGN2, useNA = "a")
roc_results_np_HGOSC_benignRRS<- lapply(c("NOTCH2_NP","CTNNB1_NP","DLL1_NP","HES1_NP"), function(col) {
  roc(response = HGSOC_BENIGN_RRS_DF$TYPE_BENIGN2, predictor = HGSOC_BENIGN_RRS_DF[[col]])})
names(roc_results_np_HGOSC_benignRRS) <- c("NOTCH2_NP","CTNNB1_NP","DLL1_NP","HES1_NP")
roc_results_np_HGOSC_benignRRS
#extract the aucs
auc_values_np_HGSOC_benignRRS <- sapply(roc_results_np_HGOSC_benignRRS, function(roc_obj) {auc(roc_obj)})
auc_values_np_HGSOC_benignRRS #extracted aucs
#roc figure 
roc_plot3RRS <- function() {
  par(pty = "s") #sets square
  plot.roc(roc_results_np_HGOSC_benignRRS[["NOTCH2_NP"]], print.auc = F, col = "#dcbeff",
           cex.main=0.8, 
           main ="Uterine lavage biomarkers HGSOC vs Non-cancer samples",
           #xlab = "1 - Specifiškumas", 
           #ylab = "Jautrumas", 
           legacy.axes = T) #title
  lines(roc_results_np_HGOSC_benignRRS[["CTNNB1_NP"]], col = "#911eb4", lwd =2) 
  lines(roc_results_np_HGOSC_benignRRS[["DLL1_NP"]], col ="#ffd8b1", lwd =2) 
  lines(roc_results_np_HGOSC_benignRRS[["HES1_NP"]], col = "#42d4f4", lwd =2) 
  legend("bottomright", legend = c( expression(italic("NOTCH2")),
                                    expression(italic("CTNNB1")),
                                    expression(italic("DLL1")), 
                                    expression(italic("HES1"))
  ),
  
  col = c("#dcbeff", "#911eb4", "#ffd8b1", "#42d4f4"), lty = 1, 
  cex = 0.7, lwd =3)
}
#plot
roc_plot3RRS()
#roc table
coords_results_np_HGSOC_benignRRS<- lapply(roc_results_np_HGOSC_benignRRS, function(roc_obj) {
  coords(roc_obj, "best", ret = c("threshold", "accuracy", "sensitivity",
                                  "specificity"),
         transpose = FALSE)
})
coords_results_np_HGSOC_benignRRS
#create df
results_np_HGSOC_benignRRS<- data.frame(
  Predictor = c("NOTCH2_NP","CTNNB1_NP","DLL1_NP","HES1_NP"),
  AUC = auc_values_np_HGSOC_benignRRS,
  do.call(rbind,coords_results_np_HGSOC_benignRRS) 
)
results_np_HGSOC_benignRRS
#lithuanize it 
results_np_HGSOC_benignRRS$Predictor <- c("NOTCH2",#fixed
                                          "CTNNB1", "DLL1", "HES1")
#nice formating of the Table metrics for ROC OC
gt_table_np_HGSOC_benRRS <- results_np_HGSOC_benignRRS %>%
  gt() %>%
  tab_header(
    title = "ROC measures for uterine lavage biomarkers", 
    subtitle = "Non-cancer vs HGSOC") %>%
  fmt_number(
    columns = everything(),
    decimals = 3
  ) %>%
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_body(columns = vars(Predictor))
  )
#show
gt_table_np_HGSOC_benRRS
#save plot
# Save the plot as a PNG file
png("C:/Users/Ieva/rprojects/outputs_all/LIQUID/NP_HGSOC_BENIGNRRS_ROC_20260618.png",
    width = 11, height = 11, res = 300, units = "cm")
roc_plot3RRS()
mtext(
  "C",
  side = 3,
  adj = -0.2,
  line = 1,
  cex = 1.2,
  font = 2
)
dev.off()
#there is no other convieneat way to save gt outputs
gtsave(gt_table_np_HGSOC_benRRS,vwidth = 10000,   
       filename = "C:/Users/Ieva/rprojects/outputs_all/LIQUID/NP_HGSOC_BENIGNRRS_ROCtable_20260618.png")
#import images
roc_image3RRS  <- image_read("C:/Users/Ieva/rprojects/outputs_all/LIQUID/NP_HGSOC_BENIGNRRS_ROC_20260618.png")
table_image3RRS <- image_read("C:/Users/Ieva/rprojects/outputs_all/LIQUID/NP_HGSOC_BENIGNRRS_ROCtable_20260618.png")

# resize table to match ROC image width
table_image3RRS <- image_resize(table_image3RRS,
                                paste0(image_info(roc_image3RRS)$width, "x"))
# combine vertically
combined3RRS <- image_append(c(roc_image3RRS, table_image3RRS), stack = TRUE)

# save
image_write(combined3RRS,
            "C:/Users/Ieva/rprojects/outputs_all/LIQUID/NP_HGSOC_BENIGNRRS_ROCcombined_20260618.png")
##ROC ENDOMETRIAL VS BENIGN / RRS########################################
#make ENDOMETRIAL VS BENIGN /RRS DF
ENDO_BENIGN_DF2<- LAVAGE_df %>%
  filter(TYPE_BENIGN2%in% c("BENIGN", "ENDOMETRIAL CANCER"))%>% #filter for right samples
  dplyr::select("TYPE_BENIGN2", "NOTCH2_NP","CTNNB1_NP","DLL1_NP","HES1_NP"
  )
ENDO_BENIGN_DF2$TYPE_BENIGN2 <- factor(ENDO_BENIGN_DF2$TYPE_BENIGN2)
ENDO_BENIGN_DF2$TYPE_BENIGN2 <- droplevels(ENDO_BENIGN_DF2$TYPE_BENIGN2) #drop unused
ENDO_BENIGN_DF2$TYPE_BENIGN2  <- relevel(ENDO_BENIGN_DF2$TYPE_BENIGN2 , ref = "BENIGN") #set control
table(ENDO_BENIGN_DF2$TYPE_BENIGN2 , useNA = "a")
roc_results_np_endo_benign2<- lapply(c("NOTCH2_NP","CTNNB1_NP","DLL1_NP","HES1_NP"), function(col) {
  roc(response = ENDO_BENIGN_DF2$TYPE_BENIGN2, predictor = ENDO_BENIGN_DF2[[col]])})
names(roc_results_np_endo_benign2) <- c("NOTCH2_NP","CTNNB1_NP","DLL1_NP","HES1_NP")
roc_results_np_endo_benign2
#extract the aucs
auc_values_np_endo_benign2 <- sapply(roc_results_np_endo_benign2, function(roc_obj) {auc(roc_obj)})
auc_values_np_endo_benign2 #extracted aucs
#roc figure 
roc_plot_ENDO <- function() {
  par(pty = "s") #sets square
  plot.roc(roc_results_np_endo_benign2[["NOTCH2_NP"]], print.auc = F, col = "#dcbeff",
           cex.main=0.8, 
           main ="Uterine lavage biomarkers EC vs Non-cancer samples",
           #xlab = "1 - Specifiškumas", 
           #ylab = "Jautrumas", 
           legacy.axes = T) #title
  lines(roc_results_np_endo_benign2[["CTNNB1_NP"]], col = "#911eb4", lwd =2) 
  lines(roc_results_np_endo_benign2[["DLL1_NP"]], col ="#ffd8b1", lwd =2) 
  lines(roc_results_np_endo_benign2[["HES1_NP"]], col = "#42d4f4", lwd =2) 
  legend("bottomright", legend = c( expression(italic("NOTCH2")),
                                    expression(italic("CTNNB1")),
                                    expression(italic("DLL1")), 
                                    expression(italic("HES1"))
  ),
  
  col = c("#dcbeff", "#911eb4", "#ffd8b1", "#42d4f4"), lty = 1, 
  cex = 0.7, lwd =3)
}
#plot
roc_plot_ENDO()
#roc table
coords_results_np_endo_benign2<- lapply(roc_results_np_endo_benign2, function(roc_obj) {
  coords(roc_obj, "best", ret = c("threshold", "accuracy", "sensitivity",
                                  "specificity"),
         transpose = FALSE)
})
coords_results_np_endo_benign2
#create df
results_np_endo_benign2<- data.frame(
  Predictor = c("NOTCH2_NP","CTNNB1_NP","DLL1_NP","HES1_NP"),
  AUC = auc_values_np_endo_benign2,
  do.call(rbind,coords_results_np_endo_benign2) 
)
results_np_endo_benign2
#change names
results_np_endo_benign2$Predictor <- c("NOTCH2",#fixed
                                       "CTNNB1", "DLL1", "HES1")
#nice formating of the Table metrics for ROC OC
gt_table_np_endo_ben2 <- results_np_endo_benign2 %>%
  gt() %>%
  tab_header(
    title = "ROC measures for uterine lavage biomarkers", 
    subtitle = "Non-cancer vs endometrial cancer") %>%
  fmt_number(
    columns = everything(),
    decimals = 3
  ) %>%
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_body(columns = vars(Predictor))
  )
#show
gt_table_np_endo_ben2
# Save the plot as a PNG file
png("C:/Users/Ieva/rprojects/outputs_all/LIQUID/NP_EC_BENIGNRRS_ROC_20260618.png",
    width = 11, height = 11, res = 300, units = "cm")
roc_plot_ENDO()
mtext(
  "D",
  side = 3,
  adj = -0.2,
  line = 1,
  cex = 1.2,
  font = 2
)
dev.off()
#there is no other convieneat way to save gt outputs
gtsave(gt_table_np_endo_ben2,vwidth = 10000,   
       filename = "C:/Users/Ieva/rprojects/outputs_all/LIQUID/NP_EC_BENIGNRRS_ROCtable_20260618.png")
#import images
roc_image4RRS  <- image_read("C:/Users/Ieva/rprojects/outputs_all/LIQUID/NP_EC_BENIGNRRS_ROC_20260618.png")
table_image4RRS <- image_read("C:/Users/Ieva/rprojects/outputs_all/LIQUID/NP_EC_BENIGNRRS_ROCtable_20260618.png")

# resize table to match ROC image width
table_image4RRS <- image_resize(table_image4RRS,
                                paste0(image_info(roc_image4RRS)$width, "x"))
# combine vertically
combined4RRS <- image_append(c(roc_image4RRS, table_image4RRS), stack = TRUE)

# save
image_write(combined4RRS,
            "C:/Users/Ieva/rprojects/outputs_all/LIQUID/NP_EC_BENIGNRRS_ROCcombined_20260618.png")
##ROC OC (HGSOC+OTHER) VS EC ########################################
#make OC VS EC DF
OCEC_DF<- LAVAGE_df %>%
  filter(TYPE_BENIGN3%in% c("OC", "ENDOMETRIAL CANCER"))%>% #filter for right samples
  dplyr::select("TYPE_BENIGN3", "NOTCH2_NP","CTNNB1_NP","DLL1_NP","HES1_NP"
  )
OCEC_DF$TYPE_BENIGN3 <- factor(OCEC_DF$TYPE_BENIGN3)
OCEC_DF$TYPE_BENIGN3 <- droplevels(OCEC_DF$TYPE_BENIGN3) #drop unused
OCEC_DF$TYPE_BENIGN3  <- relevel(OCEC_DF$TYPE_BENIGN3 , ref = "ENDOMETRIAL CANCER") #set control
table(OCEC_DF$TYPE_BENIGN3, useNA = "a")
roc_results_np_OCED<- lapply(c("NOTCH2_NP","CTNNB1_NP","DLL1_NP","HES1_NP"), function(col) {
  roc(response = OCEC_DF$TYPE_BENIGN3, predictor = OCEC_DF[[col]])})
names(roc_results_np_OCED) <- c("NOTCH2_NP","CTNNB1_NP","DLL1_NP","HES1_NP")
roc_results_np_OCED
#extract the aucs
auc_values_np_ECOC <- sapply(roc_results_np_OCED, function(roc_obj) {auc(roc_obj)})
auc_values_np_ECOC #extracted aucs
#roc figure 
roc_plot_OCEC <- function() {
  par(pty = "s") #sets square
  plot.roc(roc_results_np_OCED[["NOTCH2_NP"]], print.auc = F, col = "#dcbeff",
           cex.main=0.8, 
           main ="Uterine lavage biomarkers EC vs OC",
           #xlab = "1 - Specifiškumas", 
           #ylab = "Jautrumas", 
           legacy.axes = T) #title
  lines(roc_results_np_OCED[["CTNNB1_NP"]], col = "#911eb4", lwd =2) 
  lines(roc_results_np_OCED[["DLL1_NP"]], col ="#ffd8b1", lwd =2) 
  lines(roc_results_np_OCED[["HES1_NP"]], col = "#42d4f4", lwd =2) 
  legend("bottomright", legend = c( expression(italic("NOTCH2")),
                                    expression(italic("CTNNB1")),
                                    expression(italic("DLL1")), 
                                    expression(italic("HES1"))
  ),
  
  col = c("#dcbeff", "#911eb4", "#ffd8b1", "#42d4f4"), lty = 1, 
  cex = 0.7, lwd =3)
}
#plot
roc_plot_OCEC()
#roc table
coords_results_np_OCEC<- lapply(roc_results_np_OCED, function(roc_obj) {
  coords(roc_obj, "best", ret = c("threshold", "accuracy", "sensitivity",
                                  "specificity"),
         transpose = FALSE)
})
coords_results_np_OCEC
#create df
results_np_OCEC<- data.frame(
  Predictor = c("NOTCH2_NP","CTNNB1_NP","DLL1_NP","HES1_NP"),
  AUC = auc_values_np_ECOC,
  do.call(rbind,coords_results_np_OCEC) 
)
results_np_OCEC
#change names
results_np_OCEC$Predictor <- c("NOTCH2",#fixed
                               "CTNNB1", "DLL1", "HES1")
#nice formating of the Table metrics for ROC OC
gt_table_np_OCEC <- results_np_OCEC %>%
  gt() %>%
  tab_header(
    title = "ROC measures for uterine lavage biomarkers", 
    subtitle = "Ovarian cancer vs endometrial cancer") %>%
  fmt_number(
    columns = everything(),
    decimals = 3
  ) %>%
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_body(columns = vars(Predictor))
  )
#show
gt_table_np_OCEC
# Save the plot as a PNG file
png("C:/Users/Ieva/rprojects/outputs_all/LIQUID/NP_OCEC_ROC_20260618.png",
    width = 11, height = 11, res = 300, units = "cm")
roc_plot_OCEC()
mtext(
  "A",
  side = 3,
  adj = -0.2,
  line = 1,
  cex = 1.2,
  font = 2
)
dev.off()
#there is no other convieneat way to save gt outputs
gtsave(gt_table_np_OCEC,vwidth = 10000,   
       filename = "C:/Users/Ieva/rprojects/outputs_all/LIQUID/NP_OCEC_ROCtable_20260618.png")
#import images
roc_imageOCEC  <- image_read("C:/Users/Ieva/rprojects/outputs_all/LIQUID/NP_OCEC_ROC_20260618.png")
table_imageOCEC <- image_read("C:/Users/Ieva/rprojects/outputs_all/LIQUID/NP_OCEC_ROCtable_20260618.png")

# resize table to match ROC image width
table_imageOCEC <- image_resize(table_imageOCEC,
                                paste0(image_info(roc_imageOCEC)$width, "x"))
# combine vertically
combinedOCEC <- image_append(c(roc_imageOCEC, table_imageOCEC), stack = TRUE)

# save
image_write(combinedOCEC,
            "C:/Users/Ieva/rprojects/outputs_all/LIQUID/NP_OCEC_ROCcombined_20260618.png")
##ROC OC (HGSOC+OTHER) VS BENIGN / RRS########################################
#make OC VS NON-cancer DF
OC_BEN_RRS_DF<- LAVAGE_df %>%
  filter(TYPE_BENIGN3%in% c("OC", "BENIGN"))%>% #filter for right samples
  dplyr::select("TYPE_BENIGN3", "NOTCH2_NP","CTNNB1_NP","DLL1_NP","HES1_NP"
  )
OC_BEN_RRS_DF$TYPE_BENIGN3 <- factor(OC_BEN_RRS_DF$TYPE_BENIGN3)
OC_BEN_RRS_DF$TYPE_BENIGN3 <- droplevels(OC_BEN_RRS_DF$TYPE_BENIGN3) #drop unused
OC_BEN_RRS_DF$TYPE_BENIGN3  <- relevel(OC_BEN_RRS_DF$TYPE_BENIGN3 , ref = "BENIGN") #set control
table(OC_BEN_RRS_DF$TYPE_BENIGN3, useNA = "a")
roc_results_np_OC_BEN_RRS<- lapply(c("NOTCH2_NP","CTNNB1_NP","DLL1_NP","HES1_NP"), function(col) {
  roc(response = OC_BEN_RRS_DF$TYPE_BENIGN3, predictor = OC_BEN_RRS_DF[[col]])})
names(roc_results_np_OC_BEN_RRS) <- c("NOTCH2_NP","CTNNB1_NP","DLL1_NP","HES1_NP")
roc_results_np_OC_BEN_RRS
#extract the aucs
auc_values_np_OC_BEN_RRS <- sapply(roc_results_np_OC_BEN_RRS, function(roc_obj) {auc(roc_obj)})
auc_values_np_OC_BEN_RRS #extracted aucs
#roc figure 
roc_plot_OC_BEN_RRS <- function() {
  par(pty = "s") #sets square
  plot.roc(roc_results_np_OC_BEN_RRS[["NOTCH2_NP"]], print.auc = F, col = "#dcbeff",
           cex.main=0.8, 
           main ="Uterine lavage biomarkers OC vs Non-cancer samples",
           #xlab = "1 - Specifiškumas", 
           #ylab = "Jautrumas", 
           legacy.axes = T) #title
  lines(roc_results_np_OC_BEN_RRS[["CTNNB1_NP"]], col = "#911eb4", lwd =2) 
  lines(roc_results_np_OC_BEN_RRS[["DLL1_NP"]], col ="#ffd8b1", lwd =2) 
  lines(roc_results_np_OC_BEN_RRS[["HES1_NP"]], col = "#42d4f4", lwd =2) 
  legend("bottomright", legend = c( expression(italic("NOTCH2")),
                                    expression(italic("CTNNB1")),
                                    expression(italic("DLL1")), 
                                    expression(italic("HES1"))
  ),
  
  col = c("#dcbeff", "#911eb4", "#ffd8b1", "#42d4f4"), lty = 1, 
  cex = 0.7, lwd =3)
}
#plot
roc_plot_OC_BEN_RRS()
#roc table
coords_results_np_OC_BEN_RRS<- lapply(roc_results_np_OC_BEN_RRS, function(roc_obj) {
  coords(roc_obj, "best", ret = c("threshold", "accuracy", "sensitivity",
                                  "specificity"),
         transpose = FALSE)
})
coords_results_np_OC_BEN_RRS
#create df
results_np_OC_BEN_RRS<- data.frame(
  Predictor = c("NOTCH2_NP","CTNNB1_NP","DLL1_NP","HES1_NP"),
  AUC = auc_values_np_OC_BEN_RRS,
  do.call(rbind,coords_results_np_OC_BEN_RRS) 
)
results_np_OC_BEN_RRS
#change names
results_np_OC_BEN_RRS$Predictor <- c("NOTCH2",#fixed
                                     "CTNNB1", "DLL1", "HES1")
#nice formating of the Table metrics for ROC OC
gt_table_np_OC_BEN_RRS <- results_np_OC_BEN_RRS %>%
  gt() %>%
  tab_header(
    title = "ROC measures for uterine lavage biomarkers", 
    subtitle = "Non-cancer vs ovarian cancer") %>%
  fmt_number(
    columns = everything(),
    decimals = 3
  ) %>%
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_body(columns = vars(Predictor))
  )
#show
gt_table_np_OC_BEN_RRS
# Save the plot as a PNG file
png("C:/Users/Ieva/rprojects/outputs_all/LIQUID/NP_OC_BEN_RRS_ROC_20260618.png",
    width = 11, height = 11, res = 300, units = "cm")
roc_plot_OC_BEN_RRS()
mtext(
  "B",
  side = 3,
  adj = -0.2,
  line = 1,
  cex = 1.2,
  font = 2
)
dev.off()
#there is no other convieneat way to save gt outputs
gtsave(gt_table_np_OC_BEN_RRS,vwidth = 10000,   
       filename = "C:/Users/Ieva/rprojects/outputs_all/LIQUID/NP_OC_BEN_RRS_ROCtable_20260618.png")
#import images
roc_imageOC_BEN_RRS  <- image_read("C:/Users/Ieva/rprojects/outputs_all/LIQUID/NP_OC_BEN_RRS_ROC_20260618.png")
table_imageOC_BEN_RRS <- image_read("C:/Users/Ieva/rprojects/outputs_all/LIQUID/NP_OC_BEN_RRS_ROCtable_20260618.png")

# resize table to match ROC image width
table_imageOC_BEN_RRS <- image_resize(table_imageOC_BEN_RRS,
                                      paste0(image_info(roc_imageOC_BEN_RRS)$width, "x"))
# combine vertically
combinedOC_BEN_RRS <- image_append(c(roc_imageOC_BEN_RRS, table_imageOC_BEN_RRS), stack = TRUE)

# save
image_write(combinedOC_BEN_RRS,
            "C:/Users/Ieva/rprojects/outputs_all/LIQUID/NP_OC_BEN_RRS_ROCcombined_20260618.png")
