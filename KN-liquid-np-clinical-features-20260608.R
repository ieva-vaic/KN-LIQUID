#KN-  liquid 2026 05 08, 2026 09 04
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

#fix CA125
LAVAGE_df$CA125
LAVAGE_df <- LAVAGE_df %>%
  mutate(
    CA125 = na_if(CA125, "NA"),
    CA125 = na_if(CA125, "Neatlikta"),
    CA125 = as.numeric(CA125),
    CA125_group = case_when(
      is.na(CA125) ~ NA_character_,
      CA125 >= 35  ~ "CA125 increase",
      CA125 < 35   ~ "No CA125 increase"
    )
  )
table(LAVAGE_df$CA125_group, useNA = "ifany")#33 na
table(LAVAGE_df$CA125, useNA = "a") #33 na

#CHEK STAGE / GRADE ###########################
#make an endometrial cancer df
KN_ENDOMETRIAL <- LAVAGE_df %>%
  filter(TYPE == "ENDOMETRIAL CANCER")
#chek stage
table(KN_ENDOMETRIAL$Stage, KN_ENDOMETRIAL$TYPE) #only 2 stage II
#chek grade
table(KN_ENDOMETRIAL$Grade_simple, KN_ENDOMETRIAL$TYPE) #6 vs 5 grade g1 vs G2
#make OC df
KN_OC <- LAVAGE_df %>%
  filter(TYPE %in% c("HGSOC", "OTHER") )
#chek stage OC
table(KN_OC$Stage_simple, KN_OC$TYPE) #only 3 stage II
#chek grade OC
table(KN_OC$Grade_simple, KN_OC$TYPE) #7 vs 47

#NORMALCY stage OC##########################################
normality_resultsOC <- KN_OC %>%
  pivot_longer(
    cols = all_of(DATA),
    names_to = "variable",
    values_to = "value"
  ) %>%
  group_by(Stage_simple, variable) %>%
  summarise(
    n = sum(!is.na(value)),
    shapiro_p = if (n >= 3 & n <= 5000) shapiro.test(value)$p.value else NA_real_,
    .groups = "drop"
  ) %>%
  mutate(
    normal = ifelse(shapiro_p > 0.05, TRUE, FALSE)
  )
normality_resultsOC #all normal

#NORMALCY grade OC##########################################
normality_resultsOCgrade <- KN_OC %>%
  pivot_longer(
    cols = all_of(DATA),
    names_to = "variable",
    values_to = "value"
  ) %>%
  group_by(Grade_simple, variable) %>%
  summarise(
    n = sum(!is.na(value)),
    shapiro_p = if (n >= 3 & n <= 5000) shapiro.test(value)$p.value else NA_real_,
    .groups = "drop"
  ) %>%
  mutate(
    normal = ifelse(shapiro_p > 0.05, TRUE, FALSE)
  )
normality_resultsOCgrade #all normal

#NORMALCY grade EC##########################################
normality_resultsECgrade <- KN_ENDOMETRIAL %>%
  pivot_longer(
    cols = all_of(DATA),
    names_to = "variable",
    values_to = "value"
  ) %>%
  group_by(Grade_simple, variable) %>%
  summarise(
    n = sum(!is.na(value)),
    shapiro_p = if (n >= 3 & n <= 5000) shapiro.test(value)$p.value else NA_real_,
    .groups = "drop"
  ) %>%
  mutate(
    normal = ifelse(shapiro_p > 0.05, TRUE, FALSE)
  )
normality_resultsECgrade #all normal

#Variance stage OC#################################
#make sure stage is a factor
KN_OC$Stage_simple <- factor(KN_OC$Stage_simple)
#variance KN
variance_results_4 <- lapply(DATA, function(v) {
  # Select variable and STAGE
  df <- KN_OC[, c(v, "Stage_simple")]
  
  # Remove NA
  df <- df[!is.na(df[[v]]), ]
  
  # Only test if at least 2 groups have data
  if(length(unique(df$Stage_simple)) > 1) {
    test <- leveneTest(df[[v]] ~ df$Stage_simple, center = median)
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
variance_results_4 <- bind_rows(variance_results_4)

# View the results
variance_results_4 #all equal
#Variance grade OC#################################
#make sure stage is a factor
KN_OC$Grade_simple <- factor(KN_OC$Grade_simple)
#variance KN
variance_results_3 <- lapply(DATA, function(v) {
  # Select variable and STAGE
  df <- KN_OC[, c(v, "Grade_simple")]
  
  # Remove NA
  df <- df[!is.na(df[[v]]), ]
  
  # Only test if at least 2 groups have data
  if(length(unique(df$Grade_simple)) > 1) {
    test <- leveneTest(df[[v]] ~ df$Grade_simple, center = median)
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
variance_results_3 <- bind_rows(variance_results_3)

# View the results
variance_results_3 #all equal

#Variance grade EC#################################
#make sure stage is a factor
KN_ENDOMETRIAL$Grade_simple <- factor(KN_ENDOMETRIAL$Grade_simple)
#variance KN
variance_results_2 <- lapply(DATA, function(v) {
  # Select variable and STAGE
  df <- KN_ENDOMETRIAL[, c(v, "Grade_simple")]
  
  # Remove NA
  df <- df[!is.na(df[[v]]), ]
  
  # Only test if at least 2 groups have data
  if(length(unique(df$Grade_simple)) > 1) {
    test <- leveneTest(df[[v]] ~ df$Grade_simple, center = median)
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
variance_results_2 <- bind_rows(variance_results_2)

# View the results
variance_results_2 #all equal

#ANOVA KN STAGE##############################
#KN ANOVA:
# CTNNB1
anova_ctnnb1stage <- aov(CTNNB1_NP ~ Stage_simple, data = KN_OC)
summary(anova_ctnnb1stage) #not significant 0.818
#DLL1
anova_dll1stage <- aov(DLL1_NP ~ Stage_simple, data = KN_OC)
summary(anova_dll1stage) #not significant (but almost significant) 0.0479 
TukeyHSD(anova_dll1stage) #significant: 3 vs 2,p = 0.0459281, bet tik 2 vnt stage2
#NOTCH2
anova_notch2stage <- aov(NOTCH2_NP ~ Stage_simple, data = KN_OC) 
summary(anova_notch2stage) #not significant (but almost significant) 0.0973
TukeyHSD(anova_notch2stage)#not significant 
#HES1
anova_hes1stage <- aov(HES1_NP ~ Stage_simple, data = KN_OC) 
summary(anova_hes1stage) #not significant

#GROUPED OC STAGE###########################
#group for convenience
KN_OC$Stage_grouped <- ifelse(
  KN_OC$Stage_simple %in% c(1, 2),
  "I&II",
  "III&IV"
)

KN_OC$Stage_grouped <- factor(
  KN_OC$Stage_grouped,
  levels = c("I&II", "III&IV")
)
table(KN_OC$Stage_grouped , useNA = "a")
#NORMALCY stage OC 2 groups##########################################
normality_resultsOC2 <- KN_OC %>%
  pivot_longer(
    cols = all_of(DATA),
    names_to = "variable",
    values_to = "value"
  ) %>%
  group_by(Stage_grouped, variable) %>%
  summarise(
    n = sum(!is.na(value)),
    shapiro_p = if (n >= 3 & n <= 5000) shapiro.test(value)$p.value else NA_real_,
    .groups = "drop"
  ) %>%
  mutate(
    normal = ifelse(shapiro_p > 0.05, TRUE, FALSE)
  )
normality_resultsOC2 #all normal

#Variance stage OC grouped####################
var.test(CTNNB1_NP ~ Stage_grouped, data = KN_OC)
var.test(NOTCH2_NP ~ Stage_grouped, data = KN_OC)
var.test(HES1_NP ~ Stage_grouped, data = KN_OC)
var.test(DLL1_NP ~ Stage_grouped, data = KN_OC)
#all normal variance
#t.tests OC grouped stage###############
t.test(NOTCH2_NP ~ Stage_grouped,
       data = KN_OC,
       var.equal = TRUE)
t.test(DLL1_NP ~ Stage_grouped, #0.07831
       data = KN_OC,
       var.equal = TRUE)
t.test(CTNNB1_NP ~ Stage_grouped,
       data = KN_OC,
       var.equal = TRUE)
t.test(HES1_NP~ Stage_grouped,
       data = KN_OC,
       var.equal = TRUE) #0.03262
#plot OC grouped stage####################
#make p values
each.vs.ref_sig <- tibble::tribble(
  ~group1, ~group2, ~p.adj,   ~y.position, ~variable,
  "I&II",   "III&IV", 0.149, -2, "CTNNB1_NP",
  "I&II",   "III&IV", 0.702, -1, "NOTCH2_NP",
  "I&II",   "III&IV", 0.078, -1.5, "DLL1_NP",
  "I&II",   "III&IV", 0.033, -2, "HES1_NP",

)
#melt table for expression
GroupNP_table <- melt(KN_OC[, c(40,15:18)],
                      id.vars="Stage_grouped",
                      measure.vars=c("NOTCH2_NP",
                                     "CTNNB1_NP",
                                     "DLL1_NP",
                                     "HES1_NP"))

STAGE_OC <- ggplot(GroupNP_table, aes(x=Stage_grouped , y=value, fill = variable)) +
  geom_boxplot( outlier.shape = NA , alpha=0.3, aes(fill = Stage_grouped )) +
  geom_jitter(aes(color = Stage_grouped ), size=1, alpha=0.5) +
  ylab(label = expression("Gene expression, normalized to  " * italic("GAPDH"))) + 
  facet_wrap(.~ variable, nrow = 2, scales = "free",
             labeller = labeller(
               variable = c(
                 "CTNNB1_NP" = "CTNNB1",
                 "DLL1_NP" = "DLL1",
                 "HES1_NP" = "HES1",
                 "NOTCH2_NP" = "NOTCH2"
               ))
             ) +
  add_pvalue(each.vs.ref_sig, label = "p.adj") + #pvalue
  theme_minimal()+
  theme(
    strip.text.x = element_text(
      size = 12, face = "bold.italic"
    ),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5))+
  labs(x=NULL,
       title = "Gene expression in uterine lavage by OC stage")+
  stat_boxplot(geom ='errorbar')+
  #scale_fill_manual(values = custom_colors) +
  #scale_color_manual(values = custom_colors) +
  scale_y_continuous(labels = function(x) 
    gsub("-", "\u2212", as.character(x)))+ #add long "-" signs
  scale_fill_manual(values = c(
    "III&IV" = "#3C5488",
    "I&II" = "#4DBBD5"
  )) +
  scale_color_manual(values = c(
    "III&IV" = "#3C5488",
    "I&II" = "#4DBBD5"
  ))

STAGE_OC
#save
ggsave("C:/Users/Ieva/rprojects/outputs_all/LIQUID/stage_oc_20260709.png",
       plot = STAGE_OC,
       width = 12,
       height = 16,
       units = "cm",
       dpi = 200)
#t.tests OC grade###########################
t.test(NOTCH2_NP ~ Grade_simple,
       data = KN_OC,
       var.equal = TRUE)#0.123
t.test(DLL1_NP ~ Grade_simple, 
       data = KN_OC,
       var.equal = TRUE)# 0.8142
t.test(CTNNB1_NP ~ Grade_simple,
       data = KN_OC,
       var.equal = TRUE)#0.2212
t.test(HES1_NP~ Grade_simple,
       data = KN_OC,
       var.equal = TRUE) #0.1823
#t.tests EC grouped stage###############
t.test(NOTCH2_NP ~ Grade_simple,
       data = KN_ENDOMETRIAL,
       var.equal = TRUE)#0.4453
t.test(DLL1_NP ~ Grade_simple, 
       data = KN_ENDOMETRIAL,
       var.equal = TRUE)# 0.5434
t.test(CTNNB1_NP ~ Grade_simple,
       data = KN_ENDOMETRIAL,
       var.equal = TRUE)#0.6287
t.test(HES1_NP~ Grade_simple,
       data = KN_ENDOMETRIAL,
       var.equal = TRUE) #0.6062

#Overall normalcy EC#####################
shapiro.test(KN_ENDOMETRIAL$Age)
KN_ENDOMETRIAL$CA125_num <- as.numeric(KN_ENDOMETRIAL$CA125_num)
shapiro.test(na.omit(KN_ENDOMETRIAL$CA125_num))#too little values

shapiro.test(KN_ENDOMETRIAL$NOTCH2_NP) #not normal0.01915
shapiro.test(KN_ENDOMETRIAL$DLL1_NP)
shapiro.test(KN_ENDOMETRIAL$CTNNB1_NP)
shapiro.test(KN_ENDOMETRIAL$HES1_NP)

#Overall normalcy OC#####################
shapiro.test(KN_OC$Age)
KN_OC$CA125_num <- as.numeric(KN_OC$CA125_num)
shapiro.test(na.omit(KN_OC$CA125_num))#not normal p-value = 2.647e-09

shapiro.test(KN_OC$NOTCH2_NP) 
shapiro.test(KN_OC$DLL1_NP)
shapiro.test(KN_OC$CTNNB1_NP)# not normal p-value = 0.02548
shapiro.test(KN_OC$HES1_NP)

#age correlation EC ###############################
cor.test(KN_ENDOMETRIAL$NOTCH2_NP,
         KN_ENDOMETRIAL$Age, method = "spearman")
cor.test(KN_ENDOMETRIAL$DLL1_NP,
         KN_ENDOMETRIAL$Age, method = "pearson")
cor.test(KN_ENDOMETRIAL$HES1_NP,
         KN_ENDOMETRIAL$Age, method = "pearson") 
cor.test(KN_ENDOMETRIAL$CTNNB1_NP,
         KN_ENDOMETRIAL$Age, method = "pearson")
#age correlation EC ###############################
cor.test(KN_OC$NOTCH2_NP,
         KN_OC$Age, method = "pearson")
cor.test(KN_OC$DLL1_NP,
         KN_OC$Age, method = "pearson")
cor.test(KN_OC$HES1_NP,
         KN_OC$Age, method = "pearson") 
cor.test(KN_OC$CTNNB1_NP,
         KN_OC$Age, method = "spearman")

#ca125 correlation OC ###############################
cor.test(KN_OC$NOTCH2_NP,
         KN_OC$CA125_num, method = "spearman")
cor.test(KN_OC$DLL1_NP,
         KN_OC$CA125_num, method = "spearman")
cor.test(KN_OC$HES1_NP,
         KN_OC$CA125_num, method = "spearman") 
cor.test(KN_OC$CTNNB1_NP,
         KN_OC$CA125_num, method = "spearman")

CA_NOTCH2 <- ggplot(KN_OC, aes(x = NOTCH2_NP, y = CA125_num)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  stat_cor(method = "spearman") +
  theme_classic() +
  labs(
    x = "NOTCH2 expression",
    y = "CA125"
  )

CA_DLL1 <- ggplot(KN_OC, aes(x = DLL1_NP, y = CA125_num)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  stat_cor(method = "spearman") +
  theme_classic() +
  labs(
    x = "DLL1 expression",
    y = "CA125"
  )


CA_HES1 <- ggplot(KN_OC, aes(x = HES1_NP, y = CA125_num)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  stat_cor(method = "spearman") +
  theme_classic() +
  labs(
    x = "HES1 expression",
    y = "CA125"
  )

CA_CTNNB1 <-ggplot(KN_OC, aes(x = CTNNB1_NP, y = CA125_num)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  stat_cor(method = "spearman") +
  theme_classic() +
  labs(
    x = "CTNNB1 expression",
    y = "CA125"
  )


combined_plot <- wrap_plots(
  CA_NOTCH2,
  CA_DLL1,
  CA_HES1,
  CA_CTNNB1,
  ncol = 2
)

combined_plot

#Stage: stage 1 vs other #################################
#group for convenience
KN_OC$Stage_grouped2 <- ifelse(
  KN_OC$Stage_simple == 1,
  "I",
  "II&III&IV"
)
table(KN_OC$Stage_grouped2, useNA = "a") #10 vs 51
##NORMALCY stage OC 1 vs other groups##########################################
normality_results_stage2 <- KN_OC %>%
  pivot_longer(
    cols = all_of(DATA),
    names_to = "variable",
    values_to = "value"
  ) %>%
  group_by(Stage_grouped2, variable) %>%
  summarise(
    n = sum(!is.na(value)),
    shapiro_p = if (n >= 3 & n <= 5000) shapiro.test(value)$p.value else NA_real_,
    .groups = "drop"
  ) %>%
  mutate(
    normal = ifelse(shapiro_p > 0.05, TRUE, FALSE)
  )
normality_results_stage2 #all normal

##Variance stage OC grouped####################
var.test(CTNNB1_NP ~ Stage_grouped2, data = KN_OC)
var.test(NOTCH2_NP ~ Stage_grouped2, data = KN_OC)
var.test(HES1_NP ~ Stage_grouped2, data = KN_OC)
var.test(DLL1_NP ~ Stage_grouped2, data = KN_OC)
#all normal variance

##t.tests OC grouped stage 1 vs other ###############
t.test(NOTCH2_NP ~ Stage_grouped2,
       data = KN_OC,
       var.equal = TRUE)
t.test(DLL1_NP ~ Stage_grouped2, #0.02384
       data = KN_OC,
       var.equal = TRUE)#0.5937
t.test(CTNNB1_NP ~ Stage_grouped2,
       data = KN_OC,
       var.equal = TRUE)#0.4456
t.test(HES1_NP~ Stage_grouped2,
       data = KN_OC,
       var.equal = TRUE) #0.07291

##plot OC grouped stage####################
#make p values
each.vs.ref_sig_stage2 <- tibble::tribble(
  ~group1, ~group2, ~p.adj,   ~y.position, ~variable,
  "I",   "II&III&IV", 0.02384, -2, "CTNNB1_NP",
  "I",   "II&III&IV", 0.5937, -1, "NOTCH2_NP",
  "I",   "II&III&IV", 0.4456, -1.5, "DLL1_NP",
  "I",   "II&III&IV", 0.07291, -2, "HES1_NP",
  
)
#melt table for expression
GroupNP_table2 <- melt(KN_OC[, c(41,15:18)],
                       id.vars="Stage_grouped2",
                       measure.vars=c("NOTCH2_NP",
                                      "CTNNB1_NP",
                                      "DLL1_NP",
                                      "HES1_NP"))

STAGE_OC2 <- ggplot(GroupNP_table2, aes(x=Stage_grouped2 , y=value, fill = variable)) +
  geom_boxplot( outlier.shape = NA , alpha=0.3, aes(fill = Stage_grouped2 )) +
  geom_jitter(aes(color = Stage_grouped2 ), size=1, alpha=0.5) +
  ylab(label = expression("Gene expression, normalized to  " * italic("GAPDH"))) + 
  facet_wrap(.~ variable, nrow = 2, scales = "free",
             labeller = labeller(
               variable = c(
                 "CTNNB1_NP" = "CTNNB1",
                 "DLL1_NP" = "DLL1",
                 "HES1_NP" = "HES1",
                 "NOTCH2_NP" = "NOTCH2"
               ))
  ) +
  add_pvalue(each.vs.ref_sig_stage2, label = "p.adj") + #pvalue
  theme_minimal()+
  theme(
    strip.text.x = element_text(
      size = 12, face = "bold.italic"
    ),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5))+
  labs(x=NULL,
       title = "Gene expression in uterine lavage by OC stage")+
  stat_boxplot(geom ='errorbar')+
  #scale_fill_manual(values = custom_colors) +
  #scale_color_manual(values = custom_colors) +
  scale_y_continuous(labels = function(x) 
    gsub("-", "\u2212", as.character(x)))+ #add long "-" signs
  scale_fill_manual(values = c(
    "II&III&IV" = "#3C5488",
    "I"   = "#4DBBD5"
  )) +
  scale_color_manual(values = c(
    "II&III&IV" = "#3C5488",
    "I"   = "#4DBBD5"
  ))

STAGE_OC2
#CA125 in OC and benign only#######################
KN_OC_BEN <- LAVAGE_df %>%
  filter(TYPE != "ENDOMETRIAL CANCER")
##NORMALCY stage OC 1 vs other groups##########################################
normality_results_stageca <- KN_OC_BEN %>%
  pivot_longer(
    cols = all_of(DATA),
    names_to = "variable",
    values_to = "value"
  ) %>%
  group_by(CA125_group, variable) %>%
  summarise(
    n = sum(!is.na(value)),
    shapiro_p = if (n >= 3 & n <= 5000) shapiro.test(value)$p.value else NA_real_,
    .groups = "drop"
  ) %>%
  mutate(
    normal = ifelse(shapiro_p > 0.05, TRUE, FALSE)
  )
normality_results_stageca #all normal except CTNNB1

##Variance stage OC grouped####################
var.test(CTNNB1_NP ~ CA125_group, data = KN_OC_BEN)
var.test(NOTCH2_NP ~ CA125_group, data = KN_OC_BEN)
var.test(HES1_NP ~ CA125_group, data = KN_OC_BEN)
var.test(DLL1_NP ~ CA125_group, data = KN_OC_BEN)
#all normal variance

##t.tests OC grouped CA125 1 vs other ###############
t.test(NOTCH2_NP ~ CA125_group,
       data = KN_OC_BEN,
       var.equal = TRUE)#0.08071
t.test(DLL1_NP ~ CA125_group, 
       data = KN_OC_BEN,
       var.equal = TRUE)#0.8902
wilcox.test(
  CTNNB1_NP ~ CA125_group,
  data = KN_OC_BEN,
  exact = FALSE
)#0.9306
t.test(HES1_NP~ CA125_group,
       data = KN_OC_BEN,
       var.equal = TRUE) #0.01266


##plot OC grouped stage####################
#make p values
each.vs.ref_sig_ca <- tibble::tribble(
  ~group1, ~group2, ~p.adj,   ~y.position, ~variable,
  "CA125 increase",   "No CA125 increase", 0.931, -2, "CTNNB1_NP",
  "CA125 increase",   "No CA125 increase",0.0807, -1, "NOTCH2_NP",
  "CA125 increase",   "No CA125 increase", 0.890, -1.5, "DLL1_NP",
  "CA125 increase",   "No CA125 increase", 0.0127, -2, "HES1_NP",
  
)
#melt table for expression
GroupNP_ca <- melt(KN_OC_BEN[, c(40,15:18)],
                   id.vars="CA125_group",
                   measure.vars=c("NOTCH2_NP",
                                  "CTNNB1_NP",
                                  "DLL1_NP",
                                  "HES1_NP"))
GroupNP_ca <- GroupNP_ca[!is.na(GroupNP_ca$CA125_group), ]
ca_OCplot <- ggplot(GroupNP_ca, aes(x=CA125_group , y=value, fill = variable)) +
  geom_boxplot( outlier.shape = NA , alpha=0.3, aes(fill = CA125_group )) +
  geom_jitter(aes(color = CA125_group ), size=1, alpha=0.5) +
  ylab(label = expression("Gene expression, normalized to  " * italic("GAPDH"))) + 
  facet_wrap(.~ variable, nrow = 2, scales = "free",
             labeller = labeller(
               variable = c(
                 "CTNNB1_NP" = "CTNNB1",
                 "DLL1_NP" = "DLL1",
                 "HES1_NP" = "HES1",
                 "NOTCH2_NP" = "NOTCH2"
               ))
  ) +
  add_pvalue(each.vs.ref_sig_ca, label = "p.adj") + #pvalue
  theme_minimal()+
  theme(
    strip.text.x = element_text(
      size = 12, face = "bold.italic"
    ),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5))+
  labs(x=NULL,
       title = "Gene expression in uterine lavage by CA125 status in OC ")+
  stat_boxplot(geom ='errorbar')+
  #scale_fill_manual(values = custom_colors) +
  #scale_color_manual(values = custom_colors) +
  scale_y_continuous(labels = function(x) 
    gsub("-", "\u2212", as.character(x)))+ #add long "-" signs
  scale_fill_manual(values = c(
    "CA125 increase" = "#3C5488",
    "No CA125 increase"   = "#4DBBD5"
  )) +
  scale_color_manual(values = c(
    "CA125 increase" = "#3C5488",
    "No CA125 increase"   = "#4DBBD5"
  ))

ca_OCplot
#CA125 in all cases ##############################
LAVAGE_df_ca125_ca125 <-  LAVAGE %>%
  filter(!is.na(CA125_group))
##NORMALCY CA125##########################################
normality_results_stageca <- LAVAGE_df_ca125 %>%
  pivot_longer(
    cols = all_of(DATA),
    names_to = "variable",
    values_to = "value"
  ) %>%
  group_by(CA125_group, variable) %>%
  summarise(
    n = sum(!is.na(value)),
    shapiro_p = if (n >= 3 & n <= 5000) shapiro.test(value)$p.value else NA_real_,
    .groups = "drop"
  ) %>%
  mutate(
    normal = ifelse(shapiro_p > 0.05, TRUE, FALSE)
  )
normality_results_stageca #all normal except CTNNB1

##Variance stage OC grouped####################
var.test(CTNNB1_NP ~ CA125_group, data = LAVAGE_df_ca125)
var.test(NOTCH2_NP ~ CA125_group, data = LAVAGE_df_ca125)
var.test(HES1_NP ~ CA125_group, data = LAVAGE_df_ca125)
var.test(DLL1_NP ~ CA125_group, data = LAVAGE_df_ca125)
#all normal variance

##t.tests OC grouped CA125 1 vs other ###############
t.test(NOTCH2_NP ~ CA125_group,
       data = LAVAGE_df_ca125,
       var.equal = TRUE)#0.02086
t.test(DLL1_NP ~ CA125_group, 
       data = LAVAGE_df_ca125,
       var.equal = TRUE)#0.5112
wilcox.test(
  CTNNB1_NP ~ CA125_group,
  data = LAVAGE_df_ca125,
  exact = FALSE
)#0.7397
t.test(HES1_NP~ CA125_group,
       data = LAVAGE_df_ca125,
       var.equal = TRUE) #0.005559


##plot OC grouped stage####################
#make p values
each.vs.ref_sig_ca <- tibble::tribble(
  ~group1, ~group2, ~p.adj,   ~y.position, ~variable,
  "CA125 increase",   "No CA125 increase", 0.740, -2, "CTNNB1_NP",
  "CA125 increase",   "No CA125 increase",0.021, -1, "NOTCH2_NP",
  "CA125 increase",   "No CA125 increase", 0.511, -1.5, "DLL1_NP",
  "CA125 increase",   "No CA125 increase", 0.006, -2, "HES1_NP",
  
)
#melt table for expression
GroupNP_ca <- melt(LAVAGE_df_ca125[, c(40,15:18)],
                   id.vars="CA125_group",
                   measure.vars=c("NOTCH2_NP",
                                  "CTNNB1_NP",
                                  "DLL1_NP",
                                  "HES1_NP"))
GroupNP_ca <- GroupNP_ca[!is.na(GroupNP_ca$CA125_group), ]
ca_OCplot <- ggplot(GroupNP_ca, aes(x=CA125_group , y=value, fill = variable)) +
  geom_boxplot( outlier.shape = NA , alpha=0.3, aes(fill = CA125_group )) +
  geom_jitter(aes(color = CA125_group ), size=1, alpha=0.5) +
  ylab(label = expression("Gene expression, normalized to  " * italic("GAPDH"))) + 
  facet_wrap(.~ variable, nrow = 2, scales = "free",
             labeller = labeller(
               variable = c(
                 "CTNNB1_NP" = "CTNNB1",
                 "DLL1_NP" = "DLL1",
                 "HES1_NP" = "HES1",
                 "NOTCH2_NP" = "NOTCH2"
               ))
  ) +
  add_pvalue(each.vs.ref_sig_ca, label = "p.adj") + #pvalue
  theme_minimal()+
  theme(
    strip.text.x = element_text(
      size = 12, face = "bold.italic"
    ),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5))+
  labs(x=NULL,
       title = "Gene expression in uterine lavage by CA125 status")+
  stat_boxplot(geom ='errorbar')+
  #scale_fill_manual(values = custom_colors) +
  #scale_color_manual(values = custom_colors) +
  scale_y_continuous(labels = function(x) 
    gsub("-", "\u2212", as.character(x)))+ #add long "-" signs
  scale_fill_manual(values = c(
    "CA125 increase" = "#3C5488",
    "No CA125 increase"   = "#4DBBD5"
  )) +
  scale_color_manual(values = c(
    "CA125 increase" = "#3C5488",
    "No CA125 increase"   = "#4DBBD5"
  ))

ca_OCplot
#save
ggsave(
  filename = "C:/Users/Ieva/rprojects/outputs_all/LIQUID/ca125_20260908.png",
  plot = ca_OCplot,
  width = 6,
  height = 7,
  dpi = 300,
  bg = "white"
)