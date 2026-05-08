#KN-  liquid 2026 05 08
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
  "I&II",   "III&IV", 0.7019, -1, "NOTCH2_NP",
  "I&II",   "III&IV", 0.07831, -1.5, "DLL1_NP",
  "I&II",   "III&IV", 0.03262, -2, "HES1_NP",

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
  facet_wrap(.~ variable, nrow = 2, scales = "free") +
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
ggsave("C:/Users/Ieva/rprojects/outputs_all/LIQUID/stage_oc_20260508.png",
       plot = STAGE_OC,
       width = 15,
       height = 16,
       units = "cm",
       dpi = 400)
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

#age correlation EC ###############################
cor.test(KN_OC$NOTCH2_NP,
         KN_OC$CA125_num, method = "spearman")
cor.test(KN_OC$DLL1_NP,
         KN_OC$CA125_num, method = "spearman")
cor.test(KN_OC$HES1_NP,
         KN_OC$CA125_num, method = "spearman") 
cor.test(KN_OC$CTNNB1_NP,
         KN_OC$CA125_num, method = "spearman")
