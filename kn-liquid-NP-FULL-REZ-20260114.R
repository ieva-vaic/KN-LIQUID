#KN-  liquid 2026 01 14, 16
#FULL NP data, after additional np cases
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
#data
#read RDS
LIQUID_DF <- readRDS("C:/Users/Ieva/rprojects/OTHER DATA/KN_LIQUID/liquid_20260114.RDS")
#test normalcy of my variables#####################################
DATA_names <- c("NOTCH2_NP","CTNNB1_NP","DLL1_NP","HES1_NP",
                "NOTCH2_URINE", "CTNNB1_URINE","DLL1_URINE","HES1_URINE",
                #"NOCTH2_P","CTNNB1_P", "DLL1_P","HES1_P"
                "NOTCH2_TUMOR","CTNNB1_TUMOR","DLL1_TUMOR","HES1_TUMOR" )
DATA <- c("NOTCH2_NP","CTNNB1_NP","DLL1_NP","HES1_NP" )
results_normalcy <- lapply(DATA_names, function(v) {
  by(LIQUID_DF[[v]], LIQUID_DF$TYPE, function(x) {
    if (length(na.omit(x)) >= 3)
      shapiro.test(x)$p.value
    else NA
  })
})
names(results_normalcy) <- DATA_names
results_normalcy #Not normal NP NOTCH2 endometrial cancer and np HES1 benign
#test variance
# Loop through all variables in DATA_names
variance_results <- lapply(DATA_names, function(v) {
  # Select variable and TYPE
  df <- LIQUID_DF[, c(v, "TYPE")]
  
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
variance_results #not equal in DLL1 urine
#ANOVA###########################
#ANOVA (for normal data) - CTNNB1 and DLL1
# CTNNB1
anova_ctnnb1 <- aov(CTNNB1_NP ~ TYPE, data = LIQUID_DF)
summary(anova_ctnnb1) #significant
TukeyHSD(anova_ctnnb1) #significant: ENDOMETRIAL CANCER-BENIGN, HGSOC-ENDOMETRIAL CANCER, RSS-ENDOMETRIAL CANCER 
#DLL1
anova_dll1 <- aov(DLL1_NP ~ TYPE, data = LIQUID_DF)
summary(anova_dll1) #not significant
TukeyHSD(anova_dll1) #significant: RSS-ENDOMETRIAL CANCER 
#ANOVA (for not normal data) - NOTCH2 and HES1
#NOTCH2
kruskal.test(NOTCH2_NP ~ TYPE, data = LIQUID_DF) #significant
pairwise.wilcox.test(LIQUID_DF$NOTCH2_NP, LIQUID_DF$TYPE, p.adjust.method = "bonferroni")
#ENDOMETRIAL CANCER-BENIGN, HGSOC-ENDOMETRIAL CANCER, RSS-ENDOMETRIAL CANCER 
#HES1
kruskal.test(HES1_NP ~ TYPE, data = LIQUID_DF) #significant
pairwise.wilcox.test(LIQUID_DF$HES1_NP, LIQUID_DF$TYPE, p.adjust.method = "bonferroni") #HGSOC-ENDOMETRIAL CANCER
#boxplot full np ###########################
#melt table for expression
GroupNP_table <- melt(LIQUID_DF[, c(12,15:18)], id.vars="TYPE",  measure.vars=c("NOTCH2_NP","CTNNB1_NP","DLL1_NP","HES1_NP"))
#make p values
each.vs.ref_sig <- tibble::tribble(
    ~group1, ~group2, ~p.adj,   ~y.position, ~variable,
    "ENDOMETRIAL CANCER",   "BENIGN", 0.0016350, -2, "CTNNB1_NP",
    "ENDOMETRIAL CANCER",   "HGSOC", 0.0029426, -1, "CTNNB1_NP",
    "ENDOMETRIAL CANCER",   "RSS", 0.0004705, -1.5, "CTNNB1_NP",
    "ENDOMETRIAL CANCER",   "BENIGN", 0.0104, -2, "NOTCH2_NP",
    "ENDOMETRIAL CANCER",   "HGSOC", 0.0014, -1, "NOTCH2_NP",
    "ENDOMETRIAL CANCER",   "RSS", 0.0092, -1.5, "NOTCH2_NP",
    "ENDOMETRIAL CANCER",   "HGSOC", 0.013, -2, "HES1_NP",
    "ENDOMETRIAL CANCER",   "RSS", 0.0015864, -6, "DLL1_NP",
  )

TYPE_FULL_plot <- ggplot(GroupNP_table, aes(x=TYPE , y=value, fill = variable)) +
  geom_boxplot( outlier.shape = NA , alpha=0.3, aes(fill = TYPE )) +
  geom_jitter(aes(color = TYPE ), size=1, alpha=0.5) +
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
  labs(x=NULL)+
  stat_boxplot(geom ='errorbar')+
  #scale_fill_manual(values = custom_colors) +
  #scale_color_manual(values = custom_colors) +
  scale_y_continuous(labels = function(x) 
    gsub("-", "\u2212", as.character(x))) #add long "-" signs

TYPE_FULL_plot
#save
png("C:/Users/Ieva/rprojects/outputs_all/NP_FULL_20260114.png", width = 1200, height = 800, res = 100)  # width/height in pixels, res = resolution
TYPE_FULL_plot
dev.off()

#remove ENDOMETRIAL CANCER###################################
KN_LIQUID_no_ENDOMETRIAL <- LIQUID_DF %>%
  filter(TYPE != "ENDOMETRIAL CANCER")

#test normalcy of my variables
results_normalcy2 <- lapply(DATA_names, function(v) {
  by(KN_LIQUID_no_ENDOMETRIAL[[v]], KN_LIQUID_no_ENDOMETRIAL$TYPE, function(x) {
    if (length(na.omit(x)) >= 3)
      shapiro.test(x)$p.value
    else NA
  })
})
names(results_normalcy2) <- DATA_names
results_normalcy2 #hes1 NP not normal IN BENIGN

# Loop through all variables in DATA_names
variance_results2 <- lapply(DATA_names, function(v) {
  # Select variable and TYPE
  df <- KN_LIQUID_no_ENDOMETRIAL[, c(v, "TYPE")]
  
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
variance_results2 <- bind_rows(variance_results2)

# View the results
variance_results2 #not equal in DLL1 urine

#NP ANOVA - NO ENDOMETRIAL ###############
# CTNNB1
anova_ctnnb1 <- aov(CTNNB1_NP ~ TYPE, data = KN_LIQUID_no_ENDOMETRIAL)
summary(anova_ctnnb1) #not significant
#DLL1
anova_dll1 <- aov(DLL1_NP ~ TYPE, data = KN_LIQUID_no_ENDOMETRIAL)
summary(anova_dll1) #not significant
#NOTCH2
anova_notch2 <- aov(NOTCH2_NP ~ TYPE, data = KN_LIQUID_no_ENDOMETRIAL) #significant
summary(anova_notch2) #not significant
#ANOVA (for not normal data) - HES1
#HES1
kruskal.test(HES1_NP ~ TYPE, data = KN_LIQUID_no_ENDOMETRIAL) #significant 0.03852
pairwise.wilcox.test(KN_LIQUID_no_ENDOMETRIAL$HES1_NP,
                     KN_LIQUID_no_ENDOMETRIAL$TYPE, p.adjust.method = "bonferroni", paired = FALSE) 
#better post-hoc
dunnTest(DLL1_NP ~ TYPE, data = KN_LIQUID_no_ENDOMETRIAL, method = "bonferroni")
#not significant
#STAGE##########################
#make an endometrial cancer df
KN_ENDOMETRIAL <- LIQUID_DF %>%
  filter(TYPE == "ENDOMETRIAL CANCER")
#chek stage
table(KN_ENDOMETRIAL$Stage, KN_ENDOMETRIAL$TYPE)
table(KN_LIQUID_no_ENDOMETRIAL$Stage, KN_LIQUID_no_ENDOMETRIAL$TYPE)
table(KN_ENDOMETRIAL$Stage_simple)
table(KN_LIQUID_no_ENDOMETRIAL$Stage_simple)
#fix stage simple - make factors
KN_LIQUID_no_ENDOMETRIAL$Stage_simple <- as.factor(KN_LIQUID_no_ENDOMETRIAL$Stage_simple)

KN_ENDOMETRIAL$Stage_simple <- as.factor(KN_ENDOMETRIAL$Stage_simple)
#NORMALCY KN
normality_results <- KN_LIQUID_no_ENDOMETRIAL %>%
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

# View
normality_results #all good
#NORMALCY ENDOMETRAIL
normality_results_endometrial <- KN_ENDOMETRIAL %>%
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

# View
normality_results_endometrial #all good

#variance KN
variance_results_4 <- lapply(DATA, function(v) {
  # Select variable and STAGE
  df <- KN_LIQUID_no_ENDOMETRIAL[, c(v, "Stage_simple")]
  
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
variance_results_4 #not equal in DLL1 urine

#variance ENDOMETRIAL
variance_results_5 <- lapply(DATA, function(v) {
  # Select variable and STAGE
  df <- KN_ENDOMETRIAL[, c(v, "Stage_simple")]
  
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
variance_results_5 <- bind_rows(variance_results_5)

# View the results
variance_results_5 #notch2 not equal variances

##STAGE ANOVA KN#######################################
#remove na values
STAGE_KN_DF <- KN_LIQUID_no_ENDOMETRIAL %>%
  filter(!is.na(Stage_simple))
STAGE_ENDO_DF <- KN_ENDOMETRIAL %>%
  filter(!is.na(Stage_simple))
table(STAGE_ENDO_DF$Stage_simple) #10 VS 2 not enough to make comparissons
#KN ANOVA:
# CTNNB1
anova_ctnnb1stage <- aov(CTNNB1_NP ~ Stage_simple, data = STAGE_KN_DF)
summary(anova_ctnnb1stage) #not significant
#DLL1
anova_dll1stage <- aov(DLL1_NP ~ Stage_simple, data = STAGE_KN_DF)
summary(anova_dll1stage) #not significant (but almost significant)
TukeyHSD(anova_dll1stage) #significant: 3 vs 2, bet tik 2 vnt stage2
#NOTCH2
anova_notch2stage <- aov(NOTCH2_NP ~ Stage_simple, data = STAGE_KN_DF) 
summary(anova_notch2stage) #not significant (but almost significant)
TukeyHSD(anova_notch2stage)#not significant 
#HES1
anova_hes1stage <- aov(HES1_NP ~ Stage_simple, data = STAGE_KN_DF) 
summary(anova_hes1stage) #not significant

#GRADE #############################################
table(KN_ENDOMETRIAL$Grade, useNA = "a") #7 vs 5 more equal
table(KN_LIQUID_no_ENDOMETRIAL$Grade, useNA = "a") #needs fixing
KN_LIQUID_no_ENDOMETRIAL <- KN_LIQUID_no_ENDOMETRIAL %>%
  mutate(
    Grade_clean = case_when(
      Grade %in% c("GB", "GL") ~ NA_character_,      # set GB and GL to NA
      Grade %in% c("G1", "G1&G1", "G2&G1") ~ "G1",  # merge into one
      TRUE ~ Grade                                   # keep others as-is (G2, G3)
    )
  )
table(KN_LIQUID_no_ENDOMETRIAL$Grade_clean, useNA = "a") #7 vs 50

#NORMALCY KN grade
normality_results2 <- KN_LIQUID_no_ENDOMETRIAL %>%
  pivot_longer(
    cols = all_of(DATA),
    names_to = "variable",
    values_to = "value"
  ) %>%
  group_by(Grade_clean, variable) %>%
  summarise(
    n = sum(!is.na(value)),
    shapiro_p = if (n >= 3 & n <= 5000) shapiro.test(value)$p.value else NA_real_,
    .groups = "drop"
  ) %>%
  mutate(
    normal = ifelse(shapiro_p > 0.05, TRUE, FALSE)
  )

# View
normality_results2 #all good

#NORMALCY ENDOMETRAIL grade
normality_results_endometrial2 <- KN_ENDOMETRIAL %>%
  pivot_longer(
    cols = all_of(DATA),
    names_to = "variable",
    values_to = "value"
  ) %>%
  group_by(Grade, variable) %>%
  summarise(
    n = sum(!is.na(value)),
    shapiro_p = if (n >= 3 & n <= 5000) shapiro.test(value)$p.value else NA_real_,
    .groups = "drop"
  ) %>%
  mutate(
    normal = ifelse(shapiro_p > 0.05, TRUE, FALSE)
  )

# View
normality_results_endometrial2 #all good
#variance KN grade
variance_results_6 <- lapply(DATA, function(v) {
  # Select variable and STAGE
  df <- KN_LIQUID_no_ENDOMETRIAL[, c(v, "Grade_clean")]
  
  # Remove NA
  df <- df[!is.na(df[[v]]), ]
  
  # Only test if at least 2 groups have data
  if(length(unique(df$Grade_clean)) > 1) {
    test <- leveneTest(df[[v]] ~ df$Grade_clean, center = median)
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
variance_results_6 <- bind_rows(variance_results_6)

# View the results
variance_results_6 #all good
#variance ENDOMETRIAL grade
variance_results_7 <- lapply(DATA, function(v) {
  # Select variable and STAGE
  df <- KN_ENDOMETRIAL[, c(v, "Grade")]
  
  # Remove NA
  df <- df[!is.na(df[[v]]), ]
  
  # Only test if at least 2 groups have data
  if(length(unique(df$Grade)) > 1) {
    test <- leveneTest(df[[v]] ~ df$Grade, center = median)
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
variance_results_7 <- bind_rows(variance_results_7)

# View the results
variance_results_7 #all good

##grade t test KN##################
KN_LIQUID_no_ENDOMETRIAL$Grade_clean <- as.factor(KN_LIQUID_no_ENDOMETRIAL$Grade_clean)
t.test(NOTCH2_NP ~ Grade_clean,
       data = KN_LIQUID_no_ENDOMETRIAL,
       var.equal = TRUE)
t.test(DLL1_NP ~ Grade_clean,
       data = KN_LIQUID_no_ENDOMETRIAL,
       var.equal = TRUE)
t.test(CTNNB1_NP ~ Grade_clean,
       data = KN_LIQUID_no_ENDOMETRIAL,
       var.equal = TRUE)
t.test(HES1_NP~ Grade_clean,
       data = KN_LIQUID_no_ENDOMETRIAL,
       var.equal = TRUE) #all not significant

##grade t test ENDOMETRIAL##################
KN_ENDOMETRIAL$Grade <- as.factor(KN_ENDOMETRIAL$Grade)
t.test(NOTCH2_NP ~ Grade,
       data = KN_ENDOMETRIAL,
       var.equal = TRUE)
t.test(DLL1_NP ~ Grade,
       data = KN_ENDOMETRIAL,
       var.equal = TRUE)
t.test(CTNNB1_NP ~ Grade,
       data = KN_ENDOMETRIAL,
       var.equal = TRUE)
t.test(HES1_NP~ Grade,
       data = KN_ENDOMETRIAL,
       var.equal = TRUE) #all not significant

#AGE################################################
#make only the OC df
#make an endometrial cancer df
KN_OVARIAN_ONLY <- LIQUID_DF %>%
  filter(TYPE %in% c("HGSOC", "OTHER"))
KN_non_Cancer  <- LIQUID_DF %>%
  filter(TYPE %in% c("BENIGN", "RSS"))

shapiro.test(KN_ENDOMETRIAL$Age)
shapiro.test(KN_OVARIAN_ONLY$Age)

shapiro.test(KN_ENDOMETRIAL$NOTCH2_NP) #not normal
shapiro.test(KN_ENDOMETRIAL$DLL1_NP)
shapiro.test(KN_ENDOMETRIAL$CTNNB1_NP)
shapiro.test(KN_ENDOMETRIAL$HES1_NP)

shapiro.test(KN_OVARIAN_ONLY$NOTCH2_NP)
shapiro.test(KN_OVARIAN_ONLY$DLL1_NP)
shapiro.test(KN_OVARIAN_ONLY$CTNNB1_NP) #NOT NORMAL
shapiro.test(KN_OVARIAN_ONLY$HES1_NP)

shapiro.test(KN_non_Cancer$NOTCH2_NP)
shapiro.test(KN_non_Cancer$DLL1_NP)
shapiro.test(KN_non_Cancer$CTNNB1_NP)
shapiro.test(KN_non_Cancer$HES1_NP)

#age correlation NON KN ###############################
cor.test(KN_non_Cancer$NOTCH2_NP, KN_non_Cancer$Age, method = "pearson")#0.02814
cor.test(KN_non_Cancer$DLL1_NP, KN_non_Cancer$Age, method = "pearson")
cor.test(KN_non_Cancer$HES1_NP, KN_non_Cancer$Age, method = "pearson") #0.002786
cor.test(KN_non_Cancer$CTNNB1_NP, KN_non_Cancer$Age, method = "pearson")

#plot HES1
cor_val <- cor(KN_non_Cancer$HES1_NP,
               KN_non_Cancer$Age,
               method = "pearson", use = "complete.obs")

# Compute correlation
cor_test <- cor.test(KN_non_Cancer$HES1_NP,
                     KN_non_Cancer$Age,
                     method = "pearson")

r_val <- round(cor_test$estimate, 3)      # correlation coefficient
p_val <- signif(cor_test$p.value, 3)      # p-value

# Scatter plot with regression line and text
ggplot(KN_non_Cancer, aes(x = Age, y = HES1_NP)) +
  geom_point(color = "blue") +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  annotate("text",
           x = min(KN_non_Cancer$Age, na.rm = TRUE),
           y = max(KN_non_Cancer$HES1_NP, na.rm = TRUE),
           label = paste0("r = ", r_val, ", p = ", p_val),
           hjust = 0, vjust = 1, size = 5) +
  theme_minimal() +
  labs(title = "HES1_NP vs Age, non cancer group NP", x = "Age", y = "HES1_NP")

# plot NOTCH2
# Compute correlation
cor_test2 <- cor.test(KN_non_Cancer$NOTCH2_NP,
                      KN_non_Cancer$Age,
                     method = "pearson")

r_val2 <- round(cor_test2$estimate, 3)      # correlation coefficient
p_val2 <- signif(cor_test2$p.value, 3)      # p-value

# Scatter plot with regression line and text
ggplot(KN_non_Cancer, aes(x = Age, y = NOTCH2_NP)) +
  geom_point(color = "blue") +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  annotate("text",
           x = min(KN_non_Cancer$Age, na.rm = TRUE),
           y = max(KN_non_Cancer$DLL1_NP, na.rm = TRUE),
           label = paste0("r = ", r_val2, ", p = ", p_val2),
           hjust = 0, vjust = 1, size = 5) +
  theme_minimal() +
  labs(title = "NOTCH2_NP vs Age, non cancer NP", x = "Age", y = "NOTCH2_NP")

#age correlation ENDOMETRIAL ###############################
cor.test(KN_ENDOMETRIAL$NOTCH2_NP, KN_ENDOMETRIAL$Age, method = "spearman")
cor.test(KN_ENDOMETRIAL$DLL1_NP, KN_ENDOMETRIAL$Age, method = "pearson") 
cor.test(KN_ENDOMETRIAL$HES1_NP, KN_ENDOMETRIAL$Age, method = "pearson")
cor.test(KN_ENDOMETRIAL$CTNNB1_NP, KN_ENDOMETRIAL$Age, method = "pearson")

#age correlation OC ###############################
cor.test(KN_OVARIAN_ONLY$NOTCH2_NP, KN_OVARIAN_ONLY$Age, method = "pearson")
cor.test(KN_OVARIAN_ONLY$DLL1_NP, KN_OVARIAN_ONLY$Age, method = "pearson") 
cor.test(KN_OVARIAN_ONLY$HES1_NP, KN_OVARIAN_ONLY$Age, method = "pearson")
cor.test(KN_OVARIAN_ONLY$CTNNB1_NP, KN_OVARIAN_ONLY$Age, method = "pearson")

#SURVIVAL FULL KN###################################
#make df of only cases with survival
ALL_SURV_LIQUID <- LIQUID_DF %>%
  filter(!is.na(OS), !is.na(STATUS)) #58 observations
table(ALL_SURV_LIQUID$STATUS, useNA = "a") #29 dead, 67 NOT
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
#make only HGSOC SURV DF
HGSOC_SURV_LIQUID <- ALL_SURV_LIQUID %>%
  filter(TYPE == "HGSOC") 
table(HGSOC_SURV_LIQUID$STATUS, useNA = "a") #20 dead
#make only oc surv df
OC_SURV_LIQUID <- ALL_SURV_LIQUID %>%
  filter(TYPE %in%c("HGSOC", "OTHER")) 
table(OC_SURV_LIQUID$STATUS, useNA = "a") #23 dead
#make only ENDOMETRIAL surv df
ENDOMETRIAL_TUMOR_SURV_LIQUID <- ALL_SURV_LIQUID %>%
  filter(TYPE %in%c("ENDOMETRIAL CANCER")) 
table(ENDOMETRIAL_TUMOR_SURV_LIQUID$STATUS, useNA = "a") #2 dead, 10 not


##NP, ALL KN######################################################################
##univariable cox, KN, NP all ##################
cox_model_notch2_np <- coxph(Surv(OS, STATUS) ~ NOTCH2_NP, data = OC_SURV_LIQUID)
summary(cox_model_notch2_np)

cox_model_CTNNB1_np <- coxph(Surv(OS, STATUS) ~ CTNNB1_NP, data = OC_SURV_LIQUID)
summary(cox_model_CTNNB1_np)

cox_model_DLL1_np <- coxph(Surv(OS, STATUS) ~ DLL1_NP, data = OC_SURV_LIQUID)
summary(cox_model_DLL1_np)# P = 0.0695 .

cox_model_HES1_np <- coxph(Surv(OS, STATUS) ~ HES1_NP, data = OC_SURV_LIQUID)
summary(cox_model_HES1_np) #1.304 , hR 1.318 (1.022 - 1.664)


## Fit survival curves#################
ggsurvplot(survfit(Surv(OS, STATUS) ~ CTNNB1_NP_f, data = OC_SURV_LIQUID), 
           data = OC_SURV_LIQUID, pval = TRUE,
           title="Overall survival by CTNNB1 expression in NP, all OC cases")

ggsurvplot(survfit(Surv(OS, STATUS) ~ NOTCH2_NP_f, data = OC_SURV_LIQUID), 
           data = OC_SURV_LIQUID, pval = TRUE,
           title="Overall survival by NOTCH2 expression in NP, all OC cases")

ggsurvplot(survfit(Surv(OS, STATUS) ~ DLL1_NP_f, data = OC_SURV_LIQUID), 
           data = OC_SURV_LIQUID, pval = TRUE,
           title="Overall survival by DLL1 expression in NP, all OC cases")

ggsurvplot(survfit(Surv(OS, STATUS) ~ HES1_NP_f, data = OC_SURV_LIQUID), 
           data = OC_SURV_LIQUID, pval = TRUE,
           title="Overall survival by HES1 expression in NP, all OC cases")

##median survival NP, KN, TISSUE dataset###########################

hes1_fit <- survfit(Surv(OS, STATUS) ~ HES1_NP_f, data = OC_SURV_LIQUID)
summary(hes1_fit)$table #gives months
summary(hes1_fit, times = 12)$surv #gives survival prob at 1 year
summary(hes1_fit, times = 36)$surv #gives survival prob at 3 yrs
summary(hes1_fit, times = 60)$surv #gives survival prob at 5 yrs
#what is median survival overall?
fit_all <- survfit(Surv(OS, STATUS) ~ 1, data = OC_SURV_LIQUID)
summary(fit_all)$table #not reached
ggsurvplot(
  hes1_fit,
  data = OC_SURV_LIQUID,
  pval = TRUE,
  risk.table = TRUE, 
  surv.median.line = "hv",   # adds horizontal & vertical median lines
  title = "Overall survival by HES1 expression in NP, OC TUMOR cohort"
)
##NP, ENDOMETRIAL CANCER ##################################
## Fit survival curves ENDOMETRIAL CANCER#################
ggsurvplot(survfit(Surv(OS, STATUS) ~ CTNNB1_NP_f, data = ENDOMETRIAL_TUMOR_SURV_LIQUID), 
           data = ENDOMETRIAL_TUMOR_SURV_LIQUID, pval = TRUE,
           title="Overall survival by CTNNB1 expression in NP, ENDOMETRIAL CANCER cases")

ggsurvplot(survfit(Surv(OS, STATUS) ~ NOTCH2_NP_f, data = ENDOMETRIAL_TUMOR_SURV_LIQUID), 
           data = ENDOMETRIAL_TUMOR_SURV_LIQUID, pval = TRUE,
           title="Overall survival by NOTCH2 expression in NP, ENDOMETRIAL CANCER cases")

ggsurvplot(survfit(Surv(OS, STATUS) ~ DLL1_NP_f, data = ENDOMETRIAL_TUMOR_SURV_LIQUID), 
           data = ENDOMETRIAL_TUMOR_SURV_LIQUID, pval = TRUE,
           title="Overall survival by DLL1 expression in NP, ENDOMETRIAL CANCER cases")

ggsurvplot(survfit(Surv(OS, STATUS) ~ HES1_NP_f, data = ENDOMETRIAL_TUMOR_SURV_LIQUID), 
           data = ENDOMETRIAL_TUMOR_SURV_LIQUID, pval = TRUE,
           title="Overall survival by HES1 expression in NP, ENDOMETRIAL CANCER cases")
#ROC ENDOMETRIAL VS BENIGN ########################################
#make ENDOMETRIAL VS BENIGN DF
ENDO_BENIGN_DF<- LIQUID_DF %>%
  filter(TYPE%in% c("BENIGN", "ENDOMETRIAL CANCER"))%>% #filter for right samples
  dplyr::select("TYPE", "NOTCH2_NP","CTNNB1_NP","DLL1_NP","HES1_NP"
       )
ENDO_BENIGN_DF$TYPE <- droplevels(ENDO_BENIGN_DF$TYPE) #drop unused
ENDO_BENIGN_DF$TYPE  <- relevel(ENDO_BENIGN_DF$TYPE , ref = "BENIGN") #set control
roc_results_np_endo_benign<- lapply(c("NOTCH2_NP","CTNNB1_NP","DLL1_NP","HES1_NP"), function(col) {
  roc(response = ENDO_BENIGN_DF$TYPE, predictor = ENDO_BENIGN_DF[[col]])})
names(roc_results_np_endo_benign) <- c("NOTCH2_NP","CTNNB1_NP","DLL1_NP","HES1_NP")
roc_results_np_endo_benign
#extract the aucs
auc_values_np_endo_benign <- sapply(roc_results_np_endo_benign, function(roc_obj) {auc(roc_obj)})
auc_values_np_endo_benign #extracted aucs
#roc figure 
roc_plot <- function() {
  par(pty = "s") #sets square
  plot.roc(roc_results_np_endo_benign[["NOTCH2_NP"]], print.auc = F, col = "#dcbeff",
           cex.main=0.8, 
           main ="Gerybinių pokyčių atskyrimas nuo gimdos kūno vėžio",
           xlab = "1 - Specifiškumas", 
           ylab = "Jautrumas", 
           legacy.axes = T) #title
  lines(roc_results_np_endo_benign[["CTNNB1_NP"]], col = "#911eb4", lwd =2) 
  lines(roc_results_np_endo_benign[["DLL1_NP"]], col ="#ffd8b1", lwd =2) 
  lines(roc_results_np_endo_benign[["HES1_NP"]], col = "#42d4f4", lwd =2) 
  legend("bottomright", legend = c( expression(italic("NOTCH2")),
                                    expression(italic("CTNNB1")),
                                    expression(italic("DLL1")), 
                                    expression(italic("HES1"))
  ),
  
  col = c("#dcbeff", "#911eb4", "#ffd8b1", "#42d4f4"), lty = 1, 
  cex = 0.7, lwd =3)
}
#plot
roc_plot()
#roc table
coords_results_np_endo_benign<- lapply(roc_results_np_endo_benign, function(roc_obj) {
  coords(roc_obj, "best", ret = c("threshold", "accuracy", "sensitivity",
                                  "specificity"),
         transpose = FALSE)
})
coords_results_np_endo_benign
#create df
results_np_endo_benign<- data.frame(
  Predictor = c("NOTCH2_NP","CTNNB1_NP","DLL1_NP","HES1_NP"),
    AUC = auc_values_np_endo_benign,
  do.call(rbind,coords_results_np_endo_benign) 
)
results_np_endo_benign
#lithuanize it 
colnames(results_np_endo_benign) <- c("Biožymuo", "plotas po kreive", "slenkstinė vertė", 
                             "tikslumas", "jautrumas", "specifiškumas")
rownames(results_np_endo_benign) <- c("NOCTH2", "CTNNB1", "DLL1", "HES1")
results_np_endo_benign$Biožymuo <- c("NOCTH2", "CTNNB1", "DLL1", "HES1")
#nice formating of the Table metrics for ROC OC
gt_table_np_endo_ben <- results_np_endo_benign %>%
  gt() %>%
  tab_header(
    title = "ROC kriterijai", 
    subtitle = "Gerybinių pokyčių atskyrimas nuo gimdos kūno vėžio") %>%
  fmt_number(
    columns = everything(),
    decimals = 3
  ) %>%
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_body(columns = vars(Biožymuo))
  )
#show
gt_table_np_endo_ben
#ROC ENDOMETRIAL vs HGSOC###################################################
#make ENDOMETRIAL VS hgsoc DF
ENDO_hgsoc_DF<- LIQUID_DF %>%
  filter(TYPE%in% c("HGSOC", "ENDOMETRIAL CANCER"))%>% #filter for right samples
  dplyr::select("TYPE", "NOTCH2_NP","CTNNB1_NP","DLL1_NP","HES1_NP"
  )
ENDO_hgsoc_DF$TYPE <- droplevels(ENDO_hgsoc_DF$TYPE) #drop unused
ENDO_hgsoc_DF$TYPE  <- relevel(ENDO_hgsoc_DF$TYPE , ref = "ENDOMETRIAL CANCER") #set control
roc_results_np_endo_hgsoc<- lapply(c("NOTCH2_NP","CTNNB1_NP","DLL1_NP","HES1_NP"), function(col) {
  roc(response = ENDO_hgsoc_DF$TYPE, predictor = ENDO_hgsoc_DF[[col]])})
names(roc_results_np_endo_hgsoc) <- c("NOTCH2_NP","CTNNB1_NP","DLL1_NP","HES1_NP")
roc_results_np_endo_hgsoc
#extract the aucs
auc_values_np_endo_hgsoc <- sapply(roc_results_np_endo_hgsoc, function(roc_obj) {auc(roc_obj)})
auc_values_np_endo_hgsoc #extracted aucs
#roc figure 
roc_plot2 <- function() {
  par(pty = "s") #sets square
  plot.roc(roc_results_np_endo_hgsoc[["NOTCH2_NP"]], print.auc = F, col = "#dcbeff",
           cex.main=0.8, 
           main ="HGSOC atskyrimas nuo gimdos kūno vėžio",
           xlab = "1 - Specifiškumas", 
           ylab = "Jautrumas", 
           legacy.axes = T) #title
  lines(roc_results_np_endo_hgsoc[["CTNNB1_NP"]], col = "#911eb4", lwd =2) 
  lines(roc_results_np_endo_hgsoc[["DLL1_NP"]], col ="#ffd8b1", lwd =2) 
  lines(roc_results_np_endo_hgsoc[["HES1_NP"]], col = "#42d4f4", lwd =2) 
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
coords_results_np_endo_hgsoc<- lapply(roc_results_np_endo_hgsoc, function(roc_obj) {
  coords(roc_obj, "best", ret = c("threshold", "accuracy", "sensitivity",
                                  "specificity"),
         transpose = FALSE)
})
coords_results_np_endo_hgsoc
#create df
results_np_endo_hgsoc<- data.frame(
  Predictor = c("NOTCH2_NP","CTNNB1_NP","DLL1_NP","HES1_NP"),
  AUC = auc_values_np_endo_hgsoc,
  do.call(rbind,coords_results_np_endo_hgsoc) 
)
results_np_endo_hgsoc
#lithuanize it 
colnames(results_np_endo_hgsoc) <- c("Biožymuo", "plotas po kreive", "slenkstinė vertė", 
                                     "tikslumas", "jautrumas", "specifiškumas")
rownames(results_np_endo_hgsoc) <- c("NOCTH2", "CTNNB1", "DLL1", "HES1")
results_np_endo_hgsoc$Biožymuo <- c("NOCTH2", "CTNNB1", "DLL1", "HES1")
#nice formating of the Table metrics for ROC OC
gt_table_np_endo_hsgoc<- results_np_endo_hgsoc %>%
  gt() %>%
  tab_header(
    title = "ROC kriterijai", 
    subtitle = "HSGOC atskyrimas nuo gimdos kūno vėžio") %>%
  fmt_number(
    columns = everything(),
    decimals = 3
  ) %>%
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_body(columns = vars(Biožymuo))
  )
#show
gt_table_np_endo_hsgoc
#ROC hgsoc VS BENIGN ########################################
#make hgsoc VS BENIGN DF
HGSOC_BENIGN_DF<- LIQUID_DF %>%
  filter(TYPE%in% c("BENIGN", "HGSOC"))%>% #filter for right samples
  dplyr::select("TYPE", "NOTCH2_NP","CTNNB1_NP","DLL1_NP","HES1_NP"
  )
HGSOC_BENIGN_DF$TYPE <- droplevels(HGSOC_BENIGN_DF$TYPE) #drop unused
HGSOC_BENIGN_DF$TYPE  <- relevel(HGSOC_BENIGN_DF$TYPE , ref = "BENIGN") #set control
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
           main ="Gerybinių pokyčių atskyrimas nuo HGSOC",
           xlab = "1 - Specifiškumas", 
           ylab = "Jautrumas", 
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
#lithuanize it 
colnames(results_np_HGSOC_benign) <- c("Biožymuo", "plotas po kreive", "slenkstinė vertė", 
                                       "tikslumas", "jautrumas", "specifiškumas")
rownames(results_np_HGSOC_benign) <- c("NOCTH2", "CTNNB1", "DLL1", "HES1")
results_np_HGSOC_benign$Biožymuo <- c("NOCTH2", "CTNNB1", "DLL1", "HES1")
#nice formating of the Table metrics for ROC OC
gt_table_np_HGSOC_ben <- results_np_HGSOC_benign %>%
  gt() %>%
  tab_header(
    title = "ROC kriterijai", 
    subtitle = "Gerybinių pokyčių atskyrimas nuo HGSOC") %>%
  fmt_number(
    columns = everything(),
    decimals = 3
  ) %>%
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_body(columns = vars(Biožymuo))
  )
#show
gt_table_np_HGSOC_ben
#ROC hgsoc VS others ########################################
#make hgsoc VS OTHER DF
HGSOC_OTHER_DF<- LIQUID_DF %>%
  filter(TYPE%in% c("OTHER", "HGSOC"))%>% #filter for right samples
  dplyr::select("TYPE", "NOTCH2_NP","CTNNB1_NP","DLL1_NP","HES1_NP"
  )
HGSOC_OTHER_DF$TYPE <- droplevels(HGSOC_OTHER_DF$TYPE) #drop unused
HGSOC_OTHER_DF$TYPE  <- relevel(HGSOC_OTHER_DF$TYPE , ref = "OTHER") #set control
roc_results_np_HGOSC_OTHER<- lapply(c("NOTCH2_NP","CTNNB1_NP","DLL1_NP","HES1_NP"), function(col) {
  roc(response = HGSOC_OTHER_DF$TYPE, predictor = HGSOC_OTHER_DF[[col]])})
names(roc_results_np_HGOSC_OTHER) <- c("NOTCH2_NP","CTNNB1_NP","DLL1_NP","HES1_NP")
roc_results_np_HGOSC_OTHER
#extract the aucs
auc_values_np_HGSOC_OTHER <- sapply(roc_results_np_HGOSC_OTHER, function(roc_obj) {auc(roc_obj)})
auc_values_np_HGSOC_OTHER #extracted aucs
#roc figure 
roc_plot4 <- function() {
  par(pty = "s") #sets square
  plot.roc(roc_results_np_HGOSC_OTHER[["NOTCH2_NP"]], print.auc = F, col = "#dcbeff",
           cex.main=0.8, 
           main ="Ne HGSOC tipo KV  atskyrimas nuo HGSOC",
           xlab = "1 - Specifiškumas", 
           ylab = "Jautrumas", 
           legacy.axes = T) #title
  lines(roc_results_np_HGOSC_OTHER[["CTNNB1_NP"]], col = "#911eb4", lwd =2) 
  lines(roc_results_np_HGOSC_OTHER[["DLL1_NP"]], col ="#ffd8b1", lwd =2) 
  lines(roc_results_np_HGOSC_OTHER[["HES1_NP"]], col = "#42d4f4", lwd =2) 
  legend("bottomright", legend = c( expression(italic("NOTCH2")),
                                    expression(italic("CTNNB1")),
                                    expression(italic("DLL1")), 
                                    expression(italic("HES1"))
  ),
  
  col = c("#dcbeff", "#911eb4", "#ffd8b1", "#42d4f4"), lty = 1, 
  cex = 0.7, lwd =3)
}
#plot
roc_plot4()
#roc table
coords_results_np_HGSOC_OTHER<- lapply(roc_results_np_HGOSC_OTHER, function(roc_obj) {
  coords(roc_obj, "best", ret = c("threshold", "accuracy", "sensitivity",
                                  "specificity"),
         transpose = FALSE)
})
coords_results_np_HGSOC_OTHER
#create df
results_np_HGSOC_OTHER<- data.frame(
  Predictor = c("NOTCH2_NP","CTNNB1_NP","DLL1_NP","HES1_NP"),
  AUC = auc_values_np_HGSOC_OTHER,
  do.call(rbind,coords_results_np_HGSOC_OTHER) 
)
results_np_HGSOC_OTHER
#lithuanize it 
colnames(results_np_HGSOC_OTHER) <- c("Biožymuo", "plotas po kreive", "slenkstinė vertė", 
                                      "tikslumas", "jautrumas", "specifiškumas")
rownames(results_np_HGSOC_OTHER) <- c("NOCTH2", "CTNNB1", "DLL1", "HES1")
results_np_HGSOC_OTHER$Biožymuo <- c("NOCTH2", "CTNNB1", "DLL1", "HES1")
#nice formating of the Table metrics for ROC OC
gt_table_np_HGSOC_OTHER <- results_np_HGSOC_OTHER %>%
  gt() %>%
  tab_header(
    title = "ROC kriterijai", 
    subtitle = "Ne HGSOC tipo KV atskyrimas nuo HGSOC") %>%
  fmt_number(
    columns = everything(),
    decimals = 3
  ) %>%
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_body(columns = vars(Biožymuo))
  )
#show
gt_table_np_HGSOC_OTHER
