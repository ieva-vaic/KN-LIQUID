#KN-  liquid 2026 01 15, 16
#CASES WITH TUMOR DATA
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
library(multcomp)
library(ggpubr) 
#read RDS
LIQUID_DF <- readRDS("C:/Users/Ieva/rprojects/OTHER DATA/KN_LIQUID/liquid_20260114.RDS")
#make tumor data df#################################
#only DLL1 tumor had empty variables, thus use any other tumor column
TUMOR_df <- LIQUID_DF %>%
  filter(!is.na(NOTCH2_TUMOR))
TUMOR_df$TYPE <- droplevels(TUMOR_df$TYPE)
TUMOR_df$TYPE_tumor <- fct_recode(
  TUMOR_df$TYPE,
  BENIGN = "RSS"
)
table(TUMOR_df$TYPE_tumor, useNA = "a")
#fix grade
TUMOR_df <- TUMOR_df %>%
  mutate(
    Grade = fct_recode(
      Grade,
      G1 = "G2&G1",
      G1 = "G1&G1"
    ),
    Grade = na_if(Grade, "GL"),    # make "GL" NA
    Grade = na_if(Grade, "GB")     # make "GB" NA
  )
TUMOR_df$Grade <- droplevels(TUMOR_df$Grade)
table(TUMOR_df$Grade)
#NORMALCY TEST for groups #######################################
#chek NP in tumor tissue cohort
shapiro_results1 <- TUMOR_df[, c(36,15:18)] %>%
  pivot_longer(cols = -TYPE_tumor, names_to = "gene", values_to = "value") %>%
  group_by(TYPE_tumor, gene) %>%
  summarise(p_value = shapiro.test(value)$p.value, .groups = "drop") %>%
  filter(p_value < 0.05)
shapiro_results1 #not normal for benign dll1 and OC ctnnb1
#chek NP variance in tumor tissue cohort
variance_results1 <- TUMOR_df[, c(36, 15:18)] %>%
  pivot_longer(cols = -TYPE_tumor,
               names_to = "gene",
               values_to = "value") %>%
  group_by(gene) %>%
  summarise(
    p_value = leveneTest(value ~ TYPE_tumor)$`Pr(>F)`[1],
    .groups = "drop"
  ) 
variance_results1 #all normal
#ANOVA #######################################
#normal: NOTCH2 and HES1
#NOTCH2
res_aovNOCTH2 <- aov( NOTCH2_NP  ~ TYPE_tumor,
                data = TUMOR_df
)
shapiro.test(res_aovNOCTH2$residuals) #normality is met
summary(res_aovNOCTH2) # p = 0.437
#HES1
res_aovHES1 <- aov( HES1_NP ~ TYPE_tumor,
                      data = TUMOR_df
)
shapiro.test(res_aovHES1$residuals) #normality is met
summary(res_aovHES1) # p = 0.037 
#TUCKEY posthoc
post_test_HES1 <- glht(res_aovHES1,
                         linfct = mcp(TYPE_tumor  = "Tukey"))

summary(post_test_HES1) #hgsoc benign 0.0277
#not normal: DLL1 and CTNNB1
#DLL1
kruskalDLL1 <- kruskal.test(DLL1_NP ~ TYPE_tumor,
                         data = TUMOR_df
)
kruskalDLL1 # p =  0.7061
#CTNNB1
kruskalCTNNB1 <- kruskal.test(CTNNB1_NP ~ TYPE_tumor,
                            data = TUMOR_df
)
kruskalCTNNB1 # p = 0.5458

#BOXPLOT NP, TUMOR COHORT#################################################
#Tribble
each.vs.ref_sig <-  tibble::tribble(
    ~group1, ~group2, ~p.adj,   ~y.position, ~variable,
    "BENIGN",   "HGSOC", 0.028, -2, "HES1_NP", 
  )
#boxplot 3 groups ########################################
#melt table for expression
Group3_table <- melt(TUMOR_df[, c(36,15:18)], id.vars="TYPE_tumor",
                     measure.vars=c("NOTCH2_NP", "DLL1_NP", "HES1_NP", "CTNNB1_NP"))

custom_colors <- c("HGSOC" = "deeppink","OTHER" = "lightpink", "BENIGN" = "lightblue") 
OC_plot <- ggplot(Group3_table, aes(x=TYPE_tumor , y=value, fill = variable)) +
  geom_boxplot( outlier.shape = NA , alpha=0.3, aes(fill = TYPE_tumor )) +
  geom_jitter(aes(color = TYPE_tumor ), size=1, alpha=0.5) +
  ylab(label = expression("Relative expression, normalized to " * italic("GAPDH"))) + 
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
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  scale_y_continuous(labels = function(x) 
    gsub("-", "\u2212", as.character(x))) #add long "-" signs

OC_plot + ggtitle("NP, in TUMOR cohort")

#FC
#FC 3 groups###########################################
expression_df <- TUMOR_df[, c("TYPE_tumor", c("NOTCH2_NP", "DLL1_NP", "HES1_NP", "CTNNB1_NP"), "Laboratorinis kodas")]
rownames(expression_df) <- expression_df$`Laboratorinis kodas`
expression_df <- expression_df[, -6]
#REMOVE NA
expression_df2 <- expression_df %>%
  filter(complete.cases(.))

exp_df2 <- expression_df2 %>%
  melt(id.vars="TYPE_tumor",  measure.vars=c("NOTCH2_NP", "DLL1_NP", "HES1_NP", "CTNNB1_NP")) %>%
  rename(gene = variable) %>%
  rename(expression = value)%>%
  rename(condition = TYPE_tumor) %>%
  group_by(condition, gene)%>%
  summarise(mean_expression = mean(expression, drop_na= T)) %>%
  spread(condition, mean_expression)

mean_expression_OC2 <- exp_df2 %>%
  mutate(fold_change_HB = log2(2^(`HGSOC` - `BENIGN`))) %>% #benign vs hgsoc 
  mutate(fold_change_HO = log2(2^(`HGSOC` - `OTHER`))) %>% #others vs hgsoc 
  mutate(fold_change_OB = log2(2^(`OTHER` - `BENIGN`))) #benign vs others
mean_expression_OC2 

#ROC###################################################################
#make groupings of diseases#################
OC_HGSOC_BENIGN<- TUMOR_df[c(TUMOR_df$TYPE_tumor != "OTHER"),] #51 cases left
OC_HGSOC_BENIGN$TYPE_tumor <- relevel(factor(OC_HGSOC_BENIGN$TYPE_tumor), ref = "BENIGN")

OC_HGSOC_OTHER<- TUMOR_df[c(TUMOR_df$TYPE_tumor != "BENIGN"),] 
OC_HGSOC_OTHER$TYPE_tumor <- relevel(factor(OC_HGSOC_OTHER$TYPE_tumor), ref = "OTHER")

OC_BENIGN_OTHER<- TUMOR_df[c(TUMOR_df$TYPE_tumor != "HGSOC"),] 
OC_BENIGN_OTHER$TYPE_tumor <- relevel(factor(OC_BENIGN_OTHER$TYPE_tumor), ref = "BENIGN")

#ROC np BH########################################################
roc_results_BH <- lapply(c("NOTCH2_NP", "DLL1_NP", "HES1_NP", "CTNNB1_NP"), function(col) {
  roc(response = OC_HGSOC_BENIGN$TYPE_tumor, predictor = OC_HGSOC_BENIGN[[col]])})
names(roc_results_BH) <- c("NOTCH2_NP", "DLL1_NP", "HES1_NP", "CTNNB1_NP")
roc_results_BH
#extract the aucs
auc_values_BH <- sapply(roc_results_BH, function(roc_obj) {auc(roc_obj)})
auc_values_BH #extracted aucs

roc_plot <- function() {
  par(pty = "s") #sets square
  plot.roc(roc_results_BH[["NOTCH2_NP"]], print.auc = F, col = "#dcbeff",
           cex.main=0.8, 
           main ="Gerybinių pokyčių atskyrimas nuo HGSOC, NP, audinių imtyje",
           xlab = "1 - Specifiškumas", 
           ylab = "Jautrumas", legacy.axes = T) #title
  lines(roc_results_BH[["CTNNB1_NP"]], col = "#911eb4", lwd =2) 
  lines(roc_results_BH[["DLL1_NP"]], col ="#ffd8b1", lwd =2) 
  lines(roc_results_BH[["HES1_NP"]], col = "#42d4f4", lwd =2) 
  legend("bottomright", legend = c( expression(italic("NOTCH2")),
                                    expression(italic("CTNNB1")),
                                    expression(italic("DLL1")), 
                                    expression(italic("HES1"))
  ),
  
  col = c("#dcbeff", "#911eb4", "#ffd8b1", "#42d4f4"), lty = 1, 
  cex = 0.5, lwd =3)
}
#plot
roc_plot()

#roc table BH################################
#get roc features
coords_results_BH<- lapply(roc_results_BH, function(roc_obj) {
  coords(roc_obj, "best", ret = c("threshold", "accuracy", "sensitivity",
                                  "specificity"),
         transpose = FALSE)
})
coords_results_BH
#create df
results_tumor<- data.frame(
  Predictor = names(roc_results_BH),
  AUC = auc_values_BH,
  do.call(rbind,coords_results_BH) 
)
#lithuanize it 
colnames(results_tumor) <- c("Biožymuo", "plotas po kreive", "slenkstinė vertė", 
                             "tikslumas", "jautrumas", "specifiškumas")
#nice formating of the Table metrics for ROC OC
gt_table_tumor <- results_tumor %>%
  gt() %>%
  tab_header(
    title = "ROC kriterijai", 
    subtitle = "Gerybinių pokyčių atskyrimas nuo HGSOC, NP, audinių imtyje") %>%
  fmt_number(
    columns = everything(),
    decimals = 3
  ) %>%
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_body(columns = vars(Biožymuo))
  )
#show
gt_table_tumor

#ROC np BO########################################################
roc_results_BO <- lapply(c("NOTCH2_NP", "DLL1_NP", "HES1_NP", "CTNNB1_NP"), function(col) {
  roc(response = OC_BENIGN_OTHER$TYPE_tumor, predictor = OC_BENIGN_OTHER[[col]])})
names(roc_results_BO) <- c("NOTCH2_NP", "DLL1_NP", "HES1_NP", "CTNNB1_NP")
roc_results_BO
#extract the aucs
auc_values_BO <- sapply(roc_results_BO, function(roc_obj) {auc(roc_obj)})
auc_values_BO #extracted aucs

roc_plot2 <- function() {
  par(pty = "s") #sets square
  plot.roc(roc_results_BO[["NOTCH2_NP"]], print.auc = F, col = "#dcbeff",
           cex.main=0.8, 
           main ="Gerybinių pokyčių atskyrimas nuo ne HGSOC KV, NP, audinių imtyje",
           xlab = "1 - Specifiškumas", 
           ylab = "Jautrumas", legacy.axes = T) #title
  lines(roc_results_BO[["CTNNB1_NP"]], col = "#911eb4", lwd =2) 
  lines(roc_results_BO[["DLL1_NP"]], col ="#ffd8b1", lwd =2) 
  lines(roc_results_BO[["HES1_NP"]], col = "#42d4f4", lwd =2) 
  legend("bottomright", legend = c( expression(italic("NOTCH2")),
                                    expression(italic("CTNNB1")),
                                    expression(italic("DLL1")), 
                                    expression(italic("HES1"))
  ),
  
  col = c("#dcbeff", "#911eb4", "#ffd8b1", "#42d4f4"), lty = 1, 
  cex = 0.5, lwd =3)
}
#plot
roc_plot2()

#roc table BO################################
#get roc features
coords_results_BO<- lapply(roc_results_BO, function(roc_obj) {
  coords(roc_obj, "best", ret = c("threshold", "accuracy", "sensitivity",
                                  "specificity"),
         transpose = FALSE)
})
coords_results_BO
#create df
results_tumor2<- data.frame(
  Predictor = names(roc_results_BO),
  AUC = auc_values_BO,
  do.call(rbind,coords_results_BO) 
)
#lithuanize it 
colnames(results_tumor2) <- c("Biožymuo", "plotas po kreive", "slenkstinė vertė", 
                              "tikslumas", "jautrumas", "specifiškumas")
#nice formating of the Table metrics for ROC OC
gt_table_tumor2 <- results_tumor2 %>%
  gt() %>%
  tab_header(
    title = "ROC kriterijai", 
    subtitle = "Gerybinių pokyčių atskyrimas nuo ne HGSOC KV, NP, audinių imtyje") %>%
  fmt_number(
    columns = everything(),
    decimals = 3
  ) %>%
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_body(columns = vars(Biožymuo))
  )
#show
gt_table_tumor2

#ROC np HO########################################################
roc_results_HO <- lapply(c("NOTCH2_NP", "DLL1_NP", "HES1_NP", "CTNNB1_NP"), function(col) {
  roc(response = OC_HGSOC_OTHER$TYPE_tumor, predictor = OC_HGSOC_OTHER[[col]])})
names(roc_results_HO) <- c("NOTCH2_NP", "DLL1_NP", "HES1_NP", "CTNNB1_NP")
roc_results_HO
#extract the aucs
auc_values_HO <- sapply(roc_results_HO, function(roc_obj) {auc(roc_obj)})
auc_values_HO #extracted aucs

roc_plot2 <- function() {
  par(pty = "s") #sets square
  plot.roc(roc_results_HO[["NOTCH2_NP"]], print.auc = F, col = "#dcbeff",
           cex.main=0.8, 
           main ="HGSOC atskyrimas nuo ne HGSOC KV, NP, audinių imtyje",
           xlab = "1 - Specifiškumas", 
           ylab = "Jautrumas", legacy.axes = T) #title
  lines(roc_results_HO[["CTNNB1_NP"]], col = "#911eb4", lwd =2) 
  lines(roc_results_HO[["DLL1_NP"]], col ="#ffd8b1", lwd =2) 
  lines(roc_results_HO[["HES1_NP"]], col = "#42d4f4", lwd =2) 
  legend("bottomright", legend = c( expression(italic("NOTCH2")),
                                    expression(italic("CTNNB1")),
                                    expression(italic("DLL1")), 
                                    expression(italic("HES1"))
  ),
  
  col = c("#dcbeff", "#911eb4", "#ffd8b1", "#42d4f4"), lty = 1, 
  cex = 0.5, lwd =3)
}
#plot
roc_plot2()

#roc table HO################################
#get roc features
coords_results_HO<- lapply(roc_results_HO, function(roc_obj) {
  coords(roc_obj, "best", best.method = "closest.topleft", ret = c("threshold", "accuracy", "sensitivity",
                                                                   "specificity"),
         transpose = FALSE)
})
coords_results_HO
#create df
results_tumor3<- data.frame(
  Predictor = names(roc_results_HO),
  AUC = auc_values_HO,
  do.call(rbind,coords_results_HO) 
)
#lithuanize it 
colnames(results_tumor3) <- c("Biožymuo", "plotas po kreive", "slenkstinė vertė", 
                              "tikslumas", "jautrumas", "specifiškumas")
#nice formating of the Table metrics for ROC OC
gt_table_tumor3 <- results_tumor3 %>%
  gt() %>%
  tab_header(
    title = "ROC kriterijai", 
    subtitle = "HGSOC atskyrimas nuo ne HGSOC KV, NP, audinių imtyje") %>%
  fmt_number(
    columns = everything(),
    decimals = 3
  ) %>%
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_body(columns = vars(Biožymuo))
  )
#show
gt_table_tumor3

#CLINICAL FEATURES, NP ##########################################
#STAGE normalcy, variance ####################
TUMOR_df$Stage_simple <- as.factor(TUMOR_df$Stage_simple)
#normalcy
#NORMALCY KN
normality_results <- TUMOR_df %>%
  pivot_longer(
    cols = all_of(c("NOTCH2_NP", "DLL1_NP", "HES1_NP", "CTNNB1_NP")),
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
normality_results #CTNNB1_NP stage 1 not normal

#variance KN
variance_results <- lapply(c("NOTCH2_NP", "DLL1_NP", "HES1_NP", "CTNNB1_NP"), function(v) {
  # Select variable and STAGE
  df <- TUMOR_df[, c(v, "Stage_simple")]
  
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
variance_results <- bind_rows(variance_results)
variance_results# all  normal

##STAGE ANOVA KN#######################################
#remove na values
STAGE_KN_DF <- TUMOR_df %>%
  filter(!is.na(Stage_simple))
#KN ANOVA:
# CTNNB1 not normal
kruskalCTNNB1stage <- kruskal.test(CTNNB1_NP ~ Stage_simple,
                              data = STAGE_KN_DF
)
kruskalCTNNB1stage # p = 0.7044
#DLL1
anova_dll1stage <- aov(DLL1_NP ~ Stage_simple, data = STAGE_KN_DF)
summary(anova_dll1stage) #0.155
TukeyHSD(anova_dll1stage) 
#NOTCH2
anova_notch2stage <- aov(NOTCH2_NP ~ Stage_simple, data = STAGE_KN_DF) 
summary(anova_notch2stage) #not significant 
#HES1
anova_hes1stage <- aov(HES1_NP ~ Stage_simple, data = STAGE_KN_DF) 
summary(anova_hes1stage) # 0.0452 *
#TUCKEY posthoc
post_test_HES1stage <- glht(anova_hes1stage,
                       linfct = mcp(Stage_simple  = "Tukey"))

summary(post_test_HES1stage) #not significant


##Stage groups boxplot#####################
each.vs.ref_sig2_stage <-  tibble::tribble(
  ~group1, ~group2, ~p.adj,   ~y.position, ~variable,
  "1",   "4", 0.0452, -2, "HES1_NP", #anova
)
Stage_table <- melt(STAGE_KN_DF[, c(35,15:18)], id.vars="Stage_simple",
                    measure.vars=c("NOTCH2_NP", "DLL1_NP", "HES1_NP", "CTNNB1_NP"))

custom_colors_stage <- c("4" = "deeppink","3" = "lightpink","2" = "lightblue", "1" = "blue") 
STAGE_plot <- ggplot(Stage_table, aes(x=Stage_simple , y=value, fill = variable)) +
  geom_boxplot( outlier.shape = NA , alpha=0.3, aes(fill = Stage_simple )) +
  geom_jitter(aes(color = Stage_simple ), size=1, alpha=0.5) +
  ylab(label = expression("Relative expression, normalized " * italic("GAPDH"))) + 
  facet_wrap(.~ variable, nrow = 2, scales = "free") +
  add_pvalue(each.vs.ref_sig2_stage, label = "p.adj") + #pvalue
  theme_minimal()+
  theme(
    strip.text.x = element_text(
      size = 12, face = "bold.italic"
    ),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5))+
  labs(x=NULL)+
  stat_boxplot(geom ='errorbar')+
  scale_fill_manual(values = custom_colors_stage) +
  scale_color_manual(values = custom_colors_stage) +
  scale_y_continuous(labels = function(x) 
    gsub("-", "\u2212", as.character(x))) #add long "-" signs

STAGE_plot
#GRADE normalcy, variance########################
#remove na values
GRADE_KN_DF <- TUMOR_df %>%
  filter(!is.na(Grade))
GRADE_KN_DF
#NORMALCY KN grade
normality_results2 <- GRADE_KN_DF %>%
  pivot_longer(
    cols = all_of(c("NOTCH2_NP", "DLL1_NP", "HES1_NP", "CTNNB1_NP")),
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
normality_results2 #all good

#variance KN grade
variance_results_6 <- lapply(c("NOTCH2_NP", "DLL1_NP", "HES1_NP", "CTNNB1_NP"), function(v) {
  # Select variable and STAGE
  df <- GRADE_KN_DF[, c(v, "Grade")]
  
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
variance_results_6 <- bind_rows(variance_results_6)
variance_results_6 #all good

##grade t test KN##################
GRADE_KN_DF$Grade
t.test(NOTCH2_NP ~ Grade,
       data = GRADE_KN_DF,
       var.equal = TRUE) # 0.08866
t.test(DLL1_NP ~ Grade,
       data = GRADE_KN_DF,
       var.equal = TRUE)
t.test(CTNNB1_NP ~ Grade,
       data = GRADE_KN_DF,
       var.equal = TRUE)
t.test(HES1_NP~ Grade,
       data = GRADE_KN_DF,
       var.equal = TRUE) #all not significant

##GRADE groups BOXPLOT#####################
each.vs.ref_sig2_grade <-  tibble::tribble(
  ~group1, ~group2, ~p.adj,   ~y.position, ~variable,
  "G1",   "G3", 0.08866, -2, "NOTCH2_NP", #anova
)
GRADE_table <- melt(GRADE_KN_DF[, c(7,15:18)], id.vars="Grade",
                    measure.vars=c("NOTCH2_NP", "DLL1_NP", "HES1_NP", "CTNNB1_NP"))

custom_colors_grade <- c("G3" = "deeppink","G1" = "lightblue") 
GRADE_plot <- ggplot(GRADE_table, aes(x=Grade , y=value, fill = variable)) +
  geom_boxplot( outlier.shape = NA , alpha=0.3, aes(fill = Grade )) +
  geom_jitter(aes(color = Grade ), size=1, alpha=0.5) +
  ylab(label = expression("Relative expression, normalized " * italic("GAPDH"))) + 
  facet_wrap(.~ variable, nrow = 2, scales = "free") +
  add_pvalue(each.vs.ref_sig2_grade, label = "p.adj") + #pvalue
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

GRADE_plot

#AGE correlations #####################
age_table <- TUMOR_df[, colnames(TUMOR_df) %in% c(c("NOTCH2_NP", "DLL1_NP", "HES1_NP", "CTNNB1_NP"), "Age")]
#normalcy
sapply(age_table, function(x) shapiro.test(x)$p.value) #cttnb1 not normal
#pearson correlation for all
results_n <- lapply(age_table[, colnames(age_table) %in% c("NOTCH2_NP", "DLL1_NP", "HES1_NP", "CTNNB1_NP")], 
                    function(x) cor.test(x, age_table$Age, method = "pearson"))
results_n
#spearman for ctnnb1
cor.test(age_table$Age,age_table$CTNNB1_NP, method = "spearman")
#COMPARE BETWEEN TISSUES #################################################
##NOTCH2#######################
notch2_long <- TUMOR_df %>%
  dplyr::select(`Laboratorinis kodas`, NOTCH2_NP, NOTCH2_TUMOR)%>%
  pivot_longer(
    cols = c(NOTCH2_NP, NOTCH2_TUMOR),
    names_to = "gene",
    values_to = "expression"
  )
# Run paired test
pval <- TUMOR_df %>%
  summarise(
    p_value = t.test(NOTCH2_TUMOR, NOTCH2_NP, paired = TRUE)$p.value
  ) %>%
  pull(p_value)
pval

# Plot with p-value
ggplot(notch2_long, aes(x = gene, y = expression, group = `Laboratorinis kodas`)) +
  geom_point(size = 2) +
  geom_line(alpha = 0.5) +
  theme_classic() +
  labs(
    x = "Tissue",
    y = "Gene expression",
    title = "Paired gene expression across tissues",
    subtitle = paste0("Paired t-test p = ", signif(pval, 3))
  )

##CTNNB1########################
CTNNB1_long <- TUMOR_df %>%
  dplyr::select(`Laboratorinis kodas`, CTNNB1_NP, CTNNB1_TUMOR)%>%
  pivot_longer(
    cols = c("CTNNB1_NP", "CTNNB1_TUMOR"),
    names_to = "gene",
    values_to = "expression"
  )
# Run paired test
pval2 <- TUMOR_df %>%
  summarise(
    p_value = t.test(CTNNB1_TUMOR, CTNNB1_NP, paired = TRUE)$p.value
  ) %>%
  pull(p_value)
pval2

# Plot with p-value
ggplot(CTNNB1_long, aes(x = gene, y = expression, group = `Laboratorinis kodas`)) +
  geom_point(size = 2) +
  geom_line(alpha = 0.5) +
  theme_classic() +
  labs(
    x = "Tissue",
    y = "Gene expression",
    title = "Paired gene expression across tissues",
    subtitle = paste0("Paired t-test p = ", signif(pval2, 3))
  )

##DLL1########################
DLL1_long <- TUMOR_df %>%
  dplyr::select(`Laboratorinis kodas`, DLL1_NP, DLL1_TUMOR)%>%
  pivot_longer(
    cols = c("DLL1_NP", "DLL1_TUMOR"),
    names_to = "gene",
    values_to = "expression"
  )
# Run paired test
pval3 <- TUMOR_df %>%
  summarise(
    p_value = t.test(DLL1_TUMOR, DLL1_NP, paired = TRUE)$p.value
  ) %>%
  pull(p_value)
pval3

# Plot with p-value
ggplot(DLL1_long, aes(x = gene, y = expression, group = `Laboratorinis kodas`)) +
  geom_point(size = 2) +
  geom_line(alpha = 0.5) +
  theme_classic() +
  labs(
    x = "Tissue",
    y = "Gene expression",
    title = "Paired gene expression across tissues",
    subtitle = paste0("Paired t-test p = ", signif(pval3, 3))
  )

##HES1########################
HES1_long <- TUMOR_df %>%
  dplyr::select(`Laboratorinis kodas`, HES1_NP, HES1_TUMOR)%>%
  pivot_longer(
    cols = c("HES1_NP", "HES1_TUMOR"),
    names_to = "gene",
    values_to = "expression"
  )
# Run paired test
pval4 <- TUMOR_df %>%
  summarise(
    p_value = t.test(HES1_TUMOR, HES1_NP, paired = TRUE)$p.value
  ) %>%
  pull(p_value)
pval4

# Plot with p-value
ggplot(HES1_long, aes(x = gene, y = expression, group = `Laboratorinis kodas`)) +
  geom_point(size = 2) +
  geom_line(alpha = 0.5) +
  theme_classic() +
  labs(
    x = "Tissue",
    y = "Gene expression",
    title = "Paired gene expression across tissues",
    subtitle = paste0("Paired t-test p = ", signif(pval4, 3))
  )
#SURVIVAL NP, KN, TISSUE GROUP KN###################################
#make df of only cases with survival
ALL_SURV_LIQUID <- TUMOR_df %>%
  filter(!is.na(OS), !is.na(STATUS)) #58 observations
table(ALL_SURV_LIQUID$STATUS, useNA = "a") #19 dead, 39 NOT
#make factor of gene expresssion
DATA_names_P <- c("NOTCH2_NP","CTNNB1_NP","DLL1_NP","HES1_NP",
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
table(HGSOC_SURV_LIQUID$STATUS, useNA = "a") #17 dead
#make only oc surv df
OC_SURV_LIQUID <- ALL_SURV_LIQUID %>%
  filter(TYPE %in%c("HGSOC", "OTHER")) 
table(OC_SURV_LIQUID$STATUS, useNA = "a") #19 dead

## Fit survival curves KN NP, TUMOR COHORT#################
ggsurvplot(survfit(Surv(OS, STATUS) ~ CTNNB1_NP_f, data = OC_SURV_LIQUID), 
           data = OC_SURV_LIQUID, pval = TRUE,
           title="Overall survival by CTNNB1 expression in NP, TUMOR cohort")

ggsurvplot(survfit(Surv(OS, STATUS) ~ NOTCH2_NP_f, data = OC_SURV_LIQUID), 
           data = OC_SURV_LIQUID, pval = TRUE,
           title="Overall survival by NOTCH2 expression in NP, TUMOR cohort")

ggsurvplot(survfit(Surv(OS, STATUS) ~ DLL1_NP_f, data = OC_SURV_LIQUID), 
           data = OC_SURV_LIQUID, pval = TRUE,
           title="Overall survival by DLL1 expression in NP, TUMOR cohort")

ggsurvplot(survfit(Surv(OS, STATUS) ~ HES1_NP_f, data = OC_SURV_LIQUID), 
           data = OC_SURV_LIQUID, pval = TRUE,
           title="Overall survival by HES1 expression in NP, TUMOR cohort")
##univariable cox, KN, NP all ##################
cox_model_notch2_np <- coxph(Surv(OS, STATUS) ~ NOTCH2_NP, data = OC_SURV_LIQUID)
summary(cox_model_notch2_np)

cox_model_CTNNB1_np <- coxph(Surv(OS, STATUS) ~ CTNNB1_NP, data = OC_SURV_LIQUID)
summary(cox_model_CTNNB1_np)

cox_model_DLL1_np <- coxph(Surv(OS, STATUS) ~ DLL1_NP, data = OC_SURV_LIQUID)
summary(cox_model_DLL1_np)# P = 0.074  .

cox_model_HES1_np <- coxph(Surv(OS, STATUS) ~ HES1_NP, data = OC_SURV_LIQUID)
summary(cox_model_HES1_np) #p =  0.0662 .


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
#SURVIVAL NP, HGSOC, TISSUE GROUP ###################################
table(HGSOC_SURV_LIQUID$STATUS, useNA = "a") #17 dead

## Fit survival curves hgsoc NP, TUMOR COHORT#################
ggsurvplot(survfit(Surv(OS, STATUS) ~ CTNNB1_NP_f, data = HGSOC_SURV_LIQUID), 
           data = HGSOC_SURV_LIQUID, pval = TRUE,
           title="Overall survival by CTNNB1 expression in NP, HGSOC TUMOR cohort")

ggsurvplot(survfit(Surv(OS, STATUS) ~ NOTCH2_NP_f, data = HGSOC_SURV_LIQUID), 
           data = HGSOC_SURV_LIQUID, pval = TRUE,
           title="Overall survival by NOTCH2 expression in NP,HGSOC TUMOR cohort")

ggsurvplot(survfit(Surv(OS, STATUS) ~ DLL1_NP_f, data = HGSOC_SURV_LIQUID), 
           data = HGSOC_SURV_LIQUID, pval = TRUE,
           title="Overall survival by DLL1 expression in NP, HGSOC TUMOR cohort")

ggsurvplot(survfit(Surv(OS, STATUS) ~ HES1_NP_f, data = HGSOC_SURV_LIQUID), 
           data = HGSOC_SURV_LIQUID, pval = TRUE,
           title="Overall survival by HES1 expression in NP, HGSOC TUMOR cohort")
##univariable cox, hgsoc, NP, tumor group ##################
cox_model_notch2_np2 <- coxph(Surv(OS, STATUS) ~ NOTCH2_NP, data = HGSOC_SURV_LIQUID)
summary(cox_model_notch2_np2)

cox_model_CTNNB1_np2 <- coxph(Surv(OS, STATUS) ~ CTNNB1_NP, data = HGSOC_SURV_LIQUID)
summary(cox_model_CTNNB1_np2)

cox_model_DLL1_np2 <- coxph(Surv(OS, STATUS) ~ DLL1_NP, data = HGSOC_SURV_LIQUID)
summary(cox_model_DLL1_np2)# P = 0.023  

cox_model_HES1_np2 <- coxph(Surv(OS, STATUS) ~ HES1_NP, data = HGSOC_SURV_LIQUID)
summary(cox_model_HES1_np2) #p =  0.301


#TUMOR DATA ####################################################################
#TISSUE NORMALCY TEST for groups #######################################
#chek NP in tumor tissue cohort
shapiro_results2 <- TUMOR_df[, c(36,31:34)] %>%
  pivot_longer(cols = -TYPE_tumor, names_to = "gene", values_to = "value") %>%
  group_by(TYPE_tumor, gene) %>%
  summarise(p_value = shapiro.test(value)$p.value, .groups = "drop") %>%
  filter(p_value < 0.05)
shapiro_results2 #normal
#chek NP variance in tumor tissue cohort
variance_results2 <- TUMOR_df[, c(36,31:34)] %>%
  pivot_longer(cols = -TYPE_tumor,
               names_to = "gene",
               values_to = "value") %>%
  group_by(gene) %>%
  summarise(
    p_value = leveneTest(value ~ TYPE_tumor)$`Pr(>F)`[1],
    .groups = "drop"
  ) 
variance_results2 #ctnnb1 not good variance

#ANOVA #######################################
#ALL normal
#NOTCH2
res_aovNOCTH2t <- aov( NOTCH2_TUMOR  ~ TYPE_tumor,
                       data = TUMOR_df
)
shapiro.test(res_aovNOCTH2t$residuals) #normality is met
summary(res_aovNOCTH2t) # p = 0.000133 ***
#TUCKEY posthoc
post_test_notch2T <- glht(res_aovNOCTH2t,
                          linfct = mcp(TYPE_tumor  = "Tukey"))

summary(post_test_notch2T) #hgsoc benign 0.000201, OTHER BENIGN <0.000283 
#HES1
res_aovHES1t <- aov( HES1_TUMOR ~ TYPE_tumor,
                     data = TUMOR_df
)
shapiro.test(res_aovHES1t$residuals) #normality is met
summary(res_aovHES1t) # p =5.99e-07 ***
#TUCKEY posthoc
post_test_HES1T <- glht(res_aovHES1t,
                        linfct = mcp(TYPE_tumor  = "Tukey"))

summary(post_test_HES1T) #hgsoc benign 0.000193, OTHER HGOSC < 1e-04
#DLL1
res_aovDLL1t <- aov( DLL1_TUMOR  ~ TYPE_tumor,
                     data = TUMOR_df
)
shapiro.test(res_aovDLL1t$residuals) #normality is met
summary(res_aovDLL1t) # p = 4.31e-05 ***
#TUCKEY posthoc
post_test_DLL1T <- glht(res_aovDLL1t,
                        linfct = mcp(TYPE_tumor  = "Tukey"))

summary(post_test_DLL1T) #hgsoc benign 0.000103 , OTHER HGSOC 0.022215  
#hes1
res_aovCTNNB1t <- oneway.test(CTNNB1_TUMOR ~ TYPE_tumor, 
                              data = TUMOR_df,
                              var.equal = FALSE)  # FALSE = Welch
res_aovCTNNB1t# p = 2.983e-06
#pairwise t test posthoc
pairwise_results <- pairwise.t.test(
  x = TUMOR_df$CTNNB1_TUMOR,
  g = TUMOR_df$TYPE_tumor,
  p.adjust.method = "BH",   # Benjamini-Hochberg correction
  pool.sd = FALSE            # unequal variance → Welch
)
pairwise_results #HGSOC - BENING 6.8e-10, hgsoc - OTHER 1e-0.0014   

#Tribble
each.vs.ref_sig2 <-  tibble::tribble(
  ~group1, ~group2, ~p.adj,   ~y.position, ~variable,
  "BENIGN",   "HGSOC", 0.001 , -1, "NOTCH2_TUMOR",
  "OTHER",   "HGSOC", 0.001 , -2, "NOTCH2_TUMOR",
  "BENIGN",   "HGSOC", 0.002 , -2, "HES1_TUMOR",
  "OTHER",   "HGSOC", 0.001, -1, "HES1_TUMOR",
  "BENIGN",   "HGSOC", 0.001 , -5, "DLL1_TUMOR",
  "OTHER",   "HGSOC", 0.022, -6, "DLL1_TUMOR",
  "BENIGN",   "HGSOC", 0.001 , 0, "CTNNB1_TUMOR",
  "OTHER",   "BENIGN", 0.002, -1, "CTNNB1_TUMOR",
)
#TISSUE boxplot 3 groups ########################################
#melt table for expression
Group3_tableT <- melt(TUMOR_df[, c(36,31:34)], id.vars="TYPE_tumor",
                      measure.vars=c("NOTCH2_TUMOR", "DLL1_TUMOR", "HES1_TUMOR", "CTNNB1_TUMOR"))

custom_colors <- c("HGSOC" = "deeppink","OTHER" = "lightpink", "BENIGN" = "lightblue") 
OC_plot2 <- ggplot(Group3_tableT, aes(x=TYPE_tumor , y=value, fill = variable)) +
  geom_boxplot( outlier.shape = NA , alpha=0.3, aes(fill = TYPE_tumor )) +
  geom_jitter(aes(color = TYPE_tumor ), size=1, alpha=0.5) +
  ylab(label = expression("Relative expression, normalized to " * italic("GAPDH"))) + 
  facet_wrap(.~ variable, nrow = 2, scales = "free") +
  add_pvalue(each.vs.ref_sig2, label = "p.adj") + #pvalue
  theme_minimal()+
  theme(
    strip.text.x = element_text(
      size = 12, face = "bold.italic"
    ),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5))+
  labs(x=NULL)+
  stat_boxplot(geom ='errorbar')+
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  scale_y_continuous(labels = function(x) 
    gsub("-", "\u2212", as.character(x))) #add long "-" signs

OC_plot2 + ggtitle("TISSUE")
#tissue survival #############################################
## Fit survival curves KN tissue, TUMOR COHORT#################
ggsurvplot(survfit(Surv(OS, STATUS) ~ CTNNB1_TUMOR_f, data = OC_SURV_LIQUID), 
           data = OC_SURV_LIQUID, pval = TRUE,
           title="Overall survival by CTNNB1 expression in TISSUES, TUMOR cohort")

ggsurvplot(survfit(Surv(OS, STATUS) ~ NOTCH2_TUMOR_f, data = OC_SURV_LIQUID), 
           data = OC_SURV_LIQUID, pval = TRUE,
           title="Overall survival by NOTCH2 expression in TISSUES, TUMOR cohort")

ggsurvplot(survfit(Surv(OS, STATUS) ~ DLL1_TUMOR_f, data = OC_SURV_LIQUID), 
           data = OC_SURV_LIQUID, pval = TRUE,
           title="Overall survival by DLL1 expression in TISSUES, TUMOR cohort")

ggsurvplot(survfit(Surv(OS, STATUS) ~ HES1_TUMOR_f, data = OC_SURV_LIQUID), 
           data = OC_SURV_LIQUID, pval = TRUE,
           title="Overall survival by HES1 expression in TISSUES, TUMOR cohort") #P= 0.045
##univariable cox, KN, TISSUE , TUMOR COHORT##################
cox_model_notch2_TISSUE <- coxph(Surv(OS, STATUS) ~ NOTCH2_TUMOR, data = OC_SURV_LIQUID)
summary(cox_model_notch2_TISSUE)

cox_model_CTNNB1_TISSUE <- coxph(Surv(OS, STATUS) ~ CTNNB1_TUMOR, data = OC_SURV_LIQUID)
summary(cox_model_CTNNB1_TISSUE)

cox_model_DLL1_TISSUE <- coxph(Surv(OS, STATUS) ~ DLL1_TUMOR, data = OC_SURV_LIQUID)
summary(cox_model_DLL1_TISSUE)

cox_model_HES1_TISSUE <- coxph(Surv(OS, STATUS) ~ HES1_TUMOR, data = OC_SURV_LIQUID)
summary(cox_model_HES1_TISSUE) 


## Fit survival curves HGSOC tissue, TUMOR COHORT#################
ggsurvplot(survfit(Surv(OS, STATUS) ~ CTNNB1_TUMOR_f, data = HGSOC_SURV_LIQUID), 
           data = HGSOC_SURV_LIQUID, pval = TRUE,
           title="Overall survival by CTNNB1 expression in TISSUES,HGSOC TUMOR cohort")

ggsurvplot(survfit(Surv(OS, STATUS) ~ NOTCH2_TUMOR_f, data = HGSOC_SURV_LIQUID), 
           data = HGSOC_SURV_LIQUID, pval = TRUE,
           title="Overall survival by NOTCH2 expression in TISSUES,HGSOC TUMOR cohort")

ggsurvplot(survfit(Surv(OS, STATUS) ~ DLL1_TUMOR_f, data = HGSOC_SURV_LIQUID), 
           data = HGSOC_SURV_LIQUID, pval = TRUE,
           title="Overall survival by DLL1 expression in TISSUES,HGSOC TUMOR cohort")

ggsurvplot(survfit(Surv(OS, STATUS) ~ HES1_TUMOR_f, data = HGSOC_SURV_LIQUID), 
           data = HGSOC_SURV_LIQUID, pval = TRUE,
           title="Overall survival by HES1 expression in TISSUES,HGSOC TUMOR cohort") #P= 0.045
##univariable cox, HGSOC, TISSUE , TUMOR COHORT##################
cox_model_notch2_TISSUE_hgsoc <- coxph(Surv(OS, STATUS) ~ NOTCH2_TUMOR, data = HGSOC_SURV_LIQUID)
summary(cox_model_notch2_TISSUE_hgsoc)

cox_model_CTNNB1_TISSUE_hgsoc <- coxph(Surv(OS, STATUS) ~ CTNNB1_TUMOR, data = HGSOC_SURV_LIQUID)
summary(cox_model_CTNNB1_TISSUE_hgsoc)

cox_model_DLL1_TISSUE_hgsoc <- coxph(Surv(OS, STATUS) ~ DLL1_TUMOR, data = HGSOC_SURV_LIQUID)
summary(cox_model_DLL1_TISSUE_hgsoc)

cox_model_HES1_TISSUE_hgsoc <- coxph(Surv(OS, STATUS) ~ HES1_TUMOR, data = HGSOC_SURV_LIQUID)
summary(cox_model_HES1_TISSUE_hgsoc) 

