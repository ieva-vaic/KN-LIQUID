#KN-  liquid 2026 05 18
#POSTER TISSUE/LAVAGE
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
library(ggpubr)
library(patchwork)
library(multcomp)
library(ggpubr) 
library(ComplexHeatmap)
library(circlize)
#read RDS
LIQUID_DF_final <- readRDS("C:/Users/Ieva/rprojects/OTHER DATA/KN_LIQUID/liquid_20260415.RDS")
#leave tumor only
TUMOR_df <- LIQUID_DF_final%>%
  filter(!is.na(NOTCH2_TUMOR)) #65 cases, 4 DOES NOT HAVE np
table(TUMOR_df$TYPE, useNA = "a")
# RSS case to BENIGN 
TUMOR_df$TYPE_tumor <- fct_recode(
  TUMOR_df$TYPE,
  BENIGN = "RSS"
)
TUMOR_df$TYPE_tumor <- droplevels(TUMOR_df$TYPE_tumor)
table(TUMOR_df$TYPE_tumor, useNA = "a")
#fix grade
table(TUMOR_df$Grade)
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
#NP, tissue cohort #############################
##NORMALCY TEST for groups #######################################
#chek NP in tumor tissue cohort
shapiro_results1 <- TUMOR_df[, c(38,15:18)] %>%
  pivot_longer(cols = -TYPE_tumor, names_to = "gene", values_to = "value") %>%
  group_by(TYPE_tumor, gene) %>%
  summarise(p_value = shapiro.test(value)$p.value, .groups = "drop") %>%
  filter(p_value < 0.05)
shapiro_results1 #not normal for benign dll1 and OC ctnnb1
#chek NP variance in tumor tissue cohort
variance_results1 <- TUMOR_df[, c(38, 15:18)] %>%
  pivot_longer(cols = -TYPE_tumor,
               names_to = "gene",
               values_to = "value") %>%
  group_by(gene) %>%
  summarise(
    p_value = leveneTest(value ~ TYPE_tumor)$`Pr(>F)`[1],
    .groups = "drop"
  ) 
variance_results1 #all normal

##ANOVA #######################################
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
##BOXPLOT NP, TUMOR COHORT#################################################
#Tribble
each.vs.ref_sig <-  tibble::tribble(
  ~group1, ~group2, ~p.adj,   ~y.position, ~variable,
  "BENIGN",   "HGSOC", 0.028, -2, "HES1", 
)
##boxplot 3 groups ########################################
#melt table for expression
Group3_table <- melt(TUMOR_df[, c(38,15:18)], id.vars="TYPE_tumor",
                     measure.vars=c("NOTCH2_NP", "DLL1_NP", "HES1_NP", "CTNNB1_NP"))
Group3_table$variable <- gsub("_NP", "", Group3_table$variable) #remove _np
custom_colors <- c("HGSOC" = "#8E24AA","OTHER" = "#BA68C8", "BENIGN" = "#00ACC1") 
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

OC_plot + ggtitle("Uterine lavage gene expression")

ggsave(
  filename = "c:/Users/Ieva/rprojects/outputs_all/LIQUID/NP_boxplot_tumorcohort20260518.png",
  plot = OC_plot + ggtitle("Uterine lavage gene expression"),
  width = 8,
  height =7,
  dpi = 500
)

##ROC NP tissue cohort###################################################################
#make groupings of diseases#################
OC_HGSOC_BENIGN<- TUMOR_df[c(TUMOR_df$TYPE_tumor != "OTHER"),] #51 cases left
OC_HGSOC_BENIGN$TYPE_tumor <- relevel(factor(OC_HGSOC_BENIGN$TYPE_tumor), ref = "BENIGN")

OC_HGSOC_OTHER<- TUMOR_df[c(TUMOR_df$TYPE_tumor != "BENIGN"),] 
OC_HGSOC_OTHER$TYPE_tumor <- relevel(factor(OC_HGSOC_OTHER$TYPE_tumor), ref = "OTHER")

OC_BENIGN_OTHER<- TUMOR_df[c(TUMOR_df$TYPE_tumor != "HGSOC"),] 
OC_BENIGN_OTHER$TYPE_tumor <- relevel(factor(OC_BENIGN_OTHER$TYPE_tumor), ref = "BENIGN")

##ROC np BH########################################################
roc_results_BH <- lapply(c("NOTCH2_NP", "DLL1_NP", "HES1_NP", "CTNNB1_NP"), function(col) {
  roc(response = OC_HGSOC_BENIGN$TYPE_tumor, predictor = OC_HGSOC_BENIGN[[col]])})
names(roc_results_BH) <- c("NOTCH2_NP", "DLL1_NP", "HES1_NP", "CTNNB1_NP")
roc_results_BH
#extract the aucs
auc_values_BH <- sapply(roc_results_BH, function(roc_obj) {auc(roc_obj)})
auc_values_BH #extracted aucs
#plot the plot
roc_plot <- function() {
  par(pty = "s") #sets square
  plot.roc(roc_results_BH[["NOTCH2_NP"]], print.auc = F, col = "#dcbeff",
           cex.main=0.8, 
           main ="Uterine lavage biomarkers, HGSOC vs benign",
          # xlab = "1 - Specifiškumas", 
          # ylab = "Jautrumas",
          legacy.axes = T) #title
  lines(roc_results_BH[["CTNNB1_NP"]], col = "#911eb4", lwd =2) 
  lines(roc_results_BH[["DLL1_NP"]], col ="#ffd8b1", lwd =2) 
  lines(roc_results_BH[["HES1_NP"]], col = "#42d4f4", lwd =2) 
  legend("bottomright", legend = c( expression(italic("NOTCH2")),
                                    expression(italic("CTNNB1")),
                                    expression(italic("DLL1")), 
                                    expression(italic("HES1"))
  ),
  
  col = c("#dcbeff", "#911eb4", "#ffd8b1", "#42d4f4"), lty = 1, 
  cex = 0.8, lwd =3)
}
#plot
roc_plot()

##ROC table BH NP################################
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
results_tumor$Predictor <- c("NOTCH2", "DLL1", "HES1", "CTNNB1")
results_tumor <- results_tumor %>%
  arrange(desc(AUC))
#nice formating of the Table metrics for ROC OC
gt_table_tumor <- results_tumor %>%
  gt() %>%
  tab_header(
    title = "ROC measures for uterine lavage biomarkers", 
    subtitle = "HGSOC vs benign") %>%
  fmt_number(
    columns = everything(),
    decimals = 3
  ) %>%
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_body(columns = vars(Predictor))
  )
#show
gt_table_tumor

##Save gt BH NP########################
# Save the plot as a PNG file
png("C:/Users/Ieva/rprojects/outputs_all/LIQUID/ROC_BH_NO_tumorcohort200260518.png",
    width = 10, height = 10, res = 500, units = "cm")
roc_plot()
# mtext(
#   "A",
#   side = 3,
#   adj = -0.2,
#   line = 1,
#   cex = 1.2,
#   font = 2
# )
dev.off()
gtsave(gt_table_tumor,vwidth = 10000,   
       filename = "C:/Users/Ieva/rprojects/outputs_all/LIQUID/ROCtable_BH_NO_tumorcohort200260518.png")
#import images
roc_imageoc  <- image_read("C:/Users/Ieva/rprojects/outputs_all/LIQUID/ROC_BH_NO_tumorcohort200260518.png")
table_imageoc<- image_read("C:/Users/Ieva/rprojects/outputs_all/LIQUID/ROCtable_BH_NO_tumorcohort200260518.png")

# resize table to match ROC image width
table_imageoc <- image_resize(table_imageoc,
                                       paste0(image_info(roc_imageoc)$width, "x"))
# combine vertically
combinedoc <- image_append(c(roc_imageoc, table_imageoc), stack = TRUE)

# save
image_write(combinedoc,
            "C:/Users/Ieva/rprojects/outputs_all/LIQUID/ROCfinal_BH_NO_tumorcohort200260518.png")

#TUMOR DATA ####################################################################
##TISSUE NORMALCY TEST for groups #######################################
#chek NP in tumor tissue cohort
shapiro_results2 <- TUMOR_df[, c(38,31:34)] %>%
  pivot_longer(cols = -TYPE_tumor, names_to = "gene", values_to = "value") %>%
  group_by(TYPE_tumor, gene) %>%
  summarise(p_value = shapiro.test(value)$p.value, .groups = "drop")
shapiro_results2 #normal
#chek NP variance in tumor tissue cohort
variance_results2 <- TUMOR_df[, c(38,31:34)] %>%
  pivot_longer(cols = -TYPE_tumor,
               names_to = "gene",
               values_to = "value") %>%
  group_by(gene) %>%
  summarise(
    p_value = leveneTest(value ~ TYPE_tumor)$`Pr(>F)`[1],
    .groups = "drop"
  ) 
variance_results2 #ctnnb1 not good variance

##ANOVA tissue#######################################
#NOTCH2 - normal
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

#CTNNB1
res_aovCTNNB1t <- oneway.test(CTNNB1_TUMOR ~ TYPE_tumor, 
                              data = TUMOR_df,
                              var.equal = FALSE)  # FALSE = Welch
res_aovCTNNB1t# p = 2.983e-06
#posthoc games-howel
TUMOR_df %>%
  games_howell_test(CTNNB1_TUMOR ~ TYPE_tumor)
#HGSOC - BENING 6.76e-10, BENIGN -OTHER 3   e- 3
sprintf("%.3f", 6.76e-10)
sprintf("%.3f", 3e-3 )
#Tribble
each.vs.ref_sig2 <-  tibble::tribble(
  ~group1, ~group2, ~p.adj,   ~y.position, ~variable,
  "BENIGN",   "HGSOC", 0.001 , -1, "NOTCH2",#<
  "OTHER",   "HGSOC", 0.001 , 0, "NOTCH2",#<
  "BENIGN",   "HGSOC", 0.002 , 0, "HES1",
  "OTHER",   "HGSOC", 0.001, -1, "HES1", #<
  "BENIGN",   "HGSOC", 0.001 , -0.5, "DLL1",# <
  "OTHER",   "HGSOC", 0.022, 0.5, "DLL1",
  "BENIGN",   "HGSOC", 0.001 , 0, "CTNNB1",#<, new posthoc
  "OTHER",   "BENIGN", 0.003, 1, "CTNNB1",#new posthoc
)
each.vs.ref_sig2$p.adj <- ifelse(
  each.vs.ref_sig2$p.adj <= 0.001,
  "<0.001",
  format(each.vs.ref_sig2$p.adj, scientific = FALSE)
)
##TISSUE boxplot 3 groups ########################################
#melt table for expression
Group3_tableT <- melt(TUMOR_df[, c(38,31:34)], id.vars="TYPE_tumor",
                      measure.vars=c("NOTCH2_TUMOR", "DLL1_TUMOR", "HES1_TUMOR", "CTNNB1_TUMOR"))
Group3_tableT$variable <- c("NOTCH2", "DLL1", "HES1", "CTNNB1")
custom_colors <- c("HGSOC" = "#C2185B","OTHER" = "#F48FB1", "BENIGN" = "#4FC3F7") 
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

OC_plot2 + ggtitle("Ovarian tissue gene expression")

ggsave(
  filename = "c:/Users/Ieva/rprojects/outputs_all/LIQUID/TISSUE_boxplot_tumorcohort20260518.png",
  plot = OC_plot2 + ggtitle("Ovarian tissue gene expression"),
  width = 8,
  height =7,
  dpi = 500
)
##ROC tissue BH########################################################
roc_results_BH_TISSUE <- lapply(c("NOTCH2_TUMOR", "DLL1_TUMOR", "HES1_TUMOR", "CTNNB1_TUMOR"), function(col) {
  roc(response = OC_HGSOC_BENIGN$TYPE_tumor, predictor = OC_HGSOC_BENIGN[[col]])})
names(roc_results_BH_TISSUE) <- c("NOTCH2_TUMOR", "DLL1_TUMOR", "HES1_TUMOR", "CTNNB1_TUMOR")
roc_results_BH_TISSUE
#extract the aucs
auc_values_BH_TISSUE <- sapply(roc_results_BH_TISSUE, function(roc_obj) {auc(roc_obj)})
auc_values_BH_TISSUE #extracted aucs
#plot the plot
roc_plot_TISSUE <- function() {
  par(pty = "s") #sets square
  plot.roc(roc_results_BH_TISSUE[["NOTCH2_TUMOR"]], print.auc = F, col = "#c77dff",
           cex.main=0.8, 
           main ="Ovarian tissue biomarkers, HGSOC vs benign",
           # xlab = "1 - Specifiškumas", 
           # ylab = "Jautrumas",
           legacy.axes = T) #title
  lines(roc_results_BH_TISSUE[["CTNNB1_TUMOR"]], col = "#6a0dad", lwd =2) 
  lines(roc_results_BH_TISSUE[["DLL1_TUMOR"]], col ="#ffbf80", lwd =2) 
  lines(roc_results_BH_TISSUE[["HES1_TUMOR"]], col = "#00c2ff", lwd =2) 
  legend("bottomright", legend = c( expression(italic("NOTCH2")),
                                    expression(italic("CTNNB1")),
                                    expression(italic("DLL1")), 
                                    expression(italic("HES1"))
  ),
  
  col = c("#c77dff", "#6a0dad", "#ffbf80", "#00c2ff"), lty = 1, 
  cex = 0.8, lwd =3)
}
#plot
roc_plot_TISSUE()

##ROC table BH tissue################################
#get roc features
coords_results_BH_TISSUE<- lapply(roc_results_BH_TISSUE, function(roc_obj) {
  coords(roc_obj, "best", ret = c("threshold", "accuracy", "sensitivity",
                                  "specificity"),
         transpose = FALSE)
})
coords_results_BH_TISSUE
#create df
results_tumor_TISSUE<- data.frame(
  Predictor = names(roc_results_BH_TISSUE),
  AUC = auc_values_BH_TISSUE,
  do.call(rbind,coords_results_BH_TISSUE) 
)
#lithuanize it 
results_tumor_TISSUE$Predictor <- c("NOTCH2", "DLL1", "HES1", "CTNNB1")
results_tumor_TISSUE <- results_tumor_TISSUE %>%
  arrange(desc(AUC))
#nice formating of the Table metrics for ROC OC
gt_table_tumor_TISSUE <- results_tumor_TISSUE %>%
  gt() %>%
  tab_header(
    title = "ROC measures for ovarian tissue biomarkers", 
    subtitle = "HGSOC vs benign") %>%
  fmt_number(
    columns = everything(),
    decimals = 3
  ) %>%
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_body(columns = vars(Predictor))
  )
#show
gt_table_tumor_TISSUE

##Save gt BH tissue########################
# Save the plot as a PNG file
png("C:/Users/Ieva/rprojects/outputs_all/LIQUID/ROC_BH_TISSUE_tumorcohort200260518.png",
    width = 10, height = 10, res = 500, units = "cm")
roc_plot_TISSUE()
# mtext(
#   "A",
#   side = 3,
#   adj = -0.2,
#   line = 1,
#   cex = 1.2,
#   font = 2
# )
dev.off()
gtsave(gt_table_tumor_TISSUE,vwidth = 10000,   
       filename = "C:/Users/Ieva/rprojects/outputs_all/LIQUID/ROCtable_BH__TISSUE_tumorcohort200260518.png")
#import images
roc_imageoc_TISSUE  <- image_read("C:/Users/Ieva/rprojects/outputs_all/LIQUID/ROC_BH_TISSUE_tumorcohort200260518.png")
table_imageoc_TISSUE<- image_read("C:/Users/Ieva/rprojects/outputs_all/LIQUID/ROCtable_BH__TISSUE_tumorcohort200260518.png")

# resize table to match ROC image width
table_imageoc_TISSUE <- image_resize(table_imageoc_TISSUE,
                                     paste0(image_info(roc_imageoc_TISSUE)$width, "x"))
# combine vertically
combinedoc_TISSUE <- image_append(c(roc_imageoc_TISSUE, table_imageoc_TISSUE), stack = TRUE)

# save
image_write(combinedoc_TISSUE,
            "C:/Users/Ieva/rprojects/outputs_all/LIQUID/ROCfinal_BH__TISSUE_tumorcohort200260518.png")
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
#change gene 
notch2_long$gene <- ifelse(notch2_long$gene == "NOTCH2_NP", "Uterine lavage",
                           ifelse(notch2_long$gene == "NOTCH2_TUMOR", "Tumor tissue",
                                  notch2_long$gene))
# Plot with p-value
NOTCH2pair <- ggplot(notch2_long, aes(x = gene, y = expression, group = `Laboratorinis kodas`)) +
  geom_point(size = 2) +
  geom_line(alpha = 0.5) +
  theme_classic() +
  labs(
    x = "Tissue type",
    y = expression(italic("NOTCH2")*" relative gene expression"),
    title =expression(italic("NOTCH2")*" gene expression across tissues"),
    subtitle = paste0("Paired t-test p < 0.001 ")
  )
NOTCH2pair
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
#change gene 
CTNNB1_long$gene <- ifelse(CTNNB1_long$gene == "CTNNB1_NP", "Uterine lavage",
                           ifelse(CTNNB1_long$gene == "CTNNB1_TUMOR", "Tumor tissue",
                                  CTNNB1_long$gene))
# Plot with p-value
CTNNB1pair <- ggplot(CTNNB1_long, aes(x = gene, y = expression, group = `Laboratorinis kodas`)) +
  geom_point(size = 2) +
  geom_line(alpha = 0.5) +
  theme_classic() +
  labs(
    x = "Tissue type",
    y = expression(italic("CTNNB1")*" relative gene expression"),
    title =expression(italic("CTNNB1")*" gene expression across tissues"),
    subtitle = paste0("Paired t-test p < 0.001 ")
  )
CTNNB1pair
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
#change gene 
DLL1_long$gene <- ifelse(DLL1_long$gene == "DLL1_NP", "Uterine lavage",
                         ifelse(DLL1_long$gene == "DLL1_TUMOR", "Tumor tissue",
                                DLL1_long$gene))
# Plot with p-value
DLL1paired <- ggplot(DLL1_long, aes(x = gene, y = expression, group = `Laboratorinis kodas`)) +
  geom_point(size = 2) +
  geom_line(alpha = 0.5) +
  theme_classic() +
  labs(
    x = "Tissue type",
    y = expression(italic("DLL1")*" relative gene expression"),
    title =expression(italic("DLL1")*" gene expression across tissues"),
    subtitle = paste0("Paired t-test p = ", signif(pval3, 3))
  )
DLL1paired

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
#change gene 
HES1_long$gene <- ifelse(HES1_long$gene == "HES1_NP", "Uterine lavage",
                         ifelse(HES1_long$gene == "HES1_TUMOR", "Tumor tissue",
                                HES1_long$gene))
#plot
HES1pair <- ggplot(HES1_long, aes(x = gene, y = expression, group = `Laboratorinis kodas`)) +
  geom_point(size = 2) +
  geom_line(alpha = 0.5) +
  theme_classic() +
  labs(
    x = "Tissue type",
    y = expression(italic("HES1")*" relative gene expression"),
    title =expression(italic("HES1")*" gene expression across tissues"),
    subtitle = paste0("Paired t-test p = ", signif(pval4, 3))
  )
HES1pair

##combine into one plot#################
combined <- 
  ( HES1pair |  NOTCH2pair )/
  (CTNNB1pair | DLL1paired)
ggsave(
  filename = "c:/Users/Ieva/rprojects/outputs_all/LIQUID/TISSUEcohort_PAIRED20260518.png",
  plot = combined ,
  width = 8,
  height =8,
  dpi = 300
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
#rename genes to make more sense in a plot
OC_SURV_LIQUID_NP <- OC_SURV_LIQUID
OC_SURV_LIQUID_NP$HES1 <- factor(
  OC_SURV_LIQUID_NP$HES1_NP_f,
  levels = c("Low", "High"),
  labels = c("Low expression", "High expression")
)
OC_SURV_LIQUID_NP$HES1 <- factor(
  OC_SURV_LIQUID_NP$HES1_NP_f,
  levels = c("Low", "High"),
  labels = c("Low expression", "High expression")
)
OC_SURV_LIQUID_NP$NOTCH2 <- factor(
  OC_SURV_LIQUID_NP$NOTCH2_NP_f,
  levels = c("Low", "High"),
  labels = c("Low expression", "High expression")
)
OC_SURV_LIQUID_NP$CTNNB1 <- factor(
  OC_SURV_LIQUID_NP$CTNNB1_NP_f,
  levels = c("Low", "High"),
  labels = c("Low expression", "High expression")
)
OC_SURV_LIQUID_NP$DLL1 <- factor(
  OC_SURV_LIQUID_NP$DLL1_NP_f,
  levels = c("Low", "High"),
  labels = c("Low expression", "High expression")
)
##univariable cox, KN, NP all ##################
cox_model_notch2_np <- coxph(Surv(OS, STATUS) ~ NOTCH2_NP, data = OC_SURV_LIQUID_NP)
summary(cox_model_notch2_np)

cox_model_CTNNB1_np <- coxph(Surv(OS, STATUS) ~ CTNNB1_NP, data = OC_SURV_LIQUID_NP)
summary(cox_model_CTNNB1_np)

cox_model_DLL1_np <- coxph(Surv(OS, STATUS) ~ DLL1_NP, data = OC_SURV_LIQUID_NP)
summary(cox_model_DLL1_np)# P = 0.074  .

cox_model_HES1_np <- coxph(Surv(OS, STATUS) ~ HES1_NP, data = OC_SURV_LIQUID_NP)
summary(cox_model_HES1_np) #p =  0.0662 .
##median survival NP, KN, TISSUE dataset###########################
hes1_fit <- survfit(Surv(OS, STATUS) ~ HES1, data = OC_SURV_LIQUID_NP)
summary(hes1_fit)$table #gives months
summary(hes1_fit, times = 12)$surv #gives survival prob at 1 year
summary(hes1_fit, times = 36)$surv #gives survival prob at 3 yrs
summary(hes1_fit, times = 60)$surv #gives survival prob at 5 yrs
#what is median survival overall?
fit_all <- survfit(Surv(OS, STATUS) ~ 1, data = OC_SURV_LIQUID_NP)
summary(fit_all)$table #not reached
##SURVIVAL PLOTS UL#################################
#HES1
hes1_fit <- survfit(Surv(OS, STATUS) ~ HES1, data = OC_SURV_LIQUID_NP)
KM_HES1 <- ggsurvplot(
  hes1_fit,
  data = OC_SURV_LIQUID_NP,
  pval = TRUE,
  risk.table = TRUE, 
  surv.median.line = "hv",   # adds horizontal & vertical median lines
  xlab = "Time (months)",
  title = expression( "Overall survival by"* italic(" HES1") *" expression in uterine lavage")
)
KM_HES1
#CTNNB1
CTNNB1_fit <- survfit(Surv(OS, STATUS) ~ CTNNB1, data = OC_SURV_LIQUID_NP)
KM_CTNNB1 <- ggsurvplot(
  CTNNB1_fit,
  data = OC_SURV_LIQUID_NP,
  pval = TRUE,
  risk.table = TRUE, 
  surv.median.line = "hv",   # adds horizontal & vertical median lines
  xlab = "Time (months)",
  title = expression( "Overall survival by"* italic(" CTNNB1") *" expression in uterine lavage")
)
KM_CTNNB1
#DLL1
DLL1_fit <- survfit(Surv(OS, STATUS) ~ DLL1, data = OC_SURV_LIQUID_NP)
KM_DLL1 <- ggsurvplot(
  DLL1_fit,
  data = OC_SURV_LIQUID_NP,
  pval = TRUE,
  risk.table = TRUE, 
  surv.median.line = "hv",   # adds horizontal & vertical median lines
  xlab = "Time (months)",
  title = expression( "Overall survival by"* italic(" DLL1") *" expression in uterine lavage")
)
KM_DLL1
#NOTCH2
NOTCH2_fit <- survfit(Surv(OS, STATUS) ~ NOTCH2, data = OC_SURV_LIQUID_NP)
KM_NOTCH2 <- ggsurvplot(
  NOTCH2_fit,
  data = OC_SURV_LIQUID_NP,
  pval = TRUE,
  risk.table = TRUE, 
  surv.median.line = "hv",   # adds horizontal & vertical median lines
  xlab = "Time (months)",
  title = expression( "Overall survival by"* italic(" NOTCH2") *" expression uterine lavage")
)
KM_NOTCH2
#combine plots
combinedKM <- 
  ( KM_HES1$plot |  KM_NOTCH2$plot )/
  (KM_CTNNB1$plot | KM_DLL1$plot)


combinedKM <- arrange_ggsurvplots(
  list(KM_HES1, KM_NOTCH2, KM_CTNNB1, KM_DLL1),
  ncol = 2,
  nrow = 2
)

ggsave(
  filename = "c:/Users/Ieva/rprojects/outputs_all/LIQUID/TISSUEcohortKM_20260518.png",
  plot = combinedKM ,
  width = 13,
  height =12,
  dpi = 300
)
#TISSUE KM##################################################
#rename genes to make more sense in a plot
OC_SURV_LIQUID_TISSUE <- OC_SURV_LIQUID
OC_SURV_LIQUID_TISSUE$HES1 <- factor(
  OC_SURV_LIQUID_TISSUE$HES1_TUMOR_f,
  levels = c("Low", "High"),
  labels = c("Low expression", "High expression")
)
OC_SURV_LIQUID_TISSUE$NOTCH2 <- factor(
  OC_SURV_LIQUID_TISSUE$NOTCH2_TUMOR_f,
  levels = c("Low", "High"),
  labels = c("Low expression", "High expression")
)
OC_SURV_LIQUID_TISSUE$CTNNB1 <- factor(
  OC_SURV_LIQUID_TISSUE$CTNNB1_TUMOR_f,
  levels = c("Low", "High"),
  labels = c("Low expression", "High expression")
)
OC_SURV_LIQUID_TISSUE$DLL1 <- factor(
  OC_SURV_LIQUID_TISSUE$DLL1_TUMOR_f,
  levels = c("Low", "High"),
  labels = c("Low expression", "High expression")
)
##univariable cox, KN, NP all ##################
cox_model_notch2_tissue <- coxph(Surv(OS, STATUS) ~ NOTCH2_TUMOR, data = OC_SURV_LIQUID_TISSUE)
summary(cox_model_notch2_tissue)

cox_model_CTNNB1_tissue<- coxph(Surv(OS, STATUS) ~ CTNNB1_TUMOR, data = OC_SURV_LIQUID_TISSUE)
summary(cox_model_CTNNB1_tissue)

cox_model_DLL1_tissue <- coxph(Surv(OS, STATUS) ~ DLL1_TUMOR, data = OC_SURV_LIQUID_TISSUE)
summary(cox_model_DLL1_tissue)# P = 0.074  .

cox_model_HES1__tissue <- coxph(Surv(OS, STATUS) ~ HES1_TUMOR, data = OC_SURV_LIQUID_TISSUE)
summary(cox_model_HES1__tissue) #p =  0.0662 .
##SURVIVAL PLOTS TISSUE#################################
#HES1
hes1_fit2 <- survfit(Surv(OS, STATUS) ~ HES1, data = OC_SURV_LIQUID_TISSUE)
KM_HES1tumor <- ggsurvplot(
  hes1_fit2,
  data = OC_SURV_LIQUID_TISSUE,
  pval = TRUE,
  risk.table = TRUE, 
  surv.median.line = "hv",   # adds horizontal & vertical median lines
  xlab = "Time (months)",
  title = expression( "Overall survival by"* italic(" HES1") *" expression in ovarian tissue")
)
KM_HES1tumor
#CTNNB1
CTNNB1_fit2 <- survfit(Surv(OS, STATUS) ~ CTNNB1, data = OC_SURV_LIQUID_TISSUE)
KM_CTNNB1tumor <- ggsurvplot(
  CTNNB1_fit2,
  data = OC_SURV_LIQUID_TISSUE,
  pval = TRUE,
  risk.table = TRUE, 
  surv.median.line = "hv",   # adds horizontal & vertical median lines
  xlab = "Time (months)",
  title = expression( "Overall survival by"* italic(" CTNNB1") *" expression in ovarian tissue")
)
KM_CTNNB1tumor
#DLL1
DLL1_fit2 <- survfit(Surv(OS, STATUS) ~ DLL1, data = OC_SURV_LIQUID_TISSUE)
KM_DLL1tumor <- ggsurvplot(
  DLL1_fit2,
  data = OC_SURV_LIQUID_TISSUE,
  pval = TRUE,
  risk.table = TRUE, 
  surv.median.line = "hv",   # adds horizontal & vertical median lines
  xlab = "Time (months)",
  title = expression( "Overall survival by"* italic(" DLL1") *" expression in ovarian tissue")
)
KM_DLL1tumor
#NOTCH2
NOTCH2_fit2 <- survfit(Surv(OS, STATUS) ~ NOTCH2, data = OC_SURV_LIQUID_TISSUE)
KM_NOTCH2tumor <- ggsurvplot(
  NOTCH2_fit2,
  data = OC_SURV_LIQUID_TISSUE,
  pval = TRUE,
  risk.table = TRUE, 
  surv.median.line = "hv",   # adds horizontal & vertical median lines
  xlab = "Time (months)",
  title = expression( "Overall survival by"* italic(" NOTCH2") *" expression ovarian tissue")
)
KM_NOTCH2tumor
#combine plots

combinedKM2 <- arrange_ggsurvplots(
  list(KM_HES1tumor, KM_NOTCH2tumor, KM_CTNNB1tumor, KM_DLL1tumor),
  ncol = 2,
  nrow = 2
)

ggsave(
  filename = "c:/Users/Ieva/rprojects/outputs_all/LIQUID/TISSUEcohortKMTUMORS_20260518.png",
  plot = combinedKM2 ,
  width = 13,
  height =12,
  dpi = 300
)

#HEATMAP############################################
raiska_np <- c("NOTCH2_UL",      "CTNNB1_UL" ,     "DLL1_UL" ,       "HES1_UL")
raiska_TISSUE <- c("NOTCH2_TUMOR", "CTNNB1_TUMOR", "DLL1_TUMOR", "HES1_TUMOR")
## heatmap expression continuous data as matrix###########
Heat_data <- TUMOR_df[, c("Laboratorinis kodas",
                          raiska_np,raiska_TISSUE     
)]
colnames(Heat_data) <- gsub("_NP$", "_UL", colnames(Heat_data))
Heat_data <- as.data.frame(Heat_data)
rownames(Heat_data) <- Heat_data[, 1]
Heat_data <- Heat_data[, -1]
heatmap_raiska <- Heatmap(as.matrix(Heat_data), cluster_rows = F)
heatmap_raiska #show
## clinical data as df####################################
clinical <- TUMOR_df[, c("Laboratorinis kodas", "TYPE_tumor",  "CA125" ,"Stage_simple",
                         "Grade", "Age", "STATUS")]
clinical <- as.data.frame(clinical)
rownames(clinical) <- clinical$KN
head(clinical)
#fix type
clinical$TYPE_tumor <- dplyr::recode(clinical$TYPE,OTHER = "OTHER OC")
#fix satus na so it shows on the legend
clinical$STATUS[is.na(clinical$STATUS)] <- "NA"
clinical$STATUS <- factor(clinical$STATUS)
clinical$STATUS <- dplyr::recode(clinical$STATUS, "1" = "Alive", "2" = "Deceased")#
#fix CA125
clinical$CA125_f <- ifelse(
  clinical$CA125 > 35,
  "CA125 increase",
  "Norm"
)
clinical$CA125_f[is.na(clinical$CA125_f)] <- "NA"
#clinical$CA125_f <- dplyr::recode(clinical$CA125_f, `CA125 increase`= "CA125 padidėjimas", Norm = "Norma")
#fix grade
clinical$Grade <- replace(clinical$Grade, is.na(clinical$Grade), "NA" )
clinical$Grade <- as.character(clinical$Grade)  # ensure it's character
clinical$Grade <- case_when(
  clinical$Grade %in% c("G2", "G1&G1", "G2&G1") ~ "G2&G1",
  clinical$Grade %in% c("GL", "GB") ~ NA_character_,
  TRUE ~ clinical$Grade  # keep G1, G3, NA as is
)
#fix stage
clinical$Stage_simple[is.na(clinical$Stage_simple)] <- "NA"
clinical$Stage_simple <- factor(clinical$Stage_simple)
# clinical data annotation
col_age <- colorRamp2(
  c(40, 90),
  c("#B8E1DD", "#6D597A")
)

row_ha = rowAnnotation(
  `Tumor type` = clinical$TYPE_tumor,
  Stage = clinical$Stage_simple,
  Grade = clinical$Grade,
  Age = clinical$Age,
  CA125 = clinical$CA125_f,
  `Survival status` = clinical$STATUS,
  
  col = list(
    
    `Tumor type` = c(
      "HGSOC" = "#9D6BCE",     # soft purple
      "BENIGN" = "#FFB3C6",    # pastel coral-pink
      "OTHER OC" = "#7DC7F2"   # pastel sky blue
    ),
    
    Stage = c(
      "1" = "#B8E1DD",   # pastel mint
      "2" = "#FFF0B3",   # soft pastel yellow
      "3" = "#FFCFD2",   # coral blush
      "4" = "#CDB4FF",   # lavender
      "NA" = "grey85"
    ),
    
    Grade = c(
      "G1" = "#A8E6CF",   # mint green
      "G3" = "#9D6BCE",   # matching purple
      "NA" = "grey85"
    ),
    
    CA125 = c(
      "Norm" = "#B8E1DD",          # mint
      "CA125 increase" = "#6D597A", # muted plum
      "NA" = "grey85"
    ),
    
    `Survival status` = c(
      "Alive" = "#B8E1DD",
      "Deceased" = "#6D597A",
      "NA" = "grey85"
    ),
    
    Age = col_age
  )
)
#expression colors
col_fun = colorRamp2(
  c(2, -5, -10, -15),
  c("#CDB4FF", "#90DBF4", "#FFCFD2", "#7B2CBF")
)
#col_fun = colorRamp2(c(2, 0, -2, -4, -6, -8, -10, -12, -14, -16),
#                     c("#e7e0fe", "#cec1fd", "#8564fb", "#64b3fb","#93cafc","#325a7e", "#e088bd", "#af2745", "#9e233e", "#4f121f"))
labels <- c("\u221215", "\u221210", "\u22125", "2" ) #THIS IS NEEDED FOR MDPI AT LEAST - LONG MINUS SIGNS
# create groups from suffix
col_groups <- ifelse(grepl("_UL$", rows_order), "UL", "TUMOR")
# turn into factor with desired order
col_groups <- factor(col_groups, levels = c( "UL", "TUMOR"))
# check
col_groups
#make rows order
rows_order <-  c(  "NOTCH2_UL",     "CTNNB1_UL"  ,  "DLL1_UL" ,      "HES1_UL"   ,
                   "NOTCH2_TUMOR", "CTNNB1_TUMOR", "DLL1_TUMOR", "HES1_TUMOR" )

##final heatmap#############################
heatmap_raiska <- Heatmap(as.matrix(Heat_data), 
                          cluster_rows = F,
                          cluster_columns = F,
                          name = "Relative gene expression",  
                          right_annotation = row_ha,
                          col = col_fun, 
                          row_split = clinical$TYPE, 
                          column_names_gp = gpar(fontface = "italic"),
                          column_title = "Relative gene expression", 
                          row_names_gp = gpar(fontsize = 8), 
                          column_order = rows_order,
                          column_split = col_groups,
                          heatmap_legend_param = list( #THIS IS FOR THE LONG MINUS SIGNS
                            at = c(2, -5, -10, -15),   # Legend positions
                            labels = labels)     # Adjusted labels
)
heatmap_raiska
png("C:/Users/Ieva/rprojects/outputs_all/LIQUID/HEATMAP_tumor_cohort20260518.png",
    width = 15, height = 25, res = 300, units = "cm")
heatmap_raiska
dev.off()
