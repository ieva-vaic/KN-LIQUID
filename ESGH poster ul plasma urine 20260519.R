#KN-  liquid 2026 05 19
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
library(effectsize)
#read RDS
LIQUID_DF_final <- readRDS("C:/Users/Ieva/rprojects/OTHER DATA/KN_LIQUID/liquid_20260415.RDS")
#np genes
DATA <- c("NOTCH2_NP","CTNNB1_NP","DLL1_NP","HES1_NP" )
#leave lavage only
#notch2 the fullest data
LAVAGE_df <- LIQUID_DF_final%>%
  filter(!is.na(NOTCH2_NP)) #103 cases

#RSS+BENING grouped NP #################################################
table(LAVAGE_df$TYPE, useNA = "a") #assess the situation
#create a new grouping where RSS and BENIGN is one
LAVAGE_df <- LAVAGE_df %>%
  mutate(
    TYPE_BENIGN2 = if_else(TYPE %in% c("RSS", "BENIGN"),
                           "BENIGN",
                           TYPE)
  )
table(LAVAGE_df$TYPE_BENIGN2, useNA = "a") #now 31 benign

##test normalcy of my variables, RSS+BENIGN#####################################
results_normalcy2 <- lapply(DATA, function(v) {
  by(LAVAGE_df[[v]], LAVAGE_df$TYPE_BENIGN2, function(x) {
    if (length(na.omit(x)) >= 3)
      shapiro.test(x)$p.value
    else NA
  })
})
names(results_normalcy2) <- DATA
results_normalcy2 #Not normal NP NOTCH2 endometrial cancer p = 0.01915093

##test variance RSS+BENIGN###############################
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

##ANOVA, RSS+BENIGN###########################
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
##boxplot full np, RSS+BENIGN ###########################
#melt table for expression
GroupNP_table_BENIGN2 <- melt(LAVAGE_df[, c(38,15:18)], id.vars="TYPE_BENIGN2",
                              measure.vars=c("NOTCH2_NP","CTNNB1_NP","DLL1_NP","HES1_NP"))
#fix names
GroupNP_table_BENIGN2 <- GroupNP_table_BENIGN2 %>%
  mutate(TYPE_BENIGN2 = dplyr::recode(TYPE_BENIGN2,
                                      "ENDOMETRIAL CANCER" = "EC",
                                      "OTHER" = "OTHER OC",
                                      #"BENIGN" = "BENIGN + RSS"
                                      ))
#make p values RSS+benign
each.vs.ref_sig_BENIGN2 <- tibble::tribble(
  ~group1, ~group2, ~p.adj,   ~y.position, ~variable,
  "EC",   "BENIGN", 0.001, -2, "CTNNB1_NP",#
  "EC",   "HGSOC", 0.001, -1, "CTNNB1_NP",#
  "EC",   "BENIGN", 0.002 , -2, "NOTCH2_NP",#dunn
  "EC",   "HGSOC", 0.001, -1, "NOTCH2_NP",#dunn
  #  "EC",   "OTHER OC", 0.04504, -1.5, "NOTCH2_NP",#wilcox, dunn not significant after adj
  "EC",   "HGSOC", 0.014, -2, "HES1_NP",#
  "BENIGN",   "HGSOC", 0.034, -1.5, "HES1_NP",#
  "EC",   "HGSOC", 0.047, -6, "DLL1_NP",#
  "EC",   "BENIGN", 0.003, -7, "DLL1_NP"#
)

#make figure RSS+benign
TYPE_RSS_BENIGN_plot <- ggplot(GroupNP_table_BENIGN2, aes(x=TYPE_BENIGN2 , y=value, fill = variable)) +
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
    "EC" = "#E64164",
    "OTHER OC" = "#5B8FF9",
    "BENIGN" = "#8FAF9C",
    "HGSOC" = "#002060"
  )) +
  scale_color_manual(values = c(
    "EC" = "#E64164",
    "OTHER OC" = "#5B8FF9",
    "BENIGN" = "#8FAF9C",
    "HGSOC" = "#002060"
  ))

TYPE_RSS_BENIGN_plot
#save RSS+benign
ggsave("C:/Users/Ieva/rprojects/outputs_all/LIQUID/NP_RSS_BEN_20260519esgh.png",
       plot = TYPE_RSS_BENIGN_plot,
       width = 18,
       height = 16,
       units = "cm",
       dpi = 400)
##ROC ENDOMETRIAL VS BENIGN / RSS########################################
#make ENDOMETRIAL VS BENIGN /RSS DF
ENDO_BENIGN_DF2<- LAVAGE_df %>%
  filter(TYPE_BENIGN2%in% c("BENIGN", "ENDOMETRIAL CANCER"))%>% #filter for right samples
  dplyr::select("TYPE_BENIGN2", "NOTCH2_NP","CTNNB1_NP","DLL1_NP","HES1_NP"
  )
ENDO_BENIGN_DF2$TYPE_BENIGN2 <- factor(ENDO_BENIGN_DF2$TYPE_BENIGN2)
ENDO_BENIGN_DF2$TYPE_BENIGN2 <- droplevels(ENDO_BENIGN_DF2$TYPE_BENIGN2) #drop unused
ENDO_BENIGN_DF2$TYPE_BENIGN2  <- relevel(ENDO_BENIGN_DF2$TYPE_BENIGN2 , ref = "BENIGN") #set control
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
  plot.roc(roc_results_np_endo_benign2[["NOTCH2_NP"]], print.auc = F, col = "#E64164",
           cex.main=0.8, 
           main ="Uterine lavage biomarkers EC vs Non-cancer samples",
           #xlab = "1 - Specifiškumas", 
           #ylab = "Jautrumas", 
           legacy.axes = T) #title
  lines(roc_results_np_endo_benign2[["CTNNB1_NP"]], col = "#002060", lwd =2) 
  lines(roc_results_np_endo_benign2[["DLL1_NP"]], col ="#8FAF9C", lwd =2) 
  lines(roc_results_np_endo_benign2[["HES1_NP"]], col = "#5B8FF9", lwd =2) 
  legend("bottomright", legend = c( expression(italic("NOTCH2")),
                                    expression(italic("CTNNB1")),
                                    expression(italic("DLL1")), 
                                    expression(italic("HES1"))
  ),
  
  col = c("#E64164", "#002060", "#8FAF9C", "#5B8FF9"), lty = 1, 
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
results_np_endo_benign2$Predictor <- c("NOCTH2", "CTNNB1", "DLL1", "HES1")
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
png("C:/Users/Ieva/rprojects/outputs_all/LIQUID/NP_EC_BENIGNRSS_ROC_20260519esgh.png",
    width = 11, height = 11, res = 300, units = "cm")
roc_plot_ENDO()
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
gtsave(gt_table_np_endo_ben2,vwidth = 10000,   
       filename = "C:/Users/Ieva/rprojects/outputs_all/LIQUID/NP_EC_BENIGNRSS_ROCtable_20260519esgh.png")
#import images
roc_image4rss  <- image_read("C:/Users/Ieva/rprojects/outputs_all/LIQUID/NP_EC_BENIGNRSS_ROC_20260519esgh.png")
table_image4rss <- image_read("C:/Users/Ieva/rprojects/outputs_all/LIQUID/NP_EC_BENIGNRSS_ROCtable_20260519esgh.png")

# resize table to match ROC image width
table_image4rss <- image_resize(table_image4rss,
                                paste0(image_info(roc_image4rss)$width, "x"))
# combine vertically
combined4rss <- image_append(c(roc_image4rss, table_image4rss), stack = TRUE)

# save
image_write(combined4rss,
            "C:/Users/Ieva/rprojects/outputs_all/LIQUID/NP_EC_BENIGNRSS_ROCcombined_20260519esgh.png")

##ROC ENDOMETRIAL VS HGSOC########################################
#make ENDOMETRIAL VS HGSOC /RSS DF
ENDO_HGSOC_DF<- LAVAGE_df %>%
  filter(TYPE_BENIGN2%in% c("HGSOC", "ENDOMETRIAL CANCER"))%>% #filter for right samples
  dplyr::select("TYPE_BENIGN2", "NOTCH2_NP","CTNNB1_NP","DLL1_NP","HES1_NP"
  )
ENDO_HGSOC_DF$TYPE_BENIGN2 <- factor(ENDO_HGSOC_DF$TYPE_BENIGN2)
ENDO_HGSOC_DF$TYPE_BENIGN2 <- droplevels(ENDO_HGSOC_DF$TYPE_BENIGN2) #drop unused
ENDO_HGSOC_DF$TYPE_BENIGN2  <- relevel(ENDO_HGSOC_DF$TYPE_BENIGN2 , ref = "ENDOMETRIAL CANCER") #set control

roc_results_np_endo_hgsoc<- lapply(c("NOTCH2_NP","CTNNB1_NP","DLL1_NP","HES1_NP"),
                                   function(col) {
                                     roc(response = ENDO_HGSOC_DF$TYPE_BENIGN2, predictor = ENDO_HGSOC_DF[[col]])})
names(roc_results_np_endo_hgsoc) <- c("NOTCH2_NP","CTNNB1_NP","DLL1_NP","HES1_NP")
roc_results_np_endo_hgsoc
#extract the aucs
auc_values_np_endo_hgsoc <- sapply(roc_results_np_endo_hgsoc, 
                                   function(roc_obj) {auc(roc_obj)})
auc_values_np_endo_hgsoc #extracted aucs
#roc figure 
roc_plot_ENDOHGSOC <- function() {
  par(pty = "s") #sets square
  plot.roc(roc_results_np_endo_hgsoc[["NOTCH2_NP"]], print.auc = F, col = "#E64164",
           cex.main=0.8, 
           main ="Uterine lavage biomarkers EC vs HGSOC",
           #xlab = "1 - Specifiškumas", 
           #ylab = "Jautrumas", 
           legacy.axes = T) #title
  lines(roc_results_np_endo_hgsoc[["CTNNB1_NP"]], col = "#002060", lwd =2) 
  lines(roc_results_np_endo_hgsoc[["DLL1_NP"]], col ="#8FAF9C", lwd =2) 
  lines(roc_results_np_endo_hgsoc[["HES1_NP"]], col = "#5B8FF9", lwd =2) 
  legend("bottomright", legend = c( expression(italic("NOTCH2")),
                                    expression(italic("CTNNB1")),
                                    expression(italic("DLL1")), 
                                    expression(italic("HES1"))
  ),
  
  col = c("#E64164", "#002060", "#8FAF9C", "#5B8FF9"), lty = 1, 
  cex = 0.7, lwd =3)
}
#plot
roc_plot_ENDOHGSOC()
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
#change names
results_np_endo_hgsoc$Predictor <- c("NOCTH2", "CTNNB1", "DLL1", "HES1")
#nice formating of the Table metrics for ROC OC
gt_table_np_endo_HGSOC <- results_np_endo_hgsoc %>%
  gt() %>%
  tab_header(
    title = "ROC measures for uterine lavage biomarkers", 
    subtitle = "HGSOC vs endometrial cancer") %>%
  fmt_number(
    columns = everything(),
    decimals = 3
  ) %>%
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_body(columns = vars(Predictor))
  )
#show
gt_table_np_endo_HGSOC
# Save the plot as a PNG file
png("C:/Users/Ieva/rprojects/outputs_all/LIQUID/NP_EC_hgsocroc_20260519esgh.png",
    width = 11, height = 11, res = 300, units = "cm")
roc_plot_ENDOHGSOC()
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
gtsave(gt_table_np_endo_HGSOC,vwidth = 10000,   
       filename = "C:/Users/Ieva/rprojects/outputs_all/LIQUID/NP_EC_hgsoc_ROCtable_20260519esgh.png")
#import images
roc_image  <- image_read("C:/Users/Ieva/rprojects/outputs_all/LIQUID/NP_EC_hgsocroc_20260519esgh.png")
table_image <- image_read("C:/Users/Ieva/rprojects/outputs_all/LIQUID/NP_EC_hgsoc_ROCtable_20260519esgh.png")

# resize table to match ROC image width
table_image <- image_resize(table_image,
                            paste0(image_info(roc_image)$width, "x"))
# combine vertically
combined <- image_append(c(roc_image, table_image), stack = TRUE)

# save
image_write(combined,
            "C:/Users/Ieva/rprojects/outputs_all/LIQUID/NP_EC_hgsoc_ROCcombined_20260519esgh.png")
##ROC HGSOC VS BENIGN RSS ##########################
#make hgsoc VS BENIGN and RSS DF
HGSOC_BENIGN_RSS_DF<- LAVAGE_df %>%
  filter(TYPE_BENIGN2%in% c("BENIGN", "HGSOC"))%>% #filter for right samples
  dplyr::select("TYPE_BENIGN2", "NOTCH2_NP","CTNNB1_NP","DLL1_NP","HES1_NP"
  )
HGSOC_BENIGN_RSS_DF$TYPE_BENIGN2 <- factor(HGSOC_BENIGN_RSS_DF$TYPE_BENIGN2)
HGSOC_BENIGN_RSS_DF$TYPE_BENIGN2 <- droplevels(HGSOC_BENIGN_RSS_DF$TYPE_BENIGN2) #drop unused
HGSOC_BENIGN_RSS_DF$TYPE_BENIGN2  <- relevel(HGSOC_BENIGN_RSS_DF$TYPE_BENIGN2 , ref = "BENIGN") #set control
roc_results_np_HGOSC_benignrss<- lapply(c("NOTCH2_NP","CTNNB1_NP","DLL1_NP","HES1_NP"), function(col) {
  roc(response = HGSOC_BENIGN_RSS_DF$TYPE_BENIGN2, predictor = HGSOC_BENIGN_RSS_DF[[col]])})
names(roc_results_np_HGOSC_benignrss) <- c("NOTCH2_NP","CTNNB1_NP","DLL1_NP","HES1_NP")
roc_results_np_HGOSC_benignrss
#extract the aucs
auc_values_np_HGSOC_benignrss <- sapply(roc_results_np_HGOSC_benignrss, function(roc_obj) {auc(roc_obj)})
auc_values_np_HGSOC_benignrss #extracted aucs
#roc figure 
roc_plot3rss <- function() {
  par(pty = "s") #sets square
  plot.roc(roc_results_np_HGOSC_benignrss[["NOTCH2_NP"]], print.auc = F, col = "#E64164",
           cex.main=0.8, 
           main ="Uterine lavage biomarkers HGSOC vs Non-cancer samples",
           #xlab = "1 - Specifiškumas", 
           #ylab = "Jautrumas", 
           legacy.axes = T) #title
  lines(roc_results_np_HGOSC_benignrss[["CTNNB1_NP"]], col = "#002060", lwd =2) 
  lines(roc_results_np_HGOSC_benignrss[["DLL1_NP"]], col ="#8FAF9C", lwd =2) 
  lines(roc_results_np_HGOSC_benignrss[["HES1_NP"]], col = "#5B8FF9", lwd =2) 
  legend("bottomright", legend = c( expression(italic("NOTCH2")),
                                    expression(italic("CTNNB1")),
                                    expression(italic("DLL1")), 
                                    expression(italic("HES1"))
  ),
  
  col = c("#E64164", "#002060", "#8FAF9C", "#5B8FF9"), lty = 1, 
  cex = 0.7, lwd =3)
}
#plot
roc_plot3rss()
#roc table
coords_results_np_HGSOC_benignrss<- lapply(roc_results_np_HGOSC_benignrss, function(roc_obj) {
  coords(roc_obj, "best", ret = c("threshold", "accuracy", "sensitivity",
                                  "specificity"),
         transpose = FALSE)
})
coords_results_np_HGSOC_benignrss
#create df
results_np_HGSOC_benignrss<- data.frame(
  Predictor = c("NOTCH2_NP","CTNNB1_NP","DLL1_NP","HES1_NP"),
  AUC = auc_values_np_HGSOC_benignrss,
  do.call(rbind,coords_results_np_HGSOC_benignrss) 
)
results_np_HGSOC_benignrss
#lithuanize it 
results_np_HGSOC_benignrss$Predictor <- c("NOCTH2", "CTNNB1", "DLL1", "HES1")
#nice formating of the Table metrics for ROC OC
gt_table_np_HGSOC_benrss <- results_np_HGSOC_benignrss %>%
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
gt_table_np_HGSOC_benrss
#save plot
# Save the plot as a PNG file
png("C:/Users/Ieva/rprojects/outputs_all/LIQUID/NP_HGSOC_BENIGNRSS_ROC20260519esgh.png",
    width = 11, height = 11, res = 300, units = "cm")
roc_plot3rss()
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
gtsave(gt_table_np_HGSOC_benrss,vwidth = 10000,   
       filename = "C:/Users/Ieva/rprojects/outputs_all/LIQUID/NP_HGSOC_BENIGNRSS_ROCtable20260519esgh.png")
#import images
roc_image3rss  <- image_read("C:/Users/Ieva/rprojects/outputs_all/LIQUID/NP_HGSOC_BENIGNRSS_ROC20260519esgh.png")
table_image3rss <- image_read("C:/Users/Ieva/rprojects/outputs_all/LIQUID/NP_HGSOC_BENIGNRSS_ROCtable20260519esgh.png")

# resize table to match ROC image width
table_image3rss <- image_resize(table_image3rss,
                                paste0(image_info(roc_image3rss)$width, "x"))
# combine vertically
combined3rss <- image_append(c(roc_image3rss, table_image3rss), stack = TRUE)

# save
image_write(combined3rss,
            "C:/Users/Ieva/rprojects/outputs_all/LIQUID/NP_HGSOC_BENIGNRSS_ROCcombined20260519esgh.png")


#15 CASE DF##################################
#leave URINE/PLASMA only
#remove cases that is NA in their type: 
LIQUID_DF_15 <- LIQUID_DF_final %>%
  filter(!is.na(HES1_URINE))

##SURVIVAL PLASMA##########################################
#make df of only cases with survival
ALL_SURV_LIQUID <- LIQUID_DF_15 %>%
  filter(!is.na(OS), !is.na(STATUS)) #58 observations

table(ALL_SURV_LIQUID$STATUS, useNA = "a") #2 dead, 12 NOT
#plasma cases
table(ALL_SURV_LIQUID$STATUS)
table(ALL_SURV_LIQUID$NOTCH2_P, ALL_SURV_LIQUID$STATUS)
table(ALL_SURV_LIQUID$DLL1_P, ALL_SURV_LIQUID$STATUS)
table(ALL_SURV_LIQUID$HES1_P, ALL_SURV_LIQUID$STATUS)
table(ALL_SURV_LIQUID$CTNNB1_P, ALL_SURV_LIQUID$STATUS)

##univariable cox, KN, plasma categorical data##################
cox_model_notch2_plasma <- coxph(Surv(OS, STATUS) ~ NOTCH2_P, data = ALL_SURV_LIQUID)
summary(cox_model_notch2_plasma)

cox_model_CTNNB1_plasma <- coxph(Surv(OS, STATUS) ~ CTNNB1_P, data = ALL_SURV_LIQUID)
summary(cox_model_CTNNB1_plasma)

#cox_model_DLL1_plasma <- coxph(Surv(OS, STATUS) ~ DLL1_P, data = ALL_SURV_LIQUID)
#summary(cox_model_DLL1_plasma)# not enough data - all cases no expression

cox_model_HES1_plasma <- coxph(Surv(OS, STATUS) ~ HES1_P, data = ALL_SURV_LIQUID)
summary(cox_model_HES1_plasma) #p=0.02 9.567e+09 (0-INF)

##univariable cox, KN, plasma numbered data##################
cox_model_HES1_plasma <- coxph(Surv(OS, STATUS) ~ HES1_P_norm, data = ALL_SURV_LIQUID)
summary(cox_model_HES1_plasma) #p = 0.4  1.735(0.477-6.31), n=4

##median survival plasma, KN categorical data###########################
hes1_fit <- survfit(Surv(OS, STATUS) ~ HES1_P, data = ALL_SURV_LIQUID)
summary(hes1_fit)$table #gives months
summary(hes1_fit, times = 12)$surv #gives survival prob at 1 year
summary(hes1_fit, times = 36)$surv #gives survival prob at 3 yrs
summary(hes1_fit, times = 60)$surv #gives survival prob at 5 yrs NULL
#what is median survival overall?
fit_all <- survfit(Surv(OS, STATUS) ~ 1, data = ALL_SURV_LIQUID)
summary(fit_all)$table #not reached

##Plots KM PLASMA##########################
## Fit survival curves PLASMA
#NOTCH2
ALL_SURV_LIQUID$NOTCH2 <- factor(
  ALL_SURV_LIQUID$NOTCH2_P,
  levels = c(0, 1),
  labels = c("no expression", "expression")
)
p_notch2 <- ggsurvplot(survfit(Surv(OS, STATUS) ~ NOTCH2, data = ALL_SURV_LIQUID), 
                       data = ALL_SURV_LIQUID, pval = F,
                       xlab = "Time (months)",
                       palette = c( "#002060", "#E64164"), 
                       title =expression( "Overall survival by " * italic("NOTCH2") * " expression in plasma"))
p_notch2$plot <- p_notch2$plot + labs(subtitle = "Log-rank  p = 0.7")
print(p_notch2)
#HES1
ALL_SURV_LIQUID$HES1 <- factor(
  ALL_SURV_LIQUID$HES1_P,
  levels = c(0, 1),
  labels = c("no expression", "expression")
)
p_hes1 <-ggsurvplot(survfit(Surv(OS, STATUS) ~ HES1, data = ALL_SURV_LIQUID), 
                    data = ALL_SURV_LIQUID, pval = F,
                    xlab = "Time (months)",
                    palette = c( "#002060", "#E64164"), 
                    title =expression( "Overall survival by " * italic("HES1") * " expression in plasma"))
p_hes1$plot <- p_hes1$plot + labs(subtitle = "Log-rank  p = 0.016")
print(p_hes1)
#DLL1 - not enough data
#CTNNB1
ALL_SURV_LIQUID$CTNNB1 <- factor(
  ALL_SURV_LIQUID$CTNNB1_P,
  levels = c(0, 1),
  labels = c("no expression", "expression")
)
p_ctnnb1 <- ggsurvplot(survfit(Surv(OS, STATUS) ~ CTNNB1, data = ALL_SURV_LIQUID), 
                       data = ALL_SURV_LIQUID, pval = F,
                       xlab = "Time (months)",
                       palette = c( "#002060", "#E64164"), 
                       title =expression( "Overall survival by " * italic("CTNNB1") * " expression in plasma"))

p_ctnnb1$plot <- p_ctnnb1$plot + labs(subtitle = "p = 0.45")
print(p_ctnnb1)

combined_plot <-
  p_hes1$plot /  p_notch2$plot / p_ctnnb1$plot 

ggsave(
  filename = "c:/Users/Ieva/rprojects/outputs_all/LIQUID/OC_PLASMA_survival_combined2026010ESGH.png",
  plot = combined_plot,
  width = 6,
  height = 10,
  dpi = 300
)
#BOXPLOT URINE#######################################
##URINE HGSOC vs others######################
#urine normalcy
by(LIQUID_DF_15$NOTCH2_URINE, LIQUID_DF_15$TYPE, shapiro.test)#0.5286
#by(LIQUID_DF_15$DLL1_URINE, LIQUID_DF_15$TYPE, shapiro.test) #Not enough
by(LIQUID_DF_15$HES1_URINE, LIQUID_DF_15$TYPE, shapiro.test)#0.2706
by(LIQUID_DF_15$CTNNB1_URINE, LIQUID_DF_15$TYPE, shapiro.test)#0.9347
#urine variance
car::leveneTest(LIQUID_DF_15$NOTCH2_URINE ~ LIQUID_DF_15$TYPE, center = median)#0.8207
car::leveneTest(LIQUID_DF_15$DLL1_URINE ~ LIQUID_DF_15$TYPE, center = median) #0.8207 
car::leveneTest(LIQUID_DF_15$HES1_URINE ~ LIQUID_DF_15$TYPE, center = median)#0.1547
car::leveneTest(LIQUID_DF_15$CTNNB1_URINE ~ LIQUID_DF_15$TYPE, center = median)#0.4438
#all normal and variance normal
#t.test URINE
t.test(LIQUID_DF_15$NOTCH2_URINE ~ LIQUID_DF_15$TYPE, var.equal = TRUE)#0.6513
#t.test(LIQUID_DF_15$DLL1_URINE  ~ LIQUID_DF_15$TYPE, var.equal = FALSE) #not enough
t.test(LIQUID_DF_15$HES1_URINE  ~ LIQUID_DF_15$TYPE, var.equal = TRUE) #0.6976
t.test(LIQUID_DF_15$CTNNB1_URINE  ~ LIQUID_DF_15$TYPE, var.equal = TRUE)#0.7676
#effect size
cohens_d(CTNNB1_URINE ~ TYPE,
         data = LIQUID_DF_15,
         pooled_sd = TRUE,
         hedges.correction = TRUE)
#-0.18     | [-1.36, 1.00]
cohens_d(NOTCH2_URINE ~ TYPE,
         data = LIQUID_DF_15,
         pooled_sd = TRUE,
         hedges.correction = TRUE)
#-0.32     | [-1.68, 1.05]
cohens_d(DLL1_URINE ~ TYPE,
         data = LIQUID_DF_15,
         pooled_sd = TRUE,
         hedges.correction = TRUE)
#-0.39     | [-2.76, 2.13]
cohens_d(HES1_URINE ~ TYPE,
         data = LIQUID_DF_15,
         pooled_sd = TRUE,
         hedges.correction = TRUE)
#-0.22     | [-1.29, 0.86]
##boxplots URINE ############################
#melt table for expression
URINE_table <- melt(LIQUID_DF_15[, c(12,19:22)], id.vars="TYPE",
                    measure.vars=c("NOTCH2_URINE","CTNNB1_URINE","DLL1_URINE","HES1_URINE"))

URINE_table$TYPE <- droplevels(URINE_table$TYPE)
URINE_table <- URINE_table %>%
  mutate(TYPE = dplyr::
           recode(TYPE, "OTHER" = "OTHER OC"))
URINE_table <- URINE_table %>%
  mutate(variable = dplyr::
           recode(variable, "NOTCH2_URINE" = "NOTCH2",
                  "CTNNB1_URINE" = "CTNNB1",
                  "DLL1_URINE" = "DLL1",
                  "HES1_URINE" = "HES1"  ))
each.vs.ref_sig2_g <- tibble::tribble(
  ~group1, ~group2, ~p.adj,   ~y.position, ~variable,
  "HGSOC",   "OTHER OC", 0.6513, -3, "NOTCH2",#
  #"HGSOC",   "OTHER OC", 0.001, -1, "DLL1",#
  "HGSOC",   "OTHER OC", 0.6976 , 1, "HES1",#dunn
  "HGSOC",   "OTHER OC", 0.7676, -2, "CTNNB1"#dunn
  
)
#TUMOR boxplot
custom_colors_grade <- c("HGSOC" = "#002060","OTHER OC" = "#5B8FF9") 
URINE_plot <- ggplot(URINE_table, aes(x=TYPE , y=value, fill = variable)) +
  geom_boxplot( outlier.shape = NA , alpha=0.3, aes(fill = TYPE )) +
  geom_jitter(aes(color = TYPE ), size=1, alpha=0.5) +
  ylab(label = expression("Gene expression, normalized to  " * italic("GAPDH"))) + 
  facet_wrap(.~ variable, nrow = 2, scales = "free",
             labeller = labeller(
               variable = c(
                 "CTNNB1" = "CTNNB1 expression in urine",
                 "DLL1" = "DLL1 expression in urine",
                 "HES1" = "HES1 expression in urine",
                 "NOTCH2" = "NOTCH2 expression in urine"
               )) ) +
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

URINE_plot
#save urine 15 cases boxplot
ggsave("C:/Users/Ieva/rprojects/outputs_all/LIQUID/15_urine_20260519esgh.png",
       plot = URINE_plot,
       width = 14,
       height = 16,
       units = "cm",
       dpi = 400)
#HEATMAP############################################
HEATDATA <- LAVAGE_df %>%
  dplyr::rename(
    NOTCH2_UL = NOTCH2_NP,
    CTNNB1_UL = CTNNB1_NP,
    DLL1_UL   = DLL1_NP,
    HES1_UL   = HES1_NP
  )
raiska_np <- c("NOTCH2_UL",      "CTNNB1_UL" ,     "DLL1_UL" ,       "HES1_UL")
raiska_urine <- c("NOTCH2_URINE", "CTNNB1_URINE", "DLL1_URINE", "HES1_URINE")
raiska_p <- c("NOTCH2_P" , "CTNNB1_P" ,"DLL1_P" ,"HES1_P")
## heatmap expression continuous data as matrix###########
Heat_data <- HEATDATA[, c("Laboratorinis kodas",
                          raiska_np,raiska_urine     
)]
colnames(Heat_data) <- gsub("_NP$", "_UL", colnames(Heat_data))
Heat_data <- as.data.frame(Heat_data)
rownames(Heat_data) <- Heat_data[, 1]
Heat_data <- Heat_data[, -1]
heatmap_raiska <- Heatmap(as.matrix(Heat_data), cluster_rows = F, cluster_columns = F)
heatmap_raiska #show
## clinical data as df####################################
clinical <- HEATDATA[, c("Laboratorinis kodas", "TYPE_BENIGN2",  "CA125" ,"Stage_simple",
                         "Grade", "Age", "STATUS")]
clinical <- as.data.frame(clinical)
rownames(clinical) <- clinical$KN
head(clinical)
#fix type
clinical$TYPE_BENIGN2 <- dplyr::recode(clinical$TYPE_BENIGN2,OTHER = "OTHER OC",
                                       `ENDOMETRIAL CANCER` = "EC")
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
col_age <- colorRamp2(
  c(40, 90),
  c("#F2A7B5", "#A9B8E8")
)

row_ha = rowAnnotation(
  `Tumor type` = clinical$TYPE_BENIGN2,
  Stage = clinical$Stage_simple,
  Grade = clinical$Grade,
  Age = clinical$Age,
  CA125 = clinical$CA125_f,
  `Survival status` = clinical$STATUS,
  
  col = list(
    
    `Tumor type` = c(
      "HGSOC"    = "#7C8DB5",  # muted dusty blue
      "BENIGN"   = "#A8CBB7",  # soft sage green
      "OTHER OC" = "#A7B9F2",  # pastel periwinkle
      "EC"       = "#F2A7B5"   # dusty rose
    ),
    
    Stage = c(
      "1"  = "#CFE3F6",  # very soft blue
      "2"  = "#E6D6F5",  # lavender pastel
      "3"  = "#F7C6D0",  # soft pink
      "4"  = "#D8C6F2",  # muted violet
      "NA" = "#F2F2F2"
    ),
    
    Grade = c(
      "G1"    = "#BFDCCF",  # sage pastel
      "G2&G1" = "#D7D2F2",  # soft lavender
      "G3"    = "#A7B9F2",  # periwinkle
      "NA"    = "#F2F2F2"
    ),
    
    CA125 = c(
      "Norm"           = "#BFDCCF",
      "CA125 increase" = "#F2A7B5",
      "NA"             = "#F2F2F2"
    ),
    
    `Survival status` = c(
      "Alive"    = "#BFDCCF",
      "Deceased" = "#F2A7B5",
      "NA"       = "#F2F2F2"
    ),
    
    Age = col_age
  )
)

# expression colors (soft gradient, no harsh contrast)
col_fun = colorRamp2(
  c(2, -5, -10, -15),
  c("#8FAFE8", "#CFE3F6", "#F2A7B5", "#6E86C7")
)
#col_fun = colorRamp2(c(2, 0, -2, -4, -6, -8, -10, -12, -14, -16),
#                     c("#e7e0fe", "#cec1fd", "#8564fb", "#64b3fb","#93cafc","#325a7e", "#e088bd", "#af2745", "#9e233e", "#4f121f"))
labels <- c("\u221215", "\u221210", "\u22125", "2" ) #THIS IS NEEDED FOR MDPI AT LEAST - LONG MINUS SIGNS

#
#make rows order
rows_order <-  c( "NOTCH2_URINE",  "CTNNB1_URINE",  "DLL1_URINE",    "HES1_URINE",
                  "NOTCH2_UL",     "CTNNB1_UL"  ,  "DLL1_UL" ,      "HES1_UL"  )
##create grouping ################################
# create groups from suffix
col_groups <- ifelse(grepl("_URINE$", rows_order), "URINE",
                     ifelse(grepl("_UL$", rows_order), "UL", NA))
# turn into factor with desired order
col_groups <- factor(col_groups, levels = c("URINE", "UL"))
# check
col_groups
##final heatmap#############################

heatmap_raiska <- Heatmap(as.matrix(Heat_data), 
                          cluster_rows = F,
                          cluster_columns = F,
                          name = "Relative gene expression",  
                          right_annotation = row_ha,
                          col = col_fun, 
                          row_split = clinical$TYPE, 
                          column_order = rows_order,
                          column_split = col_groups,
                          column_names_gp = gpar(fontface = "italic"),
                          column_title = "Relative gene expression", 
                          row_names_gp = gpar(fontsize = 8), 
                          heatmap_legend_param = list( #THIS IS FOR THE LONG MINUS SIGNS
                            at = c(2, -5, -10, -15),   # Legend positions
                            labels = labels)     # Adjusted labels
)
heatmap_raiska
# png("C:/Users/Ieva/rprojects/outputs_all/LIQUID/HEATMAP_tumor_cohort20260518.png",
#     width = 15, height = 25, res = 300, units = "cm")
# heatmap_raiska
# dev.off()


#add plasma data ###################################
plasma_only <- HEATDATA[, c("Laboratorinis kodas", raiska_p)]
plasma_only$NOTCH2_P <- factor(plasma_only$NOTCH2_P, 
                               levels = c(levels(plasma_only$NOTCH2_P), "NA"))
plasma_only$HES1_P <- factor(plasma_only$HES1_P, 
                             levels = c(levels(plasma_only$HES1_P), "NA"))
plasma_only$CTNNB1_P <- factor(plasma_only$CTNNB1_P, 
                               levels = c(levels(plasma_only$CTNNB1_P), "NA"))
plasma_only$DLL1_P <- factor(plasma_only$DLL1_P, 
                             levels = c(levels(plasma_only$DLL1_P), "NA"))
plasma_only$NOTCH2_P[is.na(plasma_only$NOTCH2_P)] <- "NA"
plasma_only$CTNNB1_P[is.na(plasma_only$CTNNB1_P)] <- "NA"
plasma_only$DLL1_P[is.na(plasma_only$DLL1_P)] <- "NA"
plasma_only$HES1_P[is.na(plasma_only$HES1_P)] <- "NA"

rownames(plasma_only) <- plasma_only$`Laboratorinis kodas`
plasma_only <- plasma_only[, -1]

plasma_only$CTNNB1_P <- dplyr::recode(
  plasma_only$CTNNB1_P,
  "0" = "no expression",
  "1" = "expression"
)
plasma_only$HES1_P <- dplyr::recode(
  plasma_only$HES1_P,
  "0" = "no expression",
  "1" = "expression"
)

plasma_only$DLL1_P <- dplyr::recode(
  plasma_only$DLL1_P,
  "0" = "no expression",
  "1" = "expression"
)
plasma_only$NOTCH2_P <- dplyr::recode(
  plasma_only$NOTCH2_P,
  "0" = "no expression",
  "1" = "expression"
)
col_fun <- c("expression" = "#F2A7B5", "no expression" = "#7C8DB5", "NA" ="grey")
heatmap_plasma <- Heatmap(as.matrix(plasma_only), 
                          col = col_fun, name = "Gene expression status in plasma",
                          cluster_rows = FALSE,
                          row_split = clinical$TYPE
)
heatmap_plasma

#add together with the main plot
main_plot <- heatmap_plasma + heatmap_raiska
main_plot

png("C:/Users/Ieva/rprojects/outputs_all/LIQUID/all_liquid_heatmap20260519esgh.png",
    width = 2500, height = 2800, res = 300)
main_plot
dev.off()
