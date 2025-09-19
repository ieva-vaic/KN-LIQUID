#KN-liquid - np
#data - 20250916, changed KN-97 from the repeated expreriment
Sys.setenv(LANG = "en")
#libraries
library(tidyverse)
library(pROC)
library(glmnet)
library(gtsummary)
library(gt)
library(grid)
library(brglm2)
library(reshape2)
library(rstatix) 
library(ggprism)
library(gridExtra)
library(png)
library(htmlwidgets)
library(webshot)
library(magick)
library(multcomp)
library(ggpubr)
#open liquid np data
np_liquid <- readxl::read_xlsx("C:/Users/Ieva/rprojects/OTHER DATA/KN_LIQUID/kn_liquid_np_20250916.xlsx")
#seprate names of the genes
raiska_np <- colnames(np_liquid[19:22])
#fix rownames
np_liquid$patient_id_aud <- gsub("-np", "", np_liquid$`KN nr.`)
#remove unnecessary columns
np_liquid <- np_liquid[, c("patient_id_aud", raiska_np)]
#data - all tumor and clinical
KN_data <- readRDS("C:/Users/Ieva/rprojects/KV-AUDINIAI/FINAL_sutikrinnimas_KN_AUDINIAI_MET_EXPR/KN_data1114_essential.rds")
#biomarker groups
#expression
raiska <- colnames(KN_data[18:27])
raiska_tissue <- c("NOTCH2", "CTNNB1", "DLL1" , "HES1" )
#methylation
metilinimas <- colnames(KN_data[28:31])
biomarkers <- c(raiska, metilinimas)
#change tumor
KN_data$tumor <- recode(KN_data$tumor, OvCa = "OC", Benign ="Benign")
#methylation data must be 0 1
KN_data$HOPX <- as.numeric(KN_data$HOPX) -1
KN_data$CDX2 <- as.numeric(KN_data$CDX2) -1
KN_data$ALX4 <- as.numeric(KN_data$ALX4) -1
KN_data$ARID1A_met <- as.numeric(KN_data$ARID1A_met) -1
#fix rownames 
rownames(KN_data) <- KN_data$patient_id_aud
#add np data to the main df########
#chek
np_liquid$patient_id_aud %in% KN_data$patient_id_aud 
#add
KN_data_full <- left_join(KN_data, np_liquid, by = "patient_id_aud")
#fix factors
KN_data_full$Stage4 <- factor(KN_data_full$Stage4)
#saveRDS(KN_data_full, "C:/Users/Ieva/rprojects/OTHER DATA/KN_data_np_tissue_full20250916.RDS")
#make groupings of diseases#################
OC_HGSOC_BENIGN<- KN_data_full[c(KN_data_full$Grupė_Ieva != "Other"),] #51 cases left
OC_HGSOC_BENIGN$tumor <- relevel(factor(OC_HGSOC_BENIGN$Grupė_Ieva), ref = "Benign")

OC_HGSOC_OTHER<- KN_data_full[c(KN_data_full$Grupė_Ieva != "Benign"),] 
OC_HGSOC_OTHER$tumor <- relevel(factor(OC_HGSOC_OTHER$Grupė_Ieva), ref = "Other")

OC_BENIGN_OTHER<- KN_data_full[c(KN_data_full$Grupė_Ieva != "HGSOC"),] 
OC_BENIGN_OTHER$tumor <- relevel(factor(OC_BENIGN_OTHER$Grupė_Ieva), ref = "Benign")

#make normalcy and f tests for np data#####################
#2 GROUPS
shapiro_results1 <- KN_data_full[, c(15,32:35)] %>%
  pivot_longer(cols = -tumor, names_to = "gene", values_to = "value") %>%
  group_by(tumor, gene) %>%
  summarise(p_value = shapiro.test(value)$p.value, .groups = "drop") %>%
  filter(p_value < 0.05)
shapiro_results1 #not normal for benign dll1 and OC ctnnb1

#3 GROUPS
shapiro_results2 <- KN_data_full[, c(3,32:35)] %>%
  pivot_longer(cols = -Grupė_Ieva, names_to = "gene", values_to = "value") %>%
  group_by(Grupė_Ieva, gene) %>%
  summarise(p_value = shapiro.test(value)$p.value, .groups = "drop") %>%
  filter(p_value < 0.05)
shapiro_results2 #not notmal for OC benign dll1 and other ctnnb1

leveneTest(variable ~ group)

#var test, benign OC
var_results_oc_b <- KN_data_full[, c(15,32:35)] %>%
  pivot_longer(cols = -tumor , names_to = "variable", values_to = "value") %>%
  group_by(variable) %>%
  summarise(
    p_value = var.test(value[tumor  == unique(tumor )[1]], 
                       value[tumor  == unique(tumor )[2]])$p.value,
    .groups = "drop"
  ) 
var_results_oc_b # equal for all 4

#var test, benign hgsoc
var_results_benign_hgsoc <- KN_data_full[, c(3,32:35)] %>%
  filter(Grupė_Ieva  != "Other")%>%
  pivot_longer(cols = -Grupė_Ieva , names_to = "variable", values_to = "value") %>%
  group_by(variable) %>%
  summarise(
    p_value = var.test(value[Grupė_Ieva  == unique(Grupė_Ieva )[1]], 
                       value[Grupė_Ieva  == unique(Grupė_Ieva )[2]])$p.value,
    .groups = "drop"
  ) 
#%>%filter(p_value < 0.05)
var_results_benign_hgsoc #all equal

#var test, benign hgsoc
var_results_benign_other<- KN_data_full[, c(3,32:35)] %>%
  filter(Grupė_Ieva  != "HGSOC")%>%
  pivot_longer(cols = -Grupė_Ieva , names_to = "variable", values_to = "value") %>%
  group_by(variable) %>%
  summarise(
    p_value = var.test(value[Grupė_Ieva  == unique(Grupė_Ieva )[1]], 
                       value[Grupė_Ieva  == unique(Grupė_Ieva )[2]])$p.value,
    .groups = "drop"
  ) 
#%>%filter(p_value < 0.05)
var_results_benign_other  #all equal

#var test, other hgsoc
var_results_other_other<- KN_data_full[, c(3,32:35)] %>%
  filter(Grupė_Ieva  != "Benign")%>%
  pivot_longer(cols = -Grupė_Ieva , names_to = "variable", values_to = "value") %>%
  group_by(variable) %>%
  summarise(
    p_value = var.test(value[Grupė_Ieva  == unique(Grupė_Ieva )[1]], 
                       value[Grupė_Ieva  == unique(Grupė_Ieva )[2]])$p.value,
    .groups = "drop"
  ) 
#%>%filter(p_value < 0.05)
var_results_other_other #all equal dist. 

#levene test instead, 3 groups
levene_rez1 <- KN_data_full[, c(3,32:35)] %>%
  pivot_longer(cols = -Grupė_Ieva , names_to = "variable", values_to = "value") %>%
  group_by(variable) %>%
  summarise(
    p_value = car::leveneTest(value ~ Grupė_Ieva )$`Pr(>F)`[1],
    .groups = "drop"
  )

levene_rez1

#NORMAL DISTRIBUTION AND VARIANCE: STUJENT T.TEST
#NOT NORMAL DISTRIBUTION AND VARIANCE:  Welch’s t-test (if near normal) OR Mann-Whitney U test
#welch test is safer as it adjusts for unequal sample sizes

#ANOVA####################################################
#DLL and ctnnb1 has to use Kruskal-Wallis test
#all equal vairances, othervise would use welch anova
#NOTCH2
res_aov <- aov( NOTCH2_DELTA ~ Grupė_Ieva ,
                data = KN_data_full
)
shapiro.test(res_aov$residuals) #normality is met
summary(res_aov)
#JAG2
res_aov2 <- aov( HES1_DELTA ~ Grupė_Ieva ,
                 data = KN_data_full
)
summary(res_aov2) # p 0.37
shapiro.test(res_aov2$residuals) #normality is met
# Tukey HSD test:
post_test <- glht(res_aov2,
                  linfct = mcp(Grupė_Ieva  = "Tukey")
)

summary(post_test) #hgsoc benign 0,0277

#not normal DLL1
kruskal1 <- kruskal.test(DLL1_DELTA ~ Grupė_Ieva,
                         data = KN_data_full
)
kruskal1

#not normal CTNNB1
kruskal2 <- kruskal.test(CTNNB1_DELTA ~ Grupė_Ieva,
                         data = KN_data_full
)
kruskal2 

#Tribble all tests together 3 groups#####################
each.vs.ref_sig <- 
  each.vs.ref_sig2 <- tibble::tribble(
    ~group1, ~group2, ~p.adj,   ~y.position, ~variable,
    "Gerybiniai",   "HGSOC", 0.028, -2, "HES1", #stjudent's
  )

#boxplot 3 groups ########################################

#MELT for comarisons
#melt table for expression
Group3_table <- melt(KN_data_full[, c(3,32:35)], id.vars="Grupė_Ieva",  measure.vars=raiska_np)

#rename to lt 
Group3_table$Grupė_Ieva
Group3_table2 <- Group3_table %>%
  mutate(Grupė_Ieva = case_when(
    Grupė_Ieva == "HGSOC" ~ "HGSOC",
    Grupė_Ieva == "Benign" ~ "Gerybiniai",
    Grupė_Ieva == "Other" ~ "Kiti KV"
  )) %>%
  mutate(variable  = case_when(
    variable       == "NOTCH2_DELTA" ~ "NOTCH2",
    variable       == "CTNNB1_DELTA" ~ "CTNNB1",
    variable       == "DLL1_DELTA" ~ "DLL1",
    variable       == "HES1_DELTA" ~ "HES1"
  ))


custom_colors <- c("HGSOC" = "deeppink","Kiti KV" = "lightpink", "Gerybiniai" = "blue") 
OC_plot <- ggplot(Group3_table2, aes(x=Grupė_Ieva , y=value, fill = variable)) +
  geom_boxplot( outlier.shape = NA , alpha=0.3, aes(fill = Grupė_Ieva )) +
  geom_jitter(aes(color = Grupė_Ieva ), size=1, alpha=0.5) +
  ylab(label = expression("Santykinė genų raiška, normalizuota pagal  " * italic("GAPDH"))) + 
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

OC_plot

# Save the plot as a PNG file
png("C:/Users/Ieva/rprojects/outputs_all/np_boxplot20250918.png",
    width = 1000, height = 1100, res = 200)
OC_plot
dev.off()

#FC 3 groups###########################################
expression_df <- KN_data_full[, c("Grupė_Ieva", raiska_np, "patient_id_aud")]
rownames(expression_df) <- expression_df$patient_id_aud 
expression_df <- expression_df[, -6]
#REMOVE NA
expression_df2 <- expression_df %>%
  filter(complete.cases(.))

exp_df2 <- expression_df2 %>%
  melt(id.vars="Grupė_Ieva",  measure.vars=raiska_np) %>%
  rename(gene = variable) %>%
  rename(expression = value)%>%
  rename(condition = Grupė_Ieva) %>%
  group_by(condition, gene)%>%
  summarise(mean_expression = mean(expression, drop_na= T)) %>%
  spread(condition, mean_expression)

mean_expression_OC2 <- exp_df2 %>%
  mutate(fold_change_HB = log2(2^(`HGSOC` - `Benign`))) %>% #benign vs hgsoc 
  mutate(fold_change_HO = log2(2^(`HGSOC` - `Other`))) %>% #others vs hgsoc 
  mutate(fold_change_OB = log2(2^(`Other` / `Benign`))) #benign vs others
mean_expression_OC2 


#ROC BH########################################################
roc_results_np<- lapply(raiska_np, function(col) {
  roc(response = OC_HGSOC_BENIGN$tumor, predictor = OC_HGSOC_BENIGN[[col]])})
names(roc_results_np) <- raiska_np
roc_results_np
#extract the aucs
auc_values_tumor <- sapply(roc_results_np, function(roc_obj) {auc(roc_obj)})
auc_values_tumor #extracted aucs

#roc figure BH#############################################
roc_plot <- function() {
  par(pty = "s") #sets square
  plot.roc(roc_results_np[["NOTCH2_DELTA"]], print.auc = F, col = "#dcbeff",
           cex.main=0.8, 
           main ="Gerybinių pokyčių atskyrimas nuo HGSOC",
           xlab = "Specifiškumas", 
           ylab = "Jautrumas") #title
  lines(roc_results_np[["CTNNB1_DELTA"]], col = "#911eb4", lwd =2) 
  lines(roc_results_np[["DLL1_DELTA"]], col ="#ffd8b1", lwd =2) 
  lines(roc_results_np[["HES1_DELTA"]], col = "#42d4f4", lwd =2) 
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

# Save the plot as a PNG file
png("C:/Users/Ieva/rprojects/outputs_all/bh_roc_np20250916.png",
    width = 1000, height = 1000, res = 200)
roc_plot()
dev.off()

#roc table BH################################
#get roc features
coords_results_np<- lapply(roc_results_np, function(roc_obj) {
  coords(roc_obj, "best", ret = c("threshold", "accuracy", "sensitivity",
                                  "specificity", "precision", "npv", "tpr", "fpr"),
         transpose = FALSE)
})
coords_results_np
#create df
results_tumor<- data.frame(
  Predictor = raiska_np,
  AUC = auc_values_tumor,
  do.call(rbind,coords_results_np) 
)
#lithuanize it 
colnames(results_tumor) <- c("Biožymuo", "plotas po kreive", "slenkstinė vertė", 
                             "tikslumas", "jautrumas", "specifiškumas", 
                             "ppv", "npv", "tpr", "fpr")
rownames(results_tumor) <- c("NOCTH2", "CTNNB1", "DLL1", "HES1")
results_tumor$Biožymuo <- c("NOCTH2", "CTNNB1", "DLL1", "HES1")
#nice formating of the Table metrics for ROC OC
gt_table_tumor <- results_tumor %>%
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
gt_table_tumor

#save table
gtsave(gt_table_tumor,
       filename = "C:/Users/Ieva/rprojects/outputs_all/bh_roc_table_np_20250916.png")

#Combine the images
roc_image1<- image_read("C:/Users/Ieva/rprojects/outputs_all/bh_roc_np20250916.png")
table_image1 <- image_read("C:/Users/Ieva/rprojects/outputs_all/bh_roc_table_np_20250916.png")

combined_image1 <- image_append(c(roc_image1, table_image1), stack = F)

# Find the max width to align both
roc_info <- image_info(roc_image1)
table_info <- image_info(table_image1)
max_width <- max(roc_info$width, table_info$width)

# Pad each image to the max width
roc_image1_padded <- image_extent(roc_image1, geometry = geometry_area(max_width, roc_info$height), gravity = "center", color = "white")
table_image1_padded <- image_extent(table_image1, geometry = geometry_area(max_width, table_info$height), gravity = "center", color = "white")


combined_image1 <- image_append(c(roc_image1_padded, table_image1_padded), stack = T)

# Save the combined image
image_write(combined_image1, 
            "C:/Users/Ieva/rprojects/outputs_all/bh_roc_full_np_20250916.png")

#ROC OH########################################################
roc_results_np_oh<- lapply(raiska_np, function(col) {
  roc(response = OC_HGSOC_OTHER$tumor, predictor = OC_HGSOC_OTHER[[col]])})
names(roc_results_np_oh) <- raiska_np
roc_results_np_oh
#extract the aucs
auc_values_tumor_oh <- sapply(roc_results_np_oh, function(roc_obj) {auc(roc_obj)})
auc_values_tumor_oh #extracted aucs


#roc figure OH#############################################
roc_plot2 <- function() {
  par(pty = "s") #sets square
  plot.roc(roc_results_np_oh[["NOTCH2_DELTA"]], print.auc = F, col = "#dcbeff",
           cex.main=0.8, 
           main ="HGSOC atskyrimas nuo kitų KV",
           xlab = "Specifiškumas", 
           ylab = "Jautrumas") #title
  lines(roc_results_np_oh[["CTNNB1_DELTA"]], col = "#911eb4", lwd =2) 
  lines(roc_results_np_oh[["DLL1_DELTA"]], col ="#ffd8b1", lwd =2) 
  lines(roc_results_np_oh[["HES1_DELTA"]], col = "#42d4f4", lwd =2) 
  #SUDĖTA NE PAGAL PAVEIKLĄ BET PAGAL AUC DYDI
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

# Save the plot as a PNG file
png("C:/Users/Ieva/rprojects/outputs_all/oh_roc_np20250916.png",
    width = 1000, height = 1000, res = 200)
roc_plot2()
dev.off()

#roc table OH################################
#get roc features
coords_results_np_oh<- lapply(roc_results_np_oh, function(roc_obj) {
  coords(roc_obj, "best", method = "youden", ret = c("threshold", "accuracy", "sensitivity",
                                  "specificity", "precision", "npv", "tpr", "fpr"),
         transpose = FALSE)
})
coords_results_np_oh
#there is CTNNB1 two thresholds, choose first
# Keep only the first row of CTNNB1_DELTA
coords_results_np_oh$CTNNB1_DELTA <- coords_results_np_oh$CTNNB1_DELTA[1, , drop = FALSE]
#create df
results_tumor2<- data.frame(
  Predictor = raiska_np,
  AUC = auc_values_tumor_oh,
  do.call(rbind,coords_results_np_oh) 
)
#lithuanize it 
colnames(results_tumor2) <- c("Biožymuo", "plotas po kreive", "slenkstinė vertė", 
                              "tikslumas", "jautrumas", "specifiškumas", 
                              "ppv", "npv", "tpr", "fpr")
rownames(results_tumor2) <- c("NOCTH2", "CTNNB1", "DLL1", "HES1")
results_tumor2$Biožymuo <- c("NOCTH2", "CTNNB1", "DLL1", "HES1")
#nice formating of the Table metrics for ROC OC
gt_table_tumor2 <- results_tumor2 %>%
  gt() %>%
  tab_header(
    title = "ROC kriterijai", 
    subtitle = "HGSOC atskyrimas nuo kitų KV") %>%
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

#save table
gtsave(gt_table_tumor2,
       filename = "C:/Users/Ieva/rprojects/outputs_all/oh_roc_table_np_20250916.png")

#Combine the images
roc_image2<- image_read("C:/Users/Ieva/rprojects/outputs_all/oh_roc_np20250916.png")
table_image2 <- image_read("C:/Users/Ieva/rprojects/outputs_all/oh_roc_table_np_20250916.png")

combined_image2 <- image_append(c(roc_image2, table_image2), stack = F)

# Find the max width to align both
roc_info2 <- image_info(roc_image2)
table_info2 <- image_info(table_image2)
max_width2 <- max(roc_info2$width, table_info2$width)

# Pad each image to the max width
roc_image1_padded2 <- image_extent(roc_image2, geometry = geometry_area(max_width2, roc_info2$height), gravity = "center", color = "white")
table_image1_padded2 <- image_extent(table_image2, geometry = geometry_area(max_width2, table_info2$height), gravity = "center", color = "white")

combined_image2 <- image_append(c(roc_image1_padded2, table_image1_padded2), stack = T)

# Save the combined image
image_write(combined_image2, 
            "C:/Users/Ieva/rprojects/outputs_all/oh_roc_full_np_20250916.png")

#ROC BO########################################################
roc_results_np_bo<- lapply(raiska_np, function(col) {
  roc(response = OC_BENIGN_OTHER$tumor, predictor = OC_BENIGN_OTHER[[col]])})
names(roc_results_np_bo) <- raiska_np
roc_results_np_bo
#extract the aucs
auc_values_tumor_bo <- sapply(roc_results_np_bo, function(roc_obj) {auc(roc_obj)})
auc_values_tumor_bo #extracted aucs


#roc figure BO#############################################
roc_plot3 <- function() {
  par(pty = "s") #sets square
  plot.roc(roc_results_np_bo[["NOTCH2_DELTA"]], print.auc = F, col = "#dcbeff",
           cex.main=0.8, 
           main ="Gerybinių pokyčių atskyrimas nuo I tipo KV",
           xlab = "Specifiškumas", 
           ylab = "Jautrumas") #title
  lines(roc_results_np_bo[["CTNNB1_DELTA"]], col = "#911eb4", lwd =2) 
  lines(roc_results_np_bo[["DLL1_DELTA"]], col ="#ffd8b1", lwd =2) 
  lines(roc_results_np_bo[["HES1_DELTA"]], col = "#42d4f4", lwd =2) 
  #SUDĖTA NE PAGAL PAVEIKLĄ BET PAGAL AUC DYDI
  legend("bottomright", legend = c( expression(italic("NOTCH2")),
                                    expression(italic("CTNNB1")),
                                    expression(italic("DLL1")), 
                                    expression(italic("HES1"))
  ),
  
  col = c("#dcbeff", "#911eb4", "#ffd8b1", "#42d4f4"), lty = 1, 
  cex = 0.5, lwd =3)
}
#plot
roc_plot3()

# Save the plot as a PNG file
png("C:/Users/Ieva/rprojects/outputs_all/bo_roc_np20250916.png",
    width = 1000, height = 1000, res = 200)
roc_plot3()
dev.off()

#roc table BO################################
#get roc features
coords_results_np_bo<- lapply(roc_results_np_bo, function(roc_obj) {
  coords(roc_obj, "best", ret = c("threshold", "accuracy", "sensitivity",
                                  "specificity", "precision", "npv", "tpr", "fpr"),
         transpose = FALSE)
})
coords_results_np_bo
#create df
results_tumor3<- data.frame(
  Predictor = raiska_np,
  AUC = auc_values_tumor_bo,
  do.call(rbind,coords_results_np_bo) 
)
#lithuanize it 
colnames(results_tumor3) <- c("Biožymuo", "plotas po kreive", "slenkstinė vertė", 
                              "tikslumas", "jautrumas", "specifiškumas", 
                              "ppv", "npv", "tpr", "fpr")
rownames(results_tumor3) <- c("NOCTH2", "CTNNB1", "DLL1", "HES1")
results_tumor3$Biožymuo <- c("NOCTH2", "CTNNB1", "DLL1", "HES1")
#nice formating of the Table metrics for ROC OC
gt_table_tumor3 <- results_tumor3 %>%
  gt() %>%
  tab_header(
    title = "ROC kriterijai", 
    subtitle = "Gerybinių pokyčių atskyrimas nuo kitų KV") %>%
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


#save table
gtsave(gt_table_tumor3,
       filename = "C:/Users/Ieva/rprojects/outputs_all/bo_roc_table_np_20250916.png")

#Combine the images
roc_image3<- image_read("C:/Users/Ieva/rprojects/outputs_all/bo_roc_np20250916.png")
table_image3 <- image_read("C:/Users/Ieva/rprojects/outputs_all/bo_roc_table_np_20250916.png")

combined_image3 <- image_append(c(roc_image3, table_image3), stack = F)

# Find the max width to align both
roc_info3 <- image_info(roc_image3)
table_info3 <- image_info(table_image3)
max_width3 <- max(roc_info3$width, table_info3$width)

# Pad each image to the max width
roc_image1_padded3 <- image_extent(roc_image3, geometry = geometry_area(max_width3, roc_info3$height), gravity = "center", color = "white")
table_image1_padded3 <- image_extent(table_image3, geometry = geometry_area(max_width3, table_info3$height), gravity = "center", color = "white")


combined_image3 <- image_append(c(roc_image1_padded3, table_image1_padded3), stack = T)

# Save the combined image
image_write(combined_image3, 
            "C:/Users/Ieva/rprojects/outputs_all/bo_roc_full_np_20250916.png")

#t test 2 groups #########################################
#melt table for expression
OC_table <- melt(KN_data_full[, c(15,32:35)],
                 id.vars="tumor",  measure.vars=raiska_np)
#Stujents t
stat.test_OC <- OC_table %>%
  group_by(variable) %>%
  summarise(
    p_value = t.test(value ~ tumor, var.equal = TRUE)$p.value
  ) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) #hes* and jag
stat.test_OC

#MANN-WHITNEY - Wilcoxon test (not normal)
wilcox.test_OC <- OC_table %>%
  group_by(variable) %>%
  summarise(
    p_value = wilcox.test(value ~ tumor)$p.value
  ) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) 
wilcox.test_OC #notch and ctnnb1


#Tribble 2 groups#####################
each.vs.ref_sig2 <-  tibble::tribble(
  ~group1, ~group2, ~p.adj,   ~y.position, ~variable,
  "Gerybiniai",   "KV", 0.0435, -2, "HES1", #stjudent's
)

#boxplot 2 groups ########################################

#rename to lt 
OC_table$tumor
OC_table2 <- OC_table %>%
  mutate(tumor = case_when(
    tumor == "OC" ~ "KV",
    tumor == "Benign" ~ "Gerybiniai"
  )) %>%
  mutate(variable  = case_when(
    variable       == "NOTCH2_DELTA" ~ "NOTCH2",
    variable       == "CTNNB1_DELTA" ~ "CTNNB1",
    variable       == "DLL1_DELTA" ~ "DLL1",
    variable       == "HES1_DELTA" ~ "HES1"
  ))


custom_colors <- c("KV" = "deeppink", "Gerybiniai" = "blue") 
OC_plot2 <- ggplot(OC_table2, aes(x=tumor , y=value, fill = variable)) +
  geom_boxplot( outlier.shape = NA , alpha=0.3, aes(fill = tumor )) +
  geom_jitter(aes(color = tumor ), size=1, alpha=0.5) +
  ylab(label = expression("Santykinė genų raiška, normalizuota pagal  " * italic("GAPDH"))) + 
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

OC_plot2

# Save the plot as a PNG file
png("C:/Users/Ieva/rprojects/outputs_all/np_oc_boxplot20250916.png",
    width = 1000, height = 1100, res = 200)
OC_plot2
dev.off()
#STAGE NORMALCY#########################################
#chek out stage
table(KN_data_full$Stage4, useNA = "a")
#normalcy
shapiro_results1 <- KN_data_full[, c(10,32:35)] %>%
  filter(!is.na(Stage4), Stage4 != "II")%>% #remove na and III due to low amount of values
  pivot_longer(cols = -Stage4, names_to = "gene", values_to = "value") %>%
  group_by(Stage4, gene) %>%
  summarise(p_value = shapiro.test(value)$p.value, .groups = "drop") 
shapiro_results1 #STAGE I, ctnnb1 is not  normal

#levene test instead, 3 groups
levene_rez2 <- KN_data_full[, c(10,32:35)] %>%
  pivot_longer(cols = -Stage4 , names_to = "variable", values_to = "value") %>%
  group_by(variable) %>%
  summarise(
    p_value = car::leveneTest(value ~ Stage4 )$`Pr(>F)`[1],
    .groups = "drop"
  )

levene_rez2 #all equal

#ANOVA stage####################################################
#DLL and ctnnb1 has to use Kruskal-Wallis test
#all equal vairances, othervise would use welch anova
#NOTCH2
res_aov_Stage <- aov( NOTCH2_DELTA ~ Stage4 ,
                      data = KN_data_full
)
shapiro.test(res_aov_Stage$residuals) #normality is met
summary(res_aov_Stage)
#HES1
res_aov2_stage <- aov( HES1_DELTA ~ Stage4 ,
                       data = KN_data_full
)
summary(res_aov2_stage) # p 0.0452 
shapiro.test(res_aov2_stage$residuals) #normality is met
# Tukey HSD test:
post_test_stage <- glht(res_aov2_stage, linfct = mcp(Stage4 = "Tukey"))

summary(post_test_stage) #hgsoc benign 0,0277

#not normal DLL1
kruskal1_stage <- kruskal.test(DLL1_DELTA ~ Stage4,
                               data = KN_data_full
)
kruskal1

#not normal CTNNB1
kruskal2_stage <- kruskal.test(CTNNB1_DELTA ~ Stage4,
                               data = KN_data_full
)
kruskal2_stage 

#Tribble stage groups#####################
each.vs.ref_sig2_stage <-  tibble::tribble(
  ~group1, ~group2, ~p.adj,   ~y.position, ~variable,
  "I",   "IV", 0.0452, -2, "HES1", #stjudent's
)

#STAGE boxplot ######################################
#melt table for expression
Stage_table <- melt(KN_data_full[, c(10,32:35)], id.vars="Stage4",  measure.vars=raiska_np)

Stage_table2 <- Stage_table %>%
  mutate(variable  = case_when(
    variable       == "NOTCH2_DELTA" ~ "NOTCH2",
    variable       == "CTNNB1_DELTA" ~ "CTNNB1",
    variable       == "DLL1_DELTA" ~ "DLL1",
    variable       == "HES1_DELTA" ~ "HES1"
  )) %>%
  filter(!is.na(Stage4))
custom_colors_stage <- c("IV" = "deeppink","III" = "lightpink","II" = "lightblue", "I" = "blue") 
STAGE_plot <- ggplot(Stage_table2, aes(x=Stage4 , y=value, fill = variable)) +
  geom_boxplot( outlier.shape = NA , alpha=0.3, aes(fill = Stage4 )) +
  geom_jitter(aes(color = Stage4 ), size=1, alpha=0.5) +
  ylab(label = expression("Santykinė genų raiška, normalizuota pagal  " * italic("GAPDH"))) + 
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



#Tribble stage groups#####################
each.vs.ref_sig2_stage <-  tibble::tribble(
  ~group1, ~group2, ~p.adj,   ~y.position, ~variable,
  "I",   "IV", 0.0452, -2, "HES1", #stjudent's
)

#STAGE boxplot ######################################
Stage_table2 <- Stage_table %>%
  mutate(variable  = case_when(
    variable       == "NOTCH2_DELTA" ~ "NOTCH2",
    variable       == "CTNNB1_DELTA" ~ "CTNNB1",
    variable       == "DLL1_DELTA" ~ "DLL1",
    variable       == "HES1_DELTA" ~ "HES1"
  )) %>%
  filter(!is.na(Stage4))
custom_colors_stage <- c("IV" = "deeppink","III" = "lightpink","II" = "lightblue", "I" = "blue") 
STAGE_plot <- ggplot(Stage_table2, aes(x=Stage4 , y=value, fill = variable)) +
  geom_boxplot( outlier.shape = NA , alpha=0.3, aes(fill = Stage4 )) +
  geom_jitter(aes(color = Stage4 ), size=1, alpha=0.5) +
  ylab(label = expression("Santykinė genų raiška, normalizuota pagal  " * italic("GAPDH"))) + 
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


# Save the plot as a PNG file
png("C:/Users/Ieva/rprojects/outputs_all/stage_plot_np20250919.png",
    width = 1000, height = 1100, res = 200)
STAGE_plot
dev.off()

#HGSOC only, stage################################
#make HGSOC df
HGSOC <- KN_data_full[c(KN_data_full$Grupė_Ieva == "HGSOC"),] 
table(HGSOC$Grupė_Ieva, HGSOC$Stage4)
#normalcy
shapiro_results1 <- HGSOC[, c(10,32:35)] %>%
  filter(!is.na(Stage4), Stage4 != "II")%>% #remove na and III due to low amount of values
  pivot_longer(cols = -Stage4, names_to = "gene", values_to = "value") %>%
  group_by(Stage4, gene) %>%
  summarise(p_value = shapiro.test(value)$p.value, .groups = "drop") 
shapiro_results1 #all normal

#levene test instead, 3 groups
levene_rez3 <- HGSOC[, c(10,32:35)] %>%
  pivot_longer(cols = -Stage4 , names_to = "variable", values_to = "value") %>%
  group_by(variable) %>%
  summarise(
    p_value = car::leveneTest(value ~ Stage4 )$`Pr(>F)`[1],
    .groups = "drop"
  )

levene_rez3 #all equal

#melt table for expression
Stage_table_HGSOC <- melt(HGSOC[, c(10,32:35)], id.vars="Stage4",  measure.vars=raiska_np)

#STAGE ANOVA, HGSOC ##########################
res_df_anova <- lapply(raiska_np, function(var) {
  formula <- as.formula(paste(var, "~ Stage4"))
  fit <- aov(formula, data = HGSOC)
  data.frame(
    variable = var,
    p_value = tidy(fit)[["p.value"]][1]  # ANOVA p-value
  )
}) %>% dplyr::bind_rows()

res_df_anova

#HGSOC STAGE boxplot ######################################
Stage_table3 <- Stage_table_HGSOC %>%
  mutate(variable  = case_when(
    variable       == "NOTCH2_DELTA" ~ "NOTCH2",
    variable       == "CTNNB1_DELTA" ~ "CTNNB1",
    variable       == "DLL1_DELTA" ~ "DLL1",
    variable       == "HES1_DELTA" ~ "HES1"
  )) %>%
  filter(!is.na(Stage4))
STAGE_plot_h <- ggplot(Stage_table3, aes(x=Stage4 , y=value, fill = variable)) +
  geom_boxplot( outlier.shape = NA , alpha=0.3, aes(fill = Stage4 )) +
  geom_jitter(aes(color = Stage4 ), size=1, alpha=0.5) +
  ylab(label = expression("Santykinė genų raiška, normalizuota pagal  " * italic("GAPDH"))) + 
  facet_wrap(.~ variable, nrow = 2, scales = "free") +
  #add_pvalue(each.vs.ref_sig, label = "p.adj") + #pvalue
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

STAGE_plot_h

# Save the plot as a PNG file
png("C:/Users/Ieva/rprojects/outputs_all/hgsoc_stage_plot_np20250916.png",
    width = 1000, height = 1100, res = 200)
STAGE_plot_h
dev.off()
#GRADE#####################################
table(KN_data_full$Grade2, useNA = "a")
#normalcy
shapiro_resultsg <- KN_data_full[, c(9,32:35)] %>%
  filter(!is.na(Grade2))%>% #remove na and III due to low amount of values
  pivot_longer(cols = -Grade2, names_to = "gene", values_to = "value") %>%
  group_by(Grade2, gene) %>%
  summarise(p_value = shapiro.test(value)$p.value, .groups = "drop") 
shapiro_resultsg #all normal 
#var test, grade
var_results_g<- KN_data_full[, c(9,32:35)] %>%
  filter(!is.na(Grade2))%>%
  pivot_longer(cols = -Grade2 , names_to = "variable", values_to = "value") %>%
  group_by(variable) %>%
  summarise(
    p_value = var.test(value[Grade2  == unique(Grade2 )[1]], 
                       value[Grade2  == unique(Grade2 )[2]])$p.value,
    .groups = "drop"
  ) 
var_results_g #all normal variance

#melt table for expression
Grade_table <- melt(KN_data_full[, c(9,32:35)], id.vars="Grade2",  measure.vars=raiska_np)

#stjundents test (normal, equal variances)
t.test_grade <- Grade_table %>%
  group_by(variable) %>%
  t_test(value ~ Grade2,
         p.adjust.method = "BH", 
         var.equal = TRUE, #stjudents
         paired = FALSE, 
         #detailed=TRUE 
  )%>%
  mutate(across(c(p.adj), ~ format(., scientific = FALSE))) %>% # Format p-values to remove scientific notation
  filter(p.adj < 0.1)
t.test_grade

#Tribble 2 groups for grade
each.vs.ref_sig2_g <-  tibble::tribble(
  ~group1, ~group2, ~p.adj,   ~y.position, ~variable,
  "G1",   "G3", 0.089, -2, "NOTCH2", #stjudent's
)

#GRADE boxplot ######################################
Grade_table2 <- Grade_table %>%
  mutate(variable  = case_when(
    variable       == "NOTCH2_DELTA" ~ "NOTCH2",
    variable       == "CTNNB1_DELTA" ~ "CTNNB1",
    variable       == "DLL1_DELTA" ~ "DLL1",
    variable       == "HES1_DELTA" ~ "HES1"
  )) %>%
  filter(!is.na(Grade2))
custom_colors_grade <- c("G1" = "deeppink","G3" = "lightblue") 
GRADE_plot <- ggplot(Grade_table2, aes(x=Grade2 , y=value, fill = variable)) +
  geom_boxplot( outlier.shape = NA , alpha=0.3, aes(fill = Grade2 )) +
  geom_jitter(aes(color = Grade2 ), size=1, alpha=0.5) +
  ylab(label = expression("Santykinė genų raiška, normalizuota pagal  " * italic("GAPDH"))) + 
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

GRADE_plot

# Save the plot as a PNG file
png("C:/Users/Ieva/rprojects/outputs_all/grade_plot_np20250916.png",
    width = 1000, height = 1100, res = 200)
GRADE_plot
dev.off()

#CORRELATION WITH AGE#########################
table(KN_data_full$Amžius, useNA ="a") #all full
age_table <- KN_data_full[, colnames(KN_data_full) %in% c(raiska_np, "Amžius")]
#normalcy
sapply(age_table, function(x) shapiro.test(x)$p.value) #cttnb1 not normal
#pearson corrwlation for all
results_n <- lapply(age_table[, colnames(age_table) %in% c("NOTCH2_DELTA", "DLL1_DELTA", "HES1_DELTA")], 
                    function(x) cor.test(x, age_table$Amžius, method = "pearson"))
results_n
#spearman for ctnnb1
cor.test(age_table$Amžius,age_table$CTNNB1_DELTA, method = "spearman")
#add results to the table
results_df <- bind_rows(lapply(names(results_n), function(var) {
  broom::tidy(results_n[[var]]) %>%
    mutate(variable = var)
}))

# Spearman correlation for CTNNB1
ctn_res <- broom::tidy(cor.test(age_table$Amžius, age_table$CTNNB1_DELTA, method = "spearman")) %>%
  mutate(variable = "CTNNB1_DELTA")
#join
common_cols <- intersect(names(results_df), names(ctn_res))
all_corr_results <- bind_rows(
  results_df %>% select(all_of(common_cols)),
  ctn_res %>% select(all_of(common_cols))
)
all_corr_results
#plot
age_table_l <- age_table %>%
  pivot_longer(cols = -Amžius, names_to = "variable", values_to = "value")

age_plot <- ggplot(age_table_l, aes(x = Amžius, y = value)) +
  geom_point(alpha = 0.6) +
  facet_wrap(~ variable, scales = "free_y") +
  geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") 
age_plot


# Save the plot as a PNG file
png("C:/Users/Ieva/rprojects/outputs_all/age_plot_np20250916.png",
    width = 1000, height = 1100, res = 200)
age_plot
dev.off()

#MODEL####################################

OC_HGSOC_BENIGN<- KN_data_full[c(KN_data_full$Grupė_Ieva != "Other"),] #51 cases left
OC_HGSOC_BENIGN$tumor <- relevel(factor(OC_HGSOC_BENIGN$Grupė_Ieva), ref = "Benign")

expr_np <- OC_HGSOC_BENIGN[colnames(OC_HGSOC_BENIGN) %in% raiska_np]
logistic.model <- glm(OC_HGSOC_BENIGN$tumor ~ ., data = expr_np, family = binomial("logit"))
predicted_probs <- predict.glm(logistic.model, type='response') 
pred_data <- OC_HGSOC_BENIGN[(rownames(OC_HGSOC_BENIGN) #remove incoplete rows
                              %in% names(predicted_probs)), ]
roc_curve <- roc(pred_data$tumor, predicted_probs)
AUC<- auc(roc_curve)
AUC
plot(roc_curve)
coords <- coords(roc_curve,
                 best.method= "closest.topleft","best",
                 ret=c("threshold", "accuracy", "sensitivity",
                       "specificity", "precision", "npv",
                       "tpr", "fpr"), transpose = FALSE)
coords


#roc MODEL figure BH#############################################
roc_plot <- function() {
  par(pty = "s") #sets square
  plot.roc(roc_results_np[["NOTCH2_DELTA"]], print.auc = F, col = "#dcbeff",
           cex.main=0.8, 
           main ="Gerybinių pokyčių atskyrimas nuo HGSOC",
           xlab = "Specifiškumas", 
           ylab = "Jautrumas") #title
  lines(roc_results_np[["CTNNB1_DELTA"]], col = "#911eb4", lwd =2) 
  lines(roc_results_np[["DLL1_DELTA"]], col ="#ffd8b1", lwd =2) 
  lines(roc_results_np[["HES1_DELTA"]], col = "#42d4f4", lwd =2) 
  lines(roc_curve, col = "darkred", lwd =2)
  #SUDĖTA NE PAGAL PAVEIKLĄ BET PAGAL AUC DYDI
  legend("bottomright", legend = c( expression(italic("NOTCH2")),
                                    expression(italic("CTNNB1")),
                                    expression(italic("DLL1")), 
                                    expression(italic("HES1")),
                                    expression(italic("4 gene model"))
  ),
  
  col = c("#dcbeff", "#911eb4", "#ffd8b1", "#42d4f4", 
          "darkred"), lty = 1, 
  cex = 0.5, lwd =3)
}
#plot
roc_plot()

# Save the plot as a PNG file
png("C:/Users/Ieva/rprojects/outputs_all/bh_model_roc_np20250916.png",
    width = 1000, height = 1000, res = 200)
roc_plot()
dev.off()

#roc MODEL table BH################################
#look at the roc features
coords_results_np
#create df
# If coords_results_np is a list of dataframes or vectors:
coords_all <- do.call(rbind, c(coords_results_np, list(coords)))

# Step 2: Build final dataframe
results_np <- data.frame(
  Predictor = c(raiska_np, "Combined"),
  AUC = c(auc_values_tumor, AUC),
  coords_all,
  stringsAsFactors = FALSE
)


#lithuanize it 
colnames(results_np) <- c("Biožymuo", "plotas po kreive", "slenkstinė vertė", 
                          "tikslumas", "jautrumas", "specifiškumas", 
                          "ppv", "npv", "tpr", "fpr")
rownames(results_np) <- c("NOCTH2", "CTNNB1", "DLL1", "HES1", "4 genų kombinacija")
results_np$Biožymuo <- c("NOCTH2", "CTNNB1", "DLL1", "HES1", "4 genų kombinacija")
#nice formating of the Table metrics for ROC OC
gt_table_np <- results_np %>%
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
gt_table_np


#save table
gtsave(gt_table_np,
       filename = "C:/Users/Ieva/rprojects/outputs_all/bh_roccombined_table_np_20250916.png")

#Combine the images
roc_image1<- image_read("C:/Users/Ieva/rprojects/outputs_all/bh_model_roc_np20250916.png")
table_image1 <- image_read("C:/Users/Ieva/rprojects/outputs_all/bh_roccombined_table_np_20250916.png")

combined_image1 <- image_append(c(roc_image1, table_image1), stack = F)

# Find the max width to align both
roc_info <- image_info(roc_image1)
table_info <- image_info(table_image1)
max_width <- max(roc_info$width, table_info$width)

# Pad each image to the max width
roc_image1_padded <- image_extent(roc_image1, geometry = geometry_area(max_width, roc_info$height), gravity = "center", color = "white")
table_image1_padded <- image_extent(table_image1, geometry = geometry_area(max_width, table_info$height), gravity = "center", color = "white")


combined_image1 <- image_append(c(roc_image1_padded, table_image1_padded), stack = T)

# Save the combined image
image_write(combined_image1, 
            "C:/Users/Ieva/rprojects/outputs_all/bh_roc_model_full_np_20250916.png")

#PAIRED TEST THE TISSUE VS NP
results_list <- lapply(seq_along(raiska_np), function(i) {
  
  plasma_col <- raiska_np[i]
  tissue_col <- raiska_tissue[i]
  gene_name  <- plasma_col  # or tissue_col, just for labeling
  
  # Compute differences
  diffs <- KN_data_full[[tissue_col]] - KN_data_full[[plasma_col]]
  
  # Check normality
  normal_p <- shapiro.test(diffs)$p.value
  
  if(normal_p > 0.05) {
    # paired t-test
    test <- t.test(KN_data_full[[tissue_col]], KN_data_full[[plasma_col]], paired = TRUE)
    data.frame(
      gene = gene_name,
      test = "paired t-test",
      statistic = test$statistic,
      p_value = test$p.value,
      mean_diff = test$estimate,
      normality_p = normal_p
    )
  } else {
    # Wilcoxon signed-rank test
    test <- wilcox.test(KN_data_full[[tissue_col]], KN_data_full[[plasma_col]], paired = TRUE)
    data.frame(
      gene = gene_name,
      test = "Wilcoxon signed-rank",
      statistic = test$statistic,
      p_value = test$p.value,
      mean_diff = median(diffs),   # report median difference for Wilcoxon
      normality_p = normal_p
    )
  }
})

# Combine into a single data frame
results_df <- do.call(rbind, results_list)

# Optional: adjust p-values for multiple testing
results_df$p_adj <- p.adjust(results_df$p_value, method = "BH")

results_df

#PAIRED PLOTS #####################
# Make sure plasma and tissue columns are in the same order
plot_df <- lapply(seq_along(raiska_np), function(i) {
  plasma_col <- raiska_np[i]
  tissue_col <- raiska_tissue[i]
  gene_name  <- plasma_col  # or tissue_col
  
  data.frame(
    patient = 1:nrow(KN_data_full),
    expression = c(KN_data_full[[tissue_col]], KN_data_full[[plasma_col]]),
    source = rep(c("Tissue","NP"), each = nrow(KN_data_full)),
    gene = gene_name
  )
}) %>% bind_rows()


# Clean gene names: remove "_DELTA"
plot_df <- plot_df %>%
  mutate(gene_clean = str_replace(gene, "_DELTA", ""))

p_labels <- results_df %>%
  dplyr::select(gene, p_adj) %>%
  dplyr::mutate(
    gene_clean = str_replace(gene, "_DELTA", ""),
    label = ifelse(p_adj < 0.001, "p < 0.001", paste0("p = ", signif(p_adj, 3)))
  )

# Y-position for p-value labels
y_positions <- plot_df %>%
  group_by(gene_clean) %>%
  summarize(y_pos = max(expression, na.rm = TRUE) * 1.05, .groups = "drop")

p_labels <- left_join(p_labels, y_positions, by = "gene_clean")

# Convert gene names to expressions for italics
# Plot
paired_plot <- ggplot(plot_df, aes(x = source, y = expression, group = patient)) +
  geom_point(alpha = 0.6) +
  geom_line(alpha = 0.5, color = "gray") +
  facet_wrap(~ gene_clean, scales = "free_y") +
  geom_text(
    data = p_labels,
    aes(x = 1.5, y = y_pos, label = label),
    inherit.aes = FALSE,
    size = 3
  ) +
  theme_bw() +
  labs(y = "Expression", x = "", title = "Paired Tissue vs NP Expression")

paired_plot
# Save the plot as a PNG file
png("C:/Users/Ieva/rprojects/outputs_all/paired_plot20250928.png",
    width = 1000, height = 1100, res = 200)
paired_plot
dev.off()
