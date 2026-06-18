#KN-  liquid 2026 05 12, 13
#NP, URINE, PLASMA only the 15 cases data
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
library(coin)
library(patchwork)
library(effectsize)
library(scales)
#read RDS
LIQUID_DF_final <- readRDS("C:/Users/Ieva/rprojects/OTHER DATA/KN_LIQUID/liquid_20260415.RDS")
#leave URINE/PLASMA only
#remove cases that is NA in their type: 
LIQUID_DF_15 <- LIQUID_DF_final %>%
  filter(!is.na(HES1_URINE))
table(LIQUID_DF_15$TYPE, useNA = "a")
#BOXPLOT URINE#######################################
#URINE HGSOC vs others######################
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
#boxplots URINE ############################
#melt table for expression
URINE_table <- melt(LIQUID_DF_15[, c(12,19:22)], id.vars="TYPE",
                    measure.vars=c("NOTCH2_URINE","CTNNB1_URINE","DLL1_URINE","HES1_URINE"))
each.vs.ref_sig2_g <- tibble::tribble(
  ~group1, ~group2, ~p.adj,   ~y.position, ~variable,
  "HGSOC",   "OTHER", 0.6513, -4, "NOTCH2_URINE",#
  #"HGSOC",   "OTHER", 0.001, -1, "DLL1_URINE",#
  "HGSOC",   "OTHER", 0.6976 , 1, "HES1_URINE",#dunn
  "HGSOC",   "OTHER", 0.7676, -2, "CTNNB1_URINE"#dunn

)
#TUMOR boxplot
custom_colors_grade <- c("HGSOC" = "lightpink","OTHER" = "darkblue") 
URINE_plot <- ggplot(URINE_table, aes(x=TYPE , y=value, fill = variable)) +
  geom_boxplot( outlier.shape = NA , alpha=0.3, aes(fill = TYPE )) +
  geom_jitter(aes(color = TYPE ), size=1, alpha=0.5) +
  ylab(label = expression("Gene expression, normalized to  " * italic("GAPDH"))) + 
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

URINE_plot
#save urine 15 cases boxplot
ggsave("C:/Users/Ieva/rprojects/outputs_all/LIQUID/15_urine_20260512.png",
       plot = URINE_plot,
       width = 12,
       height = 16,
       units = "cm",
       dpi = 400)

#remove HES1 > 0 case#########################################
HES1_df <- LIQUID_DF_15[LIQUID_DF_15[,1] != "KN-108", ]
#urine normalcy DLL1 special
by(HES1_df$HES1_URINE, HES1_df$TYPE, shapiro.test)#0.2706
#urine variance
car::leveneTest(HES1_df$HES1_URINE ~ HES1_df$TYPE, center = median)#0.2289
#all normal and variance normal
#t.test URINE
t.test(HES1_df$HES1_URINE  ~ HES1_df$TYPE, var.equal = TRUE) #0.4065
t.test(HES1_df$HES1_URINE  ~ HES1_df$TYPE) #0.359 if not assuming equal variances
independence_test(HES1_URINE ~ TYPE, #also checks without assumptions
                  data = HES1_df,
                  distribution = "exact")#0.4006
#for other groups also does not affect
independence_test(CTNNB1_URINE ~ TYPE, #also checks without assumptions
                  data = LIQUID_DF_15,
                  distribution = "exact")#0.7888
independence_test(NOTCH2_URINE ~ TYPE, #also checks without assumptions
                  data = LIQUID_DF_15,
                  distribution = "exact")#0.65

#boxplots URINE 14 case############################
#melt table for expression
URINE_table2 <- melt(HES1_df[, c(12,19:22)], id.vars="TYPE",
                    measure.vars=c("NOTCH2_URINE","CTNNB1_URINE","DLL1_URINE","HES1_URINE"))
each.vs.ref_sig1_g <- tibble::tribble(
  ~group1, ~group2, ~p.adj,   ~y.position, ~variable,
  "HGSOC",   "OTHER", 0.6513, -4, "NOTCH2_URINE",#
  #"HGSOC",   "OTHER", 0.001, -1, "DLL1_URINE",#
  "HGSOC",   "OTHER", 0.4065 , -1, "HES1_URINE",#dunn
  "HGSOC",   "OTHER", 0.7676, -2, "CTNNB1_URINE"#dunn
  
)
#TUMOR boxplot
custom_colors_grade <- c("HGSOC" = "lightpink","OTHER" = "darkblue") 
URINE_plot2 <- ggplot(URINE_table2, aes(x=TYPE , y=value, fill = variable)) +
  geom_boxplot( outlier.shape = NA , alpha=0.3, aes(fill = TYPE )) +
  geom_jitter(aes(color = TYPE ), size=1, alpha=0.5) +
  ylab(label = expression("Gene expression, normalized to  " * italic("GAPDH"))) + 
  facet_wrap(.~ variable, nrow = 2, scales = "free") +
  facet_wrap(
    .~ variable,
    nrow = 2,
    scales = "free",
    labeller = labeller(
      variable = c(
        "CTNNB1_URINE" = "CTNNB1 expression in urine",
        "DLL1_URINE" = "DLL1 expression in urine",
        "HES1_URINE" = "HES1 expression in urine",
        "NOTCH2_URINE" = "NOTCH2 expression in urine"
      )
    )
  ) +
  add_pvalue(each.vs.ref_sig1_g, label = "p.adj") +
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

URINE_plot2
#save hes1, 14 cases boxplot
ggsave("C:/Users/Ieva/rprojects/outputs_all/LIQUID/14_urine_20260521.png",
       plot = URINE_plot2,
       width = 15,
       height = 16,
       units = "cm",
       dpi = 400)
#nothing much changes, so this is not used

#remove high CTNNB1 values############################
ctnnb1_df <- LIQUID_DF_15[!LIQUID_DF_15$`Laboratorinis kodas` %in% c( "KN-106", "KN-112"), ]
#urine normalcy CTNNB1 special
by(ctnnb1_df$CTNNB1_URINE, ctnnb1_df$TYPE, shapiro.test)#0.9347
#urine variance
car::leveneTest(ctnnb1_df$CTNNB1_URINE ~ ctnnb1_df$TYPE, center = median)#0.9807
#all normal and variance normal
#t.test URINE
t.test(ctnnb1_df$CTNNB1_URINE  ~ ctnnb1_df$TYPE, var.equal = TRUE) #0.1285 
#not significant still
#FISHER ####################################################
#make factor of gene expresssion
DATA_names_P <- c("NOTCH2_NP","CTNNB1_NP","DLL1_NP","HES1_NP",
                  "NOTCH2_URINE", "CTNNB1_URINE","DLL1_URINE","HES1_URINE",
                  "NOTCH2_P_norm","CTNNB1_P_norm", "DLL1_P_norm","HES1_P_norm")
#make factors in liquid df 15
table(LIQUID_DF_15$TYPE)
LIQUID_DF_15$TYPE <- droplevels(LIQUID_DF_15$TYPE)
LIQUID_DF_15 <- LIQUID_DF_15 %>%
  mutate(
    across(all_of(DATA_names_P),
           ~ factor(if_else(. > median(., na.rm = TRUE), "High", "Low")),
           .names = "{.col}_f")
  )
#fisher tests on urine#####################################
fisher.test(table(LIQUID_DF_15$TYPE, LIQUID_DF_15$NOTCH2_URINE_f))#1
fisher.test(table(LIQUID_DF_15$TYPE, LIQUID_DF_15$CTNNB1_URINE_f))#0.2657
fisher.test(table(LIQUID_DF_15$TYPE, LIQUID_DF_15$HES1_URINE_f))#0.6084
fisher.test(table(LIQUID_DF_15$TYPE, LIQUID_DF_15$DLL1_URINE_f))#1
#fisher tests on plasma#########################################
fisher.test(table(LIQUID_DF_15$TYPE, LIQUID_DF_15$NOTCH2_P))#0.3333
fisher.test(table(LIQUID_DF_15$TYPE, LIQUID_DF_15$CTNNB1_P))#0.2418
fisher.test(table(LIQUID_DF_15$TYPE, LIQUID_DF_15$HES1_P))#1
#fisher.test(table(LIQUID_DF_15$TYPE, LIQUID_DF_15$DLL1_P))#NOT ENOUGH DATA

#plot fishers plasma#######################
##ctnbb1 plasma barplot##########################
tab_ctnnb1 <- table(LIQUID_DF_15$TYPE, LIQUID_DF_15$CTNNB1_P)
colnames(tab_ctnnb1) <- c("CTNNB1-", "CTNNB1+")

ft_ctnnb1 <- fisher.test(tab_ctnnb1)

df_ctnnb1 <- as.data.frame(tab_ctnnb1)
colnames(df_ctnnb1) <- c("TYPE", "CTNNB1", "Freq")

# precise proportions (no rounding here)
df_ctnnb1 <- df_ctnnb1 %>%
  group_by(TYPE) %>%
  mutate(prop = Freq / sum(Freq))
ctnnb1_bar <- ggplot(df_ctnnb1, aes(x = TYPE, y = prop, fill = CTNNB1)) +
  
  geom_bar(stat = "identity", width = 0.6) +
  
  geom_text(
    aes(label = percent(prop, accuracy = 0.1)),
    position = position_stack(vjust = 0.5),
    size = 4
  ) +
  
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  
  scale_fill_manual(values = c(
    "CTNNB1-" = "grey80",
    "CTNNB1+" = "steelblue"
  )) +
  
  labs(
    x = NULL,
    y = "Proportion of cases",
    fill = NULL,
    title = expression(
      paste(italic("CTNNB1"), " status in plasma by OC type")
    ),
    subtitle = paste0("Fisher p = ", signif(ft_ctnnb1$p.value, 3))
  ) +
  
  theme_minimal(base_size = 14) +
  
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "top"
  )
ctnnb1_bar
#save ctnnb1 barplot plasma
ggsave("C:/Users/Ieva/rprojects/outputs_all/LIQUID/ctnnb1_plasma_20260512.png",
       plot = ctnnb1_bar,
       width = 14,
       height = 10,
       units = "cm",
       dpi = 200)
##notch2 plasma barplot####################
tab_notch2 <- table(LIQUID_DF_15$TYPE, LIQUID_DF_15$NOTCH2_P)
colnames(tab_notch2) <- c("NOTCH2-", "NOTCH2+")

ft_notch2 <- fisher.test(tab_notch2)

df_notch2 <- as.data.frame(tab_notch2)
colnames(df_notch2) <- c("TYPE", "NOTCH2", "Freq")

# precise proportions (no rounding here)
df_notch2 <- df_notch2 %>%
  group_by(TYPE) %>%
  mutate(prop = Freq / sum(Freq))
notch2_bar <- ggplot(df_notch2, aes(x = TYPE, y = prop, fill = NOTCH2)) +
  
  geom_bar(stat = "identity", width = 0.6) +
  
  geom_text(
    aes(label = ifelse(prop > 0,
                       percent(prop, accuracy = 0.1),
                       "")),
    position = position_stack(vjust = 0.5),
    size = 4
  ) +
  
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  
  scale_fill_manual(values = c(
    "NOTCH2-" = "grey80",
    "NOTCH2+" = "steelblue"
  )) +
  
  labs(
    x = NULL,
    y = "Proportion of cases",
    fill = NULL,
    title = expression(
      paste(italic("NOTCH2"), " status in plasma by OC type")
    ),
    subtitle = paste0("Fisher p = ", signif(ft_notch2$p.value, 3))
  ) +
  
  theme_minimal(base_size = 14) +
  
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "top"
  )
notch2_bar
#save notch2 barplot plasma
ggsave("C:/Users/Ieva/rprojects/outputs_all/LIQUID/notxh2_plasma_20260512.png",
       plot = notch2_bar,
       width = 14,
       height = 10,
       units = "cm",
       dpi = 200)
##hes1 plasma barplot####################
tab_hes1 <- table(LIQUID_DF_15$TYPE, LIQUID_DF_15$HES1_P)
colnames(tab_hes1) <- c("HES1-", "HES1+")

ft_hes1 <- fisher.test(tab_hes1)

df_hes1 <- as.data.frame(tab_hes1)
colnames(df_hes1) <- c("TYPE", "HES1", "Freq")

# precise proportions (no rounding here)
df_hes1 <- df_hes1 %>%
  group_by(TYPE) %>%
  mutate(prop = Freq / sum(Freq))
hes1_bar <- ggplot(df_hes1, aes(x = TYPE, y = prop, fill = HES1)) +
  
  geom_bar(stat = "identity", width = 0.6) +
  
  geom_text(
    aes(label = ifelse(prop > 0,
                       percent(prop, accuracy = 0.1),
                       "")),
    position = position_stack(vjust = 0.5),
    size = 4
  ) +
  
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  
  scale_fill_manual(values = c(
    "HES1-" = "grey80",
    "HES1+" = "steelblue"
  )) +
  
  labs(
    x = NULL,
    y = "Proportion of cases",
    fill = NULL,
    title = expression(
      paste(italic("HES1"), " status in plasma by OC type")
    ),
    subtitle = paste0("Fisher p = ", signif(ft_hes1$p.value, 3))
  ) +
  
  theme_minimal(base_size = 14) +
  
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "top"
  )
hes1_bar
#save hes1 barplot plasma
ggsave("C:/Users/Ieva/rprojects/outputs_all/LIQUID/hes1_plasma_20260512.png",
       plot = hes1_bar,
       width = 14,
       height = 10,
       units = "cm",
       dpi = 200)
#import images
ctnnb1_image   <- image_read("C:/Users/Ieva/rprojects/outputs_all/LIQUID/ctnnb1_plasma_20260512.png")
hes1_image <- image_read("C:/Users/Ieva/rprojects/outputs_all/LIQUID/hes1_plasma_20260512.png")
notch2_image <- image_read("C:/Users/Ieva/rprojects/outputs_all/LIQUID/notxh2_plasma_20260512.png")


# combine vertically
combined <- image_append(c(ctnnb1_image, hes1_image, notch2_image), stack = TRUE)

# save
image_write(combined,
            "C:/Users/Ieva/rprojects/outputs_all/LIQUID/barplots_combined_plasma_20260612.png")
#ROC#################################################
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
           main ="Urine biomarkers HGSOC vs Other OC",
           # xlab = "1 - Specifiškumas", 
           # ylab = "Jautrumas", 
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
#change names
rownames(results_URINE_15) <- c("NOTCH2", "CTNNB1", "DLL1", "HES1")
results_URINE_15$Predictor <- c("NOTCH2", "CTNNB1", "DLL1", "HES1")
#nice formating of the Table metrics for ROC OC
gt_table_URINE_15 <- results_URINE_15 %>%
  gt() %>%
  tab_header(
    title = "ROC measures for urine biomarkers", 
    subtitle = "HGOSC vs other OC") %>%
  fmt_number(
    columns = everything(),
    decimals = 3
  ) %>%
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_body(columns = vars(Predictor))
  )
#show
gt_table_URINE_15

# Save the plot as a PNG file
png("C:/Users/Ieva/rprojects/outputs_all/LIQUID/urine_OC_other_ROC_20260512.png",
    width = 11, height = 11, res = 300, units = "cm")
roc_plot2()
# mtext(
#   "B",
#   side = 3,
#   adj = -0.2,
#   line = 1,
#   cex = 1.2,
#   font = 2
# )
dev.off()
#there is no other convieneat way to save gt outputs
gtsave(gt_table_URINE_15,vwidth = 10000,   
       filename = "C:/Users/Ieva/rprojects/outputs_all/LIQUID/urine_OC_other_ROCtable_20260512.png")
#import images
roc_imageOC_urine  <- image_read("C:/Users/Ieva/rprojects/outputs_all/LIQUID/urine_OC_other_ROC_20260512.png")
table_imageOC_urine <- image_read("C:/Users/Ieva/rprojects/outputs_all/LIQUID/urine_OC_other_ROCtable_20260512.png")

# resize table to match ROC image width
table_imageOC_urine <- image_resize(table_imageOC_urine,
                                    paste0(image_info(roc_imageOC_urine)$width, "x"))
# combine vertically
combinedOC_urine <- image_append(c(roc_imageOC_urine, table_imageOC_urine), stack = TRUE)

# save
image_write(combinedOC_urine,
            "C:/Users/Ieva/rprojects/outputs_all/LIQUID/urine_OC_other_ROCcombined_20260521.png")

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
##urine cox###############################
cox_model_notch2_urine <- coxph(Surv(OS, STATUS) ~ NOTCH2_URINE_f, data = ALL_SURV_LIQUID)
summary(cox_model_notch2_urine)

cox_model_CTNNB1_urine <- coxph(Surv(OS, STATUS) ~ CTNNB1_URINE_f, data = ALL_SURV_LIQUID)
summary(cox_model_CTNNB1_urine)

cox_model_DLL1_urine <- coxph(Surv(OS, STATUS) ~ DLL1_URINE_f, data = ALL_SURV_LIQUID)
summary(cox_model_DLL1_urine)

cox_model_HES1_urine <- coxph(Surv(OS, STATUS) ~ HES1_URINE_f, data = ALL_SURV_LIQUID)
summary(cox_model_HES1_urine) 
## PLOTS KM URINE#################
#ctnnb1
ALL_SURV_LIQUID$CTNNB1 <- ALL_SURV_LIQUID$CTNNB1_URINE_f
u_ctnnb1 <- ggsurvplot(survfit(Surv(OS, STATUS) ~ CTNNB1, data = ALL_SURV_LIQUID), 
                       data = ALL_SURV_LIQUID, pval = F,
                       xlab = "Time (months)",
                       palette = c(  "#E64164", "#002060"),  
                       title =expression( "Overall survival by " * italic("CTNNB1") * " expression in urine"))
u_ctnnb1$plot <- u_ctnnb1$plot + labs(subtitle = "p = 0.3")
print(u_ctnnb1)
#notch2
ALL_SURV_LIQUID$NOTCH2 <- ALL_SURV_LIQUID$NOTCH2_URINE_f
u_notch2 <- ggsurvplot(survfit(Surv(OS, STATUS) ~ NOTCH2, data = ALL_SURV_LIQUID), 
                       data = ALL_SURV_LIQUID, pval = F,
                       xlab = "Time (months)",
                       palette = c(  "#E64164", "#002060"), 
                       title =expression( "Overall survival by " * italic("NOTCH2") * " expression in urine"))
u_notch2$plot <- u_notch2$plot + labs(subtitle = "p = 0.4")
print(u_notch2)
#DLL1
ALL_SURV_LIQUID$DLL1 <- ALL_SURV_LIQUID$DLL1_URINE_f
u_dll1 <- ggsurvplot(survfit(Surv(OS, STATUS) ~ DLL1, data = ALL_SURV_LIQUID), 
                     data = ALL_SURV_LIQUID, pval = F,
                     xlab = "Time (months)",
                     palette = c(  "#E64164", "#002060"), 
                     title =expression( "Overall survival by " * italic("DLL1") * " expression in urine"))
u_dll1$plot <- u_dll1$plot + labs(subtitle = "p = 1.00")
print(u_dll1)
#HES1
ALL_SURV_LIQUID$HES1<- ALL_SURV_LIQUID$HES1_URINE_f
u_hes1 <- ggsurvplot(survfit(Surv(OS, STATUS) ~ HES1, data = ALL_SURV_LIQUID), 
                     data = ALL_SURV_LIQUID, pval = F,
                     xlab = "Time (months)",
                     palette = c(  "#E64164", "#002060"), 
                     title =expression( "Overall survival by " * italic("HES1") * " expression in urine"))
u_hes1$plot <- u_hes1$plot + labs(subtitle = "p = 0.96")
print(u_hes1)

combined_plot_u <-
  (u_hes1$plot |  u_notch2$plot) / 
  (u_ctnnb1$plot |u_dll1 )

ggsave(
  filename = "c:/Users/Ieva/rprojects/outputs_all/LIQUID/OC_URINE_survival_combined20260521.png",
  plot = combined_plot_u,
  width = 12,
  height = 8,
  dpi = 400
)

#no dll1
combined_plot_u3 <-
  u_hes1$plot/  u_notch2$plot /u_ctnnb1$plot 
combined_plot_u3

ggsave(
  filename = "c:/Users/Ieva/rprojects/outputs_all/LIQUID/OC_URINE_survival_3genes_combined20260618.png",
  plot = combined_plot_u3,
  width = 6,
  height = 12,
  dpi = 600
)
##SURVIVAL PLASMA##########################################
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
ALL_SURV_LIQUID_P <- ALL_SURV_LIQUID <- LIQUID_DF_15 %>%
  filter(!is.na(OS), !is.na(STATUS)) #58 observations
#NOTCH2
ALL_SURV_LIQUID_P$NOTCH2 <- factor(
  ALL_SURV_LIQUID_P$NOTCH2_P,
  levels = c(0, 1),
  labels = c("no expression", "expression")
)
p_notch2 <- ggsurvplot(survfit(Surv(OS, STATUS) ~ NOTCH2, data = ALL_SURV_LIQUID_P), 
                       data = ALL_SURV_LIQUID_P, pval = F,xlab = "Time (months)",
                       palette = c( "#002060", "#E64164"), 
                       title =expression( "Overall survival by " * italic("NOTCH2") * " expression in plasma"))
p_notch2$plot <- p_notch2$plot + labs(subtitle = "Log-rank  p = 0.7")
print(p_notch2)
#HES1
ALL_SURV_LIQUID_P$HES1 <- factor(
  ALL_SURV_LIQUID_P$HES1_P,
  levels = c(0, 1),
  labels = c("no expression", "expression")
)
p_hes1 <-ggsurvplot(survfit(Surv(OS, STATUS) ~ HES1, data = ALL_SURV_LIQUID_P), 
                    data = ALL_SURV_LIQUID_P, pval = F,xlab = "Time (months)",
                    palette = c( "#002060", "#E64164"), 
                    title =expression( "Overall survival by " * italic("HES1") * " expression in plasma"))
p_hes1$plot <- p_hes1$plot + labs(subtitle = "Log-rank  p = 0.016")
print(p_hes1)
#DLL1 - not enough data
#CTNNB1
ALL_SURV_LIQUID_P$CTNNB1 <- factor(
  ALL_SURV_LIQUID_P$CTNNB1_P,
  levels = c(0, 1),
  labels = c("no expression", "expression")
)
p_ctnnb1 <- ggsurvplot(survfit(Surv(OS, STATUS) ~ CTNNB1, data = ALL_SURV_LIQUID_P), 
                       data = ALL_SURV_LIQUID_P, pval = F,xlab = "Time (months)",
                       palette = c( "#002060", "#E64164"), 
                       title =expression( "Overall survival by " * italic("CTNNB1") * " expression in plasma"))

p_ctnnb1$plot <- p_ctnnb1$plot + labs(subtitle = "p = 0.45")
print(p_ctnnb1)

combined_plot <-
  p_hes1$plot /  p_notch2$plot / p_ctnnb1$plot 

ggsave(
  filename = "c:/Users/Ieva/rprojects/outputs_all/LIQUID/OC_PLASMA_survival_combined2026021.png",
  plot = combined_plot,
  width = 6,
  height = 10,
  dpi = 300
)

#CLINICALS?##############################
table(LIQUID_DF_15$Grade_simple, useNA = "a") #g1 only 2 cases
table(LIQUID_DF_15$Stage_simple, useNA = "a") #g1 only 2 cases

LIQUID_DF_15$Stage_grouped <- ifelse(
  LIQUID_DF_15$Stage_simple %in% c(1, 2),
  "1_2",
  as.character(LIQUID_DF_15$Stage_simple)
)
table(LIQUID_DF_15$Stage_grouped, useNA = "a") #g1 only 2 cases

table(LIQUID_DF_15$Age, useNA = "a") #all have age
table(LIQUID_DF_15$CA125_num, useNA = "a") #4 is na -> 11 have data
#make factors
LIQUID_DF_15$Stage_grouped <- factor(LIQUID_DF_15$Stage_grouped)
LIQUID_DF_15$Grade_simple <- factor(LIQUID_DF_15$Grade_simple)
##plasma clinicals fishers###################################
#grade
fisher.test(table(LIQUID_DF_15$Grade_simple, LIQUID_DF_15$NOTCH2_P))#1
fisher.test(table(LIQUID_DF_15$Grade_simple, LIQUID_DF_15$CTNNB1_P))#0.3778
fisher.test(table(LIQUID_DF_15$Grade_simple, LIQUID_DF_15$HES1_P))#1
#fisher.test(table(LIQUID_DF_15$Grade_simple, LIQUID_DF_15$DLL1_P))#NOT ENOUGH DATA
#stage
fisher.test(table(LIQUID_DF_15$Stage_grouped, LIQUID_DF_15$NOTCH2_P))#1
fisher.test(table(LIQUID_DF_15$Stage_grouped, LIQUID_DF_15$CTNNB1_P))#0.1538
fisher.test(table(LIQUID_DF_15$Stage_grouped, LIQUID_DF_15$HES1_P))#0.5165
#fisher.test(table(LIQUID_DF_15$Stage_grouped, LIQUID_DF_15$DLL1_P))#NOT ENOUGH DATA
#plasma clinicals normalcy###############################
#age normalcy
#by(LIQUID_DF_15$Age, LIQUID_DF_15$NOTCH2_P, shapiro.test)#only 1 case not enough data
by(LIQUID_DF_15$Age, LIQUID_DF_15$DLL1_P, shapiro.test) #0.8039
by(LIQUID_DF_15$Age, LIQUID_DF_15$HES1_P, shapiro.test)#0.001241 - not normal
by(LIQUID_DF_15$Age, LIQUID_DF_15$CTNNB1_P, shapiro.test)#0.3839
# variance
car::leveneTest(LIQUID_DF_15$Age ~ LIQUID_DF_15$NOTCH2_P, center = median)#0.1842
car::leveneTest(LIQUID_DF_15$Age ~ LIQUID_DF_15$DLL1_P, center = median) #not enough data
car::leveneTest(LIQUID_DF_15$Age ~ LIQUID_DF_15$HES1_P, center = median)#0.1536
car::leveneTest(LIQUID_DF_15$Age ~ LIQUID_DF_15$CTNNB1_P, center = median)#0.7626
#HES1 not normal, others -normal

#CA125 normalcy
LIQUID_DF_15$CA125 <- as.numeric(LIQUID_DF_15$CA125)
by(LIQUID_DF_15$CA125 , LIQUID_DF_15$NOTCH2_P, shapiro.test)#only 1 case not enough data
by(LIQUID_DF_15$CA125, LIQUID_DF_15$DLL1_P, shapiro.test) #0.003077- not normal
by(LIQUID_DF_15$CA125, LIQUID_DF_15$HES1_P, shapiro.test)#0.01031 - not normal
by(LIQUID_DF_15$CA125, LIQUID_DF_15$CTNNB1_P, shapiro.test)#0.3888
# variance
car::leveneTest(LIQUID_DF_15$CA125 ~ LIQUID_DF_15$NOTCH2_P, center = median)#0.4777
car::leveneTest(LIQUID_DF_15$CA125 ~ LIQUID_DF_15$DLL1_P, center = median) #not enough data
car::leveneTest(LIQUID_DF_15$CA125 ~ LIQUID_DF_15$HES1_P, center = median)#0.9175
car::leveneTest(LIQUID_DF_15$CA125 ~ LIQUID_DF_15$CTNNB1_P, center = median)#0.2435
#DLL1 and HES1 not normal
#t test plasma clinical #########################
#age
t.test(LIQUID_DF_15$Age ~ LIQUID_DF_15$NOTCH2_P, var.equal = TRUE)#0.6729
t.test(LIQUID_DF_15$Age  ~ LIQUID_DF_15$DLL1_P, var.equal = TRUE) #not enough
wilcox.test(LIQUID_DF_15$Age  ~ LIQUID_DF_15$HES1_P) #0.8957
t.test(LIQUID_DF_15$Age  ~ LIQUID_DF_15$CTNNB1_P, var.equal = TRUE)#0.9599

#CA125 
t.test(LIQUID_DF_15$CA125 ~ LIQUID_DF_15$NOTCH2_P, var.equal = TRUE)#0.4743
wilcox.test(LIQUID_DF_15$CA125  ~ LIQUID_DF_15$DLL1_P) #not enough
wilcox.test(LIQUID_DF_15$CA125  ~ LIQUID_DF_15$HES1_P) #0.3758
t.test(LIQUID_DF_15$CA125  ~ LIQUID_DF_15$CTNNB1_P, var.equal = TRUE)#0.253

#urine clinicals normalcy###############################
#age/ca125 normalcy
shapiro.test(LIQUID_DF_15$Age)#0.8039
shapiro.test(LIQUID_DF_15$CA125)#0.003077- not normal
shapiro.test(LIQUID_DF_15$NOTCH2_URINE)#0.9683
shapiro.test(LIQUID_DF_15$DLL1_URINE)#0.4106
shapiro.test(LIQUID_DF_15$HES1_URINE)#0.6215
shapiro.test(LIQUID_DF_15$CTNNB1_URINE)#0.5101
#urine correlations###################################
#age
cor.test(LIQUID_DF_15$Age, LIQUID_DF_15$HES1_URINE, method = "pearson")#0.7714
cor.test(LIQUID_DF_15$Age, LIQUID_DF_15$NOTCH2_URINE, method = "pearson")#0.08791
cor.test(LIQUID_DF_15$Age, LIQUID_DF_15$CTNNB1_URINE, method = "pearson")#0.9567
cor.test(LIQUID_DF_15$Age, LIQUID_DF_15$DLL1_URINE, method = "pearson")#0.5706
#ca125
cor.test(LIQUID_DF_15$CA125, LIQUID_DF_15$HES1_URINE, method = "spearman")#0.2606
cor.test(LIQUID_DF_15$Age, LIQUID_DF_15$NOTCH2_URINE, method = "spearman")#0.09279
cor.test(LIQUID_DF_15$Age, LIQUID_DF_15$CTNNB1_URINE, method = "spearman")#0.7003
cor.test(LIQUID_DF_15$Age, LIQUID_DF_15$DLL1_URINE, method = "spearman")#1

#STAGE urine normalcy
by(LIQUID_DF_15$NOTCH2_URINE, LIQUID_DF_15$Stage_grouped, shapiro.test)#not enough data
by(LIQUID_DF_15$DLL1_URINE, LIQUID_DF_15$Stage_grouped, shapiro.test) #Not enough data
by(LIQUID_DF_15$HES1_URINE, LIQUID_DF_15$Stage_grouped, shapiro.test)#0.9978/0.5858
by(LIQUID_DF_15$CTNNB1_URINE, LIQUID_DF_15$Stage_grouped, shapiro.test)#0.1502/0.8253
#urine variance
car::leveneTest(LIQUID_DF_15$NOTCH2_URINE ~ LIQUID_DF_15$Stage_grouped, center = median)#0.801
car::leveneTest(LIQUID_DF_15$DLL1_URINE ~ LIQUID_DF_15$Stage_grouped, center = median) #2.2e-16 - not normal 
car::leveneTest(LIQUID_DF_15$HES1_URINE ~ LIQUID_DF_15$Stage_grouped, center = median)#0.2948
car::leveneTest(LIQUID_DF_15$CTNNB1_URINE ~ LIQUID_DF_15$Stage_grouped, center = median)#0.6521

#GRADE urine normalcy
by(LIQUID_DF_15$NOTCH2_URINE, LIQUID_DF_15$Grade_simple, shapiro.test)#not enough data
by(LIQUID_DF_15$DLL1_URINE, LIQUID_DF_15$Grade_simple, shapiro.test) #Not enough data
by(LIQUID_DF_15$HES1_URINE, LIQUID_DF_15$Grade_simple, shapiro.test)#Not enough data
by(LIQUID_DF_15$CTNNB1_URINE, LIQUID_DF_15$Grade_simple, shapiro.test)#Not enough data
#urine variance
car::leveneTest(LIQUID_DF_15$NOTCH2_URINE ~ LIQUID_DF_15$Grade_simple, center = median)#0.422
car::leveneTest(LIQUID_DF_15$DLL1_URINE ~ LIQUID_DF_15$Grade_simple, center = median) #not enough data
car::leveneTest(LIQUID_DF_15$HES1_URINE ~ LIQUID_DF_15$Grade_simple, center = median)#0.6458
car::leveneTest(LIQUID_DF_15$CTNNB1_URINE ~ LIQUID_DF_15$Grade_simple, center = median)#0.5364

#using wilcox.test because the data is very small and cant determine some of the normalcies / variances
wilcox.test(NOTCH2_URINE ~ Stage_grouped, data = LIQUID_DF_15)#1
wilcox.test(DLL1_URINE ~ Stage_grouped, data = LIQUID_DF_15)#1
wilcox.test(HES1_URINE ~ Stage_grouped, data = LIQUID_DF_15)#0.9495
wilcox.test(CTNNB1_URINE ~ Stage_grouped, data = LIQUID_DF_15)#0.1986
#using wilcox.test because the data is very small and cant determine some of the normalcies / variances
wilcox.test(NOTCH2_URINE ~ Grade_simple, data = LIQUID_DF_15)#1
wilcox.test(DLL1_URINE ~ Grade_simple, data = LIQUID_DF_15)#not enough data
wilcox.test(HES1_URINE ~ Grade_simple, data = LIQUID_DF_15)#0.5333
wilcox.test(CTNNB1_URINE ~ Grade_simple, data = LIQUID_DF_15)#0.8571
