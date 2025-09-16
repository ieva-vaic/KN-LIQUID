#KN-liquid - np
#data - 20250916, changed KN-97 from the repeated expreriment
#add survival

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
#libraries for heatmap
library(ComplexHeatmap)
library(circlize)
#open tissue and liquid np data ################################################
KN_data_np_tissue <- readRDS("C:/Users/Ieva/rprojects/OTHER DATA/KN_data_np_tissue_full20250916.RDS")
#open survival data###########################################################
SURVIVAL <- readxl::read_xlsx("C:/Users/Ieva/rprojects/OTHER DATA/KN-DISSERTATION FILES/KN_MIRTIES_FAILAS_20250911.xlsx")
#make smaller surv df
SURVIVAL <- SURVIVAL[, colnames(SURVIVAL) %in% c("patient_id_aud", "OS", "STATUS")]
#join
SURVIVAL$patient_id_aud %in% KN_data_np_tissue$patient_id_aud
KN_DATA_NP <- left_join(SURVIVAL, KN_data_np_tissue, by = "patient_id_aud")
#add urine####################################
urine_liquid <- readxl::read_xlsx("C:/Users/Ieva/rprojects/OTHER DATA/KN_LIQUID/KN-liquid-splapimas20250916.xlsx")
#fix patient names
urine_liquid$patient_id_aud <- gsub("-S", "", urine_liquid$`KN nr.`)
#rename urine columns
urine_liquid <- urine_liquid %>%
  rename("NOTCH2_URINE" = "NOTCH2_DELTA", 
         "CTNNB1_URINE" = "CTNNB1_DELTA",
         "DLL1_URINE" = "DLL1_DELTA",
         "HES1_URINE" = "HES1_DELTA")
urine_np <- c("NOTCH2_URINE", "CTNNB1_URINE", "DLL1_URINE", "HES1_URINE")
#make small urine df
urine_df <- urine_liquid[, colnames(urine_liquid) %in% c(urine_np, "patient_id_aud") ]
#join
urine_df$patient_id_aud %in% KN_DATA_NP$patient_id_aud
KN_DATA_NP_UR <- left_join(KN_DATA_NP, urine_df, by = "patient_id_aud")
#add plasma data ##################################################################
plasma_liquid <- readxl::read_xlsx("C:/Users/Ieva/rprojects/OTHER DATA/KN_LIQUID/KN_liquid_plasma_20250813.xlsx")
#code 0 and 1
plasma_liquid <- plasma_liquid %>%
  mutate(across(14:17, ~ ifelse(. == 0, "nenustatyta raiška", "raiška")))
plasma_liquid[, 14:17]
#fix patient names
plasma_liquid$patient_id_aud <- gsub("-P", "", plasma_liquid$`KN nr.`)
#make smaller df
plasma_liquid <- plasma_liquid[, c(14:18)]
#join
KN_DATA_NP_UR$patient_id_aud %in% plasma_liquid$patient_id_aud
KN_data_full <- left_join(KN_DATA_NP_UR,plasma_liquid, by = "patient_id_aud" )
#save 
#saveRDS(KN_data_full, "C:/Users/Ieva/rprojects/OTHER DATA/KN_LIQUID/KN_FULL_LIQUID_R250916.RDS")

#PLASMA fisher tests####################################################################
#hgsoc vs others PLASMA FISHER TEST 
OC_HGSOC_OTHER<- KN_data_full[c(KN_data_full$Grupė_Ieva != "Benign"),] 
OC_HGSOC_OTHER$Grupė_Ieva <- relevel(factor(OC_HGSOC_OTHER$Grupė_Ieva), ref = "Other")

fisher.test(table(OC_HGSOC_OTHER$Grupė_Ieva, OC_HGSOC_OTHER$HES1_P))#1
fisher.test(table(OC_HGSOC_OTHER$Grupė_Ieva, OC_HGSOC_OTHER$CTNNB1_P)) #0.2418
fisher.test(table(OC_HGSOC_OTHER$Grupė_Ieva, OC_HGSOC_OTHER$NOCTH2_P)) ##0.2418
fisher.test(table(OC_HGSOC_OTHER$Grupė_Ieva, OC_HGSOC_OTHER$DLL1_P)) #data too small

## HES1#################################################
# Save the plot as a PNG file
png("C:/Users/Ieva/rprojects/outputs_all/barplothes1plasma250916.png",
    width = 1300, height = 1100, res = 200)
# Create contingency table
tab <- table(OC_HGSOC_OTHER$Grupė_Ieva, OC_HGSOC_OTHER$HES1_P)

# Fisher's exact test
ft <- fisher.test(tab)
pval <- signif(ft$p.value, 3)  # rounded p-value

# Convert to percentages by column
tab_perc <- prop.table(tab, margin = 2) * 100

# Stacked bar plot with italic HES1 in title
bp <- barplot(tab_perc, col = c("skyblue", "salmon"),
              legend = TRUE,
              main = expression("Plasma detection of " * italic(HES1)),
              sub = paste0("Fisher's exact test p = ", pval),
              ylab = "Percentage",
              beside = FALSE)

# Add percentages on bars
for(i in 1:nrow(tab_perc)){
  text(x = bp, y = apply(tab_perc[1:i, , drop = FALSE], 2, sum) - tab_perc[i, ]/2,
       labels = paste0(round(tab_perc[i, ], 1), "%"), cex = 0.8, col = "black")
}
dev.off()


## CTNNB1################################
# Save the plot as a PNG file
png("C:/Users/Ieva/rprojects/outputs_all/barplotCTNNB1plasma250916.png",
    width = 1300, height = 1100, res = 200)
# Create contingency table
tab <- table(OC_HGSOC_OTHER$Grupė_Ieva, OC_HGSOC_OTHER$CTNNB1_P)

# Fisher's exact test
ft <- fisher.test(tab)
pval <- signif(ft$p.value, 3)  # rounded p-value

# Convert to percentages by column
tab_perc <- prop.table(tab, margin = 2) * 100

# Stacked bar plot with italic HES1 in title
bp <- barplot(tab_perc, col = c("skyblue", "salmon"),
              legend = TRUE,
              main = expression("Plasma detection of " * italic(CTNNB1)),
              sub = paste0("Fisher's exact test p = ", pval),
              ylab = "Percentage",
              beside = FALSE)

# Add percentages on bars
for(i in 1:nrow(tab_perc)){
  text(x = bp, y = apply(tab_perc[1:i, , drop = FALSE], 2, sum) - tab_perc[i, ]/2,
       labels = paste0(round(tab_perc[i, ], 1), "%"), cex = 0.8, col = "black")
}
dev.off()

## NOTCH2#############################################################
# Save the plot as a PNG file
png("C:/Users/Ieva/rprojects/outputs_all/barplotnotch2plasma250916.png",
    width = 1300, height = 1100, res = 200)
# Create contingency table
tab <- table(OC_HGSOC_OTHER$Grupė_Ieva, OC_HGSOC_OTHER$NOCTH2_P)

# Fisher's exact test
ft <- fisher.test(tab)
pval <- signif(ft$p.value, 3)  # rounded p-value

# Convert to percentages by column
tab_perc <- prop.table(tab, margin = 2) * 100

# Stacked bar plot with italic HES1 in title
bp <- barplot(tab_perc, col = c("skyblue", "salmon"),
              legend = TRUE,
              main = expression("Plasma detection of " * italic(NOTCH2)),
              sub = paste0("Fisher's exact test p = ", pval),
              ylab = "Percentage",
              beside = FALSE)

# Add percentages on bars
for(i in 1:nrow(tab_perc)){
  text(x = bp, y = apply(tab_perc[1:i, , drop = FALSE], 2, sum) - tab_perc[i, ]/2,
       labels = paste0(round(tab_perc[i, ], 1), "%"), cex = 0.8, col = "black")
}
dev.off()

