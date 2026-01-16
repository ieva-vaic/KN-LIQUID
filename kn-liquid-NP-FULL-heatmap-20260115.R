#KN-liquid - np, urine, plasma #2026 01 15
#HEATMAP NP vs TISSUE vs plasma vs urine, ALL DATA
Sys.setenv(LANG = "en")
#libraries
library(ComplexHeatmap)
library(circlize)
library(tidyverse)
#read RDS
LIQUID_DF <- readRDS("C:/Users/Ieva/rprojects/OTHER DATA/KN_LIQUID/liquid_20260114.RDS")

raiska_np <- c("NOTCH2_NP",      "CTNNB1_NP" ,     "DLL1_NP" ,       "HES1_NP")
raiska_TISSUE <- c("NOTCH2_TUMOR", "CTNNB1_TUMOR", "DLL1_TUMOR", "HES1_TUMOR")
raiska_urine <- c("NOTCH2_URINE", "CTNNB1_URINE", "DLL1_URINE", "HES1_URINE")
raiska_p <- c("NOCTH2_P" , "CTNNB1_P" ,"DLL1_P" ,"HES1_P")

# heatmap expression continuous data as matrix###########
Heat_data <- LIQUID_DF[, c("Laboratorinis kodas",
                              raiska_urine, raiska_np,raiska_TISSUE     
)]
Heat_data <- as.data.frame(Heat_data)
rownames(Heat_data) <- Heat_data[, 1]
Heat_data <- Heat_data[, -1]
heatmap_raiska <- Heatmap(as.matrix(Heat_data), cluster_rows = F)
heatmap_raiska #show

# clinical data as df####################################
clinical <- LIQUID_DF[, c("Laboratorinis kodas", "TYPE",  "CA125" ,"Stage_simple",
                             "Grade", "Age", "STATUS")]
clinical <- as.data.frame(clinical)
rownames(clinical) <- clinical$KN
head(clinical)
#fix satus na so it shows on the legend
clinical$STATUS[is.na(clinical$STATUS)] <- "NA"
clinical$STATUS <- factor(clinical$STATUS)
clinical$STATUS <- dplyr::recode(clinical$STATUS, "1" = "Gyva", "2" = "Mirusi")
clinical$CA125_f <- ifelse(
  clinical$CA125 > 35,
  "CA125 increase",
  "Norm"
)
clinical$CA125_f[is.na(clinical$CA125_f)] <- "NA"
clinical$CA125_f <- dplyr::recode(clinical$CA125_f, `CA125 increase`= "CA125 padidėjimas", Norm = "Norma")
clinical$Grade <- replace(clinical$Grade, is.na(clinical$Grade), "NA" )
clinical$Grade <- as.character(clinical$Grade)  # ensure it's character
clinical$Grade <- case_when(
  clinical$Grade %in% c("G2", "G1&G1", "G2&G1") ~ "G2&G1",
  clinical$Grade %in% c("GL", "GB") ~ NA_character_,
  TRUE ~ clinical$Grade  # keep G1, G3, NA as is
)
clinical$Stage_simple[is.na(clinical$Stage_simple)] <- "NA"
clinical$Stage_simple <- factor(clinical$Stage_simple)
clinical$TYPE <- dplyr::recode(clinical$TYPE, BENIGN= "Gerybiniai", OTHER = "Kiti KV",
                               `ENDOMETRIAL CANCER` = "Gimdos kūno vėžys",
                               RSS ="Riziką mažinanti operacija")

# clinical data annotation
col_age <- colorRamp2(c(40, 90), c("#9cd4c4", "#3c402f"))  # age colors
row_ha = rowAnnotation(
                       Navikas = clinical$TYPE, Stadija = clinical$Stage_simple,
                       `Diferenciacijos laipsnis` = clinical$Grade,
                       Amžius = clinical$Age, CA125 = clinical$CA125_f,
                       `Mirties statusas` =clinical$STATUS, 
                       col = list(Navikas = c("HGSOC" = "#a89cd4",
                                              "Gerybiniai" = "#d49cac", 
                                              "Kiti KV" = "darkblue",
                                              "Gimdos kūno vėžys" = "purple",
                                              "Riziką mažinanti operacija" = "green"
                                              ), 
                                  Stadija = c("1" = "#9cd4c4",  
                                              "2" = "#c8d49c", 
                                              "3" = "#d49cac", 
                                              "4" = "#a89cd4", 
                                              "NA" = "grey"), 
                                  `Diferenciacijos laipsnis` = c("G2&G1"="#9cd4c4",
                                                                 "G1" = "green",
                                                                 "G3" = "#a89cd4", 
                                                                 "NA" = "grey"),
                                  CA125 = c("Norma"="#9cd4c4", 
                                            "CA125 padidėjimas" = "#3c402f",
                                            "NA" = "grey"),
                                  `Mirties statusas` = c("Gyva"="#9cd4c4", 
                                                         "Mirusi" = "#3c402f",
                                                         "NA" = "grey"),
                                  Amžius = col_age
                       ))

#expression colors
col_fun = colorRamp2(c(2, -5, -10, -15), c("#8564fb",  "#64b3fb","#e088bd", "#af2745"))
#col_fun = colorRamp2(c(2, 0, -2, -4, -6, -8, -10, -12, -14, -16),
#                     c("#e7e0fe", "#cec1fd", "#8564fb", "#64b3fb","#93cafc","#325a7e", "#e088bd", "#af2745", "#9e233e", "#4f121f"))
labels <- c("\u221215", "\u221210", "\u22125", "2" ) #THIS IS NEEDED FOR MDPI AT LEAST - LONG MINUS SIGNS
#make rows order
rows_order <-  c( "NOTCH2_URINE",  "CTNNB1_URINE",  "DLL1_URINE",    "HES1_URINE",
                  "NOTCH2_NP",     "CTNNB1_NP"  ,  "DLL1_NP" ,      "HES1_NP"   ,
                  "NOTCH2_TUMOR", "CTNNB1_TUMOR", "DLL1_TUMOR", "HES1_TUMOR" )
#############################
# create groups from suffix
col_groups <- ifelse(grepl("_URINE$", rows_order), "URINE",
                     ifelse(grepl("_NP$", rows_order), "NP", "TUMOR"))
# turn into factor with desired order
col_groups <- factor(col_groups, levels = c("URINE", "NP", "TUMOR"))
# check
col_groups

#final heatmap
heatmap_raiska <- Heatmap(as.matrix(Heat_data), 
                          cluster_rows = F,
                          cluster_columns = F,
                          name = "Santykinė genų raiška",  
                          right_annotation = row_ha,
                          col = col_fun, 
                          row_split = clinical$TYPE, 
                          column_names_gp = gpar(fontface = "italic"),
                          column_title = "Santykinė genų raiška", 
                          row_names_gp = gpar(fontsize = 8), 
                          column_order = rows_order,
                          column_split = col_groups,
                          heatmap_legend_param = list( #THIS IS FOR THE LONG MINUS SIGNS
                            at = c(2, -5, -10, -15),   # Legend positions
                            labels = labels)     # Adjusted labels
)
heatmap_raiska

#add plasma data ###################################
plasma_only <- LIQUID_DF[, c("Laboratorinis kodas", raiska_p)]
plasma_only$NOCTH2_P <- factor(plasma_only$NOCTH2_P, 
                               levels = c(levels(plasma_only$NOCTH2_P), "NA"))
plasma_only$HES1_P <- factor(plasma_only$HES1_P, 
                               levels = c(levels(plasma_only$HES1_P), "NA"))
plasma_only$CTNNB1_P <- factor(plasma_only$CTNNB1_P, 
                               levels = c(levels(plasma_only$CTNNB1_P), "NA"))
plasma_only$DLL1_P <- factor(plasma_only$DLL1_P, 
                               levels = c(levels(plasma_only$DLL1_P), "NA"))
plasma_only$NOCTH2_P[is.na(plasma_only$NOCTH2_P)] <- "NA"
plasma_only$CTNNB1_P[is.na(plasma_only$CTNNB1_P)] <- "NA"
plasma_only$DLL1_P[is.na(plasma_only$DLL1_P)] <- "NA"
plasma_only$HES1_P[is.na(plasma_only$HES1_P)] <- "NA"

rownames(plasma_only) <- plasma_only$`Laboratorinis kodas`
plasma_only <- plasma_only[, -1]

plasma_only$CTNNB1_P <- dplyr::recode(
  plasma_only$CTNNB1_P,
  "0" = "nenustatyta raiška",
  "1" = "raiška"
)
plasma_only$HES1_P <- dplyr::recode(
  plasma_only$HES1_P,
  "0" = "nenustatyta raiška",
  "1" = "raiška"
)

plasma_only$DLL1_P <- dplyr::recode(
  plasma_only$DLL1_P,
  "0" = "nenustatyta raiška",
  "1" = "raiška"
)
plasma_only$NOCTH2_P <- dplyr::recode(
  plasma_only$NOCTH2_P,
  "0" = "nenustatyta raiška",
  "1" = "raiška"
)
col_fun <- c("raiška" = "magenta", "nenustatyta raiška" = "darkblue", "NA" ="grey")
heatmap_plasma <- Heatmap(as.matrix(plasma_only), 
                          col = col_fun, name = "Raiška plazmoje",
                          cluster_rows = FALSE,
                          row_split = clinical$TYPE
)
heatmap_plasma

#add together with the main plot
main_plot <- heatmap_plasma + heatmap_raiska
main_plot

png("C:/Users/Ieva/rprojects/outputs_all/all_liquid_heatmap202060115.png",
    width = 2500, height = 2800, res = 200)
main_plot
dev.off()
