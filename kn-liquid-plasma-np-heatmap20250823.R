#KN-liquid - np, urine, plasma #20250917
#HEATMAP NP vs TISSUE vs plasma vs urine
Sys.setenv(LANG = "en")
#libraries
library(ComplexHeatmap)
library(circlize)
library(tidyverse)
KN_data_full <- readRDS("C:/Users/Ieva/rprojects/OTHER DATA/KN_LIQUID/KN_FULL_LIQUID_R250916.RDS")
# fix colnames
colnames(KN_data_full)
#rename np columns
KN_data_full <- KN_data_full %>%
  rename("NOTCH2_NP" = "NOTCH2_DELTA", 
         "CTNNB1_NP" = "CTNNB1_DELTA",
         "DLL1_NP" = "DLL1_DELTA",
         "HES1_NP" = "HES1_DELTA")
#rename tissuy columns
KN_data_full <- KN_data_full %>%
  rename("NOTCH2_TISSUE" = "NOTCH2", 
         "CTNNB1_TISSUE" = "CTNNB1",
         "DLL1_TISSUE" = "DLL1",
         "HES1_TISSUE" = "HES1")
raiska_np <- c("NOTCH2_NP",      "CTNNB1_NP" ,     "DLL1_NP" ,       "HES1_NP")
raiska_TISSUE <- c("NOTCH2_TISSUE", "CTNNB1_TISSUE", "DLL1_TISSUE", "HES1_TISSUE")
raiska_urine <- c("NOTCH2_URINE", "CTNNB1_URINE", "DLL1_URINE", "HES1_URINE")
raiska_p <- c("NOCTH2_P" , "CTNNB1_P" ,"DLL1_P" ,"HES1_P")
#chek
tail(KN_data_full[, colnames(KN_data_full) %in% raiska_p])

# heatmap expression continuous data as matrix###########
Heat_data <- KN_data_full[, c("patient_id_aud",
                              raiska_urine, raiska_np,raiska_TISSUE     
)]
Heat_data <- as.data.frame(Heat_data)
rownames(Heat_data) <- Heat_data[, 1]
Heat_data <- Heat_data[, -1]
heatmap_raiska <- Heatmap(as.matrix(Heat_data))
heatmap_raiska #show

# clinical data as df
clinical <- KN_data_full[, c("patient_id_aud", "Grupė_Ieva",  "CA125" ,"Stage4",
                             "Grade2", "Amžius", "CA125_f", "Histology", "STATUS")]
clinical <- as.data.frame(clinical)
rownames(clinical) <- clinical$KN
head(clinical)
#fix na so it shows on the legend
clinical$STATUS[is.na(clinical$STATUS)] <- "NA"
clinical$STATUS <- factor(clinical$STATUS)
clinical$STATUS <- recode(clinical$STATUS, "0" = "Gyva", "1" = "Mirusi")

clinical$CA125_f <- factor(clinical$CA125_f, levels = c("CA125 increase", "Norm", "NA"), exclude = NULL)
clinical$CA125_f <- recode(clinical$CA125_f, `CA125 increase`= "CA125 padidėjimas", Norm = "Norma")
clinical$Stage4 <- replace(clinical$Stage4, is.na(clinical$Stage4), "NA" )
clinical$Grade2 <-replace(clinical$Grade2, is.na(clinical$Grade2), "NA" )
clinical$type <- factor(clinical$Grupė_Ieva)
clinical$type <- recode(clinical$type, `Benign`= "Gerybiniai", Other = "Kiti KV")
#fix histoloy names
# clinical data as df
clinical$Histology <- ifelse(clinical$Histology %in% c("Endometrial", "Endometriod"), "Endometrioidinis",
                             clinical$Histology )
clinical$Histology <- ifelse(clinical$Histology == "Endometriois", "Endometriozė",
                             clinical$Histology )
clinical$Histology <- ifelse(clinical$Histology == "Cystis", "Cista",
                             clinical$Histology )
clinical$Histology <- ifelse(clinical$Histology == "Mucinous", "Mucininis",
                             clinical$Histology )
clinical$Histology <- ifelse(clinical$Histology == "Serous", "Serozinis",
                             clinical$Histology )
clinical$Histology <- ifelse(clinical$Histology == "Clear cell", "Šviesių lastelių",
                             clinical$Histology )
clinical$Histology <- ifelse(clinical$Histology == "RSS", "Riziką mažinanti operacija",
                             clinical$Histology )
clinical$Histology <- ifelse(clinical$Histology == "Granuloza", "Granulosa",
                             clinical$Histology )
table(clinical$Histology)

# clinical data annotation
col_age <- colorRamp2(c(40, 90), c( "#9cd4c4", "#3c402f")) #age colors
row_ha = rowAnnotation(Histologija = clinical$Histology,
                       Navikas = clinical$type, Stadija = clinical$Stage4,
                       `Diferenciacijos laipsnis` = clinical$Grade2,
                       Amžius = clinical$Amžius, CA125 = clinical$CA125_f,
                       `Mirties statusas` =clinical$STATUS, 
                       col = list(Navikas = c("HGSOC" = "#a89cd4",
                                              "Gerybiniai" = "#d49cac", 
                                              "Kiti KV" = "darkblue"), 
                                  Stadija = c("I" = "#9cd4c4",  
                                              "II" = "#c8d49c", 
                                              "III" = "#d49cac", 
                                              "IV" = "#a89cd4", 
                                              "NA" = "grey"), 
                                  `Diferenciacijos laipsnis` = c("G1"="#9cd4c4", 
                                                                 "G3" = "#a89cd4", 
                                                                 "NA" = "grey"),
                                  CA125 = c("Norma"="#9cd4c4", 
                                            "CA125 padidėjimas" = "#3c402f",
                                            "NA" = "grey"),
                                  `Mirties statusas` = c("Gyva"="#9cd4c4", 
                                                         "Mirusi" = "#3c402f",
                                                         "NA" = "grey"),
                                  Amžius = col_age,
                                  Histologija = c("Šviesių lastelių" = "lightblue",
                                                  "Cista" = "lightgreen", 
                                                  "Endometrioidinis" = "green",
                                                  "Endometriozė" = "darkgreen", 
                                                  "Granulosa" = "turquoise",
                                                  "HGSOC" = "deeppink",
                                                  "Mioma" = "red", 
                                                  "Mucininis" = "yellow",
                                                  "Riziką mažinanti operacija" = "orange", 
                                                  "Serozinis" = "lightpink")
                       ))

#expression colors
col_fun = colorRamp2(c(2, -5, -10, -15), c("#8564fb",  "#64b3fb","#e088bd", "#af2745"))
#col_fun = colorRamp2(c(2, 0, -2, -4, -6, -8, -10, -12, -14, -16),
#                     c("#e7e0fe", "#cec1fd", "#8564fb", "#64b3fb","#93cafc","#325a7e", "#e088bd", "#af2745", "#9e233e", "#4f121f"))
labels <- c("\u221215", "\u221210", "\u22125", "2" ) #THIS IS NEEDED FOR MDPI AT LEAST - LONG MINUS SIGNS
#make rows order
rows_order <-  c( "NOTCH2_URINE",  "CTNNB1_URINE",  "DLL1_URINE",    "HES1_URINE",
                  "NOTCH2_NP",     "CTNNB1_NP"  ,  "DLL1_NP" ,      "HES1_NP"   ,
                  "NOTCH2_TISSUE" ,"CTNNB1_TISSUE" ,"DLL1_TISSUE" ,  "HES1_TISSUE"  )
#############################
# create groups from suffix
col_groups <- ifelse(grepl("_URINE$", rows_order), "URINE",
                     ifelse(grepl("_NP$", rows_order), "NP", "TISSUE"))
# turn into factor with desired order
col_groups <- factor(col_groups, levels = c("URINE", "NP", "TISSUE"))
# check
col_groups

#final heatmap
heatmap_raiska <- Heatmap(as.matrix(Heat_data), 
                          cluster_rows = F,
                          cluster_columns = F,
                          name = "Santykinė genų raiška",  
                          right_annotation = row_ha,
                          col = col_fun, 
                          row_split = clinical$type, 
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
plasma_only <- KN_data_full[, c("patient_id_aud", raiska_p)]
#make NA values
plasma_only$NOCTH2_P[is.na(plasma_only$NOCTH2_P)] <- "NA"
plasma_only$CTNNB1_P[is.na(plasma_only$CTNNB1_P)] <- "NA"
plasma_only$DLL1_P[is.na(plasma_only$DLL1_P)] <- "NA"
plasma_only$HES1_P[is.na(plasma_only$HES1_P)] <- "NA"

rownames(plasma_only) <- plasma_only$patient_id_aud
plasma_only <- plasma_only[, -1]
col_fun <- c("raiška" = "magenta", "nenustatyta raiška" = "darkblue", "NA" ="grey")
heatmap_plasma <- Heatmap(as.matrix(plasma_only), 
                          col = col_fun, name = "Raiška plazmoje",
                          cluster_rows = FALSE,
                          row_split = clinical$type
)
heatmap_plasma

#add together with the main plot
main_plot <- heatmap_plasma + heatmap_raiska
main_plot

png("C:/Users/Ieva/rprojects/outputs_all/all_liquid_heatmap20250917.png",
    width = 2500, height = 1500, res = 200)
main_plot
dev.off()
