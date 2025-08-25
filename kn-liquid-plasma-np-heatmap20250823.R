#HEATMAP NP vs TISSUE
Sys.setenv(LANG = "en")
#libraries
library(ComplexHeatmap)
library(circlize)
library(tidyverse)

KN_data_full <- readRDS("C:/Users/Ieva/rprojects/OTHER DATA/KN_data_np_tissue_full20250813.RDS")
#rename np columns
KN_data_full <- KN_data_full %>%
  rename("NOTCH2_NP" = "NOTCH2_DELTA", 
         "CTNNB1_NP" = "CTNNB1_DELTA",
         "DLL1_NP" = "DLL1_DELTA",
         "HES1_NP" = "HES1_DELTA")
raiska_np <- c("NOTCH2_NP",      "CTNNB1_NP" ,     "DLL1_NP" ,       "HES1_NP")
raiska_tissue <- c("NOTCH2", "CTNNB1", "DLL1", "HES1")

# heatmap data as matrix
Heat_data <- KN_data_full[, c("patient_id_aud",
                         raiska_np,    raiska_tissue     
)]
Heat_data <- as.data.frame(Heat_data)
rownames(Heat_data) <- Heat_data[, 1]
Heat_data <- Heat_data[, -1]
heatmap_raiska <- Heatmap(as.matrix(Heat_data))

# clinical data as df
clinical <- KN_data_full[, c("patient_id_aud", "Grupė_Ieva",  "CA125" ,"Stage4", "Grade2", "Amžius", "CA125_f", "Histology")]
clinical <- as.data.frame(clinical)
rownames(clinical) <- clinical$KN
head(clinical)
#fix na so it shows on the legend
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

#final heatmap
heatmap_raiska <- Heatmap(as.matrix(Heat_data), 
                          cluster_rows = F,
                          name = "Santykinė genų raiška",  
                          right_annotation = row_ha,
                          col = col_fun, 
                          row_split = clinical$type, 
                          column_names_gp = gpar(fontface = "italic"),
                          column_title = "Santykinė genų raiška", 
                          row_names_gp = gpar(fontsize = 8), 
                          #row_order = rows_order,
                          heatmap_legend_param = list( #THIS IS FOR THE LONG MINUS SIGNS
                            at = c(2, -5, -10, -15),   # Legend positions
                            labels = labels)     # Adjusted labels
)
heatmap_raiska

#save png
png("C:/Users/Ieva/rprojects/outputs_all/heatmap_np_tissue_20250813.png", width = 3000, height = 2500, res = 300) # width and height in pixels, resolution in dpi
heatmap_raiska
dev.off() # Close the PNG device

#add plasma data ###################################
plasma_liquid <- readxl::read_xlsx("C:/Users/Ieva/rprojects/OTHER DATA/KN_liquid_plasma_20250813.xlsx")
#code 0 and 1
plasma_liquid <- plasma_liquid %>%
  mutate(across(14:17, ~ ifelse(. == 0, "nenustatyta raiška", "raiška")))
plasma_liquid[, 14:17]
#fix patient names
plasma_liquid$patient_id_aud <- gsub("-P", "", plasma_liquid$`KN nr.`)
#add plasma data to the main df
KN_liquid_np_p <- merge(plasma_liquid, KN_data_full, by = "patient_id_aud", all.y = TRUE)
colnames(KN_liquid_np_p)
#make smaller merged df
colnames(KN_liquid_np_p)
KN_liquid_np_p <- KN_liquid_np_p[, c(1, 2,15:18,20:52)]
plasma_only <- KN_liquid_np_p[, c(1, 3:6)]
rownames(plasma_only) <- plasma_only[, 1]
plasma_only <- plasma_only[, -1]
col_fun <- c("raiška" = "magenta", "nenustatyta raiška" = "darkblue")
heatmap_plasma <- Heatmap(as.matrix(plasma_only), 
                          col = col_fun, name = "Raiška plazmoje",
                          cluster_rows = FALSE, row_split = clinical$type)
heatmap_plasma

#add together with the main plot
main_plot <- heatmap_plasma + heatmap_raiska
main_plot

#save png
png("C:/Users/Ieva/rprojects/outputs_all/heatmap_np_p_tissue_20250813.png", width = 3000, height = 2500, res = 300) # width and height in pixels, resolution in dpi
main_plot
dev.off() # Close the PNG device


#hgsoc vs others PLASMA FISHER TEST ##################################
OC_HGSOC_OTHER<- KN_liquid_np_p[c(KN_liquid_np_p$Grupė_Ieva != "Benign"),] 
OC_HGSOC_OTHER$Grupė_Ieva <- relevel(factor(OC_HGSOC_OTHER$Grupė_Ieva), ref = "Other")

fisher.test(table(OC_HGSOC_OTHER$Grupė_Ieva, OC_HGSOC_OTHER$HES1_P))
fisher.test(table(OC_HGSOC_OTHER$Grupė_Ieva, OC_HGSOC_OTHER$CTNNB1_P))
fisher.test(table(OC_HGSOC_OTHER$Grupė_Ieva, OC_HGSOC_OTHER$NOCTH2_P))
fisher.test(table(OC_HGSOC_OTHER$Grupė_Ieva, OC_HGSOC_OTHER$DLL1_P))

#add urine data ###################################
urine_liquid <- readxl::read_xlsx("C:/Users/Ieva/rprojects/OTHER DATA/KN-liquid-splapimas20250822.xlsx")
#fix patient names
urine_liquid$patient_id_aud <- gsub("-S", "", urine_liquid$`KN nr.`)
#rename urine columns
urine_liquid <- urine_liquid %>%
  rename("NOTCH2_URINE" = "NOTCH2_DELTA", 
         "CTNNB1_URINE" = "CTNNB1_DELTA",
         "DLL1_URINE" = "DLL1_DELTA",
         "HES1_URINE" = "HES1_DELTA")
urine_np <- c("NOTCH2_URINE",      "CTNNB1_URINE" ,     "DLL1_URINE" ,       "HES1_URINE")


#add urine data to the main df
KN_liquid_np_p_urine <- merge(urine_liquid, KN_liquid_np_p, by = "patient_id_aud", all.x = TRUE, all.y = TRUE)
colnames(KN_liquid_np_p_urine)
#remove uneeded columns form urine side of df
KN_liquid_np_p_urine <- KN_liquid_np_p_urine[, -c(40:68)]

# heatmap data as matrix
Heat_data_u <- KN_liquid_np_p_urine[, c("patient_id_aud",
                              urine_np
)]
Heat_data_u <- as.data.frame(Heat_data_u)
rownames(Heat_data_u) <- Heat_data_u[, 1]
Heat_data_u <- Heat_data_u[, -1]
heatmap_urine <- Heatmap(as.matrix(Heat_data_u), name = "Raiška šlapime")

main_plot2 <- heatmap_plasma +heatmap_urine + heatmap_raiska 

#save png
png("C:/Users/Ieva/rprojects/outputs_all/heatmap_np_p_urine_tissue_20250822.png", width = 3000, height = 2500, res = 300) # width and height in pixels, resolution in dpi
main_plot2
dev.off() # Close the PNG device

