#KN-  liquid 2026 01 12
#CLEAN data, after additional np cases
Sys.setenv(LANG = "en")
#libraries
library(tidyverse)
library(readxl)
library(openxlsx)
#upload each data sheet separate#######################
KN_clinical <- read_xlsx("C:/Users/Ieva/rprojects/OTHER DATA/KN_LIQUID/KN-LIQUID-PROJECT REZULTATAI 20260114.xlsx", sheet = 'LIQUID KLINIKA')
KN_np_full <- read_xlsx("C:/Users/Ieva/rprojects/OTHER DATA/KN_LIQUID/KN-LIQUID-PROJECT REZULTATAI 20260114.xlsx", sheet = 'LIQUID NP FULL')
KN_urine <- read_xlsx("C:/Users/Ieva/rprojects/OTHER DATA/KN_LIQUID/KN-LIQUID-PROJECT REZULTATAI 20260114.xlsx", sheet = 'LIQUID URINE')
KN_plasma <- read_xlsx("C:/Users/Ieva/rprojects/OTHER DATA/KN_LIQUID/KN-LIQUID-PROJECT REZULTATAI 20260114.xlsx", sheet = 'LIQUID PLASMA')
KN_tissue <- read_xlsx("C:/Users/Ieva/rprojects/OTHER DATA/KN_LIQUID/KN-LIQUID-PROJECT REZULTATAI 20260114.xlsx", sheet = 'LIQUID TUMOR')

#merge all liquid project data#######################
#sample each data
KN_clinical$`Laboratorinis kodas`
KN_np_full$`KN nr.`
KN_urine$`KN nr.`
KN_plasma$`KN nr.`
KN_tissue$patient_id_aud
#create a new "laboratorinis kodas" in each datasheet
KN_np_full$`Laboratorinis kodas` <- sub("-np$", "", KN_np_full$`KN nr.`)
KN_urine$`Laboratorinis kodas` <- sub("-S$", "", KN_urine$`KN nr.`)
KN_plasma$`Laboratorinis kodas` <- sub("-P$", "", KN_plasma$`KN nr.`)
KN_tissue$`Laboratorinis kodas` <- KN_tissue$patient_id_aud
#chek if "laboratorinis kodas" is equal between sheets
KN_np_full$`Laboratorinis kodas`  %in% KN_clinical$`Laboratorinis kodas`
KN_urine$`Laboratorinis kodas`  %in% KN_clinical$`Laboratorinis kodas`
KN_plasma$`Laboratorinis kodas`  %in% KN_clinical$`Laboratorinis kodas`
KN_tissue$`Laboratorinis kodas`  %in% KN_clinical$`Laboratorinis kodas`
#all in there
#joining##################################
#join 1
merged_df1 <- left_join(KN_clinical, KN_np_full, by = "Laboratorinis kodas")
#join 2
merged_df2 <- left_join(merged_df1, KN_urine, by = "Laboratorinis kodas")
#join 3
merged_df3 <- left_join(merged_df2, KN_plasma, by = "Laboratorinis kodas")
#join 4
merged_df4 <- left_join(merged_df3, KN_tissue, by = "Laboratorinis kodas")
#cleaning################################
colnames(merged_df4)
#leave these:
leave <- c("Laboratorinis kodas", "BRCA_mut_klinikinė_kodas",
           "Amžius diagnozės metu",
           "CA125" ,"Ca 125 po gydymo",
           "STAGE", "Grade","L/M", "MTS",
           "sutrumpinta_diagnozė", "Grupė_Ieva.x", "TYPE",
           "OS", "STATUS" ,
           "NOTCH2_NP", "CTNNB1_NP", "DLL1_NP","HES1_NP",
           "NOTCH2_URINE","CTNNB1_URINE" ,"DLL1_URINE", "HES1_URINE",  
           "NOCTH2_P","CTNNB1_P","DLL1_P", "HES1_P",
           "NOTCH2_TUMOR","CTNNB1_TUMOR" , "DLL1_TUMOR", "HES1_TUMOR",
           "NOTCH2_P_norm", "CTNNB1_P_norm", "DLL1_P_norm", "HES1_P_norm" )
LIQUID_DF <- merged_df4[, colnames(merged_df4)%in% leave]
#rename and fix clinicals###############
#rename columns
LIQUID_DF <- LIQUID_DF %>% rename(
  Age = `Amžius diagnozės metu`,
  BRCA_CLINICAL_MUT_STATUS = `BRCA_mut_klinikinė_kodas`,
  CA125_post_op = `Ca 125 po gydymo`,
  Stage = STAGE,
  Diagnosis = sutrumpinta_diagnozė,
  TYPE_long = Grupė_Ieva.x
)
colnames(LIQUID_DF)
#fix brca
LIQUID_DF$BRCA_CLINICAL_MUT_STATUS <- recode(
  as.character(LIQUID_DF$BRCA_CLINICAL_MUT_STATUS),
  "3" = "No data",
  "1" = "Norm",
  "2" = "Mutation"
)
LIQUID_DF$BRCA_CLINICAL_MUT_STATUS <- as.factor(LIQUID_DF$BRCA_CLINICAL_MUT_STATUS)
#remove stage information form benign / RSS
LIQUID_DF <- LIQUID_DF %>%
  mutate(
    across(
      c(Stage, Grade, `L/M`, MTS),
      ~ if_else(TYPE %in% c("RSS", "BENIGN"), NA, .)
    )
  )
#make factors
LIQUID_DF$TYPE <- as.factor(LIQUID_DF$TYPE)
LIQUID_DF$NOCTH2_P <- as.factor(LIQUID_DF$NOCTH2_P)
LIQUID_DF$CTNNB1_P <- as.factor(LIQUID_DF$CTNNB1_P)
LIQUID_DF$DLL1_P <- as.factor(LIQUID_DF$DLL1_P)
LIQUID_DF$HES1_P <- as.factor(LIQUID_DF$HES1_P)
#fix stage
LIQUID_DF <- LIQUID_DF %>%
  mutate(Stage_simple = case_when(
    grepl("^IV", Stage) ~ 4,       # check IV first
    grepl("^III", Stage) ~ 3,      # then III
    grepl("^II", Stage) ~ 2,       # then II
    grepl("^I", Stage) ~ 1,        # finally I
    TRUE ~ NA_real_                 # anything else / NA
  ))

#remove cases that is NA in their type: 
LIQUID_DF <- LIQUID_DF %>%
  filter(!`Laboratorinis kodas` %in% c("KN-098", "KN-102"))

#save xcel
write.xlsx(LIQUID_DF, "C:/Users/Ieva/rprojects/OTHER DATA/KN_LIQUID/liquid_20260114.xlsx")
#save RDS
saveRDS(LIQUID_DF, "C:/Users/Ieva/rprojects/OTHER DATA/KN_LIQUID/liquid_20260114.RDS")

