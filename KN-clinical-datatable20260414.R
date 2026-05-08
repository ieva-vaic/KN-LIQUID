#KN-  liquid 2026 04-14
#DATA TABLE
Sys.setenv(LANG = "en")
library(tidyverse)
#read RDS
LIQUID_DF <- readRDS("C:/Users/Ieva/rprojects/OTHER DATA/KN_LIQUID/liquid_20260114.RDS")
#fix grade
LIQUID_DF <- LIQUID_DF %>%
  mutate(Grade_simple = case_when(
    grepl("G1\\&G1|G2\\&G1", Grade) ~ "G1",
    grepl("GB|GL", Grade) ~ NA_character_,
    TRUE ~ Grade
  ))
#fix CA125
LIQUID_DF <- LIQUID_DF %>%
  mutate(CA125_num = case_when(
    grepl("Neatlikta", CA125) ~ NA_character_,
    TRUE ~ CA125
  ))
LIQUID_DF$CA125_num <- as.numeric(LIQUID_DF$CA125_num)
#make tumor data df#################################
#only DLL1 tumor had empty variables, thus use any other tumor column
TUMOR_df <- LIQUID_DF %>%
  filter(!is.na(NOTCH2_TUMOR)) # 65 cases

#make lavage data df#####################################
#notch2 the fullest data
LAVAGE_df <- LIQUID_DF %>%
  filter(!is.na(NOTCH2_NP)) #103 cases

#make plasma and urine data df####################################
#notch2 the fullest data
PLASMA_URINE_df <- LIQUID_DF %>%
  filter(!is.na(LIQUID_DF$NOCTH2_P)) #15 cases

#HISTOLOGY ########################################################
#tissue
table(TUMOR_df$TYPE)
prop.table(table(TUMOR_df$TYPE)) * 100
#np
table(LAVAGE_df$TYPE)
prop.table(table(LAVAGE_df$TYPE)) * 100
#urine / plasma
table(PLASMA_URINE_df$TYPE)
prop.table(table(PLASMA_URINE_df$TYPE)) * 100

#STAGE###########################################
#all tissue
table(TUMOR_df$Stage_simple)
prop.table(table(TUMOR_df$Stage_simple)) * 100
#tissue by group stage
table(TUMOR_df$Stage_simple, TUMOR_df$TYPE)
#np
table(LAVAGE_df$Stage_simple, LAVAGE_df$TYPE)
#all np
table(LAVAGE_df$Stage_simple)
prop.table(table(LAVAGE_df$Stage_simple)) * 100
#all urine / plasma
table(PLASMA_URINE_df$Stage_simple)
prop.table(table(PLASMA_URINE_df$Stage_simple)) * 100
#plasma by type stage
table(PLASMA_URINE_df$Stage_simple, PLASMA_URINE_df$TYPE)

#GRADE##############################################
#tissue
table(TUMOR_df$Grade_simple)
prop.table(table(TUMOR_df$Grade_simple)) * 100
table(TUMOR_df$Grade_simple, TUMOR_df$TYPE, useNA = "a")
#URINE
table(PLASMA_URINE_df$Grade_simple, useNA = "a")
prop.table(table(PLASMA_URINE_df$Grade_simple, useNA = "a")) * 100
table(PLASMA_URINE_df$Grade_simple, PLASMA_URINE_df$TYPE, useNA = "a")
#NP
table(LAVAGE_df$Grade_simple, useNA = "a")
prop.table(table(LAVAGE_df$Grade_simple, useNA = "a")) * 100
table(LAVAGE_df$Grade_simple, LAVAGE_df$TYPE, useNA = "a")

#Age#######################################
#test normalcy - can i use averages?
shapiro.test(LAVAGE_df$Age) #normal
shapiro.test(PLASMA_URINE_df$Age) #normal
shapiro.test(TUMOR_df$Age)#normal
#get averages and min max - tumor
round(mean(TUMOR_df$Age, na.rm = TRUE), 3)
min(TUMOR_df$Age)
max(TUMOR_df$Age)
#grouped tumor
TUMOR_df %>%
  mutate(TYPE = recode(TYPE,
                            "RSS" = "RSS_BENIGN",
                            "BENIGN" = "RSS_BENIGN")) %>%
  group_by(TYPE) %>%
  summarise(
    mean_age = round(mean(Age, na.rm = TRUE), 2),
    min_age  = min(Age, na.rm = TRUE),
    max_age  = max(Age, na.rm = TRUE),
    n        = n()
  )


#get averages and min max - URINE/ plasma
round(mean(PLASMA_URINE_df$Age, na.rm = TRUE), 3)
min(PLASMA_URINE_df$Age)
max(PLASMA_URINE_df$Age)
#grouped tumor
PLASMA_URINE_df %>%
  group_by(TYPE) %>%
  summarise(
    mean_age = round(mean(Age, na.rm = TRUE), 2),
    min_age  = min(Age, na.rm = TRUE),
    max_age  = max(Age, na.rm = TRUE),
    n        = n()
  )

#get averages and min max - LAVAGE
round(mean(LAVAGE_df$Age, na.rm = TRUE), 3)
min(LAVAGE_df$Age)
max(LAVAGE_df$Age)
#grouped tumor
LAVAGE_df %>%
  group_by(TYPE) %>%
  summarise(
    mean_age = round(mean(Age, na.rm = TRUE), 2),
    min_age  = min(Age, na.rm = TRUE),
    max_age  = max(Age, na.rm = TRUE),
    n        = n()
  )
#CA125#######################################
#test normalcy - can i use averages?
shapiro.test(LAVAGE_df$CA125_num) #not normal
shapiro.test(PLASMA_URINE_df$CA125_num) #not normal
shapiro.test(TUMOR_df$CA125_num)#not normal
#get median and ci - tumors
median(TUMOR_df$CA125_num, na.rm = TRUE)
quantile(TUMOR_df$CA125_num, 0.25, na.rm = TRUE)
quantile(TUMOR_df$CA125_num, 0.75, na.rm = TRUE)
TUMOR_df %>%
  group_by(TYPE) %>%
  summarise(sprintf(
    "%.2f (%.2f–%.2f)",
    median_CA125 = median(CA125_num, na.rm = TRUE),
    Q1 = quantile(CA125_num, 0.25, na.rm = TRUE),
    Q3 = quantile(CA125_num, 0.75, na.rm = TRUE) ),
    n = n()
  )
#get median and ci - urine / plasma
median(PLASMA_URINE_df$CA125_num, na.rm = TRUE)
quantile(PLASMA_URINE_df$CA125_num, 0.25, na.rm = TRUE)
quantile(PLASMA_URINE_df$CA125_num, 0.75, na.rm = TRUE)
PLASMA_URINE_df %>%
  group_by(TYPE) %>%
  summarise(sprintf(
    "%.2f (%.2f–%.2f)",
    median_CA125 = median(CA125_num, na.rm = TRUE),
    Q1 = quantile(CA125_num, 0.25, na.rm = TRUE),
    Q3 = quantile(CA125_num, 0.75, na.rm = TRUE) ),
    n = n()
  )

#get median and ci - lavage
median(LAVAGE_df$CA125_num, na.rm = TRUE)
quantile(LAVAGE_df$CA125_num, 0.25, na.rm = TRUE)
quantile(LAVAGE_df$CA125_num, 0.75, na.rm = TRUE)
LAVAGE_df %>%
  group_by(TYPE) %>%
  summarise(sprintf(
    "%.2f (%.2f–%.2f)",
    median_CA125 = median(CA125_num, na.rm = TRUE),
    Q1 = quantile(CA125_num, 0.25, na.rm = TRUE),
    Q3 = quantile(CA125_num, 0.75, na.rm = TRUE) ),
    n = n()
  )

#mortality #################################
#test normalcy - can i use averages?
shapiro.test(LAVAGE_df$OS) #not normal
shapiro.test(PLASMA_URINE_df$OS) #not normal
shapiro.test(TUMOR_df$OS)#not normal
#get median and ci - tumors
median(TUMOR_df$OS, na.rm = TRUE)
quantile(TUMOR_df$OS, 0.25, na.rm = TRUE)
quantile(TUMOR_df$OS, 0.75, na.rm = TRUE)
TUMOR_df %>% 
  mutate(TYPE = recode(TYPE,
        "RSS" = "RSS_BENIGN",
        "BENIGN" = "RSS_BENIGN")) %>%
  group_by(TYPE) %>%
  summarise(sprintf(
    "%.2f (%.2f–%.2f)",
    median_OS = median(OS, na.rm = TRUE),
    Q1 = quantile(OS, 0.25, na.rm = TRUE),
    Q3 = quantile(OS, 0.75, na.rm = TRUE) ),
    n = n()
  )
#STATUS for tumors
table(PLASMA_URINE_df$TYPE, PLASMA_URINE_df$STATUS, useNA = "a")
table(PLASMA_URINE_df$STATUS, useNA = "a")

#get median and ci - urine
median(PLASMA_URINE_df$OS, na.rm = TRUE)
quantile(PLASMA_URINE_df$OS, 0.25, na.rm = TRUE)
quantile(PLASMA_URINE_df$OS, 0.75, na.rm = TRUE)
PLASMA_URINE_df %>% 
  mutate(TYPE = recode(TYPE,
                       "RSS" = "RSS_BENIGN",
                       "BENIGN" = "RSS_BENIGN")) %>%
  group_by(TYPE) %>%
  summarise(sprintf(
    "%.2f (%.2f–%.2f)",
    median_OS = median(OS, na.rm = TRUE),
    Q1 = quantile(OS, 0.25, na.rm = TRUE),
    Q3 = quantile(OS, 0.75, na.rm = TRUE) ),
    n = n()
  )
#STATUS for urine
table(PLASMA_URINE_df$TYPE, PLASMA_URINE_df$STATUS, useNA = "a")
table(PLASMA_URINE_df$STATUS, useNA = "a")

#get median and ci - lavage
median(LAVAGE_df$OS, na.rm = TRUE)
quantile(LAVAGE_df$OS, 0.25, na.rm = TRUE)
quantile(LAVAGE_df$OS, 0.75, na.rm = TRUE)
LAVAGE_df %>% 
  group_by(TYPE) %>%
  summarise(sprintf(
    "%.2f (%.2f–%.2f)",
    median_OS = median(OS, na.rm = TRUE),
    Q1 = quantile(OS, 0.25, na.rm = TRUE),
    Q3 = quantile(OS, 0.75, na.rm = TRUE) ),
    n = n()
  )
#STATUS for lavage
table(LAVAGE_df$TYPE, LAVAGE_df$STATUS, useNA = "a")
table(LAVAGE_df$STATUS, useNA = "a")

#which kinds of cancer is in the "other" groups histologically?#############
#plasma samples
View(PLASMA_URINE_df[PLASMA_URINE_df$TYPE == "OTHER",])
#lavage samples
LAVAGE_OTHER <- LAVAGE_df[LAVAGE_df$TYPE == "OTHER",]
table(LAVAGE_OTHER$TYPE_long)
#lavage samples benign
LAVAGE_BEN <- LAVAGE_df[LAVAGE_df$TYPE == "BENIGN",]
table(LAVAGE_BEN$Diagnosis)
#lavage samples endometrial
LAVAGE_BEN <- LAVAGE_df[LAVAGE_df$TYPE == "ENDOMETRIAL CANCER",]
table(LAVAGE_BEN$Diagnosis)
