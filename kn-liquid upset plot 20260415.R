#upset plot for kn-liquid 20260415
library(ComplexUpset)
library(ggplot2)
Sys.setenv(LANG = "en")
library(tidyverse)
library(readxl)
library(openxlsx)
#make upset plot of sample overlap
LIQUID_DF_final <- readRDS("C:/Users/Ieva/rprojects/OTHER DATA/KN_LIQUID/liquid_20260415.RDS")
LIQUID_DF_final <-
  LIQUID_DF_final %>%
  mutate(TISSUE_SAMPLE = ifelse(!is.na(NOTCH2_TUMOR), 1, 0))%>%
  mutate(LAVAGE_SAMPLE = ifelse(!is.na(NOTCH2_NP), 1, 0))%>%
  mutate(URINE_SAMPLE = ifelse(!is.na(HES1_URINE), 1, 0))%>%
  mutate(PLASMA_SAMPLE = ifelse(!is.na(HES1_URINE), 1, 0))
LIQUID_DF_final$sample_combination
df_wide <- LIQUID_DF_final %>%
  mutate(value = 1) %>%
  pivot_wider(
    names_from = c(TYPE, sample_combination),
    values_from = value,
    values_fill = 0
  )


upset(
  LIQUID_DF_final,
  intersect = c("PLASMA_SAMPLE", "URINE_SAMPLE", "LAVAGE_SAMPLE", "TISSUE_SAMPLE"),
  
  name = "Sample type",
  
  set_sizes = (
    upset_set_size() +
      geom_text(
        aes(label = after_stat(count)),
        stat = "count",
        hjust = -0.2,
        size = 3
      )
  ),
  
  base_annotations = list(
    "Tumor type" = intersection_ratio(
      mapping = aes(fill = TYPE)
    )
  )
)

