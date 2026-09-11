#KN-  liquid 2026 09 08
#COMPARE CONCENTRATIONS
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
library(patchwork)
library(rstatix)
library(lmerTest)
library(emmeans)

conc <- readxl::read_xlsx("C:/Users/Ieva/rprojects/OTHER DATA/KN_LIQUID/kn CONCENTRATIONS 1 PAGE 20260908.xlsx")
#make summary measures
conc_summary <- conc %>%
  group_by(`SAMPLE TYPE`) %>%
  summarise(
    across(
      c(`NP ng/µL`, `NP 260/280`, `NP 260/230`),
      list(
        min = ~ min(.x, na.rm = TRUE),
        max = ~ max(.x, na.rm = TRUE),
        median = ~ median(.x, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  )

conc_summary

# Summary measures: median and IQR (Q3 − Q1)
conc_summary2 <- conc %>%
  group_by(`SAMPLE TYPE`) %>%
  summarise(
    across(
      c(`NP ng/µL`, `NP 260/280`, `NP 260/230`),
      list(
        median = ~ median(.x, na.rm = TRUE),
        IQR = ~ IQR(.x, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  )

conc_summary2
#find normalcy
# Normality within each sample type
conc %>%
  group_by(`SAMPLE TYPE`) %>%
  summarise(
    n = n(),
    p_normality = shapiro.test(`NP ng/µL`)$p.value,
    .groups = "drop"
  )#not normal

#anova of concentrations
# Overall comparison
conc %>%
  kruskal_test(`NP ng/µL` ~ `SAMPLE TYPE`)


conc_test <- conc %>%
  transmute(
    concentration = `NP ng/µL`,
    sample_type = factor(`SAMPLE TYPE`)
  ) %>%
  filter(!is.na(concentration), !is.na(sample_type))
# Pairwise post hoc comparisons
conc_test %>%
  dunn_test(
    concentration ~ sample_type,
    p.adjust.method = "holm"
  )
#open GAPDH file######################
gapdh <- readxl::read_xlsx("C:/Users/Ieva/rprojects/OTHER DATA/KN_LIQUID/GAPDH only 20260909.xlsx")
head(gapdh)
dat <- gapdh %>%
  mutate(
    Patient_ID = sub("^(KN-[0-9]+).*", "\\1", Sample),
    Patient_ID = factor(Patient_ID),
    type = factor(type)
  ) %>%
  filter(!is.na(GAPDH_MEAN), !is.na(type))
# Verify patient matching before analysis
dat %>% select(Sample, Patient_ID, type)
dat %>%
  group_by(type) %>%
  summarise(
    n = n(),
    mean_Ct = mean(GAPDH_MEAN),
    SD_Ct = sd(GAPDH_MEAN),
    median_Ct = median(GAPDH_MEAN),
    IQR_Ct = IQR(GAPDH_MEAN),
    min_Ct = min(GAPDH_MEAN),
    max_Ct = max(GAPDH_MEAN),
    .groups = "drop"
  )

ggplot(dat, aes(type, GAPDH_MEAN, fill = type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.4) +
  geom_jitter(width = 0.12, alpha = 0.7, size = 2) +
  labs(x = "Sample type", y = "GAPDH mean Ct") +
  theme_classic() +
  theme(legend.position = "none")
#compare paired
fit <- lmer(
  GAPDH_MEAN ~ type + (1 | Patient_ID),
  data = dat
)

anova(fit)  # Overall sample-type effect

# Pairwise differences in Ct, confidence intervals and adjusted p-values
emm <- emmeans(fit, ~ type)
summary(pairs(emm, adjust = "tukey"), infer = c(TRUE, TRUE))

# Check residual spread and approximate normality
plot(fit)
qqnorm(resid(fit))
qqline(resid(fit))
