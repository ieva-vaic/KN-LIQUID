#KN-  liquid 2026 04 17, 2026 05 14, 2026-06-19
#concordance of 
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
library(ggpubr)
library(patchwork)
#read RDS
LIQUID_DF_final <- readRDS("C:/Users/Ieva/rprojects/OTHER DATA/KN_LIQUID/liquid_20260415.RDS")
#leave tumor only
tumor_df <- LIQUID_DF_final%>%
  filter(!is.na(CTNNB1_TUMOR)) #65 cases
#tissue vs NP##########################################
##normalcy of genes ################################
#np
shapiro.test(tumor_df$NOTCH2_NP)#0.9807 ->0.9689
shapiro.test(tumor_df$DLL1_NP)#0.1804 -> 0.1632
shapiro.test(tumor_df$HES1_NP)#0.2917 -> 0.2694
shapiro.test(tumor_df$CTNNB1_NP)#0.002585 #not normal -> 0.001881
#tissue
shapiro.test(tumor_df$NOTCH2_TUMOR)#0.3496 -> 0.796
shapiro.test(tumor_df$DLL1_TUMOR)#0.2464 
shapiro.test(tumor_df$HES1_TUMOR)#0.2383 ->0.2623
shapiro.test(tumor_df$CTNNB1_TUMOR)#0.05054  ->0.1678
#spearmann for CTNNB1

##make correltation plot####################################
#make correlation and plots between np and tissue
hes1_corr <- cor.test(tumor_df$HES1_NP, tumor_df$HES1_TUMOR, method = "pearson")#0.04355 -> 0.07202
p.adjust(hes1_corr$p.value, method = "BH")

cor.test(tumor_df$NOTCH2_NP, tumor_df$NOTCH2_TUMOR, method = "pearson")#0.3259
cor.test(tumor_df$CTNNB1_NP, tumor_df$CTNNB1_TUMOR, method = "spearman")#0.9877
cor.test(tumor_df$DLL1_NP, tumor_df$DLL1_TUMOR, method = "pearson")#0.2198
table(tumor_df$TYPE, useNA = "a")
# genes
genes <- c("CTNNB1", "HES1", "DLL1", "NOTCH2")

# gene requiring Spearman
spearman_gene <- "CTNNB1"

make_corr_plot <- function(gene, df) {
  
  tumor_col <- paste0(gene, "_TUMOR")
  np_col <- paste0(gene, "_NP")
  
  plot_df <- df %>%
    select(
      tumor = all_of(tumor_col),
      np = all_of(np_col)
    ) %>%
    drop_na()
  
  # choose method
  method <- ifelse(gene == spearman_gene,
                   "spearman",
                   "pearson")
  
  # correlation test
  cor_res <- cor.test(
    plot_df$tumor,
    plot_df$np,
    method = method
  )
  
  r_val <- round(cor_res$estimate, 2)
  p_val <- signif(cor_res$p.value, 2)
  
  label_text <- ifelse(
    method == "spearman",
    paste0("\u03C1 = ", r_val, "\n",
           "p = ", p_val),
    paste0("r = ", r_val, "\n",
           "p = ", p_val)
  )
  
  ggplot(plot_df, aes(x = tumor, y = np)) +
    
    geom_point(
      size = 2.2,
      alpha = 0.8
    ) +
    
    geom_smooth(
      method = "lm",
      se = TRUE,
      linewidth = 0.8
    ) +
    
    annotate(
      "text",
      x = Inf,
      y = -Inf,
      label = label_text,
      hjust = 1.1,
      vjust = -0.5,
      size = 4.5
    ) +
    
    labs(
      title = gene,
      x = "Tumor expression",
      y = "Uterine lavage expression"
    ) +
    
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(
        face = "italic",
        hjust = 0.5,
        size = 15
      )
    )
}

# create plots
plot_list <- lapply(genes, make_corr_plot, df = tumor_df)

# combine
combined_plot <- wrap_plots(plot_list, ncol = 2)

combined_plot

#save
ggsave(
  filename = "c:/Users/Ieva/rprojects/outputs_all/LIQUID/np_tumor_corr_combined20260619.png",
  plot = combined_plot,
  width = 9,
  height =8,
  dpi = 350
)
##Bland-altman plot############################
genes <- c("CTNNB1", "HES1", "DLL1", "NOTCH2")

# ---- reshape ----
long_df <- bind_rows(lapply(genes, function(g) {
  
  tumor_col <- paste0(g, "_TUMOR")
  np_col    <- paste0(g, "_NP")
  
  data.frame(
    gene = g,
    mean_val = (tumor_df[[tumor_col]] + tumor_df[[np_col]]) / 2,
    diff_val = (tumor_df[[tumor_col]] - tumor_df[[np_col]])
  )
}))

# ---- stats ----
stats <- long_df %>%
  group_by(gene) %>%
  summarise(
    bias = mean(diff_val, na.rm = TRUE),
    sd = sd(diff_val, na.rm = TRUE),
    loa_upper = bias + 1.96 * sd,
    loa_lower = bias - 1.96 * sd
  ) %>%
  mutate(
    label = paste0(
      "bias = ", round(bias, 2), "\n",
      "+1.96 SD = ", round(loa_upper, 2), "\n",
      "-1.96 SD = ", round(loa_lower, 2)
    )
  )

plot_df <- long_df %>%
  left_join(stats, by = "gene")
blad_p <- ggplot(plot_df, aes(x = mean_val, y = diff_val)) +
  geom_point(alpha = 0.6) +
  
  geom_hline(aes(yintercept = bias), color = "blue", linewidth = 0.8) +
  geom_hline(aes(yintercept = loa_upper), linetype = "dashed", color = "red") +
  geom_hline(aes(yintercept = loa_lower), linetype = "dashed", color = "red") +
  
  geom_text(
    data = stats,
    aes(x = -Inf, y = Inf, label = label),
    hjust = -0.05, vjust = 1.1, size = 3.2
  ) +
  
  facet_wrap(
    ~gene,
    labeller = as_labeller(function(x) {
      paste0("italic(", x, ")")
    }, label_parsed)
  ) +
  
  labs(
    title = "Bland–Altman Plots: Tumor vs UL Expression",
    x = "Mean of Tumor and UL",
    y = "Difference (Tumor - UL)"
  ) +
  theme_minimal()
blad_p
#save
ggsave(
  filename = "c:/Users/Ieva/rprojects/outputs_all/LIQUID/np_tumor_blad-altman_combined20260619.png",
  plot = blad_p,
  width = 8,
  height =8,
  dpi = 300
)
#PLASMA / URINE/ NP/ TISSUE##########################################
#leave URINE/PLASMA only
#remove cases that is NA in their type: 
LIQUID_DF_15 <- LIQUID_DF_final %>%
  filter(!is.na(HES1_URINE))
#detetction rates plasma
table(LIQUID_DF_15$NOTCH2_P, useNA = "a")#1
table(LIQUID_DF_15$DLL1_P, useNA = "a")#0
table(LIQUID_DF_15$HES1_P, useNA = "a")#4
table(LIQUID_DF_15$CTNNB1_P, useNA = "a")#3
#detection rates urine
table(LIQUID_DF_15$NOTCH2_URINE, useNA = "a")#10
table(LIQUID_DF_15$DLL1_URINE, useNA = "a")#3
table(LIQUID_DF_15$HES1_URINE, useNA = "a")#15
table(LIQUID_DF_15$CTNNB1_URINE, useNA = "a")#13
#detection rates urine
table(LIQUID_DF_15$NOTCH2_NP, useNA = "a")#14
table(LIQUID_DF_15$DLL1_NP, useNA = "a")#14
table(LIQUID_DF_15$HES1_NP, useNA = "a")#14
table(LIQUID_DF_15$CTNNB1_NP, useNA = "a")#14
##4 way comparison plots ############################
#notch2
NOTCH2 <- c("NOTCH2_URINE", "NOTCH2_NP", "NOTCH2_P_norm","NOTCH2_TUMOR")
notch2_long <- LIQUID_DF_15  %>%
  dplyr::select(c("Laboratorinis kodas", NOTCH2)) %>%
  pivot_longer(
    cols = NOTCH2,
    names_to = "gene",
    values_to = "expression"
  )
#set order
notch2_long$gene <- factor(
  notch2_long$gene,
  levels = c("NOTCH2_P_norm", "NOTCH2_URINE", "NOTCH2_NP", "NOTCH2_TUMOR"),
  labels = c("Plasma", "Urine", "Uterine lavage", "Tumor")
)

notch2 <- ggplot(notch2_long, aes(x = gene, y = expression, group = `Laboratorinis kodas`)) +
  geom_point(size = 2) +
  geom_line(alpha = 0.5) +
  theme_classic() +
  labs(
    x = "Sample type",
    y = expression("Gene expression, normalized to"* italic(" GAPDH")),
    title = expression("Gene expression of "* italic("NOTCH2" )*" across sample types"),
  )
plot(notch2)
#hes1
hes1 <- c("HES1_URINE", "HES1_NP", "HES1_P_norm","HES1_TUMOR")
hes1_long <- LIQUID_DF_15  %>%
  dplyr::select(c("Laboratorinis kodas", hes1)) %>%
  pivot_longer(
    cols = hes1,
    names_to = "gene",
    values_to = "expression"
  )

#set order
hes1_long$gene <- factor(
  hes1_long$gene,
  levels = c("HES1_P_norm","HES1_URINE", "HES1_NP", "HES1_TUMOR") ,
  labels = c("Plasma", "Urine", "Uterine lavage", "Tumor")
)

hes1 <- ggplot(hes1_long, aes(x = gene, y = expression, group = `Laboratorinis kodas`)) +
  geom_point(size = 2) +
  geom_line(alpha = 0.5) +
  theme_classic() +
  labs(
    x = "Sample type",
    y = expression("Gene expression, normalized to"* italic(" GAPDH")),
    title = expression("Gene expression of "* italic("HES1" )*" across sample types"),
  )
plot(hes1)
#CTNNB1
CTNNB1 <- c("CTNNB1_URINE", "CTNNB1_NP", "CTNNB1_P_norm","CTNNB1_TUMOR")
CTNNB1_long <- LIQUID_DF_15  %>%
  dplyr::select(c("Laboratorinis kodas", CTNNB1)) %>%
  pivot_longer(
    cols = CTNNB1,
    names_to = "gene",
    values_to = "expression"
  )
#set order
CTNNB1_long$gene <- factor(
  CTNNB1_long$gene,
  levels = c( "CTNNB1_P_norm","CTNNB1_URINE", "CTNNB1_NP","CTNNB1_TUMOR") ,
  labels = c("Plasma", "Urine", "Uterine lavage", "Tumor")
)
ctnnb1 <- ggplot(CTNNB1_long, aes(x = gene, y = expression, group = `Laboratorinis kodas`)) +
  geom_point(size = 2) +
  geom_line(alpha = 0.5) +
  theme_classic() +
  labs(
    x = "Sample type",
    y = expression("Gene expression, normalized to"* italic(" GAPDH")),
    title = expression("Gene expression of "* italic("CTNNB1" )*" across sample types"),
  )
plot(ctnnb1)
#dll1
DLL1 <- c("DLL1_URINE", "DLL1_NP", "DLL1_P_norm","DLL1_TUMOR")
DLL1_long <- LIQUID_DF_15  %>%
  dplyr::select(c("Laboratorinis kodas", DLL1)) %>%
  pivot_longer(
    cols = DLL1,
    names_to = "gene",
    values_to = "expression"
  )

#set order
DLL1_long$gene <- factor(
  DLL1_long$gene,
  levels = c("DLL1_P_norm", "DLL1_URINE", "DLL1_NP", "DLL1_TUMOR") ,
  labels = c("Plasma", "Urine", "Uterine lavage", "Tumor")
)

dll1 <- ggplot(DLL1_long, aes(x = gene, y = expression, group = `Laboratorinis kodas`)) +
  geom_point(size = 2) +
  geom_line(alpha = 0.5) +
  theme_classic() +
  labs(
    x = "Sample type",
    y = expression("Gene expression, normalized to"* italic(" GAPDH")),
    title = expression("Gene expression of "* italic("DLL1" )*" across sample types"),
  )
plot(dll1)

combined_plot <-
 ( hes1 |  notch2 )/
  (ctnnb1 | dll1)

ggsave(
  filename = "c:/Users/Ieva/rprojects/outputs_all/LIQUID/15_sample_comparisons20260514.png",
  plot = combined_plot,
  width = 12,
  height = 8,
  dpi = 300
)

#mixed affects model anova
library(emmeans)
library(lme4)
library(lmerTest)
#notch2
model_notxh2 <- lmer(
  expression ~ gene + (1 | `Laboratorinis kodas`),
  data = notch2_long
)
summary(model_notxh2)
anova(model_notxh2)#0.0006161 
emm_notch2 <- emmeans(model_notxh2, pairwise ~ gene, adjust = "BH")

#hes1
model_hes1 <- lmer(
  expression ~ gene + (1 | `Laboratorinis kodas`),
  data = hes1_long
)
summary(model_hes1)
anova(model_hes1)#9.025e-09
emm_hes1 <- emmeans(model_hes1, pairwise ~ gene, adjust = "BH")

#DLL1
model_dll1 <- lmer(
  expression ~ gene + (1 | `Laboratorinis kodas`),
  data = DLL1_long
)
summary(model_dll1)
anova(model_dll1)#0.1289
emm_dll1 <- emmeans(model_dll1, pairwise ~ gene, adjust = "BH")

#CTNNB1
model_ctnnb1 <- lmer(
  expression ~ gene + (1 | `Laboratorinis kodas`),
  data = CTNNB1_long
)
summary(model_ctnnb1)
anova(model_ctnnb1)#0.0001466 
emm_ctnnb1 <-emmeans(model_ctnnb1, pairwise ~ gene, adjust = "BH")
emm_ctnnb1

##plots with p values###########################
#CTNNb1
pvals_ctnnb1 <- as.data.frame(emm_ctnnb1$contrasts)

pvals_ctnnb1

pvals_ctnnb1 <- pvals_ctnnb1 %>%
  separate(contrast, into = c("group1", "group2"), sep = " - ") %>%
  mutate(
    y.position = c( -1,#choose max number
                    0,
                    -1.2,
                    -1.3,
                    -0.5,
                    -1.4)
  )
pvals_ctnnb1

pvals_ctnnb1 <- pvals_ctnnb1 %>%
  filter(p.value < 0.05)
pvals_ctnnb1$p.value <- formatC(
  pvals_ctnnb1$p.value,
  format = "f",
  digits = 3
)

CTNNB1 <- ggplot(CTNNB1_long,
                 aes(x = gene,
                     y = expression,
                     group = `Laboratorinis kodas`)) +
  
  geom_point(size = 2) +
  geom_line(alpha = 0.5) +
  
  stat_pvalue_manual(
    pvals_ctnnb1,
    label = "p.value",
    tip.length = 0.01
  ) +
  
  theme_classic() +
  
  labs(
    x = "Sample type",
    y = expression("Gene expression, normalized to " * italic(GAPDH)),
    title = expression(
      "Gene expression of " * italic(CTNNB1) * " across sample types"
    )
  )

plot(CTNNB1)
#hes1
pvals_hes1 <- as.data.frame(emm_hes1$contrasts)

pvals_hes1

pvals_hes1 <- pvals_hes1 %>%
  separate(contrast, into = c("group1", "group2"), sep = " - ") %>%
  mutate(
    y.position = c( 1,#choose max number
                    0.7,
                    2.3,
                    1.5,
                    0.5,
                    1.3)
  )
pvals_hes1

pvals_hes1 <- pvals_hes1 %>%
  filter(p.value < 0.05)
pvals_hes1$p.value <- ifelse(
  pvals_hes1$p.value < 0.001,
  "< 0.001",
  formatC(pvals_hes1$p.value, format = "f", digits = 3)
)

HES1 <- ggplot(hes1_long,
               aes(x = gene,
                   y = expression,
                   group = `Laboratorinis kodas`)) +
  
  geom_point(size = 2) +
  geom_line(alpha = 0.5) +
  
  stat_pvalue_manual(
    pvals_hes1,
    label = "p.value",
    tip.length = 0.01
  ) +
  
  theme_classic() +
  
  labs(
    x = "Sample type",
    y = expression("Gene expression, normalized to " * italic(GAPDH)),
    title = expression(
      "Gene expression of " * italic(HES1) * " across sample types"
    )
  )

plot(HES1)

pvals_notch2 <- as.data.frame(emm_notch2$contrasts)

pvals_notch2

pvals_notch2 <- pvals_notch2 %>%
  separate(contrast, into = c("group1", "group2"), sep = " - ") %>%
  mutate(
    y.position = c( -1,#choose max number
                    -0.7,
                    -2,
                    -1.5,
                    -0.5,
                    -1.3)
  )
pvals_notch2

pvals_notch2 <- pvals_notch2 %>%
  filter(p.value < 0.05)
pvals_notch2$p.value <- formatC(
  pvals_notch2$p.value,
  format = "f",
  digits = 3
)
NOTCH2 <- ggplot(notch2_long,
                 aes(x = gene,
                     y = expression,
                     group = `Laboratorinis kodas`)) +
  
  geom_point(size = 2) +
  geom_line(alpha = 0.5) +
  
  stat_pvalue_manual(
    pvals_notch2,
    label = "p.value",
    tip.length = 0.01
  ) +
  
  theme_classic() +
  
  labs(
    x = "Sample type",
    y = expression("Gene expression, normalized to " * italic(GAPDH)),
    title = expression(
      "Gene expression of " * italic(NOTCH2) * " across sample types"
    )
  )

plot(NOTCH2)

compare_plots <-
  ( HES1 |  NOTCH2 )/
  (CTNNB1 | dll1)

ggsave(
  filename = "c:/Users/Ieva/rprojects/outputs_all/LIQUID/15_sample_pvals20260619.png",
  plot = compare_plots,
  width = 12,
  height = 8,
  dpi = 300
)
