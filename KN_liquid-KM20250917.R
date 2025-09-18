#KN-liquid - np
#data - 20250916, changed KN-97 from the repeated expreriment
#survival analysis
Sys.setenv(LANG = "en")
#libraries
library(survival)
library(survminer)
library(tidyverse)
library(grid)
library(gridExtra)
library(patchwork)
library(magick)
library(survivalROC)
library(purrr)
library(broom)
library(dplyr)
library(timeROC)
#set wd for plots
setwd("C:/Users/Ieva/rprojects/outputs_all/")
#data
KN_data_full <- readRDS("C:/Users/Ieva/rprojects/OTHER DATA/KN_LIQUID/KN_FULL_LIQUID_R250916.RDS")
# fix colnames
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
genes_cont <- c(raiska_np, raiska_TISSUE, raiska_urine)
raiska_np_f <- paste0(raiska_np, "_f")
#chek
tail(KN_data_full[, colnames(KN_data_full) %in% raiska_p])
#make df of only cases with survival
ALL_SURV_LIQUID <- KN_data_full%>%
  filter(!is.na(OS), !is.na(STATUS)) #58 observations
table(ALL_SURV_LIQUID$STATUS, useNA = "a") #19 dead, 39 NOT
#make factor 
ALL_SURV_LIQUID <- ALL_SURV_LIQUID %>%
  mutate(
    across(all_of(genes_cont),
           ~ factor(if_else(. > median(., na.rm = TRUE), "High", "Low")),
           .names = "{.col}_f")
  )
#make only HGSOC SURV DF
HGSOC_SURV_LIQUID <- ALL_SURV_LIQUID %>%
  filter(Grupė_Ieva == "HGSOC") #41 left#
table(HGSOC_SURV_LIQUID$STATUS, useNA = "a") #17 dead
#make only oc surv df
OC_SURV_LIQUID <- ALL_SURV_LIQUID %>%
  filter(Grupė_Ieva != "Benign") #55 left#
table(OC_SURV_LIQUID$STATUS, useNA = "a") #19 dead, 36 NOT
##COX REGRESION, UNIVARIATE, using continuous variables######################
##NP, ALL ######################################################################
univ_results_np_all <- lapply(raiska_np, function(gene) {
  formula <- as.formula(paste("Surv(OS, STATUS) ~", gene))
  cox_model <- coxph(formula, data = ALL_SURV_LIQUID)
  sum_cox <- summary(cox_model)
  
  # Extract hazard ratio, 95% CI, and p-value
  if (!is.null(sum_cox$conf.int)) {
    data.frame(
      HR = sum_cox$conf.int[,"exp(coef)"],
      lower95 = sum_cox$conf.int[,"lower .95"],
      upper95 = sum_cox$conf.int[,"upper .95"],
      pvalue = sum_cox$coefficients[,"Pr(>|z|)"]
    )
  } else {
    data.frame(HR = NA, lower95 = NA, upper95 = NA, pvalue = NA)
  }
})

# Combine results
univ_df_np_all <- do.call(rbind, univ_results_np_all)
univ_df_np_all$Gene <- raiska_np_f
univ_df_np_all <- univ_df_np_all[, c("Gene", "HR", "lower95", "upper95", "pvalue")]
print(univ_df_np_all) #HES1 DLL1

##NP, OC ######################################################################
univ_results_np_oc <- lapply(raiska_np, function(gene) {
  formula <- as.formula(paste("Surv(OS, STATUS) ~", gene))
  cox_model <- coxph(formula, data = OC_SURV_LIQUID)
  sum_cox <- summary(cox_model)
  
  # Extract hazard ratio, 95% CI, and p-value
  if (!is.null(sum_cox$conf.int)) {
    data.frame(
      HR = sum_cox$conf.int[,"exp(coef)"],
      lower95 = sum_cox$conf.int[,"lower .95"],
      upper95 = sum_cox$conf.int[,"upper .95"],
      pvalue = sum_cox$coefficients[,"Pr(>|z|)"]
    )
  } else {
    data.frame(HR = NA, lower95 = NA, upper95 = NA, pvalue = NA)
  }
})

# Combine results
univ_df_np_OC <- do.call(rbind, univ_results_np_oc)
univ_df_np_OC$Gene <- raiska_np_f
univ_df_np_OC <- univ_df_np_OC[, c("Gene", "HR", "lower95", "upper95", "pvalue")]
print(univ_df_np_OC) #HES1 DLL1 almost significant
##NP, HGSOC ######################################################################
univ_results_np_HGSOC <- lapply(raiska_np, function(gene) {
  formula <- as.formula(paste("Surv(OS, STATUS) ~", gene))
  cox_model <- coxph(formula, data = HGSOC_SURV_LIQUID)
  sum_cox <- summary(cox_model)
  
  # Extract hazard ratio, 95% CI, and p-value
  if (!is.null(sum_cox$conf.int)) {
    data.frame(
      HR = sum_cox$conf.int[,"exp(coef)"],
      lower95 = sum_cox$conf.int[,"lower .95"],
      upper95 = sum_cox$conf.int[,"upper .95"],
      pvalue = sum_cox$coefficients[,"Pr(>|z|)"]
    )
  } else {
    data.frame(HR = NA, lower95 = NA, upper95 = NA, pvalue = NA)
  }
})

# Combine results
univ_df_np_HGSOC <- do.call(rbind, univ_results_np_HGSOC)
univ_df_np_HGSOC$Gene <- raiska_np_f
univ_df_np_HGSOC <- univ_df_np_HGSOC[, c("Gene", "HR", "lower95", "upper95", "pvalue")]
print(univ_df_np_HGSOC) #HES1 DLL1 almost significant

#LONG RANK, KM PLOTS, with univariable COX HR and CI######################
##NP, ALL ##########################################
#make factor 
ALL_SURV_LIQUID <- ALL_SURV_LIQUID %>%
  mutate(
    across(all_of(genes_cont),
           ~ factor(if_else(. > median(., na.rm = TRUE), "High", "Low")),
           .names = "{.col}_f")
  )

#KM plot with univariable
plots <- list()
for (gene in raiska_np_f) {
  # Clean up gene name (remove "_f")
  gene_clean <- sub("_f$", "", gene)
  
  # Fit KM survival curve
  fit <- survfit(as.formula(paste("Surv(OS, STATUS) ~", gene)), 
                 data = ALL_SURV_LIQUID)
  
  # Extract Cox HR + CI
  hr_row <- univ_df_np_all[univ_df_np_all$Gene == gene, ]
  hr_text <- sprintf("HR = %.2f (95%% CI: %.2f–%.2f)", 
                     hr_row$HR, hr_row$lower95, hr_row$upper95)
  
  # KM log-rank p-value
  survdiff_res <- survdiff(as.formula(paste("Surv(OS, STATUS) ~", gene)), 
                           data = ALL_SURV_LIQUID)
  pval_km <- 1 - pchisq(survdiff_res$chisq, length(survdiff_res$n) - 1)
  pval_text <- sprintf("Log-rank p = %.3f", pval_km)
  
  # Subtitle (HR + CI + p-value)
  subtitle_text <- paste0(hr_text, "\n", pval_text)
  
  # Make legend labels with italic gene name
  legend_labels <- c(
    bquote(italic(.(gene_clean)) ~ " Low"),
    bquote(italic(.(gene_clean)) ~ " High")
  )
  
  # Plot
  p <- ggsurvplot(
    fit,
    data = ALL_SURV_LIQUID,
    pval = FALSE,
    risk.table = TRUE,
    legend.title = "",
    palette = c("blue", "red")
  )$plot +
    scale_color_manual(
      values = c("blue", "red"),
      labels = legend_labels
    ) +
    labs(
      title = bquote(italic(.(gene_clean))),
      subtitle = subtitle_text
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(size = 10, hjust = 0.5, lineheight = 1.1)
    )
  
  plots[[gene]] <- p
}

# Combine plots
combined_plot <- wrap_plots(plots, ncol = 4)+
  plot_annotation(
    title = "Survival analysis of ALL Ovarian cases, Uterine lavage",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
    )
  )
combined_plot
ggsave(
  filename = "KM_hr_np_liquid_250917.png",  # output file name
  plot = combined_plot,               # the patchwork plot object
  width = 16,                         # width in inches
  height = 5,                        # height in inches
  dpi = 300                           # resolution
)

##NP, OC ##########################################
#KM plot with univariable
plots2 <- list()
for (gene in raiska_np_f) {
  # Clean up gene name (remove "_f")
  gene_clean <- sub("_f$", "", gene)
  
  # Fit KM survival curve
  fit <- survfit(as.formula(paste("Surv(OS, STATUS) ~", gene)), 
                 data = OC_SURV_LIQUID)
  
  # Extract Cox HR + CI
  hr_row <- univ_df_np_OC[univ_df_np_OC$Gene == gene, ]
  hr_text <- sprintf("HR = %.2f (95%% CI: %.2f–%.2f)", 
                     hr_row$HR, hr_row$lower95, hr_row$upper95)
  
  # KM log-rank p-value
  survdiff_res <- survdiff(as.formula(paste("Surv(OS, STATUS) ~", gene)), 
                           data = OC_SURV_LIQUID)
  pval_km <- 1 - pchisq(survdiff_res$chisq, length(survdiff_res$n) - 1)
  pval_text <- sprintf("Log-rank p = %.3f", pval_km)
  
  # Subtitle (HR + CI + p-value)
  subtitle_text <- paste0(hr_text, "\n", pval_text)
  
  # Make legend labels with italic gene name
  legend_labels <- c(
    bquote(italic(.(gene_clean)) ~ " Low"),
    bquote(italic(.(gene_clean)) ~ " High")
  )
  
  # Plot
  p <- ggsurvplot(
    fit,
    data = OC_SURV_LIQUID,
    pval = FALSE,
    risk.table = TRUE,
    legend.title = "",
    palette = c("blue", "red")
  )$plot +
    scale_color_manual(
      values = c("blue", "red"),
      labels = legend_labels
    ) +
    labs(
      title = bquote(italic(.(gene_clean))),
      subtitle = subtitle_text
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(size = 10, hjust = 0.5, lineheight = 1.1)
    )
  
  plots2[[gene]] <- p
}

# Combine plots
combined_plot2 <- wrap_plots(plots2, ncol = 4)+
  plot_annotation(
    title = "Survival analysis of Ovarian Cancer cases, Uterine lavage",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
    )
  )
combined_plot2
ggsave(
  filename = "KM_hr_np_liquid_OC_250917.png",  # output file name
  plot = combined_plot2,               # the patchwork plot object
  width = 16,                         # width in inches
  height = 5,                        # height in inches
  dpi = 300                           # resolution
)
##NP, HGSOC##########################################
#KM plot with univariable
plots3 <- list()
for (gene in raiska_np_f) {
  # Clean up gene name (remove "_f")
  gene_clean <- sub("_f$", "", gene)
  
  # Fit KM survival curve
  fit <- survfit(as.formula(paste("Surv(OS, STATUS) ~", gene)), 
                 data = HGSOC_SURV_LIQUID)
  
  # Extract Cox HR + CI
  hr_row <- univ_df_np_HGSOC[univ_df_np_HGSOC$Gene == gene, ]
  hr_text <- sprintf("HR = %.2f (95%% CI: %.2f–%.2f)", 
                     hr_row$HR, hr_row$lower95, hr_row$upper95)
  
  # KM log-rank p-value
  survdiff_res <- survdiff(as.formula(paste("Surv(OS, STATUS) ~", gene)), 
                           data = HGSOC_SURV_LIQUID)
  pval_km <- 1 - pchisq(survdiff_res$chisq, length(survdiff_res$n) - 1)
  pval_text <- sprintf("Log-rank p = %.3f", pval_km)
  
  # Subtitle (HR + CI + p-value)
  subtitle_text <- paste0(hr_text, "\n", pval_text)
  
  # Make legend labels with italic gene name
  legend_labels <- c(
    bquote(italic(.(gene_clean)) ~ " Low"),
    bquote(italic(.(gene_clean)) ~ " High")
  )
  
  # Plot
  p <- ggsurvplot(
    fit,
    data = HGSOC_SURV_LIQUID,
    pval = FALSE,
    risk.table = TRUE,
    legend.title = "",
    palette = c("blue", "red")
  )$plot +
    scale_color_manual(
      values = c("blue", "red"),
      labels = legend_labels
    ) +
    labs(
      title = bquote(italic(.(gene_clean))),
      subtitle = subtitle_text
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(size = 10, hjust = 0.5, lineheight = 1.1)
    )
  
  plots3[[gene]] <- p
}

# Combine plots
combined_plot3 <- wrap_plots(plots3, ncol = 4)+
  plot_annotation(
    title = "Survival analysis of HGSOC cases, Uterine lavage",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
    )
  )
combined_plot3
ggsave(
  filename = "KM_hr_np_liquid_HGSOC_250917.png",  # output file name
  plot = combined_plot3,               # the patchwork plot object
  width = 16,                         # width in inches
  height = 5,                        # height in inches
  dpi = 300                           # resolution
)

#one model, np, all #####################################
gene_data <- ALL_SURV_LIQUID[, colnames(ALL_SURV_LIQUID) %in% raiska_np]
# Example with 4 biomarkers
cox_model_ALL <- coxph(
  Surv(OS, STATUS) ~  NOTCH2_NP + CTNNB1_NP + DLL1_NP + HES1_NP,
  data = ALL_SURV_LIQUID
)
summary(cox_model_ALL)
#get coeficients
coefs <- coef(cox_model_ALL)
# Risk score = sum( gene_expression * coefficient )
risk_scores_test <- rowSums(sweep(gene_data, 2, coefs, "*"))
# View the risk scores
print(risk_scores_test) # now I have some risk scores
#add risk scores to the clin_df_joined_test
ALL_SURV_LIQUID$RiskScore <- risk_scores_test
#create df wih survival data
surv_df_np_all <- ALL_SURV_LIQUID[, colnames(ALL_SURV_LIQUID) %in%
                                    c("OS", "STATUS", raiska_np, "RiskScore", "patient_id_aud")]

rownames(surv_df_np_all) <- surv_df_np_all$patient_id_aud

# Calculate the median risk score
median_risk <- median(surv_df_np_all$RiskScore, na.rm = TRUE) #-3.007294
# Create a new factor column based on the median value
surv_df_np_all$RiskGroup <- ifelse(surv_df_np_all$RiskScore <= median_risk,
                                   "Low Risk", "High Risk")
#Create a survival object
surv_object <- Surv(time = surv_df_np_all$OS,
                    event = surv_df_np_all$STATUS )

# Fit a Kaplan-Meier model
km_fit <- survfit(surv_object ~ RiskGroup, data = surv_df_np_all)
# Plot the Kaplan-Meier curve using ggsurvplot
test_survplot <- ggsurvplot(km_fit, data = surv_df_np_all, 
                            pval = TRUE,  # Show p-value of the log-rank test
                            risk.table = TRUE,  # Add risk table below the plot
                            title = "Kaplan-Meier kreivė: Didelės vs. mažos rizikos atvejai KV audinių imtyje",
                            xlab = "Bendras išgyvenamumo laikas",
                            ylab = "Išgyvenamumo tikimybė",
                            palette = c("turquoise", "deeppink"),  # Color palette for groups
                            legend.title = "Rizikos grupė", 
                            legend.labs = c("Žema rizika", "Didelė rizika"))
test_survplot
#save
png("C:/Users/Ieva/rprojects/outputs_all/KM_nop_all_gene_sig_20150915.png",
    width = 800, height = 600, res = 100) # width and height in pixels, resolution in dpi
test_survplot #
dev.off() # Close the PNG device

#forest 
forst_all_np <- ggforest(
  cox_model_ALL, 
  data = ALL_SURV_LIQUID,
  main = "Hazard ratios",
  cpositions = c(0.02, 0.22, 0.4), # adjust column positions
  fontsize = 1.0
)

ggsave("forestplot_np_all.png", plot = forst_all_np, width = 7, height = 5, dpi = 300)

#one model, np, OC##################
# Example with 4 biomarkers
cox_model_OC <- coxph(
  Surv(OS, STATUS) ~  NOTCH2_NP + CTNNB1_NP + DLL1_NP + HES1_NP,
  data = OC_SURV_LIQUID
)
summary(cox_model_OC)

#forest 
forst_oc_np <- ggforest(
  cox_model_OC, 
  data = OC_SURV_LIQUID,
  main = "Hazard ratios, OC",
  cpositions = c(0.02, 0.22, 0.4), # adjust column positions
  fontsize = 1.0
)

ggsave("forestplot_np_oc.png", plot = forst_oc_np, width = 7, height = 5, dpi = 300)


#one model, np, HGSOC##################
# Example with 4 biomarkers
cox_model_HGSOC <- coxph(
  Surv(OS, STATUS) ~  NOTCH2_NP + CTNNB1_NP + DLL1_NP + HES1_NP,
  data = HGSOC_SURV_LIQUID
)
summary(cox_model_HGSOC)

#forest 
forst_hgsoc_np <- ggforest(
  cox_model_HGSOC, 
  data = HGSOC_SURV_LIQUID,
  main = "Hazard ratios, HGSOC",
  cpositions = c(0.02, 0.22, 0.4), # adjust column positions
  fontsize = 1.0
)

ggsave("forestplot_np_HGSoc.png", plot = forst_hgsoc_np, width = 7, height = 5, dpi = 300)

#URINE DATA #############################
table(ALL_SURV_LIQUID$STATUS, ALL_SURV_LIQUID$NOTCH2_URINE_f,
      useNA = "a") #1 dead, 8 NOT 
table(ALL_SURV_LIQUID$STATUS, ALL_SURV_LIQUID$HES1_URINE_f,
      useNA = "a") #2 dead, 12 NOT
table(ALL_SURV_LIQUID$STATUS, ALL_SURV_LIQUID$CTNNB1_URINE_f,
      useNA = "a") #1 dead, 6 NOT
table(ALL_SURV_LIQUID$STATUS, ALL_SURV_LIQUID$DLL1_URINE_f,
      useNA = "a") #0 dead, 3 NOT

#try urine KM plots, all ######################################## 
# Fit a Kaplan-Meier model NOTCH2
km_fit_notch2_urine <- survfit(Surv(OS, STATUS) ~ NOTCH2_URINE_f, data = ALL_SURV_LIQUID)
km_fit_notch2_urine
# Plot the Kaplan-Meier curve using ggsurvplot
test_survplot_notch2_urine <- ggsurvplot(km_fit_notch2_urine, data = ALL_SURV_LIQUID, 
                            pval = TRUE,  # Show p-value of the log-rank test
                            risk.table = TRUE,  # Add risk table below the plot
                            title = "NOTCH2_URINE KV audinių imtyje",
                            xlab = "Bendras išgyvenamumo laikas",
                            ylab = "Išgyvenamumo tikimybė",
                            palette = c("turquoise", "deeppink"),  # Color palette for groups
                            legend.title = "Rizikos grupė", 
                            legend.labs = c("Maža raiška", "Didelė raiška"))
test_survplot_notch2_urine

# Fit a Kaplan-Meier model HES1
km_fit_hes1_urine <- survfit(Surv(OS, STATUS) ~ HES1_URINE_f, data = ALL_SURV_LIQUID)
km_fit_hes1_urine
# Plot the Kaplan-Meier curve using ggsurvplot
test_survplot_hes1_urine <- ggsurvplot(km_fit_hes1_urine, data = ALL_SURV_LIQUID, 
                                         pval = TRUE,  # Show p-value of the log-rank test
                                         risk.table = TRUE,  # Add risk table below the plot
                                         title = "HES1_URINE KV audinių imtyje",
                                         xlab = "Bendras išgyvenamumo laikas",
                                         ylab = "Išgyvenamumo tikimybė",
                                         palette = c("turquoise", "deeppink"),  # Color palette for groups
                                         legend.title = "Rizikos grupė", 
                                         legend.labs = c("Maža raiška", "Didelė raiška"))
test_survplot_hes1_urine

# Fit a Kaplan-Meier model CTNNB1
km_fit_cttnb1_urine <- survfit(Surv(OS, STATUS) ~ CTNNB1_URINE_f, data = ALL_SURV_LIQUID)
km_fit_cttnb1_urine
# Plot the Kaplan-Meier curve using ggsurvplot
test_survplot_ctnnb1_urine <- ggsurvplot(km_fit_cttnb1_urine, data = ALL_SURV_LIQUID, 
                                       pval = TRUE,  # Show p-value of the log-rank test
                                       risk.table = TRUE,  # Add risk table below the plot
                                       title = "CTNNB1_URINE KV audinių imtyje",
                                       xlab = "Bendras išgyvenamumo laikas",
                                       ylab = "Išgyvenamumo tikimybė",
                                       palette = c("turquoise", "deeppink"),  # Color palette for groups
                                       legend.title = "Rizikos grupė", 
                                       legend.labs = c("Maža raiška", "Didelė raiška"))
test_survplot_ctnnb1_urine

# Combine side by side
combined <- test_survplot_ctnnb1_urine$plot | test_survplot_hes1_urine$plot | test_survplot_notch2_urine$plot
# Save
ggsave("KM_combined_urine.png", combined, width = 15, height = 4, dpi = 300)

#PLASMA DATA #############################
table(ALL_SURV_LIQUID$STATUS, ALL_SURV_LIQUID$NOCTH2_P,
      useNA = "a") #2 dead, 12 NOT 
table(ALL_SURV_LIQUID$STATUS, ALL_SURV_LIQUID$HES1_P,
      useNA = "a") #2 dead, 12 NOT
table(ALL_SURV_LIQUID$STATUS, ALL_SURV_LIQUID$CTNNB1_P,
      useNA = "a") #2 dead, 11 NOT
table(ALL_SURV_LIQUID$STATUS, ALL_SURV_LIQUID$DLL1_P,
      useNA = "a") #2 dead, 12 NOT
plots_plasma <- list()

for (gene in raiska_p) {
  # KM fit
  fit <- survfit(as.formula(paste("Surv(OS, STATUS) ~", gene)),
                 data = ALL_SURV_LIQUID)
  
  # KM plot with log-rank p-value
  p <- ggsurvplot(
    fit,
    data = ALL_SURV_LIQUID,
    pval = TRUE,
    risk.table = FALSE,
    legend.title = "",
    palette = c("blue", "red"),
    title = gene# italic gene name, strip suffix
  )$plot +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  plots_plasma[[gene]] <- p
}

# Combine 4 plots into 2x2 grid
combined_plasma <- (plots_plasma[[1]] | plots_plasma[[2]]) / (plots_plasma[[3]] | plots_plasma[[4]])

# Save
ggsave("KM_plasma_combined.png", combined_plasma, width = 10, height = 8, dpi = 300)