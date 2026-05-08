#KN-  liquid 2026 05 05, 07
#FINAL PLOTS FOR LIQUID PAPER - NP ONLY - ROC COMBINATIONS
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
#read RDS
LIQUID_DF_final <- readRDS("C:/Users/Ieva/rprojects/OTHER DATA/KN_LIQUID/liquid_20260415.RDS")
#leave lavage only
#notch2 the fullest data
LAVAGE_df <- LIQUID_DF_final%>%
  filter(!is.na(NOTCH2_NP)) #103 cases
#create a new groupings
LAVAGE_df <- LAVAGE_df %>%
  mutate(
    TYPE_BENIGN2 = if_else(TYPE %in% c("RSS", "BENIGN"),
                           "BENIGN",
                           TYPE)
  )
table(LAVAGE_df$TYPE_BENIGN2, useNA = "a") #now 31 benign
LAVAGE_df <- LAVAGE_df %>%
  mutate(
    TYPE_BENIGN3 = if_else(TYPE_BENIGN2 %in% c("HGSOC", "OTHER"),
                           "OC",
                           TYPE_BENIGN2)
  )
table(LAVAGE_df$TYPE_BENIGN3, useNA = "a") #now 31 benign
DATA <- c("NOTCH2_NP","CTNNB1_NP","DLL1_NP","HES1_NP" )

#OC vs benign+RSS#################################
#HGSOC vs BENIGN DF
OC_BEN_lavage <- LAVAGE_df[c(LAVAGE_df$TYPE_BENIGN3 != "ENDOMETRIAL CANCER"),] 
OC_BEN_lavage$TYPE_BENIGN3 <- relevel(factor(OC_BEN_lavage$TYPE_BENIGN3), ref = "BENIGN")
OC_BEN_lavage <- OC_BEN_lavage[
  complete.cases(
    OC_BEN_lavage[, c(
      "TYPE_BENIGN3",
      "NOTCH2_NP",
      "CTNNB1_NP",
      "DLL1_NP",
      "HES1_NP"
    )]
  ),
]
brglm.model_2 <- glm(OC_BEN_lavage$TYPE_BENIGN3 ~ NOTCH2_NP  + CTNNB1_NP+  DLL1_NP +HES1_NP,
                     data = OC_BEN_lavage,
                     family = binomial("logit"), method = "brglm_fit") #no problems
predicted_probs_2 <- predict.glm(brglm.model_2, type='response') 
pred_data2 <- OC_BEN_lavage[(rownames(OC_BEN_lavage) #remove incoplete rows
                             %in% names(predicted_probs_2)), ]
roc_curve2 <- roc(pred_data2$TYPE_BENIGN3, predicted_probs_2)
AUC2 <- auc(roc_curve2)
coords2 <- coords(roc_curve2,
                  "best", ret=c("accuracy", "sensitivity", "specificity"), transpose = FALSE)
coords2$AUC <- AUC2
coords2
#on the same data
#Duos###################################################
# Generate all pairwise combinations
gene_pairs <- combn(DATA, 2, simplify = FALSE)

# Store results
results_list <- list()
roc_list <- list()
# Loop through gene pairs
for(i in seq_along(gene_pairs)) {
  
  pair <- gene_pairs[[i]]
  
  vars_needed <- c("TYPE_BENIGN3", pair)
  
  # Remove missing values
  model_df <- OC_BEN_lavage[
    complete.cases(OC_BEN_lavage[, vars_needed]),
  ]
  
  # Create formula dynamically
  formula_text <- paste(
    "TYPE_BENIGN3 ~",
    paste(pair, collapse = " + ")
  )
  
  model_formula <- as.formula(formula_text)
  
  # Fit model
  fit <- glm(
    model_formula,
    data = model_df,
    family = binomial("logit"),
    method = "brglm_fit"
  )
  
  # Predict
  probs <- predict(fit, type = "response")
  
  # ROC + AUC
  roc_obj <- roc(model_df$TYPE_BENIGN3, probs)
  
  # Save ROC object
  roc_name <- paste(pair, collapse = "_")
  roc_list[[roc_name]] <- roc_obj
  
  auc_val <- auc(roc_obj)
  
  best_coords <- coords(
    roc_obj,
    "best",
    ret = c("accuracy", "sensitivity", "specificity"),
    transpose = FALSE
  )
  
  # Store results
  results_list[[i]] <- data.frame(
    Gene1 = pair[1],
    Gene2 = pair[2],
    AUC = as.numeric(auc_val),
    Accuracy = best_coords$accuracy,
    Sensitivity = best_coords$sensitivity,
    Specificity = best_coords$specificity
  )
}

# Combine all resultDATA# Combine all results
results_df <- do.call(rbind, results_list)

# Sort by AUC
results_df <- results_df[order(-results_df$AUC), ]

results_df
#Trios###########################################
# Generate all 3-gene combinations
gene_trios <- combn(DATA, 3, simplify = FALSE)

# Store results
results_list2 <- list()
roc_list2 <- list()

# Loop through trios
for(i in seq_along(gene_trios)) {
  
  trio <- gene_trios[[i]]
  
  vars_needed <- c("TYPE_BENIGN3", trio)
  
  # Remove rows with missing values
  model_df <- OC_BEN_lavage[
    complete.cases(OC_BEN_lavage[, vars_needed]),
  ]
  
  # Build formula
  formula_text <- paste(
    "TYPE_BENIGN3 ~",
    paste(trio, collapse = " + ")
  )
  
  model_formula <- as.formula(formula_text)
  
  # Fit model
  fit <- glm(
    model_formula,
    data = model_df,
    family = binomial("logit"),
    method = "brglm_fit"
  )
  
  # Predictions
  probs <- predict(fit, type = "response")
  
  # ROC/AUC
  roc_obj <- roc(model_df$TYPE_BENIGN3, probs)
  
  # Save ROC object
  roc_name <- paste(trio, collapse = "_")
  roc_list2[[roc_name]] <- roc_obj
  
  auc_val <- auc(roc_obj)
  
  best_coords <- coords(
    roc_obj,
    "best",
    ret = c("accuracy", "sensitivity", "specificity"),
    transpose = FALSE
  )
  
  # Store results
  results_list2[[i]] <- data.frame(
    Gene1 = trio[1],
    Gene2 = trio[2],
    Gene3 = trio[3],
    AUC = as.numeric(auc_val),
    Accuracy = best_coords$accuracy,
    Sensitivity = best_coords$sensitivity,
    Specificity = best_coords$specificity
  )
}

# Combine all results
results_df2 <- do.call(rbind, results_list2)

# Sort by AUC
results_df2 <- results_df2[order(-results_df2$AUC), ]

results_df2

##ROC PLOT OC BEN #########################
roc_plotOC_BEN3 <- function() {
  par(pty = "s") #sets square
  plot.roc(roc_list2[["NOTCH2_NP_CTNNB1_NP_DLL1_NP"]], print.auc = F, col = "#dcbeff",
           cex.main=0.8, 
           main ="Uterine lavage biomarkers OC vs Non-cancer samples",
           # xlab = "1 - Specifiškumas", #lithuanian version
           # ylab = "Jautrumas", 
           legacy.axes = T) #title
  lines(roc_list2[["NOTCH2_NP_CTNNB1_NP_HES1_NP"]], col = "#911eb4", lwd =2) 
  lines(roc_list2[["NOTCH2_NP_DLL1_NP_HES1_NP"]], col ="#ffd8b1", lwd =2) 
  lines(roc_list2[["CTNNB1_NP_DLL1_NP_HES1_NP"]], col = "#42d4f4", lwd =2) 
  
  
  lines(roc_list[["NOTCH2_NP_CTNNB1_NP"]], col = "#3cb44b", lwd =2) 
  lines(roc_list[["NOTCH2_NP_DLL1_NP"]], col ="#e6194B", lwd =2) 
  lines(roc_list[["NOTCH2_NP_HES1_NP"]], col = "#4363d8", lwd =2) 
  lines(roc_list[["CTNNB1_NP_DLL1_NP"]], col = "#f58231", lwd =2) 
  lines(roc_list[["CTNNB1_NP_HES1_NP"]], col ="#a9a9a9", lwd =2) 
  lines(roc_list[["DLL1_NP_HES1_NP"]], col = "#800000", lwd =2) 
  
  lines(roc_curve2, col = "black", lwd =2, lty = 2) 
  
  legend("bottomright", legend = c( expression(italic("NOTCH2 + CTNNB1 + DLL1")),
                                    expression(italic("NOTCH2 + CTNNB1 + HES1")),
                                    expression(italic("NOTCH2 + DLL1 + HES1")), 
                                    expression(italic("CTNNB1 + DLL1 + HES1")),
                                    
                                    expression(italic("NOTCH2 + CTNNB1")),
                                    expression(italic("NOTCH2 + DLL1")), 
                                    expression(italic("NOTCH2 + HES1")),
                                    expression(italic("CTNNB1 + DLL1")),
                                    expression(italic("CTNNB1 + HES1")), 
                                    expression(italic("DLL1 + HES1")),
                                    
                                    expression(italic("DLL1 + HES1 + CTNNB1 + NOTCH2"))
                                    
                                    
  ),
  
  col = c("#dcbeff", "#911eb4", "#ffd8b1", "#42d4f4",
          "#3cb44b", "#e6194B","#4363d8", "#f58231", "#a9a9a9", "#800000", "black"  ), lty = 1, 
  cex = 0.7, lwd =3)
}
#plot
roc_plotOC_BEN3()

# Save the plot as a PNG file
png("C:/Users/Ieva/rprojects/outputs_all/LIQUID/non-cancer-oc-combs200260507.png",
    width = 15, height = 15, res = 300, units = "cm")
roc_plotOC_BEN3()
# mtext(
#   "A",
#   side = 3,
#   adj = -0.2,
#   line = 1,
#   cex = 1.2,
#   font = 2
# )
dev.off()
##GT TABLE##########################################
#make coords
coords2
results_df
results_df2

## ---- Single gene model ----
single_df <- data.frame(
  Model = "NOTCH2_NP + CTNNB1_NP + DLL1_NP + HES1_NP",
  Genes = 4,
  AUC = as.numeric(coords2$AUC),
  Accuracy = coords2$accuracy,
  Sensitivity = coords2$sensitivity,
  Specificity = coords2$specificity
)
## ---- Pair models ----
pairs_df <- results_df %>%
  mutate(
    Model = paste(Gene1, Gene2, sep = " + "),
    Genes = 2
  ) %>%
  select(Model, Genes, AUC, Accuracy, Sensitivity, Specificity)

## ---- Trio models ----
trios_df <- results_df2 %>%
  mutate(
    Model = paste(Gene1, Gene2, Gene3, sep = " + "),
    Genes = 3
  ) %>%
  select(Model, Genes, AUC, Accuracy, Sensitivity, Specificity)

## ---- Combine everything ----
final_results <- bind_rows(
  single_df,
  pairs_df,
  trios_df
) %>%
  mutate(Model = gsub("_NP", "", Model)) %>%
  arrange(desc(AUC))

## ---- GT table ----
gt <- final_results %>%
  gt() %>%
  fmt_number(
    columns = c(AUC, Accuracy, Sensitivity, Specificity),
    decimals = 3
  ) %>%
  cols_label(
    Model = "Gene Combination",
    Genes = "N Genes",
    AUC = "AUC",
    Accuracy = "Accuracy",
    Sensitivity = "Sensitivity",
    Specificity = "Specificity"
  ) %>%
  tab_header(
    title = "OC vs Non-cancer Logistic Regression Models"
  )%>%
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_body(columns = vars(Model))
  )

gt
#SAVE gt
gtsave(gt,vwidth = 10000,   
       filename = "C:/Users/Ieva/rprojects/outputs_all/LIQUID/non-cancer-oc-combstable200260507.png")
#import images
roc_imageoc_be_combs  <- image_read("C:/Users/Ieva/rprojects/outputs_all/LIQUID/non-cancer-oc-combs200260507.png")
table_imageoc_be_combs <- image_read("C:/Users/Ieva/rprojects/outputs_all/LIQUID/non-cancer-oc-combstable200260507.png")

# resize table to match ROC image width
table_imageoc_be_combs <- image_resize(table_imageoc_be_combs,
                                      paste0(image_info(roc_imageoc_be_combs)$width, "x"))
# combine vertically
combinedoc_be_combs <- image_append(c(roc_imageoc_be_combs, table_imageoc_be_combs), stack = TRUE)

# save
image_write(combinedoc_be_combs,
            "C:/Users/Ieva/rprojects/outputs_all/LIQUID/non-cancer-oc-combinations200260507.png")

#EC vs non-cancer#################################################
#HGSOC vs BENIGN DF
EC_BEN_lavage <- LAVAGE_df[c(LAVAGE_df$TYPE_BENIGN3 != "OC"),] 
EC_BEN_lavage$TYPE_BENIGN3 <- factor(EC_BEN_lavage$TYPE_BENIGN3)
EC_BEN_lavage$TYPE_BENIGN3 <- relevel(factor(EC_BEN_lavage$TYPE_BENIGN3), ref = "BENIGN")
EC_BEN_lavage <- EC_BEN_lavage[
  complete.cases(
    EC_BEN_lavage[, c(
      "TYPE_BENIGN3",
      "NOTCH2_NP",
      "CTNNB1_NP",
      "DLL1_NP",
      "HES1_NP"
    )]
  ),
]
brglm.model_3<- glm(EC_BEN_lavage$TYPE_BENIGN3 ~ NOTCH2_NP  + CTNNB1_NP+  DLL1_NP +HES1_NP,
                    data = EC_BEN_lavage,
                    family = binomial("logit"), method = "brglm_fit") #no problems
predicted_probs_3 <- predict.glm(brglm.model_3, type='response') 
pred_data3 <- EC_BEN_lavage[(rownames(EC_BEN_lavage) #remove incoplete rows
                             %in% names(predicted_probs_3)), ]
roc_curve3 <- roc(pred_data3$TYPE_BENIGN3, predicted_probs_3)
AUC3 <- auc(roc_curve3)
coords3 <- coords(roc_curve3,
                  "best", ret=c("accuracy", "sensitivity", "specificity"), transpose = FALSE)
coords3$AUC <- AUC3
coords3
#on the same data
##Duos###################################################
# Generate all pairwise combinations
gene_pairs <- combn(DATA, 2, simplify = FALSE)

# Store results
results_list_EC <- list()
roc_list_ec <- list()
# Loop through gene pairs
for(i in seq_along(gene_pairs)) {
  
  pair <- gene_pairs[[i]]
  
  vars_needed <- c("TYPE_BENIGN3", pair)
  
  # Remove missing values
  model_df <- EC_BEN_lavage[
    complete.cases(EC_BEN_lavage[, vars_needed]),
  ]
  
  # Create formula dynamically
  formula_text <- paste(
    "TYPE_BENIGN3 ~",
    paste(pair, collapse = " + ")
  )
  
  model_formula <- as.formula(formula_text)
  
  # Fit model
  fit <- glm(
    model_formula,
    data = model_df,
    family = binomial("logit"),
    method = "brglm_fit"
  )
  
  # Predict
  probs <- predict(fit, type = "response")
  
  # ROC + AUC
  roc_obj <- roc(model_df$TYPE_BENIGN3, probs)
  
  # Save ROC object
  roc_name <- paste(pair, collapse = "_")
  roc_list_ec[[roc_name]] <- roc_obj
  
  auc_val <- auc(roc_obj)
  
  best_coords <- coords(
    roc_obj,
    "best",
    ret = c("accuracy", "sensitivity", "specificity"),
    transpose = FALSE
  )
  
  # Store results
  results_list_EC[[i]] <- data.frame(
    Gene1 = pair[1],
    Gene2 = pair[2],
    AUC = as.numeric(auc_val),
    Accuracy = best_coords$accuracy,
    Sensitivity = best_coords$sensitivity,
    Specificity = best_coords$specificity
  )
}

# Combine all resultDATA# Combine all results
results_df_EC <- do.call(rbind, results_list_EC)

# Sort by AUC
results_df_EC <- results_df_EC[order(-results_df_EC$AUC), ]

results_df_EC

##Trios###########################################
# Generate all 3-gene combinations
gene_trios <- combn(DATA, 3, simplify = FALSE)

# Store results
results_list2EC <- list()
roc_list2EC <- list()

# Loop through trios
for(i in seq_along(gene_trios)) {
  
  trio <- gene_trios[[i]]
  
  vars_needed <- c("TYPE_BENIGN3", trio)
  
  # Remove rows with missing values
  model_df <- EC_BEN_lavage[
    complete.cases(EC_BEN_lavage[, vars_needed]),
  ]
  
  # Build formula
  formula_text <- paste(
    "TYPE_BENIGN3 ~",
    paste(trio, collapse = " + ")
  )
  
  model_formula <- as.formula(formula_text)
  
  # Fit model
  fit <- glm(
    model_formula,
    data = model_df,
    family = binomial("logit"),
    method = "brglm_fit"
  )
  
  # Predictions
  probs <- predict(fit, type = "response")
  
  # ROC/AUC
  roc_obj <- roc(model_df$TYPE_BENIGN3, probs)
  
  # Save ROC object
  roc_name <- paste(trio, collapse = "_")
  roc_list2EC[[roc_name]] <- roc_obj
  
  auc_val <- auc(roc_obj)
  
  best_coords <- coords(
    roc_obj,
    "best",
    ret = c("accuracy", "sensitivity", "specificity"),
    transpose = FALSE
  )
  
  # Store results
  results_list2EC[[i]] <- data.frame(
    Gene1 = trio[1],
    Gene2 = trio[2],
    Gene3 = trio[3],
    AUC = as.numeric(auc_val),
    Accuracy = best_coords$accuracy,
    Sensitivity = best_coords$sensitivity,
    Specificity = best_coords$specificity
  )
}

# Combine all results
results_df2EC <- do.call(rbind, results_list2EC)

# Sort by AUC
results_df2EC <- results_df2EC[order(-results_df2EC$AUC), ]

results_df2EC

##ROC PLOT OC BEN #########################
roc_plotEC_BEN3 <- function() {
  par(pty = "s") #sets square
  plot.roc(roc_list2EC[["NOTCH2_NP_CTNNB1_NP_DLL1_NP"]], print.auc = F, col = "#dcbeff",
           cex.main=0.8, 
           main ="Uterine lavage biomarkers EC vs Non-cancer samples",
           # xlab = "1 - Specifiškumas", #lithuanian version
           # ylab = "Jautrumas", 
           legacy.axes = T) #title
  lines(roc_list2EC[["NOTCH2_NP_CTNNB1_NP_HES1_NP"]], col = "#911eb4", lwd =2) 
  lines(roc_list2EC[["NOTCH2_NP_DLL1_NP_HES1_NP"]], col ="#ffd8b1", lwd =2) 
  lines(roc_list2EC[["CTNNB1_NP_DLL1_NP_HES1_NP"]], col = "#42d4f4", lwd =2) 
  
  
  lines(roc_list_ec[["NOTCH2_NP_CTNNB1_NP"]], col = "#3cb44b", lwd =2) 
  lines(roc_list_ec[["NOTCH2_NP_DLL1_NP"]], col ="#e6194B", lwd =2) 
  lines(roc_list_ec[["NOTCH2_NP_HES1_NP"]], col = "#4363d8", lwd =2) 
  lines(roc_list_ec[["CTNNB1_NP_DLL1_NP"]], col = "#f58231", lwd =2) 
  lines(roc_list_ec[["CTNNB1_NP_HES1_NP"]], col ="#a9a9a9", lwd =2) 
  lines(roc_list_ec[["DLL1_NP_HES1_NP"]], col = "#800000", lwd =2) 
  
  lines(roc_curve3, col = "black", lwd =2, lty = 2) 
  
  legend("bottomright", legend = c( expression(italic("NOTCH2 + CTNNB1 + DLL1")),
                                    expression(italic("NOTCH2 + CTNNB1 + HES1")),
                                    expression(italic("NOTCH2 + DLL1 + HES1")), 
                                    expression(italic("CTNNB1 + DLL1 + HES1")),
                                    
                                    expression(italic("NOTCH2 + CTNNB1")),
                                    expression(italic("NOTCH2 + DLL1")), 
                                    expression(italic("NOTCH2 + HES1")),
                                    expression(italic("CTNNB1 + DLL1")),
                                    expression(italic("CTNNB1 + HES1")), 
                                    expression(italic("DLL1 + HES1")),
                                    
                                    expression(italic("DLL1 + HES1 + CTNNB1 + NOTCH2"))
                                    
                                    
  ),
  
  col = c("#dcbeff", "#911eb4", "#ffd8b1", "#42d4f4",
          "#3cb44b", "#e6194B","#4363d8", "#f58231", "#a9a9a9", "#800000", "black"  ), lty = 1, 
  cex = 0.7, lwd =3)
}
#plot
roc_plotEC_BEN3()

# Save the plot as a PNG file
png("C:/Users/Ieva/rprojects/outputs_all/LIQUID/non-cancer-ec-combs200260508.png",
    width = 15, height = 15, res = 300, units = "cm")
roc_plotEC_BEN3()
# mtext(
#   "A",
#   side = 3,
#   adj = -0.2,
#   line = 1,
#   cex = 1.2,
#   font = 2
# )
dev.off()

#make coords
coords3
results_df_EC
results_df2EC

## GT TABLE###################################
## ---- Single gene model ----
single_dfEC <- data.frame(
  Model = "NOTCH2_NP + CTNNB1_NP + DLL1_NP + HES1_NP",
  Genes = 4,
  AUC = as.numeric(coords3$AUC),
  Accuracy = coords3$accuracy,
  Sensitivity = coords3$sensitivity,
  Specificity = coords3$specificity
)
## ---- Pair models ----
pairs_dfEC <- results_df_EC %>%
  mutate(
    Model = paste(Gene1, Gene2, sep = " + "),
    Genes = 2
  ) %>%
  select(Model, Genes, AUC, Accuracy, Sensitivity, Specificity)

## ---- Trio models ----
trios_dfEC <- results_df2EC %>%
  mutate(
    Model = paste(Gene1, Gene2, Gene3, sep = " + "),
    Genes = 3
  ) %>%
  select(Model, Genes, AUC, Accuracy, Sensitivity, Specificity)

## ---- Combine everything ----
final_resultsEC <- bind_rows(
  single_dfEC,
  pairs_dfEC,
  trios_dfEC
) %>%
  mutate(Model = gsub("_NP", "", Model)) %>%
  arrange(desc(AUC))

## ---- GT table ----
gt_EC <- final_resultsEC %>%
  gt() %>%
  fmt_number(
    columns = c(AUC, Accuracy, Sensitivity, Specificity),
    decimals = 3
  ) %>%
  cols_label(
    Model = "Gene Combination",
    Genes = "N Genes",
    AUC = "AUC",
    Accuracy = "Accuracy",
    Sensitivity = "Sensitivity",
    Specificity = "Specificity"
  ) %>%
  tab_header(
    title = "EC vs Non-cancer Logistic Regression Models"
  )%>%
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_body(columns = vars(Model))
  )

gt_EC
#SAVE gt
gtsave(gt_EC,vwidth = 10000,   
       filename = "C:/Users/Ieva/rprojects/outputs_all/LIQUID/non-cancer-ec-combstable200260509.png")
#import images
roc_imageoc_be_combs  <- image_read("C:/Users/Ieva/rprojects/outputs_all/LIQUID/non-cancer-ec-combs200260508.png")
table_imageoc_be_combs <- image_read("C:/Users/Ieva/rprojects/outputs_all/LIQUID/non-cancer-ec-combstable200260509.png")

# resize table to match ROC image width
table_imageoc_be_combs <- image_resize(table_imageoc_be_combs,
                                       paste0(image_info(roc_imageoc_be_combs)$width, "x"))
# combine vertically
combinedoc_be_combs <- image_append(c(roc_imageoc_be_combs, table_imageoc_be_combs), stack = TRUE)

# save
image_write(combinedoc_be_combs,
            "C:/Users/Ieva/rprojects/outputs_all/LIQUID/non-cancer-ec-combinations200260508.png")
#EC vs OC#################################################
#HGSOC vs BENIGN DF
EC_OC_lavage <- LAVAGE_df[c(LAVAGE_df$TYPE_BENIGN3 != "BENIGN"),] 
EC_OC_lavage$TYPE_BENIGN3 <- factor(EC_OC_lavage$TYPE_BENIGN3)
EC_OC_lavage$TYPE_BENIGN3 <- relevel(factor(EC_OC_lavage$TYPE_BENIGN3), ref = "ENDOMETRIAL CANCER")
EC_OC_lavage <- EC_OC_lavage[
  complete.cases(
    EC_OC_lavage[, c(
      "TYPE_BENIGN3",
      "NOTCH2_NP",
      "CTNNB1_NP",
      "DLL1_NP",
      "HES1_NP"
    )]
  ),
]
brglm.model_4<- glm(EC_OC_lavage$TYPE_BENIGN3 ~ NOTCH2_NP  + CTNNB1_NP+  DLL1_NP +HES1_NP,
                    data = EC_OC_lavage,
                    family = binomial("logit"), method = "brglm_fit") #no problems
predicted_probs_4 <- predict.glm(brglm.model_4, type='response') 
pred_data4 <- EC_OC_lavage[(rownames(EC_OC_lavage) #remove incoplete rows
                            %in% names(predicted_probs_4)), ]
roc_curve4 <- roc(pred_data4$TYPE_BENIGN3, predicted_probs_4)
AUC4 <- auc(roc_curve4)
coords4 <- coords(roc_curve4,
                  "best", ret=c("accuracy", "sensitivity", "specificity"),
                  transpose = FALSE)
coords4$AUC <- AUC4
coords4
#on the same data
##Duos###################################################
# Generate all pairwise combinations
gene_pairs <- combn(DATA, 2, simplify = FALSE)

# Store results
results_list_EC_OC <- list()
roc_list_EC_OC <- list()
# Loop through gene pairs
for(i in seq_along(gene_pairs)) {
  
  pair <- gene_pairs[[i]]
  
  vars_needed <- c("TYPE_BENIGN3", pair)
  
  # Remove missing values
  model_df <- EC_OC_lavage[
    complete.cases(EC_OC_lavage[, vars_needed]),
  ]
  
  # Create formula dynamically
  formula_text <- paste(
    "TYPE_BENIGN3 ~",
    paste(pair, collapse = " + ")
  )
  
  model_formula <- as.formula(formula_text)
  
  # Fit model
  fit <- glm(
    model_formula,
    data = model_df,
    family = binomial("logit"),
    method = "brglm_fit"
  )
  
  # Predict
  probs <- predict(fit, type = "response")
  
  # ROC + AUC
  roc_obj <- roc(model_df$TYPE_BENIGN3, probs)
  
  # Save ROC object
  roc_name <- paste(pair, collapse = "_")
  roc_list_EC_OC[[roc_name]] <- roc_obj
  
  auc_val <- auc(roc_obj)
  
  best_coords <- coords(
    roc_obj,
    "best",
    ret = c("accuracy", "sensitivity", "specificity"),
    transpose = FALSE
  )
  
  # Store results
  results_list_EC_OC[[i]] <- data.frame(
    Gene1 = pair[1],
    Gene2 = pair[2],
    AUC = as.numeric(auc_val),
    Accuracy = best_coords$accuracy,
    Sensitivity = best_coords$sensitivity,
    Specificity = best_coords$specificity
  )
}

# Combine all resultDATA# Combine all results
results_df_EC_OC <- do.call(rbind, results_list_EC_OC)

# Sort by AUC
results_df_EC_OC <- results_df_EC_OC[order(-results_df_EC_OC$AUC), ]

results_df_EC_OC

##Trios###########################################
# Generate all 3-gene combinations
gene_trios <- combn(DATA, 3, simplify = FALSE)

# Store results
results_list2EC_OC <- list()
roc_list2EC_OC <- list()

# Loop through trios
for(i in seq_along(gene_trios)) {
  
  trio <- gene_trios[[i]]
  
  vars_needed <- c("TYPE_BENIGN3", trio)
  
  # Remove rows with missing values
  model_df <- EC_OC_lavage[
    complete.cases(EC_OC_lavage[, vars_needed]),
  ]
  
  # Build formula
  formula_text <- paste(
    "TYPE_BENIGN3 ~",
    paste(trio, collapse = " + ")
  )
  
  model_formula <- as.formula(formula_text)
  
  # Fit model
  fit <- glm(
    model_formula,
    data = model_df,
    family = binomial("logit"),
    method = "brglm_fit"
  )
  
  # Predictions
  probs <- predict(fit, type = "response")
  
  # ROC/AUC
  roc_obj <- roc(model_df$TYPE_BENIGN3, probs)
  
  # Save ROC object
  roc_name <- paste(trio, collapse = "_")
  roc_list2EC_OC[[roc_name]] <- roc_obj
  
  auc_val <- auc(roc_obj)
  
  best_coords <- coords(
    roc_obj,
    "best",
    ret = c("accuracy", "sensitivity", "specificity"),
    transpose = FALSE
  )
  
  # Store results
  results_list2EC_OC[[i]] <- data.frame(
    Gene1 = trio[1],
    Gene2 = trio[2],
    Gene3 = trio[3],
    AUC = as.numeric(auc_val),
    Accuracy = best_coords$accuracy,
    Sensitivity = best_coords$sensitivity,
    Specificity = best_coords$specificity
  )
}

# Combine all results
results_df2EC_OC <- do.call(rbind, results_list2EC)

# Sort by AUC
results_df2EC_OC <- results_df2EC_OC[order(-results_df2EC_OC$AUC), ]

results_df2EC_OC

##ROC PLOT OC BEN #########################
roc_plotEC_OC <- function() {
  par(pty = "s") #sets square
  plot.roc(roc_list2EC_OC[["NOTCH2_NP_CTNNB1_NP_DLL1_NP"]], print.auc = F, col = "#dcbeff",
           cex.main=0.8, 
           main ="Uterine lavage biomarkers EC vs OC samples",
           # xlab = "1 - Specifiškumas", #lithuanian version
           # ylab = "Jautrumas", 
           legacy.axes = T) #title
  lines(roc_list2EC_OC[["NOTCH2_NP_CTNNB1_NP_HES1_NP"]], col = "#911eb4", lwd =2) 
  lines(roc_list2EC_OC[["NOTCH2_NP_DLL1_NP_HES1_NP"]], col ="#ffd8b1", lwd =2) 
  lines(roc_list2EC_OC[["CTNNB1_NP_DLL1_NP_HES1_NP"]], col = "#42d4f4", lwd =2) 
  
  
  lines(roc_list_EC_OC[["NOTCH2_NP_CTNNB1_NP"]], col = "#3cb44b", lwd =2) 
  lines(roc_list_EC_OC[["NOTCH2_NP_DLL1_NP"]], col ="#e6194B", lwd =2) 
  lines(roc_list_EC_OC[["NOTCH2_NP_HES1_NP"]], col = "#4363d8", lwd =2) 
  lines(roc_list_EC_OC[["CTNNB1_NP_DLL1_NP"]], col = "#f58231", lwd =2) 
  lines(roc_list_EC_OC[["CTNNB1_NP_HES1_NP"]], col ="#a9a9a9", lwd =2) 
  lines(roc_list_EC_OC[["DLL1_NP_HES1_NP"]], col = "#800000", lwd =2) 
  
  lines(roc_curve4, col = "black", lwd =2, lty = 2) 
  
  legend("bottomright", legend = c( expression(italic("NOTCH2 + CTNNB1 + DLL1")),
                                    expression(italic("NOTCH2 + CTNNB1 + HES1")),
                                    expression(italic("NOTCH2 + DLL1 + HES1")), 
                                    expression(italic("CTNNB1 + DLL1 + HES1")),
                                    
                                    expression(italic("NOTCH2 + CTNNB1")),
                                    expression(italic("NOTCH2 + DLL1")), 
                                    expression(italic("NOTCH2 + HES1")),
                                    expression(italic("CTNNB1 + DLL1")),
                                    expression(italic("CTNNB1 + HES1")), 
                                    expression(italic("DLL1 + HES1")),
                                    
                                    expression(italic("DLL1 + HES1 + CTNNB1 + NOTCH2"))
                                    
                                    
  ),
  
  col = c("#dcbeff", "#911eb4", "#ffd8b1", "#42d4f4",
          "#3cb44b", "#e6194B","#4363d8", "#f58231", "#a9a9a9", "#800000", "black"  ), lty = 1, 
  cex = 0.7, lwd =3)
}
#plot
roc_plotEC_OC()

# Save the plot as a PNG file
png("C:/Users/Ieva/rprojects/outputs_all/LIQUID/oc-ec-combs200260508.png",
    width = 15, height = 15, res = 300, units = "cm")
roc_plotEC_OC()
# mtext(
#   "A",
#   side = 3,
#   adj = -0.2,
#   line = 1,
#   cex = 1.2,
#   font = 2
# )
dev.off()

#make coords
coords4
results_df_EC_OC
results_df2EC_OC

## GT TABLE###################################
## ---- Single gene model ----
single_dfEC_OC <- data.frame(
  Model = "NOTCH2_NP + CTNNB1_NP + DLL1_NP + HES1_NP",
  Genes = 4,
  AUC = as.numeric(coords4$AUC),
  Accuracy = coords4$accuracy,
  Sensitivity = coords4$sensitivity,
  Specificity = coords4$specificity
)
## ---- Pair models ----
pairs_dfEC_OC <- results_df_EC_OC %>%
  mutate(
    Model = paste(Gene1, Gene2, sep = " + "),
    Genes = 2
  ) %>%
  select(Model, Genes, AUC, Accuracy, Sensitivity, Specificity)

## ---- Trio models ----
trios_dfEC_OC <- results_df2EC_OC %>%
  mutate(
    Model = paste(Gene1, Gene2, Gene3, sep = " + "),
    Genes = 3
  ) %>%
  select(Model, Genes, AUC, Accuracy, Sensitivity, Specificity)

## ---- Combine everything ----
final_resultsEC_OC <- bind_rows(
  single_dfEC_OC,
  pairs_dfEC_OC,
  trios_dfEC_OC
) %>%
  mutate(Model = gsub("_NP", "", Model)) %>%
  arrange(desc(AUC))

## ---- GT table ----
gt_EC_OC <- final_resultsEC_OC %>%
  gt() %>%
  fmt_number(
    columns = c(AUC, Accuracy, Sensitivity, Specificity),
    decimals = 3
  ) %>%
  cols_label(
    Model = "Gene Combination",
    Genes = "N Genes",
    AUC = "AUC",
    Accuracy = "Accuracy",
    Sensitivity = "Sensitivity",
    Specificity = "Specificity"
  ) %>%
  tab_header(
    title = "EC vs OC Logistic Regression Models"
  )%>%
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_body(columns = vars(Model))
  )

gt_EC_OC
#SAVE gt
gtsave(gt_EC_OC,vwidth = 10000,   
       filename = "C:/Users/Ieva/rprojects/outputs_all/LIQUID/oceccombstable200260509.png")
#import images
roc_imageoc_be_combs  <- image_read("C:/Users/Ieva/rprojects/outputs_all/LIQUID/oc-ec-combs200260508.png")
table_imageoc_be_combs <- image_read("C:/Users/Ieva/rprojects/outputs_all/LIQUID/oceccombstable200260509.png")

# resize table to match ROC image width
table_imageoc_be_combs <- image_resize(table_imageoc_be_combs,
                                       paste0(image_info(roc_imageoc_be_combs)$width, "x"))
# combine vertically
combinedoc_be_combs <- image_append(c(roc_imageoc_be_combs, table_imageoc_be_combs), stack = TRUE)

# save
image_write(combinedoc_be_combs,
            "C:/Users/Ieva/rprojects/outputs_all/LIQUID/oc-ec-combinations200260508.png")
