#KN-liquid - urine 250916
Sys.setenv(LANG = "en")
#libraries
library(tidyverse)
library(pROC)
library(glmnet)
library(gtsummary)
library(gt)
library(grid)
library(brglm2)
library(reshape2)
library(rstatix) 
library(ggprism)
library(gridExtra)
library(png)
library(htmlwidgets)
library(webshot)
library(magick)

#add urine data ###################################
KN_data_full <- readRDS( "C:/Users/Ieva/rprojects/OTHER DATA/KN_LIQUID/KN_FULL_LIQUID_R250916.RDS")
urine_liquid <- KN_data_full[51:65, c(1:19, 38:41)] #leave only URINE DATA
urine_np <- c("NOTCH2_URINE", "CTNNB1_URINE" , "DLL1_URINE" ,  "HES1_URINE") #names

#SHAPIRO TEST##########################################
#2 GROUPS # except DLL1, too litle data
shapiro_results1 <- urine_liquid[, c(5, 20, 21, 23)] %>%
  pivot_longer(cols = -Grupė_Ieva , names_to = "gene", values_to = "value") %>%
  group_by(Grupė_Ieva , gene) %>%
  summarise(p_value = shapiro.test(value)$p.value, .groups = "drop") 
shapiro_results1 #normal
#VARTEST#####################################################
var_results_2 <- urine_liquid[,  c(5, 20, 21, 23)] %>%
  pivot_longer(cols = -Grupė_Ieva , names_to = "variable", values_to = "value") %>%
  group_by(variable) %>%
  summarise(
    p_value = var.test(value[Grupė_Ieva  == unique(Grupė_Ieva )[1]], 
                       value[Grupė_Ieva  == unique(Grupė_Ieva )[2]])$p.value,
    .groups = "drop"
  ) 
var_results_2 #all equal

#PAIRWISE STJUDENTS T TEST ###################################
#melt table for expression
Group2_table <- melt(urine_liquid[, c(5, 20, 21, 22, 23)], id.vars="Grupė_Ieva",  measure.vars=urine_np)

#stjundents test (normal, equal variances)
t.test_2groups <- Group2_table %>%
  group_by(variable) %>%
  t_test(value ~ Grupė_Ieva,
         p.adjust.method = "BH", 
         var.equal = TRUE, #stjudents
         paired = FALSE, 
         #detailed=TRUE 
  )
t.test_2groups
#BOXPLOT##################################
#rename to lt 
custom_colors <- c("HGSOC" = "deeppink","Other" = "lightpink") 
OC_plot <- ggplot(Group2_table, aes(x=Grupė_Ieva , y=value, fill = variable)) +
  geom_boxplot( outlier.shape = NA , alpha=0.3, aes(fill = Grupė_Ieva )) +
  geom_jitter(aes(color = Grupė_Ieva ), size=1, alpha=0.5) +
  ylab(label = expression("Santykinė genų raiška šlapime, normalizuota pagal  " * italic("GAPDH"))) + 
  facet_wrap(.~ variable, nrow = 2, scales = "free") +
  #add_pvalue(each.vs.ref_sig, label = "p.adj") + #pvalue
  theme_minimal()+
  theme(
    strip.text.x = element_text(
      size = 12, face = "bold.italic"
    ),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5))+
  labs(x=NULL)+
  stat_boxplot(geom ='errorbar')+
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  scale_y_continuous(labels = function(x) 
    gsub("-", "\u2212", as.character(x))) #add long "-" signs

OC_plot

# Save the plot as a PNG file
png("C:/Users/Ieva/rprojects/outputs_all/URINE_boxplot20250916.png",
    width = 1000, height = 1100, res = 200)
OC_plot
dev.off()

#CLINICAL FEATURES###########################################
#chek clinical features
table(urine_liquid$Grade2, urine_liquid$NOTCH2_URINE) #only 1 G1
table(urine_liquid$Stage2, urine_liquid$NOTCH2_URINE) #only 2 STAGE 1&2
table(urine_liquid$CA125_f, urine_liquid$NOTCH2_URINE) #all >35
#not enough cases to compare further

#CORRELATION WITH AGE#########################
age_table <- urine_liquid[, colnames(urine_liquid) %in% c(urine_np, "Amžius")]
#normalcy
sapply(age_table, function(x) shapiro.test(x)$p.value) #cttnb1 not normal
#pearson corrwlation for all
results_n <- lapply(age_table[, colnames(age_table) %in% urine_np], 
                    function(x) cor.test(x, age_table$Amžius, method = "pearson"))
results_n
#spearman for ctnnb1
cor.test(age_table$Amžius,age_table$CTNNB1_URINE, method = "spearman")
#plot
age_table_l <- age_table %>%
  pivot_longer(cols = -Amžius, names_to = "variable", values_to = "value")

age_plot <- ggplot(age_table_l, aes(x = Amžius, y = value)) +
  geom_point(alpha = 0.6) +
  facet_wrap(~ variable, scales = "free_y") +
  geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") 
age_plot

#add p values
results_df <- bind_rows(
  lapply(names(results_n), function(nm) {
    out <- broom::tidy(results_n[[nm]])
    out$variable <- nm
    out
  })
)

results_df <- results_df %>%
  mutate(
    label = paste0("rho = ", round(estimate, 2),
                   ", p = ", signif(p.value, 3))
  )
#plot
age_plot_p <- age_plot +
  geom_text(
    data = results_df,
    aes(x = max(age_table_l$Amžius, na.rm = TRUE),
        y = Inf,
        label = label),
    inherit.aes = FALSE,
    hjust = 1.1, vjust = 1.5,
    size = 3
  )
age_plot_p
# Save the plot as a PNG file
png("C:/Users/Ieva/rprojects/outputs_all/age_plot_urine20250916.png",
    width = 1000, height = 1100, res = 200)
age_plot_p
dev.off()
