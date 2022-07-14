library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(tibble)
library(readxl)
library(survival)
library(ggfortify)
library(survminer)
library(rlang)
library(janitor)

# IMPORT
MutationData <- read_excel("1- Raw data/Raw mutations data/ovary_mutations_final_27-08.xlsx")
MutationData <-  MutationData %>%
  mutate(change = 1) %>%
mutate(mut = "_mut") %>%
  select(-3,-4) %>% distinct() %>%
  unite(col = Gene,Gene,mut,sep = "") %>%
  pivot_wider(names_from = Gene,values_from = change,values_fill = 0)
clinical_data <- read_excel("1- Raw data/Raw clinical data/Ov_car_final_clinical.xlsx",na = "NA",
                            col_types = c("text","skip","skip","skip","skip","skip",
                                          "skip","text","numeric","text","numeric"))
clinical_data <- clean_names(clinical_data)

full_data_m <- MutationData %>%
  left_join(clinical_data,by = c("SAMPLE_ID"="sample_id"))
# Add functions
full_data_m$Disease_Free <- ifelse(full_data_m$disease_free_status == "DiseaseFree",1,0)
full_data_m$Overall_Free <- ifelse(full_data_m$overall_survival_status == "alive",1,0)
genes <- str_remove(names(MutationData[,2:16]),"_mut")
result_table <- data.frame(genes = genes)
result_table$Disease_free_survival_mut <- NA
result_table$Overall_survival_mut <- NA

# LOOP
for (i in 1:length(genes)) {
  
  gene <- genes[i]
    
  full_data_m_gene <- full_data_m %>%
    select( 1,1+i,17:22)
  names(full_data_m_gene)[2] <- "gene_group"
  km_trt_fit5 <- survfit(Surv(disease_free_months,Disease_Free ) ~ gene_group,data=full_data_m_gene)
  pvalue_df <- surv_pvalue(km_trt_fit5)$pval
  result_table$Disease_free_survival_mut[i] <- pvalue_df
  km_trt_fit8 <- survfit(Surv(overall_survival_months,Overall_Free) ~ gene_group,data=full_data_m_gene)
  pvalue_o <- surv_pvalue(km_trt_fit8 )$pval
  result_table$Overall_survival_mut[i] <- pvalue_o
  
  if(pvalue_df < 0.05){
    
    SurvPlot_df <- ggsurvplot(
      km_trt_fit5,  data = full_data_m_gene ,
      size = 1,                # change line size
      palette = c( "#E7B800","#2E9FDF"),# custom color palettes
      conf.int = TRUE,         # Add confidence interval
      pval = T,             # Add p-value
      risk.table = TRUE,       # Add risk table
      risk.table.col = "strata",# Risk table color by groups
      legend.labs =
        c("No alteration","Alteration"),   # Change legend labels
      risk.table.height = 0.25,# Useful to change when you have multiple groups
      ggtheme = theme_bw(),     # Change ggplot2 theme
      xlab = "Disease Free Survival (Months)", ylab = "Disease Free Survival",title = paste0(gene ,"ovarian cancer Disease Free survival" )
    ) 
    ggsave(filename = paste0(gene,"ovarian cancer Disease Free survival.png"),plot = print(SurvPlot_df),
           height = 6,width = 6,units = 'in')
    
  }
  
  if (pvalue_o < 0.05) {
    
    SurvPlot_o <- ggsurvplot(
      km_trt_fit8,
      data = full_data_m_gene, size = 1,                # change line size
      palette = c("#E7B800","#2E9FDF"),# custom color palettes
      conf.int = TRUE,         # Add confidence interval
      pval = T,             # Add p-value
      risk.table = TRUE,       # Add risk table
      risk.table.col = "strata",# Risk table color by groups
      legend.labs =
        c("No alteration","Alteration"),   # Change legend labels
      risk.table.height = 0.25,# Useful to change when you have multiple groups
      ggtheme = theme_bw(),     # Change ggplot2 theme
      xlab = "Overall Survival (Months)", ylab = "Overall Survival", title = paste0(gene," - ovarian cancer (Overall survival)")
    ) 
    ggsave(filename = paste0(gene," - ovarian cancer - overall survival.png"),plot = print(SurvPlot_o),
           height = 6,width = 6,units = 'in')
    
  }
  
}

write_csv(result_table,"2 - Mutational analysis/Survival analysis/result_survival_testis_antigen.csv")
