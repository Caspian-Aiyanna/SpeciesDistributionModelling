library(terra)
library(h2o)
library(dplyr)
library(readr)
library(ggplot2)
library(stringr)

# Species folders
species_dirs <- list(
  Species1 = "D:/Harin/Projects/SDM/SDM/All/H2o_Results/A/thinned_E3A",
  Species2 = "D:/Harin/Projects/SDM/SDM/All/H2o_Results/A/thinned_E4A",
  Species3 = "D:/Harin/Projects/SDM/SDM/All/H2o_Results/A/thinned_E5A"
)

h2o.init()

for (sp in names(species_dirs)) {
  species_dir <- species_dirs[[sp]]
  
  #Find the leader ensemble model only
  model_path <- list.files(
    species_dir, 
    pattern = "StackedEnsemble_BestOfFamily|StackedEnsemble_AllModels", 
    full.names = TRUE
  )[1]
  
  if (is.na(model_path)) {
    cat("No StackedEnsemble found for:", sp, " — skipping.\n")
    next
  }
  
  cat("Visualizing ensemble for:", sp, "\n")
  ensemble <- h2o.loadModel(model_path)
  
  if (is.null(ensemble@model$metalearner)) {
    cat("ℹModel is not an ensemble for:", sp, " — skipping plot.\n")
    next
  }
  
  meta_name <- ensemble@model$metalearner$name
  metalearner <- h2o.getModel(meta_name)
  varimp <- h2o.varimp(metalearner)
  
  if (is.null(varimp) || nrow(varimp) == 0) {
    cat("No base learner contributions found for:", sp, "\n")
    next
  }
  
  #Shorten names
  varimp$ShortName <- str_extract(varimp$variable,
                                  "^[^_]+(_grid_[0-9]+)?|GBM_[0-9]+|XRT_[0-9]+|GLM_[0-9]+|DRF_[0-9]+|DeepLearning_[0-9]+")
  varimp$ShortName <- ifelse(is.na(varimp$ShortName), varimp$variable, varimp$ShortName)
  
  #Normalize and group by ShortName to sum duplicate rows
  varimp$percentage <- (varimp$percentage / sum(varimp$percentage)) * 100
  varimp_summary <- varimp %>%
    group_by(ShortName) %>%
    summarise(percentage = sum(percentage)) %>%
    arrange(desc(percentage))
  
  #Plot
  p <- ggplot(varimp_summary, aes(x = reorder(ShortName, percentage), y = percentage, fill = percentage)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = sprintf("%.1f%%", percentage)),
              hjust = -0.1, size = 4, color = "black") +
    scale_fill_gradient(low = "#9ecae1", high = "#08519c") +
    coord_flip() +
    labs(title = paste("Base Learner Contributions -", sp),
         x = "Base Learner",
         y = "Contribution (%)",
         fill = "Contribution %") +
    theme_minimal(base_size = 16) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      plot.title = element_text(face = "bold")
    )
  #Export HD JPG
  out_file <- file.path(species_dir, paste0("BaseLearner_Contributions_", sp, ".jpg"))
  ggsave(out_file, plot = p, width = 10, height = 6, dpi = 500)
  
  cat("Saved base learner plot for:", sp, "\n")
}

h2o.shutdown(prompt = FALSE)



#Variable importance per each base learner model----
library(h2o)
library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(patchwork)

h2o.init()

ensemble_paths <- list(
  E3A = "D:/Harin/Projects/SDM/SDM/All/H2o_Results/A/thinned_E3A/StackedEnsemble_BestOfFamily_1_AutoML_9_20250701_140429",
  E4A = "D:/Harin/Projects/SDM/SDM/All/H2o_Results/A/thinned_E4A/StackedEnsemble_AllModels_1_AutoML_10_20250701_141911"
)

data_paths <- list(
  E3A = "D:/Harin/Projects/SDM/SDM/All/Data/A/thinned_E3A.csv",
  E4A = "D:/Harin/Projects/SDM/SDM/All/Data/A/thinned_E4A.csv"
)

for (sp in names(ensemble_paths)) {
  cat("🔍 Processing:", sp, "\n")
  
  ensemble <- h2o.loadModel(ensemble_paths[[sp]])
  data <- h2o.importFile(data_paths[[sp]])
  
  combined_varimp <- data.frame(variable = character(), relative_importance = numeric())
  
  base_models <- ensemble@model$base_models
  
  for (base_name in base_models) {
    base_model <- h2o.getModel(base_name)
    if (grepl("Metalearner", base_name)) next
    
    cat("   📊 Base learner:", base_name, "\n")
    
    varimp <- h2o.varimp(base_model)
    if (!is.null(varimp)) {
      varimp_df <- as.data.frame(varimp)[, c("variable", "relative_importance")]
      combined_varimp <- bind_rows(combined_varimp, varimp_df)
    }
  }
  
  if (nrow(combined_varimp) > 0) {
    # Sum and normalize across base learners
    combined_summary <- combined_varimp %>%
      group_by(variable) %>%
      summarise(relative_importance = sum(relative_importance)) %>%
      mutate(relative_importance = relative_importance / sum(relative_importance) * 100) %>%
      arrange(desc(relative_importance))
    
    p <- ggplot(combined_summary, aes(x = reorder(variable, relative_importance), y = relative_importance)) +
      geom_col(fill = "#2c7fb8") +
      coord_flip() +
      labs(title = paste("Combined Variable Importance -", sp),
           x = "Variable", y = "Relative Importance (%)") +
      theme_minimal()
    
    ggsave(file.path(dirname(ensemble_paths[[sp]]), paste0("Combined_VarImp_", sp, ".jpg")),
           plot = p, width = 8, height = 6, dpi = 300)
  }
}

h2o.shutdown(prompt = FALSE)

