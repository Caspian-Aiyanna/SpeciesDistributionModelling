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
