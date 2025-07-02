####Discription---- What our final training loop does
For each species CSV:
We run one call to h2o.automl().
AutoML builds many models internally: GBMs, GLMs, XGBoosts, DRFs, etc.
Plus ensembles like StackedEnsemble_BestOfFamily and StackedEnsemble_AllModels or even even if stand alone models performs the best in the leaderboard.

Then: leader <- aml@leader
takes only the single “best model” (the “leader” of that AutoML run).
And we save just this one: model_path <- h2o.saveModel(leader, path=species_dir, force=TRUE)
So we only save the leader model, not the whole leaderboard.
Each species’ best model may be a different type!AutoML picks the leader based on validation AUC — sometimes that’s a StackedEnsemble_BestOfFamily, sometimes AllModels, sometimes a base GBM.
 exactly how h2o.automl() works!
The other models exist in the AutoML leaderboard but they’re usually used internally for stacking.
Key point
For projections, we only need to load the leader for each species
Why does each species (elephant) get a different leader model?
-We have three separate individuals with different telemetry CSVs.
-Even though they’re all inside the same national park (~12,000 ha):
-Each elephant’s movement, habitat preference, or micro-range may vary.
-Their data points may differ in sample size, spatial spread, or seasonal coverage.

------How does this affect H2O AutoML?
AutoML always selects the best model for that exact dataset.
If one elephant has a broader spread, a StackedEnsemble may outperform.
If another has fewer or more clustered points, a simple GBM may be more stable and become the leader.

So the leader model type can absolutely vary because:
 Different sample sizes → affects variance and bias tradeoff.
 Different spatial autocorrelation → affects variable importance and predictability.
 Different background/pseudo-absence generation → affects how well ensembles generalize.

Imagine:
Elephant A ranges widely → more environmental gradient → ensemble helps blend models.

Elephant B stays near the river → simpler structure → GBM wins.

Elephant C has patchy points → AllModels ensemble works best.

🗝️ Key takeaway
it means AutoML pipeline is adapting to real biological variation in our tracking data.
This is exactly why using automl() is robust: it tests multiple algorithms and picks what generalizes best and this variation is ecologically realistic.

####Current day models script----
library(terra)
library(h2o)
library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)

dirs <- list(
  A = list(env_dir = "Envi/A", sp_dir = "Data/A", res_dir = "H2o_Results/A"),
  B = list(env_dir = "Envi/B", sp_dir = "Data/B", res_dir = "H2o_Results/B")
)

h2o.init()

for(run in names(dirs)){
  env_dir <- dirs[[run]]$env_dir
  sp_dir  <- dirs[[run]]$sp_dir
  res_dir <- dirs[[run]]$res_dir
  dir.create(res_dir, showWarnings=FALSE, recursive=TRUE)
  
  env_files <- list.files(env_dir, pattern="\\.tif$", full.names=TRUE)
  env_names <- tools::file_path_sans_ext(basename(env_files))
  env_vect  <- rast(env_files)
  names(env_vect) <- make.names(env_names, unique=TRUE)
  all_vars  <- names(env_vect)
  
  metrics_list <- list()
  varimp_list  <- list()
  model_list   <- list()
  
  for(f in list.files(sp_dir, pattern="\\.csv$", full.names=TRUE)){
    sp <- read_csv(f, show_col_types=FALSE)
    pres <- sp %>% select(lon, lat) %>% mutate(pa=1)
    bg_pts <- spatSample(env_vect[[1]], size=nrow(pres), method="random", as.points=TRUE)
    bg_coords <- crds(bg_pts)
    bg_df <- data.frame(lon=bg_coords[,1], lat=bg_coords[,2], pa=0)
    df_sp <- bind_rows(pres, bg_df)
    pts <- vect(df_sp, geom=c("lon","lat"), crs=crs(env_vect))
    vals <- terra::extract(env_vect, pts)[,-1]
    df <- bind_cols(df_sp, as.data.frame(vals)) %>% na.omit()
    hf <- as.h2o(df)
    hf["pa"] <- as.factor(hf["pa"])
    aml <- h2o.automl(x=all_vars, y="pa", training_frame=hf, max_models=30, seed=123)
    leader <- aml@leader
    dset <- tools::file_path_sans_ext(basename(f))
    
    species_dir <- file.path(res_dir, dset)
    dir.create(species_dir, recursive=TRUE, showWarnings=FALSE)
    model_path <- h2o.saveModel(leader, path=species_dir, force=TRUE)
    
    model_list[[dset]] <- leader
    
    perf <- h2o.performance(leader, newdata=hf)
    metrics_list[[dset]] <- data.frame(dataset=dset,
                                       rmse=h2o.rmse(perf),
                                       auc=h2o.auc(perf),
                                       logloss=h2o.logloss(perf))
    varimp <- tryCatch(as.data.frame(h2o.varimp(leader)), error=function(e) NULL)
    if(!is.null(varimp) && nrow(varimp)>0){
      varimp$dataset <- dset
      varimp_list[[dset]] <- varimp
      pdf(file.path(species_dir, paste0("varimp_", dset, ".pdf")))
      h2o.varimp_plot(leader, num_of_features=10)
      dev.off()
      pdf(file.path(species_dir, paste0("partial_", dset, ".pdf")))
      h2o.partialPlot(object=leader, data=hf, cols=head(all_vars,5))
      dev.off()
    }
    
    # Save the current prediction to species folder
    pred_rast <- terra::app(env_vect, fun=function(vals){
      df_pred <- as.data.frame(vals)
      hf2 <- as.h2o(df_pred)
      p <- h2o.predict(leader, hf2)
      as.vector(p[["p1"]])
    })
    writeRaster(pred_rast, filename=file.path(species_dir, paste0("prediction_", dset, ".tif")), overwrite=TRUE)
  }
  
  metrics_all <- bind_rows(metrics_list)
  varimp_all  <- bind_rows(varimp_list)
  
  ggsave(file.path(res_dir, "metrics_rmse.png"),
         ggplot(metrics_all, aes(x=dataset, y=rmse, fill=dataset)) +
           geom_col(show.legend=FALSE) +
           labs(x="Dataset", y="RMSE") +
           theme(axis.text.x=element_text(angle=45,hjust=1)),
         width=8, height=4)
  
  ggsave(file.path(res_dir, "metrics_auc.png"),
         ggplot(metrics_all, aes(x=dataset, y=auc, fill=dataset)) +
           geom_col(show.legend=FALSE) +
           labs(x="Dataset", y="AUC") +
           theme(axis.text.x=element_text(angle=45,hjust=1)),
         width=8, height=4)
  
  if(nrow(varimp_all)>0){
    ggsave(file.path(res_dir, "varimp_heatmap.png"),
           ggplot(varimp_all %>% select(dataset, variable, relative_importance),
                  aes(x=variable, y=dataset, fill=relative_importance)) +
             geom_tile() +
             labs(x="Environmental Variable", y="Dataset", fill="Rel. Importance") +
             theme(axis.text.x=element_text(angle=90,hjust=1)),
           width=10, height=6)
  }
}

h2o.shutdown(prompt=FALSE)

## Future Predictions----
#resample future rasters to the extents of A
library(terra)

ref_files <- list.files("D:/Harin/Projects/SDM/SDM/All/Envi/A", pattern="\\.tif$", full.names=TRUE)
ref_stack <- rast(ref_files)
future_dir <- "D:/Harin/Projects/SDM/SDM/All/Envi/future"
time_periods <- c("2021-2040", "2041-2060", "2061-2080", "2081-2100")
ssps <- c("ssp126", "ssp245", "ssp370", "ssp585")

for(tp in time_periods){
  for(ssp in ssps){
    in_path <- file.path(future_dir, tp, ssp)
    out_path <- file.path(future_dir, tp, ssp, "resampled")
    dir.create(out_path, recursive=TRUE, showWarnings=FALSE)
    
    r_files <- list.files(in_path, pattern="\\.tif$", full.names=TRUE)
    for(f in r_files){
      r <- rast(f)
      r_resampled <- resample(r, ref_stack[[1]], method="bilinear")
      writeRaster(r_resampled, filename=file.path(out_path, basename(f)), overwrite=TRUE)
    }
  }
}

## Future Predictions — Multi-Species, One Scenario----
library(terra)
library(h2o)
library(tools)

species_list <- c("Species1", "Species2", "Species3")  # Replace with your real names
future_env <- "Envi/future/2021-2040/ssp126/resampled"
output_base <- "H2o_Results/future_test/2021-2040/ssp126"
env_files <- list.files(future_env, pattern="\\.tif$", full.names=TRUE)
stopifnot(length(env_files) > 0)
env_names <- make.names(file_path_sans_ext(basename(env_files)), unique=TRUE)
env_stack <- rast(env_files)
names(env_stack) <- env_names

vals_df <- as.data.frame(env_stack, na.rm=FALSE)
cat("Loaded raster rows:", nrow(vals_df), " Variables:", ncol(vals_df), "\n")

h2o.init()
hf_future <- as.h2o(vals_df)

for(sp in species_list) {
  #Exact model path for this species
  model_dir <- file.path("H2o_Results/A", sp)
  model_file <- list.files(model_dir, pattern="StackedEnsemble", full.names=TRUE)[1]
  
  cat("\nProjecting:", sp, "\n")
  cat("Model file:", basename(model_file), "\n")
  
  model <- h2o.loadModel(model_file)
  pred <- h2o.predict(model, hf_future)
  pred_vals <- as.vector(as.data.frame(pred[["p1"]])[,1])
  
  cat("Prediction length:", length(pred_vals), " | Raster cells:", ncell(env_stack[[1]]), "\n")
  
  if(length(pred_vals) == ncell(env_stack[[1]])){
    r_out <- setValues(env_stack[[1]], pred_vals)
    plot(r_out, main=paste("Projected Suitability -", sp))
    out_dir <- file.path(output_base, sp)
    dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
    writeRaster(r_out, filename=file.path(out_dir, paste0("prediction_", sp, ".tif")), overwrite=TRUE)
    cat("Saved:", sp, "\n")
  } else {
    cat("Skipped:", sp, "- length mismatch\n")
  }
}

h2o.shutdown(prompt=FALSE)

library(terra)
library(h2o)
library(tools)

setwd("D:/Harin/Projects/SDM/SDM/All")

model_path <- "H2o_Results/A/<species_name>/StackedEnsemble_BestOfFamily_1_AutoML_1_20250627_135706"
future_env <- "Envi/future/2021-2040/ssp126/resampled"
output_dir <- "Projections/2021-2040/ssp126/<species_name>"
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)

env_files <- list.files(future_env, pattern="\\.tif$", full.names=TRUE)
stopifnot(length(env_files) > 0)

env_names <- make.names(file_path_sans_ext(basename(env_files)), unique=TRUE)
env_stack <- rast(env_files)
names(env_stack) <- env_names

h2o.init()

vals_df <- as.data.frame(env_stack, na.rm=FALSE)
hf_future <- as.h2o(vals_df)

cat("Loaded raster rows:", nrow(vals_df), " Variables:", ncol(vals_df), "\n")

model <- h2o.loadModel(model_path)

cat("Projecting with model:", basename(model_path), "\n")
pred <- h2o.predict(model, hf_future)
pred_vals <- as.vector(as.data.frame(pred[["p1"]])[,1])

cat("Prediction length:", length(pred_vals), "\n")
cat("Raster cells:", ncell(env_stack[[1]]), "\n")

if(length(pred_vals) == ncell(env_stack[[1]])){
  r_out <- setValues(env_stack[[1]], pred_vals)
  plot(r_out, main="Projected Suitability")
  writeRaster(r_out, filename=file.path(output_dir, paste0("prediction_", basename(model_path), ".tif")), overwrite=TRUE)
  cat("Prediction raster saved.\n")
} else {
  cat("Prediction skipped: length mismatch\n")
}

h2o.shutdown(prompt=FALSE)


#looooooooooooppppp-----
library(terra)
library(h2o)
library(tools)
species_list <- c("thinned_E3A", "thinned_E4A")

#Time periods and SSPs
time_periods <- c("2021-2040", "2041-2060", "2061-2080", "2081-2100")
ssps <- c("ssp126", "ssp245", "ssp370", "ssp370", "ssp585")

future_dir <- "Envi/future"
output_base <- "H2o_Results/future_test"

h2o.init()

for(tp in time_periods) {
  for(ssp in ssps) {
    env_path <- file.path(future_dir, tp, ssp, "resampled")
    env_files <- list.files(env_path, pattern="\\.tif$", full.names=TRUE)
    
    if(length(env_files) == 0) {
      cat("No rasters for:", tp, ssp, "\n")
      next
    }
    
    env_names <- make.names(file_path_sans_ext(basename(env_files)), unique=TRUE)
    env_stack <- rast(env_files)
    names(env_stack) <- env_names
    
    vals_df <- as.data.frame(env_stack, na.rm=FALSE)
    cat("\n\nLoaded:", tp, ssp, "| Rows:", nrow(vals_df), " Vars:", ncol(vals_df), "\n")
    
    hf_future <- as.h2o(vals_df)
    
    for(sp in species_list) {
      model_dir <- file.path("H2o_Results/A", sp)
      model_file <- list.files(model_dir, pattern="StackedEnsemble", full.names=TRUE)[1]
      
      if (is.na(model_file)) {
        cat("No model found for:", sp, "\n")
        next
      }
      
      cat("Projecting:", sp, " | Model:", basename(model_file), "\n")
      
      model <- h2o.loadModel(model_file)
      pred <- h2o.predict(model, hf_future)
      pred_vals <- as.vector(as.data.frame(pred[["p1"]])[,1])
      
      cat("Prediction:", length(pred_vals), " | Raster cells:", ncell(env_stack[[1]]), "\n")
      
      if(length(pred_vals) == ncell(env_stack[[1]])) {
        r_out <- setValues(env_stack[[1]], pred_vals)
        out_dir <- file.path(output_base, tp, ssp, sp)
        dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
        writeRaster(r_out, filename=file.path(out_dir, paste0("prediction_", sp, ".tif")), overwrite=TRUE)
        cat("Saved:", sp, tp, ssp, "\n")
      } else {
        cat("Skipped:", sp, "- length mismatch\n")
      }
    }
  }
}
h2o.shutdown(prompt=FALSE)
