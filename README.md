# Species Distribution Modeling with H2O AutoML

This repository contains an R pipeline for training, saving, and projecting species distribution models (SDMs) using H2O’s AutoML.

---

## What it does

1. **Prepares species data**  
   - Reads species CSVs (e.g., different individuals like elephants).
   - Generates pseudo-absences (equal number, randomly sampled).

2. **Extracts environmental variables**  
   - Extracts raster values for all points (presence + background).

3. **Trains models with H2O AutoML**  
   - Runs `h2o.automl()` to test multiple algorithms: GBMs, GLMs, XGBoosts, DRFs, Stacked Ensembles.
   - Selects the best-performing model (leader) based on AUC.

4. **Saves only the leader model**  
   - Each species gets its own leader model — could be different types.
   - Leader model saved for future projections.

5. **Generates current predictions**  
   - Maps habitat suitability for current climate.

6. **Projects to future climate scenarios**  
   - Resamples future rasters to match training extents.
   - Runs each leader model across all time periods & SSPs.

---

## Why this design is ecologically sound

- **Individual variation**: Each animal may have unique movement or habitat use.
- **Pseudo-absences**: Balances data, standard for presence-only data.
- **AutoML adaptiveness**: Avoids algorithm bias by testing many models.
- **Leader-only storage**: Keeps storage clean & deployment practical.
- **Future raster resampling**: Ensures predictions align spatially with training data.

---

## How to use

1. Drop your CSVs in the species folders.
2. Organize future rasters under `Envi/future/<time_period>/<ssp>/`.
3. Run the script — it loops through all species, scenarios, and periods.
4. Outputs:
   - Saved leader models.
   - Current suitability maps.
   - Future scenario projections.

---

## Key takeaway

This pipeline ensures your SDMs:
- Respect biological realism.
- Generalize well.
- Are fully reproducible and adaptable for conservation and climate impact assessments.

---
