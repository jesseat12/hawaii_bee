# Jesse Tabor
# March 25, 2021
# Hawaii non-native bee study using BIOMOD2

# This script creates species distribution models (SDMs).
# We create the SDMs using species occurrence data (longitude, latitude) from gbif.org
# and current bioclimatic variables downloaded from Worldclim.org. We project these models 
# onto current climate scenarios in Hawaii to predict suitable habitat for non-native bees. 
# We also project models onto future climate scenarios in Hawaii to predict future suitable habitat.

# Set where the program saves and looks for data.
setwd()
list.files()

# Load libraries into the R environment.
library(rgbif)
library(biomod2)
library(ggplot2)
library(gridExtra)
library(ade4)
library(rgdal)
library(raster)
library(maptools)
library(rasterVis)
library(latticeExtra)
library(lattice)
library(sp)
library(dplyr)
library(devtools)
print("libraries loaded")

# Create a folder in your working directory and write your plots into that folder.
dir.create("plots")

# Load species occurence data into R
bee <- read.table("species.csv",
                  header=TRUE, 
                  sep=",",
                  row.names=NULL)
print("species occurences loaded")

# Load your bioclimatic variables into R.
bio_1 <-raster("bioclim/bio_1.grd")
bio_2 <-raster("bioclim/bio_2.grd")
bio_3 <-raster("bioclim/bio_3.grd")
bio_4 <-raster("bioclim/bio_4.grd")
bio_5 <-raster("bioclim/bio_5.grd")
bio_6 <-raster("bioclim/bio_6.grd")
bio_7 <-raster("bioclim/bio_7.grd")
bio_8 <-raster("bioclim/bio_8.grd")
bio_9 <-raster("bioclim/bio_9.grd")
bio_10 <-raster("bioclim/bio_10.grd")
bio_11 <-raster("bioclim/bio_11.grd")
bio_12 <-raster("bioclim/bio_12.grd")
bio_13 <-raster("bioclim/bio_13.grd")
bio_14 <-raster("bioclim/bio_14.grd")
bio_15 <-raster("bioclim/bio_15.grd")
bio_16 <-raster("bioclim/bio_16.grd")
bio_17 <-raster("bioclim/bio_17.grd")
bio_18 <-raster("bioclim/bio_18.grd")
bio_19 <-raster("bioclim/bio_19.grd")
print("area1 range raster loaded")

# Stack bioclim rasters in raster stack.
bioclim <- stack(c(bio_1, bio_10, bio_11,
                      bio_12, bio_13, bio_14,
                      bio_15, bio_16, bio_17,
                      bio_18, bio_19, bio_2,
                      bio_3, bio_4, bio_5,
                      bio_6, bio_7, bio_8,
                      bio_9))
print("area1 range raster stacked")

# ------------------------------------------------------------------------------Principal component analysis (PCA)-------------------------------------------------------------------------------------------------

# To reduce multicollinearity within the environmental variables, a principal component analysis is conducted to highlight the relationship between the target species 
# occurrences and the specific environmental combinations. Use the PCA plot to select a set of variables that are not too colinear (two variables pointing in
# orthogonal directions are independent, two variables pointing in the same or opposite directions are highly dependent, positively or negatively, respectively),
# and significantly contribute to the overall environmental variation (the longer the arrow, the more important the variable).

# Obtain identifiers of the cells where species occurs (lat/long).
points_bee <- data.frame(bee)[1:999, c("Longitude","Latitude")]
print("identifiers ready")

# Extract the value from each cell and put it in one central location.
bee_cell_id <- cellFromXY(subset(bioclim,1),points_bee)
print("extracted cell id")

# Convert raster object into a data frame and remove non-defined area from the dataset, which gives values to every point.
bioclim_df <- na.omit(as.data.frame(bioclim))
print("removed non-defined area complete")

# Perform a principal component analysis (PCA) on environment to avoid model overfitting by decreasing the amount of variables used.
pca_NA <- dudi.pca(bioclim_df, scannf = F, nf = 2)
print("pca complete")

# Write PCA plot to plots folder and determine which variables to use in SDMs.
jpeg("/plots/species_pca.jpeg", type="cairo")
s.corcircle(pca_NA$co,clabel = 1)
dev.off()

# ------------------------------------------------------------------------------BIOMOD_Modeling()------------------------------------------------------------------------------------------------------------------

# Create our initial models that ensemble models will be built from.

# Stack your subset area1 range rasters into a new rasterstack.
bioclim <- stack(c(bio_2,bio_4,bio_11,bio_12))
print("subset stack area1 rasters complete")

# Rearrange input data (species occurrence data and range raster) to make sure they can be used within biomod2. 
# The function allows to select pseudo-absences or background data in the case that true absences data are not available,
# or to add pseudo-absence data to an existing set of absence.
bee_data <- BIOMOD_FormatingData(resp.var = rep(1,nrow(bee)),
                                 expl.var = bioclim,
                                 resp.xy = bee[,c('Longitude','Latitude')],
                                 resp.name = "name.model",
                                 PA.nb.rep = 3,PA.nb.absences = 10000,
                                 PA.strategy = 'random')
print("BIOMOD_FormattingData complete")

# Check to see if BIOMOD_FormatingData() function worked.
bee_data

# Plot bee_data and save to the plots folder.
jpeg("/plots/bee_data.jpeg", type="cairo")
plot(bee_data)
dev.off()

# Parameterize and/or tune biomod single models options with BIOMOD_ModelingOptions() function.
bee_opt <- BIOMOD_ModelingOptions(GLM = list(type = 'quadratic',
                                             interaction.level = 1),
                                  GAM = list(algo ='GAM_mgcv'),
                                  GBM = list(n.trees = 1000),
                                  MARS = list(type = 'quadratic'))
print("BIOMOD_ModelingOptions complete")

# Next we run 9 different models on bee_data paired individually with 3 different pseudo-absence datasets
# (created with BIOMOD_FormatingData() function). We do this 4 times. This creates 108 total models.
# This function allows to calibrate and evaluate a range of species distribution models techniques
# run over a given species. Calibrations are made on the whole sample or a random subpart. The
# predictive power of the different models is estimated using a range of evaluation metrics.
bee_models <- BIOMOD_Modeling(data = bee_data, models = c("GLM","GAM","GBM","MARS","RF","ANN","FDA","CTA","MAXENT.Phillips"),
                              models.options = bee_opt, 
                              NbRunEval = 4, DataSplit = 80,
                              VarImport = 3, do.full.models = F,
                              modeling.id = "ex.2")
print("BIOMOD_Modeling complete")

# Get model evaluation scores.
bee_models_scores <- get_evaluations(bee_models)

# dim() function: Get or set the number of rows, columns, and layers of a Raster* object.
dim(bee_models_scores)

# Retrieve or set the dimnames of an object.
dimnames(bee_models_scores)

# model_scores_graph() function: This function is a graphic tool to represent evaluation scores of models
# produced with biomod2 according to 2 different evaluation methods. Models can be grouped in several ways
# (by algo, by CV run, ...) to highlight potential differences in models quality due to chosen models,
# or cross validation sampling bias. Points represent mean of evaluation score for a given condition and
# the lines represent the associated standard deviations.

# Model scores. Save the plot in the plots folder.
jpeg("/plots/area1_model_scores.jpeg", type="cairo")
models_scores_graph(bee_models, by = "models",metrics = c("ROC","TSS"),
                    xlim = c(0.1,1), ylim = c(0.1,1))
dev.off()

# Get variable importance in the models.
(bee_models_var_import <- get_variables_importance(bee_models))

# Cross-validation RUN scores. Save the plot in the plots folder.
jpeg("/plots/area1_run_scores.jpeg", type="cairo")
models_scores_graph(bee_models, by = "cv_run",metrics = c("ROC","TSS"),
                    xlim = c(0.1,1), ylim = c(0.1,1))
dev.off()

# Pseudo-absence scores. Save the plot in the plots folder.
jpeg("/plots/area1_pa_scores.jpeg", type="cairo")
models_scores_graph(bee_models, by = "data_set",metrics = c("ROC","TSS"),
                    xlim = c(0.1,1), ylim = c(0.1,1))
dev.off()

# Calculate the MEAN of variable importance by algorighm. This function tells us the importance of the individual
# variables in the model.
apply(bee_models_var_import, c(1,2), mean)

# Analyze how each environmental variable influence the species probability of presence by creating response curve plots.
# These are graphical visualizations of the response curve of each variable. How does each environmental variable influence
# probability of presence? Each line corresponds to a diffrent model. To do this we first have to load the produced models.
bee_glm <- BIOMOD_LoadModels(bee_models, models = 'GLM')
bee_gbm <- BIOMOD_LoadModels(bee_models, models = 'GBM')
bee_rf <- BIOMOD_LoadModels(bee_models, models = 'RF')
bee_cta <- BIOMOD_LoadModels(bee_models, models = 'CTA')
bee_ann <- BIOMOD_LoadModels(bee_models, models = 'ANN')
bee_fda <- BIOMOD_LoadModels(bee_models, models = 'FDA')
bee_mars <- BIOMOD_LoadModels(bee_models, models = 'MARS')
bee_gam <- BIOMOD_LoadModels(bee_models, models = 'GAM')
bee_phi <- BIOMOD_LoadModels(bee_models, models = 'MAXENT.Phillips')

# biomod2::response.plot2() function creates the plot. Save the plot in the plots folder.
jpeg("/plots/glm.jpeg", type="cairo")
glm_eval_strip <- biomod2::response.plot2(models = bee_glm, Data = get_formal_data(bee_models,'expl.var'),
                                          show.variables = get_formal_data(bee_models,'expl.var.names'),
                                          do.bivariate = FALSE, fixed.var.metric = 'median', legend = FALSE,
                                          display_title = FALSE, data_species = get_formal_data(bee_models, 'resp.var'))
dev.off()

# biomod2::response.plot2() function creates the plot. Save the plot in the plots folder.
jpeg("/plots/gbm.jpeg", type="cairo")
gbm_eval_strip <- biomod2::response.plot2(models = bee_gbm, Data = get_formal_data(bee_models, 'expl.var'),
                                          show.variables = get_formal_data(bee_models, 'expl.var.names'),
                                          do.bivariate = FALSE, fixed.var.metric = 'median', legend = FALSE,
                                          display_title = FALSE, data_species = get_formal_data(bee_models, 'resp.var'))
dev.off()

# biomod2::response.plot2() function creates the plot. Save the plot in the plots folder.
jpeg("/plots/rf.jpeg", type="cairo")
rf_eval_strip <- biomod2::response.plot2(models = bee_rf, Data = get_formal_data(bee_models, 'expl.var'),
                                         show.variables = get_formal_data(bee_models, 'expl.var.names'),
                                         do.bivariate = FALSE, fixed.var.metric = 'median', legend = FALSE,
                                         display_title = FALSE, data_species = get_formal_data(bee_models, 'resp.var'))
dev.off()

# biomod2::response.plot2() function creates the plot. Save the plot in the plots folder.
jpeg("/plots/gam.jpeg", type="cairo")
gam_eval_strip <- biomod2::response.plot2(models = bee_gam, Data = get_formal_data(bee_models, 'expl.var'),
                                          show.variables = get_formal_data(bee_models, 'expl.var.names'),
                                          do.bivariate = FALSE, fixed.var.metric = 'median', legend = FALSE,
                                          display_title = FALSE, data_species = get_formal_data(bee_models, 'resp.var'))
dev.off()

# biomod2::response.plot2() function creates the plot. Save the plot in the plots folder.
jpeg("/plots/cta.jpeg", type="cairo")
cta_eval_strip <- biomod2::response.plot2(models = bee_cta, Data = get_formal_data(bee_models, 'expl.var'),
                                          show.variables = get_formal_data(bee_models, 'expl.var.names'),
                                          do.bivariate = FALSE, fixed.var.metric = 'median', legend = FALSE,
                                          display_title = FALSE, data_species = get_formal_data(bee_models, 'resp.var'))
dev.off()

# biomod2::response.plot2() function creates the plot. Save the plot in the plots folder.
jpeg("/plots/ann.jpeg", type="cairo")
ann_eval_strip <- biomod2::response.plot2(models = bee_ann, Data = get_formal_data(bee_models, 'expl.var'),
                                          show.variables = get_formal_data(bee_models, 'expl.var.names'),
                                          do.bivariate = FALSE, fixed.var.metric = 'median', legend = FALSE,
                                          display_title = FALSE, data_species = get_formal_data(bee_models, 'resp.var'))
dev.off()

# biomod2::response.plot2() function creates the plot. Save the plot in the plots folder.
jpeg("/plots/fda.jpeg", type="cairo")
fda_eval_strip <- biomod2::response.plot2(models = bee_fda, Data = get_formal_data(bee_models, 'expl.var'),
                                          show.variables = get_formal_data(bee_models, 'expl.var.names'),
                                          do.bivariate = FALSE, fixed.var.metric = 'median', legend = FALSE,
                                          display_title = FALSE, data_species = get_formal_data(bee_models, 'resp.var'))
dev.off()

# biomod2::response.plot2() function creates the plot. Save the plot in the plots folder.
jpeg("/plots/mars.jpeg", type="cairo")
mars_eval_strip <- biomod2::response.plot2(models = bee_mars, Data = get_formal_data(bee_models, 'expl.var'),
                                           show.variables = get_formal_data(bee_models, 'expl.var.names'),
                                           do.bivariate = FALSE, fixed.var.metric = 'median', legend = FALSE,
                                           display_title = FALSE, data_species = get_formal_data(bee_models, 'resp.var'))
dev.off()

# biomod2::response.plot2() function creates the plot. Save the plot in the plots folder.
jpeg("/plots/maxent.jpeg", type="cairo")
phi_eval_strip <- biomod2::response.plot2(models = bee_phi, Data = get_formal_data(bee_models, 'expl.var'),
                                          show.variables = get_formal_data(bee_models, 'expl.var.names'),
                                          do.bivariate = FALSE, fixed.var.metric = 'median', legend = FALSE,
                                          display_title = FALSE, data_species = get_formal_data(bee_models, 'resp.var'))
dev.off()
print("response plots complete")

# Congratulations you have built 108 SDMs calibrated from the area1 range and
# analyzed the importance of algorithms, pseudo-absence, runs, and variables  in model performance. 
# Now we we will create an ensemble by averaging  all 108 SDMs into 1 final ensemble model. 

# ------------------------------------------------------------------------------Ensemble modeling------------------------------------------------------------------------------------------------------------------

#BIOMOD_EnsembleModeling combines models and make ensemble predictions built with BIOMOD_Modeling.
#The ensemble predictions can also be evaluated against the original data given to BIOMOD_Modeling.
#Biomod2 proposes a range of options to build ensemble models and predictions and to assess the
#modeling uncertainty. The created ensemble models can then be used to project distributions over
#space and time as classical biomod2 models.

# In this case to reduce the number of outputs we only concider two "ensembling" options:
# Committee averaging and weighted mean. We also produce coefficient of variation that tell us the extent the models agree or disagree. 
# In this case we made a decision to mix  all models (i.e. all algorithms, all pseudo-absence sampling, all cross-validation runs) 
# to produce our ensemble models. TSS is used as the evaluation reference for committee building and defining weights. 
# In this case only models with a TSS greater than or equal to 0.8 are kept to build the final ensemble model. 
bee_ensemble_models <- BIOMOD_EnsembleModeling(modeling.output = bee_models,
                                               em.by = 'all',
                                               eval.metric = 'TSS',
                                               eval.metric.quality.threshold = 0.8,
                                               models.eval.meth = c('KAPPA','TSS','ROC'),
                                               prob.mean = FALSE,
                                               prob.cv = TRUE, 
                                               committee.averaging = TRUE,
                                               prob.mean.weight = TRUE,
                                               VarImport = 0)

# Now check the scores for the ensemble models.
(bee_ensemble_models_scores <- get_evaluations(bee_ensemble_models))
print("area1 ensemble modeling complete")

# ------------------------------------------------------------------------------Current area1 range projection----------------------------------------------------------------------------------------------------

# For all the models currently implemented, BIOMOD_Projection() function is able to project potential distributions of
# species in other areas, other resolutions or other time scales. We now project all 108 current individual area1 range SDMs.
bee_models_proj_current <- BIOMOD_Projection(modeling.output = bee_models,
                                             new.env = bioclim,
                                             proj.name = "current", 
                                             binary.meth = 'TSS', 
                                             do.stack = FALSE)
print("area1 BIOMOD_Projection complete")

# This function use projections of individual models and ensemble models from BIOMOD_EnsembleModeling
# to build an ensemble of species projections over space and time. We now project current area1 range ensemble SDMs.
bee_ensemble_models_proj_current <- BIOMOD_EnsembleForecasting(EM.output = bee_ensemble_models,
                                                               projection.output = bee_models_proj_current,
                                                               binary.meth = "TSS",
                                                               output.format = ".img",
                                                               do.stack = FALSE)
print("area1 BIOMOD_EnsembleForcasting complete")

# Stack the area1 range ensemble SDMs.
stk_bee_ensemble_models_proj_current <- get_predictions(bee_ensemble_models_proj_current)

# Keep committee averaging and weighted mean ensamble models only.
stk_bee_ensemble_models_proj_current <- subset(stk_bee_ensemble_models_proj_current, grep("EMca|EMwmean", names(stk_bee_ensemble_models_proj_current)))

# Simplify the layer names for plotting conveniences.
names(stk_bee_ensemble_models_proj_current) <-sapply(strsplit(names(stk_bee_ensemble_models_proj_current), "_"), getElement, 2)

# Plot the committee averaging and weighted mean current area1 range ensemble models and write to plots folder.
jpeg("/plots/area1_ensemble_projection.jpeg", type="cairo")
levelplot(stk_bee_ensemble_models_proj_current, main = "A. mellifera current ensemble area1 projection",
          col.regions = colorRampPalette(c("grey90","yellow4","green4"))(100))
dev.off()
print("current area1 range projection complete")

# ------------------------------------------------------------------------------Project current area1 range onto current Hawaii environment-----------------------------------------------------------------------

# We project the "area1 niche" that the species fills in its "area1 range" onto the "Hawaiian environment" to predict suitable habitat in "Hawaii" for that species.

# Load our Hawaii range bioclimatic variables into R.
bio_1 <-raster("haw_raster/bio_1.grd")
bio_2 <-raster("haw_raster/bio_2.grd")
bio_3 <-raster("haw_raster/bio_3.grd")
bio_4 <-raster("haw_raster/bio_4.grd")
bio_5 <-raster("haw_raster/bio_5.grd")
bio_6 <-raster("haw_raster/bio_6.grd")
bio_7 <-raster("haw_raster/bio_7.grd")
bio_8 <-raster("haw_raster/bio_8.grd")
bio_9 <-raster("haw_raster/bio_9.grd")
bio_10 <-raster("haw_raster/bio_10.grd")
bio_11 <-raster("haw_raster/bio_11.grd")
bio_12 <-raster("haw_raster/bio_12.grd")
bio_13 <-raster("haw_raster/bio_13.grd")
bio_14 <-raster("haw_raster/bio_14.grd")
bio_15 <-raster("haw_raster/bio_15.grd")
bio_16 <-raster("haw_raster/bio_16.grd")
bio_17 <-raster("haw_raster/bio_17.grd")
bio_18 <-raster("haw_raster/bio_18.grd")
bio_19 <-raster("haw_raster/bio_19.grd")

# Stack Hawaii range rasters in raster stack.
bioclim_current <- stack(c(bio_2,bio_4,bio_11,bio_12))
print("Hawaii raster stacked")

# For all the models currently implemented, BIOMOD_Projection() function is able to project potential distributions of
# species in other areas, other resolutions or other time scales. We now project all 108 current individual Hawaii range SDMs.
species_predict <- BIOMOD_Projection(modeling.output = bee_models,
                                   new.env = bioclim_current,
                                   proj.name = "species_predict", 
                                   binary.meth = 'TSS', 
                                   do.stack = FALSE)
print("Hawaii BIOMOD_Projecion complete")

# This function use projections of individual models and ensemble models from BIOMOD_EnsembleModeling
# to build an ensemble of species projections over space and time. We now project current Hawaii range ensemble SDMs.
bee_ensemble_models_predict <- BIOMOD_EnsembleForecasting(EM.output = bee_ensemble_models,
                                                          projection.output = species_predict,
                                                          binary.meth = "TSS",
                                                          output.format = ".img",
                                                          do.stack = FALSE)
print("Hawaii BIOMOD_EnsembleForcasting complete")

# Stack the ensemble models.
stk_bee_ensemble_models_predict <- get_predictions(bee_ensemble_models_predict)

# Keep committee averaging and weighted mean ensemble models only.
stk_bee_ensemble_models_predict <- subset(stk_bee_ensemble_models_predict, grep("EMca|EMwmean", names(stk_bee_ensemble_models_predict)))

# Simplify the layer names for plotting conveniences.
names(stk_bee_ensemble_models_predict) <-sapply(strsplit(names(stk_bee_ensemble_models_predict), "_"), getElement, 2)

# Plot the committee averaging and weighted mean current Hawaii range ensemble models and write to plots folder.
jpeg("/plots/area1_hawaii_projection.jpeg", type="cairo")
levelplot(stk_bee_ensemble_models_predict, main = "Current climate Hawaii ensemble projection calibrated witih area1 range",
          col.regions = colorRampPalette(c("grey90","yellow4","green4"))(100))
dev.off()
print("Hawaii range hawaii projection complete")



# ------------------------------------------------------------------------------Project current area1 range onto Future Hawaii climate projections RCP 8.5--------------------------------------------------------

# We project the "area1 niche" that the species fills in its "area1 range" onto the "future Hawaiian environment" to predict suitable habitat in "future Hawaii climate" for that species.
# We use representative concentration pathway 8.5 (RCP8.5) for our future climate scenario.
# Since our future climate scenario variables are the extent of the entire world, we need to clip them to Hawaii before we can use them.

#Load 2070 bioclim variables.
bioclim_world_2070_BC85 <- stack(c(bio_2 = "2070_8.5_ensemble/bio_2.tif",
                                   bio_4 = "2070_8.5_ensemble/bio_4.tif",
                                   bio_11 = "2070_8.5_ensemble/bio_11.tif",
                                   bio_12 = "2070_8.5_ensemble/bio_12.tif"))
print("2070 8.5 rasters loaded")

# Load a Hawaii state shapefile to clip bioclim raster. 
hawaii <- readOGR("/Users/JTAdmin/Desktop/policaris/new_coast_n83/hawaii.shp")

#crop of our area (Hawaii)
bioclim_2070_BC85 <- crop(bioclim_world_2070_BC85, hawaii)
bioclim_2070_BC85 <- mask(bioclim_2070_BC85, hawaii)
bioclim_2070_BC85 <- stack(bioclim_2070_BC85)
print("2070 8.5 rasters clipped to Hawaii shapefile")

# For all the models currently implemented, BIOMOD_Projection() function is able to project potential distributions of
# species in other areas, other resolutions or other time scales. We now project all 108 2070 RCP8.5 individual Hawaii range SDMs.
bee_models_proj_2070_BC85 <- BIOMOD_Projection(modeling.output = bee_models,
                                               new.env = bioclim_2070_BC85,
                                               proj.name = "2070_BC85",
                                               binary.meth = "TSS",
                                               do.stack = FALSE)
print("2070 8.5 BIOMOD_Projection complete")

# This function use projections of individual models and ensemble models from BIOMOD_EnsembleModeling
# to build an ensemble of species projections over space and time. We now project 2070 RCP8.5 Hawaii range ensemble SDMs.
bee_ensemble_models_proj_2070_BC85 <- BIOMOD_EnsembleForecasting(EM.output = bee_ensemble_models,
                                                                 projection.output = bee_models_proj_2070_BC85,
                                                                 binary.meth = "TSS",
                                                                 output.format = ".img",
                                                                 do.stack = FALSE)
print("2070 8.5 BIOMOD_EnsembleForcasting complete")

# Stack the ensemble models.
stk_bee_ef_2070_BC85 <- get_predictions(bee_ensemble_models_proj_2070_BC85)

# Keep committee averaging and weighted mean projections only ensemble models.
stk_bee_ef_2070_BC85 <- subset(stk_bee_ef_2070_BC85, grep("EMca|EMwmean", names(stk_bee_ef_2070_BC85)))

# Simplify the layer names for plotting conveniences.
names(stk_bee_ef_2070_BC85) <-sapply(strsplit(names(stk_bee_ef_2070_BC85), "_"), getElement, 2)

# Plot the committee averaging and weighted mean 2070 8.5 Hawaii range ensemble models and write to plots folder.
jpeg("/plots/area1_hawaii_projection_2070_85.jpeg", type="cairo")
levelplot(stk_bee_ef_2070_BC85, main = "2070 RCP 8.5 Hawaii ensemble projection calibrated witih area1 range",
          col.regions = colorRampPalette(c("grey90","yellow4","green4"))(100))
dev.off()

# ------------------------------------------------------------------------------RCP 8.5 Species range change-------------------------------------------------------------------------------------------------------

# BIOMOD_RangeSize(): This function allows to estimate the proportion and relative number of pixels (or habitat) lost, gained
# and stable for the time slice considered in species-climate modeling under future scenarios.

# Load and stack binary projections for current, 2050, and 2070.
bee_bin_proj_current <- stack(c(ca = "/species/proj_species_predict/individual_projections/name.model_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.img",
                                wm = "/species/proj_species_predict/individual_projections/name.model_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.img"))
bee_bin_proj_2070_BC85 <- stack(c(ca = "/species/proj_2070_BC85/individual_projections/name.model_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.img",
                                  wm = "/species/proj_2070_BC85/individual_projections/name.model_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.img"))
print("binary projections stacked")

# Calculate the range change from current current to 2070.
# SRC current -> 2070
SRC_current_2070_BC85 <- BIOMOD_RangeSize(bee_bin_proj_current, bee_bin_proj_2070_BC85)
SRC_current_2070_BC85$Compt.By.Models

# Plot the predicted changes.
bee_src_map <- stack(SRC_current_2070_BC85$Diff.By.Pixel)
names(bee_src_map) <- c("ca cur-2070", "wm cur-2070")
my.at <- seq(-2.5,1.5,1)
myColorkey <- list(at = my.at, labels = list(labels = c("lost","pres","abs","gain"),at = my.at[-1]-0.5))

# Write plot to plots folder.
jpeg("/plots/species_range_change_85.jpeg", type="cairo")
rasterVis::levelplot(bee_src_map, main = "RCP 8.5 range change", colorkey = myColorkey, layout = c(1,2))
dev.off()
print("8.5 range change plot complete")

print("Hawaii non-native bee study using BIOMOD2 complete")

 

