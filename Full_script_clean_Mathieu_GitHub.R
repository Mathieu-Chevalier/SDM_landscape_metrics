## clean up memory
rm(list=ls()); gc()

### Load libraries
require(rgdal)
require(raster)
require(plyr)
require(dplyr)
require(usdm)
require(viridis)
require(biomod2)
require(scales)
require(ggside)
require(reshape2)
require(factoextra)
require(VoCC)
require(landscapetools)
require(landscapemetrics) 
require(ggplot2)
require(RColorBrewer)
require(ggpubr)
require(rio)
require(tidyr)
require(gridExtra)
require(ade4)
require(ggpmisc)
require(gbm)

### Set working directory
workingDirectory <- "Y:/Objective 4 - Population_Modelling/SDM_new/Mathieu/"
setwd(workingDirectory)

### Define extent
EXT <- extent(-16, 5.014622, 32, 61)

####################################################
################### 1- Load data ################### 
####################################################

### Load data
observations_points <-readOGR("Data/observations", stringsAsFactors=FALSE)
master_raster <- crop(raster("Data/master_rasters/output/master_raster_final.tif"), EXT) # resolution of 5 arcmin (i.e. ~10km at the equator)

### Point selection
observations_points <- observations_points[which(observations_points$PRES_ABS %in% "P"), ] # Appropriate Presence data (no need for absences)
observations_points <- observations_points[which(observations_points$YEAR >= 2000), ] # Years between 2000 and 2019
observations_points <- observations_points[which(observations_points$POS_ACC %in% c("< 10m","< 100m","< 1km")),] # Position accuracy
observations_points <- observations_points[which(observations_points$TIDAL_ZONE == "Intertidal"),] # Select points in tidal zone 

### Spatial selection
observations_points$PRES_ABS <- 1
Rast.obs <- rasterize(observations_points, master_raster, observations_points$PRES_ABS, fun="max")
r_max <- mask(Rast.obs, master_raster, maskvalue=NA)
spObservationPoints <- rasterToPoints(r_max, spatial=TRUE)
proj4string(spObservationPoints) <- CRS("+init=epsg:4326")
names(spObservationPoints)[1] <- "PRES_ABS"
saveRDS(spObservationPoints, file="Data/observations/selected_presence_absence.Rdata")

###########################################################
################### 2- Prepare predictors #################
###########################################################

#######################################################################
#--------- Raster operations for present and future climatic conditions

master_raster <- crop(raster("Data/master_rasters/output/master_raster_final.tif"), EXT) 
selec.time <- c("present", "future")

for(time.step in 1:length(selec.time)){
  
  tVariables <- import(paste0("Data/config_files_predictors/1_predictors_preparation_", selec.time[time.step], ".csv"))
  
  for (i in 1:nrow(tVariables)){
    
    # Get some settings
    variableName <- tVariables[i,"shortName"]
    variableDataset <- tVariables[i,"dataset"]
    variableFolder <- tVariables[i,"folder"]
    variableInputFile <- tVariables[i,"input_file"]
    variableOutputFile <- tVariables[i,"output_file"]
    requires_resampling <- tVariables[i,"requires_resampling"]==1
    number_of_fill_gaps_iterations <- tVariables[i,"fill_gaps_iterations"]
    
    # Load raster
    raster_variable <- raster(paste0("Data/", variableFolder, "/", variableInputFile))
    raster_variable <- crop(raster_variable, master_raster)
    
    # Specific treatments for some variables
    if(variableDataset=="WorldClim")if(grepl("air_temp", variableName)) raster_variable <- raster_variable/10 # divide by 10 for terrestrial temperature
    if(variableDataset=="fetch_burrows_100m") raster_variable <- aggregate(raster_variable, fact=90, fun=mean) # upscale resolution
    if(requires_resampling) raster_variable <- resample(raster_variable, master_raster, method="ngb") 
    if (number_of_fill_gaps_iterations!=0){
      for (j in 1:number_of_fill_gaps_iterations){
        raster_variable <- focal(raster_variable, w=matrix(1, nrow=3, ncol=3), fun=mean, NAonly=TRUE, na.rm=TRUE) # focal to replace missing values along the coast
      }
    }
    
    raster_variable <- mask(raster_variable, master_raster, maskvalue=NA)
    writeRaster(raster_variable, paste0("Data/rasters/ready_for_modeling/", selec.time[time.step], "/", variableOutputFile), overwrite=TRUE)
    
    print(i)
    
  }
}

###################################
#--------- Extract predictor values

rasters.config <- import("Data/config_files_predictors/2_predictor_values_extraction.csv")
new.Names <- rasters.config[,"shortName"]
rastlist <- rasters.config[,3]
allrasters <- stack(paste0("Data/rasters/ready_for_modeling/present/", rastlist))

observations_points <- readRDS("Data/observations/selected_presence_absence.Rdata")
output <- raster::extract(allrasters, observations_points)
colnames(output) <- new.Names
output <- cbind(observations_points, output)
saveRDS(output, file="Data/observations/selected_presence_absence_with_predictors.Rdata")

###################################
#--------- Check for collinearity

### Load presence data
Presence <- readRDS("Data/observations/selected_presence_absence_with_predictors.Rdata") 
Presence <- cbind(coordinates(Presence), Presence@data)

### Discard some variables we don't want (based on expert knowledge)
cor(Presence[,-c(1:3, 4, 10, 15, 14, 5, 8, 11)], use="na.or.complete")

### Compute correlations and remove too collinear variables based on VIF
v1 <- vifstep(Presence[,-c(1:3, 4, 10, 15, 14, 5, 8, 11)], th=7)
Presence.pred <- exclude(Presence, v1)
Presence.pred <- Presence.pred[,-7] # Exclude current velocity cause it has a low influence and there are no future values
saveRDS(Presence.pred, file="Data/observations/selected_presence_absence_with_predictors_selected.Rdata")

############################################################################
#--------- Measure climatic exposition between current and future conditions

rasters.config <- import("Data/config_files_predictors/2_predictor_values_extraction.csv")
Presence.pred <- readRDS("Data/observations/selected_presence_absence_with_predictors_selected.Rdata")

### Get list of current and future climatic layers
list.pres <- list.files("Data/rasters/ready_for_modeling/present/")
list.fut <- list.files("Data/rasters/ready_for_modeling/future/")

### Focus on 2050 and RCP45
list.fut <- list.fut[grep("2050_RCP45", list.fut)]
tmp.list.fut <- gsub("_2050_RCP45", "", list.fut)

### Keep only variables that change in the future
list.pres <- list.pres[which(list.pres %in% tmp.list.fut)]

### Keep only variables that will be included in the modeling framework
keep <- rasters.config[which(rasters.config$shortName %in% colnames(Presence.pred)), "raster"]

### Load current and future climatic data
rast.pres <- stack(paste0("Data/rasters/ready_for_modeling/present/", list.pres[which(list.pres %in% keep)]))
rast.fut <- stack(paste0("Data/rasters/ready_for_modeling/future/", list.fut[which(list.pres %in% keep)]))

### Build a PCA using climatic conditions over the study area for the current period
tmp.df <- na.omit(as.data.frame(rast.pres, xy=T))
df.pres <- na.omit(as.data.frame(rast.pres))
df.fut <- na.omit(as.data.frame(rast.fut))
PCA <- dudi.pca(df.pres, scannf=F, nf=2)
PCA$eig/sum(PCA$eig) # percentage of variance explained by each axis (56% and 25% for the first and second axes, respectively)
proj.pca <- suprow(PCA, df.fut)

### Plot the PCA
tmp.plot <- fviz_pca_biplot(PCA, label ="var", alpha.ind = 0.4, col.ind="contrib")+
  geom_point(data=proj.pca$lisup, aes(x=Axis1, y=Axis2), col="darkgreen", alpha=0.2, shape=4)+
  scale_colour_gradient("Contribution", low="blue", high="red")
tiff("Outputs/PCA_exposure.tiff", width=2500, height=2500, res=250)
print(tmp.plot)
dev.off()

### Measure change between current and future conditions using Euclidean distance-based metric within the two-dimensional space defined by the PCA
euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))
Dist.vec <- numeric()
for(i in 1:nrow(tmp.df)) Dist.vec[i] <- euc.dist(PCA$li[i,], proj.pca$lisup[i,])

### Create a raster out of the metric
rast.chge <- rasterize(tmp.df[,1:2],  rast.pres, Dist.vec)

### Plot the exposure
tmp <- cbind(tmp.df[,1:2], Dist.vec)
colnames(tmp)[3] <- "Exposure"
tmp.plot <- ggplot(tmp, aes(x=x, y=y, fill=Exposure))+
  geom_raster()+
  scale_fill_gradient(low="green", high="red")+
  theme_void()+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),strip.text.x = element_text(size=8), strip.background = element_blank())+
  coord_equal()
tiff("Outputs/Exposure.tiff", width=2500, height=2500, res=250)
print(tmp.plot)
dev.off()

### Save all rasters used to compute future changes
rast.pres <- crop(rast.pres, extent(rast.chge))
rast.fut <- crop(rast.fut, extent(rast.chge))
all.rast <- list(rast.pres, rast.fut, rast.chge)
names(all.rast) <- c("Present", "Future_2050_RCP45", "Change(distance)")
saveRDS(all.rast, file="Outputs/clim_rasters_pres-fut-chge.Rdata")

### Create rasters of absolute relative differences between future and present for each variable (to check which variable is changing the most and where)
norm.airT.pres <- (df.pres[,1] - min(c(df.pres[,1], df.fut[,1])))/(max(c(df.pres[,1], df.fut[,1]))-min(c(df.pres[,1], df.fut[,1])))
norm.airT.fut <- (df.fut[,1] - min(c(df.pres[,1], df.fut[,1])))/(max(c(df.pres[,1], df.fut[,1]))-min(c(df.pres[,1], df.fut[,1])))
chge.airT <- abs((((norm.airT.fut+1) - (norm.airT.pres+1))/(norm.airT.pres+1))*100)
rast.chge.airT <- rasterize(tmp.df[,1:2],  rast.pres, chge.airT)

norm.NOC.pres <- (df.pres[,2] - min(c(df.pres[,2], df.fut[,2])))/(max(c(df.pres[,2], df.fut[,2]))-min(c(df.pres[,2], df.fut[,2])))
norm.NOC.fut <- (df.fut[,2] - min(c(df.pres[,2], df.fut[,2])))/(max(c(df.pres[,2], df.fut[,2]))-min(c(df.pres[,2], df.fut[,2])))
chge.NOC <- abs((((norm.NOC.fut+1) - (norm.NOC.pres+1))/(norm.NOC.pres+1))*100)
rast.chge.NOC <- rasterize(tmp.df[,1:2],  rast.pres, chge.NOC)

norm.sal.pres <- (df.pres[,3] - min(c(df.pres[,3], df.fut[,3])))/(max(c(df.pres[,3], df.fut[,3]))-min(c(df.pres[,3], df.fut[,3])))
norm.sal.fut <- (df.fut[,3] - min(c(df.pres[,3], df.fut[,3])))/(max(c(df.pres[,3], df.fut[,3]))-min(c(df.pres[,3], df.fut[,3])))
chge.sal <- abs((((norm.sal.fut+1) - (norm.sal.pres+1))/(norm.sal.pres+1))*100)
rast.chge.sal <- rasterize(tmp.df[,1:2],  rast.pres, chge.sal)

norm.SST.pres <- (df.pres[,4] - min(c(df.pres[,4], df.fut[,4])))/(max(c(df.pres[,4], df.fut[,4]))-min(c(df.pres[,4], df.fut[,4])))
norm.SST.fut <- (df.fut[,4] - min(c(df.pres[,4], df.fut[,4])))/(max(c(df.pres[,4], df.fut[,4]))-min(c(df.pres[,4], df.fut[,4])))
chge.SST <- abs((((norm.SST.fut+1) - (norm.SST.pres+1))/(norm.SST.pres+1))*100)
rast.chge.SST <- rasterize(tmp.df[,1:2],  rast.pres, chge.SST)

rast.chge.var <- stack(rast.chge.airT, rast.chge.NOC, rast.chge.sal, rast.chge.SST)
names(rast.chge.var) <- colnames(df.pres)
plot(rast.chge.var, zlim=c(0,12))
writeRaster(rast.chge.var, filename = "outputs/raster_files/raster_change_each_variable.tif", overwrite=T)

######################################################
################### 3- Model fitting #################
######################################################

#####################################
#--------- Formatting data for biomod

### Get presence data
Presence <- readRDS("Data/observations/selected_presence_absence_with_predictors.Rdata") 
Presence <- cbind(coordinates(Presence), Presence@data)
Presence.pred <- readRDS("Data/observations/selected_presence_absence_with_predictors_selected.Rdata")
df.pres <- cbind(Presence[,1:3], Presence.pred) 

### Get environmental values in non-presence pixels
rasters.config <- import("Data/config_files_predictors/2_predictor_values_extraction.csv")
list.pres <- list.files("Data/rasters/ready_for_modeling/present/")
tmp.pred <- stack(paste0("Data/rasters/ready_for_modeling/present/", list.pres))
mat <- match(colnames(Presence.pred), rasters.config[,1])
rast.keep <- gsub(".tif", "", rasters.config[mat,3])
tmp.pred <- tmp.pred[[rast.keep]]
Predictors <- crop(tmp.pred, EXT)
names(Predictors) <- colnames(Presence.pred)

### Generate absence data
TMP <- raster::extract(Predictors, df.pres[,1:2], cellnumbers=T)
Predictors.abs <- Predictors
Predictors.abs[TMP[,1]] <- NA
df.abs <- na.omit(as.data.frame(Predictors.abs, xy=T))
df.abs$PRES_ABS <- NA
df.abs <- df.abs[,c(1:2, 9, 3:8)]

### Combine presence and absence data
colnames(df.pres)[1:2] <- c("x", "y")
df.pres.abs <- rbind(df.pres, df.abs)

### Run the models
Reehab.data <- BIOMOD_FormatingData(resp.var = df.pres.abs$PRES_ABS,
                                    expl.var = Predictors,
                                    resp.xy = df.pres.abs[,c("x", "y")],
                                    resp.name = "PRES.ABS",
                                    PA.nb.rep = 10,
                                    PA.nb.absences = length(which(Z$PRES_ABS==1)), # Equal number of presence and absence
                                    PA.strategy = "random")
saveRDS(Reehab.data, file="Outputs/Biomod_data.Rdata")

###########################
#--------- Biomod modeling

### Set options
reehab.mod.options <- BIOMOD_ModelingOptions(GLM = list(type = 'quadratic'),
                                             RF = list(n.trees = 500),
                                             GBM = list(n.trees = 500),
                                             GAM = list(algo = 'GAM_mgcv'))

### Run models
reehab.models <- BIOMOD_Modeling(data = Reehab.data,
                                 models = c("GLM", "GBM", "RF", "GAM"), # GBM is equal to BRT
                                 models.options = reehab.mod.options,
                                 NbRunEval = 10,
                                 DataSplit = 70,
                                 VarImport = 3,
                                 do.full.models = FALSE,
                                 modeling.id = "test")
saveRDS(reehab.models, file="Outputs/Biomod_individual_runs.Rdata")

### Ensemble modeling
reehab.ensembles <- BIOMOD_EnsembleModeling(modeling.output = reehab.models, 
                                            em.by = "all",
                                            eval.metric = "TSS",
                                            eval.metric.quality.threshold = 0.5,
                                            models.eval.meth = c("KAPPA", "TSS", "ROC"),
                                            prob.mean = FALSE,
                                            prob.cv = TRUE,
                                            committee.averaging = FALSE,
                                            prob.mean.weight = TRUE,
                                            VarImport = 0)
saveRDS(reehab.ensembles, file="Outputs/Biomod_ensemble_run.Rdata")

######################################################
################### 4- Model projection ##############
######################################################

reehab.models <- readRDS("Outputs/Biomod_individual_runs.Rdata")
reehab.ensembles <- readRDS("Outputs/Biomod_ensemble_run.Rdata")

#####################################
#-------- Current climatic conditions

### Individual projections
reehab.models.proj.current <- BIOMOD_Projection(modeling.output = reehab.models,
                                                new.env = Predictors,
                                                proj.name = "current",
                                                binary.meth = "TSS",
                                                output.format = ".img",
                                                do.stack = FALSE)
saveRDS(reehab.models.proj.current, file="Outputs/Individual_projections_current.Rdata")

### Ensemble projections
reehab.ensemble.proj.current <- BIOMOD_EnsembleForecasting(EM.output = reehab.ensembles,
                                                           projection.output = reehab.models.proj.current,
                                                           binary.meth = "TSS",
                                                           output.format = ".img",
                                                           do.stack = FALSE)
saveRDS(reehab.ensemble.proj.current, file="Outputs/Ensemble_projections_current.Rdata")

#####################################
#--------- Future climatic conditions

### Change predictors that are changing in the future with projection values
Predictors.fut <- Predictors
Predictors.fut[["AT_MIN"]] <- raster("Data/rasters/ready_for_modeling/future/air_temp_min_2050_RCP45.tif")
Predictors.fut[["SAL_MEAN"]] <- raster("Data/rasters/ready_for_modeling/future/sal_mean_2050_RCP45.tif")
Predictors.fut[["HS_P90"]] <- raster("Data/rasters/ready_for_modeling/future/NOC_hs_p90_2050_RCP45.tif")
Predictors.fut[["SST_MAX"]] <- raster("Data/rasters/ready_for_modeling/future/SST_max_2050_RCP45.tif")

### Individual projections
reehab.models.proj.futur <- BIOMOD_Projection(modeling.output = reehab.models,
                                              new.env = Predictors.fut,
                                              proj.name = "futur",
                                              binary.meth = "TSS",
                                              output.format = ".img",
                                              do.stack = FALSE)
saveRDS(reehab.models.proj.futur, file="Outputs/Individual_projections_futur.Rdata")

### Ensemble projections
reehab.ensemble.proj.futur <- BIOMOD_EnsembleForecasting(EM.output = reehab.ensembles,
                                                         projection.output = reehab.models.proj.futur,
                                                         binary.meth = "TSS",
                                                         output.format = ".img",
                                                         do.stack = FALSE)
saveRDS(reehab.ensemble.proj.futur, file="Outputs/Ensemble_projections_futur.Rdata")

##############################################################
################### 5- Process model outputs #################
##############################################################

####################################################
#--------- Model evaluations and variable importance

#-------- Individual models

reehab.models <- readRDS("Outputs/Biomod_individual_runs.Rdata")

### Get evalations
reehab.model.score <- get_evaluations(reehab.models)
apply(reehab.model.score, c(1, 2), sd)
a <- apply(reehab.model.score, c(1, 3, 2), function(x)length(which(x<0.5)))
a[1,,1] # Models we do not keep (TSS below 0.5)

### Plot evaluations
plot1 <- models_scores_graph(reehab.models, by = "models", metrics = c("ROC","TSS"), plot=F) +
  theme_bw() + 
  theme(legend.title=element_blank())
plot2 <- models_scores_graph(reehab.models, by = "cv_run", metrics = c("ROC","TSS"), plot=F) +
  theme_bw() + 
  theme(legend.title=element_blank())
plot3 <- models_scores_graph(reehab.models, by = "data_set", metrics = c("ROC","TSS"), plot=F) +
  theme_bw() + 
  theme(legend.title=element_blank())

svg("Outputs/Evaluation_TSS_ROC.svg", width=16, height=5)
# tiff("Outputs/Evaluation_TSS_ROC.tiff", width=2500, height=2500, res=250)
grid.arrange(plot1, plot2, plot3, ncol=3)
dev.off()

### Variable importance
reehab.models.var.import <- get_variables_importance(reehab.models)
tab.var.imp <- apply(reehab.models.var.import, c(1, 2), mean)
write.table(tab.var.imp, file="Outputs/Variable_importance.txt")

#-------- Evaluation ensemble model

reehab.ensembles <- readRDS("Outputs/Biomod_ensemble_run.Rdata")
ensemble.models.scores <- get_evaluations(reehab.ensembles)

###############################
#--------- Draw response curves

reehab.models <- readRDS("Outputs/Biomod_individual_runs.Rdata")

### Define some settings to plots response curves
pred.df <- get_formal_data(reehab.models, 'expl.var')
var.names <- colnames(pred.df)
means <- apply(pred.df, 2, mean)
tmp.new.data <- data.frame(sapply(means, function(x)rep(x, 100)))

### Get model names
mod.list <- list.files("PRES.ABS/models/test/")
if(length(grep("Data", mod.list)) !=0) mod.list <- mod.list[-grep("Data", mod.list)]

### Define names for algorithms, runs and PA datasets
tmp <- strsplit(mod.list, "_")
algo.names <- unlist(lapply(tmp, function(x)x[4]))
run.names <- unlist(lapply(tmp, function(x)x[3]))
PA.names <- unlist(lapply(tmp, function(x)x[2]))

### Loop to perform predictions
out.all <- NULL # Define an array to store predictions

for(i in 1:length(mod.list)){
  
  load(paste0("PRES.ABS/models/test/", mod.list[i]))
  tmp.mod <- get_formal_model(get(mod.list[i]))
  
  for(j in 1:length(var.names)){
    
    # Make predictions
    foc.var <- pred.df[,var.names[j]]
    new.data <- tmp.new.data
    var.new <- seq(min(foc.var), max(foc.var), length=100)
    new.data[,j] <- var.new
    
    if(algo.names[i] %in% c("GLM", "GAM")){
      tmp.pred <- predict(tmp.mod, newdata=new.data, type="response")
    } else if(algo.names[i] == "GBM"){
      tmp.perf <- gbm.perf(tmp.mod, method = "cv", plot.it = F)
      tmp.pred <- predict(tmp.mod, newdata=new.data, type="response", n.trees=tmp.perf)
    } else tmp.pred <- predict(tmp.mod, newdata=new.data, type="prob")[,2]
    
    # Store results
    tmp1 <- cbind(Occ.prob=tmp.pred, Env.val=var.new)
    tmp2 <- data.frame(cbind(Algorithm=rep(algo.names[i], each=100),
                             RUN=rep(run.names[i], each=100),
                             PA=rep(PA.names[i], each=100),
                             Var.name=var.names[j]))
    
    out.all <- rbind(out.all,cbind(tmp1,tmp2))
    
  }
  
}
saveRDS(out.all, file="Outputs/Individual_predictions_response_curves.Rdata")

### The plot
# sort(apply(tab.var.imp,1,mean))
out.all <-readRDS("Outputs/Individual_predictions_response_curves.Rdata")
out.all$Var.name <- factor(out.all$Var.name, levels = c("FETCH", "SST_MAX", "SAL_MEAN", "AT_MIN", "HS_P90", "TIDE", "CUR_MEAN"))
out.all$Algorithm <- factor(out.all$Algorithm, levels = c("GLM", "GAM", "GBM", "RF"))
out.all[which(out.all$Var.name == "FETCH"), "Env.val"] <- out.all[which(out.all$Var.name == "FETCH"), "Env.val"]/100

tmp.plot <- ggplot(out.all, aes(x=Env.val, y=Occ.prob)) +
  stat_summary(geom="ribbon", fun.data=mean_cl_normal, fill="grey70")+
  stat_summary(geom="line", fun=mean, size=1, col="darkgreen")+
  # geom_rug(sides="b",alpha = 1/2, position = "jitter")+
  facet_grid(Algorithm ~ Var.name, scale="free_x")+
  theme_bw()+
  labs(x="Values",y="Occurence probability")+
  theme_bw()

svg("Outputs/Average_response_curves.svg", width=14, height=12)
# tiff("Outputs/Average_response_curves.tiff", width=3500, height=2500, res=300)
print(tmp.plot)
dev.off()

###############################################################
#--------- Projections and uncertainties of the ensemble model

pal1 <- rev(colorRampPalette(brewer.pal(11,"RdYlGn"))(64))

#------------- Present and future HS and associated uncertainty

### Current projections
raster_files.all <- list.files("PRES.ABS/proj_current/individual_projections/", pattern = ".img")
raster_files.rm.1 <- list.files("PRES.ABS/proj_current/individual_projections/", pattern = ".img.")
raster_files.rm.2 <- list.files("PRES.ABS/proj_current/individual_projections/", pattern = "bin.img")
raster_files.keep <- raster_files.all[-which(raster_files.all %in% c(raster_files.rm.1, raster_files.rm.2))]
raster_files.ensemble <- raster_files.keep[grep("mergedAlgo_mergedRun_mergedData", raster_files.keep)]
stack.rast <- stack(paste0("PRES.ABS/proj_current/individual_projections/", raster_files.ensemble))
stack.rast <- stack.rast/1000
names(stack.rast) <- c("Coefficient of variation", "Weighted mean")
HS.cur <- stack.rast$Weighted.mean
tmp.df.rast <- as.data.frame(stack.rast, xy=T, na.rm=T)
df.pres <- melt(tmp.df.rast, id=c("x", "y"))

### Future projections
raster_files.all <- list.files("PRES.ABS/proj_futur/individual_projections/", pattern = ".img")
raster_files.rm.1 <- list.files("PRES.ABS/proj_futur/individual_projections/", pattern = ".img.")
raster_files.rm.2 <- list.files("PRES.ABS/proj_futur/individual_projections/", pattern = "bin.img")
raster_files.keep <- raster_files.all[-which(raster_files.all %in% c(raster_files.rm.1, raster_files.rm.2))]
raster_files.ensemble <- raster_files.keep[grep("mergedAlgo_mergedRun_mergedData", raster_files.keep)]
stack.rast <- stack(paste0("PRES.ABS/proj_futur/individual_projections/", raster_files.ensemble))
stack.rast <- stack.rast/1000
names(stack.rast) <- c("Coefficient of variation", "Weighted mean")
HS.fut <- stack.rast$Weighted.mean
tmp.df.rast <- as.data.frame(stack.rast, xy=T, na.rm=T)
df.fut <- melt(tmp.df.rast, id=c("x", "y"))

### Plot results
plot1 <- ggplot(df.pres[which(df.pres$variable == "Coefficient.of.variation"),], aes(x=x, y=y, fill=value))+
  geom_raster()+
  scale_fill_gradient2("Coefficient \nof variation", low="white", high="red")+
  theme_void()+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),
        strip.text.x = element_text(size=8),
        strip.background = element_blank())+
  coord_equal()+
  ggtitle("Present")

plot2 <- ggplot(df.pres[which(df.pres$variable == "Weighted.mean"),], aes(x=x, y=y, fill=value))+
  geom_raster()+
  scale_fill_gradientn("Habitat\nsuitability",colours=pal1, lim=c(0,1))+
  theme_void()+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),
        strip.text.x = element_text(size=8),
        strip.background = element_blank())+
  coord_equal()

plot3 <- ggplot(df.fut[which(df.fut$variable == "Coefficient.of.variation"),], aes(x=x, y=y, fill=value))+
  geom_raster()+
  scale_fill_gradient2("Coefficient \nof variation", low="white", high="red")+
  theme_void()+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),
        strip.text.x = element_text(size=8),
        strip.background = element_blank())+
  coord_equal()+
  ggtitle("Futur")

plot4 <- ggplot(df.fut[which(df.fut$variable == "Weighted.mean"),], aes(x=x, y=y, fill=value))+
  geom_raster()+
  scale_fill_gradientn("Habitat\nsuitability",colours=pal1, lim=c(0,1))+
  theme_void()+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),
        strip.text.x = element_text(size=8),
        strip.background = element_blank())+
  coord_equal()

svg("Outputs/Predicted_HS_and_uncertainty.svg", width=9, height=9)
# tiff("Outputs/Predicted_HS_and_uncertainty.tiff", width=2500, height=2500, res=250)
grid.arrange(plot1, plot2, plot3, plot4, ncol=2)
dev.off()

#---------- Differences between the two projections

### Difference in HS
Diff <- df.fut[which(df.fut$variable == "Weighted.mean"),"value"] - df.pres[which(df.pres$variable == "Weighted.mean"),"value"]
df.dif <- cbind(tmp.df.rast[,1:2], Diff)
colnames(df.dif) <- c("x", "y", "Difference")
pal2 <- colorRampPalette(brewer.pal(11,"RdYlGn"))(64)
write.table(df.dif, "Outputs/Table_difference_HS.txt")

tmp.plot <- ggplot(df.dif, aes(x=x, y=y, fill=Difference))+
  geom_raster()+
  scale_fill_gradientn("Change in HS",colours=pal2)+
  theme_void()+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),
        strip.text.x = element_text(size=8),
        strip.background = element_blank())+
  coord_equal()

svg("Outputs/Change_HS.svg", width=9, height=9)
# tiff("Outputs/Change_HS.tiff", width=2500, height=2500, res=250)
print(tmp.plot)
dev.off()

###############################################
#--------- Change in species range

PA.cur <- raster("PRES.ABS/proj_current/individual_projections/PRES.ABS_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData_ROCbin.img")
PA.fut <- raster("PRES.ABS/proj_futur/individual_projections/PRES.ABS_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData_ROCbin.img")
Rg.chge <- BIOMOD_RangeSize(PA.cur, PA.fut)
Rg.chge$Compt.By.Models
tmp <- Rg.chge$Diff.By.Pixel
df.chge <- na.omit(as.data.frame(tmp, xy=T))
colnames(df.chge)[3] <- "value"
df.chge$value <- as.factor(df.chge$value)
levels(df.chge$value) <- c("Lost", "Stable presence", "Stable absence", "Gain")

tmp.plot <- ggplot(df.chge, aes(x=x, y=y, fill=factor(value))) +
  geom_raster()+
  scale_fill_brewer(palette = "Set1") +
  coord_equal()+
  theme_void() +
  labs(fill="")

svg("Outputs/Range_change.svg", width=9, height=9)
# tiff("Outputs/Range_change.tiff", width=2500, height=2500, res=250)
print(tmp.plot)
dev.off()

###############################################
#--------- Range shift & percentages of (1) pixels gained, (2) pixels lost, (3) stable presence and (4) stable absence

PA.cur <- raster("PRES.ABS/proj_current/individual_projections/PRES.ABS_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData_ROCbin.img")
PA.fut <- raster("PRES.ABS/proj_futur/individual_projections/PRES.ABS_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData_ROCbin.img")

### Get data
df <- na.omit(cbind(as.data.frame(PA.cur, xy=T), as.data.frame(PA.fut)))
colnames(df)[3:4] <- c("Current", "Futur")
table(df$Current, df$Futur)

### Compute current leading and trailing edges
cur <- df[which(df$Current==1),]
optimum.cur <- median(cur$y)
qu.cur <- quantile(cur$y, probs = c(0.025, 0.975))
trailing.cur <- qu.cur[1]
leading.cur <- qu.cur[2]

### Compute future leading and trailing edges
fut <- df[which(df$Futur==1),]
optimum.fut <- median(fut$y)
qu.fut <- quantile(fut$y, probs = c(0.025, 0.975))
trailing.fut <- qu.fut[1]
leading.fut <- qu.fut[2]

### Global range change
Nb.pix.pres.cur <- nrow(df[which(df$Current==1),])
Nb.pix.pres.fut <- nrow(df[which(df$Futur==1),])
Nb.pix.range <- abs(((Nb.pix.pres.fut - Nb.pix.pres.cur)/Nb.pix.pres.cur)*100) # ~27.5%

### Percentage of pixels gained
Nb.pix.pres.cur <- nrow(df[which(df$Current==1),])
Nb.pix.abs.cur.pres.fut <- nrow(df[which(df$Current==0 & df$Futur==1),])
Nb.pix.gain <- 100-abs(((Nb.pix.abs.cur.pres.fut - Nb.pix.pres.cur)/Nb.pix.pres.cur)*100) # ~35.8%

### Percentage of pixels lost
Nb.pix.pres.cur <- nrow(df[which(df$Current==1),])
Nb.pix.pres.cur.abs.fut <- nrow(df[which(df$Current==1 & df$Futur==0),])
Nb.pix.lost <- 100 - abs(((Nb.pix.pres.cur - Nb.pix.pres.cur.abs.fut)/Nb.pix.pres.cur)*100) #~8.3%

### Percentage of stable presence pixels 
Nb.pix.pres.cur.pres.fut <- nrow(df[which(df$Current==1 & df$Futur==1),])
Nb.pix.stableP <- 100 - Nb.pix.lost #~91.7%

### Percentage of stable absence pixels 
Nb.pix.abs.cur <- nrow(df[which(df$Current==0),])
Nb.pix.abs.cur.abs.fut <- nrow(df[which(df$Current==0 & df$Futur==0),])
Nb.pix.stableA <- 100 - Nb.pix.gain #~72.5%

### Put all metrics in a table
tab.shift <- as.data.frame(cbind(c(optimum.cur, trailing.cur, leading.cur), c(optimum.fut, trailing.fut, leading.fut)))
colnames(tab.shift) <- c("Current", "Futur")
rownames(tab.shift) <- c("Optimum", "Trailing edge", "Leading edge")
tab.shift$difference <- ((tab.shift$Futur - tab.shift$Current)/tab.shift$Futur)*100
tab.shift <- rbind(tab.shift, c(Nb.pix.pres.cur, Nb.pix.abs.cur.pres.fut, Nb.pix.gain)) # pixel gained
tab.shift <- rbind(tab.shift, c(Nb.pix.pres.cur, Nb.pix.pres.cur.abs.fut, Nb.pix.lost)) # pixel lost
tab.shift <- rbind(tab.shift, c(Nb.pix.pres.cur, Nb.pix.pres.cur.pres.fut, Nb.pix.stableP)) # Stable presence
tab.shift <- rbind(tab.shift, c(Nb.pix.abs.cur, Nb.pix.abs.cur.abs.fut, Nb.pix.stableA)) # Stable absence
tab.shift <- rbind(tab.shift, c(Nb.pix.pres.cur, Nb.pix.pres.fut, Nb.pix.range)) # Global range change
rownames(tab.shift)[4:8] <- c("Pixel gained (%)", "Pixel lost (%)", "Stable presence pixel (%)", "Stable absence pixel (%)", "Range change (%)")

### Write down results
write.table(tab.shift, file="Outputs/Estimated_range_shift.txt")

###############################################
#--------- Write down rasters

PA.cur <- raster("PRES.ABS/proj_current/individual_projections/PRES.ABS_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData_ROCbin.img")
PA.fut <- raster("PRES.ABS/proj_futur/individual_projections/PRES.ABS_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData_ROCbin.img")
HS.dif <- HS.fut-HS.cur
Expo <- readRDS("Outputs/clim_rasters_pres-fut-chge.Rdata")
Expo <- Expo$'Change(distance)'
Rast.all <- stack(Expo, PA.cur, PA.fut, Rg.chge$Diff.By.Pixel, HS.cur, HS.fut, HS.dif)
names(Rast.all) <- c("Climatic.exposure", "PA.current", "PA.futur", "Range.change", "HS.current", "HS.futur", "HS.difference") 
saveRDS(Rast.all, file="Outputs/All.predicted.rasters.Rdata")

dir.create("outputs/raster_files")
writeRaster(Expo, filename = "outputs/raster_files/Climatic.exposure.tif", overwrite=T)
writeRaster(PA.cur, filename = "outputs/raster_files/PA.current.tif", overwrite=T)
writeRaster(PA.fut, filename = "outputs/raster_files/PA.futur.tif", overwrite=T)
writeRaster(Rg.chge$Diff.By.Pixel, filename = "outputs/raster_files/Range.change.tif", overwrite=T)
writeRaster(HS.cur, filename = "outputs/raster_files/HS.current.tif", overwrite=T)
writeRaster(HS.fut, filename = "outputs/raster_files/HS.futur.tif", overwrite=T)
writeRaster(HS.dif, filename = "outputs/raster_files/HS.difference.tif", overwrite=T)

############################################################################
################### 6- Spatial pattern distribution ranges #################
############################################################################

###################################################################################
#--------- Look at changes in spatial patterns between current and future projections

# list_lsm() # The different landscape metrics available

#---------- Check out some landscape metrics

### New projection to give equal area to all pixels
New.proj <- "+proj=laea +lat_0=44.3 +lon_0=-3.2 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

### Get present and future PA data
PA.cur <- raster("PRES.ABS/proj_current/individual_projections/PRES.ABS_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData_ROCbin.img")
PA.fut <- raster("PRES.ABS/proj_futur/individual_projections/PRES.ABS_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData_ROCbin.img")
Both <- stack(PA.cur, PA.fut)
Rast <- projectRaster(Both, method = "ngb", res=5000, crs = New.proj) 

### Compute some landscape metrics
land_metrics <- calculate_lsm(Rast, what = c("patch", "class"))

### Summary of some landscape metrics for presence and absence patches
land_metrics %>%
  filter(level == "class") %>%
  filter(metric %in% c("np", "ca","area_mn","area_sd","enn_mn","enn_sd","clumpy", "cohesion", "ai")) %>%
  spread(metric, value) %>%
  select(-layer,-level, -id) %>%
  mutate_if(is.numeric, round, 1)  %>%
  mutate(Layer=rep(c("Current", "Future"),2))

#---------- Change in patch area for different class  of patches (1-1, 0-0, 1-0, 0-1)

### Compute patch area
tmp <- show_lsm(Rast, what = "lsm_p_area", directions = 8, class = "all", label_lsm = T)
Area.cur <- tmp[[1]]$data
Area.cur <- Area.cur[order(Area.cur$x, Area.cur$y),]
Area.cur <- Area.cur[which(!is.na(Area.cur$id)),]
Area.fut <- tmp[[2]]$data
Area.fut <- Area.fut[order(Area.fut$x, Area.fut$y),]
Area.fut <- Area.fut[which(!is.na(Area.fut$id)),]

### Define transitions
pres <- which(Area.cur$class.get_patches==1 & Area.fut$class.get_patches==1)
abs <- which(Area.cur$class.get_patches==0 & Area.fut$class.get_patches==0)
pres.abs <- which(Area.cur$class.get_patches==1 & Area.fut$class.get_patches==0)
abs.pres <- which(Area.cur$class.get_patches==0 & Area.fut$class.get_patches==1)

### Re-arrange a bit the data
TMP.cur <- c(Area.cur[pres,"value"], Area.cur[abs,"value"], Area.cur[pres.abs,"value"], Area.cur[abs.pres,"value"])
TMP.fut <- c(Area.fut[pres,"value"], Area.fut[abs,"value"], Area.fut[pres.abs,"value"], Area.fut[abs.pres,"value"])
TMP.period <- rep(c("Present", "Futur"), each=length(TMP.cur))
TMP.class <- rep(c(rep("Pres", length(pres)), rep("Abs", length(abs)), rep("pres.abs", length(pres.abs)), rep("abs.pres", length(abs.pres))),2)
TMP.y <- rep(Area.cur[,"y"], 2)
A <- as.data.frame(cbind(c(TMP.cur, TMP.fut), TMP.period, TMP.class))
A$y <- TMP.y
colnames(A) <- c("value", "period", "class", "y")
A$value <- log(as.numeric(A$value))
A$period <- as.factor(A$period)
A$class <- as.factor(A$class)
A$period <- factor(A$period, levels = c("Present", "Futur"))
A$class  <- revalue(A$class, c("Abs"="Absence", "Pres"="Presence", "pres.abs"="Presence -> Absence", "abs.pres"="Absence -> Presence"))

### the plot
tmp.plot <- ggplot(A, aes(y=value, x=class, fill=period))+
  geom_boxplot(position = position_dodge(0.9), outlier.shape = NA)+
  theme_bw()+
  labs(y="Patch area", x="")

tiff("Outputs/Boxplot_patch_area_classes.tiff", width=2500, height=2500, res=250)
print(tmp.plot)
dev.off()

#---------- Combine patch and pixel information into one DF

Rast.expo <- readRDS("Outputs/clim_rasters_pres-fut-chge.Rdata")
Rast.expo <- Rast.expo$'Change(distance)'

### Current and future suitability
HS.cur <- raster("PRES.ABS/proj_current/individual_projections/PRES.ABS_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData.img")
HS.cur <- HS.cur/1000
HS.fut <- raster("PRES.ABS/proj_futur/individual_projections/PRES.ABS_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData.img")
HS.fut <- HS.fut/1000
Rast.HS <- stack(HS.cur, HS.fut)

### Current and future PA
PA.cur <- raster("PRES.ABS/proj_current/individual_projections/PRES.ABS_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData_ROCbin.img")
PA.fut <- raster("PRES.ABS/proj_futur/individual_projections/PRES.ABS_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData_ROCbin.img")
Rast.PA <- stack(PA.cur, PA.fut)

### Change projection to have equal pixel sizes for all rasters
Rast.expo <- projectRaster(Rast.expo, method = "ngb", res=5000, crs = New.proj) 
Rast.HS <- projectRaster(Rast.HS, method = "ngb", res=5000, crs = New.proj) 
Rast.PA <- projectRaster(Rast.PA, method = "ngb", res=5000, crs = New.proj) 

### combine all rasters and make a data frame out of it
Rast <- stack(Rast.expo, Rast.HS, Rast.PA)
names(Rast) <- c("Exposition", "HS.cur", "HS.fut", "PA.cur", "PA.fut")
DF <- na.omit(as.data.frame(Rast, xy=T))
DF <- DF[order(DF$x, DF$y),] # order coordinates 

### Compute landscape metrics 
func.name <- as.data.frame(list_lsm()[,5])[,1]
metric.name <- as.data.frame(list_lsm()[,1])[,1]
out.metric.patch <- NULL
out.metric.other <- NULL

for(i in 1:length(metric.name)){
  
  ### Deal with patch level metrics
  if(grepl("_p_", func.name[i])){
    
    ### Get the metric on PA rasters
    tmp.met <- show_lsm(Rast.PA, what = func.name[i], directions = 8, class = "all", label_lsm = T)
    
    #--- Present
    tmp.pres <- tmp.met[[1]]$data
    tmp1 <- na.omit(tmp.pres[which(tmp.pres$class.get_patches==1), c(1,2,3,9)]) # for the ones
    tmp2 <- na.omit(tmp.pres[which(tmp.pres$class.get_patches==0), c(1,2,3,9)]) # for the zeroes
    df.pres <- rbind(tmp1, tmp2)
    df.pres <- df.pres[order(df.pres$x, df.pres$y),] # order coordinates 
    
    #--- Futur
    tmp.fut <- tmp.met[[2]]$data
    tmp1 <- na.omit(tmp.fut[which(tmp.fut$class.get_patches==1), c(1,2,3,9)])
    tmp2 <- na.omit(tmp.fut[which(tmp.fut$class.get_patches==0), c(1,2,3,9)])
    df.fut <- rbind(tmp1, tmp2)
    df.fut <- df.fut[order(df.fut$x, df.fut$y),]
    
    ### Assemble the metric
    tmp1 <- cbind(df.pres[,4], df.fut[,4])
    foc.name <- paste(c("Cur", "Fut"), metric.name[i], sep="_")
    colnames(tmp1) <- foc.name
    if(i == 1) tmp1 <- cbind(ID.current=df.pres$id, ID.futur=df.fut$id, tmp1)
    out.metric.patch <- cbind(out.metric.patch, tmp1)
    
  } else {   ### Deal with other level metrics
  
    tmp <- get(func.name[i])
    metric.foc <- tmp(Rast.PA)
    out.metric.other <- rbind(out.metric.other, metric.foc)
    
  }

  print(i)
}

### Combine and save metrics
DF.patch <- cbind(DF, out.metric.patch)
saveRDS(DF.patch, file="Outputs/Table.metrics.patch.Rdata") # Save outputs
saveRDS(out.metric.other, file="Outputs/Table.metrics.other.Rdata") # Save outputs

###  Data re-arragement
DF.all <- readRDS("Outputs/Table.metrics.patch.Rdata")
DF.all$HS.dif <- DF.all$HS.fut - DF.all$HS.cur # compute difference in suitability
write.table(DF.all, "Outputs/Patch_level_metrics.table")
DF.resume <- DF.all %>% # Aggregate measures based on current ID of patches
  group_by(ID.current) %>% 
  summarise_all(mean) 
DF.resume <- bind_cols(DF.resume, N=as.numeric(table(DF.all$ID.current)))
DF.resume$PA.cur <- as.factor(DF.resume$PA.cur)
levels(DF.resume$PA.cur) <- c("Absence", "Presence")

### Get back to original coordinates
coordinates(DF.resume) <- c("x", "y")
proj4string(DF.resume) <- crs(Rast) 
DF.resume <- as.data.frame(spTransform(DF.resume, crs(PA.cur))) # WGS 84

#----------- Make some plots/tables

### Changes in class metrics over time
DF.other <- readRDS("Outputs/Table.metrics.other.Rdata")
write.table(DF.other, "Outputs/Lanscape_class_level_metrics.txt")
DF.class <- DF.other[-which(DF.other$level=="landscape"),] 
DF.class$class <- factor(DF.class$class) 
levels(DF.class$class) <- c("Absence", "Presence")
DF.class$layer <- factor(DF.class$layer) 
levels(DF.class$layer) <- c("Current", "Future")
tmp.plot <- ggplot(DF.class, aes(x=factor(layer), y=value, color=factor(class), group=factor(class)))+ 
  geom_point()+
  geom_line()+
  theme_bw()+
  facet_wrap(~metric, scales="free")+
  labs(y="Metric value", x="Time", color="Patch status")

tiff("Outputs/temporal_change_class_metrics.tiff", width=3200, height=2000, res=250)
print(tmp.plot)
dev.off()

### proportional changes in class and landscape metrics
tab.prop <- NULL
DF.other$pst <- paste(DF.other$level, DF.other$metric, sep="_")
lev.metric <- levels(factor(DF.other$pst))
for(i in 1:length(lev.metric)){
  tmp <- DF.other[which(DF.other$pst == lev.metric[i]),]
  if(levels(factor(tmp$level)) == "class"){
    tmp0 <- pull(tmp[which(tmp$class==0), "value"])
    prop.chge.0 <- ((tmp0[2]-tmp0[1])/tmp0[1])*100
    tmp1 <- pull(tmp[which(tmp$class==1), "value"])
    prop.chge.1 <- ((tmp1[2]-tmp1[1])/tmp1[1])*100
    tab.prop <- rbind(tab.prop, cbind(metric=lev.metric[i], class=c(0,1), prop.chge=c(prop.chge.0, prop.chge.1)))
  } else {
    tmp.glob <- pull(tmp, "value")
    prop.chge.glob <- ((tmp.glob[2]-tmp.glob[1])/tmp.glob[1])*100
    tab.prop <- rbind(tab.prop, cbind(metric=lev.metric[i], class="global", prop.chge=prop.chge.glob))
  }
}
write.table(tab.prop, "Outputs/Proportional_change_landscape_metrics.txt")

### Relationship with latitude (all patches)
Pred.tab <- NULL
mod <- lm(y ~ HS.dif + I(HS.dif^2), weights=N, data=DF.resume)
new <- data.frame(HS.dif = seq(min(DF.resume$HS.dif), max(DF.resume$HS.dif), length.out=100))
pred <- predict(mod, new, se.fit=T)
Pred.tab <- as.data.frame(pred$fit)
Pred.tab$pred.min <- pred$fit - 1.96*pred$se.fit
Pred.tab$pred.max <- pred$fit + 1.96*pred$se.fit
colnames(Pred.tab)[1] <- "pred"
Pred.tab$HS.dif <- new$HS.dif
Pred.tab$pred.min <- ifelse(Pred.tab$pred.min<30, 30, Pred.tab$pred.min)

tmp.plot <- ggplot(DF.resume, aes(x=HS.dif, y=y, size=sqrt(N)))+ 
  geom_point(data=DF.resume[which(DF.resume$PA.cur=="Absence"),], alpha=0.6, col="#d95f02")+
  geom_point(data=DF.resume[which(DF.resume$PA.cur=="Presence"),], alpha=0.6, col="#1b9e77")+
  geom_ribbon(data=Pred.tab, aes(x=HS.dif, ymin=pred.min, ymax=pred.max), fill="gray", inherit.aes=F, alpha=0.5) +
  geom_line(data=Pred.tab, aes(x=HS.dif, y=pred), col="blue", inherit.aes=F, size=1) +
  stat_poly_eq(aes(label = paste(..adj.rr.label.., ..p.value.label.., sep = "~~~"), weight = N), 
               label.x.npc = "right", label.y.npc = "top",
               formula = y~x+I(x^2), parse = TRUE, size = 5)+
  theme_bw()+
  theme(legend.position = "none", axis.text = element_text(size=12), axis.title = element_text(size=14))+
  geom_hline(yintercept = 46.2, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  labs(y="Latitude", x="Difference in habitat suitability")

svg("Outputs/Relationship_latitude_avgHSdif_based_on_patch_identity.svg", width=9, height=9)
# tiff("Outputs/Relationship_latitude_avgHSdif_based_on_patch_identity.tiff", width=2200, height=1800, res=250)
print(tmp.plot)
dev.off()

### Map a subset of patches along the latitudinal gradient 
dif.HS <- read.table("Outputs/Table_difference_HS.txt", h=T)
tmp <- DF.all 
tmp$Cur_area <- log(tmp$Cur_area)
tmp <- tmp[which(tmp$Cur_area>13),]
cols <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99")
map_highlight_patch <- ggplot(data=tmp, aes(x=x, y=y, color=factor(ID.current)))+ 
  geom_point(size=0.1)+
  theme_void()+
  coord_equal()+
  scale_color_manual(values=cols)+
  labs(color="Patch ID")+
  guides(colour = guide_legend(override.aes = list(size=10)))
r <- rasterize(SpatialPoints(tmp[,1:2]), Rast.PA, field=tmp$ID.current)
writeRaster(r, filename=file.path("Outputs/Raster_patch_ID.tif"), overwrite=TRUE)

svg("Outputs/Map_11_selected_patches.svg", width=9, height=9)
# tiff("Outputs/Map_11_selected_patches.tiff", width=2200, height=1800, res=250)
print(map_highlight_patch)
dev.off()

Z1 <- DF.resume[which(DF.resume$ID.current %in% tmp$ID.current),]
new_plot_reg <- tmp.plot + geom_point(data=Z1[which(Z1$PA.cur=="Absence"),], alpha=0.6, fill="#d95f02", shape=21, col="Black", stroke=2)+
  geom_point(data=Z1[which(Z1$PA.cur=="Presence"),], alpha=0.6, fill="#1b9e77", col="Black", stroke=2)+
  geom_text(data=Z1, aes(label=Z1$ID.current, y=y+1.3), size=5)

svg("Outputs/Relationship_latitude_avgHSdif_based_on_patch_identity_highlighted_pixels.svg", width=9, height=9)
print(new_plot_reg)
dev.off()

#----------- Identify points of particular interest to exemplify what happens within the species range

dif.HS <- read.table("Outputs/Table_difference_HS.txt", h=T)
tmp <- DF.all 

#########################
# Example of colonization
#########################

### Process data
tmp0 <- tmp[which(tmp$PA.cur == 0),]
z <- tapply(tmp0$ID.futur, tmp0$ID.current, function(x)nlevels(factor(x)))
id.patch <- names(which(z!=1))
tmp0 <- tmp0[which(tmp0$ID.current %in% id.patch),]
z <- tapply(tmp0$ID.futur, tmp0$ID.current, function(x)nlevels(factor(x)))

tmp1 <- subset(tmp0, ID.current==2) 
tmp1$index <- paste(tmp1$x, tmp1$y, sep="_") # necessary to remove pixels outside the selected one
coordinates(tmp1) <- c("x", "y")
proj4string(tmp1) <- crs(Rast)
tmp1 <- as.data.frame(spTransform(tmp1, crs(PA.cur))) # WGS 84

tmp2 <- tmp0[which(tmp0$ID.futur %in% levels(factor(tmp1$ID.futur))),]
tmp2$index <- paste(tmp2$x, tmp2$y, sep="_")
tmp2 <- tmp2[which(tmp2$index %in% tmp1$index),] 
coordinates(tmp2) <- c("x", "y")
proj4string(tmp2) <- crs(Rast)
tmp2 <- as.data.frame(spTransform(tmp2, crs(PA.cur))) # WGS 84

### Write down rasters
rast1 <- rasterize(SpatialPoints(tmp1[c("x", "y")]), PA.cur, field=tmp1$PA.cur)
rast2 <- rasterize(SpatialPoints(tmp2[c("x", "y")]), PA.cur, field=tmp2$PA.fut)
stack.rast <- stack(rast1, rast2)
names(stack.rast) <- c("Current", "Futur")
writeRaster(stack.rast, filename = "outputs/raster_files/Raster_colonisation_no_merge.tif", overwrite=T)

### Select only pixels transitioning for plotting purposes
tmp2 <- tmp2[which(tmp2$PA.fut == 1),]

plot.colo.1 <- ggplot(dif.HS, aes(x=x, y=y))+ #
  geom_raster(alpha=0.4, fill="grey")+
  geom_point(data=tmp1, aes(x=x, y=y, fill=NULL, color="#EDF285"), size=0.1)+
  geom_point(data=tmp2, aes(x=x, y=y, fill=NULL, color="#612FE0"), size=0.1)+
  theme_void()+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),
        strip.text.x = element_text(size=8),
        strip.background = element_blank())+
  coord_equal()+
  scale_color_manual(values = c("#612FE0", "#EDF285"), labels=c("A->P", "A->A"))+
  labs(color="Patch transition")+
  guides(color = guide_legend(override.aes = list(size=8)))

svg("Outputs/Example_colonisation_no_merge.svg", width=9, height=9)
print(plot.colo.1)
dev.off()

### Compute landscape metrics for this patch
tmp1 <- subset(tmp0, ID.current==2)
tmp1$index <- paste(tmp1$x, tmp1$y, sep="_")

tmp2 <- tmp0[which(tmp0$ID.futur %in% levels(factor(tmp1$ID.futur))),]
tmp2$index <- paste(tmp2$x, tmp2$y, sep="_")
tmp2 <- tmp2[which(tmp2$index %in% tmp1$index),]

rast1 <- rasterize(SpatialPoints(tmp1[c("x", "y")]), Rast.PA, field=tmp1$PA.cur)
rast2 <- rasterize(SpatialPoints(tmp2[c("x", "y")]), Rast.PA, field=tmp2$PA.fut)
stack.rast <- stack(rast1, rast2)

func.name <- as.data.frame(list_lsm()[,5])[,1]
metric.name <- as.data.frame(list_lsm()[,1])[,1]
out.metric.other <- NULL
for(i in 1:length(metric.name)){
  func <- get(func.name[i])
  metric.foc <- func(stack.rast)
  out.metric.other <- rbind(out.metric.other, metric.foc)
}
write.table(out.metric.other, "Outputs/Landscape_metrics_colonisation_example.txt")

##########################################
# Example of expansion of occurrence patch
##########################################

### Process data
tmp0 <- tmp[which(tmp$PA.fut == 1),]
z <- tapply(tmp0$ID.current, tmp0$ID.futur, function(x)nlevels(factor(x)))
id.patch <- names(which(z!=1))
tmp0 <- tmp0[which(tmp0$ID.futur %in% id.patch),]
z <- tapply(tmp0$ID.current, tmp0$ID.futur, function(x)nlevels(factor(x)))

tmp1 <- subset(tmp0, ID.futur==79) 
tmp1$index <- paste(tmp1$x, tmp1$y, sep="_")
coordinates(tmp1) <- c("x", "y")
proj4string(tmp1) <- crs(Rast)
tmp1 <- as.data.frame(spTransform(tmp1, crs(PA.cur))) # WGS 84

tmp2 <- tmp0[which(tmp0$ID.current %in% levels(factor(tmp1$ID.current))),]
tmp2$index <- paste(tmp2$x, tmp2$y, sep="_")
tmp2 <- tmp2[which(tmp2$index %in% tmp1$index),]
coordinates(tmp2) <- c("x", "y")
proj4string(tmp2) <- crs(Rast)
tmp2 <- as.data.frame(spTransform(tmp2, crs(PA.cur))) # WGS 84

### Write down rasters
rast1 <- rasterize(SpatialPoints(tmp1[c("x", "y")]), PA.cur, field=tmp1$PA.fut) # futur
rast2 <- rasterize(SpatialPoints(tmp2[c("x", "y")]), PA.cur, field=tmp2$PA.cur) # current
stack.rast <- stack(rast2, rast1)
names(stack.rast) <- c("Current", "Futur")
writeRaster(stack.rast, filename = "outputs/raster_files/Raster_expansion_occ_patch.tif", overwrite=T)

### Select only pixels transitioning for plotting purposes
tmp2 <- tmp2[which(tmp2$PA.cur == 0),] 

plot.exp.occ <- ggplot(dif.HS, aes(x=x, y=y))+
  geom_raster(alpha=0.2, fill="grey")+
  geom_point(data=tmp1, aes(x=x, y=y, fill=NULL, color="#612FE0"), size=0.1)+
  geom_point(data=tmp2, aes(x=x, y=y, fill=NULL, color="#24E0A2"), size=0.1)+
  theme_void()+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),
        strip.text.x = element_text(size=8),
        strip.background = element_blank())+
  coord_equal()+
  scale_color_manual(values = c("#612FE0", "#24E0A2"), labels=c("A->P", "P->P"))+
  labs(color="Patch transition")+
  guides(color = guide_legend(override.aes = list(size=8)))

svg("Outputs/Example_expansion_occ_patch.svg", width=9, height=9)
print(plot.exp.occ)
dev.off()

### Compute landscape metrics for this patch
tmp1 <- subset(tmp0, ID.futur==79) 
tmp1$index <- paste(tmp1$x, tmp1$y, sep="_")

tmp2 <- tmp0[which(tmp0$ID.current %in% levels(factor(tmp1$ID.current))),]
tmp2$index <- paste(tmp2$x, tmp2$y, sep="_")
tmp2 <- tmp2[which(tmp2$index %in% tmp1$index),]

rast1 <- rasterize(SpatialPoints(tmp1[c("x", "y")]), Rast.PA, field=tmp1$PA.fut)
rast2 <- rasterize(SpatialPoints(tmp2[c("x", "y")]), Rast.PA, field=tmp2$PA.cur)
stack.rast <- stack(rast2, rast1)

func.name <- as.data.frame(list_lsm()[,5])[,1]
metric.name <- as.data.frame(list_lsm()[,1])[,1]
out.metric.other <- NULL
for(i in 1:length(metric.name)){
  func <- get(func.name[i])
  metric.foc <- func(stack.rast)
  out.metric.other <- rbind(out.metric.other, metric.foc)
}
write.table(out.metric.other, "Outputs/Landscape_metrics_expansion.occ_example.txt")

#########################################
# Example of expansion of absence patches
#########################################

### Process data
tmp0 <- tmp[which(tmp$PA.fut == 0),]
z <- tapply(tmp0$ID.current, tmp0$ID.futur, function(x)nlevels(factor(x)))
id.patch <- names(which(z!=1))
tmp0 <- tmp0[which(tmp0$ID.futur %in% id.patch),]
z <- tapply(tmp0$ID.current, tmp0$ID.futur, function(x)nlevels(factor(x)))

tmp1 <- subset(tmp0, ID.futur==20) 
coordinates(tmp1) <- c("x", "y")
proj4string(tmp1) <- crs(Rast)
tmp1 <- as.data.frame(spTransform(tmp1, crs(PA.cur))) # WGS 84

tmp2 <- tmp0[which(tmp0$ID.current %in% levels(factor(tmp1$ID.current))),]
coordinates(tmp2) <- c("x", "y")
proj4string(tmp2) <- crs(Rast)
tmp2 <- as.data.frame(spTransform(tmp2, crs(PA.cur))) # WGS 84

### Write down rasters
rast1 <- rasterize(SpatialPoints(tmp1[c("x", "y")]), PA.cur, field=tmp1$PA.fut) # futur
rast2 <- rasterize(SpatialPoints(tmp2[c("x", "y")]), PA.cur, field=tmp1$PA.cur) # current
stack.rast <- stack(rast2, rast1)
names(stack.rast) <- c("Current", "Futur")
writeRaster(stack.rast, filename = "outputs/raster_files/Raster_expansion_abs_patch.tif", overwrite=T)

### Select only pixels transitioning for plotting purposes
tmp2 <- tmp2[which(tmp2$PA.cur == 1),]

plot.exp.abs <- ggplot(dif.HS, aes(x=x, y=y))+
  geom_raster(alpha=0.2, fill="grey")+
  geom_point(data=tmp1, aes(x=x, y=y, fill=NULL, color="#EDF285"), size=0.1)+
  geom_point(data=tmp2, aes(x=x, y=y, fill=NULL, color="#E05C3A"), size=0.1)+
  theme_void()+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),
        strip.text.x = element_text(size=8),
        strip.background = element_blank())+
  coord_equal()+
  scale_color_manual(values = c("#E05C3A", "#EDF285"), labels=c("P->A", "A->A"))+
  labs(color="Patch transition")+
  guides(color = guide_legend(override.aes = list(size=8)))

svg("Outputs/Example_expansion_abs_patch.svg", width=9, height=9)
print(plot.exp.abs)
dev.off()

### Compute landscape metrics for this patch
tmp1 <- subset(tmp0, ID.futur==20) 
tmp2 <- tmp0[which(tmp0$ID.current %in% levels(factor(tmp1$ID.current))),]

rast1 <- rasterize(SpatialPoints(tmp1[c("x", "y")]), Rast.PA, field=tmp1$PA.fut)
rast2 <- rasterize(SpatialPoints(tmp2[c("x", "y")]), Rast.PA, field=tmp2$PA.cur)
stack.rast <- stack(rast2, rast1)

func.name <- as.data.frame(list_lsm()[,5])[,1]
metric.name <- as.data.frame(list_lsm()[,1])[,1]
out.metric.other <- NULL
for(i in 1:length(metric.name)){
  func <- get(func.name[i])
  metric.foc <- func(stack.rast)
  out.metric.other <- rbind(out.metric.other, metric.foc)
}
write.table(out.metric.other, "Outputs/Landscape_metrics_expansion.abs_example.txt")

##############################################
# Example of fragmentation of presence patches
##############################################

### Process data
tmp0 <- tmp[which(tmp$PA.cur == 1),]
z <- tapply(tmp0$ID.futur, tmp0$ID.current, function(x)nlevels(factor(x)))
id.patch <- names(which(z!=1))
tmp0 <- tmp0[which(tmp0$ID.current %in% id.patch),]
z <- tapply(tmp0$ID.futur, tmp0$ID.current, function(x)nlevels(factor(x)))

tmp1 <- subset(tmp0, ID.current==74) 
coordinates(tmp1) <- c("x", "y")
proj4string(tmp1) <- crs(Rast)
tmp1 <- as.data.frame(spTransform(tmp1, crs(PA.cur))) # WGS 84

tmp2 <- tmp0[which(tmp0$ID.futur %in% levels(factor(tmp1$ID.futur))),]
coordinates(tmp2) <- c("x", "y")
proj4string(tmp2) <- crs(Rast)
tmp2 <- as.data.frame(spTransform(tmp2, crs(PA.cur))) # WGS 84
tmp2$index <- paste(tmp2$x, tmp2$y, sep="_")

### Write down rasters
rast1 <- rasterize(SpatialPoints(tmp1[c("x", "y")]), PA.cur, field=tmp1$PA.cur) # futur
rast2 <- rasterize(SpatialPoints(tmp2[c("x", "y")]), PA.cur, field=tmp1$PA.fut) # current
stack.rast <- stack(rast1, rast2)
names(stack.rast) <- c("Current", "Futur")
writeRaster(stack.rast, filename = "outputs/raster_files/Raster_fragmentation_presence_patch.tif", overwrite=T)

### Select only pixels transitioning for plotting purposes
tmp2 <- tmp2[which(tmp2$PA.fut == 0),]

plot.frag <- ggplot(dif.HS, aes(x=x, y=y))+
  geom_raster(alpha=0.2, fill="grey")+
  geom_point(data=tmp1, aes(x=x, y=y, fill=NULL, color="#24E0A2"), size=0.1)+
  geom_point(data=tmp2, aes(x=x, y=y, fill=NULL, color="#E05C3A"), size=0.1)+
  theme_void()+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),
        strip.text.x = element_text(size=8),
        strip.background = element_blank())+
  coord_equal()+
  scale_color_manual(values = c("#24E0A2", "#E05C3A"), labels=c("P->P", "P->A"))+
  labs(color="Patch status")+
  guides(color = guide_legend(override.aes = list(size=8)))

svg("Outputs/Example_fragmentation_presence_patch.svg", width=9, height=9)
print(plot.frag)
dev.off()

### Compute landscape metrics for this patch
tmp1 <- subset(tmp0, ID.current==74) 
tmp2 <- tmp0[which(tmp0$ID.futur %in% levels(factor(tmp1$ID.futur))),]

rast1 <- rasterize(SpatialPoints(tmp1[c("x", "y")]), Rast.PA, field=tmp1$PA.cur)
rast2 <- rasterize(SpatialPoints(tmp2[c("x", "y")]), Rast.PA, field=tmp2$PA.fut)
stack.rast <- stack(rast1, rast2)

func.name <- as.data.frame(list_lsm()[,5])[,1]
metric.name <- as.data.frame(list_lsm()[,1])[,1]
out.metric.other <- NULL
for(i in 1:length(metric.name)){
  func <- get(func.name[i])
  metric.foc <- func(stack.rast)
  out.metric.other <- rbind(out.metric.other, metric.foc)
}
write.table(out.metric.other, "Outputs/Landscape_metrics_extinction_example.txt")

#----------- Histograms ENN and patch area

Current <- projectRaster(PA.cur, method = "ngb", res=5000, crs = New.proj) 
Future <- projectRaster(PA.fut, method = "ngb", res=5000, crs = New.proj) 

###########
# patch ENN
###########

tmp1 <- show_lsm(Current, what = "lsm_p_enn", directions = 8, label_lsm = T)
tmp2 <- show_lsm(Future, what = "lsm_p_enn", directions = 8, label_lsm = T)
A <- tmp1$data
B <- tmp2$data

### Raster maps
A$Period <- "Current"
B$Period <- "Futur"
A <- A[which(!is.na(A$id)),]
B <- B[which(!is.na(B$id)),]
All <- rbind(A, B)

### Save rasters
rast1 <- rasterize(SpatialPoints(A[c("x", "y")]), Current, field=log(A$value)) # futur
rast2 <- rasterize(SpatialPoints(B[c("x", "y")]), Current, field=log(B$value)) # current
stack.rast <- stack(rast1, rast2)
names(stack.rast) <- c("Current", "Futur")
writeRaster(stack.rast, filename = "outputs/raster_files/Raster_isolated_patch.tif", overwrite=T)

### Plot the spatial pattern
Fig.isolated.patch <- ggplot(All, aes(x=x, y=y, fill=log(value)))+
  geom_raster()+
  scale_fill_viridis("ENN")+
  theme_void()+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),
        strip.text.x = element_text(size=8),
        strip.background = element_blank())+
  coord_equal()+
  facet_wrap(~Period+class)
svg("Outputs/isolated.patch.svg", width=9, height=9)
print(Fig.isolated.patch)
dev.off()

### Density plots of ENN
z1 <- cbind(Val=A$value, class=A$class, Period="Current")
z2 <- cbind(Val=B$value, class=A$class, Period="Future")
z0 <- data.frame(na.omit(rbind(z1, z2)))
z0$Val <- as.numeric(as.character(z0$Val))
z0$Val <- log(z0$Val)

hist_ENN <- ggplot(z0, aes(x=Val, fill=Period))+ 
  scale_fill_manual(values=c("#ff7f00", "#6a3d9a"))+
  scale_color_manual(values=c("#ff7f00", "#6a3d9a"))+
  geom_density(data=subset(z0, Period == 'Current'), alpha = 0.4, bw=0.2) +
  geom_density(data=subset(z0, Period == 'Future'),  alpha = 0.4, bw=0.2) +
  geom_rug(aes(color=Period))+
  theme_classic()+
  facet_wrap(~class)+
  labs(y="Density", x="Euclidean Nearest Neighbour (log-scale)")
svg("Outputs/Histogram_ENN.svg", width=9, height=9)
print(hist_ENN)
dev.off()

############
# patch Area
############

tmp1 <- show_lsm(Current, what = "lsm_p_area", directions = 8, label_lsm = T)
tmp2 <- show_lsm(Future, what = "lsm_p_area", directions = 8, label_lsm = T)
A <- tmp1$data
B <- tmp2$data
z1 <- cbind(Val=A$value, class=A$class, Period="Current")
z2 <- cbind(Val=B$value, class=A$class, Period="Future")
z0 <- data.frame(na.omit(rbind(z1, z2)))
z0$Val <- as.numeric(as.character(z0$Val))
z0$Val <- log(z0$Val)

hist_area <- ggplot(z0, aes(x=Val, fill=Period))+ 
  scale_fill_manual(values=c("#ff7f00", "#6a3d9a"))+
  scale_color_manual(values=c("#ff7f00", "#6a3d9a"))+
  geom_density(data=subset(z0, Period == 'Current'), alpha = 0.4, bw=0.25, adjust=1) +
  geom_density(data=subset(z0, Period == 'Future'),  alpha = 0.4, bw=0.25, adjust=1) +
  geom_rug(aes(color=Period))+
  theme_classic()+
  facet_wrap(~class)+
  labs(y="Density", x="Patch area")
svg("Outputs/Histogram_AREA.svg", width=9, height=9)
print(hist_area)
dev.off()

#----------- Patch transition matrix (A. Boye)

DF.all <- readRDS("Outputs/Table.metrics.patch.Rdata")
DF.all$HS.dif <- DF.all$HS.fut - DF.all$HS.cur # compute difference in suitability
DF.all <- unite(DF.all, "transition", c("PA.cur", "PA.fut"), sep = "->", remove = FALSE)

fact_order_cur <- DF.all %>%
  dplyr::group_by(ID.current) %>%
  summarise(y=max(y), PA.cur=mean(PA.cur)) %>%
  ungroup() %>%
  arrange(y) %>%
  mutate(color = if_else(PA.cur == 0, muted("red"), muted("green")))

fact_order_fut <- DF.all %>%
  dplyr::group_by(ID.futur) %>%
  summarise(y=max(y), PA.fut=mean(PA.fut)) %>%
  ungroup() %>%
  arrange(y) %>%
  mutate(color = if_else(PA.fut == 0, muted("red"), muted("green")))

transition_mat <- DF.all %>% 
  group_by(ID.current,ID.futur) %>%
  summarise(mean_hs_diff=mean(HS.dif)) %>%
  ungroup()

transition_mat <- table(DF.all$ID.current, DF.all$ID.futur) %>%
  as.data.frame() %>%
  rename(ID.current = Var1, ID.futur = Var2) %>%
  mutate_at(vars(c(ID.current,ID.futur)), as.character) %>%
  mutate_at(vars(c(ID.current,ID.futur)), as.numeric) %>%
  left_join(., transition_mat) %>%
  mutate(ID.current = factor(ID.current, levels=fact_order_cur$ID.current, ordered = TRUE)) %>%
  mutate(ID.futur = factor(ID.futur, levels=fact_order_fut$ID.futur, ordered = TRUE)) %>%
  drop_na()

transition_mat_summarised_x <- DF.all %>% 
  group_by(ID.current) %>%
  summarise(freq = n(),mean_hs_diff=mean(HS.dif)) %>%
  ungroup() %>%
  mutate(ID.current = factor(ID.current, levels=fact_order_cur$ID.current, ordered = TRUE))

transition_mat_summarised_y <- DF.all %>% 
  group_by(ID.futur) %>%
  summarise(freq=n(), mean_hs_diff=mean(HS.dif)) %>%
  ungroup() %>%
  mutate(ID.futur = factor(ID.futur, levels=fact_order_fut$ID.futur, ordered = TRUE))

transition.mat <- ggplot(transition_mat, aes(x=ID.current, y=ID.futur, fill = mean_hs_diff)) +
  geom_tile(col = "black", size = 0.8, alpha = 0.8) +
  scale_fill_gradient2(high = muted("green"),midpoint=0, name = "Mean difference\nin habitat suitability") +
  scale_x_discrete(labels = transition_mat_summarised_x$ID.current) +
  scale_y_discrete(labels = transition_mat_summarised_y$ID.futur) +
  theme_bw() +
  theme(axis.text.x = element_text(colour = fact_order_cur$color),
        axis.text.y = element_text(colour = fact_order_fut$color),
        panel.grid.major.x = element_line(colour = fact_order_cur$color, size = 0.1),
        panel.grid.major.y = element_line(colour = fact_order_fut$color, size = 0.1),
        legend.key = element_rect(fill = NA), ggside.panel.scale = .025) +
  xlab("Current patches\n(ordered by poleward edge latitude; green label = 1, red = 0)") +
  ylab("Future patches\n(ordered by poleward edge latitude; green label = 1, red = 0)") +
  geom_ysidepoint(data = transition_mat_summarised_y, aes(x="", y = ID.futur,size = log(freq), fill = mean_hs_diff), shape = 21, alpha = 0.8) +
  geom_xsidepoint(data = transition_mat_summarised_x, aes(y="", x = ID.current,size = log(freq), fill = mean_hs_diff), shape = 21, alpha = 0.8) +
  # Remove ticks and labels from the panels
  scale_ysidex_discrete(breaks = NULL, labels = NULL) +
  scale_xsidey_discrete(breaks = NULL, labels = NULL) +
  scale_size_continuous(range = c(0.1,10)) +
  # Clip or not the points within the panels
  coord_cartesian(clip = 'off') +
  guides(size = guide_legend(title = "Number of pixels (log)", override.aes=list(fill = "grey")))
svg("Outputs/Patch_transition_matrix_aurel.svg", width=20, height=10)
print(transition.mat)
dev.off()

#----------- Patch transition vectors (M. Marzloff)

DF.all$transition <- as.factor(DF.all$transition)
DF.all$transition_num <- as.numeric(DF.all$transition)
DF.all$ID.current <- factor(DF.all$ID.current) # levels = fact_order_cur$ID.current, ordered = T)
DF.all$ID.futur <- factor(DF.all$ID.futur) # levels = fact_order_fut$ID.futur, ordered = T)
AlltransitionPerPatch <- as.data.frame(xtabs(~ ID.current+ID.futur+transition, data = DF.all))
AlltransitionPerPatch[(AlltransitionPerPatch$Freq!=0),]

### Present patch point of view - Changes
TransitionPresentPatch <- (xtabs(~ ID.current+transition, data = DF.all))
colnames(TransitionPresentPatch) <- levels(DF.all$transition)
rownames(TransitionPresentPatch) <- levels(DF.all$ID.current)
TransitionPresentPatch <- cbind(TransitionPresentPatch, round(TransitionPresentPatch/rowSums(TransitionPresentPatch), digits = 2), rowSums(TransitionPresentPatch))
TransitionPresentPatch[TransitionPresentPatch == 0] <- NA
data_plot <- expand.grid(X=rownames(TransitionPresentPatch), Y=levels(DF.all$transition)) # Turn to data.frame for GGPLOT2
data_plot $Z <- c(TransitionPresentPatch[,5:8])
data_plot $CurArea <- rep(TransitionPresentPatch[,9], rep= 4)
data_plot$X <- factor(data_plot$X, levels=fact_order_cur$ID.current, ordered = TRUE) # Define northwards order of patches

### Future patch point of view - Changes
TransitionFuturePatch <- xtabs(~ ID.futur+transition, data = DF.all)
colnames(TransitionFuturePatch) <- levels(DF.all$transition)
rownames(TransitionFuturePatch) <- levels(DF.all$ID.futur)
TransitionFuturePatch <- cbind(TransitionFuturePatch, round(TransitionFuturePatch/rowSums(TransitionFuturePatch), digits = 2), rowSums(TransitionFuturePatch))
TransitionFuturePatch[TransitionFuturePatch == 0] <- N
data_plot = expand.grid(X=rownames(TransitionFuturePatch), Y=levels(DF.all$transition) )
data_plot $Z = c( TransitionFuturePatch[,5:8])
data_plot $CurArea = rep(TransitionFuturePatch[,9], rep= 4)
data_plot$X <- factor(data_plot$X, levels=fact_order_fut$ID.futur, ordered = TRUE)
data_plot$Y <- factor(data_plot$Y, levels = (levels(data_plot$Y)[c(1,3,2,4)]))

### The plots
patch_pres <- ggplot(data=data_plot, aes(x =X, y = Y, size= Z , fill = log(CurArea))) + 
  geom_point(shape = 21)+
  scale_size_binned(name = "Relative Contribution", breaks = c(0,0.25,0.75,1))+
  scale_fill_gradientn(name = 'Patch size', colours = colorRamps::matlab.like2(100), limits=c(0, 8), labels = c('S','M','L','XL','XXL')) +
  scale_x_discrete(breaks = c(min(data_plot$X),max(data_plot$X)), labels = c('South','North')) +
  scale_y_discrete(labels = levels(data_plot$Y)) +
  xlab("")+
  ylab("")+
  coord_flip()+
  ggtitle('Changes occuring to current patches', subtitle = 'Column-focused interpretation of Transition matrix')

patch_fut <- ggplot(data=data_plot, aes(x =X, y = Y, size= Z , fill = log(CurArea))) + 
  geom_point( shape = 21)+
  scale_size_binned(name = "Relative Contribution", breaks = c(0,0.25,0.75,1))+
  scale_fill_gradientn(name = 'Patch size', colours = colorRamps::matlab.like2(100), limits=c(0, 8), labels = c('S','M','L','XL','XXL')) +
  scale_x_discrete(breaks = c(min(data_plot$X),max(data_plot$X)), labels = c('South','North')) +
  scale_y_discrete(labels = levels(data_plot$Y)) +
  xlab("")+
  ylab("")+
  coord_flip()+
  ggtitle('Changes leading to future patches', subtitle = 'Row-focused interpretation of transition matrix')

svg("Outputs/Patch_transition_vector_martin.svg", width=18, height=10)
ggarrange(patch_pres, patch_fut, ncol=2)
dev.off()
