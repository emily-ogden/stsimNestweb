# zzz:
# Load spatial data to disambiguating habitat suitability (strata, site)
# Create a datasheet for Site
# Fix project-scoep site datasheet

# Set environment variable TZ when running on AWS EC2 instance
Sys.setenv(TZ='UTC')

# Settings
options(stringsAsFactors=FALSE, SHAPE_RESTORE_SHX=T, useFancyQuotes = F, digits=10)

# Load libraries
library(rsyncrosim)
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(terra))
suppressPackageStartupMessages(library(MuMIn))
suppressPackageStartupMessages(library(VGAM))
suppressPackageStartupMessages(library(unmarked))

# Setup ----
progressBar(type = "message", message = "Preparing inputs...")

## Function definitions ----
# Define function to facilitate recoding a vector using a look-up table
lookup <- function(x, old, new){
  dplyr::recode(x, !!!set_names(new, old))
}

## Connect to SyncroSim ----
myScenario <- scenario()

# Load relevant datasheets
# stsim
RunControl <- datasheet(myScenario, "stsim_RunControl")
Stratum <- datasheet(myScenario, "stsim_Stratum")
SecondaryStratum <- datasheet(myScenario, "stsim_SecondaryStratum")
StateClass <- datasheet(myScenario, "stsim_StateClass")
InitialConditionsSpatial <- datasheet(myScenario, "stsim_InitialConditionsSpatial")
OutputSpatialState <- datasheet(myScenario, "stsim_OutputSpatialState")

#stsimsf
OutputSpatialStockGroup <- datasheet(myScenario, "stsimsf_OutputSpatialStockGroup")
StockType <- datasheet(myScenario, "stsimsf_StockType")
StockGroup <- datasheet(myScenario, "stsimsf_StockGroup")
StockTypeGroupMembership <- datasheet(myScenario, "stsimsf_StockTypeGroupMembership")

# stsimNetweb
SiteType <- datasheet(myScenario, "stsimNestweb_SiteType")
SpeciesID <- datasheet(myScenario, "stsimNestweb_Species", includeKey = TRUE) %>% 
  pull(SpeciesID, name = Name)
Site <- datasheet(myScenario, "stsimNestweb_SiteValue")
OutputOptions <- datasheet(myScenario, "stsimNestweb_OutputOptions")
OutputOptionsSpatial <- datasheet(myScenario, "stsimNestweb_OutputOptionsSpatial")
HabitatModel <- datasheet(myScenario, "stsimNestweb_HabitatModel")
OutputHabitatAmount <- data.frame()


## Setup Parameters ----
# Iterations 
iterations <- seq(RunControl$MinimumIteration, RunControl$MaximumIteration)

# Timesteps
timestepsTabular <- if(OutputOptions$SummaryOutputHA){
  seq(RunControl$MinimumTimestep, RunControl$MaximumTimestep, by = OutputOptions$SummaryOutputHATimesteps) %>% 
    c(RunControl$MaximumTimestep) %>% 
    unique()
} else c()

timestepsSpatial <- if(OutputOptionsSpatial$RasterOutputHA){
  seq(RunControl$MinimumTimestep, RunControl$MaximumTimestep, by = OutputOptionsSpatial$RasterOutputHATimesteps) %>% 
    c(RunControl$MaximumTimestep) %>% 
    unique()
} else c()

timesteps <- c(timestepsTabular, timestepsSpatial) %>% 
  unique() %>% 
  sort()

# Species
species <- HabitatModel$Name

# Square meter to hectare conversion
scaleFactor <- 0.0001

## Load spatial data ----
StrataData <- data.frame(
  StratumID = rast(InitialConditionsSpatial$StratumFileName)[] %>% as.vector(),
  SecondaryStratumID = rast(InitialConditionsSpatial$SecondaryStratumFileName)[] %>% as.vector(),
  Site = rast(Site$FileName)[] %>% as.vector() %>% lookup(SiteType$ID, SiteType$Name))

# Load template raster
templateRaster <- rast(InitialConditionsSpatial$StratumFileName)
templateRaster[!is.na(templateRaster)] <- 1
names(templateRaster) <- "habitatSuitability"

# Pixel resolution
cellResolution <- templateRaster %>% res()

# Pixel count
cellCount <- freq(templateRaster, value = 1)$count

# Area of a pixel (square meter)
cellArea <- cellResolution[1]^2

# Total area
totalArea <- cellCount * cellArea



## Setup files and folders ----

# Create temp folder, ensure it is empty
tempDir <- ssimEnvironment()$TempDirectory

# Generate filenames for potential outputs

## Handle empty values ----

# Main Code Here ----
progressBar(type = "message", message = "Running main code...")
progressBar(type = "begin", totalSteps = length(iterations) * length(timesteps) * length(species))

# Build parameter sampling table
modelNames <- map_chr(HabitatModel$ModelFileName, load)
for(m in HabitatModel$ModelFileName) load(m)
models <- modelNames %>% 
  map(get) %>% 
  set_names(HabitatModel$Name)
rm(list = modelNames)

# Parameterize the sampling distribution for each parameter and model
parameterTable <- imap_dfr(
  models,
  ~{
    # Connect clean parameter names to variable names in model fit
    parameters <- c('Aspen Cover' = 'scale(Perc_At)',
                    'Diameter' = 'scale(Median_DBH)')
    means <- coef(.x)[parameters]
    stdErrors <- .x$coefArray[, 'Std. Error', parameters] %>% 
      apply(MARGIN = 2, FUN = mean, na.rm = TRUE) # MARGIN = 2 is used to average within columns 
    
    tibble(
      species = .y,
      parameter = names(parameters),
      mean = means,
      sd = stdErrors)}) # Double check this assumption



for(iteration in iterations){
  # Sample parameters
  parameterTable <- parameterTable %>% 
    mutate(value = map2_dbl(mean, sd, rnorm, n = 1))
  
  for(timestep in timesteps){
    # Load aspen cover and diameter
    aspenCover <- datasheetRaster(
      ssimObject = myScenario, 
      datasheet = "stsimsf_OutputSpatialStockGroup", 
      iteration = iteration, 
      timestep = timestep,
      filterColumn = "StockGroupID",
      filterValue = "Aspen Cover (%) [Type]")
    
    # Convert aspen raster to proportion (rather than percentage)
    aspenCover[] <- aspenCover[]/100
    
    diameter <- datasheetRaster(
      ssimObject = myScenario, 
      datasheet = "stsimsf_OutputSpatialStockGroup", 
      iteration = iteration, 
      timestep = timestep,
      filterColumn = "StockGroupID",
      filterValue = "Diameter (cm) [Type]")
    
    # Create dataframe of habitat suitability model inputs
    habitatSuitabilityDf <- data.frame(Perc_At = aspenCover[], 
                                       Median_DBH = diameter[],
                                       edge_near = 0,
                                       Num_2BI = 0,
                                       Mean_decay = 0,
                                       dist_to_cut = 0,
                                       cut_harvest0 = "N",
                                       Site = StrataData$Site)
    
    for(aSpecies in species){
      # Predict habitat suitability
      model <- models[[aSpecies]]
      habitatSuitabilityDf$pred <- predict(model, newdata = habitatSuitabilityDf, type = "response", allow.new.levels = TRUE)
      
      # Output raster
      if(timestep %in% timestepsSpatial) {
        outputFilename <- file.path(tempDir, str_c("hs.sp", SpeciesID[aSpecies], ".it", iteration, ".ts", timestep, ".tif")) %>% 
          normalizePath(mustWork = FALSE)
        
        rast(templateRaster, vals = habitatSuitabilityDf$pred) %>% 
          writeRaster(outputFilename, overwrite = TRUE)}
      
      # Calculate tabular output and append to 
      if(timestep %in% timestepsTabular) {
        OutputHabitatAmount <- bind_rows(
          OutputHabitatAmount, 
          habitatSuitabilityDf %>% 
            dplyr::select(Amount = pred) %>% 
            bind_cols(StrataData) %>% 
            filter(!is.na(Amount)) %>% 
            group_by(StratumID, SecondaryStratumID, Site) %>%
            summarise(Amount = mean(Amount), .groups = "drop") %>%
            mutate(
              Timestep = timestep,
              Iteration = iteration,
              Species = aSpecies
            ))}
      
      # Increment progress bar
      progressBar()
    }
  }
}

# Save spatial outputs
OutputSpatialHabitat <- tibble(FileName = list.files(tempDir, ".tif", full.names = TRUE) %>% normalizePath()) %>%
  mutate(
    temp = basename(FileName),
    Iteration = temp %>% str_extract("it\\d+") %>% str_replace("it", "") %>% as.numeric(),
    Timestep = temp %>% str_extract("ts\\d+") %>% str_replace("ts", "") %>% as.numeric(),
    Species = temp %>% str_extract("sp\\d+") %>% str_replace("sp", "") %>% as.numeric(),
    Species = lookup(Species, SpeciesID, names(SpeciesID))) %>% 
  dplyr::select(-temp) %>% 
  as.data.frame()

saveDatasheet(myScenario, OutputSpatialHabitat, "stsimNestweb_OutputSpatialHabitat")

OutputHabitatAmount <- OutputHabitatAmount %>% 
  mutate(
    StratumID = StratumID %>% lookup(Stratum$ID, Stratum$Name),
    SecondaryStratumID = SecondaryStratumID %>% lookup(SecondaryStratum$ID, SecondaryStratum$Name))
saveDatasheet(myScenario, OutputHabitatAmount, "stsimNestweb_OutputHabitatAmount")


# Clean up ----

# Wrap up SyncroSim progress bar
progressBar("end")