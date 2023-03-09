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

# Setup ----
progressBar(type = "message", message = "Preparing inputs...")

## Connect to SyncroSim ----
myScenario <- scenario()

# Load relevant datasheets
OutputOptionsSpatial <- datasheet(myScenario, "stsimNestweb_OutputOptionsSpatial")
OutputSpatialHabitat <- datasheet(myScenario, "stsimNestweb_OutputSpatialHabitat")
SpeciesID <- datasheet(myScenario, "stsimNestweb_Species", includeKey = TRUE) %>% 
  pull(SpeciesID, name = Name)

## Setup Parameters ----
# Timesteps 
timesteps <- OutputSpatialHabitat$Timestep %>% 
  unique() %>% 
  sort()

# Species
species <- OutputSpatialHabitat$Species %>% 
  unique() %>% 
  sort()

## Setup files and folders ----

# Create temp folder, ensure it is empty
tempDir <- ssimEnvironment()$TempDirectory

spatialOutputDir <- file.path(tempDir, "SpatialOutputs") %>% normalizePath(mustWork = FALSE)

unlink(spatialOutputDir, recursive = TRUE, force = TRUE)

dir.create(spatialOutputDir)

# Main Code Here ----
if(OutputOptionsSpatial$RasterOutputHAAverage) {
  progressBar(type = "message", message = "Running main code...")
  progressBar(type = "begin", totalSteps = length(iterations) * length(timesteps) * length(species))
  
  for(timestep in timesteps){
    for(aSpecies in species){
      outputFilename <- file.path(spatialOutputDir, str_c("hsa.sp", SpeciesID[aSpecies], ".ts", timestep, ".tif")) %>% 
        normalizePath(mustWork = FALSE)
      
      habitatSuitability <- datasheetRaster(
        ssimObject = myScenario, 
        datasheet = "stsimNestweb_OutputSpatialHabitat", 
        timestep = timestep,
        #iteration = 1,
        filterColumn = "Species",
        filterValue = aSpecies)
    }
  }

}
