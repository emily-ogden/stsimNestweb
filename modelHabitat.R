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


## Connect to SyncroSim ----
myScenario <- scenario()

# Load relevant datasheets
# stsim
RunControl <- datasheet(myScenario, "stsim_RunControl")
Stratum <- datasheet(myScenario, "stsim_Stratum")
SecondaryStrtaum <- datasheet(myScenario, "stsim_SecondaryStratum")
StateClass <- datasheet(myScenario, "stsim_StateClass")
InitialConditionsSpatial <- datasheet(myScenario, "stsim_InitialConditionsSpatial")
OutputSpatialState <- datasheet(myScenario, "stsim_OutputSpatialState")

#stsimsf
OutputSpatialStockGroup <- datasheet(myScenario, "stsimsf_OutputSpatialStockGroup")
StockType <- datasheet(myScenario, "stsimsf_StockType")
StockGroup <- datasheet(myScenario, "stsimsf_StockGroup")
StockGroupMembership <- datasheet(myScenario, "stsimsf_StockGroupMembership")

# stsimNetweb
OutputOptions <- datasheet(myScenario, "stsimNestweb_OutputOptions")
OutputOptionsSpatial <- datasheet(myScenario, "stsimNestweb_OutputOptionsSpatial")
HabitatModel <- datasheet(myScenario, "stsimNestweb_HabitatModel")

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

# Square meter to hectare conversion
scaleFactor <- 0.0001

# Load template raster
templateRaster <- rast(InitialConditionsSpatial$StratumFileName)
templateRaster[!is.na(templateRaster)] <- 1

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

## Function definitions ----

# Main Code Here zzz ----
progressBar(type = "message", message = "Running main code...")

# Build parameter sampling table
models <- map_chr(HabitatModel$ModelFileName, load) %>% 
  map(get) %>% 
  set_names(HabitatModel$Name)

# enframe(models, name = "Species", value = "Model") %>% 
#   mutate()

# imap_dfr(
#   models, 
#   ~{
#     tibble(
#       Species = .y,
#       Parameter = c("Aspen Cover", "Diameter"),
#       Mean = c(0),
#       SD = c(0)
#     )
#     })

# piwo_mods$coefArray %>% as_tibble

for(iteration in iterations){
  # Sample parameters
  for(timestep in timesteps){
    # Get relevant spatial data
    for(aSpecies in species){
      # Predict and accumulate outputs
      
    }
  }
}

# Clean up ----

# Wrap up SyncroSim progress bar
progressBar("end")