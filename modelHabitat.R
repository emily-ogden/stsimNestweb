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
RunControl <- datasheet(myScenario, "stsim_RunControl") # 



## Setup Parameters ----
# zzz: read cell size and total area from raster
# Cell size (90m x 90m) to Ha conversion
scaleFactor <- 0.81

# Total area (Ha) of analysis area
totalArea <- 336475.6200

## Setup files and folders ----

# Create temp folder, ensure it is empty
tempDir <- ssimEnvironment()$TempDirectory

# Generate filenames for potential outputs

## Handle empty values ----

## Function definitions ----

# Main Code Here zzz ----
progressBar(type = "message", message = "Running main code...")

# Clean up ----

# Wrap up SyncroSim progress bar
progressBar("end")