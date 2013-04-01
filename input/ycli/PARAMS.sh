#! /bin/bash

# -----------------------------------------------------------
# MAP_SIM : the simulation maps WITHOUT modes subtracted.
#           These maps are used in reference calculation, 
#           for auto power, we use MAP_SIM + "_temperature";
#           for cross power, we use MAP_SIM + "temperature" 
#           instead of GBT, and MAP_SIM + "_delta" instead
#           of WiggleZ. These suffixes are added automatically
#           in pipeline.
# MAP_SSM : the simulation maps WITH modes subtracted 
#           these maps are used for calculate the transfer
#           functions
# MAP_CLN : the cleaned maps, used for power spectrum. 
#           and their noise_inv maps are used as the weight
#           in reference calculation and transfer function
#           calculation.
#           the pipeline will automatically choose the "_combined"
#           maps, except for the auto power spectrum.
# MAP_WIG : WiggleZ maps. the pipeline will automatically
#           choose their "_delta_binned_data", "_delta_mock",
#           and "_montecarlo"


#export MAP_SIM="sim_15hr_oldmap_ideal"
export MAP_SIM="None"
export MAP_SSM="None"
export MAP_CLN="None"
export MAP_WGZ="None"
export MAP_MULTIPLIER=1
export MODE="15"
export HOUR="15"
# -----------------------------------------------------------

export GOODMODES="5"
export SIMROOT="/cita/d/data-pen3/ycli/map_simulation/"
export SIMFILE="None"
export MAPFILE="sim_15hr_oldmap_str_beam_plus_fg"
export SUBREAL=""
export MAPGBT_MY="/mnt/scratch-3week/ycli/"
export NOISEWEIGHT="sim_15hr_oldmap_str_beam_plus_fg"

# -----------------------------------------------------------
# For Parkes
export BEAMS="None"
