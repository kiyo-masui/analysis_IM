#! /bin/bash

export PIPELINE_DIR="/home/ycli/pipeline"
export PIPELINE_DIR_P=$PIPELINE_DIR"/powerspectrum"
cd $PIPELINE_DIR
source $PIPELINE_DIR/PARAMS.sh

#===== Cross Power ======

# cross power for 15hr
MAP_SIM="/mnt/raid-project/gmrt/ycli/map_sim/15hr/"
MAP_SSM="/mnt/raid-project/gmrt/ycli/foreground_cleand/15hr_IE_legendre_modes_0gwj/"
MAP_CLN="/mnt/raid-project/gmrt/ycli/foreground_cleand/15hr_IE_legendre_modes_0gwj/"
MAP_WGZ="WiggleZ_15hr"
MODE="20"
HOUR="15"
python $PIPELINE_DIR/manager.py $PIPELINE_DIR_P/reference_cros.pipe
python $PIPELINE_DIR/manager.py $PIPELINE_DIR_P/reference_cros_beam.pipe
python $PIPELINE_DIR/manager.py $PIPELINE_DIR_P/transfer_cros.pipe
python $PIPELINE_DIR/manager.py $PIPELINE_DIR_P/power_cros.pipe
#python $PIPELINE_DIR/manager.py $PIPELINE_DIR/powespectrum_P/sim_power_cros.pipe

MAP_SIM="/mnt/raid-project/gmrt/ycli/map_sim/15hr/"
MAP_SSM="/mnt/raid-project/gmrt/ycli/foreground_cleand/15hr_EI_legendre_modes_0gwj/"
MAP_CLN="/mnt/raid-project/gmrt/ycli/foreground_cleand/15hr_EI_legendre_modes_0gwj/"
MAP_WGZ="WiggleZ_15hr"
MODE="20"
HOUR="15"
python $PIPELINE_DIR/manager.py $PIPELINE_DIR_P/reference_cros.pipe
python $PIPELINE_DIR/manager.py $PIPELINE_DIR_P/reference_cros_beam.pipe
python $PIPELINE_DIR/manager.py $PIPELINE_DIR_P/transfer_cros.pipe
python $PIPELINE_DIR/manager.py $PIPELINE_DIR_P/power_cros.pipe
#python $PIPELINE_DIR/manager.py $PIPELINE_DIR/powespectrum_P/sim_power_cros.pipe

echo 'Done'
