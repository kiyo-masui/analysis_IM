#! /bin/bash

export PIPELINE_DIR="/home/ycli/pipeline"
export PIPELINE_DIR_P=$PIPELINE_DIR"/powerspectrum"
cd $PIPELINE_DIR
source $PIPELINE_DIR/PARAMS.sh

#===== Cross Power ======

# cross power for 1hr commen beam
MAP_SIM="/mnt/raid-project/gmrt/ycli/map_sim/1hr/"
MAP_SSM="/mnt/raid-project/gmrt/ycli/foreground_cleand/1hr_II_legendre_modes_0gwj_conv/"
MAP_CLN="/mnt/raid-project/gmrt/ycli/foreground_cleand/1hr_II_legendre_modes_0gwj_conv/"
MAP_WGZ="WiggleZ_1hr"
MODE="20"
HOUR="1"
python $PIPELINE_DIR/manager.py $PIPELINE_DIR_P/reference_cros.pipe
python $PIPELINE_DIR/manager.py $PIPELINE_DIR_P/reference_cros_degradebeam.pipe
python $PIPELINE_DIR/manager.py $PIPELINE_DIR_P/transfer_cros.pipe
python $PIPELINE_DIR/manager.py $PIPELINE_DIR_P/power_cros.pipe
#python $PIPELINE_DIR/manager.py $PIPELINE_DIR/powespectrum_P/sim_power_cros.pipe

MAP_SIM="/mnt/raid-project/gmrt/ycli/map_sim/1hr/"
MAP_SSM="/mnt/raid-project/gmrt/ycli/foreground_cleand/1hr_IQ_legendre_modes_0gwj_conv/"
MAP_CLN="/mnt/raid-project/gmrt/ycli/foreground_cleand/1hr_IQ_legendre_modes_0gwj_conv/"
MAP_WGZ="WiggleZ_1hr"
MODE="20"
HOUR="1"
python $PIPELINE_DIR/manager.py $PIPELINE_DIR_P/reference_cros.pipe
python $PIPELINE_DIR/manager.py $PIPELINE_DIR_P/reference_cros_degradebeam.pipe
python $PIPELINE_DIR/manager.py $PIPELINE_DIR_P/transfer_cros.pipe
python $PIPELINE_DIR/manager.py $PIPELINE_DIR_P/power_cros.pipe
#python $PIPELINE_DIR/manager.py $PIPELINE_DIR/powespectrum_P/sim_power_cros.pipe

MAP_SIM="/mnt/raid-project/gmrt/ycli/map_sim/1hr/"
MAP_SSM="/mnt/raid-project/gmrt/ycli/foreground_cleand/1hr_IU_legendre_modes_0gwj_conv/"
MAP_CLN="/mnt/raid-project/gmrt/ycli/foreground_cleand/1hr_IU_legendre_modes_0gwj_conv/"
MAP_WGZ="WiggleZ_1hr"
MODE="20"
HOUR="1"
python $PIPELINE_DIR/manager.py $PIPELINE_DIR_P/reference_cros.pipe
python $PIPELINE_DIR/manager.py $PIPELINE_DIR_P/reference_cros_degradebeam.pipe
python $PIPELINE_DIR/manager.py $PIPELINE_DIR_P/transfer_cros.pipe
python $PIPELINE_DIR/manager.py $PIPELINE_DIR_P/power_cros.pipe
#python $PIPELINE_DIR/manager.py $PIPELINE_DIR/powespectrum_P/sim_power_cros.pipe

