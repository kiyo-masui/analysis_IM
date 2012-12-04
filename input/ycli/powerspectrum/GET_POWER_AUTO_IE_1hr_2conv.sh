#! /bin/bash
cd $PBS_O_WORKDIR

export PIPELINE_DIR="/home/ycli/pipeline"
export PIPELINE_DIR_P="/home/ycli/pipeline/powerspectrum"
cd $PIPELINE_DIR
source $PIPELINE_DIR/PARAMS.sh

#===== Auto  Power ======
MAP_SIM="/mnt/raid-project/gmrt/ycli/map_sim/1hr_2/"

MAP_SSM="/mnt/raid-project/gmrt/ycli/foreground_cleand/1hr_IE_legendre_modes_0gwj_2conv/"
MAP_CLN="/mnt/raid-project/gmrt/ycli/foreground_cleand/1hr_IE_legendre_modes_0gwj_2conv/"

MODE="20"
HOUR="1"
python $PIPELINE_DIR/manager.py $PIPELINE_DIR_P/power_noise_auto.pipe
#python $PIPELINE_DIR/manager.py $PIPELINE_DIR_P/reference_auto.pipe
#python $PIPELINE_DIR/manager.py $PIPELINE_DIR_P/reference_auto_degradebeam.pipe
#python $PIPELINE_DIR/manager.py $PIPELINE_DIR_P/transfer_auto.pipe
#python $PIPELINE_DIR/manager.py $PIPELINE_DIR_P/power_auto.pipe

#MODE="40"
#HOUR="1"
#mpirun -np 30 -pernode python $PIPELINE_DIR/mpimanager.py $PIPELINE_DIR_P/power_noise_auto.pipe
#mpirun -np 30 -pernode python $PIPELINE_DIR/mpimanager.py $PIPELINE_DIR_P/reference_auto.pipe
#mpirun -np 30 -pernode python $PIPELINE_DIR/mpimanager.py $PIPELINE_DIR_P/reference_auto_degradebeam.pipe
#mpirun -np 30 -pernode python $PIPELINE_DIR/mpimanager.py $PIPELINE_DIR_P/transfer_auto.pipe
#mpirun -np 30 -pernode python $PIPELINE_DIR/mpimanager.py $PIPELINE_DIR_P/power_auto.pipe

#MODE="80"
#HOUR="1"
#mpirun -np 30 -pernode python $PIPELINE_DIR/mpimanager.py $PIPELINE_DIR_P/power_noise_auto.pipe
#mpirun -np 30 -pernode python $PIPELINE_DIR/mpimanager.py $PIPELINE_DIR_P/reference_auto.pipe
#mpirun -np 30 -pernode python $PIPELINE_DIR/mpimanager.py $PIPELINE_DIR_P/reference_auto_degradebeam.pipe
#mpirun -np 30 -pernode python $PIPELINE_DIR/mpimanager.py $PIPELINE_DIR_P/transfer_auto.pipe
#mpirun -np 30 -pernode python $PIPELINE_DIR/mpimanager.py $PIPELINE_DIR_P/power_auto.pipe



#===== Auto  Power no conv ======
#MAP_SIM="/mnt/raid-project/gmrt/ycli/map_sim/1hr/"
#MAP_SSM="/mnt/raid-project/gmrt/ycli/foreground_cleand/1hr_ABCD_legendre_modes_0gwj/"
#MAP_CLN="/mnt/raid-project/gmrt/ycli/foreground_cleand/1hr_ABCD_legendre_modes_0gwj/"

#mpirun -np 30 -pernode python $PIPELINE_DIR/mpimanager.py $PIPELINE_DIR_P/power_noise_auto.pipe
#mpirun -np 30 -pernode python $PIPELINE_DIR/mpimanager.py $PIPELINE_DIR_P/reference_auto.pipe
#mpirun -np 30 -pernode python $PIPELINE_DIR/mpimanager.py $PIPELINE_DIR_P/reference_auto_beam.pipe
#mpirun -np 30 -pernode python $PIPELINE_DIR/mpimanager.py $PIPELINE_DIR_P/transfer_auto.pipe
#mpirun -np 30 -pernode python $PIPELINE_DIR/mpimanager.py $PIPELINE_DIR_P/power_auto.pipe

