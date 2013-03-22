#! /bin/bash
cd $PBS_O_WORKDIR

export PIPELINE_DIR="/home/ycli/pipeline"
export PIPELINE_DIR_P="/home/ycli/pipeline/powerspectrum"
cd $PIPELINE_DIR
source $PIPELINE_DIR/PARAMS.sh

export EXP="-x MAP_SIM -x MAP_SSM -x MAP_CLN -x MAP_MULTIPLIER -x MODE -x HOUR"

#===== Auto  Power ======
MAP_SIM="/mnt/scratch-gl/ycli/simulation_map/1hr_14/"
MAP_SSM="/mnt/raid-project/gmrt/ycli/foreground_cleand/1hr_IE_legendre_modes_0gwj_14conv_new/"
MAP_CLN="/mnt/raid-project/gmrt/ycli/foreground_cleand/1hr_IE_legendre_modes_0gwj_14conv_new/"
MAP_MULTIPLIER=1.

#MODE="20"
##HOUR="1"
#mpirun $EXP -np 45 -npernode 15 -hostfile $PIPELINE_DIR/HOSTFILE python $PIPELINE_DIR/mpimanager.py $PIPELINE_DIR_P/power_noise_auto.pipe
#mpirun $EXP -np 45 -npernode 15 -hostfile $PIPELINE_DIR/HOSTFILE python $PIPELINE_DIR/mpimanager.py $PIPELINE_DIR_P/reference_auto.pipe
#mpirun $EXP -np 45 -npernode 15 -hostfile $PIPELINE_DIR/HOSTFILE python $PIPELINE_DIR/mpimanager.py $PIPELINE_DIR_P/reference_auto_degradebeam.pipe
##mpirun -np 20 -npernode 20 python $PIPELINE_DIR/mpimanager.py $PIPELINE_DIR_P/transfer_auto.pipe
#mpirun $EXP -np 45 -npernode 15 -hostfile $PIPELINE_DIR/HOSTFILE python $PIPELINE_DIR/mpimanager.py $PIPELINE_DIR_P/power_auto.pipe
#
MODE="10"
HOUR="1"
mpirun $EXP -np 45 -npernode 15 -hostfile $PIPELINE_DIR/HOSTFILE python $PIPELINE_DIR/mpimanager.py $PIPELINE_DIR_P/power_noise_auto.pipe
mpirun $EXP -np 45 -npernode 15 -hostfile $PIPELINE_DIR/HOSTFILE python $PIPELINE_DIR/mpimanager.py $PIPELINE_DIR_P/power_auto.pipe

#MODE="15"
#HOUR="1"
#mpirun $EXP -np 45 -npernode 15 -hostfile $PIPELINE_DIR/HOSTFILE python $PIPELINE_DIR/mpimanager.py $PIPELINE_DIR_P/power_noise_auto.pipe
#mpirun $EXP -np 45 -npernode 15 -hostfile $PIPELINE_DIR/HOSTFILE python $PIPELINE_DIR/mpimanager.py $PIPELINE_DIR_P/power_auto.pipe
##mpirun $EXP -np 45 -npernode 15 -hostfile $PIPELINE_DIR/HOSTFILE python $PIPELINE_DIR/mpimanager.py $PIPELINE_DIR_P/reference_auto.pipe
##mpirun $EXP -np 45 -npernode 15 -hostfile $PIPELINE_DIR/HOSTFILE python $PIPELINE_DIR/mpimanager.py $PIPELINE_DIR_P/reference_auto_degradebeam.pipe
#
#MODE="20"
#HOUR="1"
#mpirun $EXP -np 45 -npernode 15 -hostfile $PIPELINE_DIR/HOSTFILE python $PIPELINE_DIR/mpimanager.py $PIPELINE_DIR_P/power_noise_auto.pipe
#mpirun $EXP -np 45 -npernode 15 -hostfile $PIPELINE_DIR/HOSTFILE python $PIPELINE_DIR/mpimanager.py $PIPELINE_DIR_P/power_auto.pipe
#
#MODE="25"
#HOUR="1"
#mpirun $EXP -np 45 -npernode 15 -hostfile $PIPELINE_DIR/HOSTFILE python $PIPELINE_DIR/mpimanager.py $PIPELINE_DIR_P/power_noise_auto.pipe
#mpirun $EXP -np 45 -npernode 15 -hostfile $PIPELINE_DIR/HOSTFILE python $PIPELINE_DIR/mpimanager.py $PIPELINE_DIR_P/power_auto.pipe
#
#MODE="30"
#HOUR="1"
#mpirun $EXP -np 45 -npernode 15 -hostfile $PIPELINE_DIR/HOSTFILE python $PIPELINE_DIR/mpimanager.py $PIPELINE_DIR_P/power_noise_auto.pipe
#mpirun $EXP -np 45 -npernode 15 -hostfile $PIPELINE_DIR/HOSTFILE python $PIPELINE_DIR/mpimanager.py $PIPELINE_DIR_P/power_auto.pipe
#
#MODE="35"
#HOUR="1"
#mpirun $EXP -np 45 -npernode 15 -hostfile $PIPELINE_DIR/HOSTFILE python $PIPELINE_DIR/mpimanager.py $PIPELINE_DIR_P/power_noise_auto.pipe
#mpirun $EXP -np 45 -npernode 15 -hostfile $PIPELINE_DIR/HOSTFILE python $PIPELINE_DIR/mpimanager.py $PIPELINE_DIR_P/power_auto.pipe
#
#MODE="40"
#HOUR="1"
#mpirun $EXP -np 45 -npernode 15 -hostfile $PIPELINE_DIR/HOSTFILE python $PIPELINE_DIR/mpimanager.py $PIPELINE_DIR_P/power_noise_auto.pipe
#mpirun $EXP -np 45 -npernode 15 -hostfile $PIPELINE_DIR/HOSTFILE python $PIPELINE_DIR/mpimanager.py $PIPELINE_DIR_P/power_auto.pipe
#
#
