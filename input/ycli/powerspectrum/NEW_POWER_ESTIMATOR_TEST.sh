#! /bin/bash
cd $PBS_O_WORKDIR

export PIPELINE_DIR="/home/ycli/pipeline"
export PIPELINE_DIR_P="/home/ycli/pipeline/powerspectrum"
cd $PIPELINE_DIR
source $PIPELINE_DIR/PARAMS.sh

MAP_SIM="/mnt/scratch-gl/ycli/simulation_map/15hr_14/"
MAP_CLN="/mnt/scratch-gl/ycli/cln_result/15hr_ABCD_legendre_modes_0gwj_14conv_new/Emap_clean_themselves/"
MAP_SSM="/mnt/scratch-gl/ycli/cln_result/15hr_ABCD_legendre_modes_0gwj_14conv_new/simmapmode_simmap_beam_%03d_subreal/"
HOUR="15"

export EXP="-x MAP_SIM -x MAP_SSM -x MAP_CLN -x MAP_MULTIPLIER -x MODE -x HOUR"

#===== Auto  Power ======
MODE="10"
mpirun $EXP -np 45 -npernode 15 -hostfile $PIPELINE_DIR/HOSTFILE python $PIPELINE_DIR/mpimanager.py $PIPELINE_DIR_P/new_power_auto.pipe
MODE="15"
mpirun $EXP -np 45 -npernode 15 -hostfile $PIPELINE_DIR/HOSTFILE python $PIPELINE_DIR/mpimanager.py $PIPELINE_DIR_P/new_power_auto.pipe
MODE="20"
mpirun $EXP -np 45 -npernode 15 -hostfile $PIPELINE_DIR/HOSTFILE python $PIPELINE_DIR/mpimanager.py $PIPELINE_DIR_P/new_power_auto.pipe
MODE="25"
mpirun $EXP -np 45 -npernode 15 -hostfile $PIPELINE_DIR/HOSTFILE python $PIPELINE_DIR/mpimanager.py $PIPELINE_DIR_P/new_power_auto.pipe
MODE="30"
mpirun $EXP -np 45 -npernode 15 -hostfile $PIPELINE_DIR/HOSTFILE python $PIPELINE_DIR/mpimanager.py $PIPELINE_DIR_P/new_power_auto.pipe
MODE="35"
mpirun $EXP -np 45 -npernode 15 -hostfile $PIPELINE_DIR/HOSTFILE python $PIPELINE_DIR/mpimanager.py $PIPELINE_DIR_P/new_power_auto.pipe
MODE="40"
mpirun $EXP -np 45 -npernode 15 -hostfile $PIPELINE_DIR/HOSTFILE python $PIPELINE_DIR/mpimanager.py $PIPELINE_DIR_P/new_power_auto.pipe
