#! /bin/bash
export PIPELINE_DIR="/home/ycli/pipeline"
#export PIPELINE_DIR="/cita/h/home-2/ycli/analysis_IM/input/ycli/"
cd $PIPELINE_DIR
source $PIPELINE_DIR/PARAMS.sh

#===== Foreground Clean 15hr IxABCDQUV======
#export MAPFILE="/mnt/raid-project/gmrt/tcv/maps/1hr_41-18_avg_fdgp_new/"
export MAPFILE="/mnt/scratch-gl/ycli/maps/1hr_41-18_avg_fdgp/"
export GOODMODES="0"
#export SIMROOT="/mnt/raid-project/gmrt/ycli/map_sim/1hr_ideal/"
#export SIMROOT="/mnt/scratch-gl/ycli/simulation_map/sim_ideal/"
#export SIMROOT="/mnt/scratch-gl/ycli/simulation_map/1hr_11/"
#export SIMROOT="/mnt/scratch-gl/ycli/simulation_map/15hr_14/"
export SIMROOT="/mnt/scratch-gl/ycli/simulation_map/1hr_14/"
export SUBREAL=True

export SIMFILE=1

export EXP="-x MAPFILE -x GOODMODES -x SIMROOT -x SUBREAL -x SIMFILE"

export pipeline=foreground_clean/create_simulation.pipe
mpirun $EXP -np 45 -npernode 15 -hostfile $PIPELINE_DIR/HOSTFILE python $PIPELINE_DIR/mpimanager.py $PIPELINE_DIR/$pipeline


##===== Foreground Clean 1hr IxABCDQUV======
#export MAPFILE="/mnt/scratch-3week/ycli/1hr_41-18_fdg/"
#export NOISEWEIGHT="/mnt/scratch-3week/ycli/1hr_41-18_fdg/"
#export GOODMODES="0"
##export SIMROOT="/mnt/raid-project/gmrt/ycli/map_sim/1hr_ideal/"
##export SIMROOT="/mnt/scratch-gl/ycli/simulation_map/sim_ideal/"
##export SIMROOT="/mnt/scratch-gl/ycli/simulation_map/1hr_11/"
#export SIMROOT="/mnt/scratch-gl/ycli/simulation_map/1hr_20/"
#export SUBREAL=True
#
#export SIMFILE=1
#
#export EXP="-x MAPFILE -x NOISEWEIGHT -x GOODMODES -x SIMROOT -x SUBREAL -x SIMFILE"
#
#export pipeline=foreground_clean/create_simulation.pipe
#mpirun $EXP -np 45 -npernode 15 -hostfile $PIPELINE_DIR/HOSTFILE python $PIPELINE_DIR/mpimanager.py $PIPELINE_DIR/$pipeline
#
