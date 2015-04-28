#! /bin/bash
export PIPELINE_DIR="/home/ycli/pipeline"
#export PIPELINE_DIR="/cita/h/home-2/ycli/analysis_IM/input/ycli/"
cd $PIPELINE_DIR
source $PIPELINE_DIR/PARAMS.sh


#===== Foreground Clean 1hr IxABCDQUV======
export MAPFILE="/mnt/scratch-3week/ycli/1hr_41-18_fdg/"
export NOISEWEIGHT="/mnt/scratch-3week/ycli/1hr_41-18_fdg/"
export GOODMODES="0"
#export SIMROOT="/mnt/raid-project/gmrt/ycli/map_sim/1hr_2/"
export SIMROOT="/mnt/scratch-gl/ycli/simulation_map/1hr_20/"
export SUBREAL=True

export pipeline=foreground_clean/cleaned_AQU_1hr_extend_2conv_mpi.pipe
python $PIPELINE_DIR/manager.py $PIPELINE_DIR/$pipeline
#export pipeline=foreground_clean/cleaned_BQU_1hr_extend_2conv_mpi.pipe
#python $PIPELINE_DIR/manager.py $PIPELINE_DIR/$pipeline
#export pipeline=foreground_clean/cleaned_CQU_1hr_extend_2conv_mpi.pipe
#python $PIPELINE_DIR/manager.py $PIPELINE_DIR/$pipeline
#export pipeline=foreground_clean/cleaned_DQU_1hr_extend_2conv_mpi.pipe
#python $PIPELINE_DIR/manager.py $PIPELINE_DIR/$pipeline
#
#export SIMFILE=1
#
#export EXP="-x MAPFILE -x NOISEWEIGHT -x GOODMODES -x SIMROOT -x SUBREAL -x SIMFILE"
#
#export pipeline=foreground_clean/cleaned_AQU_1hr_extend_2conv_mpi.pipe
#mpirun $EXP -np 60 -npernode 20 -hostfile $PIPELINE_DIR/HOSTFILE python $PIPELINE_DIR/mpimanager.py $PIPELINE_DIR/$pipeline
#
#export pipeline=foreground_clean/cleaned_BQU_1hr_extend_2conv_mpi.pipe
#mpirun $EXP -np 60 -npernode 20 -hostfile $PIPELINE_DIR/HOSTFILE python $PIPELINE_DIR/mpimanager.py $PIPELINE_DIR/$pipeline
#
#export pipeline=foreground_clean/cleaned_CQU_1hr_extend_2conv_mpi.pipe
#mpirun $EXP -np 60 -npernode 20 -hostfile $PIPELINE_DIR/HOSTFILE python $PIPELINE_DIR/mpimanager.py $PIPELINE_DIR/$pipeline
#
#export pipeline=foreground_clean/cleaned_DQU_1hr_extend_2conv_mpi.pipe
#mpirun $EXP -np 60 -npernode 20 -hostfile $PIPELINE_DIR/HOSTFILE python $PIPELINE_DIR/mpimanager.py $PIPELINE_DIR/$pipeline
#
#export pipeline=foreground_clean/combine_ABCD_2conv.pipe
#python $PIPELINE_DIR/manager.py $PIPELINE_DIR/$pipeline
#mpirun $EXP -np 60 -npernode 20 -hostfile $PIPELINE_DIR/HOSTFILE python $PIPELINE_DIR/mpimanager.py $PIPELINE_DIR/$pipeline
#
