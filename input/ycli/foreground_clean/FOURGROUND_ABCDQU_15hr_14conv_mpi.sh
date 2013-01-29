#! /bin/bash
export PIPELINE_DIR="/home/ycli/pipeline"
#export PIPELINE_DIR="/cita/h/home-2/ycli/analysis_IM/input/ycli/"
cd $PIPELINE_DIR
source $PIPELINE_DIR/PARAMS.sh


#===== Foreground Clean 1hr IxABCDQUV======
export MAPFILE="/mnt/raid-project/gmrt/tcv/maps/15hr_41-80_avg_fdgp/"
export GOODMODES="0"
export SIMROOT="/mnt/scratch-gl/ycli/simulation_map/15hr_14/"
export SUBREAL=True

#export pipeline=foreground_clean/cleaned_AQU_15hr_extend_14conv_mpi.pipe
#python $PIPELINE_DIR/manager.py $PIPELINE_DIR/$pipeline
#export pipeline=foreground_clean/cleaned_BQU_15hr_extend_14conv_mpi.pipe
#python $PIPELINE_DIR/manager.py $PIPELINE_DIR/$pipeline
#export pipeline=foreground_clean/cleaned_CQU_15hr_extend_14conv_mpi.pipe
#python $PIPELINE_DIR/manager.py $PIPELINE_DIR/$pipeline
#export pipeline=foreground_clean/cleaned_DQU_15hr_extend_14conv_mpi.pipe
#python $PIPELINE_DIR/manager.py $PIPELINE_DIR/$pipeline

export SIMFILE=1

export EXP="-x MAPFILE -x NOISEWEIGHT -x GOODMODES -x SIMROOT -x SUBREAL -x SIMFILE"

#export pipeline=foreground_clean/cleaned_AQU_15hr_extend_14conv_mpi.pipe
#mpirun $EXP -np 60 -npernode 20 -hostfile $PIPELINE_DIR/HOSTFILE python $PIPELINE_DIR/mpimanager.py $PIPELINE_DIR/$pipeline
#
#export pipeline=foreground_clean/cleaned_BQU_15hr_extend_14conv_mpi.pipe
#mpirun $EXP -np 60 -npernode 20 -hostfile $PIPELINE_DIR/HOSTFILE python $PIPELINE_DIR/mpimanager.py $PIPELINE_DIR/$pipeline
#
#export pipeline=foreground_clean/cleaned_CQU_15hr_extend_14conv_mpi.pipe
#mpirun $EXP -np 60 -npernode 20 -hostfile $PIPELINE_DIR/HOSTFILE python $PIPELINE_DIR/mpimanager.py $PIPELINE_DIR/$pipeline
#
#export pipeline=foreground_clean/cleaned_DQU_15hr_extend_14conv_mpi.pipe
#mpirun $EXP -np 60 -npernode 20 -hostfile $PIPELINE_DIR/HOSTFILE python $PIPELINE_DIR/mpimanager.py $PIPELINE_DIR/$pipeline

export pipeline=foreground_clean/combine_ABCD_14conv.pipe
python $PIPELINE_DIR/manager.py $PIPELINE_DIR/$pipeline
mpirun $EXP -np 60 -npernode 20 -hostfile $PIPELINE_DIR/HOSTFILE python $PIPELINE_DIR/mpimanager.py $PIPELINE_DIR/$pipeline

