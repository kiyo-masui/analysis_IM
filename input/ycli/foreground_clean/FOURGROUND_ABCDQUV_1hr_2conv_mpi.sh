#! /bin/bash
export PIPELINE_DIR="/home/ycli/pipeline"
#export PIPELINE_DIR="/cita/h/home-2/ycli/analysis_IM/input/ycli/"
cd $PIPELINE_DIR
source $PIPELINE_DIR/PARAMS.sh


#===== Foreground Clean 1hr IxABCDQUV======
export MAPFILE="/mnt/scratch-3week/ycli/1hr_41-18_fdg/"
export NOISEWEIGHT="/mnt/scratch-3week/ycli/1hr_41-18_fdg/"
export GOODMODES="0"
export SIMROOT="/mnt/raid-project/gmrt/ycli/map_sim/1hr_2/"
export SUBREAL=True


#python $PIPELINE_DIR/manager.py $PIPELINE_DIR/$pipeline

export SIMFILE=1

export EXP="-x MAPFILE -x NOISEWEIGHT -x GOODMODES -x SIMROOT -x SUBREAL -x SIMFILE"

export pipeline=foreground_clean/cleaned_AQUV_1hr_extend_2conv_mpi.pipe
mpirun $EXP -np 45 -npernode 15 -hostfile $PIPELINE_DIR/HOSTFILE python $PIPELINE_DIR/mpimanager.py $PIPELINE_DIR/$pipeline

#export pipeline=foreground_clean/cleaned_BQUV_1hr_extend_2conv_mpi.pipe
#mpirun $EXP -np 45 -npernode 15 -hostfile $PIPELINE_DIR/HOSTFILE python $PIPELINE_DIR/mpimanager.py $PIPELINE_DIR/$pipeline
#
#export pipeline=foreground_clean/cleaned_CQUV_1hr_extend_2conv_mpi.pipe
#mpirun $EXP -np 45 -npernode 15 -hostfile $PIPELINE_DIR/HOSTFILE python $PIPELINE_DIR/mpimanager.py $PIPELINE_DIR/$pipeline
#
#export pipeline=foreground_clean/cleaned_DQUV_1hr_extend_2conv_mpi.pipe
#mpirun $EXP -np 45 -npernode 15 -hostfile $PIPELINE_DIR/HOSTFILE python $PIPELINE_DIR/mpimanager.py $PIPELINE_DIR/$pipeline
#
