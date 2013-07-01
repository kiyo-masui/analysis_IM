#! /bin/bash
cd $PBS_O_WORKDIR

export PIPELINE_DIR="/home/ycli/pipeline"
export PIPELINE_DIR_P="/home/ycli/pipeline/powerspectrum"
cd $PIPELINE_DIR
source $PIPELINE_DIR/PARAMS.sh

HOUR="15"
MAP_SIM="/mnt/scratch-gl/ycli/simulation_map/15hr_14/"
#MAP_CLN="/mnt/raid-project/gmrt/ycli/foreground_cleand/15hr_IE_14conv/mapmode_map/"
#MAP_SSM="/mnt/raid-project/gmrt/ycli/foreground_cleand/15hr_IE_14conv/simmapmode_simmap_beam_%03d_subreal/"
#MAP_SSM="/mnt/raid-project/gmrt/ycli/foreground_cleand/15hr_IE_14conv/simmapmode_simmap_beam_%03d/"
#MAP_SSM="/mnt/raid-project/gmrt/ycli/foreground_cleand/15hr_IE_14conv/simmapmode_simmap_beam_%03d_x0.40/"
#NAME="15hr_IE_14conv"

#MAP_WGZ="/mnt/raid-project/gmrt/eswitzer/wiggleZ/binned/15hr/"
#MAP_CLN="/mnt/scratch-gl/ycli/cln_result/15hr_ABCD_legendre_modes_0gwj_14conv_new/mapmode_map/"
#MAP_SSM="/mnt/scratch-gl/ycli/cln_result/15hr_ABCD_legendre_modes_0gwj_14conv_new/simmapmode_simmap_beam_%03d_subreal/"
#NAME="15hr_ABCD_14conv"

#MAP_WGZ="/mnt/raid-project/gmrt/eswitzer/wiggleZ/binned/15hr/"
#MAP_CLN="/mnt/raid-project/gmrt/eswitzer/GBT/cleaned_maps/GBT_15hr_map_autopaper/"
#MAP_SSM="/mnt/scratch-gl/ycli/cln_result/15hr_ABCD_legendre_modes_0gwj_14conv_new/simmapmode_simmap_beam_%03d_subreal/"
#
#NAME="15hr_ABCD_14conv_test_eric"

MAP_WGZ="/mnt/raid-project/gmrt/eswitzer/wiggleZ/binned/15hr_delta/"
MAP_CLN="/mnt/scratch-gl/ycli/cln_result/15hr_ABCD_14conv/Emap_clean_themselves/"
MAP_SSM="/mnt/scratch-gl/ycli/cln_result/15hr_ABCD_14conv/simmapmode_simmap_beam_%03d_subreal/"

NAME="15hr_ABCD_14conv_test_me"

#HOUR="1"
#MAP_SIM="/mnt/scratch-gl/ycli/simulation_map/1hr_14/"
#MAP_CLN="/mnt/raid-project/gmrt/ycli/foreground_cleand/1hr_IE_14conv/mapmode_map/"
#MAP_SSM="/mnt/raid-project/gmrt/ycli/foreground_cleand/1hr_IE_14conv/simmapmode_simmap_beam_%03d_subreal/"
#NAME="1hr_IE_14conv"




export EXP="-x MAP_SIM -x MAP_SSM -x MAP_CLN -x MAP_WGZ -x MAP_MULTIPLIER -x MODE -x HOUR -x NAME"

#===== Auto  Power ======
#MODE="10"
#mpirun $EXP -np 45 -npernode 15 -hostfile $PIPELINE_DIR/HOSTFILE python $PIPELINE_DIR/mpimanager.py $PIPELINE_DIR_P/new_power_auto.pipe
#MODE="15"
#mpirun $EXP -np 45 -npernode 15 -hostfile $PIPELINE_DIR/HOSTFILE python $PIPELINE_DIR/mpimanager.py $PIPELINE_DIR_P/new_power_auto.pipe
MODE="20"
mpirun $EXP -np 45 -npernode 15 -hostfile $PIPELINE_DIR/HOSTFILE python $PIPELINE_DIR/mpimanager.py $PIPELINE_DIR_P/new_power_auto.pipe
#MODE="25"
#mpirun $EXP -np 45 -npernode 15 -hostfile $PIPELINE_DIR/HOSTFILE python $PIPELINE_DIR/mpimanager.py $PIPELINE_DIR_P/new_power_auto.pipe
#MODE="30"
#mpirun $EXP -np 45 -npernode 15 -hostfile $PIPELINE_DIR/HOSTFILE python $PIPELINE_DIR/mpimanager.py $PIPELINE_DIR_P/new_power_auto.pipe
#MODE="35"
#mpirun $EXP -np 45 -npernode 15 -hostfile $PIPELINE_DIR/HOSTFILE python $PIPELINE_DIR/mpimanager.py $PIPELINE_DIR_P/new_power_auto.pipe
#MODE="40"
#mpirun $EXP -np 45 -npernode 15 -hostfile $PIPELINE_DIR/HOSTFILE python $PIPELINE_DIR/mpimanager.py $PIPELINE_DIR_P/new_power_auto.pipe
