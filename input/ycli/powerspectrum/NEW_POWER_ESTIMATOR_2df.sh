#! /bin/bash
cd $PBS_O_WORKDIR

export PIPELINE_DIR="/home/ycli/pipeline"
export PIPELINE_DIR_P="/home/ycli/pipeline/powerspectrum"
cd $PIPELINE_DIR
source $PIPELINE_DIR/PARAMS.sh

HOUR="29"
#MAP_2DF="/mnt/scratch-gl/ycli/2df_catalog/map/map_2929.5_large/"
#MAP_2DF="/mnt/scratch-gl/ycli/2df_catalog/map/map_2929.5_selection/"
#MAP_2DF="/mnt/scratch-gl/ycli/2df_catalog/map/map_2929.5_separable/"
#NAME="29RA"

#MAP_2DF="/mnt/scratch-gl/ycli/2df_catalog/map/map_2929.5_full_separable/"
#MAP_2DF="/mnt/scratch-gl/ycli/2df_catalog/map/map_2929.5_full_selection/"
MAP_2DF="/mnt/scratch-gl/ycli/2df_catalog/map/map_2929.5_full_selection_10000mock/"
#MAP_2DF="/mnt/scratch-gl/ycli/2df_catalog/map/map_2929.5_full_oneseed_selection/"
NAME="FULL"

#HOUR="15"
#MAP_2DF="/mnt/raid-project/gmrt/eswitzer/wiggleZ_v1/complete_binned_delta/15hr/"
#MAP_2DF="/mnt/raid-project/gmrt/eswitzer/wiggleZ_v1/binned_delta/15hr/"
#NAME="15hr_wigglez"

export EXP="-x MAP_2DF -x HOUR -x NAME"

#===== Auto  Power ======
mpirun $EXP -np 20  python $PIPELINE_DIR/mpimanager.py $PIPELINE_DIR_P/new_power_2df.pipe
#mpirun $EXP -np 20 -npernode 5 -hostfile $PIPELINE_DIR/HOSTFILE python $PIPELINE_DIR/mpimanager.py $PIPELINE_DIR_P/new_power_2df.pipe
