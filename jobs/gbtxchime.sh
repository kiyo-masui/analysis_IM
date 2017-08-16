#!/bin/bash
#PBS -l nodes=1:m32g:ppn=8
#PBS -l walltime=0:40:00
#PBS -N full522

cd $SCRATCH/parkes/analysis_IM

source ../setup2

export MAP_GBT='/scratch2/p/pen/pberger/ch_pipeline_runs/pass2_k/chimexgbt/regridded_maps/regridded_gbtmap_gbtvec_map_I_'
export MAP_CHIME='/scratch2/p/pen/pberger/ch_pipeline_runs/pass2_k/chimexgbt/regridded_maps/regridded_chimemap_gbtvec_map_I_'

export PKS_FIELD='182'
export OPT_FIELD='p1820'
export MODE=10
export NAME='full312'
#mpirun -npernode 5 python pipeline/manager.py pipeline/new_power_px2.pipe

export PKS_FIELD='n18'
export OPT_FIELD='n1800'
export MODE=10
export NAME='full522'
mpirun -npernode 5 python pipeline/manager.py pipeline/gbt_power.pipe
#mpirun -npernode 5 python pipeline/manager.py pipeline/red_power.pipe

export PKS_FIELD='165'
export OPT_FIELD='p1650'
#export MODE=15
export NAME='ful510'
#mpirun -npernode 5 python pipeline/manager.py pipeline/new_power_px2.pipe
#mpirun -npernode 5 python pipeline/manager.py pipeline/red_power.pipe

export PKS_FIELD='199'
export OPT_FIELD='p1990'
#export MODE=20
export NAME='full511'
#mpirun -npernode 5 python pipeline/manager.py pipeline/new_power_px2.pipe
#mpirun -npernode 5 python pipeline/manager.py pipeline/red_power.pipe

export PKS_FIELD='33'
export OPT_FIELD='p3300'
export NAME='full512'
#mpirun -npernode 5 python pipeline/manager.py pipeline/new_power_px2.pipe
#mpirun -npernode 5 python pipeline/manager.py pipeline/red_power.pipe

export PKS_FIELD='225'
export OPT_FIELD='p2250'
#export MODE=25
export NAME='full435'
#mpirun -npernode 5 python pipeline/manager.py pipeline/new_power_px2.pipe

