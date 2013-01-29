#! /bin/bash

export PIPELINE_DIR="/cita/h/home-2/ycli/analysis_IM/input/ycli/"
cd $PIPELINE_DIR
source $PIPELINE_DIR/PARAMS.sh

#pipeline=parkes/parkes_pre_test.pipe
#pipeline=parkes/parkes_pre_test2.pipe
#pipeline=parkes/parkes_pre_test3.pipe

pipeline=parkes/parkes_pre.pipe
#pipeline=parkes/parkes_dirtymap.pipe
#pipeline=parkes/parkes_cleanmap.pipe

#===== Parkes Analysis ======
#pipeline=parkes/parkes_cleanmap_25.pipe
python $PIPELINE_DIR/manager.py $PIPELINE_DIR/$pipeline 
#export BEAMS=4
#python $PIPELINE_DIR/manager.py $PIPELINE_DIR/$pipeline 
#export BEAMS=5
#python $PIPELINE_DIR/manager.py $PIPELINE_DIR/$pipeline 
#export BEAMS=6
#python $PIPELINE_DIR/manager.py $PIPELINE_DIR/$pipeline 

