#! /bin/bash

export PIPELINE_DIR="/cita/h/home-2/ycli/analysis_IM/input/ycli/"
cd $PIPELINE_DIR
source $PIPELINE_DIR/PARAMS.sh

#pipeline=parkes.pipe
#pipeline=parkes/parkes_dirtymap.pipe
#pipeline=parkes/parkes_cleanmap.pipe
pipeline=parkes/parkes_cleanmap2.pipe

#===== Parkes Analysis ======
python $PIPELINE_DIR/manager.py $PIPELINE_DIR/$pipeline 
export BEAMS=0
python $PIPELINE_DIR/manager.py $PIPELINE_DIR/$pipeline 
export BEAMS=1
python $PIPELINE_DIR/manager.py $PIPELINE_DIR/$pipeline 
export BEAMS=2
python $PIPELINE_DIR/manager.py $PIPELINE_DIR/$pipeline 
export BEAMS=3
python $PIPELINE_DIR/manager.py $PIPELINE_DIR/$pipeline 

