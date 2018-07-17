#!/bin/bash


echo "$SHARED_SCRATCH/jdb12/pythia_bbbar/"
echo "JOBID $1"
source /projects/geurts/jb31/ROOT/bin/thisroot.sh

# export ROOTSYS=/projects/geurts/jb31/ROOT
# export PATH=$ROOTSYS/bin:$PATH
# export LD_LIBRARY_PATH=$ROOTSYS/lib/root
#export LD_LIBRARY_PATH=/projects/geurts/jb31/ROOTbuild/pythia6:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/projects/geurts/jb31/ROOT/pythia6/:$LD_LIBRARY_PATH

/home/jdb12/work/PythiaBBbar/GENERATOR 5 200000 >& /dev/null
