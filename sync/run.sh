#!/bin/bash
export PARAVIEW=/path/to/paraview
export PATH=$PARAVIEW/bin:$PATH
export LD_LIBRARY_PATH=$PARAVIEW/lib:$LD_LIBRARY_PATH

export PYTHONPATH=$PARAVIEW/lib/python3.8/site-packages:$PYTHONPATH
export PYTHONPATH=$PARAVIEW/lib/python3.8/site-packages/paraview/:$PYTHONPATH

export LD_LIBRARY_PATH=/path/to/pythonlib
export PYTHONPATH=/path/to/pythonlib:$PATH:$LD_LIBRARY_PATH:`pwd`'/'

echo tgv > SESSION.NAME
echo `pwd`'/' >> SESSION.NAME
rm -rf perf
mkdir perf 
rm -rf fig
mkdir fig

mpirun -n 4 nek5000