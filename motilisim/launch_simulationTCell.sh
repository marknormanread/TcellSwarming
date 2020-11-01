#!/bin/bash
# 
# Launches a simulation T cell experiment
# $> launch_simulationTCell.sh to/directory
#
# Assumes that "to/directory" contains parameters.xml, and this is also where 
# data will be written.
#
# Mark N. Read, 2018

# Pull directory from commmand line using $1
OUTDIR=$1/rep$2
java -Xms3000m -Xmx3000m -classpath "lib/*:j3dlibs/*" core.SimulationTCell -p $1/parameters.xml -o $OUTDIR -s $2
#python plot_reenactment.py $1
#./makemovie_supply_dir.sh $OUTDIR
