#!/bin/bash
#
# Launches a simulation reenactment. Call as:
# $> launch_reenactment.sh dir/to/reenact
#
# Mark N. Read, 2018

# Pull directory from commmand line using $1
java -Xms2000m -Xmx2000m -classpath "lib/*:j3dlibs/*" reenact.GUI_Reenact -d $1 -m

# Make movie of the stills generated.
MOVIELOC=$(pwd)
cd $1
$MOVIELOC/makemovie.sh

tar -zcvf stills.tar.gz stills  # Save some space. 
