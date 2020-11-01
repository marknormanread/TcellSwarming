#!/bin/bash

# This script will create a Jar file holding all of the class files required to run the Treg_2D simulation.
# The Jar file is called 'Treg_2D.jar' and will be held in the Treg_2D directory. 
# 
# Compilation requires the creation of a temporary directory, and this is cleaned up at the end of the script.
# The script means that bin directories no longer need to be held in the svn repository, and eclipse does not
# require customisation to handle svn folders that exist in the bin directories. Instead the jar file can be
# included in the svn repository. Remember to call this script before commiting changes to svn though!
# 

if [ -d "bin" ] ; then					# check that the bin directory exists. abort if it does not. 
	# do nothing.
	echo "bin directory exists." 
else
	echo "bin directory does not exist. Jar file cannot be created without it. Aborting!"
	exit 1
fi

cd bin
find . -name "*.svn*" -exec rm -rvf {} \;	# delete all svn files from the bin directory, svn should only be used to track source. 

jar cvf ../simulation.jar *	# create a jar containing everything in this temporary direcetory, and place the jar one level up. 
cd ..				# cd back out of the temporary dir. 
mv simulation.jar lib/.

echo "jar compilation complete." 
