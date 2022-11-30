#!/bin/bash

#===================================================================#
# allrun script for testcase as part of test routine 
# run settlingTest
# Christoph Goniva - August 2011
#===================================================================#



#- define variables
casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"

#- include functions
source $CFDEM_SRC_DIR/lagrangian/cfdemParticle/etc/functions.sh
#- define variables
runOctave="true"
showPlot="false"
cleanFig="true"
cleanup="true"

casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"
ctrDctPath="$casePath/CFD/system/controlDict"
caseName=${PWD##*/}
testHarnessPath="$CFDEM_TEST_HARNESS_PATH"
casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"
separateDEM="run.asx"

#- clean the old data
if [ $cleanFig == "true" ]; then 
    rm $casePath/CFD/octave/*.png
fi

 cd $casePath/CFD
 blockMesh
 cd ../..
 #- run parallel CFD-DEM in new terminal
 . $casePath/parCFDDEMrun.sh
   
# Postprocessing
postProcPath="$casePath/CFD/postProcessing/singleGraph/"
endTime=$(ls -t $postProcPath | head -1 )
echo "Latest timestep:" $endTime

postProcFile="\"..\/postProcessing\/singleGraph\/$endTime\/line_mag(U)_p_rho_T_voidfraction.xy\""

#- octave
if [ $runOctave == "true" ]
    then
        #- change path
        cd octave
        sed -i "s/^path\ .*/path\ =\ $postProcFile;/" main.m

        #- run octave
        octave --no-gui main.m

        #- show plot
        eog cfdemSolverRhoPimple_nozzleFlow.png

        #- copy log file to test harness
        cp ../../$logfileName $testHarnessPath
        cp *.png $testHarnessPath
fi

#- clean up case
if [ $cleanup == "true" ]; then
    keepDEMrestart="true"
    cleanCFDEMcase $casePath/CFD $keepDEMrestart
fi

#------------------------------------------------------------------------------

