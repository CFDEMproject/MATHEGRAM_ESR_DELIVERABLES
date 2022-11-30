#!/bin/bash

#===================================================================#
# allrun script for testcase as part of test routine 
# run settlingTest CFD part
# Christoph Goniva - May. 2011
#===================================================================#

#- source CFDEM env vars
. ~/.bashrc

#- include functions
source $CFDEM_SRC_DIR/lagrangian/cfdemParticle/etc/functions.sh

#--------------------------------------------------------------------------------#
#- define variables
casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"
logpath=$casePath
headerText="run_parallel_cfdemSolverRhoPimple_parallelPlates_CFDDEM"
logfileName="log_$headerText"
solverName="cfdemSolverBuoyantPimple" #"cfdemSolverPisoSTM" #
nrProcs="2"
machineFileName="none"   # yourMachinefileName | none
debugMode="off"          # on | off| strict
separateDEM="in.liggghts_run"
testHarnessPath="$CFDEM_TEST_HARNESS_PATH"
runOctave="true"
postproc="false"
cleanup="true"
#--------------------------------------------------------------------------------#

#- call function to run a parallel CFD-DEM case
parCFDDEMrun $logpath $logfileName $casePath $headerText $solverName $nrProcs $machineFileName $debugMode $separateDEM

#------------------------------#
if [ $runOctave == "true" ]
  then

    #- prepare file to be read by octave 
    cd ../DEM/post/
    files=$(find . -name 'dump*')
    for i in ${files}; do
       sed -i "1,9 s/^/%/" $i
    done

    cd ../../CFD/octave/

    #- change path
    cd octave

    #- rmove old graph
    rm *.png

    #- run octave
    octave-cli plotTprofile.m

    #- show plots 
    eog cfdemSolver*.png
    #------------------------------#

fi

#-------------------------------------------------------#
if [ $postproc == "true" ]
  then

    #- keep terminal open (if started in new terminal)
    echo "simulation finished? ...press enter to proceed"
    read

    #- get VTK data from liggghts dump file
    cd $casePath/DEM/post
    python -i $CFDEM_LPP_DIR/lpp.py  dump.liggghts_run

    #- get VTK data from CFD sim
    cd $casePath/CFD
    foamToVTK                                                   #- serial run of foamToVTK
    #source $CFDEM_SRC_DIR/lagrangian/cfdemParticle/etc/functions.sh                       #- include functions
    #pseudoParallelRun "foamToVTK" $nrPostProcProcessors          #- pseudo parallel run of foamToVTK

    #- start paraview
    paraview

    #- keep terminal open (if started in new terminal)
    echo "...press enter to clean up case"
    echo "press Ctr+C to keep data"
    read
fi

#- clean up case
if [ $cleanup == "true" ]
  then
    #- clean up case
    keepDEMrestart="false"
    cleanCFDEMcase $casePath/CFD $keepDEMrestart
fi


