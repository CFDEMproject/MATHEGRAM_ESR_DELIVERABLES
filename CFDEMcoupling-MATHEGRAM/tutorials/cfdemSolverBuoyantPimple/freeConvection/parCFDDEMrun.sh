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
headerText="run_parallel_cfdemSolverPisoScalar_freeConvection_CFDDEM"
logfileName="log_$headerText"
solverName="cfdemSolverBuoyantPimple"
nrProcs="4"
machineFileName="none"   # yourMachinefileName | none
debugMode="off"          # on | off| strict
testHarnessPath="$CFDEM_TEST_HARNESS_PATH"
runOctave="true"
postproc="false"
cleanup="true"
#--------------------------------------------------------------------------------#

#- call function to run a parallel CFD-DEM case
parCFDDEMrun $logpath $logfileName $casePath $headerText $solverName $nrProcs $machineFileName $debugMode

#------------------------------#

#------------------------------#
if [ $runOctave == "true" ]
  then

    #- change path
    cd octave

    #- rmove old graph
    rm *.png

    #- run octave
    octave-cli localNu.m

    #- show plots 
    eog cfdemSolver*.png
    #------------------------------#

fi


#- clean up case
if [ $cleanup == "true" ]
  then
    keepDEMrestart="false"
    cleanCFDEMcase $casePath/CFD $keepDEMrestart
fi


