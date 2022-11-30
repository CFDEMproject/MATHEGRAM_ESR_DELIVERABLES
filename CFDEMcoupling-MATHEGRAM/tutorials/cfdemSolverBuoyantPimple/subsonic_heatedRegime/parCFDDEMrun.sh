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
headerText="HeatedFlowRegime"
logfileName="log_$headerText"
solverName="cfdemSolverBuoyantPimple"
nrProcs="2"
machineFileName="none"   # yourMachinefileName | none
debugMode="off"          # on | off| strict
separateDEM="in.liggghts_run"
testHarnessPath="$CFDEM_TEST_HARNESS_PATH"
doPost="true"
cleanup="true"
#--------------------------------------------------------------------------------#

#- call function to run a parallel CFD-DEM case
parCFDDEMrun $logpath $logfileName $casePath $headerText $solverName $nrProcs $machineFileName $debugMode $separateDEM

