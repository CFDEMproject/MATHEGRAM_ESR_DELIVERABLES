/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
                                Copyright (C) 1991-2009 OpenCFD Ltd.
                                Copyright (C) 2009-2012 JKU, Linz
                                Copyright (C) 2012-     DCS Computing GmbH,Linz
-------------------------------------------------------------------------------
License
    This file is part of CFDEMcoupling.

    CFDEMcoupling is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    CFDEMcoupling is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with CFDEMcoupling.  If not, see <http://www.gnu.org/licenses/>.

Application
    cfdemSolverBuoyantPimple

Description
    Transient solver for compressible flow.
    The code is an evolution of the solver buoyantPimpleFoam in 
    OpenFOAM(R) 6, with additional functionality for CFD-DEM coupling.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fluidThermo.H" 
#include "turbulentFluidThermoModel.H"
#include "radiationModel.H"


#include "OFversion.H"
//#if defined(version30)
//    #include "turbulentTransportModel.H"
    #include "pimpleControl.H"
//#else
    #include "turbulenceModel.H"
//#endif
#if defined(versionv1606plus) || defined(version40)
    #include "fvOptions.H"
#else
    #include "fvIOoptionList.H"
#endif
#include "fixedFluxPressureFvPatchScalarField.H"

#if defined(MS)
    #if defined(superquadrics_flag)
        #include "cfdemCloudRotationSuperquadric.H"
    #else
        #include "cfdemCloudMS.H"
    #endif
#else
    #include "cfdemCloud.H"
#endif

#include "implicitCouple.H"
#include "clockModel.H"
#include "smoothingModel.H"
#include "forceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCase.H" 
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #if defined(version30)
        //pimpleControl pimple(mesh);
        #include "createTimeControls.H"
    #endif
    #include "createFields.H"
    #include "createFieldRefs.H"
    #include "createFvOptions.H"
    #include "initContinuityErrs.H"

    // create cfdemCloud
    #include "checkImCoupleM.H"

    #if defined(MS)
        #if defined(superquadrics_flag)
            cfdemCloudRotationSuperquadric particleCloud(mesh);
        #else
            cfdemCloudMS particleCloud(mesh);
        #endif
    #else
        cfdemCloud particleCloud(mesh);
    #endif

    #include "checkModelType.H"

    // get ref to tempExchange models
    labelList tempExchangeModels(0);
    label id(particleCloud.registryM().getProperty("LaEuScalarTemp_index"));
    if(id>=0) tempExchangeModels.append(id);
    id = particleCloud.registryM().getProperty("LaEuScalarRadiation_index");
    if(id>=0) tempExchangeModels.append(id);
    
    // check ordering tempExchange models
    if
    (
        tempExchangeModels.size() > 1 
     && tempExchangeModels[0] > tempExchangeModels[1]
    ) 
    {
        FatalError << "Please use correct order of forceModels:" 
           << " LaEuScalarTemp before LaEuScalarRadiation."
           << abort(FatalError);
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
    Info<< "\nStarting time loop\n" << endl;
    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #if defined(version30)
            #include "readTimeControls.H"
            #include "compressibleCourantNo.H"
            #include "setDeltaT.H"
        #else
            #include "readPIMPLEControls.H"
            #include "CourantNo.H"
        #endif

        // do particle stuff
        particleCloud.clockM().start(1,"Global");
        particleCloud.clockM().start(2,"Coupling");
        bool hasEvolved = particleCloud.evolve(voidfraction,Us,U);

        if(hasEvolved && particleCloud.solveFlow())
        {
            particleCloud.smoothingM().smoothenAbsolutField
            (
                particleCloud.forceM(0).impParticleForces()
            );
        }

        Ksl = 
           particleCloud.momCoupleM
           (
               particleCloud.registryM().getProperty("implicitCouple_index")
           ).impMomSource();
        Ksl.correctBoundaryConditions();

        //Force Checks
        #include "forceCheckIm.H"

        particleCloud.clockM().stop("Coupling");

        particleCloud.clockM().start(26,"Flow");

        // get scalar source from DEM
        forAll(tempExchangeModels,i)
        {
            particleCloud.forceM(tempExchangeModels[i]).manipulateScalarField(Qsource);
            Qsource.correctBoundaryConditions();
            particleCloud.forceM(tempExchangeModels[i]).commToDEM();
        }

       	if(particleCloud.solveFlow())
        {

	       #include "rhoEqn.H"
	    
	    // Pressure-velocity PIMPLE corrector
	    while (pimple.loop())
            {
               #include "UEqn.H"
               #include "EEqn.H"   

                // --- Pressure corrector loop
                #if defined(version30)
                    while (pimple.correct())
                #else
                    for (int corr=0; corr<nCorr; corr++)
                #endif
                {
		    #include "pEqn.H"
                } 

               if (pimple.turbCorr())
               {
                  turbulence->correct();
	       }
            } // end pimple loop

            rho = thermo.rho();

        }// end solveFlow
        else
        {
            Info << "skipping flow solution." << endl;
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

        particleCloud.clockM().stop("Flow");
        particleCloud.clockM().stop("Global");
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
