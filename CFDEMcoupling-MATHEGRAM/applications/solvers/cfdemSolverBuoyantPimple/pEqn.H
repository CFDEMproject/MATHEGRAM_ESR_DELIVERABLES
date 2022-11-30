rho = thermo.rho();

volScalarField rhoVF = rho*voidfraction; 

// Thermodynamic density needs to be updated by psi*d(p) after the
// pressure solution
const volScalarField psip0(psi*p);

volScalarField rAU(1.0/UEqn.A());
//surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rhoVF*rAU));

volScalarField rhorAUvoidfraction("rhorAUvoidfraction",rho*rAU*voidfraction);
surfaceScalarField rhorAUvoidfractionf
(
    "rhorAUvoidfractionf", fvc::interpolate(rhorAUvoidfraction)
);

//surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p_rgh));

surfaceScalarField phig(-rhorAUvoidfractionf*ghf*fvc::snGrad(rhoVF)*mesh.magSf());
if (modelType=="A")
   phig *= fvc::interpolate(voidfraction);

surfaceScalarField phiUs(fvc::interpolate(rhorAUvoidfraction*Ksl*Us) & mesh.Sf());

surfaceScalarField phiHbyA
(
    "phiHbyA",
    (
        fvc::flux(rho*HbyA*voidfraction)
//      + MRF.zeroFilter(rhorAUf*fvc::ddtCorr(rhoVF, U, phi))
    )
  + phig
  + phiUs
);

MRF.makeRelative(fvc::interpolate(rho), phiHbyA);

// Update the pressure BCs to ensure flux consistency
constrainPressure(p_rgh, rhoVF, U, phiHbyA, rhorAUvoidfractionf, MRF);


fvScalarMatrix p_rghDDtEqn
(
    fvc::ddt(rhoVF)
  + psi*voidfraction*correction(fvm::ddt(p_rgh))
  + fvc::div(phiHbyA)
  ==
    fvOptions(psi, p_rgh, rho.name())
);

if (modelType=="A")
   rhorAUvoidfraction *= voidfraction;
if (modelType=="A")
   rhorAUvoidfractionf *= fvc::interpolate(voidfraction);

while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix p_rghEqn
    (
        p_rghDDtEqn
      - fvm::laplacian(rhorAUvoidfractionf, p_rgh)
    );

    if (thermo.incompressible())
    {
        p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));
    }

    p_rghEqn.solve(mesh.solver(p_rgh.select(pimple.finalInnerIter())));

    if (pimple.finalNonOrthogonalIter())
    {
        // Calculate the conservative fluxes
        phi = phiHbyA + p_rghEqn.flux();
        // Explicitly relax pressure for momentum corrector
        p_rgh.relax();

        // Correct the momentum source with the pressure gradient flux
        // calculated from the relaxed pressure
        U.correctBoundaryConditions();
	if (modelType=="A")
	{
            U 
          = 
	    HbyA 
	  + rAU*( 
                   fvc::reconstruct
                   (
                        (
                            phig 
	                  + fvc::interpolate(voidfraction)
                           *p_rghEqn.flux()
                        )/rhorAUvoidfractionf
                    )
	          + Ksl*Us
               );
	}
	else
	{
            U 
          = 
	    HbyA 
	  + rAU*(
                    fvc::reconstruct
                    (
                        (phig + p_rghEqn.flux())/rhorAUvoidfractionf
                    )
	          + Ksl*Us
                );
	}
        fvOptions.correct(U);
        K = 0.5*magSqr(U);
    }
}

p = p_rgh + voidfraction*rho*gh;

// Thermodynamic density update
thermo.correctRho(psi*p - psip0);

if (thermo.dpdt())
{
    dpdt = fvc::ddt(voidfraction,p);
}

#include "rhoEqn.H"
#include "compressibleContinuityErrs.H"
