Info<< "Reading thermophysical properties\n" << endl;
// rhoThermo replaced with fluidThermo
autoPtr<fluidThermo> pThermo(fluidThermo::New(mesh));
fluidThermo& thermo = pThermo();
thermo.validate(args.executable(), "h", "e");

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.rho()
);

volScalarField& p = thermo.p();

Info<< "Reading physical velocity field U" << endl;
Info<< "Note: only if voidfraction at boundary is 1, U is superficial velocity!!!\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "\nReading momentum exchange field Ksl\n" << endl;
volScalarField Ksl
(
    IOobject
    (
        "Ksl",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "\nReading voidfraction field voidfraction = (Vgas/Vparticle)\n" << endl;
volScalarField voidfraction
(
    IOobject
    (
        "voidfraction",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading particle velocity field Us\n" << endl;
volVectorField Us
(
    IOobject
    (
        "Us",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading/calculating face flux field phi\n" << endl;
surfaceScalarField phi
(
    IOobject
    (
       "phi",
       runTime.timeName(),
        mesh,
       IOobject::READ_IF_PRESENT,
       IOobject::AUTO_WRITE
    ),
    linearInterpolate(rho*voidfraction*U) & mesh.Sf()
);

Info<< "Creating turbulence model\n" << endl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);


#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Force p_rgh to be consistent with p
p_rgh = p - rho*voidfraction*gh;

label pRefCell = 0;
scalar pRefValue = 0.0;

if (thermo.incompressible())
{
    setRefCell
    (
        p,
        p_rgh,
        pimple.dict(),
        pRefCell,
        pRefValue
    );
}

mesh.setFluxRequired(p_rgh.name());

Info<< "Creating field dpdt\n" << endl;
volScalarField dpdt
(
    IOobject
    (
        "dpdt",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("dpdt", p.dimensions()/dimTime, 0)
);

Info<< "Creating field kinetic energy K\n" << endl;
volScalarField K("K", 0.5*magSqr(U));

Info<< "\nCreating fluid-particle heat flux field\n" << endl;
volScalarField Qsource
(
    IOobject
    (
        "Qsource",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh//,
    //dimensionedScalar("0", dimensionSet(0, 0, -1, 1, 0), 0.0)
);

#include "createMRF.H"
#include "createRadiationModel.H"

