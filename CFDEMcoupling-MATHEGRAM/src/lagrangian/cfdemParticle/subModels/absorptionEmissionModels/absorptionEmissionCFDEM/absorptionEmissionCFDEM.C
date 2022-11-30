/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
                                Copyright 2012-     DCS Computing GmbH, Linz
-------------------------------------------------------------------------------
License
    This file is part of CFDEMcoupling.

    CFDEMcoupling is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    CFDEMcoupling is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with CFDEMcoupling; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
    This code is designed to realize coupled CFD-DEM simulations using ASPHERIX
    and OpenFOAM(R). Note: this code is not part of OpenFOAM(R) (see DISCLAIMER).
\*---------------------------------------------------------------------------*/

#include "absorptionEmissionCFDEM.H"
#include "physicoChemicalConstants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(absorptionEmissionCFDEM, 0);

        addToRunTimeSelectionTable
        (
            absorptionEmissionModel,
            absorptionEmissionCFDEM,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::absorptionEmissionCFDEM::absorptionEmissionCFDEM
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    absorptionEmissionModel(dict, mesh),
    coeffsDict_(dict.subDict(typeName + "Coeffs")),
    apName_
    (
       coeffsDict_.lookupOrDefault<word>("AbsorptionCoeffFieldName","ap")
    ),
    EpName_(coeffsDict_.lookupOrDefault<word>("EmissionFieldName","Ep")),
    useParticleAbsorptionCoefficient_
    (
       coeffsDict_.lookupOrDefault<Switch>
       (
           "useParticleAbsorptionCoefficient",true
       )
    ),
    useParticleEmission_
    (
       coeffsDict_.lookupOrDefault<Switch>("useParticleEmission",true)
    )
{
    Info << "absorptionEmissionCFDEM::absorptionEmissionCFDEM "
         << "useParticleAbsorptionCoefficient_=" 
         << useParticleAbsorptionCoefficient_ << endl;
    Info << "absorptionEmissionCFDEM::absorptionEmissionCFDEM "
         << "useParticleEmission_=" << useParticleEmission_ << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::absorptionEmissionCFDEM::~absorptionEmissionCFDEM()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::radiation::absorptionEmissionCFDEM::aDisp(const label) const
{
    tmp<volScalarField> ta
    (
        new volScalarField
        (
            IOobject
            (
                "a",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("a", dimless/dimLength, 0.0)
        )
    );

    if(useParticleAbsorptionCoefficient_)
    {
        const volScalarField& ap_=mesh_.lookupObject<volScalarField>(apName_);
        ta.ref().primitiveFieldRef()=ap_.primitiveField();
    }
    return ta;
}

Foam::tmp<Foam::volScalarField>
Foam::radiation::absorptionEmissionCFDEM::eDisp(const label bandI) const
{
    tmp<volScalarField> te
    (
        new volScalarField
        (
            IOobject
            (
                "e",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("e", dimless/dimLength, 0.0)
        )
    );

    Warning << "Emission coefficient for dispersed phase"
            << " is not implemented!" << endl;

    return te;
}

Foam::tmp<Foam::volScalarField>
Foam::radiation::absorptionEmissionCFDEM::EDisp(const label bandI) const
{
    tmp<volScalarField> tE
    (
        new volScalarField
        (
            IOobject
            (
                "E",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("E", dimMass/dimLength/pow3(dimTime), 0.0)
        )
    );

    if(useParticleEmission_)
    {
        const volScalarField& Ep_=mesh_.lookupObject<volScalarField>(EpName_);
        tE.ref().primitiveFieldRef()=Ep_.primitiveField();
    }

    // Total emission is 4 times the projected emission
    return 4*tE;
}


// ************************************************************************* //
