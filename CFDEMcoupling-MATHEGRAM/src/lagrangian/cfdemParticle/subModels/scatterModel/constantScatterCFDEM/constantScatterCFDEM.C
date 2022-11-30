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

#include "OFversion.H"
#include "constantScatterCFDEM.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(constantScatterCFDEM, 0);

        addToRunTimeSelectionTable
        (
            scatterModel,
            constantScatterCFDEM,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::constantScatterCFDEM::constantScatterCFDEM
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    scatterModel(dict, mesh),
    coeffsDict_(dict.optionalSubDict(typeName + "Coeffs")),
    sigma_("sigma", dimless/dimLength, coeffsDict_),
    C_("C", dimless, coeffsDict_), 
    sigmapName_
    (
       coeffsDict_.lookupOrDefault<word>("ScatteringCoeffFieldName","sigmap")
    ),
    useParticleScatteringCoefficient_
    (
       coeffsDict_.lookupOrDefault<Switch>
       (
           "useParticleScatteringCoefficient",true
        )
    )
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //
Foam::radiation::constantScatterCFDEM::~constantScatterCFDEM()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
Foam::tmp<Foam::volScalarField>
Foam::radiation::constantScatterCFDEM::sigmaEff() const
{
    tmp<volScalarField> tsigma
    (
        new volScalarField
        (
            IOobject
            (
                "sigma",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            3.0*sigma_
        )
    );

    if(useParticleScatteringCoefficient_)
    {
        const volScalarField& sigmap_ = 
            mesh_.lookupObject<volScalarField>(sigmapName_) * (3.0-C_);
        tsigma.ref().primitiveFieldRef() += sigmap_.primitiveField();
    }

   return tsigma; 
}

// ************************************************************************* //
