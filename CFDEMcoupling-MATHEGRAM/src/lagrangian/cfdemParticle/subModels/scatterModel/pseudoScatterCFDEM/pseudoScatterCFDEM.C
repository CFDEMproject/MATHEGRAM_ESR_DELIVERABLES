/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "pseudoScatterCFDEM.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(pseudoScatterCFDEM, 0);

        addToRunTimeSelectionTable
        (
            scatterModel,
            pseudoScatterCFDEM,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::pseudoScatterCFDEM::pseudoScatterCFDEM
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
           "useParticleScatteringCoefficient",
           true
        )
    )
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //
Foam::radiation::pseudoScatterCFDEM::~pseudoScatterCFDEM()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
Foam::tmp<Foam::volScalarField>
Foam::radiation::pseudoScatterCFDEM::sigmaEff() const
{

    const volScalarField& voidfraction =
         mesh_.lookupObject<volScalarField>("voidfraction"); 
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

    // Fix with voidfraction
    tsigma.ref().primitiveFieldRef() *= voidfraction.primitiveField();

    if(useParticleScatteringCoefficient_)
    {
        const volScalarField& sigmap_ =
           mesh_.lookupObject<volScalarField>(sigmapName_)*(3.0-C_);
        tsigma.ref().primitiveFieldRef() += sigmap_.primitiveField();
    }

   return tsigma; 
}

// ************************************************************************* //
