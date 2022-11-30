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
    and OpenFOAM(R). Note: this code is not part of OpenFOAM(R) (see DISCLAIMER)
\*---------------------------------------------------------------------------*/

#include "P1CFDEM.H"
#include "fvCFD.H"
#include "absorptionEmissionModel.H"
#include "scatterModel.H"
#include "constants.H"
#include "addToRunTimeSelectionTable.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
   namespace radiation
   {
       defineTypeNameAndDebug(P1CFDEM, 0);
       addToRadiationRunTimeSelectionTables(P1CFDEM);
   }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::P1CFDEM::P1CFDEM(const volScalarField& T)
:
    radiationModel(typeName, T),
    G_
    (
        IOobject
        (
            "G",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    qr_
    (
        IOobject
        (
            "qr",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("qr", dimMass/pow3(dimTime), 0)
    ),
    radHeatFlux
    (
        IOobject
        (
            "radHeatFlux",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector
        (
             "radHeatFlux", 
             dimMass/pow3(dimTime), 
             Foam::vector(0,0,0)
        )
    ),
    a_
    (
        IOobject
        (
            "a",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("a", dimless/dimLength, 0)
    ),
    e_
    (
        IOobject
        (
            "e",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("e", dimless/dimLength, 0)
    ),
    E_
    (
        IOobject
        (
            "E",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("E", dimMass/dimLength/pow3(dimTime), 0)
    ), 
    gamma_ 
    (
        IOobject
        (
            "gammaRad",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("gammaRad", dimLength, 0)
    )
{
    // change all calculated BC to zeroGradient
    forAll(gamma_.boundaryFieldRef(), patchI)
    {
       if(gamma_.boundaryFieldRef()[patchI].type() == "calculated") 
       {
             gamma_.boundaryFieldRef().set
             (
                  patchI, 
                  fvPatchField<scalar>::New
                  (
                        "zeroGradient",
                        mesh_.boundary()[patchI],
                        gamma_
                  )
             );
       }
    }
}


Foam::radiation::P1CFDEM::P1CFDEM
(
     const dictionary& dict, 
     const volScalarField& T
)
:
    radiationModel(typeName, dict, T),
    G_
    (
        IOobject
        (
            "G",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    qr_
    (
        IOobject
        (
            "qr",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("qr", dimMass/pow3(dimTime), 0)
    ),
    radHeatFlux
    (
        IOobject
        (
            "radHeatFlux",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector
        (
            "radHeatFlux", 
            dimMass/pow3(dimTime), 
            Foam::vector(0,0,0)
        )
    ),
    a_
    (
        IOobject
        (
            "a",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("a", dimless/dimLength, 0)
    ),
    e_
    (
        IOobject
        (
            "e",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("e", dimless/dimLength, 0)
    ),
    E_
    (
        IOobject
        (
            "E",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("E", dimMass/dimLength/pow3(dimTime), 0)
    ), 
    gamma_
    (
        IOobject
        (
            "gammaRad",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("gammaRad", dimLength, 0)
    ) 
{
    // change all calculated BC to zeroGradient
    forAll(gamma_.boundaryFieldRef(), patchI)
    {
       if(gamma_.boundaryFieldRef()[patchI].type() == "calculated") 
       {
             gamma_.boundaryFieldRef().set
             (
                  patchI, 
                  fvPatchField<scalar>::New
                  (
                       "zeroGradient",
                       mesh_.boundary()[patchI], 
                       gamma_
                  )
             );
       }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::P1CFDEM::~P1CFDEM()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::radiation::P1CFDEM::read()
{
    if (radiationModel::read())
    {
        return true;
    }
    else
    {
        return false;
    }
}


void Foam::radiation::P1CFDEM::calculate()
{
    a_ = absorptionEmission_->a();
    e_ = absorptionEmission_->e();
    E_ = absorptionEmission_->E();

    const volScalarField sigmaEff(scatter_->sigmaEff());

    const dimensionedScalar a0 ("a0", a_.dimensions(), ROOTVSMALL);

    // Construct diffusion
    gamma_ = 1.0/(3.0*a_ + sigmaEff + a0);
    gamma_.correctBoundaryConditions(); 

    // Solve G transport equation
    solve
    (
        fvm::laplacian(gamma_, G_)
      - fvm::Sp(a_, G_)
     ==
      - 4.0*(e_*physicoChemical::sigma*pow4(T_) ) - E_
    ); 

    radHeatFlux = G_*fvc::grad(gamma_) - fvc::grad(gamma_*G_); 

    volScalarField::Boundary& qrBf = qr_.boundaryFieldRef();

    // Calculate radiative heat flux on boundaries.
    forAll(mesh_.boundaryMesh(), patchi)
    {
        if (!G_.boundaryField()[patchi].coupled())
        {
            qrBf[patchi] =
                -gamma_.boundaryField()[patchi]
                *G_.boundaryField()[patchi].snGrad();
        }
    }
}


Foam::tmp<Foam::volScalarField> Foam::radiation::P1CFDEM::Rp() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "Rp",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            4.0*absorptionEmission_->eCont()*physicoChemical::sigma
        )
    );

}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::radiation::P1CFDEM::Ru() const
{
    const volScalarField::Internal& G = G_();
    const volScalarField::Internal E = absorptionEmission_->ECont()()();
    const volScalarField::Internal a = absorptionEmission_->aCont()()();

    return a*G - E;
}


// ************************************************************************* //
