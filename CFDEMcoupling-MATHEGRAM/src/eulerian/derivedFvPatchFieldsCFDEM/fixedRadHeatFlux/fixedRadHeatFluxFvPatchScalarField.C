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

#include "zeroGradientFvPatchFields.H"
#include "fixedRadHeatFluxFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "radiationModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fixedRadHeatFluxFvPatchScalarField::
fixedRadHeatFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF), 
    radHeatFlux_()
{}


Foam::fixedRadHeatFluxFvPatchScalarField::
fixedRadHeatFluxFvPatchScalarField
(
    const fixedRadHeatFluxFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper), 
    radHeatFlux_(ptf.radHeatFlux_, false)
{}


Foam::fixedRadHeatFluxFvPatchScalarField::
fixedRadHeatFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF, dict), 
    radHeatFlux_()
{
    radHeatFlux_ = Function1<scalar>::New("radHeatFlux", dict); 
}


Foam::fixedRadHeatFluxFvPatchScalarField::
fixedRadHeatFluxFvPatchScalarField
(
    const fixedRadHeatFluxFvPatchScalarField& tppsf
)
:
    fixedGradientFvPatchScalarField(tppsf), 
    radHeatFlux_(tppsf.radHeatFlux_, false)
{}


Foam::fixedRadHeatFluxFvPatchScalarField::
fixedRadHeatFluxFvPatchScalarField
(
    const fixedRadHeatFluxFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(tppsf, iF), 
    radHeatFlux_(tppsf.radHeatFlux_, false)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fixedRadHeatFluxFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag+1;

    // Diffusion coefficient - created by radiation model's ::updateCoeffs()
    const scalarField& gamma =
        patch().lookupPatchField<volScalarField, scalar>("gammaRad");
    const scalar radHeatFlux = 
        radHeatFlux_->value(this->db().time().timeOutputValue());
    gradient() = radHeatFlux/gamma;  

    UPstream::msgType() = oldTag;

    fixedGradientFvPatchScalarField::updateCoeffs();
}


void Foam::fixedRadHeatFluxFvPatchScalarField::write(Ostream& os) const
{
    fixedGradientFvPatchScalarField::write(os);
    radHeatFlux_->writeData(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        fixedRadHeatFluxFvPatchScalarField
    );
}

// ************************************************************************* //
