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

Class
    Foam::fixedRadHeatFluxFvPatchScalarField

Description
    This boundary condition provides an imposed radiative heat flux boundary 
    condition for P1 approximation.

Usage
    \table
        Property     | Description                | Required    | Default value
        radHeatFlux  | radiative heat flux [W/m2] | yes         |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            fixedRadHeatFlux;
        radHeatFlux     1e3; 
        gradient        uniform 0.0; 
    }
    \endverbatim

See also
    Foam::zeroGradientFvPatchScalarField
    Foam::fixedGradientFvPatchScalarField
    

SourceFiles
    fixedRadHeatFluxFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef fixedRadHeatFluxFvPatchScalarField_H
#define fixedRadHeatFluxFvPatchScalarField_H

#include "fixedGradientFvPatchFields.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class fixedRadHeatFluxFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class fixedRadHeatFluxFvPatchScalarField
:
    public fixedGradientFvPatchScalarField
{
    protected: 

    // Protected data

        //- radHeatFlux Assigned heat flux magnitude [W/m^2]
        autoPtr<Function1<scalar>> radHeatFlux_; 

    public:

       //- Runtime type information
       TypeName("fixedRadHeatFlux");


       // Constructors

        //- Construct from patch and internal field
        fixedRadHeatFluxFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fixedRadHeatFluxFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given fixedRadHeatFluxFvPatchScalarField
        // onto a new patch
        fixedRadHeatFluxFvPatchScalarField
        (
            const fixedRadHeatFluxFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy constructor
        fixedRadHeatFluxFvPatchScalarField
        (
            const fixedRadHeatFluxFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new fixedRadHeatFluxFvPatchScalarField(*this)
            );
        }

        //- Copy constructor setting internal field reference
        fixedRadHeatFluxFvPatchScalarField
        (
            const fixedRadHeatFluxFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new fixedRadHeatFluxFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
