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
    This code is designed to realize coupled CFD-DEM simulations using LIGGGHTS
    and OpenFOAM(R). Note: this code is not part of OpenFOAM(R) (see DISCLAIMER).
    two way LaEu Scalar Exchange between DEM and CFD radiative heat transfer 
    model. 
Class
    LaEuScalarRadiation
SourceFiles
    LaEuScalarRadiation.C
\*---------------------------------------------------------------------------*/

#ifndef LaEuScalarRadiation_H
#define LaEuScalarRadiation_H

#include "forceModel.H"
#include "averagingModel.H"
#include "interpolationCellPoint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class LaEuScalarRadiation Declaration
\*---------------------------------------------------------------------------*/

class LaEuScalarRadiation
:
    public forceModel
{
protected:

    dictionary propsDict_;

    //- Absorption coefficient field for dispersed phase
    mutable volScalarField ap_;

    //- emission field for dispersed phase
    mutable volScalarField Ep_;

    //- Scattering field for dispersed phase
    mutable volScalarField sigmap_;

    //- Incident radiation intensity
    const volScalarField& G_;

    //- Absorption efficiency (=emissivity of material for large particles)
    scalar Qabs_;

    //- Scattering efficiency (=reflectivity of material for large particles)
    scalar Qsca_;

    //- Effects of dependent scattering/shadowing
    Switch useBrewsterCorrection_;

    //- Smoothing of fields
    Switch useSmoothing_;

    word voidfractionFieldName_;

    const volScalarField& voidfraction_; // ref to voidfraction field

    word partTempName_;

    mutable double **partTemp_;                  // Lagrangian array

    word partHeatFluxName_;

    mutable double **partHeatFlux_;              // Lagrangian array

    void allocateMyArrays() const;

public:

    //- Runtime type information
    TypeName("LaEuScalarRadiation");

    // Constructors

        //- Construct from components
        LaEuScalarRadiation
        (
            const dictionary& dict,
            cfdemCloud& sm,
            word name = word("")
        );

    // Destructor

        ~LaEuScalarRadiation();


    // Member Functions
        void setForce() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
