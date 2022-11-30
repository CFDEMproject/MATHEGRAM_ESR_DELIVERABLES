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
    Foam::radiationModels::P1CFDEM

Description
    P1 approximation for particle-fluid systems. 
    Evolved from OF5.x P1 class. 

SourceFiles
    P1CFDEM.C

\*---------------------------------------------------------------------------*/

#ifndef radiationModel_P1CFDEM_H
#define radiationModel_P1CFDEM_H

#include "radiationModel.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                           Class P1CFDEM Declaration
\*---------------------------------------------------------------------------*/

class P1CFDEM
:
    public radiationModel
{
    // Private Data

        //- Incident radiation / [W/m^2]
        volScalarField G_;

        //- Total radiative heat flux at the boubdaries [W/m^2]
        volScalarField qr_;

        //- Radiative heat flux throughout the bed [w/m^2]
        volVectorField radHeatFlux; 

        //- Absorption coefficient
        volScalarField a_;

        //- Emission coefficient
        volScalarField e_;

        //- Emission contribution
        volScalarField E_;

        //-Radiative diffusion coefficient
        volScalarField gamma_; 

public:

    //- Runtime type information
    TypeName("P1CFDEM");


    // Constructors

        //- Construct from components
        P1CFDEM(const volScalarField& T);

        //- Construct from components
        P1CFDEM(const dictionary& dict, const volScalarField& T);

        //- Disallow default bitwise copy construction
        P1CFDEM(const P1CFDEM&) = delete;


    //- Destructor
    virtual ~P1CFDEM();


    // Member Functions

        // Edit

            //- Solve radiation equation(s)
            void calculate();

            //- Read radiation properties dictionary
            bool read();


        // Access

            //- Source term component (for power of T^4)
            virtual tmp<volScalarField> Rp() const;

            //- Source term component (constant)
            virtual tmp<volScalarField::Internal> Ru() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const P1CFDEM&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
