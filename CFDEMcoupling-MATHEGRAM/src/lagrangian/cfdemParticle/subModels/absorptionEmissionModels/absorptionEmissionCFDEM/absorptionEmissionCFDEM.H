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
    Foam::radiation::absorptionEmissionCFDEM

Description
    Creates absorption/emission data based on CFDEM cloud data. The class 
    is based on OF3.0.x absorptionEmission for lagrange particles.

SourceFiles
    absorptionEmissionCFDEM.C

\*---------------------------------------------------------------------------*/

#ifndef radiationAbsorptionEmissionCFDEM_H
#define radiationAbsorptionEmissionCFDEM_H

#include "absorptionEmissionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                  Class absorptionEmissionCFDEM Declaration
\*---------------------------------------------------------------------------*/

class absorptionEmissionCFDEM
:
    public absorptionEmissionModel
{
    // Private data

        //- Coefficients dictionary
        dictionary coeffsDict_;

        //- Absorption coefficient field name for dispersed phase
        word apName_;

        //- emission field name for dispersed phase
        word EpName_;

        Switch useParticleAbsorptionCoefficient_;

        Switch useParticleEmission_;


public:

    //- Runtime type information
    TypeName("absorptionEmissionCFDEM");


    // Constructors

        //- Construct from components
        absorptionEmissionCFDEM(const dictionary& dict, const fvMesh& mesh);


    //- Destructor
    virtual ~absorptionEmissionCFDEM();


    // Member Operators

        // Access

            // Absorption coefficient

                //- Absorption coefficient for dispersed phase
                virtual tmp<volScalarField> aDisp(const label bandI = 0) const;


            // Emission coefficient

                //- Emission coefficient for dispersed phase
                virtual tmp<volScalarField> eDisp(const label bandI = 0) const;


            // Emission contribution

                //- Return emission contribution for dispersed phase
                virtual tmp<volScalarField> EDisp(const label bandI = 0) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ///

#endif 
