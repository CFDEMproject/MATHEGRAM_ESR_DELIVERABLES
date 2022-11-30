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
    Foam::radiation::constantScatterCFDEM

Description
    Creates constant radiation scatter coefficient based on CFDEM cloud data. 
    The class is based on OF5.x constantScatter class.

SourceFiles
    scatterModel.C

\*---------------------------------------------------------------------------*/

#ifndef constantScatterCFDEM_H
#define constantScatterCFDEM_H

#include "scatterModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{


/*---------------------------------------------------------------------------*\
                       Class constantScatterCFDEM Declaration
\*---------------------------------------------------------------------------*/

class constantScatterCFDEM
:
    public scatterModel
{

    // Private data

        //- Coefficients dictionary
        dictionary coeffsDict_;

        //- Scattering coefficient for air phase/ [1/m]
        dimensionedScalar sigma_;

        //- Linear-anisotropic phase function coefficient / []
        //  -1 < C < 1
        //  - = backward scattering
        //  0 = isotropic scattering (reasonable default value)
        //  + = forward scattering
        dimensionedScalar C_;

        word sigmapName_; 

        Switch useParticleScatteringCoefficient_;

public:

    //- Runtime type information
    TypeName("constantScatterCFDEM");


    // Constructors

        //- Construct from components
        constantScatterCFDEM(const dictionary& dict, const fvMesh& mesh);


    //- Destructor
    virtual ~constantScatterCFDEM();


    // Member Functions

        //- Return scatter coefficient
        tmp<volScalarField> sigmaEff() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace radiation
} // End namespace Foam

#endif

// ************************************************************************* //
