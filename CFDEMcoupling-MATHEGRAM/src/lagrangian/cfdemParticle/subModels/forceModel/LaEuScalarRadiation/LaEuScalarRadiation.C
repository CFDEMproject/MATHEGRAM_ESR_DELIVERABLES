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
\*---------------------------------------------------------------------------*/

#include "error.H"

#include "LaEuScalarRadiation.H"
#include "addToRunTimeSelectionTable.H"
#include "dataExchangeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(LaEuScalarRadiation, 0);

addToRunTimeSelectionTable
(
    forceModel,
    LaEuScalarRadiation,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
LaEuScalarRadiation::LaEuScalarRadiation
(
    const dictionary& dict,
    cfdemCloud& sm,
    word name
)
:
    forceModel(dict,sm),
    propsDict_(dict.subDict(name == "" ? typeName + "Props" : name + "Props")),
    ap_
    (
        IOobject
        (
            "ap",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        sm.mesh(),
        dimensionedScalar("zero", dimless/dimLength, 0.) 
    ),
    Ep_
    (
        IOobject
        (
            "Ep",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        sm.mesh(),
        dimensionedScalar("zero", dimMass/dimLength/pow3(dimTime), 0.) 
    ),
    sigmap_
    (
        IOobject
        (
            "sigmap",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        sm.mesh(),
        dimensionedScalar("zero", dimless/dimLength, 0.)
    ),
    G_(sm.mesh().lookupObject<volScalarField>("G")),
    Qabs_(readScalar(propsDict_.lookup("Qabs"))),
    Qsca_(readScalar(propsDict_.lookup("Qsca"))),
    useBrewsterCorrection_
    (
        propsDict_.lookupOrDefault<Switch>("useBrewsterCorrection",true)
    ),
    useSmoothing_(propsDict_.lookupOrDefault<Switch>("useSmoothing",true)),
    voidfractionFieldName_
    (
        propsDict_.lookupOrDefault<word>
       (
          "voidfractionFieldName", 
          "voidfraction"
       )
    ),
    voidfraction_
    (
       sm.mesh().lookupObject<volScalarField>(voidfractionFieldName_)
    ),
    partTempName_(propsDict_.lookupOrDefault<word>("partTempName","Temp")),
    partTemp_(NULL),
    partHeatFluxName_
    (
      propsDict_.lookupOrDefault<word>
      (
         "partHeatFluxName", 
         "convectiveHeatFlux"
      )
    ),
    partHeatFlux_(NULL)
{
    allocateMyArrays();

    // init force sub model
    setForceSubModels(propsDict_);

    // define switches which can be read from dict
    // activate search for verbose switch
    forceSubM(0).setSwitchesList(3,true); 
    // activate search for interpolate switch
    forceSubM(0).setSwitchesList(4,true);

    // read those switches defined above, if provided in dict
    for (int iFSub=0;iFSub<nrForceSubModels();iFSub++)
        forceSubM(iFSub).readSwitches();

    particleCloud_.checkCG(true);

    Warning << "\nLaEuScalarRadiation does yet not work "
            << "in combination with LaEuScalarTemp "
            << "(use either of both)\n" << endl;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

LaEuScalarRadiation::~LaEuScalarRadiation()
{
    delete partTemp_;
    delete partHeatFlux_;
}

// * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * //
void LaEuScalarRadiation::allocateMyArrays() const
{
    // get memory for 2d arrays
    double initVal=0.0;
    // field/initVal/with/lenghtFromLigghts
    particleCloud_.dataExchangeM().allocateArray(partTemp_,initVal,1);  
    particleCloud_.dataExchangeM().allocateArray(partHeatFlux_,initVal,1);
}
// * * * * * *  * * * * * * public Member Functions  * * * * * * * * * * * * //

void LaEuScalarRadiation::setForce() const
{
    Info << "\nLaEuScalarRadiation::setForce()\n" << endl;
    // realloc the arrays
    allocateMyArrays();

    // get DEM data
    particleCloud_.dataExchangeM().getData(partTempName_,"scalar-atom",partTemp_);

    // reset ap_
    forAll(ap_,cellI)
    {
        ap_[cellI] = 0.;
        Ep_[cellI] = 0.;
        sigmap_[cellI] = 0.;
    }
    ap_.correctBoundaryConditions();
    Ep_.correctBoundaryConditions();

    // set per particle absorbtion and emission
    scalar cellI(0);
    vector position(0,0,0);
    scalar Vc(0);
    scalar dscaled(0);
    scalar dparcel(0);
    scalar numberParticlesInParcel(0);
    scalar dt(particleCloud_.mesh().time().deltaT().value());
    scalar AsProj(0);
    scalar T(0);
    scalar ap(0);
    scalar Ep(0);
    scalar sigmap(0);
    scalar G(0);

    #include "resetGInterpolator.H"

    for(int index = 0;index < particleCloud_.numberOfParticles(); ++index)
    {
        cellI = particleCloud_.cellIDs()[index][0];
        if(cellI >= 0)
        {
            if(forceSubM(0).interpolation())
            {
                position = particleCloud_.position(index);
                G = GInterpolator_().interpolate(position,cellI);
            }
            else
            {
                G = G_[cellI];
            }

            // cell volume
            Vc = particleCloud_.mesh().V()[cellI];

            // projected area of particles in parcel
            dscaled = 2*particleCloud_.radius(index);
            dparcel = dscaled;
            //caution: this fct will scale ds!
            forceSubM(0).scaleDia(dscaled,index); 
            numberParticlesInParcel = dparcel/dscaled;
            numberParticlesInParcel *= 
               numberParticlesInParcel*numberParticlesInParcel;
            AsProj = dscaled*dscaled*M_PI*numberParticlesInParcel/4.;

            // calc absorbtion
            ap = AsProj*Qabs_; // [m2]
            if(useBrewsterCorrection_)
                ap /= voidfraction_[cellI];
            ap_[cellI] += ap/Vc; // [1/m]

            // calc scattering 
            sigmap = AsProj*Qsca_; // [m2]
            if(useBrewsterCorrection_)
                sigmap /= voidfraction_[cellI];
            sigmap_[cellI] += sigmap/Vc; // [1/m]

            // calc emission
            T = partTemp_[index][0];
            // [m2 *K4 *kg/s3/K4] = [m2*kg/s3] = [W]
            Ep = ap*T*T*T*T*Foam::constant::physicoChemical::sigma.value(); 
            Ep_[cellI] += Ep/Vc;    // [kg/s3/m] = [W/m3]

            // calc convective heat flux [W]
            scalar partHeatFlux     = G*ap - 4*Ep; // [m2*kg/s3] = [W]
            partHeatFlux_[index][0] += partHeatFlux;

            if(forceSubM(0).verbose() && index >=0 && index <2)
            {
                const volScalarField& Tfluid = 
                    particleCloud_.mesh().lookupObject<volScalarField> ("T");
                Pout << "dscaled = " << dscaled << endl;
                Pout << "dparcel = " << dparcel << endl;
                Pout << "AsProj = " << AsProj << endl;
                Pout << "T = " << T << endl;
                Pout << "Tfluid[cellI] = " << Tfluid[cellI] << endl;
                Pout << "G = " << G << endl;
                Pout << "ap = " << ap << endl;
                Pout << "G*ap = " << G*ap << endl;
                Pout << "Ep = " << Ep << endl;
                Pout << "partTemp_[index][0] = " 
                     << partTemp_[index][0] << endl;
                Pout << "partHeatFlux_[index][0] = " 
                     << partHeatFlux_[index][0] << endl;
            }
        }
    }
    // give DEM data
    particleCloud_.dataExchangeM().giveData
    (
        partHeatFluxName_,
        "scalar-atom", 
        partHeatFlux_
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
