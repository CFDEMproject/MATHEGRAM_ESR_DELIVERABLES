/* ----------------------------------------------------------------------
    This is the

    ██╗     ██╗ ██████╗  ██████╗  ██████╗ ██╗  ██╗████████╗███████╗
    ██║     ██║██╔════╝ ██╔════╝ ██╔════╝ ██║  ██║╚══██╔══╝██╔════╝
    ██║     ██║██║  ███╗██║  ███╗██║  ███╗███████║   ██║   ███████╗
    ██║     ██║██║   ██║██║   ██║██║   ██║██╔══██║   ██║   ╚════██║
    ███████╗██║╚██████╔╝╚██████╔╝╚██████╔╝██║  ██║   ██║   ███████║
    ╚══════╝╚═╝ ╚═════╝  ╚═════╝  ╚═════╝ ╚═╝  ╚═╝   ╚═╝   ╚══════╝®

    DEM simulation engine, released by
    DCS Computing Gmbh, Linz, Austria
    http://www.dcs-computing.com, office@dcs-computing.com

    LIGGGHTS® is part of CFDEM®project:
    http://www.liggghts.com | http://www.cfdem.com

    Core developer and main author:
    Christoph Kloss, christoph.kloss@dcs-computing.com

    LIGGGHTS® is open-source, distributed under the terms of the GNU Public
    License, version 2 or later. It is distributed in the hope that it will
    be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. You should have
    received a copy of the GNU General Public License along with LIGGGHTS®.
    If not, see http://www.gnu.org/licenses . See also top-level README
    and LICENSE files.

    LIGGGHTS® and CFDEM® are registered trade marks of DCS Computing GmbH,
    the producer of the LIGGGHTS® software and the CFDEM®coupling software
    See http://www.cfdem.com/terms-trademark-policy for details.

-------------------------------------------------------------------------
    Contributing author and copyright for this file:
    This file is from LAMMPS
    LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
    http://lammps.sandia.gov, Sandia National Laboratories
    Steve Plimpton, sjplimp@sandia.gov

    Copyright (2003) Sandia Corporation.  Under the terms of Contract
    DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
    certain rights in this software.  This software is distributed under
    the GNU General Public License.

-------------------------------------------------------------------------

    Copyright 2022      Jelena Macak (DCS Computing GmbH, Linz; TU Graz)
    Copyright 2022      DCS Computing GmbH, Linz
------------------------------------------------------------------------- */

#include <cmath>
#include <stdlib.h>
#include <string.h>
#include "fix_property_thermal.h"
#include "fix_property_atom.h"
#include "atom.h"
#include "update.h"
#include "error.h"
#include "force.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixPropertyThermal::FixPropertyThermal(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), 
  nvalues_conductivity(0), 
  nvalues_capacity(0), 
  fix_variable_thermal_conductivity_(0), 
  fix_variable_thermal_capacity_(0), 
  variable_conductivity_flag_(0), 
  variable_capacity_flag_(0)
{
  if (narg < 5) error->all(FLERR,"Illegal fix property thermal command");

  int iarg = 3; 
  int i = 1; 

  variable_conductivity_flag_ = variable_capacity_flag_ = 0; 

  while (iarg < narg) {
      if (strcmp(arg[iarg],"thermalConductivity") == 0) {
          variable_conductivity_flag_ = 1; 
          i = 1; 
          while ( (iarg+i < narg) && (strcmp(arg[iarg+i],"thermalCapacity") != 0) ) {
              i++; 
          }
          nvalues_conductivity = i-1; 
          defaultvalues_conductivity = new double[nvalues_conductivity];
          for (int j = 1; j < i; j++ ) {
              defaultvalues_conductivity[j-1] = force->numeric(FLERR,arg[iarg+j]);
          }
          iarg += i; 
      }
      else if (strcmp(arg[iarg],"thermalCapacity") == 0) {
          variable_capacity_flag_ = 1; 
          i = 1;
          while ( (iarg+i < narg) && (strcmp(arg[iarg+i],"thermalConductivity") != 0) ) {
              i++; 
          }
          nvalues_capacity = i-1;
          defaultvalues_capacity = new double[nvalues_capacity];
          for (int j = 1; j < i; j++ ) {
              defaultvalues_capacity[j-1] = force->numeric(FLERR,arg[iarg+j]);
          }
          iarg += i; 
      }
      else
      {
          error->fix_error(FLERR,this,"Expecting 'thermalConductivity' or 'thermalCapacity'"); 
      }
  }
}

/* ---------------------------------------------------------------------- */

FixPropertyThermal::~FixPropertyThermal()
{
}

/* ---------------------------------------------------------------------- */

int FixPropertyThermal::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixPropertyThermal::init()
{

  f_Temp = static_cast<FixPropertyAtom*>(modify->find_fix_property("Temp","property/atom","scalar",1,0,style));

  fix_variable_thermal_conductivity_ = static_cast<FixPropertyAtom*>(modify->find_fix_property("variableThermalConductivity","property/atom","scalar",0,0,this->style,false));
  if(!fix_variable_thermal_conductivity_ && variable_conductivity_flag_)
  {
    const char* fixarg[10];
    fixarg[0]="variableThermalConductivity";
    fixarg[1]="all";
    fixarg[2]="property/atom";
    fixarg[3]="variableThermalConductivity";
    fixarg[4]="scalar";
    fixarg[5]="no";
    fixarg[6]="yes";
    fixarg[7]="yes";
    fixarg[8]="0.";
    fix_variable_thermal_conductivity_ = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
  }

  fix_variable_thermal_capacity_ = static_cast<FixPropertyAtom*>(modify->find_fix_property("variableThermalCapacity","property/atom","scalar",0,0,this->style,false));
  if(!fix_variable_thermal_capacity_ && variable_capacity_flag_)
  {
    const char* fixarg[10];
    fixarg[0]="variableThermalCapacity";
    fixarg[1]="all";
    fixarg[2]="property/atom";
    fixarg[3]="variableThermalCapacity";
    fixarg[4]="scalar";
    fixarg[5]="no";
    fixarg[6]="yes";
    fixarg[7]="yes";
    fixarg[8]="0.";
    fix_variable_thermal_capacity_ = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
  }

}

/* ---------------------------------------------------------------------- */

void FixPropertyThermal::pre_force(int vflag)
{
   int nlocal = atom->nlocal;
   
   if (variable_conductivity_flag_) {
     for (int i = 0; i < nlocal; i++) { 
         fix_variable_thermal_conductivity_->vector_atom[i] = 
             calc_polynomial(i, nvalues_conductivity, defaultvalues_conductivity);
     }
   }

   if (variable_capacity_flag_) {
     for (int i = 0; i < nlocal; i++) {
         fix_variable_thermal_capacity_->vector_atom[i] =
             calc_polynomial(i, nvalues_capacity, defaultvalues_capacity);
     }
   }

}

double FixPropertyThermal::calc_polynomial(int i, int nvalues, double *defaultvalues)
{
   double Temp_p = f_Temp->vector_atom[i]; 
   double value = defaultvalues[nvalues-1]; 
   for (int n = nvalues-2; n >= 0; n--) {
       value *= Temp_p; 
       value += defaultvalues[n]; 
   }
   
   return value;  
}

/* ---------------------------------------------------------------------- */
