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
    (if not contributing author is listed, this file has been contributed
    by the core developer)

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz

    Copyright 2022      Aman Rastogi (University of Surrey)
    Copyright 2022-     DCS Computing GmbH, Linz
------------------------------------------------------------------------- */

#include "fix_reaction.h"

#include "atom.h"
#include "fix_property_atom.h"
#include "global_properties.h"
#include "fix_scalar_transport_equation.h"
#include "force.h"
#include "group.h"
#include "math_extra.h"
#include "modify.h"
#include "pair_gran.h"
#include <stdlib.h>
#include "update.h"
#include <cmath>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixReaction::FixReaction(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg){

  if ((!atom->radius_flag)||(!atom->rmass_flag)) error->all(FLERR,"Fix reaction needs per particle radius and mass");

  if (narg < 9)
    error->fix_error(FLERR,this,"not enough arguments");

  int iarg = 3;
  bool hasargs = true;
  while(iarg < narg && hasargs)
  {
  hasargs = false;
   if(strcmp(arg[iarg],"initial_concentration")==0){
     //error->fix_error(FLERR,this,"expecting keyword 'initial_concentration'");
     printf("processing initial conc");
     C0 = atof(arg[iarg+1]);
     printf("C0 is equal to %f\n", C0);
     iarg += 2;
     hasargs= true;
    }else if(strcmp(arg[iarg],"pre_exponent") == 0) {
      pre_exponent = atof(arg[iarg+1]);
      printf("pre_exponent is equal to %f\n", pre_exponent);
      iarg += 2;
      hasargs = true;
    }else if(strcmp(arg[iarg],"activation_energy") == 0) {
      
      activation_energy = atof(arg[iarg+1]);
      printf("Ea is equal to %f\n", activation_energy);
      iarg += 2;
      hasargs = true;
    }else if(strcmp(arg[iarg],"heat_reaction") == 0) {
      
      if (strcmp(arg[iarg+1],"yes") == 0){
      heat_reaction_flag = 1;  
      heat_reaction = -1*atof(arg[iarg+2]);
      printf("Heat of Reaction is equal to %f\n", heat_reaction);// heat of reaction in J per mol
      iarg += 3;
      hasargs = true;}
      else if (strcmp(arg[iarg+1],"no") == 0){
      heat_reaction_flag = 0;  
      iarg += 2;
      hasargs = true;
      }
    }else if(strcmp(arg[iarg],"maximum_conversion") == 0) {
      
      max_conversion = atof(arg[iarg+1]);
      printf("Maximum conversion is %f\n", max_conversion);
      iarg += 2;
      hasargs = true;
    }else if (strcmp(arg[iarg],"model") == 0){
      
      model = atof(arg[iarg+1]);
      printf("Model is type %f \n", model);
      exponent = atof(arg[iarg+2]);
      iarg += 3;
      hasargs = true;
      }
  
    
  }
  
  fix_concentration = NULL;
  fix_dconcentration = NULL;
  fix_temp = NULL;
  fix_heatSource = NULL;

  peratom_flag = 1;      
  size_peratom_cols = 0; 
  peratom_freq = 1;

  scalar_flag = 1; 
  global_freq = 1; 

  //cpl = NULL;
}

/* ---------------------------------------------------------------------- */

void FixReaction::post_create()
{
   // register concentration	
   if(!fix_concentration)
  {
    const char* fixarg[9];
    fixarg[0]="concentration";   // This part is similar to defining a fix statement from the input script. so this fix defines the variable concentration. 
    fixarg[1]="all";
    fixarg[2]="property/atom";
    fixarg[3]="concentration";
    fixarg[4]="scalar";
    fixarg[5]="yes";
    fixarg[6]="yes";
    fixarg[7]="no";    
    char arg8[30];
    sprintf(arg8,"%e",C0);
    fixarg[8]=arg8;
    
    fix_concentration= modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);  // This is the line that is adding the new atom/property.  
  }
  
    // register dconcentration	
   if(!fix_dconcentration)
  {
    const char* fixarg[9];
    fixarg[0]="dconcentration";   // This part is similar to defining a fix statement from the input script. so this fix defines the variable concentration. 
    fixarg[1]="all";
    fixarg[2]="property/atom";
    fixarg[3]="dconcentration";
    fixarg[4]="scalar";
    fixarg[5]="yes";
    fixarg[6]="yes";
    fixarg[7]="no";    
    fixarg[8]="0.";
    
    fix_dconcentration= modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);  // This is the line that is adding the new atom/property.  
  }
  
  fix_temp = static_cast<FixPropertyAtom*>(modify->find_fix_property("Temp","property/atom","scalar",0,0,style));
  fix_heatSource = static_cast<FixPropertyAtom*>(modify->find_fix_property("heatSource","property/atom","scalar",0,0,style));
  fix_concentration = static_cast<FixPropertyAtom*>(modify->find_fix_property("concentration","property/atom","scalar",0,0,style));
  fix_dconcentration = static_cast<FixPropertyAtom*>(modify->find_fix_property("dconcentration","property/atom","scalar",0,0,style));
  if(!fix_concentration || !fix_temp)
    error->one(FLERR,"internal error");
}

/* ---------------------------------------------------------------------- */

void FixReaction::updatePtrs()
{
  concentration = fix_concentration->vector_atom;
  dconcentration = fix_dconcentration->vector_atom;
  heatSource = fix_heatSource->vector_atom;
  //vector_atom = Temp; 
  Temp = fix_temp->vector_atom;

}

/* ---------------------------------------------------------------------- */

void FixReaction::init()
{
  
  if (!atom->radius_flag || !atom->rmass_flag)
    error->fix_error(FLERR,this,"must use a granular atom style ");

    // check if a fix of this style already exists
  //if(modify->n_fixes_style(style) > 1)
   // error->fix_error(FLERR,this,"cannot have more than one fix of this style");

  if(!force->pair_match("gran", 0))
    error->fix_error(FLERR,this,"needs a granular pair style to be used");

  pair_gran = static_cast<PairGran*>(force->pair_match("gran", 0));
  history_flag = pair_gran->is_history();
  //int max_type = atom->get_properties()->max_type();
  
  
  //fix_temp = static_cast<FixPropertyAtom*>(modify->find_fix_property("Temp","property/atom","scalar",0,0,style));
  //fix_heatSource = static_cast<FixPropertyAtom*>(modify->find_fix_property("heatSource","property/atom","scalar",0,0,style));
  //fix_concentration = static_cast<FixPropertyAtom*>(modify->find_fix_property("concentration","property/atom","scalar",0,0,style));
  
  //fix_dconcentration = static_cast<FixPropertyAtom*>(modify->find_fix_property("dconcentration","property/atom","scalar",0,0,style));
  
  if(!fix_concentration || !fix_temp)
    error->one(FLERR,"internal error");

  updatePtrs();
}

/* ---------------------------------------------------------------------- */

int FixReaction::setmask()
{
  int mask = 0;
  mask |= PRE_FINAL_INTEGRATE;
  //mask |= POST_INTEGRATE;
  //mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */
void FixReaction::pre_final_integrate()
{
  
  updatePtrs();
  
  //double *radius = atom->radius;
  double dt = update->dt;
  double *rmass = atom->rmass;
  double dalpha =0;
  double falpha =0;
  double modelcalc=0;
  int nlocal = atom->nlocal;
  for (int i=0; i<nlocal; i++)
   {
   
   //Concentration is mass of reactant/mass of particle 
   //C0 is initial concentration
   
   double alpha = (C0 - concentration[i])/(C0*(1-max_conversion)) ; //Extent of reaction
   
   if(model==1.0){
   //power (alpha)
   falpha = pow(alpha,exponent);
   }
   else if(model==2.0){
   //power (1-alpha)
   falpha = pow(1.0-alpha,exponent);
   }
   else if(model==3.0){
   //Avrami Erofeyev
   if (alpha == 0){
   
   alpha = 10e-12; //to avoid log(1)
   }
   falpha = exponent*(1-alpha)*pow((-1*log(1-alpha)),(1-(1/exponent)));
   }
   
   
   
   // dalpha = A*exp(E/RT)*f(alpha)*dt //Explicit integration for a single rate equation
   dconcentration[i] = -1*pre_exponent*exp(-1*activation_energy/(8.314*Temp[i]))*falpha*dt*(C0*(1-max_conversion));
   concentration[i] += dconcentration[i];
   
   
   if (heat_reaction_flag == 1){
   
   heatSource[i] = dconcentration[i]*rmass[i]*heat_reaction/(dt); // heat of reaction is in J/kg // Qdot = mass of reactant converted in time dt * heat of reaction/dt
   }
   }
   
  updatePtrs();
}



/* ---------------------------------------------------------------------- */
/*
void FixReaction::post_integrate(int vflag)
{
  
  updatePtrs();
  
  //double *radius = atom->radius;
  double dt = update->dt;
  
  
  int nlocal = atom->nlocal;
  for (int i=0; i<nlocal; i++)
   {
   
   dconcentration[i] = -1*pre_exponent*exp(activation_energy/(8.314*Temp[i]))*(pow(1-concentration[i],(1/3)))*dt;
   
   //char conc[30];
   //sprintf(conc,"%e",dconcentration[i]);
   //error->warningAll(FLERR, conc);
   //fprintf(universe->ulogfile, conc);
   concentration[i] += dconcentration[i];
   concentration[i] = 2.0;
   }
   
  updatePtrs();
  fix_concentration->do_forward_comm();
}
*/
/* ---------------------------------------------------------------------- */

