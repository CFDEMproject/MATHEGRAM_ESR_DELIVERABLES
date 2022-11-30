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
#ifdef FIX_CLASS

FixStyle(reaction,FixReaction)


#else


#ifndef LMP_FIX_REACTION_H
#define LMP_FIX_REACTION_H

#include "fix.h"
#include "fix_heat_gran.h"




namespace LAMMPS_NS {

  class FixReaction : public Fix {

    friend class FixMultisphere;
    friend class Multisphere;

  public:
    FixReaction(class LAMMPS *, int, char **);
    ~FixReaction(){};
    virtual void post_create();
    virtual void pre_delete(bool unfixflag){ UNUSED(unfixflag); };

    //void post_integrate(int vflag);
    //void end_of_step();
    void pre_final_integrate();

    virtual int setmask();
    virtual void init();

    // per default these three methods throw errors.
    //virtual void cpl_evaluate(class ComputePairGranLocal *);
    //virtual void register_compute_pair_local(class ComputePairGranLocal *);
    //virtual void unregister_compute_pair_local(class ComputePairGranLocal *);
    virtual void updatePtrs();

  protected:
    class FixPropertyAtom* fix_concentration;
    class FixPropertyAtom* fix_dconcentration;
    class FixPropertyAtom* fix_temp;
    class FixPropertyAtom* fix_heatSource;

    double *concentration;
    double *dconcentration;          
    double C0;          
    double *Temp; 
    double *heatSource;
    double pre_exponent;
    double activation_energy;
    double heat_reaction;
    double model;
    double exponent;
    int heat_reaction_flag;
    double max_conversion;
    
    class PairGran *pair_gran;
    int history_flag;
  };

}
#endif
#endif
