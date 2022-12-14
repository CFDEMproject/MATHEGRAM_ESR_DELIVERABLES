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

    Implementation of temperature-sensitive thermal conductivity
    Copyright 2022      Jelena Macak (DCS Computing GmbH, Linz; TU Graz)
    Copyright 2022      DCS Computing GmbH, Linz
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(heat/gran/conduction,FixHeatGranCond)
FixStyle(heat/gran,FixHeatGranCond)

#else

#ifndef LMP_FIX_HEATGRAN_CONDUCTION_H
#define LMP_FIX_HEATGRAN_CONDUCTION_H

#include "fix_heat_gran.h"

namespace LAMMPS_NS {

  class FixHeatGranCond : public FixHeatGran {
  public:
    FixHeatGranCond(class LAMMPS *, int, char **);
    ~FixHeatGranCond();
    virtual void post_create();
    virtual void pre_delete(bool);

    int setmask();
    void init();
    virtual void pre_force(int vflag);
    virtual void post_force(int vflag);

    virtual void cpl_evaluate(class ComputePairGranLocal *);
    void register_compute_pair_local(ComputePairGranLocal *);
    void unregister_compute_pair_local(ComputePairGranLocal *);

    virtual void updatePtrs();

  protected:
    int iarg_;

    template <int,int> void post_force_eval(int,int);

    class FixPropertyGlobal* fix_conductivity_;
    double *conductivity_;

    bool store_contact_data_;
    class FixPropertyAtom* fix_conduction_contact_area_;
    class FixPropertyAtom* fix_n_conduction_contacts_;
    class FixPropertyAtom* fix_wall_heattransfer_coeff_;
    class FixPropertyAtom* fix_wall_temperature_;
    class FixPropertyAtom* fix_variable_thermal_conductivity_; 
    double *conduction_contact_area_;
    double *n_conduction_contacts_;
    double *wall_heattransfer_coeff_;
    double *wall_temp_;

    // model for contact area calculation
    int area_calculation_mode_;

    double fixed_contact_area_;

    // for heat transfer area correction
    int area_correction_flag_;
    double const* const* deltan_ratio_;

    // for temperature dependent conductivity
    int variable_conductivity_flag_; 
  };

}

#endif
#endif

