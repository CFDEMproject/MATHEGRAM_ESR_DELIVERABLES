# MATHEGRAM_ESR_DELIVERABLES

Collection of source code and data deliverables of MATHEGRAM ESRs

## Description

Numerical models for investigation of the fundamental mechanisms of heat
generation and transfer in granular materials

### Heat generation due to chemical reaction

Chemical reactions, by generating or consuming heat, affect the thermomechanical
behaviour of granular materials. A DEM model incorporating a solid state model
for a single reactant is developed. Reaction progression is modelled through
alpha, extent of reaction.  

- [Reaction example case](https://github.com/CFDEMproject/MATHEGRAM_ESR_DELIVERABLES/tree/main/LIGGGHTS-MATHEGRAM/examples/LIGGGHTS/Tutorials_public/reaction)

### Heat generation due to biological activity

Heat generation in biological material is mainly triggered by moisture content.
A DEM model for swelling is developed. 

- [Swelling example case](https://github.com/CFDEMproject/MATHEGRAM_ESR_DELIVERABLES/tree/main/LIGGGHTS-MATHEGRAM/examples/LIGGGHTS/Tutorials_public/swelling)

### Thermal radiation model

At high temperatures, radiation becomes a dominant mode of heat transfer. A
CFD-DEM model using P1 approximation is developed. 

Examples: 
- [Single particle cooling](https://github.com/CFDEMproject/MATHEGRAM_ESR_DELIVERABLES/CFDEMcoupling-MATHEGRAM/tutorials/cfdemSolverBuoyantPimple/MajorTomToGroundControl)
- [Imposed radiative heat
  flux](https://github.com/CFDEMproject/MATHEGRAM_ESR_DELIVERABLES/CFDEMcoupling-MATHEGRAM/tutorials/cfdemSolverBuoyantPimple/radiativeHeatFlux)
- [Small pebble bed under
  vacuum (only radiation)](https://github.com/CFDEMproject/MATHEGRAM_ESR_DELIVERABLES/CFDEMcoupling-MATHEGRAM/tutorials/cfdemSolverBuoyantPimple/pebbleBed_radiation)

### Temperature-dependant particle properties

Particle properties can change significantly in relation to temperature. A DEM
model for temperature-dependant heat capacity, and temperature-dependant
thermal conductivity is implemented. 

- [Small pebble bed under vacuum (radiation and conduction example case)](https://github.com/CFDEMproject/MATHEGRAM_ESR_DELIVERABLES/tree/main/CFDEMcoupling-MATHEGRAM/tutorials/cfdemSolverBuoyantPimple/pebbleBed_full)

### Non-constant density model

Strong heat exchange affects the fluid density. A CFD-DEM model for compressible
fluid flow is developed. 

Examples: 
- [Free
  convection](https://github.com/CFDEMproject/MATHEGRAM_ESR_DELIVERABLES/tree/main/CFDEMcoupling-MATHEGRAM/tutorials/cfdemSolverBuoyantPimple/freeConvection)
- [Friction
  regime](https://github.com/CFDEMproject/MATHEGRAM_ESR_DELIVERABLES/tree/main/CFDEMcoupling-MATHEGRAM/tutorials/cfdemSolverBuoyantPimple/subsonic_frictionRegime)
- [Heat
  regime](https://github.com/CFDEMproject/MATHEGRAM_ESR_DELIVERABLES/tree/main/CFDEMcoupling-MATHEGRAM/tutorials/cfdemSolverBuoyantPimple/subsonic_heatedRegime)
- [Nozzle
  flow](https://github.com/CFDEMproject/MATHEGRAM_ESR_DELIVERABLES/tree/main/CFDEMcoupling-MATHEGRAM/tutorials/cfdemSolverBuoyantPimple/subsonic_nozzleFlow)

One-dimensional compressible flows can be assessed using a GNU Octave 
[verification
tool](https://github.com/CFDEMproject/MATHEGRAM_ESR_DELIVERABLES/tree/main/CFDEMcoupling-MATHEGRAM/tutorials/cfdemSolverBuoyantPimple/subsonic_functions).

### Heat generation due to deformation

A FVM [contact deformation model](https://github.com/dllrun/foam-extend-general-contact) 
is developed for heat generation due to contact between solid particles. 

## Installation

Requires [LIGGGHTS-PUBLIC
3.8.0](https://github.com/CFDEMproject/LIGGGHTS-PUBLIC) and
[CFDEMcoupling-PUBLIC 3.8.0](https://github.com/CFDEMproject/CFDEMcoupling-PUBLIC). 
Details on installation are given on the "www.cfdem.com" website.

## Credit

| Feature                                    | Author                 |
| ------------------------------------------ | ---------------------- |
| FVM deformation model                      | Ranjan Dhakal (ESR2)   |
| DEM chemical reaction model                | Aman Rastogi (ESR3)    |
| DEM swelling model                         | Domenica Braile (ESR4) |
| DEM model with temperature dep. properties | Jelena Mačak (ESR5)    |
| CFD-DEM thermal radiation model            | Jelena Mačak (ESR5)    |
| CFD-DEM model with non-constant density    | Jelena Mačak (ESR5)    |

## Acknowledgment

Marie SKŁODOWSKA-CURIE Innovative Training Network **MATHEGRAM**, the People
Programme (Marie SKŁODOWSKA-CURIE Actions) of the European Union's Horizon
2020 Programme H2020 under REA grant agreement No.813202.
