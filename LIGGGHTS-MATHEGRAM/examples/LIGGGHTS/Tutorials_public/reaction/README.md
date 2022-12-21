## LIGGGHTS_reaction

Solid state reaction models for a single reactant. Reaction progression is modelled through alpha, extent of reaction. More details can be found in [1]. 

Most of the reaction models in [1] fall into one of these categories:
  1. (alpha)^n - Type 1
  2. (1-alpha)^n - Type 2
  3. n(1-alpha)(-ln(1-alpha))^(1-1/n) - Type 3

**fix command**

|**fix id** group reaction initial_concentration C0 pre_exponent A activation_energy Ea maximum_conversion Cmax model type exponent keywords args|
|:---|

**keywords**
  1. _model_ #default=0
   - _type_ values = 0 or 1 or 2
   - exponent values should be positive      
  2. _heat_reaction_ #calculate heat reaction flag, default=no 
   - values = yes or no
   - if yes, then H_R (heat of reaction in J/g) has to be provided	

**Example:**

| **fix	react** all reaction initial_concentration 0.8 pre_exponent 100000000 activation_energy 150000 heat_reaction yes 4000000 maximum_conversion 0.5 model 3 2 |
|:-------------------|

**Notes:**
Model exponents need to be positive
heat/gran needs to be defined prior to this fix so that temperature and heatsource properties are created. This fix resets the heatSource at every timestep, hence cannot be paired with another heatSource setting.
All units are SI units.



**Output:**
Two scalars: f_concentration[0], f_dconcentration[0]


[1]	A. Khawam and D. R. Flanagan, “Solid-state kinetic models: basics and mathematical fundamentals,” The journal of physical chemistry B, vol. 110, no. 35, pp. 17315–17328, 2006.

