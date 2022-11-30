function res = makeDimlessData(data)
% INPUT: 
%  data - stucture with fields p, U, T, rho, M 
% OUTPUT: 
%  res - structure containing dimensionless fields: 
%    res.p    - pressure/outlet pressure             [Pa]
%    res.U    - superficial velocity/outlet veocity  [m/s]
%    res.T    - temperature/outlet temperature          [K]
%    res.rho  - density/outlet density              [kg/m^3]
%    res.M    - Mach number/outlet Mach number      [-]

   fieldNames = {'p', 'U', 'T', 'rho', 'M'}; 
   
   for i = 1:numel(fieldNames)
       outletVal = data.(fieldNames{i})(end); 
       res.(fieldNames{i}) = data.(fieldNames{i})/outletVal; 
   endfor  

endfunction
