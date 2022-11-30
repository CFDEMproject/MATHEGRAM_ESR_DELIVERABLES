function res = calcError(data1, data2)
% INPUT: 
%  data1, data2 - structures containing fields p, U, T, rho.
% OUTPUT: 
%  res.p - error in pressure
%  res.U - error in velocity
%  res.T - error in temperature
%  res.rho - error in density 

   fieldNames = {'p', 'U', 'T', 'rho'}; 

   for i =  1:numel(fieldNames)
       res.(fieldNames{i}) = data1.(fieldNames{i})./data2.(fieldNames{i}) - 1; 
   endfor

endfunction
