function res = avgError(errData)
% INPUT: 
%  errData - structures with errors in fields p, U, T, rho. 
% OUTPUT: 
%  res - average error

   err_p = mean(errData.p);
   err_U = mean(errData.U);
   err_T = mean(errData.T);
   err_rho = mean(errData.rho);

   res = sqrt(err_p.^2 + err_U.^2 + err_T.^2 + err_rho.^2);  

endfunction
