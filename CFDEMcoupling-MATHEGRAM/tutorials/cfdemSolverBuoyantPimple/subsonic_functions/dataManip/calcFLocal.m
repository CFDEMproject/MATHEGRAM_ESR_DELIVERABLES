function f = calcFLocal(simData)
% OUTPUT:  
%  calculates local friction factor based on the
%  local simulation values
% INPUT:
%  simData - structure contaning simulation data

   delta_p = simData.p(1:end-1) - simData.p(2:end); 
   delta_L = simData.L_dp(2:end) - simData.L_dp(1:end-1); 

   mean_U = 0.5*(simData.U(1:end-1)+simData.U(2:end)); 
   mean_rho = 0.5*(simData.rho(1:end-1)+simData.rho(2:end));
   mean_phip = 0.5*(simData.phip(1:end-1)+simData.phip(2:end)); 

   f = delta_p ...
       .* (1-mean_phip).^2 ...
       ./ ( ...
               mean_U.^2 ...
            .* mean_rho ...
            .* delta_L ...
          );
           
endfunction
