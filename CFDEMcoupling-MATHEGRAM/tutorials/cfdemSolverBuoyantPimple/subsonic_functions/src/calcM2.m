function M2 = calcM2(rho, U, p, k, phip)
% OUTPUT:  
%  calculates squared Mach number M^2
% INPUT:
%  rho - density, 
%  U - superficial velocity, 
%  p - pressure, 
%  k - isentropic exponent, 
%  phip - solid fraction.

    M2 = rho .* U.^2 ./ (k .* (1-phip).^2 .* p);

endfunction
