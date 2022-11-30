function res = calckM2plus1( M, k )
% OUTPUT:
%  calcultes frequently used relation:
%  (k-1)/2 * M^2 + 1. 
% INPUT:
%  M - Mach number,
%  k - the isentropic exponent. 

    res =   0.5 .* (k-1) .* M.^2 +  1;

endfunction
