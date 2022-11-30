function [res] = theoreticalRegimes(Min, Mout, k, regime)
% OUTPUT: 
%  dimensionless state vector normalized by the outlet state
%  res(1) = coeff*L/dp ... dimensionless length
%  res(2) = p/pout ... pressure
%  res(3) = U/Uout ... velocity
%  res(4) = T/Tout ... temperature
%  res(5) = rho/rhoout ... density

% INPUT: 
%  Min - inlet Mach number
%  Mout - outlet Mach number
%  k - isentropic exponent
%  regime - 'friction', or 'heat', or 'isentropic'

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

switch ( regime )

   case { "friction" } 
   % constant friction coefficient
   % no heat exchange, constant solid fraction

     TinTout = calckM2plus1(Mout,k) ./ calckM2plus1(Min,k);

     pinpout = (Mout./Min) .* sqrt(TinTout); 

     UinUout = (Min./Mout) .* sqrt(TinTout); 

     % dimless length: fL/dp
     coeffLdp = ( Mout.^2 - Min.^2 ) ... 
          ./ ( 2 .* k * Min.^2 .* Mout.^2) ... 
          + ( (k+1) ./ (4.*k) ) ...
          .* log( TinTout .*  Min.^2 ./ Mout.^2 );

   case { "heat" }
   % constant heating/cooling coefficient
   % frictionless, constant solid fraction

     pinpout = ( k .* Mout.^2 + 1 ) ./ ( k .* Min.^2 + 1 ); 

     UinUout = ( Min.^2 ./ Mout.^2 ) .* pinpout; 

     TinTout = ( Min.^2 ./ Mout.^2 ) .* pinpout.^2; 

     % dimless length: QL/dp
     coeffLdp = log( ( Mout.^2 ./ Min.^2 ) .* pinpout.^( (k+1)./k ) ); 
      
   case { "isentropic" }
   % linear change in solid fraction: "nozzle" or "diffuser" configuration
   % frictionless, no heat exchange

     TinTout = calckM2plus1(Mout,k) ./ calckM2plus1(Min,k);

     pinpout = TinTout .^ ( k./(k-1) ); 

     UinUout = ( 1 ./ TinTout ) .^ ( 1./(k-1) );

     % dimless length: CLdp = Aout/Ain - 1
     coeffLdp = (Min./Mout) .* sqrt( TinTout.^( (k+1)./(k-1) ) ) - 1; 

   otherwise

     error...
     ( ...
          "Unfamiliar regime.\n \
          Please choose: 'friction', 'heat', or 'isentropic' " ... 
     );

endswitch

     rhoinrhoout = 1 ./ UinUout;

     res = [coeffLdp; pinpout; UinUout; TinTout; rhoinrhoout]; 

endfunction
