% tests Runge-Kutta--based solver against analytical solution 
addpath('../src'); 

% flow regimes: 
caseList = ...
{ ...
  "Friction", ...
  "Heating", ...
  "Cooling", ...
  "Nozzle", ...
  "Diffuser" ...
}; 

% boundary conditions:
Uin     = [   5,  10,  10, 150,  10];  % velocity [m/s]
pin     = [1000, 100, 100, 100, 100];  % pressure [Pa]
Tin     = [ 900, 300, 300, 800, 300];  % temperature[K]
phipin  = [ 0.5, 0.5, 0.5, 0,   0.4];  % inlet solid fraction [-]
phipout = [ 0.5, 0.5, 0.5, 0.5,   0];  % outlet solid fraction [-]
L_dp    = [ 500, 100, 100, 20,  100];  % pipe length

R = 287.1;            % ideal gas constant [J/kgK]
k = 1.4;              % isentropic exponent [-]
rhoin = pin./(R*Tin); % inlet density [kg/m^3]

% flow regime coefficients:  
% friction:
f     = [1,    0,     0, 0, 0]; 
% heating/cooling:
Q     = [0, 0.01, -0.01, 0, 0]; 
% linear variation in solid fraction: 
C0    = [0, 0, 0, 1-phipin(4), 1-phipin(5)]; 
C     = [0, 0, 0, (1-phipout(4)-C0(4))/L_dp(4), (1-phipout(5)-C0(5))/L_dp(5)]; 
Cnorm = [0, 0, 0, C(4)/C0(4), C(5)/C0(5)]; 

fields = {'p', 'U', 'T', 'rho'}; 
yLabelList = {'p/p_{out}', 'U/U_{out}', 'T/T_{out}', '\rho/\rho_{out}'}; 
regimeList = ...
{ ...
  "friction", ...
  "heat", ...
  "heat", ...
  "isentropic", ...
  "isentropic" ...
};

for j=1:numel(Uin)

   yin = [pin(j) Uin(j) Tin(j) rhoin(j)]; % inlet state 
   Ldp = linspace(0, L_dp(j), 100);       % distance from inlet
   Min = sqrt(calcM2(rhoin(j), Uin(j), pin(j), k, phipin(j))); % inlet Mach 

   % other inputs for numerical solution
   otherInputs.k = k;
   otherInputs.f = f(j);
   otherInputs.Q = Q(j);
   otherInputs.C = C(j);
   otherInputs.C0 = C0(j);
   otherInputs.Cnorm = Cnorm(j); 
   otherInputs.phip = min([phipin(j), phipout(j)]);

   [Linlet, y_num] = ...
      ode45(@(x,y) flowDerivatives(x, y, otherInputs), Ldp, yin);

   % outlet state 
   pout = y_num(end,1); 
   Uout = y_num(end,2); 
   Tout = y_num(end,3); 
   rhoout = y_num(end,4); 
   % outlet Mach number
   Mout = sqrt(calcM2(rhoout, Uout, pout, k, phipout(j))); 

   % analytical solution
   y_sol = theoreticalRegimes(Min, linspace(Min, Mout, 100), k, regimeList{j}); 

   fig = figure(j);
   printf("\n%s\n\n", caseList{j})
   for i=1:4
     subplot(2,2,i)
     num = y_num(1,i)./y_num(:,i); % non-dimensional numerical sol
     % regime coefficient 
     X = otherInputs.f + otherInputs.Q + otherInputs.Cnorm; 
     plot ...
     ( ...
         Ldp(1:10:end), num(1:10:end), 'o', ...
         y_sol(1,:)/X, y_sol(i+1,:), '-' ... 
     )
     % relative error calculation
     err = abs(num'-y_sol(i+1,:))./y_sol(i+1,:); 
     maxerr = max(err); 
     meanerr = mean(err); 
     printf...
     (...
         "Field %s: max err %f%%, mean err %f%% \n", ...
          fields{i}, 100*maxerr, 100*meanerr ...
      )
     xlabel('x/dp')
     ylabel(yLabelList{i})
     if i==1
        title(caseList{j})
        legend('numerical', 'analytic', 'location', 'northwest')
        legend boxoff
     endif
   endfor
   nameFig = strcat(caseList{j}, ".png"); 
   saveas(fig, nameFig); 
endfor
