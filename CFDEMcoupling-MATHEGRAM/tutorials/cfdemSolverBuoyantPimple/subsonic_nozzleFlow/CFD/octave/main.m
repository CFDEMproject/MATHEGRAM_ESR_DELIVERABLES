% this is a script for plotting comparison between analytic and numerical
% solution for isentropic flow with variable solid fraction

% FUNCTIONS % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
addpath('../../../subsonic_functions/src/');
addpath('../../../subsonic_functions/dataManip/');
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% VARIABLES % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% constant variables
dp = 1e-3;  % particle diameter   [m]
R = 287.1;  % ideal gas constant  [J/(kg.K)]
k = 1.4;    % isentropic exponent [-]

% path to simulation data
path = "../postProcessing/singleGraph/0.0025/line_mag(U)_p_rho_T_voidfraction.xy";

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% SIMULATION DATA % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
simData = makeSimulationData(path);
simData.L_dp = simData.L/dp; % dimensionless bed length [-]
simData.M = ...
   sqrt ...
   ( ...
      calcM2 ...
      ( ...
          simData.rho, ...
          simData.U, ...
          simData.p, ...
          k, ...
          simData.phip ...
      )...
    );


% NUMERICAL SOLUTION  % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
otherInputs.k = k;
otherInputs.Q = 0;

% Inlet state
p_init = simData.p(1);
U_init = simData.U(1);
T_init = simData.T(1);
rho_init = simData.rho(1);

phipin = simData.phip(1); 
phipout = simData.phip(end); 

% coefficients which describe linear variation in solid fraction
otherInputs.C0 = 1-phipin; 
otherInputs.C = (1-phipout-otherInputs.C0)/simData.L_dp(end); 
otherInputs.Cnorm = otherInputs.C / otherInputs.C0; 

otherInputs.phip = 0.; 

phip = calcPhipDistribution(otherInputs.C, otherInputs.C0, simData.L_dp);

% calculate local friction factor
local_f = calcFLocal(simData);
% force same-length vector
local_f = [local_f; local_f(end)]; 

fields = {'p', 'U', 'T', 'rho'};
numberOfElements = length(simData.p);

% inlet state
y_init = [p_init, U_init, T_init, rho_init];
% dimensionless bed length
xSpan = simData.L/dp; 

for j=1:numel(fields)
   numData.(fields{j})(1,1) = y_init(j);
endfor
numData.L_dp = xSpan;

for i=1:(numberOfElements-1)

   otherInputs.f = local_f(i); 
   x_init = [simData.L_dp(i), simData.L_dp(i+1)]; 
   [x_num, y_num] = ...
       ode45(@(x,y) flowDerivatives(x, y, otherInputs), x_init, y_init);
   y_init = y_num(end, :); 
   numData.L_dp(i+1) = x_num(end);

   for j = 1:numel(fields)
         numData.(fields{j})(i+1,1) = y_num(end, j) ;
   endfor

endfor

numData.M = ...
  sqrt ...
  ( ...
     calcM2 ...
     ( ...
         numData.rho, ...
         numData.U, ...
         numData.p, ...
         k, ... 
         phip ...
      ) ...
   );

dimlessSimData = makeDimlessData(simData);
dimlessNumData = makeDimlessData(numData);

% ERRORS  % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% CFD-DEM vs RK
simNumErr = calcError(simData, numData);
avg_simNumErr = avgError(simNumErr);

printf("\nAvg Err: %f%%\n", avg_simNumErr*100)

% PLOTTING  % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
fig = figure(1);
fields = {'p', 'U', 'T','rho'};
fname = {'p', 'U', 'T', '\rho'};
for j = 1:numel(fields)

    subplot(2,2,j)
    plot...
    ( ...
       simData.L_dp, dimlessSimData.(fields{j}), ...
       'r-', 'linewidth', 2, ...
       numData.L_dp, dimlessNumData.(fields{j}), ...
       'k--', 'linewidth', 2 ...
    )
    xlabel('L/d_p [-]', 'interpreter', 'latex')
    ylabel ...
    ( ...
        strcat( fname{j}, '/', fname{j}, '_{out} [-]' ), ...
        'interpreter', 'tex' ...
    )
    if j==1
        leg = legend('CFD-DEM', 'RK-4');
        set(leg, 'location', 'southwest', 'orientation', 'horizontal')
        legend boxoff
    endif
endfor

saveas(fig, 'cfdemSolverRhoPimple_nozzleFlow.png')

