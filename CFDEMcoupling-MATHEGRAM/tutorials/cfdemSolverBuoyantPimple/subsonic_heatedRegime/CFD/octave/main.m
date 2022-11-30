% compares the CFD-DEM simulation output 
% to analytic solution and 4th order Runge-Kutta solution

% FUNCTIONS % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
addpath('../../../subsonic_functions/src/');
addpath('../../../subsonic_functions/dataManip/');
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% VARIABLES % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% constant variables
dp = 1e-3;  % particle diameter   [m]
R = 287.1;  % ideal gas constant  [J/(kg.K)]
k = 1.4;    % isentropic exponent [-]
mu = 1e-5;  % dynamic viscosity   [Pa.s]
Pr = 0.72;  % Prandtl number      [-]
qp = 0.001;     % heat source to particles [W/m^3]
lambda = 0.014; % fluid conducitivity [W/(m.K)]
cp = 1007;      % fluid heat capacity [J/(kg.K)]

% path to simulation data
path = "../postProcessing/singleGraph/0.2/line_mag(U)_p_rho_T_voidfraction_Qsource.xy";

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

% NUMERICAL INTEGRATION % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

otherInputs.k = k; 
otherInputs.phip = mean(simData.phip); 
otherInputs.C = 0;
otherInputs.C0 = 0;
otherInputs.Cnorm = 0;

% Outlet state
pIn = simData.p(1); 
UIn = simData.U(1); 
TIn = simData.T(1); 
rhoIn = simData.rho(1); 

% calculate local friction factor
local_f = calcFLocal(simData);
% calculate heat transfer coeff
QSim = simData.Qs .* dp ./ (simData.rho .* simData.U .* simData.T .* cp);

% inlet state
y0 = [pIn UIn TIn rhoIn]; 
% dimensionless bed length
xSpan = simData.L/dp;

fields = {'p', 'U', 'T', 'rho'};
numberElements = length(simData.p);

for j=1:numel(fields)
   numData.(fields{j})(1,1) = y0(j);
endfor
numData.L_dp = xSpan; 


for i=1:(numberElements-1)

    otherInputs.Q = QSim(i); 
    otherInputs.f = local_f(i); 
    x0 = [xSpan(i) xSpan(i+1)]; 
    [x_num, y_num] = ...
        ode45(@(x,y) flowDerivatives(x, y, otherInputs), x0, y0);

    y0 = y_num(end,:);
    for j = 1:numel(fields)
          numData.(fields{j})(i+1,1) = y0(j) ;
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
         otherInputs.phip ...
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
fields = {'p', 'U', 'T','rho','M'};
fname = {'p', 'U', 'T', '\\rho', 'M'};
for j = 1:numel(fields)
    subplot(3,2,j)
    plot ...
    ( ...
        simData.L_dp, dimlessSimData.(fields{j}), ...
        'r-', 'linewidth', 2, ...
        numData.L_dp, dimlessNumData.(fields{j}), ...
        'k--', 'linewidth', 2 ...
    )
    xlabel('L/d_p [-]', 'interpreter', 'latex')
    ylabel...
    ( ...
       strcat( fname{j}, '/', fname{j}, '_{out} [-]'), ...
       'interpreter', 'latex' ...
    )
    if j==1
        leg = legend('CFD-DEM', 'RK-4');
        set(leg, 'location', 'north', 'orientation', 'horizontal')
        legend boxoff
    endif
endfor
subplot(3,2,6)
plot(simData.L_dp, QSim, 'r-', 'linewidth', 2)
xlabel('L/d_p [-]', 'interpreter', 'latex')
ylabel('Q [-]', 'interpreter', 'latex')
saveas(fig, "cfdemSolverRhoPimple_heatedFlow.png")
