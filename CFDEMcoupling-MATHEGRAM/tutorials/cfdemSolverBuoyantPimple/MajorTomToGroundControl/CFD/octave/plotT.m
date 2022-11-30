% m file
close all;
clear;
clc;

% input data
rho = 1000;        % particle density [kg/m^3]
cp = 1000;         % heat capacity [J/(kgK)]
dp = 0.05;         % particle diameter [m]
Qabs = 1;          % particle absorption efficiency [-]

Tmin = 500;        % minimal temperature [K]
Tmax = 1000;       % maximal temperature [K]


% simulation data
path = '../../DEM/tempFile';
data = load(path); 
Tsim =  data(:,2);
tsim = data(:,1);

% analytical solution

%% function
function time = calcSingleParticleTime(Tp, T0, otherInputs)
% Analytical estimation of the time needed for a particle with initial
% temperature Tp, to attain temperature of surroundings T0. 
% Goes together with 'Major Tom' test-case
% input: 
%  Tp - initial particle temperature
%  T0 - temperature of the enclosure
%  otherInputs.rho - particle density [kg/m^3]
%  otherInputs.c - particle heat capacity [J/(kgK)]
%  otherInputs.dp - particle diameter [m]
%  otherInputs.Qabs - absorption efficiency [-]

   sigma = 5.67e-8; % Stefan-Boltzmann constant
   tau =  otherInputs.rho.*otherInputs.c.*otherInputs.dp ...
       ./ (24*sigma*otherInputs.Qabs.*T0.^3);

   time = tau*(2*atan2(Tp,T0)-log(abs(Tp-T0)./(Tp+T0)));
endfunction
%% end function

otherInputs.rho = rho; 
otherInputs.c = cp; 
otherInputs.dp = dp; 
otherInputs.Qabs = Qabs; 

% calculate analytical solution
Tcooling = Tmax:-1:Tmin;
timeC = ...
   calcSingleParticleTime(Tcooling, Tmin, otherInputs) ...
 - calcSingleParticleTime(Tmax, Tmin, otherInputs);
Theating = Tmin:1:Tmax;
timeH = ...
    calcSingleParticleTime(Theating, Tmax, otherInputs) ...
  - calcSingleParticleTime(Tmin, Tmax, otherInputs);

% plot data
figure(1)
step = 2; 
plot ...
( ...
   timeC, Tcooling, 'b-', 'linewidth', 2, ...
   timeH, Theating, 'r-', 'linewidth', 2, ...
   tsim(1:step:end), Tsim(1:step:end), 'ko', 'markersize', 2, 'markerfacecolor', 'k' ...
)
legend('cooling', 'heating', 'simulation')
legend boxoff
xlabel("time, [s]")
ylabel("temperature, [K]")
title('simulation vs analytical')
print -color "cfdemSolverRhoPimple_Tsolid.png"

