%    cfdemSolverBuoyantPimple results 
% vs experimental data (Corvaro and Paroncini, Appl. Therm. Eng., 2008)

Th = 297;                      % heating strip temperature, [K]
Tc = 285;                      % aluminum walls temperature, [K]
l = 0.05;                      % cavity side length, [m]
k = 0.025;                     % air thermal conductvity, [W/(mK)]

constant = l/(k*(Th-Tc)); 

% experimental data (Corvaro and Paroncini, 2008) 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
expData = load('./data/d04Ra171.csv'); 
x_exp = expData(:,1);  % length, [cm]
Nu_exp = expData(:,2); 
Nu_exp_avg = sum(expData(:,2))/length(expData(:,2))

% simulation data 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
path = "../postProcessing/localWallHeatFlux1/10/wallHeatFlux_heatExchanger.raw";
data = load(path);
x = data(:,1)*100;       % length, [cm]
Nu = data(:,4)*constant; % Nu number, [-]
Nu_sim_avg=sum(data(:,4))*constant/length(data(:,4))
% without edge values
Nu_sim_mid = sum(data(2:end-1,4)*constant)/length(data(2:end-1,4))


% plotting
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
fig = figure(1)
plot...
( ...
     x, Nu, 'rx-', ...
     x_exp, Nu_exp, 'bo' ...
)
xlabel('length, [cm]')
ylabel('Nu, [-]')
title("Local Nu for asymmetric configuration and Ra = 1.71 e5")
a = ...
   strcat ...
   (...
      'numerical, Nu_{avg} = ', ...
      num2str(Nu_sim_avg), ...
      ', w/o edges: ', ...
      num2str(Nu_sim_mid) ...
   );
b = ...
  strcat...
  ( ...
     'experimental, Nu_{avg} = ', ....
      num2str(Nu_exp_avg)
  );
leg = legend(a,b)
set(leg, 'Interpreter', 'tex')
saveas(fig, 'cfdemSolverBuoyantPimple_localNu.png')
