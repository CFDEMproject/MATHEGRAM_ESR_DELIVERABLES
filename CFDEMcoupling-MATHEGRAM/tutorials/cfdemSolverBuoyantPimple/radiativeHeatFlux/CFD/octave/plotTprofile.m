%    P1 simulation results
%    vs analtical solution (Heaslet and Warming, Int. J. Heat Mass Transfer, 1965)

Thot4=1000^4;       % theoretical hot plate temperature^4 [K^4]
Tcool4=500^4;       % cold plate temperature^4 [K^4]
len = 0.02;         % distance between plates [m]
tau = 10;           % optical thickness [-]

% exact solution
dataLit = load('./data/tau10.csv'); 
xE = dataLit(:,1); 
psiE = dataLit(:,2); 

% simulation results
path = "../../DEM/post/dmpTemp10000"; 
data = load(path); 
x = data(:,3)/len; 
T = data(:,4); 
psi = (T.^4 - Tcool4)/(Thot4-Tcool4); 

fig = figure(1); 
step = 10; 
plot ...
( ...
  xE, psiE, 'k-', ...
  x(1:step:end), psi(1:step:end), ...
    'o', 'markersize', 3, 'color', [0.64, 0.64, 0.82] ... 
)
axis([0 1 0 1])
legend...
(...
   'exact', ...
   'simulation' ...
)
xlabel('distance x/L')
ylabel('normalized temperature')

saveas(fig, "cfdemSolverBuoyantPimple_temperatureDistribution.png")
