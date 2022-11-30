%    P1 simulation results 
% vs view-factor results (Tausendschoen and Radl, Int. J. Heat Mass Transfer, 2021)
% vs analtical solution (Heaslet and Warming, Int. J. Heat Mass Transfer, 1965)

T1 = 1073.5;        % hot wall temperature, [K]
T2 = 436.5;         % cold wall temperature, [K]
L = 0.42;           % domain length, [m]
dp = 0.06;          % particle diameter, [m]
epsilon = 0.8;      % graphite emissivity, [-]
Qabs = 0.8;         % absorption efficiency = graphite emissivity, [-]
Qext = 1;           % extinction efficiency for large particles, [-]
Qsca = Qext-Qabs;   % scattering efficiency, [-]
Np = 332;           % number of particles
sigma = 5.67e-8;    % Stefan-Boltzmann constant, [W/(m^2K^4)]

Vp = dp^3*pi/6;             % particle volume, [m^3]
V = L^3;                    % enclosure volume, [m^3]
phip = Np*Vp/V;             % solid fraction, [-]
A = Np*dp^2*pi/4;           % projected area, [m^2]
betap = Qext*A/(1-phip)/V;  % extinction coefficient, [1/m]
tau = betap*L;              % opacity, [-]


% exact solution (Heaslet and Warming, 1965)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
psiB = (4/3)/(1.42089+tau);        % dimless heat flux, black walls 
psi = psiB/(1+psiB*(2/epsilon-2)); % dimless heat flux, gray wall

% calc dimless temp. profile
tauL = linspace(0, tau, 10);     
PHI = (-(3/4)*(tauL-tau) + (1/epsilon-1/2)) ./ (2/epsilon-1+(3/4)*tau);
% calc temp. profile
T_exact = (PHI*(T1^4-T2^4)+T2^4).^(1/4);
L_exact = tauL/tau*L; 

% view-factor solution (Tausendschoen and Radl, 2021)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
dataVF = load('./data/vf_radiation.csv'); 
L_VF = dataVF(:,1)*dp;          % dimensionalise length
T_VF = dataVF(:,2)*(T1-T2)+T2;  % dimensionalise temperature

% P1-results 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
path = '../../DEM/post/dump270000.liggghts_run'; 
dataP1 = load(path); 
L_P1 = dataP1(:,1); 
T_P1 = dataP1(:,4); 

% plotting
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

fig = figure(1); 
plot ...
( ...
   L_P1, T_P1, 'o', 'markersize', 3, 'color', [0.64, 0.64, 0.82] , ...
   L_VF, T_VF, 'x', 'markersize', 4, 'linewidth', 2, 'color', 'b', ...
   L_exact, T_exact, '--', 'color', 'b' ...
)
xlabel('length, [m]')
ylabel('temperature, [K]')
legend('P1', 'view-factor', 'analytic')
saveas(fig,'cfdemSolverBuoyantPimple_pebbleBedRadiation.png')

