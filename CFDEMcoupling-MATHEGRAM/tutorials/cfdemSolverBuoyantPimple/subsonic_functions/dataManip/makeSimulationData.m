function res = makeSimulationData(path)
% INPUT: 
%  path - path to data sampled via singleGraph openFoam utility
% OUTPUT: 
%  res - structure containing dimensioned fields: 
%    res.L    - length               [m]
%    res.p    - pressure             [Pa]
%    res.U    - superficial velocity [m/s]
%    res.T    - temperature          [K]
%    res.rho  - density              [kg/m^3]
%    res.phip - solid fraction       [-]

   data = load(path); 
   res.L = data(:,1)-data(1,1);  % length                [m]
   res.p = data(:,3);            % pressure              [Pa]
   res.T = data(:,5);            % temperture            [K]
   res.rho = data(:,4);          % density               [kg/m^3]
 
   I = ones(size(data(:,1)));    % unit vector
   phiv = data(:,6);             % void fraction         [-]
   
   res.phip = I - phiv;          % solid fraction        [-]
   
   res.U = data(:,2).*phiv;      % superificial velocity [m/s]

   if columns(data) >= 7
      res.Qs = data(:,7);
   endif

   if columns(data) >= 8
      res.Tp = data(:,8);
   endif


   
endfunction
