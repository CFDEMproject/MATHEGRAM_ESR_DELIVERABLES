inputData = fopen("L50H20pc100.data", "w");
fdisp(inputData, "LIGGGHTS data file \n");

% user-defined properties
lenPerParticle = 50;           % bed length per particle
heightPerParticle = 20;        % bed height in particle diameters
percent = 1.0;                 % percentage of outlet cross section fullness

% fixed properties
dp = 1e-3;                     % particle diameter

% calculation of bed dimensions
atoms_in_x = 6; 
atoms_in_y = round(percent*heightPerParticle); 
inclination = atoms_in_y/lenPerParticle; 
atoms_in_z = lenPerParticle+atoms_in_x; 
emptyLength = atoms_in_x*dp; 


%distributionFigure = figure(1); 
atomNumber = 0; 
for iz = 1:atoms_in_z
   if round(inclination*iz)<=0
      localHeight = 1; 
   elseif round(inclination*iz)<=atoms_in_y
      localHeight = round(inclination*iz);
   else
      localHeight = atoms_in_y;  
   endif
   for iy = 1:localHeight
      for ix = 1:atoms_in_x
      atomNumber += 1; 
      data.no(atomNumber) = atomNumber; 
      data.x(atomNumber) = 0.5*dp + (ix-1)*dp; 
      data.y(atomNumber) = 0.5*dp + (iy-1)*dp; 
      data.z(atomNumber) = emptyLength+0.5*dp+(iz-1)*dp;  
%      plot(data.z(atomNumber), data.y(atomNumber), 'mo', 'markerfacecolor', 'm')
%      hold on 
      endfor
   endfor
endfor
%saveas(distributionFigure, "distributionFigure.png")

fdisp(inputData, strcat(num2str(data.no(end)), " atoms\n"));
fdisp(inputData, "1 atom types\n");

xDim = atoms_in_x*dp; 
yDim = atoms_in_y*dp; 
zDim = emptyLength+atoms_in_z*dp; 
xEndStr = num2str(xDim);
yEndStr = num2str(yDim); 
zStartStr = num2str(emptyLength); 
zEndStr = num2str(emptyLength+zDim);

x = ["0.0   " xEndStr " xlo " "xhi"];
y = ["0.0   " yEndStr " ylo " "yhi"];
z = [zStartStr " " zEndStr " zlo " "zhi"];

fdisp(inputData, x);
fdisp(inputData, y);
fdisp(inputData, z);

fdisp(inputData, "\n");
fdisp(inputData, "Atoms\n");

rho=2500;  % particle density
rhoStr=num2str(rho);
dpStr=num2str(dp);

for i=1:data.no(end)
   string = ...
   [ ...
     num2str(data.no(i)), ...
      "\t1\t ", ...
      num2str(dp), ... 
      "\t", ...
      num2str(rho), ...
      "\t", ...
      num2str(data.x(i)), ...
      "\t", ...
      num2str(data.y(i)), ...
      "\t", ...
      num2str(data.z(i)),
   ];
   fdisp(inputData, string)
endfor


