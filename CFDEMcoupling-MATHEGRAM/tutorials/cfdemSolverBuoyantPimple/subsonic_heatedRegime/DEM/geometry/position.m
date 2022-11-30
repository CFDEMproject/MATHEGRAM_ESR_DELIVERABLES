inputData = fopen("L20.data", "w");
fdisp(inputData, "LIGGGHTS data file \n");

lenPerParticle = 20;
dp = 1e-3;                    % particle diameter
numPerCell = 2;               % number particles per cfd cell
Nx = 3*numPerCell;            % number of particles in x-direction
Ny = Nx;                      % number of particles in y-direction
Nz = 6/Nx;                    % number of particles in one layer in z-direction

numOfAtoms=Nx*Ny*lenPerParticle/Nz; % total number of particles 
fdisp(inputData, strcat(num2str(numOfAtoms), " atoms\n"));
fdisp(inputData, "1 atom types\n");

side = 6*dp;   
emptyLength = 6*dp; 
bedLength = lenPerParticle*dp;

% packed bed boundaries
sideStr = num2str(side);
startStr = num2str(emptyLength);
endStr = num2str(emptyLength+bedLength);

x = ["0.0   " sideStr " xlo " "xhi"];
y = ["0.0   " sideStr " ylo " "yhi"];
z = [startStr " " endStr " zlo " "zhi"];

fdisp(inputData, x);
fdisp(inputData, y);
fdisp(inputData, z);
fdisp(inputData, "\n");
fdisp(inputData, "Atoms\n");

rho=1000;  % particle density
rhoStr=num2str(rho);
dpStr=num2str(dp);

xstep = side/Nx; 
xLen = xstep/2 : xstep : side-xstep/2; 
ystep = side/Ny; 
yLen = ystep/2 : ystep : side-ystep/2;
zstep= bedLength/(lenPerParticle/Nz);
zLen = emptyLength + zstep/2 : zstep : emptyLength + bedLength - zstep/2; 

num=1;
for iz = 1:(lenPerParticle/Nz)
    for iy=1:Ny
        for  ix=1:Nx
            x = num2str(xLen(ix)); 
    	    y = num2str(xLen(iy));
	    z = num2str(zLen(iz));
            Num = num2str(num);
	    num = num+1;
            string=[Num " 1 " dpStr " " rhoStr " " x " " y " " z];
	    fdisp(inputData,string);
        end
    end
end
