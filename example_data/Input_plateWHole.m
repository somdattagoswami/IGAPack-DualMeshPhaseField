function [geometry, Mater, Fract, Integ] = Input_plateWHole
% Input File for Rectangular Plate
% use C1 cubics with adaptivity

geometry = struct;
geometry.dim = 2;
geometry.nstress = 3; % Number of stress components.
geometry.patchBoundaries = {1,2,3,1; 2,3,3,1; 3,4,3,1; 4,8,2,4; [3,8],7,[2,1],[4,3];
    [2,7],6,[2,1],[4,3]; [1,6],5,[2,1],[4,3]; 8,10,2,3; 10,9,2,4; 9,12,2,2;
    [10,12],11,[4,4],[2,4]; 11,13,3,3; 13,16,4,2; 16,15,4,2; [7,13,15],14,[2,2,4],[3,2,4];
    15,17,3,3; 17,20,2,2; 20,19,4,2; [6,17,19],18,[2,4,4],[3,2,4]; 19,21,3,3; 12,22,3,4;
    [16,22],23,[3,1],[4,3]; [20,23],24,[3,1],[4,3]; [21,24],25,[2,1],[4,3]; 22,26,2,4;
    [23,26],27,[2,1],[4,3]; [24,27],28,[2,1],[4,3]; [25,28],29,[2,1],[4,3];
    26,30,2,4; [27,30],31,[2,1],[4,3]; [28,31],32,[2,1],[4,3]; [29,32],33,[2,1],[4,3]};

% Considering linear B-Splines in u and v directions
geometry.p = 3; % Degree of the polynomial in U Direction
geometry.q = 3; % Degree of the polynomial in V Direction
geometry.numPatches = 33; % Total number of patches
geometry.numElemU = 3; % Initial mesh
geometry.numElemV = 3;
geometry.toler = 1e-2;
geometry.ngaussX = geometry.p+1;
geometry.ngaussY = geometry.q+1;
geometry.maxERefLevel = 5; % Maximum level of Refinement for Elastic Mesh
geometry.maxPhRefLevel = 6; % Maximum level of Refinement for Phase Field Mesh
geometry.threshPhi = 0.5; % Threshold for Refinement 
geometry.B = 1e2;% Parameter for initial history function

% Material properties
Mater = struct;
Mater.E = 20.8*1e3; % Young's Modulus based on (N/mm2)
Mater.nu = 0.3; % Poisson's Ratio
% Mater.C = Mater.E/(1-Mater.nu^2)*[1, Mater.nu, 0; Mater.nu, 1, 0; 0, 0, (1-Mater.nu)/2]; % Plane Stress
Mater.C = (Mater.E/((1+Mater.nu)*(1-2*Mater.nu)))*[ 1-Mater.nu,Mater.nu,0;Mater.nu,1-Mater.nu,0;0,0,0.5-Mater.nu]; % Plane Strain
Mater.lamda = Mater.nu*Mater.E/((1+Mater.nu)*(1-2*Mater.nu)); % Lame Constant
Mater.mu = Mater.E/(2*(1+Mater.nu)); % Lame Constant

% Properties for Fracture
Fract = struct;
Fract.cenerg = 1; % Critical energy release for unstable crack (Gc)
Fract.constl = 0.025; % L0 : Length parameter which controls the spread of the damage

% Time Integration Paramenters
Integ = struct;
Integ.nstep = 2000; % Number of Displacement steps
Integ.tfacto = 0; % Total increment factor for displacement increment
Integ.dfacto1 = -1e-3; % Displacement increment per time steps upto numStepsLimit
Integ.dfacto2 = -1e-5; % Displacement increment per time steps
Integ.numStepsLimit = 192;
Integ.nprint = 1; % Printing the results after how many steps

end