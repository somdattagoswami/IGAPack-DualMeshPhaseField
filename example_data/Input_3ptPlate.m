function [geometry, Mater, Fract, Integ] = Input_3ptPlate
% Input File for Rectangular Plate
% use C1 cubics with adaptivity

geometry = struct;
geometry.dim = 2;
geometry.nstress = 3; % Number of stress components.
geometry.patchBoundaries = {1,2,3,1; 2,3,2,4; 3,4,1,3};

% Considering linear B-Splines in u and v directions
geometry.p = 3; % Degree of the polynomial in U Direction
geometry.q = 3; % Degree of the polynomial in V Direction
geometry.numPatches = 4; % Total number of patches
geometry.L = 8; % Length of the plate in the U direction
geometry.W = 2; % Height of the plate in the V direction
geometry.notchHt = 0.5;
geometry.toler = 1e-5;
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
Fract.cenerg = 0.5; % Critical energy release for unstable crack (Gc)
Fract.constl = 0.03; % L0 : Length parameter which controls the spread of the damage

% Time Integration Paramenters
Integ = struct;
Integ.nstep = 10000; % Number of Displacement steps
Integ.tfacto = 0; % Total increment factor for displacement increment
Integ.dfacto1 = -1e-4; % Displacement increment per time steps upto numStepsLimit
Integ.dfacto2 = -1e-5; % Displacement increment per time steps
Integ.numStepsLimit = 360;
Integ.nprint = 1; % Printing the results after how many steps

end