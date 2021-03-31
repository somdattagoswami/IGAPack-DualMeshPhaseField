function [geometry, Mater, Fract, Integ] = Input_cubeWHole
% Input File for Rectangular Plate
% use C1 cubics with adaptivity

geometry = struct;
geometry.dim = 3;
geometry.nstress = 6; % Number of stress components.
geometry.patchBoundaries = {1,2,2,4;1,3,3,1;[2,3],4,[3,2],[1,4]};

% Considering linear B-Splines in u and v directions
geometry.p = 3; % Degree of the polynomial in U Direction
geometry.q = 3; % Degree of the polynomial in V Direction
geometry.r = 3;
geometry.L = 1; % Length of the plate in the U direction
geometry.W = 1; % Width of the plate in the V direction
geometry.H = 1; % Height of the plate in W Direction
geometry.numPatches = 4; % Total number of patches
geometry.numElemU = 2; % Initial mesh
geometry.numElemV = 2;
geometry.numElemW = 2;
geometry.ngaussX = geometry.p+1;
geometry.ngaussY = geometry.q+1;
geometry.ngaussZ = geometry.r+1;
geometry.maxERefLevel = 3; % Maximum level of Refinement for Elastic Mesh
geometry.maxPhRefLevel = 4; % Maximum level of Refinement for Phase Field Mesh
geometry.threshPhi = 0.5; % Threshold for Refinement 
geometry.B = 50;% Parameter for initial history function
geometry.toler = 1e-2;

% Material properties
Mater = struct;
Mater.E = 20.8*1e3; % Young's Modulus based on (N/mm2)
Mater.nu = 0.3; % Poisson's Ratio
% Mater.C = Mater.E/(1-Mater.nu^2)*[1, Mater.nu, 0; Mater.nu, 1, 0; 0, 0, (1-Mater.nu)/2]; % Plane Stress
Mater.C = zeros(6,6);
Mater.C(1:3,1:3)= Mater.E/(1+Mater.nu)/(1-2*Mater.nu)*[1-Mater.nu,Mater.nu,Mater.nu; Mater.nu,1-Mater.nu,Mater.nu; Mater.nu,Mater.nu,1-Mater.nu];
Mater.C(4:6,4:6)= Mater.E/(1+Mater.nu)*eye(3)/2;
Mater.lamda = Mater.nu*Mater.E/((1+Mater.nu)*(1-2*Mater.nu)); % Lame Constant
Mater.mu = Mater.E/(2*(1+Mater.nu)); % Lame Constant

% Properties for Fracture
Fract = struct;
Fract.cenerg = 0.5; % Critical energy release for unstable crack (Gc)
Fract.constl = 0.1; % L0 : Length parameter which controls the spread of the damage

% Time Integration Paramenters
Integ = struct;
Integ.nstep = 1000; % Number of time steps
Integ.tfacto = 0; % Total increment factor for displacement increment
Integ.dfacto1 = 1e-4; % Displacement increment per time steps
Integ.dfacto2 = 5e-5; % Displacement increment per time steps
Integ.numStepsLimit = 30;
Integ.nprnt = 1; % Print frequency to output the results to file

end