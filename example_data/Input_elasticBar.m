function [geometry,Mater,Fract] = Input_elasticBar
% Input File for 1D elastic bar
% use C1 cubics with adaptivity
geometry = struct;
geometry.dim = 1;

% Considering linear B-Splines in u direction
geometry.p = 3; % Degree of the polynomial in U Direction
geometry.L = 1; % Length of the plate in the U direction
geometry.numPatches = 1; % Total number of patches
geometry.ngaussX = geometry.p+1;
geometry.numElemU = 12; % Initial mesh
geometry.maxERefLevel = 6; % Maximum level of Refinement for Elastic Mesh
geometry.maxPhRefLevel = 8; % Maximum level of Refinement for Phase Field Mesh
geometry.threshPhi = 0.5; % Threshold for Refinement
geometry.B = 1e3;% Parameter for initial history function
geometry.toler = 7*1e-12;

% Material properties
Mater = struct;
Mater.E = 1; % Young's Modulus based on (N/mm2)
Mater.nu = 0.3; % Poisson's Ratio

% Properties for Fracture
Fract = struct;
Fract.cenerg = 2.7; % Critical energy release for unstable crack (Gc)
Fract.constl = 0.0125/2; % L0 : Length parameter which controls the spread of the damage

end