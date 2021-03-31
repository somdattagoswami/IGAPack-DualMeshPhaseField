% Script make_c0_cube_with_hole.m
% This requires the nurbs toolbox to be in the path

addpath ../nurbs/inst
close all
clear 
clc
% Start with C0 plate with a hole

% Initialize geometry on coarsest mesh
rad = 1;
side_fac = 1;

knotU = [0, 0, 0, 0.5, 0.5, 1, 1, 1];
knotV = [0, 0, 1, 1];

coefs(1:3,1,1) = rad.*[-1;0;0];
coefs(1:3,2,1) = rad.*[-0.853553390593274; 0.353553390593274; 0];
coefs(1:3,3,1) = rad.*[-0.603553390593274; 0.603553390593274; 0];
coefs(1:3,4,1) = rad.*[-0.353553390593274;0.853553390593274;0];
coefs(1:3,5,1) = rad.*[0;1;0];

coefs(1:3,1,2) = side_fac.*[-2;0;0];
coefs(1:3,2,2) = side_fac.*[-2;1;0];
coefs(1:3,3,2) = side_fac.*[-2;2;0];
coefs(1:3,4,2) = side_fac.*[-1;2;0];
coefs(1:3,5,2) = side_fac.*[0;2;0];


coefs(4,1,1) = 1;
coefs(4,2,1) = 0.853553390593274;
coefs(4,3,1) = 0.853553390593274;
coefs(4,4,1) = 0.853553390593274;
coefs(4,5,1) = 1;

coefs(4,1,2) = 1;
coefs(4,2,2) = 1;
coefs(4,3,2) = 1;
coefs(4,4,2) = 1;
coefs(4,5,2) = 1;
srf1 = nrbmak(coefs,{knotU,knotV});

vol1 = nrbrevolve(srf1,[0,0,0],[0,1,0],pi/2);
vol1 = nrbkntins(vol1,{[0.5,0.5],[],[]});
vol1.coefs(:,3,1,2) = [-2;0;2;1];
vol1.coefs(:,3,2,2) = [-3.414213562373095/2;1.707106781186548/2;3.414213562373095/2;0.853553390593274];
vol1.coefs(:,3,3,2) = [-2;2;2;1];

% Split in 4 patches
coefs = vol1.coefs;
knotU = [0,0,0,1,1,1];
knotV = [0,0,0,1,1,1];
knotW = [0,0,1,1];
coefsSW1 = coefs(:,1:3,1:3,:);
coefsSE1 = coefs(:,3:5,1:3,:);
coefsNW1 = coefs(:,1:3,3:5,:);
coefsNE1 = coefs(:,3:5,3:5,:);

solidSW1 = nrbmak(coefsSW1,{knotU,knotV,knotW});
solidSE1 = nrbmak(coefsSE1,{knotU,knotV,knotW});
solidNW1 = nrbmak(coefsNW1,{knotU,knotV,knotW});
solidNE1 = nrbmak(coefsNE1,{knotU,knotV,knotW});

% Making the inner hole patch to make the structure rigid
knotU = [0, 0, 0, 0.5, 0.5, 1, 1, 1];
knotV = [0, 0, 1, 1];

inRad = 0.5;
coefs = zeros(4,5,2);
coefs(1:3,1,1) = inRad.*[-1;0;0];
coefs(1:3,2,1) = inRad.*[-0.853553390593274; 0.353553390593274; 0];
coefs(1:3,3,1) = inRad.*[-0.603553390593274; 0.603553390593274; 0];
coefs(1:3,4,1) = inRad.*[-0.353553390593274;0.853553390593274;0];
coefs(1:3,5,1) = inRad.*[0;1;0];

coefs(1:3,1,2) = rad.*[-1;0;0];
coefs(1:3,2,2) = rad.*[-0.853553390593274; 0.353553390593274; 0];
coefs(1:3,3,2) = rad.*[-0.603553390593274; 0.603553390593274; 0];
coefs(1:3,4,2) = rad.*[-0.353553390593274;0.853553390593274;0];
coefs(1:3,5,2) = rad.*[0;1;0];

coefs(4,1,1) = 1;
coefs(4,2,1) = 0.853553390593274;
coefs(4,3,1) = 0.853553390593274;
coefs(4,4,1) = 0.853553390593274;
coefs(4,5,1) = 1;

coefs(4,1,2) = 1;
coefs(4,2,2) = 0.853553390593274;
coefs(4,3,2) = 0.853553390593274;
coefs(4,4,2) = 0.853553390593274;
coefs(4,5,2) = 1;
srf2 = nrbmak(coefs,{knotU,knotV});

vol2 = nrbrevolve(srf2,[0,0,0],[0,1,0],pi/2);
vol2 = nrbkntins(vol2,{[0.5,0.5],[],[]});

%split in 4 patches
coefs = vol2.coefs;
knotU = [0,0,0,1,1,1];
knotV = [0,0,0,1,1,1];
knotW = [0,0,1,1];
coefsSW2 = coefs(:,1:3,1:3,:);
coefsSE2 = coefs(:,3:5,1:3,:);
coefsNW2 = coefs(:,1:3,3:5,:);
coefsNE2 = coefs(:,3:5,3:5,:);

solidSW2 = nrbmak(coefsSW2,{knotU,knotV,knotW});
solidSE2 = nrbmak(coefsSE2,{knotU,knotV,knotW});
solidNW2 = nrbmak(coefsNW2,{knotU,knotV,knotW});
solidNE2 = nrbmak(coefsNE2,{knotU,knotV,knotW});

savefile = 'cube_with_hole_C0.mat';
save(savefile,'solidNE1','solidNW1','solidSE1','solidSW1','solidNE2','solidNW2','solidSE2','solidSW2');

% figure
% nrbctrlplot(solidSW1)
% hold on
% nrbctrlplot(solidSE1)
% hold on
% nrbctrlplot(solidNW1)
% hold on
% nrbctrlplot(solidNE1)
% hold on
% nrbctrlplot(solidSW2)
% hold on
% nrbctrlplot(solidSE2)
% hold on
% nrbctrlplot(solidNW2)
% hold on
% nrbctrlplot(solidNE2)