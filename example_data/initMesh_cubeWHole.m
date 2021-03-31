function [PHTelem,controlPts,meshInfo] = initMesh_cubeWHole(geometry)

numberElements = 0;
load('cube_with_hole_C0.mat')

[PHTelem{1},controlPts{1},dimBasis(1),octupleList{1}] = genGeometry(solidSW1,geometry);
numberElements(1) = 8*size(octupleList{1}(:,1),1);

[PHTelem{2},controlPts{2},dimBasis(2),octupleList{2}] = genGeometry(solidSE1,geometry);
numberElements(2) = 8*size(octupleList{2}(:,1),1);

[PHTelem{3},controlPts{3},dimBasis(3),octupleList{3}] = genGeometry(solidNW1,geometry);
numberElements(3) = 8*size(octupleList{3}(:,1),1);

[PHTelem{4},controlPts{4},dimBasis(4),octupleList{4}] = genGeometry(solidNE1,geometry);
numberElements(4) = 8*size(octupleList{4}(:,1),1);

meshInfo = struct;
meshInfo.dimBasis = dimBasis;
meshInfo.octupleList = octupleList;
meshInfo.numElements = numberElements;

end

function[PHTelem,controlPts,dimBasis,octupleList] = genGeometry(solid,geometry)

numberElemU = 1;
numberElemV = 1;
numberElemW = 1;

% Refine into numberElemU by numberElemV knotspans
knotU = linspace(0,1,numberElemU+1);
knotV = linspace(0,1,numberElemV+1);
knotW = linspace(0,1,numberElemW+1);

geometry.numElmtU = length(unique(knotU))-1;
geometry.numElmtV = length(unique(knotV))-1;
geometry.numElmtW = length(unique(knotW))-1;

% Increase polynomial order
p = geometry.p;
q = geometry.q;
r = geometry.r;
solid = nrbdegelev(solid,[p,q,r]-(solid.order-1));
solid = nrbkntins(solid,{knotU(2:end-1) knotV(2:end-1) knotW(2:end-1)});

% Repeat the knots to get C1 continuity
solid = nrbkntins(solid,{knotU(2:end-1) knotV(2:end-1) knotW(2:end-1)});

[PHTelem,controlPts,dimBasis,octupleList] = genControlPts3D(solid,geometry);

end