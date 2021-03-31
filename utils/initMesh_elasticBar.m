function [PHTelem,controlPts,meshInfo] = initMesh_elasticBar(geometry)

p = geometry.p;
L = geometry.L;
numberElemU = geometry.numElemU;

% Initialize the PHT geometry on coarse mesh
controlPts = cell(geometry.numPatches,1);
PHTelem = cell(geometry.numPatches, 1);
dimBasis = zeros(1,geometry. numPatches);
tupleList = cell(geometry.numPatches,1);
% Divide the patches along the x direction
xVertices = linspace(-L,L,geometry.numPatches+1);
numberElements = 0;
for patchIndex = 1:geometry.numPatches
    
    % Set the dimensions of the patch
    patchMinX = xVertices(patchIndex);
    patchMaxX = xVertices(patchIndex+1);
    
    % Initialize geometry on coarsest mesh
    coefs(1:4,1) = [patchMinX;0;0;1];
    coefs(1:4,2) = [patchMaxX;0;0;1];
    
    knotU = [0 0 1 1];
    nurbs = nrbmak(coefs,knotU);
    p_init = nurbs.order(1)-1;
    % Refine into numberElemU knotspans
    knotU = linspace(0,1,numberElemU+1);
    numberElementsU = length(unique(knotU))-1;
    
    % Increase polynomial order
    nurbs = nrbdegelev(nurbs,p-p_init);
    nurbs = nrbkntins(nurbs,knotU(2:end-1));
    % Repeat the knots to get C1 continuity
    nurbs = nrbkntins(nurbs,knotU(2:end-1));
    [controlPts{patchIndex},PHTelem{patchIndex},dimBasis(patchIndex),tupleList{patchIndex}] = genControlPts1D(nurbs,p,numberElementsU);
    numberElements = numberElements + 2*size(tupleList{patchIndex},1);
end
meshInfo.dimBasis = dimBasis;
meshInfo.tupleList = tupleList;
meshInfo.numElements = numberElements;
end