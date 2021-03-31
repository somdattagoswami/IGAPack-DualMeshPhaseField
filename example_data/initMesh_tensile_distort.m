function [PHTelem,controlPts,meshInfo] = initMesh_tensile_distort(geometry)

meshInfo = struct;
numPatches = geometry.numPatches;
p = geometry.p;
q = geometry.q;
numberElemU = geometry.numElemU;
numberElemV = geometry.numElemV;

% Initialize the PHT geometry on coarse mesh
controlPts = cell(numPatches,1);
PHTelem = cell(numPatches, 1);
dimBasis = zeros(1, numPatches);
quadList = cell(numPatches,1);

numberElements = 0;
patchCords = [1,0,0,0,0.5,0.5,0.5,0.7,0; 2,0.7,0,0.5,0.5,1,0.7,1,0; 3,0.5,0.5,0.7,1,1,1,1,0.7; 4,0,0.5,0,1,0.7,1,0.5,0.5];

for iPatch = 1:size(patchCords,1)
    
    patchIndex = patchCords(iPatch,1);
    %set the dimensions of the patch
    side1 = nrbline([patchCords(iPatch,2), patchCords(iPatch,3)],[patchCords(iPatch,8), patchCords(iPatch,9)]);
    side2 = nrbline([patchCords(iPatch,2), patchCords(iPatch,3)],[patchCords(iPatch,4), patchCords(iPatch,5)]);
    side3 = nrbline([patchCords(iPatch,4), patchCords(iPatch,5)],[patchCords(iPatch,6), patchCords(iPatch,7)]);
    side4 = nrbline([patchCords(iPatch,8), patchCords(iPatch,9)],[patchCords(iPatch,6), patchCords(iPatch,7)]);
    nurbs = nrbcoons(side1,side3,side2,side4);
    
    p_init = nurbs.order(1)-1;
    q_init = nurbs.order(2)-1;
    
    % Refine into numberElemU by numberElemV knotspans
    knotU = linspace(0,1,numberElemU+1);
    knotV = linspace(0,1,numberElemV+1);
    %     nurbs = nrbkntins(nrb,{knotU(2:end-1) knotV(2:end-1)});
    numberElementsU = length(unique(knotU))-1;
    numberElementsV = length(unique(knotV))-1;
    
    % Increase polynomial order
    nurbs = nrbdegelev(nurbs,[p-p_init,q-q_init]);
    nurbs = nrbkntins(nurbs,{knotU(2:end-1) knotV(2:end-1)});
    % Repeat the knots to get C1 continuity
    nurbs = nrbkntins(nurbs,{knotU(2:end-1) knotV(2:end-1)});
    [controlPts{patchIndex},PHTelem{patchIndex},dimBasis(patchIndex),quadList{patchIndex}] = genControlPtsNoRep(nurbs,p,q,numberElementsU,numberElementsV);
    
    % Refine the area between (0.45...0.55) in the y-direction
    [quadRef] = refineMeshCoordinates(PHTelem{patchIndex},controlPts{patchIndex},quadList{patchIndex},p,q,0.45+0.01,0.55-0.01,0.45+0.01,0.55-0.01);
    [quadList{patchIndex},PHTelem{patchIndex}, controlPts{patchIndex},dimBasis(patchIndex)] = refineMeshIso(quadRef,quadList{patchIndex},PHTelem{patchIndex},controlPts{patchIndex},p,q,dimBasis(patchIndex));
    
    % Refine the area between (0.475...0.525) in the y-direction
    [quadRef] = refineMeshCoordinates(PHTelem{patchIndex},controlPts{patchIndex},quadList{patchIndex},p,q,0.475+0.001,0.525-0.001,0.475+0.001,0.525-0.001);
    [quadList{patchIndex},PHTelem{patchIndex}, controlPts{patchIndex},dimBasis(patchIndex)] = refineMeshIso(quadRef,quadList{patchIndex},PHTelem{patchIndex},controlPts{patchIndex},p,q,dimBasis(patchIndex));
    
    % Refine the area between (0.4875...0.5125) in the y-direction
    [quadRef] = refineMeshCoordinates(PHTelem{patchIndex},controlPts{patchIndex},quadList{patchIndex},p,q,0.4875+0.0001,0.5125-0.0001,0.4875+0.0001,0.5125-0.0001);
    [quadList{patchIndex},PHTelem{patchIndex},controlPts{patchIndex},dimBasis(patchIndex)] = refineMeshIso(quadRef, quadList{patchIndex}, PHTelem{patchIndex}, controlPts{patchIndex},p,q,dimBasis(patchIndex));
    numberElements = numberElements + 4*size(quadList{patchIndex},1);
end

meshInfo.dimBasis = dimBasis;
meshInfo.quadList = quadList;
meshInfo.numElements = numberElements;

end