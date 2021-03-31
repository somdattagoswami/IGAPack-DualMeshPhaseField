function [PHTelem,controlPts,meshInfo] = initMesh_3ptPlate(geometry)

meshInfo = struct;
numPatches = geometry.numPatches;
p = geometry.p;
q = geometry.q;
L = geometry.L;
W = geometry.W;
notchHt = geometry.notchHt;
controlPts = cell(numPatches,1);
PHTelem = cell(numPatches,1);
quadList = cell(numPatches,1);
dimBasis = zeros(1, numPatches);

% Divide the patches along the x direction
xCord = [0,(0.5*L-0.1),0.5*L,0; 0,0.5*L,0.5*L,0; 0.5*L,L,L,0.5*L; (0.5*L+0.1),L,L,0.5*L];
yCord = [0,0,notchHt,notchHt; notchHt,notchHt,W,W; notchHt,notchHt,W,W;0,0,notchHt,notchHt];

numberElements = 0;
for patchIndex = 1:numPatches
    
     % Initialize geometry on coarsest mesh
    coefs(1:3,1,1) = [xCord(patchIndex,1); yCord(patchIndex,1);0];
    coefs(1:3,1,2) = [xCord(patchIndex,4); yCord(patchIndex,4);0];
    coefs(1:3,2,1) = [xCord(patchIndex,2); yCord(patchIndex,2);0];
    coefs(1:3,2,2) = [xCord(patchIndex,3); yCord(patchIndex,3);0];
    coefs(4,1,1) = 1;
    coefs(4,1,2) = 1;
    coefs(4,2,1) = 1;
    coefs(4,2,2) = 1;
    
    knotU = [0 0 1 1];
    knotV = [0 0 1 1];
    
    nurbs = nrbmak(coefs,{knotU,knotV});
    p_init = nurbs.order(1)-1;
    q_init = nurbs.order(2)-1;
    
    if patchIndex == 1 || patchIndex == 4
        numberElemU = 12;
        numberElemV = 2;
    else
        numberElemU = 12;
        numberElemV = 6;
    end
    
    % Refine into numberElemU by numberElemV knotspans
    knotU = linspace(0,1,numberElemU+1);
    knotV = linspace(0,1,numberElemV+1);
    %     nurbs = nrbkntins(nrb,{knotU(2:end-1) knotV(2:end-1)});
    numberElementsU = length(unique(knotU))-1;
    numberElementsV = length(unique(knotV))-1;
%     numberElements = numberElements + numberElementsU*numberElementsV;
    
    % Increase polynomial order
    nurbs = nrbdegelev(nurbs,[p-p_init,q-q_init]);
    nurbs = nrbkntins(nurbs,{knotU(2:end-1) knotV(2:end-1)});
    
    % Repeat the knots to get C1 continuity
    nurbs = nrbkntins(nurbs,{knotU(2:end-1) knotV(2:end-1)});
    [controlPts{patchIndex},PHTelem{patchIndex},dimBasis(patchIndex),quadList{patchIndex}] = genControlPtsNoRep(nurbs,p,q,numberElementsU,numberElementsV);
    
    % Refine the area between (0.45...0.55) in the y-direction
    [quadRef] = refineMeshCoordinates(PHTelem{patchIndex},controlPts{patchIndex},quadList{patchIndex},p,q,3.95+0.01,4.1-0.01,0.45+0.01,0.55-0.01);
    [quadList{patchIndex},PHTelem{patchIndex}, controlPts{patchIndex},dimBasis(patchIndex)] = refineMeshIso(quadRef,quadList{patchIndex},PHTelem{patchIndex},controlPts{patchIndex},p,q,dimBasis(patchIndex));
    
    % Refine the area between (0.495...0.505) in the y-direction
    [quadRef] = refineMeshCoordinates(PHTelem{patchIndex},controlPts{patchIndex},quadList{patchIndex},p,q,3.975+0.001,4.025-0.001,0.49+0.001,0.51-0.001);
    [quadList{patchIndex},PHTelem{patchIndex}, controlPts{patchIndex},dimBasis(patchIndex)] = refineMeshIso(quadRef,quadList{patchIndex},PHTelem{patchIndex},controlPts{patchIndex},p,q,dimBasis(patchIndex));
    
    % Refine the area between (0.4995...0.5005) in the y-direction
    [quadRef] = refineMeshCoordinates(PHTelem{patchIndex},controlPts{patchIndex},quadList{patchIndex},p,q,3.9875+0.0001,4.0125-0.0001,0.4995+0.0001,0.5005-0.0001);
    [quadList{patchIndex},PHTelem{patchIndex},controlPts{patchIndex},dimBasis(patchIndex)] = refineMeshIso(quadRef, quadList{patchIndex}, PHTelem{patchIndex}, controlPts{patchIndex},p,q,dimBasis(patchIndex));
    
    numberElements = numberElements + 4*size(quadList{patchIndex},1);
end

meshInfo.dimBasis = dimBasis;
meshInfo.quadList = quadList;
meshInfo.numElements = numberElements;

end