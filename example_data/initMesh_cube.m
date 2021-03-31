function [PHTelem,controlPts,meshInfo] = initMesh_cube(geometry)

meshInfo = struct;
L = geometry.L;
W = geometry.W;
H = geometry.H;
numberElemU = geometry.numElemU;
numberElemV = geometry.numElemV;
numberElemW = geometry.numElemW;

% Initialize the PHT geometry on coarse mesh
controlPts = cell(geometry.numPatches,1);
PHTelem = cell(geometry.numPatches,1);
dimBasis = zeros(1,geometry.numPatches);
numElements = zeros(1,geometry.numPatches);
octupleList = cell(geometry.numPatches,1);

% Divide the patches along the x direction
xVertices = linspace(0,L,3);
zVertices = linspace(0,H,3);

patchCounter = 0;
patchIndexSet = [1,2,4,3];

for patchIndexZ = 1:2
    for patchIndexX = 1:2
        
        % Set the dimensions of the patch
        patchCounter = patchCounter + 1;
        patchIndex = patchIndexSet(patchCounter);
        patchMinX = xVertices(patchIndexX);
        patchMaxX = xVertices(patchIndexX+1);
        patchMinY = 0;
        patchMaxY = W;
        patchMinZ = zVertices(patchIndexZ);
        patchMaxZ = zVertices(patchIndexZ+1);
        
        % Initialize geometry on coarsest mesh
        coefs(1:3,1,1,1) = [patchMinX; patchMinY; patchMinZ];
        coefs(1:3,1,2,1) = [patchMinX; patchMaxY; patchMinZ];
        coefs(1:3,2,1,1) = [patchMaxX; patchMinY; patchMinZ];
        coefs(1:3,2,2,1) = [patchMaxX; patchMaxY; patchMinZ];
        coefs(1:3,1,1,2) = [patchMinX; patchMinY; patchMaxZ];
        coefs(1:3,1,2,2) = [patchMinX; patchMaxY; patchMaxZ];
        coefs(1:3,2,1,2) = [patchMaxX; patchMinY; patchMaxZ];
        coefs(1:3,2,2,2) = [patchMaxX; patchMaxY; patchMaxZ];
        coefs(4,1,1,1) = 1;
        coefs(4,1,2,1) = 1;
        coefs(4,2,1,1) = 1;
        coefs(4,2,2,1) = 1;
        coefs(4,1,1,2) = 1;
        coefs(4,1,2,2) = 1;
        coefs(4,2,1,2) = 1;
        coefs(4,2,2,2) = 1;
        
        knotU = [0 0 1 1];
        knotV = [0 0 1 1];
        knotW = [0 0 1 1];
        nurbs = nrbmak(coefs,{knotU,knotV,knotW});
        
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
        nurbs = nrbdegelev(nurbs,[p,q,r]-(nurbs.order-1));
        nurbs = nrbkntins(nurbs,{knotU(2:end-1) knotV(2:end-1) knotW(2:end-1)});
        
        % Repeat the knots to get C1 continuity
        nurbs = nrbkntins(nurbs,{knotU(2:end-1) knotV(2:end-1) knotW(2:end-1)});
        
        [PHTelem{patchIndex},controlPts{patchIndex},dimBasis(patchIndex),octupleList{patchIndex}] = genControlPts3D(nurbs,geometry);
        
        %pre-refine along the expected crack propagation path
        %if patchIndexX==2
        % Refine the area between (0.4...0.6) in the y-direction
        [octupleRef] = refineMeshCoordinates3D(PHTelem{patchIndex},controlPts{patchIndex},octupleList{patchIndex},geometry,L/2,0.55-0.01,0,W,0.45+0.01,0.55-0.01);
        [octupleList{patchIndex},PHTelem{patchIndex}, controlPts{patchIndex},dimBasis(patchIndex)] = refineMesh3D(octupleRef,octupleList{patchIndex},PHTelem{patchIndex},controlPts{patchIndex},geometry,dimBasis(patchIndex));
        
         % Refine the area between (0.45...0.55) in the y-direction
        [octupleRef] = refineMeshCoordinates3D(PHTelem{patchIndex},controlPts{patchIndex},octupleList{patchIndex},geometry,L/2,0.525-0.001,0,W,0.475+0.001,0.525-0.001);
        [octupleList{patchIndex},PHTelem{patchIndex}, controlPts{patchIndex},dimBasis(patchIndex)] = refineMesh3D(octupleRef,octupleList{patchIndex},PHTelem{patchIndex},controlPts{patchIndex},geometry,dimBasis(patchIndex));

         % Refine the area between (0.475...0.525) in the y-direction
        [octupleRef] = refineMeshCoordinates3D(PHTelem{patchIndex},controlPts{patchIndex},octupleList{patchIndex},geometry,L/2,0.5125-0.0001,0,W,0.4875+0.0001,0.5125-0.0001);
        [octupleList{patchIndex},PHTelem{patchIndex}, controlPts{patchIndex},dimBasis(patchIndex)] = refineMesh3D(octupleRef,octupleList{patchIndex},PHTelem{patchIndex},controlPts{patchIndex},geometry,dimBasis(patchIndex));

        numElements(patchIndex) = 8*size(octupleList{patchIndex},1);
    end
end

meshInfo.dimBasis = dimBasis;
meshInfo.octupleList = octupleList;
meshInfo.numElements = numElements;

end