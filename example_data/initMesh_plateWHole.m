function [PHTelem,controlPts,meshInfo] = initMesh_plateWHole(geometry)

meshInfo = struct;
numPatches = geometry.numPatches;
p = geometry.p;
q = geometry.q;
nurbsList = cell(1,numPatches);
controlPts = cell(numPatches,1);
PHTelem = cell(numPatches,1);
quadList = cell(numPatches,1);
numElemU = geometry.numElemU;
numElemV = geometry.numElemV;
dimBasis = zeros(1, numPatches);
% Data for rectangular patches
% patchCordsRect = [PatchIndex, PatchXMin, PatchXMax, PatchYMin, PatchYMax]
patchCordsRect = [1,0,1,0,1; 2,0,1,1,3.75; 3,0,1,3.75,5.75; 4,0,1,5.75,8;
    5,1,4,0,1; 6,1,4,1,3.75; 7,1,4,3.75,5.75; 8,1,4,5.75,8;
    21,4,8,0,1;
    22,8,10,5.75,8; 23,8,10,3.75,5.75; 24,8,10,1,3.75; 25,8,10,0,1;
    26,10,19,5.75,8; 27,10,19,3.75,5.75; 28,10,19,1,3.75; 29,10,19,0,1;
    30,19,20,5.75,8; 31,19,20,3.75,5.75; 32,19,20,1,3.75; 33,19,20,0,1;];

% patchCordsTri = [PatchIndex for the 1st Quad in the square, PatchXMin, PatchXMax, PatchYMin, PatchYMax, X-Centre of the circle, Y-Centre of the circle]
patchCordsTri = [9,4,8,5.75,8,6,6.75; 13,4,8,3.75,5.75,6,4.75; 17,4,8,1,3.75,6,2.75];
numberElements = 0;

for iPatch = 1:size(patchCordsRect,1)
    coefs = zeros(4,2,2);
    patchIndex = patchCordsRect(iPatch,1);
    dimBasis(patchIndex) = (p+1)*(q+1); %coarsest 1x1 mesh has dimension (p+1)*(q+1)
    
    % Set the dimensions of the patch
    patchMinX = patchCordsRect(iPatch,2);
    patchMaxX = patchCordsRect(iPatch,3);
    patchMinY = patchCordsRect(iPatch,4);
    patchMaxY = patchCordsRect(iPatch,5);
    
    % Initialize geometry on coarsest mesh
    coefs(1:3,1,1) = [patchMinX; patchMinY; 0];
    coefs(1:3,1,2) = [patchMinX; patchMaxY; 0];
    coefs(1:3,2,1) = [patchMaxX; patchMinY; 0];
    coefs(1:3,2,2) = [patchMaxX; patchMaxY; 0];
    coefs(4,1,1) = 1;
    coefs(4,1,2) = 1;
    coefs(4,2,1) = 1;
    coefs(4,2,2) = 1;
    knotU = [0 0 1 1];
    knotV = [0 0 1 1];
    nurbs = nrbmak(coefs,{knotU,knotV});
    nurbsList{patchIndex} = nurbs;
    [controlPts{patchIndex},PHTelem{patchIndex},dimBasis(patchIndex),quadList{patchIndex}] = initElements(nurbs,p,q,numElemU,numElemV);
    numberElements = numberElements + 4*size(quadList{patchIndex},1);
end

for iPatch = 1:size(patchCordsTri,1)
    
    patchMinX = patchCordsTri(iPatch,2);
    patchMaxX = patchCordsTri(iPatch,3);
    patchMinY = patchCordsTri(iPatch,4);
    patchMaxY = patchCordsTri(iPatch,5);
    centreX = patchCordsTri(iPatch,6);
    centreY = patchCordsTri(iPatch,7);
    rad = 0.25;
    
    patchIndex = patchCordsTri(iPatch,1);
    startAngle = pi/4;
    endAngle = 3*pi/4;
    arcSide = nrbcirc(rad, [centreX, centreY, 0], startAngle, endAngle);
    arcSide.coefs = fliplr(arcSide.coefs);
    startPtX = arcSide.coefs(1,1);
    startPtY = arcSide.coefs(2,1);
    side1 = nrbline([startPtX, startPtY],[patchMinX,patchMaxY]);
    side2 = nrbline([patchMinX, patchMaxY],[patchMaxX, patchMaxY]);
    endPtX = arcSide.coefs(1,end);
    endPtY = arcSide.coefs(2,end);
    side3 = nrbline([endPtX, endPtY],[patchMaxX, patchMaxY]);
    nurbs = nrbcoons(arcSide,side2,side1,side3);
    % flip the u-direction to maintain consistent direction across patch
    % boundaries
    if iPatch==2
        nurbs = nrbreverse(nurbs,[1]);
    end
    nurbsList{patchIndex} = nurbs;
    [controlPts{patchIndex},PHTelem{patchIndex},dimBasis(patchIndex),quadList{patchIndex}] = initElements(nurbs,p,q,numElemU,numElemV);
    numberElements = numberElements + 4*size(quadList{patchIndex},1);
    
    patchIndex = patchCordsTri(iPatch,1)+ 1;
    startAngle = 3*pi/4;
    endAngle = 5*pi/4;
    arcSide = nrbcirc(rad, [centreX, centreY, 0], startAngle, endAngle);
    arcSide.coefs = fliplr(arcSide.coefs);
    startPtX = arcSide.coefs(1,1);
    startPtY = arcSide.coefs(2,1);
    side1 = nrbline([startPtX, startPtY],[patchMinX,patchMinY]);
    side2 = nrbline([patchMinX,patchMinY], [patchMinX, patchMaxY]);
    endPtX = arcSide.coefs(1,end);
    endPtY = arcSide.coefs(2,end);
    side3 = nrbline([endPtX, endPtY],[patchMinX, patchMaxY]);
    nurbs = nrbcoons(arcSide,side2,side1,side3);
    nurbsList{patchIndex} = nurbs;
    [controlPts{patchIndex},PHTelem{patchIndex},dimBasis(patchIndex),quadList{patchIndex}] = initElements(nurbs,p,q,numElemU,numElemV);
    numberElements = numberElements + 4*size(quadList{patchIndex},1);
    
    patchIndex = patchCordsTri(iPatch,1)+ 2;
    startAngle = 5*pi/4;
    endAngle = 7*pi/4;
    arcSide = nrbcirc(rad, [centreX, centreY, 0], startAngle, endAngle);
    arcSide.coefs = fliplr(arcSide.coefs);
    startPtX = arcSide.coefs(1,1);
    startPtY = arcSide.coefs(2,1);
    side1 = nrbline([startPtX, startPtY],[patchMaxX,patchMinY]);
    side2 = nrbline([patchMaxX,patchMinY] , [patchMinX, patchMinY]);
    endPtX = arcSide.coefs(1,end);
    endPtY = arcSide.coefs(2,end);
    side3 = nrbline([endPtX, endPtY],[patchMinX, patchMinY]);
    nurbs = nrbcoons(arcSide,side2,side1,side3);
    %flip the u-direction to maintain consistent direction across patch
    %boundaries
    if (iPatch == 2) || (iPatch == 3)
        nurbs = nrbreverse(nurbs,[1]);
    end
    nurbsList{patchIndex} = nurbs;
    [controlPts{patchIndex},PHTelem{patchIndex},dimBasis(patchIndex),quadList{patchIndex}] = initElements(nurbs,p,q,numElemU,numElemV);
    numberElements = numberElements + 4*size(quadList{patchIndex},1);
    
    patchIndex = patchCordsTri(iPatch,1)+ 3;
    startAngle = 7*pi/4;
    endAngle = pi/4;
    arcSide = nrbcirc(rad, [centreX, centreY, 0], startAngle, endAngle);
    arcSide.coefs = fliplr(arcSide.coefs);
    startPtX = arcSide.coefs(1,1);
    startPtY = arcSide.coefs(2,1);
    side1 = nrbline([startPtX, startPtY],[patchMaxX,patchMaxY]);
    side2 = nrbline([patchMaxX,patchMaxY], [patchMaxX, patchMinY]);
    endPtX = arcSide.coefs(1,end);
    endPtY = arcSide.coefs(2,end);
    side3 = nrbline([endPtX, endPtY],[patchMaxX, patchMinY]);
    nurbs = nrbcoons(arcSide,side2,side1,side3);
    %flip the u-direction to maintain consistent direction across patch
    %boundaries
    nurbs = nrbreverse(nurbs,[1]);
    nurbsList{patchIndex} = nurbs;
    [controlPts{patchIndex},PHTelem{patchIndex},dimBasis(patchIndex),quadList{patchIndex}] = initElements(nurbs,p,q,numElemU,numElemV);
    numberElements = numberElements + 4*size(quadList{patchIndex},1);
end

meshInfo.dimBasis = dimBasis;
meshInfo.quadList = quadList;
meshInfo.numElements = numberElements;

end