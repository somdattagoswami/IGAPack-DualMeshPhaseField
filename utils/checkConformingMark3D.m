function [PHTelem,controlPts,dimBasis,octupleList,markRef] = checkConformingMark3D(PHTelem,controlPts,dimBasis,geometry,octupleList,markRef)
% Checks that the patches are conforming and if needed makes them conforming through mesh refinement
% patchBoundaries format:
% patchA, patchB, faceA, faceB
% patchA should be 1st
% edge format: 1-front, 2-right, 3-back, 4-left, 5-down, 6-up
patchBoundaries = geometry.patchBoundaries;
numBoundaries = size(patchBoundaries,1);
numPatches = length(PHTelem);

keepChecking = 1;
while keepChecking
    keepChecking = 0;
    %create/set nodesGlobal entries in all patches to be equal to local nodes
    %entries
    for patchIndex = 1:numPatches
        for elemIndex = 1:length(PHTelem{patchIndex})
            PHTelem{patchIndex}(elemIndex).nodesGlobal = PHTelem{patchIndex}(elemIndex).nodes;
        end
    end
    
    for boundaryIndex = 1:numBoundaries
        %get the nodes on the boundary edge in patchA and patchB
        patchAList = patchBoundaries{boundaryIndex,1};
        patchB = patchBoundaries{boundaryIndex,2};
        faceAList = patchBoundaries{boundaryIndex,3};
        faceBList = patchBoundaries{boundaryIndex,4};
        
        for indexPatch=1:length(patchAList)
            patchA = patchAList(indexPatch);
            faceA = faceAList(indexPatch);
            faceB = faceBList(indexPatch);
            
            octupleListA = octupleList{patchA};
            octupleListB = octupleList{patchB};
            
            [elementsA] = sortFaceElem(PHTelem{patchA},faceA);
            [elementsB] = sortFaceElem(PHTelem{patchB},faceB);
            
            [octupleRefA,octupleRefB] = makeConforming3D(PHTelem{patchA},PHTelem{patchB},elementsA,elementsB,faceA,faceB,octupleListA,octupleListB);
            indexOctupleA = find(octupleRefA > 0);
            indexOctupleB = find(octupleRefB > 0);
            
            if ~isempty(indexOctupleA)
                numNewPatches = length(indexOctupleA);
                disp(['In patch ', num2str(patchA), ' refining ',num2str(numNewPatches), ' octuples to keep conformity with patch ', num2str(patchB)])
                [octupleList{patchA},PHTelem{patchA},controlPts{patchA},dimBasis(patchA),markRefTemp] = refineMeshMark3D(octupleRefA,octupleListA,PHTelem{patchA},controlPts{patchA},geometry,dimBasis(patchA));
                keepChecking = 1;
                indexRef = find(markRefTemp > 0);
                if indexRef
                    markRef{patchA}(indexRef) = 1;
                end
            end
            
            if ~isempty(indexOctupleB)
                numNewPatches = length(indexOctupleB);
                disp(['In patch ', num2str(patchB), ' refining ',num2str(numNewPatches), ' octuples to keep conformity with patch ', num2str(patchA)])
                [octupleList{patchB},PHTelem{patchB},controlPts{patchB},dimBasis(patchB),markRefTemp] = refineMeshMark3D(octupleRefB,octupleListB,PHTelem{patchB},controlPts{patchB},geometry,dimBasis(patchB));
                keepChecking = 1;
                indexRef = find(markRefTemp > 0);
                if indexRef
                    markRef{patchB}(indexRef) = 1;
                end
            end
        end
    end
end