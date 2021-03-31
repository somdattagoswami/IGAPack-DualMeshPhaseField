function [PHTelem,controlPts,dimBasis,fieldDataPatch,numberElements,markRef] = checkConformingProjElemMark3D(PHTelem,controlPts,dimBasis,geometry,fieldDataPatch,numberElements,markRef)
%checks that the patches are conforming and if needed makes them conforming through mesh refinement
%patchBoundaries format:
% patchA, patchB, edgeA, edgeB
% patchA should be 1
%edge format: 1-down, 2-right, 3-up, 4-left
%transfer field data at patch level

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
    
    % loop over the boundaries and match patches if needed
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
            
            [elementsA] = sortFaceElem(PHTelem{patchA},faceA);
            [elementsB] = sortFaceElem(PHTelem{patchB},faceB);
            
            [elemRefA,elemRefB] = makeConformingElem3D(PHTelem{patchA},PHTelem{patchB},elementsA,elementsB,faceA,faceB);
            indexElemA = find(elemRefA > 0);
            indexElemB = find(elemRefB > 0);
            
            if ~isempty(indexElemA)
                keepChecking = 1;
                numNewElems = length(indexElemA);
                disp(['In patch ', num2str(patchA), ' refining ',num2str(numNewElems), ' elements to keep conformity with patch ', num2str(patchB)])
                [PHTelem{patchA},controlPts{patchA},dimBasis(patchA),fieldDataPatch{patchA},numberElements(patchA),markRefTemp] = refineElemProjGradedIso3D(elemRefA,PHTelem{patchA},controlPts{patchA},geometry,dimBasis(patchA),fieldDataPatch{patchA},numberElements(patchA));
                indexRef = find(markRefTemp > 0);
                if indexRef
                    markRef{patchA}(indexRef) = 1;
                end
            end
            
            if ~isempty(indexElemB)
                keepChecking = 1;
                numNewElems = length(indexElemB);
                disp(['In patch ', num2str(patchB), ' refining ',num2str(numNewElems), ' elements to keep conformity with patch ', num2str(patchA)])
                [PHTelem{patchB},controlPts{patchB},dimBasis(patchB),fieldDataPatch{patchB},numberElements(patchB),markRefTemp] = refineElemProjGradedIso3D(elemRefB,PHTelem{patchB},controlPts{patchB},geometry,dimBasis(patchB),fieldDataPatch{patchB},numberElements(patchB));
                indexRef = find(markRefTemp > 0);
                if indexRef
                    markRef{patchB}(indexRef) = 1;
                end
            end
        end
    end
end