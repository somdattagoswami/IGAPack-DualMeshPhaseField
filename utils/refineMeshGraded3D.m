function [patchList,PHTelem,controlPts,dimBasis,markRef] = refineMeshGraded3D(patchRef,patchList,PHTelem,controlPts,geometry,dimBasis)
%refines the patches marked by patchRef
%also refine the neighbor patches if the level_neighbor<level_refined_quad

numElem = length(PHTelem);
markRef = zeros(length(PHTelem),1);
indexPatch = find(patchRef > 0);
numNewPatches = length(indexPatch);
newPatchList = zeros(size(patchList,1)+numNewPatches,8);

keepChecking = 1;
% Check the refinement level for the neighbor octuple and mark it for refinement if needed
extraPatchCounter = 0;
while keepChecking
    keepChecking = 0;
    patchRefOrig = patchRef;
    for i=1:length(patchRefOrig)
        if patchRefOrig(i)
            curLevel = PHTelem(patchList(i,1)).level;
            for iPatchElem = 1:8
                curElem = patchList(i,iPatchElem);
                if ~isempty(PHTelem(curElem).neighbor_front)
                    neighborFront = PHTelem(curElem).neighbor_front(1);
                    if curLevel>PHTelem(neighborFront).level                    
                        [patchIndex, ~] = find(patchList==neighborFront);
                        if patchRef(patchIndex)==0
                            patchRef(patchIndex) = 1;
                            extraPatchCounter = extraPatchCounter + 1;
                            keepChecking = 1;
                        end
                    end
                end
                if ~isempty(PHTelem(curElem).neighbor_down)
                    neighborDown = PHTelem(curElem).neighbor_down(1);
                    if curLevel>PHTelem(neighborDown).level
                        [patchIndex, ~] = find(patchList==neighborDown);
                        if patchRef(patchIndex)==0
                            patchRef(patchIndex) = 1;
                            extraPatchCounter = extraPatchCounter + 1;
                            keepChecking = 1;
                        end
                    end
                end
                if ~isempty(PHTelem(curElem).neighbor_left)
                    neighborLeft = PHTelem(curElem).neighbor_left(1);
                    if curLevel>PHTelem(neighborLeft).level
                        [patchIndex, ~] = find(patchList==neighborLeft);
                        if patchRef(patchIndex)==0
                            patchRef(patchIndex) = 1;
                            extraPatchCounter = extraPatchCounter + 1;
                            keepChecking = 1;
                        end
                    end
                end
                if ~isempty(PHTelem(curElem).neighbor_down_left)
                    neighborDownLeft = PHTelem(curElem).neighbor_down_left(1);
                    if curLevel>PHTelem(neighborDownLeft).level
                        [patchIndex, ~] = find(patchList==neighborDownLeft);
                        if patchRef(patchIndex)==0
                            patchRef(patchIndex) = 1;
                            extraPatchCounter = extraPatchCounter + 1;
                            keepChecking = 1;
                        end
                    end
                end
                if ~isempty(PHTelem(curElem).neighbor_left_front)
                    neighborLeftFront = PHTelem(curElem).neighbor_left_front(1);
                    if curLevel>PHTelem(neighborLeftFront).level
                        [patchIndex, ~] = find(patchList==neighborLeftFront);
                        if patchRef(patchIndex)==0
                            patchRef(patchIndex) = 1;
                            extraPatchCounter = extraPatchCounter + 1;
                            keepChecking = 1;
                        end
                    end
                end
                
                if ~isempty(PHTelem(curElem).neighbor_down_right)
                    neighborDownRight = PHTelem(curElem).neighbor_down_right(1);
                    if curLevel>PHTelem(neighborDownRight).level
                        [patchIndex, ~] = find(patchList==neighborDownRight);
                        if patchRef(patchIndex)==0
                            patchRef(patchIndex) = 1;
                            extraPatchCounter = extraPatchCounter + 1;
                            keepChecking = 1;
                        end
                    end
                end
                if ~isempty(PHTelem(curElem).neighbor_left_back)
                    %disp('West')
                    neighborLeftBack = PHTelem(curElem).neighbor_left_back(1);
                    if curLevel>PHTelem(neighborLeftBack).level
                        [patchIndex, ~] = find(patchList==neighborLeftBack);
                        if patchRef(patchIndex)==0
                            patchRef(patchIndex) = 1;
                            extraPatchCounter = extraPatchCounter + 1;
                            keepChecking = 1;
                        end
                    end
                end
                if ~isempty(PHTelem(curElem).neighbor_down_front)
                    neighborDownFront = PHTelem(curElem).neighbor_down_front(1);
                    if curLevel>PHTelem(neighborDownFront).level
                        [patchIndex, ~] = find(patchList==neighborDownFront);
                        if patchRef(patchIndex)==0
                            patchRef(patchIndex) = 1;
                            extraPatchCounter = extraPatchCounter + 1;
                            keepChecking = 1;
                        end
                    end
                end
                if ~isempty(PHTelem(curElem).neighbor_down_back)
                    neighborDownBack = PHTelem(curElem).neighbor_down_back(1);
                    if curLevel>PHTelem(neighborDownBack).level
                        [patchIndex, ~] = find(patchList==neighborDownBack);
                        if patchRef(patchIndex)==0
                            patchRef(patchIndex) = 1;
                            extraPatchCounter = extraPatchCounter + 1;
                            keepChecking = 1;
                        end
                    end
                end
                if ~isempty(PHTelem(curElem).neighbor_back)
                    %disp('North')
                    neighborBack = PHTelem(curElem).neighbor_back(1);
                    if curLevel>PHTelem(neighborBack).level
                        [patchIndex, ~] = find(patchList==neighborBack);
                        if patchRef(patchIndex)==0
                            patchRef(patchIndex) = 1;
                            extraPatchCounter = extraPatchCounter + 1;
                            keepChecking = 1;
                        end
                    end
                end
                
                if ~isempty(PHTelem(curElem).neighbor_up)
                    %disp('North')
                    neighborUp = PHTelem(curElem).neighbor_up(1);
                    if curLevel>PHTelem(neighborUp).level
                        [patchIndex, ~] = find(patchList==neighborUp);
                        if patchRef(patchIndex)==0
                            patchRef(patchIndex) = 1;
                            extraPatchCounter = extraPatchCounter + 1;
                            keepChecking = 1;
                        end
                        %pause
                    end
                end
                if ~isempty(PHTelem(curElem).neighbor_right)
                    %disp('East')
                    neighborRight = PHTelem(curElem).neighbor_right(1);
                    if curLevel>PHTelem(neighborRight).level
                        [patchIndex, ~] = find(patchList==neighborRight);
                        if patchRef(patchIndex)==0
                            patchRef(patchIndex) = 1;
                            extraPatchCounter = extraPatchCounter + 1;
                            keepChecking = 1;
                        end
                    end
                end
                
                if ~isempty(PHTelem(curElem).neighbor_up_left)
                    %disp('North')
                    neighborUpLeft = PHTelem(curElem).neighbor_up_left(1);
                    if curLevel>PHTelem(neighborUpLeft).level
                        [patchIndex, ~] = find(patchList==neighborUpLeft);
                        if patchRef(patchIndex)==0
                            patchRef(patchIndex) = 1;
                            extraPatchCounter = extraPatchCounter + 1;
                            keepChecking = 1;
                        end
                        %pause
                    end
                end
                if ~isempty(PHTelem(curElem).neighbor_right_front)
                    %disp('East')
                    neighborRightFront = PHTelem(curElem).neighbor_right_front(1);
                    if curLevel>PHTelem(neighborRightFront).level
                        [patchIndex, ~] = find(patchList==neighborRightFront);
                        if patchRef(patchIndex)==0
                            patchRef(patchIndex) = 1;
                            extraPatchCounter = extraPatchCounter + 1;
                            keepChecking = 1;
                        end
                    end
                end
                
                if ~isempty(PHTelem(curElem).neighbor_up_right)
                    %disp('North')
                    neighborUpRight = PHTelem(curElem).neighbor_up_right(1);
                    if curLevel>PHTelem(neighborUpRight).level
                        [patchIndex, ~] = find(patchList==neighborUpRight);
                        if patchRef(patchIndex)==0
                            patchRef(patchIndex) = 1;
                            extraPatchCounter = extraPatchCounter + 1;
                            keepChecking = 1;
                        end
                        %pause
                    end
                end
                if ~isempty(PHTelem(curElem).neighbor_right_back)
                    %disp('East')
                    neighborRightBack = PHTelem(curElem).neighbor_right_back(1);
                    if curLevel>PHTelem(neighborRightBack).level
                        [patchIndex, ~] = find(patchList==neighborRightBack);
                        if patchRef(patchIndex)==0
                            patchRef(patchIndex) = 1;
                            extraPatchCounter = extraPatchCounter + 1;
                            keepChecking = 1;
                        end
                    end
                end
                if ~isempty(PHTelem(curElem).neighbor_up_front)
                    %disp('East')
                    neighborUpFront = PHTelem(curElem).neighbor_up_front(1);
                    if curLevel>PHTelem(neighborUpFront).level
                        [patchIndex, ~] = find(patchList==neighborUpFront);
                        if patchRef(patchIndex)==0
                            patchRef(patchIndex) = 1;
                            extraPatchCounter = extraPatchCounter + 1;
                            keepChecking = 1;
                        end
                    end
                end
                if ~isempty(PHTelem(curElem).neighbor_up_back)
                    %disp('East')
                    neighborUpBack = PHTelem(curElem).neighbor_up_back(1);
                    if curLevel>PHTelem(neighborUpBack).level
                        [patchIndex, ~] = find(patchList==neighborUpBack);
                        if patchRef(patchIndex)==0
                            patchRef(patchIndex) = 1;
                            extraPatchCounter = extraPatchCounter + 1;
                            keepChecking = 1;
                        end
                    end
                end
            end
        end
    end
    disp(['Extra patches refined: ',num2str(extraPatchCounter)]);
end

% Sort patchList by level of elements to prevent coarse->fine level errors
levelList = zeros(1,size(patchList,1));
for i=1:size(patchList,1)
    levelList(i) = PHTelem(patchList(i,1)).level;
end

[~,sortedIndex] = sort(levelList);
patchList = patchList(sortedIndex,:);
patchRef = patchRef(sortedIndex);
tempPHTelem = PHTelem;
patchCounter = 0;
for i=1:length(patchRef)
    if patchRef(i)
        curPatchElem = patchList(i,:);
        markRef(curPatchElem) = 1;
        for iCount = 1:8
            curElem = curPatchElem(iCount);
            left = tempPHTelem(curElem).neighbor_left;
            right = tempPHTelem(curElem).neighbor_right;
            down = tempPHTelem(curElem).neighbor_down;
            up = tempPHTelem(curElem).neighbor_up;
            front = tempPHTelem(curElem).neighbor_front;
            back = tempPHTelem(curElem).neighbor_back;
            upleft = tempPHTelem(curElem).neighbor_up_left;
            downleft = tempPHTelem(curElem).neighbor_down_left;
            upright = tempPHTelem(curElem).neighbor_up_right;
            downright = tempPHTelem(curElem).neighbor_down_right;
            upfront = tempPHTelem(curElem).neighbor_up_front;
            downfront = tempPHTelem(curElem).neighbor_down_front;
            upback = tempPHTelem(curElem).neighbor_up_back;
            downback = tempPHTelem(curElem).neighbor_down_back;
            leftfront = tempPHTelem(curElem).neighbor_left_front;
            rightfront = tempPHTelem(curElem).neighbor_right_front;
            leftback = tempPHTelem(curElem).neighbor_left_back;
            rightback = tempPHTelem(curElem).neighbor_right_back;
            
            neighbor = [left,right,down,up,front,back,upleft,downleft,upright,downright,upfront,downfront,upback,downback,leftfront,rightfront,leftback,rightback];
            markRef(neighbor) = 1;
        end
        [PHTelem,controlPts,dimBasis] = crossInsertIso3D(PHTelem,controlPts,curPatchElem,dimBasis,geometry);
        newPatchList(patchCounter+1,:) = numElem+1:numElem+8;
        newPatchList(patchCounter+2,:) = numElem+9:numElem+16;
        newPatchList(patchCounter+3,:) = numElem+17:numElem+24;
        newPatchList(patchCounter+4,:) = numElem+25:numElem+32;
        newPatchList(patchCounter+5,:) = numElem+33:numElem+40;
        newPatchList(patchCounter+6,:) = numElem+41:numElem+48;
        newPatchList(patchCounter+7,:) = numElem+49:numElem+56;
        newPatchList(patchCounter+8,:) = numElem+57:numElem+64;
        numElem = numElem + 64;
        patchCounter = patchCounter + 8;
    else
        newPatchList(patchCounter+1,:) = patchList(i,:);
        patchCounter = patchCounter + 1;
    end
end
patchList = newPatchList;
end