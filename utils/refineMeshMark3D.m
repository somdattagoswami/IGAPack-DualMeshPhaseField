function [octupleList,PHTelem,controlPts,dimBasis,markRef] = refineMeshMark3D(octupleRef,octupleList,PHTelem,controlPts,geometry,dimBasis)
% Refines the patches marked by patchRef

numElem = length(PHTelem);
indexOctuple = find(octupleRef > 0);
numNewOctuples = length(indexOctuple);
newOctupleList = zeros(size(octupleList,1)+numNewOctuples,8);
markRef = zeros(length(PHTelem),1);
%sort patchList by level of elements to prevent coarse->fine level errors
levelList = zeros(1,size(octupleList,1));
for i=1:size(octupleList,1)
    levelList(i) = PHTelem(octupleList(i,1)).level;
end
    
[~,sortedIndex] = sort(levelList);
octupleList = octupleList(sortedIndex, :);
octupleRef = octupleRef(sortedIndex);
tempPHTelem = PHTelem;
octupleCounter = 0;
octupleIndexCounter = 0;
for i=1:length(octupleRef)
    if octupleRef(i)
        octupleIndexCounter = octupleIndexCounter + 1;
        curPatchElem = octupleList(i,:);
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
            markRef(neighbor) = 1;% The neighbours of the elemnent that has been refined is also marks since there is C1 continuity, so there is sharing of nodes between neighbouring elements. MarkRef would be used as a reference for the elemnts that have been changed incalculation of derivatives and also the transfer function
            
        end        
        [PHTelem,dimBasis,controlPts] = crossInsert3D(PHTelem,controlPts,curPatchElem,dimBasis,geometry);
        newOctupleList(octupleCounter+1,:) = numElem+1:numElem+8;
        newOctupleList(octupleCounter+2,:) = numElem+9:numElem+16;
        newOctupleList(octupleCounter+3,:) = numElem+17:numElem+24;
        newOctupleList(octupleCounter+4,:) = numElem+25:numElem+32;
        newOctupleList(octupleCounter+5,:) = numElem+33:numElem+40;
        newOctupleList(octupleCounter+6,:) = numElem+41:numElem+48;
        newOctupleList(octupleCounter+7,:) = numElem+49:numElem+56;
        newOctupleList(octupleCounter+8,:) = numElem+57:numElem+64;
        numElem = numElem + 64;
        octupleCounter = octupleCounter + 8;
    else
        newOctupleList(octupleCounter+1,:) = octupleList(i,:);
        octupleCounter = octupleCounter + 1;
    end                    
end
 octupleList = newOctupleList;
