function [tupleList,PHTelem,controlPts,dimBasis,markRef] = refineMeshIso1D(tupleRef,tupleList,PHTelem,controlPts,p,dimBasis)
% Refines the patches marked by patchRef

numElem = length(PHTelem);
markRef = zeros(length(PHTelem),1);
indexTuple = find(tupleRef > 0);
numNewTuples = length(indexTuple);
newTupleList = zeros(size(tupleList,1)+numNewTuples,2);

%sort patchList by level of elements to prevent coarse->fine level errors
levelList = zeros(1,size(tupleList,1));
for i=1:size(tupleList,1)
    levelList(i) = PHTelem(tupleList(i,1)).level;
end
    
[~,sortedIndex] = sort(levelList);
tupleList = tupleList(sortedIndex, :);
tupleRef = tupleRef(sortedIndex);
tempPHTelem = PHTelem;
tupleCounter = 0;
for i=1:length(tupleRef)
    if tupleRef(i)
        curTupleElem = tupleList(i,:);                
        markRef(curTupleElem) = 1;
        for iCount = 1:2
            curElem = curTupleElem(iCount);
            east = tempPHTelem(curElem).neighbor_right;
            west = tempPHTelem(curElem).neighbor_left;
            neighbor = [east,west];
            markRef(neighbor) = 1;
        end
        [PHTelem,controlPts,dimBasis] = crossInsert1D(PHTelem,controlPts,curTupleElem,dimBasis,p);
        newTupleList(tupleCounter+1,:) = numElem+1:numElem+2;
        newTupleList(tupleCounter+2,:) = numElem+3:numElem+4;        
        numElem = numElem + 4;
        tupleCounter = tupleCounter + 2;
    else
        newTupleList(tupleCounter+1,:) = tupleList(i,:);
        tupleCounter = tupleCounter + 1;
    end                    
end
 tupleList = newTupleList;
end