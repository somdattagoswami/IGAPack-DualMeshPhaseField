function [PHTelem,controlPts,dimBasis,fieldData,numberElements,markRef] = refineElemProjGradedIso1D(elemRef,PHTelem,controlPts,p,dimBasis,fieldData,numberElements)
% Refines the elements marked by elemRef also refine the neighbor patches if the level_neighbor<level_refined_quad
% Project the coefficients included in the field data

numElem = length(PHTelem);
markRef = zeros(length(PHTelem),1); % This array has a marking of 0 for elemnts that are not refined the ones that are not neighbouring to the elemnets that are refined.
% Check the refinement level for the neighbor quad and mark it for refinement if needed

% Sort elements by level to prevent coarse->fine level errors
% We must refine the coarsest marked elements first
levelList = zeros(1,numElem);
for i=1:numElem
    levelList(i) = PHTelem(i).level;
end

[~,sortedIndex] = sort(levelList);
elemList = 1:length(elemRef);
elemList = elemList(sortedIndex);
elemRef = elemRef(sortedIndex);
tempPHTelem = PHTelem;
for i=1:length(elemRef)
    if elemRef(i)
        curElem = elemList(i);
        if isempty(PHTelem(curElem).children)
            markRef(curElem) = 1;% The value of markRef for this element that is refined is marked 1
            east = tempPHTelem(curElem).neighbor_right;
            west = tempPHTelem(curElem).neighbor_left;
            neighbor = [east,west];
            markRef(neighbor) = 1; % The neighbours of the elemnent that has been refined is also marks since there is C1 continuity, so there is sharing of nodes between neighbouring elements. MarkRef would be used as a reference for the elemnts that have been changed incalculation of derivatives and also the transfer function
            [PHTelem, controlPts,dimBasis,fieldData] = crossInsertProjIso1D(PHTelem, controlPts,curElem,dimBasis,p,fieldData);
            % Three new elements created (1 deactivated parent and 4 children)
            numberElements = numberElements + 1;            
        else
            warning('Something wrong')
            pause
        end
        
    end
end
