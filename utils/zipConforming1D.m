function[PHTelem,meshInfo] = zipConforming1D(PHTelem,meshInfo)
% Connects two conforming patches by changing the nodesGlobal entry
numPatches = length(PHTelem);
for patchIndex = 1:numPatches
    for elemIndex = 1:length(PHTelem{patchIndex})
        PHTelem{patchIndex}(elemIndex).nodesGlobal = PHTelem{patchIndex}(elemIndex).nodes;
    end
end

overlapCounter = 0;
meshInfo.sizeBasis = sum(meshInfo.dimBasis) - overlapCounter;

end