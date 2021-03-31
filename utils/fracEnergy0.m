function[fenerg] = fracEnergy0(PHTelem,geometry)

mgauss = geometry.ngaussX*geometry.ngaussY*geometry.ngaussZ;

fenerg = cell(1,geometry.numPatches); % Initilaizing the Elastic Strain Energy Matrix
for indexPatch = 1:length(PHTelem)
    fenerg(indexPatch) = {cell(length(PHTelem{indexPatch}),1)};
    for i=1:length(PHTelem{indexPatch})
        if isempty(PHTelem{indexPatch}(i).children)
            fenerg{indexPatch}{i} = zeros(1,mgauss);
        end
    end
end

