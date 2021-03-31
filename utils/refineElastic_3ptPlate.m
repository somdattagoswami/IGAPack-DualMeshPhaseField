function [EPHTelem,EcontrolPts,EmeshInfo,EmarkRef] = refineElastic_3ptPlate(EelemRef,EmeshInfo,EPHTelem,EcontrolPts,geometry)

EmarkRef = cell(1,geometry.numPatches);
p = geometry.p;
q = geometry.q;

for indexPatch = 1:geometry.numPatches
    EmarkRef(indexPatch) = {zeros(length(EPHTelem{indexPatch}),1)};
    if sum(EelemRef{indexPatch}) > 0
        [EPHTelem{indexPatch},EcontrolPts{indexPatch},EmeshInfo.dimBasis(indexPatch),EmeshInfo.numElements,EmarkRef{indexPatch}] = ...
            refineElemGradedIso(EelemRef{indexPatch},EPHTelem{indexPatch},EcontrolPts{indexPatch},p,q,EmeshInfo.dimBasis(indexPatch),EmeshInfo.numElements);
    end
end
[EPHTelem,EcontrolPts,EmeshInfo.dimBasis,EmeshInfo.numElements,EmarkRef] = checkConformingElemIsoMark(EPHTelem,EcontrolPts,EmeshInfo.dimBasis,geometry.patchBoundaries,p,q,EmeshInfo.numElements,EmarkRef);
[EPHTelem,EmeshInfo] = zipConforming(EPHTelem,geometry,EmeshInfo);

plot1 = subplot(2,2,1);
cla(plot1)
plotMesh(EPHTelem,EcontrolPts,geometry,0)
axis equal

end