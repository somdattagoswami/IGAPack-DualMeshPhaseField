function [EPHTelem,EcontrolPts,EmeshInfo,EmarkRef] = refineElastic_plateWHole(quadRef,EmeshInfo,EPHTelem,EcontrolPts,geometry)

EmarkRef = cell(1,geometry.numPatches);
p = geometry.p;
q = geometry.q;

for indexPatch = 1:geometry.numPatches
    EmarkRef(indexPatch) = {zeros(length(EPHTelem{indexPatch}),1)};
    if sum(quadRef{indexPatch}) > 0
        [EmeshInfo.quadList{indexPatch},EPHTelem{indexPatch},EcontrolPts{indexPatch},EmeshInfo.dimBasis(indexPatch),EmarkRef{indexPatch}] = ...
            refineMeshGradedIso(quadRef{indexPatch},EmeshInfo.quadList{indexPatch},EPHTelem{indexPatch},EcontrolPts{indexPatch},p,q,EmeshInfo.dimBasis(indexPatch));
    end
end
[EPHTelem,EcontrolPts,EmeshInfo.dimBasis,EmeshInfo.quadList,EmarkRef] = checkConformingMark_plateWHole(EPHTelem,EcontrolPts,EmeshInfo.dimBasis,geometry.patchBoundaries,p,q,EmeshInfo.quadList,EmarkRef);
[EPHTelem,EmeshInfo] = zipConforming(EPHTelem,geometry,EmeshInfo);

plot1 = subplot(2,2,1);
cla(plot1)
plotMesh(EPHTelem,EcontrolPts,geometry,0)
axis equal

end