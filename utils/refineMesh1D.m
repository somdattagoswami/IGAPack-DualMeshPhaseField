function [PhmarkRef,EmarkRef,PhPHTelem,PhcontrolPts,Phidisp,EPHTelem,EcontrolPts,EmeshInfo,PhmeshInfo] = refineMesh1D(elemRef,tupleRef,PhPHTelem,EPHTelem,EcontrolPts,PhcontrolPts,PhmeshInfo,EmeshInfo,geometry,Phidisp)

PhmarkRef = cell(1,geometry.numPatches);
EmarkRef = cell(1,geometry.numPatches);
p = geometry.p;
[solPhiPatch] = transferFieldGlob2Loc(PhPHTelem,PhmeshInfo,Phidisp);

if (sum(elemRef{1})>0)
    
    % Refine and Update the mesh
    
    disp(['In patch ', num2str(1), ' refining ', num2str(sum(elemRef{1})), ' elements...'])
    [PhPHTelem{1},PhcontrolPts{1},PhdimBasis(1),solPhiPatch{1},PhmeshInfo.NumElements,PhmarkRef{1}] = ...
        refineElemProjGradedIso1D(elemRef{1},PhPHTelem{1},PhcontrolPts{1},geometry.p,PhmeshInfo.dimBasis(1),...
        solPhiPatch{1},PhmeshInfo.numElements);
    
end
PhmeshInfo.dimBasis = PhdimBasis;

[PhPHTelem,PhmeshInfo] = zipConforming1D(PhPHTelem,PhmeshInfo);
Phidisp = transferFieldLoc2Glob(PhPHTelem,PhmeshInfo,solPhiPatch);

if sum([tupleRef{:}])
    [EtupleList{1},EPHTelem{1},EcontrolPts{1},EdimBasis(1),EmarkRef{1}] = refineMeshIso1D(tupleRef{1},EmeshInfo.tupleList{1},...
        EPHTelem{1},EcontrolPts{1},p,EmeshInfo.dimBasis(1));
    EmeshInfo.dimBasis = EdimBasis;
    EmeshInfo.tupleList = EtupleList;
    [EPHTelem,EmeshInfo] = zipConforming1D(EPHTelem,EmeshInfo);
    
end
end