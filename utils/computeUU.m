function [stiffUU,PhmarkRef,EoctupleRef,Edisp] = computeUU(geometry,EmeshInfo,Fract,Mater,Edirichlet,tfacto,Edisp,Phidisp,EshapeTrans,PhRefElem,EBasis,PhBasis)

addfun = @plus;
stiffUU = sparse(geometry.dim*EmeshInfo.sizeBasis,geometry.dim*EmeshInfo.sizeBasis);
PhmarkRef = cell(1,geometry.numPatches);
EoctupleRef = cell(1,geometry.numPatches);
parfor indexPatch = 1:geometry.numPatches
    
    [stiffUUCell,PhmarkRef{indexPatch},EoctupleRef{indexPatch}] = gStiffnessUU3D(EmeshInfo.numElements(indexPatch),...
        EmeshInfo.octupleList{indexPatch},EBasis.elemLevel{indexPatch},EBasis.sctrxElem{indexPatch},EmeshInfo.sizeBasis,...
        EBasis.dgdx{indexPatch},EBasis.volume{indexPatch},EshapeTrans{indexPatch},PhBasis.elemLevel{indexPatch},...
        PhBasis.sctrxElem{indexPatch},Phidisp,Fract,Mater,geometry,PhRefElem{indexPatch});
    stiffUU = addfun(stiffUU,stiffUUCell);
    
end
clear stiffUUCell
toc

disp('Imposing boundary conditions and solving...')
tic
Edisp = applyBoundary3D(Edirichlet,tfacto,stiffUU,Edisp);
end