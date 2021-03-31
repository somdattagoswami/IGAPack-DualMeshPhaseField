function [Phidisp,fenerg] = computePhiPhi(geometry,Fract,Mater,Edisp,EBasis,PhBasis,fenerg,PhdgdxTrans,ERefElem,PhmeshInfo)

addfun = @plus;
stiffPhiPhi = sparse(PhmeshInfo.sizeBasis,PhmeshInfo.sizeBasis);
RHSPhi = zeros(PhmeshInfo.sizeBasis,1);
parfor indexPatch = 1:geometry.numPatches
    [stiffPhiPhiCell,RHSPhiCell,fenerg{indexPatch}] = gStiffnessPhiPhi3D(PhmeshInfo.numElements(indexPatch),PhBasis.elemLevel{indexPatch},...
        PhBasis.sctrxElem{indexPatch},PhmeshInfo.sizeBasis,PhdgdxTrans{indexPatch},PhBasis.dgdx{indexPatch},PhBasis.shape{indexPatch},...
        PhBasis.volume{indexPatch},Fract,Mater,geometry,fenerg{indexPatch},Edisp,ERefElem{indexPatch},EBasis.sctrxElem{indexPatch});
    stiffPhiPhi = addfun(stiffPhiPhi, stiffPhiPhiCell);
    RHSPhi = addfun(RHSPhi,RHSPhiCell);
end
Phidisp = stiffPhiPhi\RHSPhi;
end