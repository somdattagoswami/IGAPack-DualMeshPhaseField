function [dgdx,ERefElem] = Elastic2PhaseRefinement1D(PhPHTelem,EPHTelem,EcontrolPts,dgdx,ERefElem,PhmarkRef,geometry)

p = geometry.p;
nGauss = geometry.ngaussX;

for indexPatch = 1:geometry.numPatches
    for i = 1:length(PhmarkRef{indexPatch})
        if PhmarkRef{indexPatch}(i)
            ElasElmtNum = unique(ERefElem{indexPatch}{i});
            if isempty(PhPHTelem{indexPatch}(i).children)
                xmin = PhPHTelem{indexPatch}(i).vertex(1);
                xmax = PhPHTelem{indexPatch}(i).vertex(2);
                               
                [~,dgdx{indexPatch}{i},ERefElem{indexPatch}{i}] = transferCalculation1D(ElasElmtNum,EPHTelem{indexPatch},EcontrolPts{indexPatch},p,xmin,xmax,nGauss);
            else
                childElmt = PhPHTelem{indexPatch}(i).children;
                for ichild = 1:2
                    elmtNum = childElmt(ichild);
                    
                    xmin = PhPHTelem{indexPatch}(elmtNum).vertex(1);
                    xmax = PhPHTelem{indexPatch}(elmtNum).vertex(2);            
                    [~,dgdx{indexPatch}{elmtNum},ERefElem{indexPatch}{elmtNum}] = transferCalculation1D(ElasElmtNum,EPHTelem{indexPatch},EcontrolPts{indexPatch},p,xmin,xmax,nGauss);
                end
                
                dgdx{indexPatch}{i} = [];
                ERefElem{indexPatch}{i} = [];
            end            
        end
    end
end
end

