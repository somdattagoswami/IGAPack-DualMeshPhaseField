function [PhTrans] = Elastic2PhaseRefinement(PhPHTelem,EPHTelem,EcontrolPts,PhTrans,PhmarkRef,geometry)

p = geometry.p;
q = geometry.q;
nGauss = geometry.ngaussX;
numPatches = geometry.numPatches;

for indexPatch = 1:numPatches
    for i = 1:length(PhmarkRef{indexPatch})
        if PhmarkRef{indexPatch}(i)
            ElasElmtNum = unique(PhTrans.ERefelem{indexPatch}{i});
            if isempty(PhPHTelem{indexPatch}(i).children)
                xmin = PhPHTelem{indexPatch}(i).vertex(1);
                xmax = PhPHTelem{indexPatch}(i).vertex(3);
                ymin = PhPHTelem{indexPatch}(i).vertex(2);
                ymax = PhPHTelem{indexPatch}(i).vertex(4);
                
                [PhTrans.shape{indexPatch}{i},PhTrans.dgdx{indexPatch}{i},PhTrans.ERefelem{indexPatch}{i}] = transferCalculation(ElasElmtNum,EPHTelem{indexPatch},EcontrolPts{indexPatch},p,q,xmin,xmax,ymin,ymax,nGauss);
            else
                childElmt = PhPHTelem{indexPatch}(i).children;
                for ichild = 1:4
                    elmtNum = childElmt(ichild);
                    
                    xmin = PhPHTelem{indexPatch}(elmtNum).vertex(1);
                    xmax = PhPHTelem{indexPatch}(elmtNum).vertex(3);
                    ymin = PhPHTelem{indexPatch}(elmtNum).vertex(2);
                    ymax = PhPHTelem{indexPatch}(elmtNum).vertex(4);                    
                    [PhTrans.shape{indexPatch}{elmtNum},PhTrans.dgdx{indexPatch}{elmtNum},PhTrans.ERefelem{indexPatch}{elmtNum}] = transferCalculation(ElasElmtNum,EPHTelem{indexPatch},EcontrolPts{indexPatch},p,q,xmin,xmax,ymin,ymax,nGauss);
                end
                PhTrans.PhTrans.shape{indexPatch}{i} = [];
                PhTrans.dgdx{indexPatch}{i} = [];
                PhTrans.ERefelem{indexPatch}{i} = [];
            end            
        end
    end
end
end

