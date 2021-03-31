function [ETrans] = Phase2ElasticRefinement(EPHTelem,PhPHTelem,PhcontrolPts,ETrans,EmarkRef,geometry)

p = geometry.p;
q = geometry.q;
nGauss = geometry.ngaussX;
numPatches = geometry.numPatches;

for indexPatch = 1:numPatches
    for i = 1:length(EmarkRef{indexPatch})
        if EmarkRef{indexPatch}(i)
            PhaseElmtNum = unique(ETrans.PhRefelem{indexPatch}{i});
            if isempty(EPHTelem{indexPatch}(i).children)
                xmin = EPHTelem{indexPatch}(i).vertex(1);
                xmax = EPHTelem{indexPatch}(i).vertex(3);
                ymin = EPHTelem{indexPatch}(i).vertex(2);
                ymax = EPHTelem{indexPatch}(i).vertex(4);                
                [ETrans.shape{indexPatch}{i},~,ETrans.PhRefelem{indexPatch}{i}] = transferCalculation(PhaseElmtNum,PhPHTelem{indexPatch},PhcontrolPts{indexPatch},p,q,xmin,xmax,ymin,ymax,nGauss);
            else
                childElmt = EPHTelem{indexPatch}(i).children;
                for ichild = 1:4
                    elmtNum = childElmt(ichild);
                    xmin = EPHTelem{indexPatch}(elmtNum).vertex(1);
                    xmax = EPHTelem{indexPatch}(elmtNum).vertex(3);
                    ymin = EPHTelem{indexPatch}(elmtNum).vertex(2);
                    ymax = EPHTelem{indexPatch}(elmtNum).vertex(4);                    
                    [ETrans.shape{indexPatch}{elmtNum},~,ETrans.PhRefelem{indexPatch}{elmtNum}] = transferCalculation(PhaseElmtNum,PhPHTelem{indexPatch},PhcontrolPts{indexPatch},p,q,xmin,xmax,ymin,ymax,nGauss);
                end
                ETrans.shape{indexPatch}{i} = [];
                ETrans.PhRefelem{indexPatch}{i} = [];
            end            
        end
    end
end
end

