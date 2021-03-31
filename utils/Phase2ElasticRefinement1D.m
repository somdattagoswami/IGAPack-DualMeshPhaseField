function [shape,PhRefElem] =Phase2ElasticRefinement1D(EPHTelem,PhPHTelem,PhcontrolPts,shape,PhRefElem,EmarkRef,geometry)

p = geometry.p;
nGauss = geometry.ngaussX;

for indexPatch = 1:geometry.numPatches
    for i = 1:length(EmarkRef{indexPatch})
        if EmarkRef{indexPatch}(i)
            PhaseElmtNum = unique(PhRefElem{indexPatch}{i});
            if isempty(EPHTelem{indexPatch}(i).children)
                xmin = EPHTelem{indexPatch}(i).vertex(1);
                xmax = EPHTelem{indexPatch}(i).vertex(2);
                 
                [shape{indexPatch}{i},~,PhRefElem{indexPatch}{i}] = transferCalculation1D(PhaseElmtNum,PhPHTelem{indexPatch},PhcontrolPts{indexPatch},p,xmin,xmax,nGauss);
            else
                childElmt = EPHTelem{indexPatch}(i).children;
                for ichild = 1:2
                    elmtNum = childElmt(ichild);
                    xmin = EPHTelem{indexPatch}(elmtNum).vertex(1);
                    xmax = EPHTelem{indexPatch}(elmtNum).vertex(2);
                           
                    [shape{indexPatch}{elmtNum},~,PhRefElem{indexPatch}{elmtNum}] = transferCalculation1D(PhaseElmtNum,PhPHTelem{indexPatch},PhcontrolPts{indexPatch},p,xmin,xmax,nGauss);
                end
                shape{indexPatch}{i} = [];
                PhRefElem{indexPatch}{i} = [];
            end            
        end
    end
end
end

