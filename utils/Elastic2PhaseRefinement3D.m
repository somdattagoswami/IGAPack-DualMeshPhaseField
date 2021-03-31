function [dgdx,ERefElem] = Elastic2PhaseRefinement3D(PhPHTelem,EPHTelem,EcontrolPts,dgdx,ERefElem,PhmarkRef,geometry)

numPatches = length(PhPHTelem);
for indexPatch = 1:numPatches
    for i = 1:length(PhmarkRef{indexPatch})
        if PhmarkRef{indexPatch}(i)
            ElasElmtNum = unique(ERefElem{indexPatch}{i});
            if isempty(PhPHTelem{indexPatch}(i).children)
                xmin = PhPHTelem{indexPatch}(i).vertex(1);
                xmax = PhPHTelem{indexPatch}(i).vertex(4);
                ymin = PhPHTelem{indexPatch}(i).vertex(2);
                ymax = PhPHTelem{indexPatch}(i).vertex(5);
                zmin = PhPHTelem{indexPatch}(i).vertex(3);
                zmax = PhPHTelem{indexPatch}(i).vertex(6);
                
                [~,dgdx{indexPatch}{i},ERefElem{indexPatch}{i}] = transferCalculation3D(ElasElmtNum,EPHTelem{indexPatch},EcontrolPts{indexPatch},geometry,xmin,xmax,ymin,ymax,zmin,zmax);
            else
                childElmt = PhPHTelem{indexPatch}(i).children;
                for ichild = 1:8
                    elmtNum = childElmt(ichild);
                    
                    xmin = PhPHTelem{indexPatch}(elmtNum).vertex(1);
                    xmax = PhPHTelem{indexPatch}(elmtNum).vertex(4);
                    ymin = PhPHTelem{indexPatch}(elmtNum).vertex(2);
                    ymax = PhPHTelem{indexPatch}(elmtNum).vertex(5);
                    zmin = PhPHTelem{indexPatch}(elmtNum).vertex(3);
                    zmax = PhPHTelem{indexPatch}(elmtNum).vertex(6);          
                    [~,dgdx{indexPatch}{elmtNum},ERefElem{indexPatch}{elmtNum}] = transferCalculation3D(ElasElmtNum,EPHTelem{indexPatch},EcontrolPts{indexPatch},geometry,xmin,xmax,ymin,ymax,zmin,zmax);
                end
                dgdx{indexPatch}{i} = [];
                ERefElem{indexPatch}{i} = [];
            end            
        end
    end
end
end

