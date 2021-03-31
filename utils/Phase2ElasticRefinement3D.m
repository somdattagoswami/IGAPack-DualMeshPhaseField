function [shape,PhRefElem] =Phase2ElasticRefinement3D(EPHTelem,PhPHTelem,PhcontrolPts,shape,PhRefElem,EmarkRef,geometry)

numPatches = length(EPHTelem);
for indexPatch = 1:numPatches
    for i = 1:length(EmarkRef{indexPatch})
        if EmarkRef{indexPatch}(i)
            PhaseElmtNum = unique(PhRefElem{indexPatch}{i});
            if isempty(EPHTelem{indexPatch}(i).children)
                xmin = EPHTelem{indexPatch}(i).vertex(1);
                xmax = EPHTelem{indexPatch}(i).vertex(4);
                ymin = EPHTelem{indexPatch}(i).vertex(2);
                ymax = EPHTelem{indexPatch}(i).vertex(5); 
                zmin = EPHTelem{indexPatch}(i).vertex(3);
                zmax = EPHTelem{indexPatch}(i).vertex(6);  
                [shape{indexPatch}{i},~,PhRefElem{indexPatch}{i}] = transferCalculation3D(PhaseElmtNum,PhPHTelem{indexPatch},PhcontrolPts{indexPatch},geometry,xmin,xmax,ymin,ymax,zmin,zmax);
            else
                childElmt = EPHTelem{indexPatch}(i).children;
                for ichild = 1:8
                    elmtNum = childElmt(ichild);
                    
                    xmin = EPHTelem{indexPatch}(elmtNum).vertex(1);
                    xmax = EPHTelem{indexPatch}(elmtNum).vertex(4);
                    ymin = EPHTelem{indexPatch}(elmtNum).vertex(2);
                    ymax = EPHTelem{indexPatch}(elmtNum).vertex(5);
                    zmin = EPHTelem{indexPatch}(elmtNum).vertex(3);
                    zmax = EPHTelem{indexPatch}(elmtNum).vertex(6);              
                    [shape{indexPatch}{elmtNum},~,PhRefElem{indexPatch}{elmtNum}] = transferCalculation3D(PhaseElmtNum,PhPHTelem{indexPatch},PhcontrolPts{indexPatch},geometry,xmin,xmax,ymin,ymax,zmin,zmax);
                end
                shape{indexPatch}{i} = [];
                PhRefElem{indexPatch}{i} = [];
            end            
        end
    end
end
end

