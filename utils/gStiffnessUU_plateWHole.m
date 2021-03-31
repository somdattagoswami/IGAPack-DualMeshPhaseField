function [stiffkUU,PhmarkRef,quadRef] = gStiffnessUU_plateWHole(EPHTelem,EmeshInfo,Ebasis,Phbasis,ETrans,PhPHTelem,Phidisp,Mater,geometry)
% Assemble the stiffness matrix and RHS
% Uses separate meshes for elasticity and phase field

dim = geometry.dim;
nstress = geometry.nstress;
numPatches = geometry.numPatches;
ngaussX = geometry.ngaussX;
ngaussY = geometry.ngaussY;

PhmarkRef = cell(numPatches,1);
EmarkRef = cell(numPatches,1);

indexCounter = 0;
for indexPatch = 1:length(EPHTelem)
    PhmarkRef(indexPatch) = {zeros(1,length(PhPHTelem{indexPatch}))};
    EmarkRef(indexPatch) = {zeros(1,length(EPHTelem{indexPatch}))};
    for i=1:length(EPHTelem{indexPatch})
        if isempty(EPHTelem{indexPatch}(i).children)
            
            nument = size(EPHTelem{indexPatch}(i).C,1);
            sctrx = EPHTelem{indexPatch}(i).nodesGlobal(1:nument);
            localkUU = zeros(dim*nument,dim*nument);
            kgauss = 0;
            for ii=1:ngaussX
                for jj=1:ngaussY
                    kgauss = kgauss+1;
                    PhelmtNum = ETrans.PhRefelem{indexPatch}{i}(kgauss);
                    Phnument = size(PhPHTelem{indexPatch}(PhelmtNum).C,1);
                    Phsctrx = PhPHTelem{indexPatch}(PhelmtNum).nodesGlobal(1:Phnument);
                    ephi = Phidisp(Phsctrx);
                    phigp = (1-Phbasis.isComp{indexPatch}{PhelmtNum}(kgauss))*ETrans.shape{indexPatch}{i}(kgauss,:)*ephi;
                    
                    if (phigp >= geometry.threshPhi) && (PhPHTelem{indexPatch}(PhelmtNum).level < geometry.maxPhRefLevel)
                        PhmarkRef{indexPatch}(PhelmtNum)=1;
                        if (PhPHTelem{indexPatch}(PhelmtNum).level - EPHTelem{indexPatch}(i).level) >= (geometry.maxPhRefLevel - geometry.maxERefLevel)
                            EmarkRef{indexPatch}(i) = 1;
                        end
                    end
                    
                    [Bu,~,D] = strainGrad(Ebasis.dgdx{indexPatch}{i},nument,nstress,dim,kgauss,Mater.C);
                    
                    % Calculation of kUU
                    localkUU = localkUU+((1.0-phigp).^2).*(Bu'*D).*Ebasis.volume{indexPatch}{i}(kgauss);
                end %jgaus
            end %igaus
            dsctrx = reshape([2*sctrx-1;2*sctrx],1,dim*nument);
            II(indexCounter+1:indexCounter+dim^2*nument^2) = repmat(dsctrx,1,dim*nument);
            JJ(indexCounter+1:indexCounter+dim^2*nument^2) = reshape(repmat(dsctrx,dim*nument,1),1,dim^2*nument^2);
            S(indexCounter+1:indexCounter+dim^2*nument^2) = reshape(localkUU,1,dim^2*nument^2);
            indexCounter = indexCounter +dim^2*nument^2;
        end
    end
end
stiffkUU = sparse(II,JJ,S,dim*EmeshInfo.sizeBasis,dim*EmeshInfo.sizeBasis);
clear II JJ S

quadRef = cell(1,numPatches); % flag which quads to be refined
for indexPatch = 1:numPatches
    quadRef{indexPatch} = zeros(1,size(EmeshInfo.quadList{indexPatch},1));
    for iQuad = 1:size(EmeshInfo.quadList{indexPatch},1)
        curLevel = EPHTelem{indexPatch}(EmeshInfo.quadList{indexPatch}(iQuad,1)).level;
        quadRef{indexPatch}(iQuad) = sum(EmarkRef{indexPatch}(EmeshInfo.quadList{indexPatch}(iQuad,:)));
        if quadRef{indexPatch}(iQuad)>0 && curLevel < geometry.maxERefLevel
            quadRef{indexPatch}(iQuad) = 1;
        else
            quadRef{indexPatch}(iQuad) = 0;
        end
    end
end
end