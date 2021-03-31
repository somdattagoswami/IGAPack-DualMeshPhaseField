function [stiffkUU,solU,PhelemRef,tupleRef] = gStiffnessUU1D(EPHTelem,EmeshInfo,EBasis,EshapeTrans,PhPHTelem,solPhi,PhRefElem,Fract,Mater,geometry,Edirichlet)

% Assembles the stiffness matrix and RHS (Galerkin method)
% Uses GIFT mapping Supports multipatches
% Gauss points
% checks for compressed elements (condition 3 in Ambati et al)
% assemble the stiffness matrix and RHS
% uses separate meshes for elasticity and phase field

dim = 1;
numPatches = length(EPHTelem);
PhelemRef = cell(numPatches,1);
EelemRef = cell(numPatches,1);
rhs = zeros(dim*EmeshInfo.sizeBasis,1);
indexCounter = 0;
for indexPatch = 1:length(EPHTelem)
    PhelemRef(indexPatch) = {zeros(1,length(PhPHTelem{indexPatch}))};
    EelemRef(indexPatch) = {zeros(1,length(EPHTelem{indexPatch}))};
    for i=1:length(EPHTelem{indexPatch})
        if isempty(EPHTelem{indexPatch}(i).children)
            nument = size(EPHTelem{indexPatch}(i).C,1);
            sctrx = EPHTelem{indexPatch}(i).nodesGlobal(1:nument);
            localkUU = zeros(dim*nument,dim*nument);
            localrhs = zeros(nument,1); %local RHS vector
            kgauss = 0;
            for ii=1:geometry.ngaussX
                kgauss = kgauss+1;
                PhelmtNum = PhRefElem{indexPatch}{i}(kgauss);
                Phnument = size(PhPHTelem{indexPatch}(PhelmtNum).C,1);
                Phsctrx = PhPHTelem{indexPatch}(PhelmtNum).nodesGlobal(1:Phnument);
                ephi = solPhi(Phsctrx);
                phigp = EshapeTrans{indexPatch}{i}(kgauss,:)*ephi;
                
                if (phigp>= geometry.threshPhi) && (PhPHTelem{indexPatch}(PhelmtNum).level<geometry.maxPhRefLevel)
                    PhelemRef{indexPatch}(PhelmtNum)=1;
                    if (PhPHTelem{indexPatch}(PhelmtNum).level - EPHTelem{indexPatch}(i).level) >= (geometry.maxPhRefLevel - geometry.maxERefLevel)
                        EelemRef{indexPatch}(i) = 1;
                    end
                end
                Bu = EBasis.dgdx{indexPatch}{i}(kgauss,:);
                D = Mater.E*Bu;
                coord = EBasis.gaussCord{indexPatch}{i}(kgauss,1);
                P = sin(pi*coord);
                % Calculation of kUU
                localkUU = localkUU+((1.0-phigp).^2).*(Bu'*D).*EBasis.volume{indexPatch}{i}(kgauss);
                localrhs = localrhs + EBasis.shape{indexPatch}{i}(kgauss,:)'*P*EBasis.volume{indexPatch}{i}(kgauss);
            end %igaus
            II(indexCounter+1:indexCounter+dim^2*nument^2) = repmat(sctrx,1,dim*nument);
            JJ(indexCounter+1:indexCounter+dim^2*nument^2) = reshape(repmat(sctrx,dim*nument,1),1,dim^2*nument^2);
            S(indexCounter+1:indexCounter+dim^2*nument^2) = reshape(localkUU,1,dim^2*nument^2);
            indexCounter = indexCounter +dim^2*nument^2;
            rhs(sctrx) = rhs(sctrx) + localrhs;
        end
    end
end
stiffkUU = sparse(II,JJ,S,dim*EmeshInfo.sizeBasis,dim*EmeshInfo.sizeBasis);
clear II JJ S

tupleRef = cell(1,numPatches); %flag which quads to be refined
for indexPatch = 1:numPatches
    tupleRef{indexPatch} = zeros(1,size(EmeshInfo.tupleList{indexPatch},1));
    for iQuad = 1:size(EmeshInfo.tupleList{indexPatch},1)
        curLevel = EPHTelem{indexPatch}(EmeshInfo.tupleList{indexPatch}(iQuad,1)).level;
        tupleRef{indexPatch}(iQuad) = sum(EelemRef{indexPatch}(EmeshInfo.tupleList{indexPatch}(iQuad,:)));
        if tupleRef{indexPatch}(iQuad)>0 && curLevel < geometry.maxERefLevel
            tupleRef{indexPatch}(iQuad) = 1;
        else
            tupleRef{indexPatch}(iQuad) = 0;
        end
    end
end

disp('Imposing boundary conditions and solving.')
solU = applyBoundary1D(stiffkUU,rhs,Edirichlet.X,Edirichlet.ValX);
end