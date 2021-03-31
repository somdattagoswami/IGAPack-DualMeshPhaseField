function [stiffkPhiPhi,Phi] = gStiffnessPhiPhi(PHTelem,PhmeshInfo,geometry,Fract,Mater,PhBasis,fenerg)
% Assembles the stiffness matrix and RHS (Galerkin method)

dim = geometry.dim;
nstress = geometry.nstress;
numPatches = geometry.numPatches;
ngaussX = geometry.ngaussX;
ngaussY = geometry.ngaussY;

Phi = zeros(PhmeshInfo.sizeBasis,1);
indexCounter = 0;

for indexPatch = 1:numPatches
    for i=1:length(PHTelem{indexPatch})
        if isempty(PHTelem{indexPatch}(i).children)            
            nument = size(PHTelem{indexPatch}(i).C,1);
            sctrx = PHTelem{indexPatch}(i).nodesGlobal(1:nument);
            localkPhiPhi = zeros(nument,nument);
            localfPhi = zeros(nument,1);
            kgauss = 0;
            for ii=1:ngaussX
                for jj=1:ngaussY
                    kgauss = kgauss+1;  
                    [~,Bphi,~]=strainGrad(PhBasis.dgdx{indexPatch}{i},nument,nstress,dim,kgauss,Mater.C);
                    
                    senerg = fenerg{indexPatch}{i}(kgauss);
                  
                    % Calculation of kPhiPhi                    
                    localkPhiPhi = localkPhiPhi + Fract.cenerg*Fract.constl* (Bphi'*Bphi).*PhBasis.volume{indexPatch}{i}(kgauss);
                    localkPhiPhi = localkPhiPhi + (((Fract.cenerg/Fract.constl) + 2.0*senerg))...
                        .*(PhBasis.shape{indexPatch}{i}(kgauss,:)'*PhBasis.shape{indexPatch}{i}(kgauss,:)).*PhBasis.volume{indexPatch}{i}(kgauss);   
                    localfPhi = localfPhi +2*PhBasis.shape{indexPatch}{i}(kgauss,:)'*senerg*PhBasis.volume{indexPatch}{i}(kgauss);
                end %jgaus
            end %igaus
            II(indexCounter+1:indexCounter + nument^2) = repmat(sctrx,1,nument);
            JJ(indexCounter+1:indexCounter + nument^2) = reshape(repmat(sctrx,nument,1),1,nument^2);
            S(indexCounter+1:indexCounter + nument^2) = reshape(localkPhiPhi,1,nument^2);
            Phi(sctrx) = Phi(sctrx) + localfPhi;
            indexCounter = indexCounter + nument^2;
        end
    end
end
stiffkPhiPhi = sparse(II,JJ,S,PhmeshInfo.sizeBasis,PhmeshInfo.sizeBasis);
end

