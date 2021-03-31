function [stiffkPhiPhi,Phi] = gStiffnessPhiPhi1D(PHTelem,PhmeshInfo,PhBasis,fenerg,geometry,Fract)
% Assembles the stiffness matrix and RHS (Galerkin method)
% Uses GIFT mapping Supports multipatches
% %Gauss points
% assemble the stiffness matrix and RHS
ngaussX = geometry.ngaussX;
Phi = zeros(PhmeshInfo.dimBasis,1);
indexCounter = 0;
for indexPatch = 1:length(PHTelem)
    for i=1:length(PHTelem{indexPatch})
        if isempty(PHTelem{indexPatch}(i).children)            
            nument = size(PHTelem{indexPatch}(i).C,1);
            sctrx = PHTelem{indexPatch}(i).nodesGlobal(1:nument);
            localkPhiPhi = zeros(nument,nument);
            localfPhi = zeros(nument,1);
            kgauss = 0;
            for ii=1:ngaussX                
                    kgauss = kgauss+1;
                    Bphi = PhBasis.dgdx{indexPatch}{i}(kgauss,:);
                    senerg = fenerg{indexPatch}{i}(kgauss);
                    ddphi = 2;
                    % Calculation of kPhiPhi                    
                    localkPhiPhi = localkPhiPhi + Fract.cenerg*Fract.constl* (Bphi'*Bphi).*PhBasis.volume{indexPatch}{i}(kgauss);
                    localkPhiPhi = localkPhiPhi + (((Fract.cenerg/Fract.constl) + ddphi*senerg))...
                        .*(PhBasis.shape{indexPatch}{i}(kgauss,:)'*PhBasis.shape{indexPatch}{i}(kgauss,:)).*PhBasis.volume{indexPatch}{i}(kgauss);   
                    localfPhi = localfPhi +ddphi*PhBasis.shape{indexPatch}{i}(kgauss,:)'*senerg*PhBasis.volume{indexPatch}{i}(kgauss);                
            end %igaus
            II(indexCounter+1:indexCounter + nument^2) = repmat(sctrx,1,nument);
            JJ(indexCounter+1:indexCounter + nument^2) = reshape(repmat(sctrx,nument,1),1,nument^2);
            S(indexCounter+1:indexCounter + nument^2) = reshape(localkPhiPhi,1,nument^2);
            Phi(sctrx) = Phi(sctrx) + localfPhi;
            indexCounter = indexCounter + nument^2;
        end
    end
end
stiffkPhiPhi = sparse(II,JJ,S,PhmeshInfo.dimBasis,PhmeshInfo.dimBasis);
end