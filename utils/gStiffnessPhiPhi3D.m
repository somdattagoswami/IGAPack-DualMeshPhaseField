function [stiffPhiPhi,RHSPhi,fenerg] = gStiffnessPhiPhi3D(PhnumberElements,PhelemLevel,PhsctrxElem,PhsizeBasis,PhdgdxTrans,Phdgdx,Phshape,Phvolume,Fract,Mater,geometry,fenerg,Edisp,ERefElem,EsctrxElem)
% Assembles the Phi part of the stiffness matrix and RHS (Galerkin method)
% Gauss points
RHSPhi = zeros(PhsizeBasis,1);
% Calculate the length of the the triplet array
indexSize = zeros(1,PhnumberElements);
for iElem = 1:PhnumberElements
    sctrx = PhsctrxElem{iElem};
    nument = length(sctrx);
    indexSize(iElem) = nument^2;
end

%initialize the triplet arrays
lengthTriplet = sum(indexSize);
II = zeros(1,lengthTriplet);
JJ = zeros(1,lengthTriplet);
S = zeros(1,lengthTriplet);

iCounter = 0;

for iElem = 1:length(PhelemLevel)
    if ~isempty(PhelemLevel{iElem})
        sctrx = PhsctrxElem{iElem};
        nument = length(sctrx);
        localkPhiPhi = zeros(nument,nument);
        localPhi = zeros(nument,1);        
        Edgdx = PhdgdxTrans{iElem};
        dgdx = Phdgdx{iElem};
        volume = Phvolume{iElem};
        shape = Phshape{iElem};
        BuPhi = zeros(geometry.dim,geometry.dim*nument);
        for ii=1:geometry.ngaussX
            for jj=1:geometry.ngaussY
                for kk = 1:geometry.ngaussZ
                    kgauss = kk+geometry.ngaussY*(jj-1)+geometry.ngaussX*geometry.ngaussY*(ii-1); 
                    EelmtNum = ERefElem{iElem}(kgauss);
                    Esctrx = EsctrxElem{EelmtNum};
                    Edsctrx = reshape([3*Esctrx-2;3*Esctrx-1;3*Esctrx],1,geometry.dim*nument);
                    dispElmt = Edisp(Edsctrx);
                    
                    % Calculation of Buphi
                    BuPhi(1,1:3:3*nument-2) = Edgdx(1,:,kgauss);
                    BuPhi(2,2:3:3*nument-1) = Edgdx(2,:,kgauss);
                    BuPhi(3,3:3:3*nument) =  Edgdx(3,:,kgauss);
                    
                    strainElmt = BuPhi*dispElmt;
                    strainPos = (sum(strainElmt) + abs(sum(strainElmt)))/2;
                    senerg = 0.5*Mater.lamda*(strainPos)^2 + Mater.mu*sum(strainPos.^2);
                    senerg = max(senerg,fenerg{iElem}(kgauss));
                    fenerg{iElem}(kgauss) = senerg;
                    
                    
                    % Calculation of kPhiPhi
                    localkPhiPhi = localkPhiPhi + Fract.cenerg*Fract.constl* (dgdx(:,:,kgauss)'*dgdx(:,:,kgauss)).*volume(kgauss);
                    localkPhiPhi = localkPhiPhi + (((Fract.cenerg/Fract.constl) + 2.0*senerg))...
                        .*(shape(:,kgauss)*shape(:,kgauss)').*volume(kgauss);
                    localPhi = localPhi +2*shape(:,kgauss)*senerg*volume(kgauss);
                end % kguass
            end %jgauss
        end %igauss
        
        RHSPhi(sctrx) = RHSPhi(sctrx) + localPhi;
        
        II(iCounter+1:iCounter + nument^2) = repmat(sctrx,1,nument);
        JJ(iCounter+1:iCounter + nument^2) = reshape(repmat(sctrx,nument,1),1,nument^2);
        S(iCounter+1:iCounter + nument^2) = reshape(localkPhiPhi,1,nument^2);
        iCounter = iCounter + nument^2;
    end
end
stiffPhiPhi = sparse(II,JJ,S,PhsizeBasis,PhsizeBasis);
end