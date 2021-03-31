function [stiff,PhmarkRef,octupleRef] = gStiffnessUU3D(EnumberElements,EoctupleList,EelemLevel,EsctrxElem,EsizeBasis,Edgdx,Evolume,EshapeTrans,PhelemLevel,PhsctrxElem,Phidisp,Fract,Mater,geometry,PhRefElem)
% Assembles the stiffness matrix and RHS (Galerkin method)
% Gauss points

dim = geometry.dim;
PhmarkRef = zeros(1,length(PhelemLevel));
EmarkRef = zeros(1,length(EelemLevel));
iCounter = 0;

%calculate the length of the the triplet array
indexSize = zeros(1,EnumberElements);
for iElem = 1:EnumberElements
    sctrx = EsctrxElem{iElem};
    nument = length(sctrx);
    indexSize(iElem) = dim^2*nument^2;
end

%initialize the triplet arrays
lengthTriplet = sum(indexSize);
II = zeros(1,lengthTriplet);
JJ = zeros(1,lengthTriplet);
S = zeros(1,lengthTriplet);

for iElem = 1:length(EelemLevel)
    if ~isempty(EelemLevel{iElem})
        sctrx = EsctrxElem{iElem};
        nument = length(sctrx);
        localkUU = zeros(dim*nument,dim*nument);
        dgdx = Edgdx{iElem};
        volume = Evolume{iElem};        
        for ii=1:geometry.ngaussX
            for jj=1:geometry.ngaussY
                for kk = 1:geometry.ngaussZ
                    Bu = zeros(geometry.nstress,geometry.dim*nument);
                    kgauss = kk+geometry.ngaussY*(jj-1)+geometry.ngaussX*geometry.ngaussY*(ii-1);                    
                    PhelmtNum = PhRefElem{iElem}(kgauss);
                    Phsctrx = PhsctrxElem{PhelmtNum};
                    ephi = Phidisp(Phsctrx);
                    phigp = EshapeTrans{iElem}(:,kgauss)'*ephi;
                    
                    if (phigp>= geometry.threshPhi) && (PhelemLevel{PhelmtNum} < geometry.maxPhRefLevel)
                        PhmarkRef(PhelmtNum)=1;
                        if (PhelemLevel{PhelmtNum} - EelemLevel{iElem}) >= (geometry.maxPhRefLevel - geometry.maxERefLevel)
                            EmarkRef(iElem) = 1;
                        end
                    end
                    
                    Bu(1,1:3:3*nument-2) = dgdx(1,:,kgauss);
                    Bu(2,2:3:3*nument-1) = dgdx(2,:,kgauss);
                    Bu(3,3:3:3*nument) =  dgdx(3,:,kgauss);
                    
                    Bu(4,1:3:3*nument-2) = dgdx(2,:,kgauss);
                    Bu(4,2:3:3*nument-1) = dgdx(1,:,kgauss);
                    
                    Bu(5,2:3:3*nument-1) = dgdx(3,:,kgauss);
                    Bu(5,3:3:3*nument) = dgdx(2,:,kgauss);
                    
                    Bu(6,1:3:3*nument-2) = dgdx(3,:,kgauss);
                    Bu(6,3:3:3*nument) = dgdx(1,:,kgauss);
                    
                    
                    % Calculation of kUU
                    localkUU=localkUU+((1.0-phigp).^2).*(Bu'*Mater.C*Bu).*volume(kgauss);
                end
            end %jgaus
        end %igaus
        dsctrx = reshape([3*sctrx-2;3*sctrx-1;3*sctrx],1,dim*nument);
        II(iCounter+1:iCounter + dim^2*nument^2) = repmat(dsctrx,1,dim*nument);
        JJ(iCounter+1:iCounter + dim^2*nument^2) = reshape(repmat(dsctrx,dim*nument,1),1,dim^2*nument^2);
        S(iCounter+1:iCounter + dim^2*nument^2) = reshape(localkUU,1,dim^2*nument^2);
        iCounter = iCounter + dim^2*nument^2;
    end
end
stiff = sparse(II,JJ,S,dim*EsizeBasis,dim*EsizeBasis);

octupleRef = zeros(1,size(EoctupleList,1));% Flag which quads to be refined
if sum(EoctupleList(:))
    for iOct = 1:size(EoctupleList,1)
        curLevel = EelemLevel{EoctupleList(iOct,1)};
        octupleRef(iOct) = sum(EmarkRef(EoctupleList(iOct,:)));
        if octupleRef(iOct)> 0 && curLevel < geometry.maxERefLevel
            octupleRef(iOct) = 1;
        else
            octupleRef(iOct) = 0;
        end
    end
end
end