function [fenerg] = history_cubeWHole(Basis,Fract,PHTelem,geometry)
% Computes the initial history field value based on the distance of the tip
% of the crack to other points on the domain
nGaussX = geometry.p+1;
nGaussY = geometry.q+1;
nGaussZ = geometry.r+1;
constB = geometry.B;
fenerg = cell(geometry.numPatches,1); % Initilaizing the Elastic Strain Energy Matrix

for indexPatch = 1:length(PHTelem)
    fenerg(indexPatch) = {cell(1,length(PHTelem{indexPatch}))};
    for i=1:length(PHTelem{indexPatch})
        if isempty(PHTelem{indexPatch}(i).children)            
            fenerg{indexPatch}{i} = zeros(1,nGaussX*nGaussY*nGaussZ); 

            kgauss = 0;
            for ii = 1:nGaussX
                for jj = 1:nGaussY
                    for kk = 1:nGaussZ
                        kgauss = kgauss + 1;                        
                        dis = sqrt((Basis.gaussCord{indexPatch}{i}(1,kgauss)-0.0)^2+(Basis.gaussCord{indexPatch}{i}(2,kgauss)-0.0)^2+(Basis.gaussCord{indexPatch}{i}(3,kgauss)-1.0)^2);
                        if dis <= Fract.constl/2
                            fenerg{indexPatch}{i}(kgauss) = constB*Fract.cenerg*(1-dis/(Fract.constl/2))/(2*Fract.constl);
                        end
                        
                    end
                end
            end            
        end
    end
end
end