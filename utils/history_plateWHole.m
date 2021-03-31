function [fenerg] = history_plateWHole(Basis,Fract,PHTelem,geometry)
% Computes the initial history field value based on the distance of the tip
% of the crack to other points on the domain
gaussCord = Basis.gaussCord;
nGaussX = geometry.ngaussX;
nGaussY = geometry.ngaussY;
constB = geometry.B;
numPatches = length(PHTelem);
fenerg = cell(1,numPatches);

for indexPatch = 1:length(PHTelem)
    fenerg(indexPatch) = {cell(1,length(PHTelem{indexPatch}))};
    for i=1:length(PHTelem{indexPatch})
        if isempty(PHTelem{indexPatch}(i).children)
            fenerg{indexPatch}{i} = zeros(1,nGaussX*nGaussY);
            kgauss = 0;
            for ii=1:nGaussX
                for jj=1:nGaussY
                    kgauss = kgauss + 1;
                    if (gaussCord{indexPatch}{i}(kgauss,2) > 1.0)
                        dis = sqrt((gaussCord{indexPatch}{i}(kgauss,1) - 4.0)^2+(gaussCord{indexPatch}{i}(kgauss,2) - 1.0)^2);
                        
                        if dis <= Fract.constl/2
                            fenerg{indexPatch}{i}(kgauss) = constB*Fract.cenerg*(1-dis/(Fract.constl/2))/(2*Fract.constl);
%                             hold on
%                             plot(gaussCord{indexPatch}{i}(kgauss,1), gaussCord{indexPatch}{i}(kgauss,2), '.r')
                        end
                    elseif (gaussCord{indexPatch}{i}(kgauss,2) <= 1.0)
                        dis = abs((gaussCord{indexPatch}{i}(kgauss,1) - 4.0));
                        if dis <= Fract.constl/2
                            fenerg{indexPatch}{i}(kgauss) = constB*Fract.cenerg*(1-dis/(Fract.constl/2))/(2*Fract.constl);
%                             hold on
%                             plot(gaussCord{indexPatch}{i}(kgauss,1), gaussCord{indexPatch}{i}(kgauss,2), '.r')
                        end
                    end
                end
            end
        end
    end
end
end