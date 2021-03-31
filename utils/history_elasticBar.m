function [fenerg] = history_elasticBar(Basis,Fract,PHTelem,geometry)
% Computes the initial history field value based on the distance of the tip
% of the crack to other points on the domain
gaussCord = Basis.gaussCord;
ngaussX = geometry.ngaussX;
constB = geometry.B;
fenerg = cell(1,geometry.numPatches);

for indexPatch = 1:length(PHTelem)
    fenerg(indexPatch) = {cell(1,length(PHTelem{indexPatch}))};
    for i=1:length(PHTelem{indexPatch})
        if isempty(PHTelem{indexPatch}(i).children)
            fenerg{indexPatch}{i} = zeros(1,ngaussX);
            kgauss = 0;
            for ii=1:ngaussX
                kgauss = kgauss + 1;
                dis = abs((gaussCord{indexPatch}{i}(kgauss,1)-0));
                if dis <= Fract.constl/2
                    fenerg{indexPatch}{i}(kgauss) = constB*Fract.cenerg*(1-dis/(Fract.constl/2))/(2*Fract.constl);
                end
            end
        end
    end
end
end









