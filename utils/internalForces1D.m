function [fenerg] = internalForces1D(PhPHTelem,PhdgdxTrans,ERefElem,EPHTelem,Edisp,Mater,fenerg,geometry)
% Compute the decomposd stress

ngaussX = nGauss;
elementCounter = 0;
for indexPatch = 1:length(PhPHTelem)
    for i=1:length(PhPHTelem{indexPatch})
        if isempty(PhPHTelem{indexPatch}(i).children)
            elementCounter =  elementCounter+1;            
            % Calculate strains and stresses at integration points
            kgauss=0;
            for ii=1:geometry.ngaussX                
                    kgauss=kgauss+1;
                    EelmtNum = ERefElem{indexPatch}{i}(kgauss);
                    Enument = size(EPHTelem{indexPatch}(EelmtNum).C,1);
                    Esctrx = EPHTelem{indexPatch}(EelmtNum).nodesGlobal(1:Enument);                    
                    dispElmt = Edisp(Esctrx);                    
                    Bu = PhdgdxTrans{indexPatch}{i}(kgauss,:);                    
                    strainElmt = Bu*dispElmt;
                    stress =  Mater.E*strainElmt;
                    senerg = 0.5*stress* strainElmt;                                     
                    senerg = max(senerg,fenerg{indexPatch}{i}(kgauss));
                    fenerg{indexPatch}{i}(kgauss) = senerg;  
            end %jgaus
        end
    end
end