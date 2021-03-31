function [fenerg] = internalForces(PhPHTelem,PhTrans,EPHTelem,Edisp,Mater,fenerg,geometry)
% Compute the decomposd stress

dim = geometry.dim;
nstress = geometry.nstress;
numPatches = geometry.numPatches;
ngaussX = geometry.ngaussX;
ngaussY = geometry.ngaussY;

for indexPatch = 1:numPatches
    for i=1:length(PhPHTelem{indexPatch})
        if isempty(PhPHTelem{indexPatch}(i).children)               
            % Calculate strains and stresses at integration points
            kgauss=0;
            for ii=1:ngaussX
                for jj=1:ngaussY
                    kgauss=kgauss+1;
                    EelmtNum = PhTrans.ERefelem{indexPatch}{i}(kgauss);
                    Enument = size(EPHTelem{indexPatch}(EelmtNum).C,1);
                    Esctrx = EPHTelem{indexPatch}(EelmtNum).nodesGlobal(1:Enument);
                    Edsctrx = reshape([2*Esctrx-1;2*Esctrx],1,dim*Enument);
                    dispElmt = Edisp(Edsctrx);
                    [Bu,~,~]=strainGrad(PhTrans.dgdx{indexPatch}{i},Enument,nstress,dim,kgauss,Mater.C);
                    strainElmt = Bu*dispElmt;
                    [positive_elast]=localDStress2(strainElmt(1),strainElmt(2),strainElmt(3)/2,Mater.lamda,Mater.mu);                    
                    senerg = max(positive_elast,fenerg{indexPatch}{i}(kgauss));
                    fenerg{indexPatch}{i}(kgauss) = senerg;                    
                end %igaus
            end %jgaus
        end
    end
end

