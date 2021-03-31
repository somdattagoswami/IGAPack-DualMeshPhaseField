function [dgdx,ERefelem] = Elastic2Phase1D(PhPHTelem,EPHTelem,EcontrolPts,geometry)
% In the Script, the data from the Elastic mesh is transferred to the Phase field Mesh
% The code is executed on each element in the Phase field Mesh and the
% corresponding element in the Elastic mesh is found by seaching
% through the parent of the element in the Phase field mesh. Note: The parents
% of the elemnts in both the meshes are same.
maxRefLevel = geometry.maxERefLevel;
p = geometry.p;
nGauss = geometry.ngaussX;
[Gpts,~]=GaussQuad(nGauss);
toleq = 1e-10;

dgdx = cell(length(PhPHTelem),1);% Preallocating memory to store the value of first derivative at all gauss points on Phase Field Mesh
ERefelem = cell(length(PhPHTelem),1); % Stores the element no. from the Elastic Mesh corresponding to each gauss point of each element in the Phase Field Mesh

for indexPatch = 1:length(PhPHTelem)
    dgdx(indexPatch) = {cell(1,length(PhPHTelem{indexPatch}),1)};
    ERefelem(indexPatch) = {cell(1,length(PhPHTelem{indexPatch}),1)};
    for i=1:length(PhPHTelem{indexPatch})% For each element in the Phase mesh, we search for the corresponding element in the Elastic Mesh and then evaluate the basis function
        if isempty(PhPHTelem{indexPatch}(i).children)
            
            umin = PhPHTelem{indexPatch}(i).vertex(1);
            umax = PhPHTelem{indexPatch}(i).vertex(2);
            
            indexParent = i; % For each element in the Phase mesh, we search for the corresponding element in the Elastic Mesh and then evaluate the basis function
            % Parents are the elements at elevl 1 and they are the same for both the meshes
            while PhPHTelem{indexPatch}(indexParent).level >1
                indexParent = PhPHTelem{indexPatch}(indexParent).parent; % The Parent index in the Parent mesh
            end
            leafcount = 0;
            leaf = [];
            branchNew = [];
            if isempty(EPHTelem{indexPatch}(indexParent).children)
                leaf = indexParent;
                vertex = EPHTelem{indexPatch}(indexParent).vertex;
            else
                branch = EPHTelem{indexPatch}(indexParent).children;
                for ilevel = 2 : maxRefLevel
                    branchCount = 0;
                    for  ibranch = 1 : length(branch)
                        if isempty(EPHTelem{indexPatch}(branch(ibranch)).children)
                            leafcount = leafcount + 1;
                            leaf(leafcount,1) = branch(ibranch);
                            vertex(leafcount,:) = EPHTelem{indexPatch}(branch(ibranch)).vertex;
                        else
                            branchCount = branchCount + 1;
                            branchNew(branchCount,:) = EPHTelem{indexPatch}(branch(ibranch)).children;
                        end
                    end
                    if isempty(branchNew)
                        break;
                    else
                        branch = reshape(branchNew',1,[]);
                        branchNew = [];
                    end
                end
            end
            
            kgauss=0;
            for ii=1:nGauss
                kgauss=kgauss+1;
                u_hat = Gpts(ii);  
                % Mapping u_hat on [-1,1] to xi on[xmin,xmax]
                xi = u_hat*(umax-umin)/2+(umax+umin)/2;                
                
                
                for iElmt = 1:size(vertex,1)
                    if (xi <= vertex(iElmt,2)+toleq) && (xi >= vertex(iElmt,1)-toleq)
                        EelmtNum = leaf(iElmt);
                        Pxmin = vertex(iElmt,1);
                        Pxmax = vertex(iElmt,2);                        
                        break;
                    end
                end
                uu_hat = (2*xi - Pxmin - Pxmax)/(Pxmax-Pxmin);                
                
                nument = size(EPHTelem{indexPatch}(EelmtNum).C,1);
                nodes = EPHTelem{indexPatch}(EelmtNum).nodes(1:nument);
                cpts = EcontrolPts{indexPatch}(nodes,1);
                wgts = EcontrolPts{indexPatch}(nodes,2);
                [~,dR] = phtBasisIso1D(uu_hat,EPHTelem{indexPatch}(EelmtNum).C,p,wgts);
                dR(1,:) = dR(1,:)*2/(Pxmax-Pxmin);
                dxdxi = dR*cpts;
                dR = dxdxi\dR;
                dgdx{indexPatch}{i}(kgauss,1,:) = dR;
                ERefelem{indexPatch}{i}(kgauss) = EelmtNum;
            end %igaus
        end
    end
end
end %endfunction