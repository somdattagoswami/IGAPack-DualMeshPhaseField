function [PhTrans] = Elastic2Phase(PhPHTelem,EPHTelem,EcontrolPts,geometry)
% In the Script, the data from the Elastic mesh is transferred to the Phase field Mesh
% The code is executed on each element in the Phase field Mesh and the
% corresponding element in the Elastic mesh is found by seaching
% through the parent of the element in the Phase field mesh. Note: The parents
% of the elemnts in both the meshes are same.

p = geometry.p;
q = geometry.q;

ngaussX = geometry.ngaussX;
ngaussY = geometry.ngaussY;
maxRefLevel = geometry.maxERefLevel;

[Gpts,~]=GaussQuad(ngaussX);
toleq = 1e-10;
PhTrans = struct;
PhTrans.shape = cell(length(PhPHTelem),1); % Preallocating memory to store the value of shape functions at all gauss points on Phase Field Mesh
PhTrans.dgdx = cell(length(PhPHTelem),1);% Preallocating memory to store the value of first derivative at all gauss points on Phase Field Mesh
PhTrans.ERefelem = cell(length(PhPHTelem),1); % Stores the element no. from the Elastic Mesh corresponding to each gauss point of each element in the Phase Field Mesh

for indexPatch = 1:length(PhPHTelem)
    PhTrans.shape(indexPatch) = {cell(1,length(PhPHTelem{indexPatch}))};
    PhTrans.dgdx(indexPatch) = {cell(1,length(PhPHTelem{indexPatch}))};
    PhTrans.ERefelem(indexPatch) = {cell(1,length(PhPHTelem{indexPatch}))};
    for i=1:length(PhPHTelem{indexPatch})% For each element in the Phase mesh, we search for the corresponding element in the Elastic Mesh and then evaluate the basis function
        if isempty(PhPHTelem{indexPatch}(i).children)
            
            xmin = PhPHTelem{indexPatch}(i).vertex(1);
            xmax = PhPHTelem{indexPatch}(i).vertex(3);
            ymin = PhPHTelem{indexPatch}(i).vertex(2);
            ymax = PhPHTelem{indexPatch}(i).vertex(4);
            
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
            for ii=1:ngaussX
                for jj=1:ngaussY
                    kgauss=kgauss+1;
                    u_hat = Gpts(ii);
                    v_hat = Gpts(jj);
                    
                    % Mapping (u_hat,v_hat) on [-1,1]x[-1,1] to (xi,eta) on[xmin,xmax]x[ymin,ymax]
                    xi = u_hat*(xmax-xmin)/2+(xmax+xmin)/2;
                    eta = v_hat*(ymax-ymin)/2+(ymax+ymin)/2;
                    
                    
                    for iElmt = 1:size(vertex,1)
                        if (xi <= vertex(iElmt,3)+toleq) && (xi >= vertex(iElmt,1)-toleq) && (eta <= vertex(iElmt,4)+toleq) && (eta >= vertex(iElmt,2)-toleq)
                            EelmtNum = leaf(iElmt);
                            Pxmin = vertex(iElmt,1);
                            Pxmax = vertex(iElmt,3);
                            Pymin = vertex(iElmt,2);
                            Pymax = vertex(iElmt,4);
                            break;
                        end
                    end
                    uu_hat = (2*xi - Pxmin - Pxmax)/(Pxmax-Pxmin);
                    vv_hat = (2*eta - Pymin - Pymax)/(Pymax-Pymin);
                    
                    nument = size(EPHTelem{indexPatch}(EelmtNum).C,1);
                    nodes = EPHTelem{indexPatch}(EelmtNum).nodes(1:nument);
                    cpts = EcontrolPts{indexPatch}(nodes,1:2);
                    wgts = EcontrolPts{indexPatch}(nodes,3);
                    PhTrans.shape{indexPatch}{i}(kgauss,:) = phtBasis(uu_hat,vv_hat,EPHTelem{indexPatch}(EelmtNum).C,p,q);
                    [~,dR] = phtBasisIso(uu_hat,vv_hat,EPHTelem{indexPatch}(EelmtNum).C,p,q,wgts);
                    dR(1,:) = dR(1,:)*2/(Pxmax-Pxmin);
                    dR(2,:) = dR(2,:)*2/(Pymax-Pymin);
                    dxdxi = dR*cpts;
                    dR = dxdxi\dR;
                    PhTrans.dgdx{indexPatch}{i}(kgauss,:,:) = dR;
                    PhTrans.ERefelem{indexPatch}{i}(kgauss) = EelmtNum;
                end %igauss
            end %jgauss
        end
    end
end
end %endfunction