function [shape,PhRefelem] = Phase2Elastic1D(EPHTelem,PhPHTelem,geometry)
% In the Script the data from the Phase Field mesh is transferred to the Elastic Mesh
% The code is executed on each element in the elastic Mesh and the
% corresponding element in the Phase field mesh is found by seaching
% through the parent of the element in the Elastic mesh. Note: The parent
% of the elemnts in both the meshes are same.
maxRefLevel = geometry.maxPhRefLevel;
p = geometry.p;
nGauss = geometry.ngaussX;

[Gpts,~]=GaussQuad(nGauss);
toleq = 1e-10;

shape = cell(length(EPHTelem),1);% Preallocating memory to store the value of shape function at all gauss points of Elastic Mesh
PhRefelem = cell(length(EPHTelem),1);% Stores the element no. from the Phase Field Mesh corresponding to each gauss point of each element in the Elastic Mesh

for indexPatch = 1:length(EPHTelem)
    
    shape(indexPatch) = {cell(length(EPHTelem{indexPatch}),1)};
    PhRefelem(indexPatch) = {cell(length(EPHTelem{indexPatch}),1)};
    
    for i=1:length(EPHTelem{indexPatch})% For each element in the Elastic mesh, we search for the corresponding element in the Phase Field Mesh and then evaluate the basis function
        if isempty(EPHTelem{indexPatch}(i).children)
            
            umin = EPHTelem{indexPatch}(i).vertex(1);
            umax = EPHTelem{indexPatch}(i).vertex(2);
                        
            indexParent = i; % For element number 'i': we need to compute element index of its parent whose level is 1. We dont go below 1 as upto Level 1, both the meshes are same.
            % Intially the Parent is assumed same as the element number. If the level of the parent index is 1, then this loop is not evaluated.
            % If not, this loop is evaluated over level, to reach Level1
            while EPHTelem{indexPatch}(indexParent).level >1
                indexParent = EPHTelem{indexPatch}(indexParent).parent; % The indexParent is updated for each level by assigning to to the corresponding parent index
            end
            % The element number denoted by indexParent is same for
            % both the meshes.
            leafCount = 0; % Counts the elements without children in the Phase Field mesh for the parent which has been found.
            leaf = [];% Stores the element numbers of elements without children in the Phase Field mesh for the parent which has been found.
            branchNew = [];% Stores the element number of the elemnts at each (level+1) which have children
            if isempty(PhPHTelem{indexPatch}(indexParent).children)% If the parent has no children, then leaf is the indexParent itself and vertices are the vertices of the element. This is executed only when the refinement level of the element is 1 and no further refinement has taken place on the element.
                leaf = indexParent;
                vertex = PhPHTelem{indexPatch}(indexParent).vertex;
            else
                branch = PhPHTelem{indexPatch}(indexParent).children;% This contains all the children of the element at each level
                for ilevel = 2 : maxRefLevel % This loop is evaluated from 2 as the level for index Parent is 1 and already branch contains the children. The loop is evaluated till the max refinement level of the Phase Field mesh
                    branchCount = 0;
                    for  ibranch = 1 : length(branch)% The length of branch and branch is updated at each level
                        if isempty(PhPHTelem{indexPatch}(branch(ibranch)).children)% This loop is evaluated if a particular branch element in branch has no children. Then its vertices and the element number is stored at it has reached its own max refinement
                            leafCount = leafCount + 1;
                            leaf(leafCount,1) = branch(ibranch);
                            vertex(leafCount,:) = PhPHTelem{indexPatch}(branch(ibranch)).vertex;
                        else % The 'else' section is evaluated when a particular branch element has children and then the children are stored in branchNew
                            branchCount = branchCount + 1;
                            branchNew(branchCount,:) = PhPHTelem{indexPatch}(branch(ibranch)).children;
                        end
                    end
                    if isempty(branchNew) % Branch new would be empty when all the children have been noted for all the elements within indexParent
                        break;
                    else
                        branch = reshape(branchNew',1,[]);% The branch is updated at each level by branch new which contains the element number that has children and have not reached its maximum level of Refinement.
                        branchNew = [];
                    end
                end
            end
            
            kgauss = 0;
            for ii=1:nGauss                
                    kgauss = kgauss+1;
                    u_hat = Gpts(ii);                    
                    
                    % Mapping (u_hat,v_hat) on [-1,1]x[-1,1] to (xi,eta) on[xmin,xmax]x[ymin,ymax]
                    xi = u_hat*(umax-umin)/2+(umax+umin)/2;
                    
                    for iElmt = 1:size(vertex,1)% The [xi,eta] computed have to be searched within the elements stored in leaf
                        if (xi <= vertex(iElmt,2)+toleq) && (xi >= vertex(iElmt,1)-toleq)% Since the comparison is between floating point numbers, eps has been used
                            PhelmtNum = leaf(iElmt);
                            Pxmin = vertex(iElmt,1);
                            Pxmax = vertex(iElmt,2);                            
                            break;
                        end
                    end
                    % Mapping (xi,eta) on [Pxmin,Pxmax]x[Pymin,Pymax] to (uu_hat,vv_hat) on [-1,1]
                    uu_hat = (2*xi - Pxmin - Pxmax)/(Pxmax-Pxmin);                    
                    shape{indexPatch}{i}(kgauss,:) = phtBasis1D(uu_hat,PhPHTelem{indexPatch}(PhelmtNum).C,p);
                    PhRefelem{indexPatch}{i}(kgauss) = PhelmtNum;% The element number from the Phase Field mesh corresponding to that on the elastic mesh is stored
                
            end %igaus
        end
    end
end
end