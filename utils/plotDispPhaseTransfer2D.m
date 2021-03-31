function plotDispPhaseTransfer2D(PhPHTelem,EPHTelem,PhcontrolPts,EcontrolPts,Phidisp,Edisp,PhmeshInfo,geometry,Mater)

% Evaluats the displacements and phase field values at the element vertices
% Displaying the displacements
% uses the PHT element structure
% for multipatch geometries
PhNumElements = PhmeshInfo.numElements;
p = geometry.p;
q = geometry.q;
vIEN = zeros(PhNumElements,4);
physcoord = zeros(4*PhNumElements,2);
dispcoord = zeros(4*PhNumElements,2);
wcoord = zeros(4*PhNumElements,1);
sigmacoord = zeros(4*PhNumElements,1);
dim = 2;

fudge = 0;
nx = 2;
ny = 2;
px = linspace(-1+fudge,1-fudge, nx);
py = linspace(-1+fudge,1-fudge, ny);

%1D bernstein polynomials evaluated at the Gauss points on the master element
[B_u,~] = bernstein_basis(px,p);
[B_v,~] = bernstein_basis(py,q);

R = zeros(nx, ny, (p+1)*(q+1));

% The derivatives of the 2D Bernstein polynomials at Gauss points on the master element
basisCounter = 0;
for j=1:q+1
    for i=1:p+1
        basisCounter = basisCounter + 1;
        R(:,:,basisCounter) = B_u(:,i)*B_v(:,j)';
    end
end


elementCounter = 0;
for indexPatch = 1:length(PhPHTelem)
    for i=1:length(PhPHTelem{indexPatch})
        if isempty(PhPHTelem{indexPatch}(i).children)
            elementCounter =  elementCounter+1;
            
            vIEN(elementCounter,:) = [(elementCounter-1)*4+1:(elementCounter-1)*4+4];
            ximin = PhPHTelem{indexPatch}(i).vertex(1);
            ximax = PhPHTelem{indexPatch}(i).vertex(3);
            etamin = PhPHTelem{indexPatch}(i).vertex(2);
            etamax = PhPHTelem{indexPatch}(i).vertex(4);
            
            coordt = cell(ny,nx);
            dispmatx = zeros(ny,nx);
            dispmaty = zeros(ny,nx);
            stressvect = cell(ny,nx);
            stressVM = cell(ny,nx);
            wmatx = zeros(ny,nx);
            
            nument = size(PhPHTelem{indexPatch}(i).C,1);
            nodes = PhPHTelem{indexPatch}(i).nodes(1:nument);
            sctrx = PhPHTelem{indexPatch}(i).nodesGlobal(1:nument);
            cpts = PhcontrolPts{indexPatch}(nodes,1:2);
            wgts = PhcontrolPts{indexPatch}(nodes,3);
            
            indexParent = i;
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
                for ilevel = 2 : geometry.maxERefLevel
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
            
            for ii=1:nx
                for jj=1:ny
                    % To get the values of strain and displacement from the Elastic Mesh to the Phase Field Mesh
                    u_hat = px(ii);
                    v_hat = py(jj);
                    xi = u_hat*(ximax-ximin)/2+(ximax+ximin)/2;
                    eta = v_hat*(etamax-etamin)/2+(etamax+etamin)/2;
                    
                    for iElmt = 1:size(vertex,1)
                        if (xi <= vertex(iElmt,3)+eps) && (xi >= vertex(iElmt,1)-eps) && (eta <= vertex(iElmt,4)+eps) && (eta >= vertex(iElmt,2)-eps)
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
                    
                    Enument = size(EPHTelem{indexPatch}(EelmtNum).C,1);
                    nodes = EPHTelem{indexPatch}(EelmtNum).nodes(1:Enument);
                    Ecpts = EcontrolPts{indexPatch}(nodes,1:2);
                    Ewgts = EcontrolPts{indexPatch}(nodes, 3);
                    Esctrx = EPHTelem{indexPatch}(EelmtNum).nodesGlobal(1:Enument);
                    Edsctrx = reshape([2*Esctrx-1;2*Esctrx],1,dim*Enument);
                    ER = phtBasis(uu_hat,vv_hat,EPHTelem{indexPatch}(EelmtNum).C,p,q);
                    [~,EdR] = phtBasisIso(uu_hat,vv_hat,EPHTelem{indexPatch}(EelmtNum).C,p,q,Ewgts);
                    EdR(1,:) = EdR(1,:)*2/(Pxmax-Pxmin);
                    EdR(2,:) = EdR(2,:)*2/(Pymax-Pymin);
                    Edxdxi = EdR*Ecpts;
                    EdR = Edxdxi\EdR;
                    Bu = zeros(3,dim*Enument);
                    
                    for inode=1:Enument
                        % Calculation of Bu
                        Bu(1,2*inode-1) = EdR(1,inode);
                        Bu(2,2*inode) = EdR(2,inode);
                        Bu(3,2*inode-1) = EdR(2,inode);
                        Bu(3,2*inode) = EdR(1,inode);
                    end                    
                    
                    % Calculate tdisp values
                    dispmatx(jj,ii) = dispmatx(jj,ii) + ER*Edisp(2*Esctrx-1);
                    dispmaty(jj,ii) = dispmaty(jj,ii) + ER*Edisp(2*Esctrx);
                    
                    % Calculate the stress values
                    stressvect{jj,ii} = Mater.C*Bu*Edisp(Edsctrx);
                    stressVM{jj,ii} = sqrt(stressvect{jj,ii}(1)^2 - stressvect{jj,ii}(1)*stressvect{jj,ii}(2) + stressvect{jj,ii}(2)^2 +3*stressvect{jj,ii}(3)^2);
                    
                    % Obtaining the values for Phase Field at the predefined points
                    RR = (PhPHTelem{indexPatch}(i).C)*squeeze(R(ii,jj,:));
                    RR = RR .* wgts;
                    w_sum = sum(RR);
                    RR = RR/w_sum;
                    
                    coord = RR'*cpts;
                    coordt{jj,ii} = coord;
                    % Calculate phase field values
                    wmatx(jj,ii) = wmatx(jj,ii) + RR'*Phidisp(sctrx);
                    
                end
            end
            physcoord((elementCounter-1)*4+1,:) = coordt{1,1};
            physcoord((elementCounter-1)*4+2,:) = coordt{1,2};
            physcoord((elementCounter-1)*4+3,:) = coordt{2,2};
            physcoord((elementCounter-1)*4+4,:) = coordt{2,1};
            
            dispcoord((elementCounter-1)*4+1,:) = [dispmatx(1,1) dispmaty(1,1)];
            dispcoord((elementCounter-1)*4+2,:) = [dispmatx(1,2) dispmaty(1,2)];
            dispcoord((elementCounter-1)*4+3,:) = [dispmatx(2,2) dispmaty(2,2)];
            dispcoord((elementCounter-1)*4+4,:) = [dispmatx(2,1) dispmaty(2,1)];
            
            sigmacoord((elementCounter-1)*4+1, 1) = stressVM{1,1};
            sigmacoord((elementCounter-1)*4+2, 1) = stressVM{1,2}';
            sigmacoord((elementCounter-1)*4+3, 1) = stressVM{2,2}';
            sigmacoord((elementCounter-1)*4+4, 1) = stressVM{2,1}';
            
            wcoord((elementCounter-1)*4+1) = wmatx(1,1);
            wcoord((elementCounter-1)*4+2) = wmatx(1,2);
            wcoord((elementCounter-1)*4+3) = wmatx(2,2);
            wcoord((elementCounter-1)*4+4) = wmatx(2,1);
            
        end
    end
end

% Magnifcation factor for the tdisp plot
factor = 10;

plot3 = subplot(2,2,3);
cla(plot3)
trisurf(vIEN,physcoord(:,1)+dispcoord(:,1)*factor, physcoord(:,2)+dispcoord(:,2)*factor, zeros(size(physcoord,1),1), log10(sigmacoord), 'EdgeColor','none','facecolor','interp')
view(0,90)
title('Displacements and \sigma_{VM}')
colorbar('vert')
colormap('parula')
hold on
axis tight
axis equal


plot4 = subplot(2,2,4);
cla(plot4)
factor = 0;
trisurf(vIEN,physcoord(:,1)+dispcoord(:,1)*factor, physcoord(:,2)+dispcoord(:,2)*factor, zeros(size(physcoord,1),1), wcoord, 'facecolor','interp', 'EdgeColor','none')
view(0,90)
title('Phase Field')
colorbar('vert')
colormap('parula')
hold on
axis tight
axis equal
drawnow

numVertices=4*size(vIEN,1);
end
