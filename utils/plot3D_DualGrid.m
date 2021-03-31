function plot3D_DualGrid(PhPHTelem,EPHTelem,Phidisp,Edisp,PhcontrolPts,EcontrolPts,geometry,vtuFile,PhmeshInfo,fudge)
% Evaluate the displacements, Strain and Potential values at the element vertex
maxRefLevel = geometry.maxERefLevel;
toleq = 1e-10;
if nargin < 10
   fudge = 0.0;
end
nx = 2;
ny = 2;
nz = 2;
p = geometry.p;
q = geometry.q;
r = geometry.r;
px = linspace(-1+fudge,1-fudge,nx);
py = linspace(-1+fudge,1-fudge,ny);
pz = linspace(-1+fudge,1-fudge,nz);
noGpEle = nx*ny*nz;% Number of Gauss points per element

noElems = sum(PhmeshInfo.numElements);
physcoord = zeros(3,noGpEle,noElems);  % Global coords of Gauss points
dispcoord     = zeros(3,noGpEle,noElems);  % Displacements of Gauss points
epsiloncoord = zeros(6,noGpEle,noElems);% Strain of Gauss points
wcoord = zeros(1,noGpEle,noElems);  % Phase Field at Gauss points

% 1D Bernstein polynomials evaluated at the Gauss points on the master element
[B_u,~] = bernstein_basis(px,p);
[B_v,~] = bernstein_basis(py,q);
[B_w,~] = bernstein_basis(pz,r);

R = zeros(nx,ny,nz,(p+1)*(q+1)*(r+1));

basisCounter = 0;
for k=1:r+1
    for j=1:q+1
        for i=1:p+1
            basisCounter = basisCounter + 1;
            for kk=1:nz
                for jj=1:ny
                    for ii=1:nx
                        R(ii,jj,kk,basisCounter) = B_u(ii,i)*B_v(jj,j)*B_w(kk,k);
                    end
                end
            end
        end
    end
end

elementCounter = 1;
for indexPatch = 1:length(PhPHTelem)
    for i=1:length(PhPHTelem{indexPatch})
        if isempty(PhPHTelem{indexPatch}(i).children)
            
            ximin = PhPHTelem{indexPatch}(i).vertex(1);
            ximax = PhPHTelem{indexPatch}(i).vertex(4);
            etamin = PhPHTelem{indexPatch}(i).vertex(2);
            etamax = PhPHTelem{indexPatch}(i).vertex(5);
            zetamin = PhPHTelem{indexPatch}(i).vertex(3);
            zetamax = PhPHTelem{indexPatch}(i).vertex(6);
            
            nument = size(PhPHTelem{indexPatch}(i).C,1);
            nodes = PhPHTelem{indexPatch}(i).nodes(1:nument);
            sctrx = PhPHTelem{indexPatch}(i).nodesGlobal(1:nument);
            
            cpts = PhcontrolPts{indexPatch}(nodes,1:3);
            wgts = PhcontrolPts{indexPatch}(nodes,4);
            
            indexParent = i;
            while PhPHTelem{indexPatch}(indexParent).level > 1
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
            
            gp = 1;
            for kk=1:nz
                for jj=1:ny
                    for ii=1:nx
                        u_hat = px(ii);
                        v_hat = py(jj);
                        w_hat = pz(kk);
                        % Mapping (u_hat,v_hat) on [-1,1]x[-1,1] to (xi,eta) on[xmin,xmax]x[ymin,ymax]
                        xi = u_hat*(ximax-ximin)/2+(ximax+ximin)/2;
                        eta = v_hat*(etamax-etamin)/2+(etamax+etamin)/2;
                        zeta = w_hat*(zetamax-zetamin)/2+(zetamax+zetamin)/2;
                        
                        for iElmt = 1:size(vertex,1)
                            if (xi <= vertex(iElmt,4)+toleq) && (xi >= vertex(iElmt,1)-toleq) && (eta <= vertex(iElmt,5)+toleq) && (eta >= vertex(iElmt,2)-toleq) && (zeta <= vertex(iElmt,6)+toleq) && (zeta>=vertex(iElmt,3)-toleq)
                                EelmtNum = leaf(iElmt);
                                Pxmin = vertex(iElmt,1);
                                Pxmax = vertex(iElmt,4);
                                Pymin = vertex(iElmt,2);
                                Pymax = vertex(iElmt,5);
                                Pzmin = vertex(iElmt,3);
                                Pzmax = vertex(iElmt,6);
                                break;
                            end
                        end
                        uu_hat = (2*xi - Pxmin - Pxmax)/(Pxmax-Pxmin);
                        vv_hat = (2*eta - Pymin - Pymax)/(Pymax-Pymin);
                        ww_hat = (2*zeta - Pzmin - Pzmax)/(Pzmax-Pzmin);
                        Enument = size(EPHTelem{indexPatch}(EelmtNum).C,1);
                        nodes = EPHTelem{indexPatch}(EelmtNum).nodes(1:Enument);
                        Ecpts = EcontrolPts{indexPatch}(nodes,1:3);
                        Ewgts = EcontrolPts{indexPatch}(nodes,4);
                        Esctrx = EPHTelem{indexPatch}(EelmtNum).nodesGlobal(1:Enument);                        
                    
                         ER = phtBasis3D(uu_hat,vv_hat,ww_hat,EPHTelem{indexPatch}(EelmtNum).C,p,q,r);
                        [~,EdR] = phtBasisIso3D(uu_hat,vv_hat,ww_hat,EPHTelem{indexPatch}(EelmtNum).C,p,q,r,Ewgts);
                        EdR(1,:) = EdR(1,:)*2/(Pxmax-Pxmin);
                        EdR(2,:) = EdR(2,:)*2/(Pymax-Pymin);
                        EdR(3,:) = EdR(3,:)*2/(Pzmax-Pzmin);
                        Edxdxi = EdR*Ecpts;
                        EdR = Edxdxi\EdR;  
                        
                        cdRdx = EdR(1,:);
                        cdRdy = EdR(2,:);
                        cdRdz = EdR(3,:);
                        
                        % Calculate displacement values
                        dispcoord(1,gp,elementCounter) = ER*Edisp(3*Esctrx-2);
                        dispcoord(2,gp,elementCounter) = ER*Edisp(3*Esctrx-1);
                        dispcoord(3,gp,elementCounter) = ER*Edisp(3*Esctrx);
                        
                        % u_disp = u(:,gp,elementCounter)
                        
                        strain11 = cdRdx*Edisp(3*Esctrx-2);
                        strain22 = cdRdy*Edisp(3*Esctrx-1);
                        strain33 = cdRdz*Edisp(3*Esctrx);
                        strain12 = cdRdy*Edisp(3*Esctrx-2) + cdRdx*Edisp(3*Esctrx-2);
                        strain23 = cdRdy*Edisp(3*Esctrx) + cdRdz*Edisp(3*Esctrx-1);
                        strain31 = cdRdz*Edisp(3*Esctrx-2) + cdRdx*Edisp(3*Esctrx);
                        
                        RR = (PhPHTelem{indexPatch}(i).C)*squeeze(R(ii,jj,kk,:));                        
                        RR = RR .* wgts;
                        w_sum = sum(RR);
                        RR = RR/w_sum;
                        coord = RR'*cpts;
                        
                        physcoord(1,gp,elementCounter) = coord(1);
                        physcoord(2,gp,elementCounter) = coord(2);
                        physcoord(3,gp,elementCounter) = coord(3);
                        epsiloncoord(:,gp,elementCounter) = [strain11;strain22;strain33;strain12;strain23;strain31];
                        wcoord(1,gp,elementCounter)= RR'*Phidisp(sctrx);
                        gp = gp +1 ;
                    end
                end
            end
            elementCounter = elementCounter + 1;
        end
    end
end
msh_to_vtu_3dVM(physcoord,dispcoord,epsiloncoord,wcoord,[nx ny nz],vtuFile);
end