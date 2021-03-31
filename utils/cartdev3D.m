function Basis = cartdev3D(PHTelem,controlPts,geometry)
% Compute the derivatives using the PHTelem structure at Gauss points

p = geometry.p;
q = geometry.q;
r = geometry.r;

ngaussX = geometry.p+1;
ngaussY = geometry.q+1;
ngaussZ = geometry.r+1;
mgauss = ngaussX*ngaussY*ngaussZ;

[GptsX,GWtsX]=GaussQuad(ngaussX);
[GptsY,GWtsY]=GaussQuad(ngaussY);
[GptsZ,GWtsZ]=GaussQuad(ngaussZ);

%1D bernstein polynomials evaluated at the Gauss points on the master element
[B_u,dB_u,~] = bernstein_basis2nd(GptsX,geometry.p);
[B_v,dB_v,~] = bernstein_basis2nd(GptsY,geometry.q);
[B_w,dB_w,~] = bernstein_basis2nd(GptsZ,geometry.r);

dBdu = zeros(ngaussX,ngaussY,ngaussZ,(p+1)*(q+1)*(r+1));
dBdv = zeros(ngaussX,ngaussY,ngaussZ,(p+1)*(q+1)*(r+1));
dBdw = zeros(ngaussX,ngaussY,ngaussZ,(p+1)*(q+1)*(r+1));
R = zeros(ngaussX,ngaussY,ngaussZ,(p+1)*(q+1)*(r+1));

% The derivatives of the 2D Bernstein polynomials at Gauss points on the master element
basisCounter = 0;
for k=1:r+1
    for j=1:q+1
        for i=1:p+1
            basisCounter = basisCounter + 1;
            for kk=1:ngaussZ
                for jj=1:ngaussY
                    for ii=1:ngaussX
                        R(ii,jj,kk,basisCounter) = B_u(ii,i)*B_v(jj,j)*B_w(kk,k);
                        dBdu(ii,jj,kk,basisCounter) = dB_u(ii,i)*B_v(jj,j)*B_w(kk,k);                        
                        dBdv(ii,jj,kk,basisCounter) = B_u(ii,i)*dB_v(jj,j)*B_w(kk,k);                        
                        dBdw(ii,jj,kk,basisCounter) = B_u(ii,i)*B_v(jj,j)*dB_w(kk,k); 
                    end
                end
            end
        end
    end
end

% Initialize the volume, dgdx and shape arrays
Basis = struct;
Basis.volume =cell(1,geometry.numPatches); 
Basis.dgdx = cell(1,geometry.numPatches,1);
Basis.shape = cell(1,geometry.numPatches);
Basis.gaussCord = cell(1,geometry.numPatches);
Basis.sctrxElem = cell(1,geometry.numPatches);
Basis.elemLevel = cell(1,geometry.numPatches);
Basis.RefElem = cell(1,geometry.numPatches);

for indexPatch = 1:length(PHTelem)
    
    Basis.volume(indexPatch) = {cell(length(PHTelem{indexPatch}),1)};
    Basis.dgdx(indexPatch) = {cell(length(PHTelem{indexPatch}),1)};
    Basis.shape(indexPatch) = {cell(length(PHTelem{indexPatch}),1)};
    Basis.gaussCord(indexPatch) = {cell(length(PHTelem{indexPatch}),1)};
    Basis.sctrxElem(indexPatch) = {cell(length(PHTelem{indexPatch}),1)};
    Basis.RefElem(indexPatch) = {cell(length(PHTelem{indexPatch}),1)};
    Basis.elemLevel{indexPatch} = cell(length(PHTelem{indexPatch}),1);
    
    for i=1:length(PHTelem{indexPatch})
        if isempty(PHTelem{indexPatch}(i).children)
            Basis.elemLevel{indexPatch}{i} = PHTelem{indexPatch}(i).level;
            Basis.RefElem{indexPatch}{i}(1:ngaussX*ngaussY*ngaussZ) = i;
            xmin = PHTelem{indexPatch}(i).vertex(1);
            xmax = PHTelem{indexPatch}(i).vertex(4);
            ymin = PHTelem{indexPatch}(i).vertex(2);
            ymax = PHTelem{indexPatch}(i).vertex(5);
            zmin = PHTelem{indexPatch}(i).vertex(3);
            zmax = PHTelem{indexPatch}(i).vertex(6);
            
            % The jacobian of the transformation from [-1,1]x[-1,1] to [xmin,xmax]x[ymin,ymax]
            scalefac = (xmax - xmin)*(ymax - ymin)*(zmax-zmin)/8;
            nument = size(PHTelem{indexPatch}(i).C,1);
            nodes = PHTelem{indexPatch}(i).nodes(1:nument);
            cpts = controlPts{indexPatch}(nodes,1:3);
            wgts = controlPts{indexPatch}(nodes,4);
            Basis.sctrxElem{indexPatch}{i} = PHTelem{indexPatch}(i).nodesGlobal(1:nument);
            kgauss = 0;
            for kk = 1:ngaussZ
                for jj=1:ngaussY
                    for ii=1:ngaussX
                        kgauss = kgauss + 1;
                        dRdx = (PHTelem{indexPatch}(i).C)*squeeze(dBdu(ii,jj,kk,:))*2/(xmax-xmin);
                        dRdy = (PHTelem{indexPatch}(i).C)*squeeze(dBdv(ii,jj,kk,:))*2/(ymax-ymin);
                        dRdz = (PHTelem{indexPatch}(i).C)*squeeze(dBdw(ii,jj,kk,:))*2/(zmax-zmin);
                        RR = (PHTelem{indexPatch}(i).C)*squeeze(R(ii,jj,kk,:));
                        
                        RR = RR .* wgts;
                        dRdx = dRdx .* wgts;
                        dRdy = dRdy .* wgts;
                        dRdz = dRdz .* wgts;
                        
                        w_sum = sum(RR);
                        dw_xi = sum(dRdx);
                        dw_eta = sum(dRdy);
                        dw_zeta = sum(dRdz);

                        dRdx = dRdx/w_sum - RR*dw_xi/w_sum^2;
                        dRdy = dRdy/w_sum - RR*dw_eta/w_sum^2;
                        dRdz = dRdz/w_sum - RR*dw_zeta/w_sum^2;
                        RR = RR/w_sum;
                        
                        % multiply by the jacobian of the transformation from reference
                        % space to the parameter space
                        
                        dR  = [dRdx';dRdy';dRdz'];                         
                        dxdxi = dR*cpts;
                        coord = RR'*cpts;
                        Basis.gaussCord{indexPatch}{i}(:,kgauss) = coord;
                        
                        % Solve for first derivatives in global coordinates
                        dR = dxdxi\dR;
                        J = abs(det(dxdxi));
                        Basis.volume{indexPatch}{i}(kgauss) = J*scalefac*GWtsX(ii).*GWtsY(jj).*GWtsZ(kk);
                        Basis.shape{indexPatch}{i}(:,kgauss) = RR;
                        Basis.dgdx{indexPatch}{i}(:,:,kgauss) = dR;            
                    end
                end
            end
        end
    end
end
end

