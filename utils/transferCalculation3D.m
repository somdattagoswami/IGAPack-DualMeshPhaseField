function[shape,dgdx,Refelem] = transferCalculation3D(ElmtList,PHTelem,controlPts,geometry,xmin,xmax,ymin,ymax,zmin,zmax)

ngaussX = geometry.ngaussX;
ngaussY = geometry.ngaussY;
ngaussZ = geometry.ngaussZ;
[Gpts,~]=GaussQuad(ngaussX);
toleq = 1e-10;
counter = 1;
for iCount = 1:length(ElmtList)
    if isempty(PHTelem(ElmtList(iCount)).children)% This loop is for the elements that have been refined in the other mesh and has the reference elemnt number affected in the current mesh
        leaf(counter) = ElmtList(iCount);% Leaf contains the elemnt number ogf the most refined element. The elemt number stored in leaf does not have children
        vertex(counter,:) = PHTelem(ElmtList(iCount)).vertex;
        counter = counter + 1;
    else
        leafTemp = PHTelem(ElmtList(iCount)).children;% Children have to be checked since in the current miter, there might be refinement in other Mesh
        for  ileaf = 1 : 8 % Coz each element has 8 children
            leaf(counter) = leafTemp(ileaf);
            vertex(counter,:) = PHTelem(leafTemp(ileaf)).vertex;
            counter = counter + 1;
        end
    end
end
kgauss=0;
p = geometry.p;
q = geometry.q;
r = geometry.r;
nument = (p+1)*(q+1)*(r+1);
shape = zeros(nument,ngaussX*ngaussY*ngaussZ);
dgdx = zeros(3,nument,ngaussX*ngaussY*ngaussZ);
Refelem = zeros(1,ngaussX*ngaussY);

for ii = 1:ngaussX
    for jj = 1:ngaussY
        for kk = 1:ngaussZ
        kgauss=kgauss+1;
        u_hat = Gpts(ii);
        v_hat = Gpts(jj);
        w_hat = Gpts(kk);
        % Mapping (u_hat,v_hat) on [-1,1]x[-1,1] to (xi,eta) on[xmin,xmax]x[ymin,ymax]
        xi = u_hat*(xmax-xmin)/2+(xmax+xmin)/2;
        eta = v_hat*(ymax-ymin)/2+(ymax+ymin)/2;
        zeta = w_hat*(zmax-zmin)/2+(zmax+zmin)/2;
        
        for iElmt = 1:counter-1
            if (xi <= vertex(iElmt,4)+toleq) && (xi >= vertex(iElmt,1)-toleq) && (eta <= vertex(iElmt,5)+toleq) && (eta >= vertex(iElmt,2)-toleq) && (zeta <= vertex(iElmt,6)+toleq) && (zeta>=vertex(iElmt,3)-toleq)
                elmtNum = leaf(iElmt);
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
        nument = size(PHTelem(elmtNum).C,1);
        nodes = PHTelem(elmtNum).nodes(1:nument);
        cpts = controlPts(nodes,1:3);
        wgts = controlPts(nodes,4);
        shape(:,kgauss) = phtBasis3D(uu_hat,vv_hat,ww_hat,PHTelem(elmtNum).C,p,q,r);
        [~,dR] = phtBasisIso3D(uu_hat,vv_hat,ww_hat,PHTelem(elmtNum).C,p,q,r,wgts);
        dR(1,:) = dR(1,:)*2/(Pxmax-Pxmin);
        dR(2,:) = dR(2,:)*2/(Pymax-Pymin);
        dR(3,:) = dR(3,:)*2/(Pzmax-Pzmin);
        dxdxi = dR*cpts;
        dR = dxdxi\dR;
        dgdx(:,:,kgauss) = dR;
        Refelem(kgauss) = elmtNum;
        end
    end %jgaus
end %igaus