function[shape,dgdx,Refelem] = transferCalculation1D(ElmtList,PHTelem,controlPts,p,xmin,xmax,nGauss)

dim = 1;
ngaussX = nGauss;
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
        for  ileaf = 1 : 2 % Coz each element has 2 children
            leaf(counter) = leafTemp(ileaf);
            vertex(counter,:) = PHTelem(leafTemp(ileaf)).vertex;
            counter = counter + 1;
        end
    end
end
kgauss=0;
shape = zeros(ngaussX,(p+1));
dgdx = zeros(ngaussX,dim,(p+1));
Refelem = zeros(1,ngaussX);
for ii=1:ngaussX
    kgauss=kgauss+1;
    u_hat = Gpts(ii);
    
    % Mapping (u_hat,v_hat) on [-1,1]x[-1,1] to (xi,eta) on[xmin,xmax]x[ymin,ymax]
    xi = u_hat*(xmax-xmin)/2+(xmax+xmin)/2;
    
    for iElmt = 1:size(vertex,1)
        if (xi <= vertex(iElmt,2)+toleq) && (xi >= vertex(iElmt,1)-toleq)
            elmtNum = leaf(iElmt);
            Pxmin = vertex(iElmt,1);
            Pxmax = vertex(iElmt,2);
            break;
        end
    end
    uu_hat = (2*xi - Pxmin - Pxmax)/(Pxmax-Pxmin);
    
    nument = size(PHTelem(elmtNum).C,1);
    nodes = PHTelem(elmtNum).nodes(1:nument);
    cpts = controlPts(nodes,1);
    wgts = controlPts(nodes,2);
    shape(kgauss,:) = phtBasis1D(uu_hat,PHTelem(elmtNum).C,p);
    [~,dR] = phtBasisIso1D(uu_hat,PHTelem(elmtNum).C,p,wgts);
    dR(1,:) = dR(1,:)*2/(Pxmax-Pxmin);
    dxdxi = dR*cpts;
    dR = dxdxi\dR;
    dgdx(kgauss,:,:) = dR;
    Refelem(kgauss) = elmtNum;
end %igaus
end
