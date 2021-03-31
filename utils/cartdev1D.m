function Basis = cartdev1D(PHTelem,controlPts,geometry)

[Gpts,GWts]=GaussQuad(geometry.ngaussX);
[B_u,dB_u] = bernstein_basis(Gpts,geometry.p);

% Initialize the volume, dgdx and shape arrays
Basis = struct;
Basis.volume =cell(1,geometry.numPatches); 
Basis.dgdx = cell(1,geometry.numPatches,1);
Basis.shape = cell(1,geometry.numPatches);
Basis.gaussCord = cell(1,geometry.numPatches);

for indexPatch = 1:length(PHTelem)
    
    Basis.volume(indexPatch) = {cell(length(PHTelem{indexPatch}),1)};
    Basis.dgdx(indexPatch) = {cell(length(PHTelem{indexPatch}),1)};
    Basis.shape(indexPatch) = {cell(length(PHTelem{indexPatch}),1)};
    Basis.gaussCord(indexPatch) = {cell(length(PHTelem{indexPatch}),1)};
    for i=1:length(PHTelem{indexPatch})
        if isempty(PHTelem{indexPatch}(i).children)
            umin = PHTelem{indexPatch}(i).vertex(1);
            umax = PHTelem{indexPatch}(i).vertex(2);
            scalefac = (umax - umin)/2;
            nument = size(PHTelem{indexPatch}(i).C,1);
            nodes = PHTelem{indexPatch}(i).nodes(1:nument);
            cpts = controlPts{indexPatch}(nodes,1);
            wgts = controlPts{indexPatch}(nodes,2);
            kgauss = 0;
            for igaussX=1:geometry.ngaussX
                kgauss = kgauss + 1;
                RR = (PHTelem{indexPatch}(i).C)*B_u(igaussX,:)';
                dRdx = (PHTelem{indexPatch}(i).C)*dB_u(igaussX,:)'*2/(umax-umin);
                RR = RR .* wgts;
                dRdx = dRdx .* wgts;
                w_sum = sum(RR);
                dw_xi = sum(dRdx);
                dRdx = dRdx/w_sum - RR*dw_xi/w_sum^2;
                RR = RR/w_sum;
                dR  = dRdx';
                dxdxi = dR*cpts;
                coord = RR'*cpts;
                Basis.gaussCord{indexPatch}{i}(kgauss,:) = coord(:);
                % Solve for first derivatives in global coordinates
                dR = dxdxi\dR;
                J = det(dxdxi);
                
                Basis.volume{indexPatch}{i}(kgauss) = J*scalefac*GWts(igaussX);
                Basis.shape{indexPatch}{i}(kgauss,:) = RR(:);
                Basis.dgdx{indexPatch}{i}(kgauss,:) = dR(1,:);
            end
        end
    end
end
end