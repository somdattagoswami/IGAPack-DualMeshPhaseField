function [sol0] = applyBoundary2D(dirichlet,tfacto,stiffness,tdisp,RHS)
% Apply boundary conditions while maintaining symmetry and solve the linear system
% Apply full displacement
dim = 2;
fixedPts = length(dirichlet.XY);
bcdof = [];
bcval = [];
if nargin<5
    RHS = zeros(size(tdisp));
end

for ivfix = 1:fixedPts
    lnode = dirichlet.XY(ivfix);    
    for idofn =1:dim
        if(dirichlet.ValXY(ivfix,idofn) == 1)
            itotv =(lnode-1)*dim +idofn;
            dispFull = dirichlet.restrainedPts(ivfix,idofn).*tfacto;
            bcdof = [bcdof,itotv];
            bcval = [bcval,dispFull];
        end
    end
end

dof_all = 1:length(tdisp);
dof_int = setdiff(dof_all, bcdof);
stiffness2 = stiffness(dof_int, dof_int);
RHS=RHS-stiffness(:,bcdof)*bcval';
RHS2 = RHS(dof_int);

sol0_int = stiffness2\RHS2;
sol0 = zeros(length(tdisp),1);
sol0(bcdof) = bcval;
sol0(dof_int) = sol0_int;

end %endfunction