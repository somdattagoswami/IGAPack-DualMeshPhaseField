% Script for a elastic bar with a crack at the center
% The bar is subjected to sinosoidal loading
% Uses a coarser elastic mesh and finer phase mesh
% Implement history function and hybrid staggered solver
close all
clear

addpath ./utils
addpath ./example_data
addpath ./nurbs/inst
scrsz = get(groot, 'ScreenSize');
hFig = figure('Position',[1 scrsz(4)/6 3*scrsz(3)/5 3*scrsz(4)/4]);

[geometry,Mater,Fract] = Input_elasticBar;
[EPHTelem,PhPHTelem,EcontrolPts,PhcontrolPts,EmeshInfo,PhmeshInfo] = modelDualMesh(geometry,'elasticBar');
Edirichlet = initialBC_elasticBar(EPHTelem{1});
EBasis = cartdev1D(EPHTelem,EcontrolPts,geometry);
PhBasis = EBasis;
Fract.constl = 2*Fract.constl;
[fenerg]=history_elasticBar(PhBasis,Fract,PhPHTelem,geometry); % Uses the data of the fracture mesh
Fract.constl = Fract.constl/2;

Phidisp = zeros(PhmeshInfo.sizeBasis,1); % Solution Vector for phase field on Phase Field Mesh
Edisp = zeros(geometry.dim*EmeshInfo.sizeBasis,1); % Solution Vector for the displacement on Elastic Mesh

exactSolU = @(x)exact_sol_crack(x);
exactPhi = @(x,Fract)exact_sol_phi(x,Fract);

[PhdgdxTrans,ERefElem] = Elastic2Phase1D(PhPHTelem,EPHTelem,EcontrolPts,geometry);% The Transfer function which trafers data from the Elastic mesh to the Phase Field Mesh. It also contains a reference for
[EshapeTrans,PhRefElem] = Phase2Elastic1D(EPHTelem,PhPHTelem,geometry);

% Begin inner iteration
normInnerStep = Inf;
miter = 0;
convergence=[];

while (normInnerStep > geometry.toler)
    
    disp('Assembling the stiffness matrix...')
    tic
    [stiffUU,Edisp,elemRef,tupleRef] = gStiffnessUU1D(EPHTelem,EmeshInfo,EBasis,EshapeTrans,PhPHTelem,Phidisp,PhRefElem,Fract,Mater,geometry,Edirichlet);
    toc
    
    disp(['Number of degrees of freedom: ', num2str(length(Edisp))])
    [l2relerr] = calcErrorNorms_1D(Edisp,EPHTelem,EcontrolPts,geometry,exactSolU);
    
    if (miter==0) || (length(Edisp) > convergence(miter,1))
        miter = miter + 1
    end
    
    convergence(miter,1:2) = [length(Edisp),l2relerr];
    tic
    disp('Updating phase field...')
    
    [stiffPhiPhi,RHSPhi] = gStiffnessPhiPhi1D(PhPHTelem,PhmeshInfo,PhBasis,fenerg,geometry,Fract);
    solPhiOld = Phidisp;
    Phidisp = stiffPhiPhi\RHSPhi;
    toc
    
    plot1 = subplot(1,2,1);
    cla(plot1)
    plotErrorSol1DIso(PhPHTelem,PhcontrolPts,Phidisp,geometry.p,exactPhi,Fract)
    set(get(gca,'ylabel'),'String','$\phi(x)$','FontSize',16','FontWeight','bold','FontName','Times','Interpreter','latex')
    set(get(gca,'xlabel'),'String','Axial co-ordinate of the bar','FontSize',16','FontWeight','bold','FontName','Times','Interpreter','tex')

    plot2 = subplot(1,2,2);
    cla(plot2)
    plotErrorSol1DIso(EPHTelem,EcontrolPts,Edisp,geometry.p,exactSolU)
    set(get(gca,'ylabel'),'String','Displacement','FontSize',16','FontWeight','bold','FontName','Times','Interpreter','tex')
    set(get(gca,'xlabel'),'String','Axial co-ordinate of the bar','FontSize',16','FontWeight','bold','FontName','Times','Interpreter','tex')
    saveas(hFig, ['Iteration', num2str(miter),'.png'])
    
    normInnerStep = norm(Phidisp-solPhiOld)/sqrt(PhmeshInfo.sizeBasis)
    refFlag = 0;    
    
    if sum([elemRef{:}])
        refFlag = 1;
        [PhmarkRef,EmarkRef,PhPHTelem,PhcontrolPts,Phidisp,EPHTelem,EcontrolPts,EmeshInfo,PhmeshInfo] = ...
            refineMesh1D(elemRef,tupleRef,PhPHTelem,EPHTelem,EcontrolPts,PhcontrolPts,PhmeshInfo,EmeshInfo,geometry,Phidisp);
        PhBasis = cartdev1D(PhPHTelem,PhcontrolPts,geometry);
        [fenerg]=history_elasticBar(PhBasis,Fract,PhPHTelem,geometry); % Uses the data of the fracture mesh
        Edirichlet = initialBC_elasticBar(EPHTelem{1});
        EBasis = cartdev1D(EPHTelem,EcontrolPts,geometry);
    end
        
    if refFlag        
        Edisp = zeros(geometry.dim*EmeshInfo.sizeBasis,1); % Solution Vector
        [PhmarkRef,EmarkRef] = markRefinement(PhmarkRef,EmarkRef,ERefElem,PhRefElem);% marks the elemnts on either meshes that have to checked which have been refined on either of the meshes
        [PhdgdxTrans,ERefElem] = Elastic2PhaseRefinement1D(PhPHTelem,EPHTelem,EcontrolPts,PhdgdxTrans,ERefElem,PhmarkRef,geometry);
        [EshapeTrans,PhRefElem] = Phase2ElasticRefinement1D(EPHTelem,PhPHTelem,PhcontrolPts,EshapeTrans,PhRefElem,EmarkRef,geometry);
        normInnerStep = Inf; % Iterate at least one more time after refinement
    end
    toc
end
figure
loglog(convergence(:,1),convergence(:,2),'-o','LineWidth',3)
set(get(gca,'ylabel'),'String','Relative $\mathcal{L}_2$ norm','FontSize',16','FontWeight','bold','FontName','Times','Interpreter','latex')
set(get(gca,'xlabel'),'String','Degrees of freedom','FontSize',16','FontWeight','bold','FontName','Times','Interpreter','tex')
saveas(gcf,'Convergence.png')