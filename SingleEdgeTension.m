% Script for a two dimensional plate under tension
% Uses a coarser elastic mesh and finer phase mesh

close all
clear

addpath('./utils')
addpath('./example_data')
addpath('./nurbs/inst')

savePlot = input('Do you want to save the mesh refinement and the phase field plots: \n Yes : 1 \n No  : 2 \n Your choice: ','s');

file_name = 'FD-tensile.txt';
output = fopen(file_name,'w');
fprintf(output,'%14.6e %14.6e\n',0,0);
fclose(output);
scrsz = get(groot, 'ScreenSize');
hFig = figure('Position',[1 scrsz(4)/6 3*scrsz(3)/5 3*scrsz(4)/4]);

[geometry,Mater,Fract,Integ] = Input_tensilePlate;
[EPHTelem,PhPHTelem,EcontrolPts,PhcontrolPts,EmeshInfo,PhmeshInfo] = modelDualMesh(geometry,'singleEdgeTension');

EBasis = cartdev(EPHTelem,EcontrolPts,geometry);
PhBasis = cartdev(PhPHTelem,PhcontrolPts,geometry);
Edirichlet = initialBC_tensile(EPHTelem,geometry);
fenerg = history_tensile(PhBasis,Fract,PhPHTelem,geometry); % Uses the data of the fracture mesh

% The Transfer function which trafers data from the Elastic mesh to the Phase Field Mesh.
PhTrans = Elastic2Phase(PhPHTelem,EPHTelem,EcontrolPts,geometry);
ETrans = Phase2Elastic(EPHTelem,PhPHTelem,geometry);

Phidisp = zeros(PhmeshInfo.sizeBasis,1); % Solution Vector for phase field on Phase Field Mesh
Edisp = zeros(geometry.dim*EmeshInfo.sizeBasis,1); % Solution Vector for the displacement on Elastic Mesh

for istep = 1:Integ.nstep
    
    istep
    
    if (istep < Integ.numStepsLimit)
        Integ.tfacto = Integ.tfacto + Integ.dfacto1;
    else
        Integ.tfacto = Integ.tfacto + Integ.dfacto2;
    end
    
    % Begin Inner Iteration
    normInnerStep = Inf;
    miter = 0;
    while (normInnerStep > geometry.toler)
        
        disp('Assembling the stiffness matrix.')
        tic
        [stiffUU,markRef,quadRef] = gStiffnessUU(EPHTelem,EmeshInfo,EBasis,ETrans,PhPHTelem,Phidisp,Mater,geometry);
        disp('Imposing boundary conditions and solving.')
        Edisp = applyBoundary2D(Edirichlet,Integ.tfacto,stiffUU,Edisp);
        toc
        
        tic
        disp('Update and internal forces.')
        fenerg = internalForces(PhPHTelem,PhTrans,EPHTelem,Edisp,Mater,fenerg,geometry);
        disp('Updating phase field.')
        [stiffPhiPhi,RHSPhi] = gStiffnessPhiPhi(PhPHTelem,PhmeshInfo,geometry,Fract,Mater,PhBasis,fenerg);
        solPhiOld = Phidisp;
        Phidisp = stiffPhiPhi\RHSPhi;
        miter = miter + 1
        normInnerStep = norm(Phidisp-solPhiOld)/sqrt(PhmeshInfo.sizeBasis)
        toc
        
        tic
        refFlag = 0;
        clear stiffPhiPhi
        if sum([markRef{:}])
            solPhiPatch = transferFieldGlob2Loc(PhPHTelem,PhmeshInfo,Phidisp);
            refFlag = 1;
            [PhPHTelem,PhcontrolPts,PhmeshInfo,Phidisp,PhmarkRef] = refinePhase(markRef,geometry,solPhiPatch,PhcontrolPts,PhPHTelem,PhmeshInfo);
            title(['Phase Field Mesh for Loadstep', num2str(istep) ,' and Iteration ', num2str(miter),'.']);
            PhBasis = cartdevRefinement(PhPHTelem,PhcontrolPts,geometry,PhBasis,PhmarkRef);
            fenerg = history_tensile(PhBasis,Fract,PhPHTelem,geometry);
            clear stiffUU
        end
        
        if refFlag
            if sum([quadRef{:}])
                [EPHTelem,EcontrolPts,EmeshInfo,EmarkRef] = refineElastic(quadRef,EmeshInfo,EPHTelem,EcontrolPts,geometry);
                title(['Elastic Mesh for Loadstep ',num2str(istep),' and Iteration ',num2str(miter),'.']);
                Edirichlet = initialBC_tensile(EPHTelem,geometry);
                EBasis = cartdevRefinement(EPHTelem,EcontrolPts,geometry,EBasis,EmarkRef);
            else
                EmarkRef = cell(1,geometry.numPatches);
                for indexPatch = 1:geometry.numPatches
                    EmarkRef(indexPatch) = {zeros(length(EPHTelem{indexPatch}),1)};
                end
            end
            Edisp = zeros(geometry.dim*EmeshInfo.sizeBasis,1); % Solution Vector
            [PhmarkRef,EmarkRef] = markRefinement(PhmarkRef,EmarkRef,PhTrans.ERefelem,ETrans.PhRefelem);
            PhTrans = Elastic2PhaseRefinement(PhPHTelem,EPHTelem,EcontrolPts,PhTrans,PhmarkRef,geometry);
            ETrans = Phase2ElasticRefinement(EPHTelem,PhPHTelem,PhcontrolPts,ETrans,EmarkRef,geometry);
            normInnerStep = Inf; % Iterate at least one more time after refinement
            clear PhmarkRef EmarkRef
        end
        toc
    end
    disp('Print data for Force-Displacement curves.')
    compTreacTension(stiffUU,Edisp,Integ.tfacto,Edirichlet,file_name);
    tic
    if(mod(istep,Integ.nprint) == 0)% Print Results
        fprintf('Done step: %5d\n',istep);
        plotDispPhaseTransfer2D(PhPHTelem,EPHTelem,PhcontrolPts,EcontrolPts,Phidisp,Edisp,PhmeshInfo,geometry,Mater);
        plot1 = subplot(2,2,1);
        title(['Elastic Mesh for Loadstep ',num2str(istep),'.']);
        plot2 = subplot(2,2,2);
        title(['Phase Field Mesh for Loadstep ', num2str(istep) ,'.']);
        if savePlot == '1'
            saveas(hFig, ['Loadstep', num2str(istep),'.png'])
        end
    end %if
    toc
end %istep
plotFD(file_name)
saveas(gcf,'ForceDisp_Tension.png')