% Script for a cube under tensile loading
% Uses a coarser elastic mesh and finer phase mesh
% Run the make_c0_cube_with_hole.m file from the example_data folder to
% create the initial geometry.

close all
clear

addpath ./utils
addpath ./example_data
addpath ./nurbs/inst

file_name = 'FD-CubeWithHole.txt';
output = fopen(file_name,'w');
fprintf(output,'%14.6e %14.6e\n',0,0);
fclose(output);

% Importing the Input File
[geometry,Mater,Fract,Integ] = Input_cubeWHole;
addfun = @plus;
[EPHTelem,PhPHTelem,EcontrolPts,PhcontrolPts,EmeshInfo,PhmeshInfo] = modelDualMesh(geometry,'cubeWHole');
[Edirichlet] = initialBC_cubeWHole(EPHTelem,geometry);
EBasis = cartdev3D(EPHTelem,EcontrolPts,geometry);
PhBasis = EBasis;

Phidisp = zeros(PhmeshInfo.sizeBasis,1); % Solution Vector for phase field on Phase Field Mesh
Edisp = zeros(geometry.dim*EmeshInfo.sizeBasis,1); % Solution Vector for the displacement on Elastic Mesh

Fract.constl = 16*Fract.constl;
[fenerg] = history_cubeWHole(PhBasis,Fract,PhPHTelem,geometry);
Fract.constl = Fract.constl/16;

PhdgdxTrans = EBasis.dgdx;
EshapeTrans = PhBasis.shape;
PhRefElem = EBasis.RefElem;
ERefElem = EBasis.RefElem;

for istep = 1:Integ.nstep+1
    
    istep
    
    if (istep < Integ.numStepsLimit)
        Integ.tfacto = Integ.tfacto + Integ.dfacto1;
    else
        Integ.tfacto = Integ.tfacto + Integ.dfacto2;
    end
    
    % Begin inner iteration
    normInnerStep = Inf;
    miter = 0;
    while (normInnerStep > geometry.toler)
        tic
        disp('Assembling the stiffness matrix.')
        [stiffUU,PhmarkRef,EoctupleRef,Edisp] = computeUU(geometry,EmeshInfo,Fract,Mater,...
            Edirichlet,Integ.tfacto,Edisp,Phidisp,EshapeTrans,PhRefElem,EBasis,PhBasis);
        toc
        
        solPhiOld = Phidisp;
        [Phidisp,fenerg] = computePhiPhi(geometry,Fract,Mater,Edisp,EBasis,PhBasis,fenerg,PhdgdxTrans,ERefElem,PhmeshInfo);
        miter = miter + 1
        normInnerStep = norm(Phidisp - solPhiOld)
        clear solPhiOld
        toc
        refFlag = 0;
        
        if sum([PhmarkRef{:}])
            refFlag = 1;
            clear stiffUU
            [PhPHTelem,EPHTelem,PhcontrolPts,EcontrolPts,PhmeshInfo,EmeshInfo,Phidisp,PhmarkRef,EmarkRef] = ...
                refineMesh(PhPHTelem,PhmeshInfo,Phidisp,geometry,PhcontrolPts,EPHTelem,EmeshInfo,EcontrolPts,EoctupleRef,PhmarkRef);
            PhBasis = cartdev3D(PhPHTelem,PhcontrolPts,geometry);
            fenerg0 = fracEnergy0(PhPHTelem,geometry);
            if (istep == 1) && (miter < 3)
                Fract.constl = 16*Fract.constl;
                [fenerg] = history_cubeWHole(PhBasis,Fract,PhPHTelem,geometry);
                Fract.constl = Fract.constl/16;
            else
                fenerg = fenerg0;
            end
            [Edirichlet] = initialBC_cubeWHole(EPHTelem,geometry);
            EBasis = cartdev3D(EPHTelem,EcontrolPts,geometry);
        end
        if refFlag
            
            Edisp = zeros(geometry.dim*EmeshInfo.sizeBasis,1); % Solution Vector
            [PhmarkRef,EmarkRef] = markRefinement(PhmarkRef,EmarkRef,ERefElem,PhRefElem);% Marks the elemnts on either meshes that have to checked which have been refined on either of the meshes
            [PhdgdxTrans,ERefElem] = Elastic2PhaseRefinement3D(PhPHTelem,EPHTelem,EcontrolPts,PhdgdxTrans,ERefElem,PhmarkRef,geometry);
            [EshapeTrans,PhRefElem] = Phase2ElasticRefinement3D(EPHTelem,PhPHTelem,PhcontrolPts,EshapeTrans,PhRefElem,EmarkRef,geometry);
            clear PhmarkRef EmarkRef
            
            normInnerStep = Inf; % Iterate at least one more time after refinement
        end
    end
    
    disp('Print data for force-tdisp curves')
    tic
    compTreacCube(stiffUU,Edisp,Integ.tfacto,Edirichlet,file_name)
    if(mod(istep,Integ.nprnt) == 0)% Print results
        numVCtrlElmt  = 8;
        vtuFile = ['time_',num2str(istep),'.vtu'];
        fudge = 1e-3;
        plot3D_DualGrid(PhPHTelem,EPHTelem,Phidisp,Edisp,PhcontrolPts,EcontrolPts,geometry,vtuFile,PhmeshInfo,fudge)
        fprintf('Done step: %5d\n',istep);
    end
    toc
end
plotFD(file_name)
saveas(gcf,'ForceDisp_cubeWHole.png')