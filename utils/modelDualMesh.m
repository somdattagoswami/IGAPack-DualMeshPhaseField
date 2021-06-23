function [EPHTelem,PhPHTelem,EcontrolPts,PhcontrolPts,EmeshInfo,PhmeshInfo] = modelDualMesh(geometry,example)
% Generates the initial coarse mesh for the elastic field and phase field
% Non-uniform refinement for cubic PHT splines

if example == "elasticBar"
    [PHTelem,controlPts,meshInfo] = initMesh_elasticBar(geometry);
    [PHTelem,meshInfo] = zipConforming1D(PHTelem,meshInfo);
    
elseif example == "singleEdgeTension" || example == "singleEdgeShear"
    [PHTelem,controlPts,meshInfo] = initMesh_tensile(geometry);
    [PHTelem,controlPts,meshInfo] = checkConforming(PHTelem,controlPts,meshInfo,geometry);
    [PHTelem,meshInfo] = zipConforming(PHTelem,geometry,meshInfo);
    
elseif example == "singleEdgeTension_distort"
    [PHTelem,controlPts,meshInfo] = initMesh_tensile_distort(geometry);
    [PHTelem,controlPts,meshInfo] = checkConforming(PHTelem,controlPts,meshInfo,geometry);
    [PHTelem,meshInfo] = zipConforming(PHTelem,geometry,meshInfo);
    
elseif example == "plateWHoles"
    [PHTelem,controlPts,meshInfo] = initMesh_plateWHole(geometry);
    [PHTelem,controlPts,meshInfo] = checkConforming(PHTelem,controlPts,meshInfo,geometry);
    [PHTelem,meshInfo] = zipConforming(PHTelem,geometry,meshInfo);
    
elseif example == "threePointBending"
    [PHTelem,controlPts,meshInfo] = initMesh_3ptPlate(geometry);
    [PHTelem,controlPts,meshInfo] = checkConforming(PHTelem,controlPts,meshInfo,geometry);
    [PHTelem,meshInfo] = zipConforming(PHTelem,geometry,meshInfo);
    
elseif example == "tensileCube"
    [PHTelem,controlPts,meshInfo] = initMesh_cube(geometry);
    [PHTelem,controlPts,meshInfo] = checkConforming3D(PHTelem,controlPts,meshInfo,geometry);
    [PHTelem,meshInfo] = zipConforming3D(PHTelem,meshInfo,geometry);
    
elseif example == "cubeWHole"
    [PHTelem,controlPts,meshInfo] = initMesh_cubeWHole(geometry);
    [PHTelem,controlPts,meshInfo] = checkConforming3D(PHTelem,controlPts,meshInfo,geometry);
    [PHTelem,meshInfo] = zipConforming3D(PHTelem,meshInfo,geometry);
end

% The Initial evaluation mesh is the parent Mesh and which is same for
% Elastic field and the Phase field
EPHTelem = PHTelem;
EcontrolPts = controlPts;
EmeshInfo = meshInfo;

PhPHTelem = PHTelem;
PhcontrolPts = controlPts;
PhmeshInfo = meshInfo;

if geometry.dim == 2
    
    plot1 = subplot(2,2,1);
    cla(plot1)
    plotMesh(EPHTelem,EcontrolPts,geometry,0)
    axis equal
    title('Intial Elastic Mesh');
    
    plot2 = subplot(2,2,2);
    cla(plot2)
    plotMesh(PhPHTelem,PhcontrolPts,geometry,0)
    axis equal
    title('Intial Phase Field Mesh');
end

end
