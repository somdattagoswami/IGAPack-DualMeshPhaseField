function [dirichlet] = initialBC_cubeWHole(PHTelem,geometry)
% Initialize the boundary conditions
p = geometry.p;
q = geometry.q;
r = geometry.r;

% Detect which nodes have support on the left and right boundary
bcdof_left = [];
bcdof_front = [];
bcdof_right = [];
bcdof_down = [];
bcdof_back = [];
bcdof_up = [];

% Define side node indices
down_nodes = 1:(p+1)*(q+1);
right_nodes = (p+1):(p+1):(p+1)*(q+1)*(r+1);
up_nodes = 1+(p+1)*(q+1)*r:(p+1)*(q+1)*(r+1);
left_nodes = 1:(p+1):(1+(p+1)*(q+1)*(r+1)-p);
front_nodes = [];

numPatches = length(PHTelem);

for i=1:r+1
    front_nodes = [front_nodes, (p+1)*(q+1)*(i-1)+1:(p+1)*((q+1)*(i-1)+1)];
end

back_nodes = front_nodes + (p+1)*q;

for indexPatch=1:numPatches
    for i=1:length(PHTelem{indexPatch})
        if isempty(PHTelem{indexPatch}(i).children)
            if isempty(PHTelem{indexPatch}(i).neighbor_front) && indexPatch <= 2
                %x fixed in patches 1 and 2, front side                
                bcdof_front = [bcdof_front, PHTelem{indexPatch}(i).nodesGlobal(front_nodes)];
                
            end
            if isempty(PHTelem{indexPatch}(i).neighbor_down)
                bcdof_down = [bcdof_down, PHTelem{indexPatch}(i).nodesGlobal(down_nodes)];
            end
            if isempty(PHTelem{indexPatch}(i).neighbor_up) && (indexPatch == 2)
                bcdof_up = [bcdof_up, PHTelem{indexPatch}(i).nodesGlobal(up_nodes)];
            end
            
            if isempty(PHTelem{indexPatch}(i).neighbor_left) && ((indexPatch == 1) || (indexPatch == 3))
                %z fixed in patches 1 and 3, left side
                bcdof_left = [bcdof_left, PHTelem{indexPatch}(i).nodesGlobal(left_nodes)];
            end
            if isempty(PHTelem{indexPatch}(i).neighbor_right) && ((indexPatch == 2) || (indexPatch == 4))
                %y fixed in patches 2 and 4, right side
                bcdof_right = [bcdof_right, PHTelem{indexPatch}(i).nodesGlobal(right_nodes)];
            end
            if isempty(PHTelem{indexPatch}(i).neighbor_back) && (indexPatch >= 3)
                %axis of revolution, fixed in the y and z directions in
                %patches 3 and 4
                bcdof_back = [bcdof_back, PHTelem{indexPatch}(i).nodesGlobal(back_nodes)];
            end
        end
    end
end

%remove duplicated entries
dirichlet.Front = unique(bcdof_front,'stable');
dirichlet.Right = unique(bcdof_right,'stable');
dirichlet.Left = unique(bcdof_left,'stable');
dirichlet.Back = unique(bcdof_back,'stable');
dirichlet.Top = unique(bcdof_up,'stable');
dirichlet.Down = unique(bcdof_down,'stable');


dirichlet.XYZ =[dirichlet.Front,dirichlet.Right,dirichlet.Back,dirichlet.Left,dirichlet.Top,dirichlet.Down]; % Contains the node numbers whose DOF/s have to be restrained.
dirichlet.ValXYZ = zeros(length(dirichlet.XYZ),3); %  0 = Unrestrained 1 = Restrained. X and Y DOFs of the nodes are considered.
dirichlet.ValXYZ(length(dirichlet.Front)+length(dirichlet.Right)+1:end,3) = 1; % Restrained Z
dirichlet.ValXYZ(length(dirichlet.Front)+1:length(dirichlet.Front)+length(dirichlet.Right)+length(dirichlet.Back),1) = 1; % Restrained X
dirichlet.ValXYZ(1:length(dirichlet.Front),2) = 1; % Restrained Y
dirichlet.ValXYZ(length(dirichlet.Front)+length(dirichlet.Right)+length(dirichlet.Back)+length(dirichlet.Left)+length(dirichlet.Top)+1:end,1) = 1;
dirichlet.ValXYZ(length(dirichlet.Front)+length(dirichlet.Right)+length(dirichlet.Back)+length(dirichlet.Left)+length(dirichlet.Top)+1:end,2) = 1;


dirichlet.restrainedPts = zeros(size(dirichlet.ValXYZ));% Prescribed values of Contrained Dofs.
dirichlet.restrainedPts(length(dirichlet.Front)+length(dirichlet.Right)+length(dirichlet.Back)+length(dirichlet.Left)+1:...
    length(dirichlet.Front)+length(dirichlet.Right)+length(dirichlet.Back)+length(dirichlet.Left)+length(dirichlet.Top),3) = 1;
dirichlet.reactForce = dirichlet.restrainedPts;

end