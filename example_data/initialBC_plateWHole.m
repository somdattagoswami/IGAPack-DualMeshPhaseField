function [dirichlet,cornerElem] = initialBC_plateWHole(PHTelem,geometry)
% Initialize the boundary conditions

p = geometry.p;
q = geometry.q;

down_nodes = 1:p+1;
up_nodes = 1+(p+1)*q:(p+1)*(q+1);
cornerElem = [];
simpleSupport = [];

% Set the fixed boundary degrees of freedom
for patchIndex = 5
    for i=1:length(PHTelem{patchIndex})
        if isempty(PHTelem{patchIndex}(i).children)
            if isempty(PHTelem{patchIndex}(i).neighbor_down) && isempty(PHTelem{patchIndex}(i).neighbor_left)
                bottomEdge = PHTelem{patchIndex}(i).nodesGlobal(down_nodes);
                simpleSupport = bottomEdge(1);
                cornerElem = [cornerElem; patchIndex, i];
                break;
            end
        end
    end
end

for patchIndex = 1
    for i=1:length(PHTelem{patchIndex})
        if isempty(PHTelem{patchIndex}(i).children)
            if isempty(PHTelem{patchIndex}(i).neighbor_down) && isempty(PHTelem{patchIndex}(i).neighbor_right)
                cornerElem = [cornerElem; patchIndex, i];
                break;
            end
        end
    end
end

% Set the fixed boundary degrees of freedom
for patchIndex = 29
    for i=1:length(PHTelem{patchIndex})
        if isempty(PHTelem{patchIndex}(i).children)
            if isempty(PHTelem{patchIndex}(i).neighbor_down) && isempty(PHTelem{patchIndex}(i).neighbor_right)
                bottomEdge = PHTelem{patchIndex}(i).nodesGlobal(down_nodes);
                rollerSupport = bottomEdge(end);
                cornerElem = [cornerElem; patchIndex, i];
                break;
            end
        end
    end
end

% Set the fixed boundary degrees of freedom
for patchIndex = 33
    for i=1:length(PHTelem{patchIndex})
        if isempty(PHTelem{patchIndex}(i).children)
            if isempty(PHTelem{patchIndex}(i).neighbor_down) && isempty(PHTelem{patchIndex}(i).neighbor_left)
                cornerElem = [cornerElem; patchIndex, i];
                break;
            end
        end
    end
end

% Set the loading boundary degrees of freedom
for patchIndex = 26
    for i=1:length(PHTelem{patchIndex})
        if isempty(PHTelem{patchIndex}(i).children)
            if isempty(PHTelem{patchIndex}(i).neighbor_up)&& isempty(PHTelem{patchIndex}(i).neighbor_left)
                topEdge = PHTelem{patchIndex}(i).nodesGlobal(up_nodes);
                loadingPt = topEdge(1);
                break;
            end
        end
    end
end

for patchIndex = [1,22,26]
    for i=1:length(PHTelem{patchIndex})
        if isempty(PHTelem{patchIndex}(i).children)
            cornerElem = [cornerElem; patchIndex, i];
        end
    end
end


dirichlet.XY =[simpleSupport,rollerSupport,loadingPt]; % Contains the node numbers whose DOF/s have to be restrained.
dirichlet.ValXY = zeros(length(dirichlet.XY'),2); %  0 = Unrestrained 1 = Restrained. X and Y DOFs of the nodes are considered.
dirichlet.ValXY (1,1:2) = 1; % Restrained X and Y
dirichlet.ValXY (2:3,2) = 1; % Restrained Y

dirichlet.restrainedPts = zeros(size(dirichlet.ValXY));% Prescribed values of Contrained Dofs.
dirichlet.restrainedPts(3,2) = 1;
dirichlet.reactForce = dirichlet.restrainedPts;
end