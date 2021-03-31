function [octupleRef] = refineMeshCoordinates3D(PHTelem,controlPts,octupleList,geometry,xmin,xmax,ymin,ymax,zmin,zmax)
%marks for refinement all the elements on the PHTelem mesh which are inside the
%rectangle (xmin,xmax) x (ymin,ymax)

% Plots the elements stored in PHTelem structure array
p = geometry.p;
q = geometry.q;
r = geometry.r;

numPts = 2; %number of evaluation points to use on each edge
uref = linspace(-1,1,numPts);
vref = linspace(-1,1,numPts);
wref = linspace(-1,1,numPts);

% 1D bernstein polynomials evaluated at the Gauss points on the master element
[B_u] = bernstein_basis(uref,p);
[B_v] = bernstein_basis(vref,q);
[B_w] = bernstein_basis(wref,r);

R = zeros(numPts,numPts,numPts,(p+1)*(q+1)*(r+1));
basisCounter = 0;

for k=1:r+1
    for j=1:q+1
        for i=1:p+1
            basisCounter = basisCounter + 1;
            for kk=1:numPts
                for jj=1:numPts
                    for ii=1:numPts
                        R(ii,jj,kk,basisCounter) = B_u(ii,i)*B_v(jj,j)*B_w(kk,k);
                    end
                end
            end
        end
    end
end

numOctuple = size(octupleList,1);

% Initialize the marked patch array to zeros (i.e. no refinement)
octupleRef = zeros(1, numOctuple);


for i=1:numOctuple
    for j=1:8
        % Determine the edges of the element in the physical space
        coord_up_left = zeros(numPts,3);
        coord_down_left = zeros(numPts,3);
        coord_up_right = zeros(numPts,3);
        coord_down_right = zeros(numPts,3);
        coord_up_front = zeros(numPts,3);
        coord_down_front = zeros(numPts,3);
        coord_up_back = zeros(numPts,3);
        coord_down_back = zeros(numPts,3);
        coord_left_front = zeros(numPts,3);
        coord_right_front = zeros(numPts,3);
        coord_left_back = zeros(numPts,3);
        coord_right_back = zeros(numPts,3);
        
        
        
        curElem = octupleList(i,j);
        nument = min(size(PHTelem(curElem).C,1),length(PHTelem(i).nodes));
        nodes = PHTelem(curElem).nodes(1:nument);
        cpts = controlPts(nodes,1:3);
                
        for jj=1:numPts
                      
            UpLeft = (PHTelem(curElem).C(1:nument,:))*squeeze(R(1,jj,numPts,:));
            DownLeft = (PHTelem(curElem).C(1:nument,:))*squeeze(R(1,jj,1,:));
            UpRight = (PHTelem(curElem).C(1:nument,:))*squeeze(R(numPts,jj,numPts,:));
            DownRight = (PHTelem(curElem).C(1:nument,:))*squeeze(R(numPts,jj,1,:));
            UpFront = (PHTelem(curElem).C(1:nument,:))*squeeze(R(jj,1,numPts,:));
            DownFront = (PHTelem(curElem).C(1:nument,:))*squeeze(R(jj,1,1,:));
            UpBack=(PHTelem(curElem).C(1:nument,:))*squeeze(R(jj,numPts,numPts,:));
            DownBack=(PHTelem(curElem).C(1:nument,:))*squeeze(R(jj,numPts,1,:));
            LeftFront=(PHTelem(curElem).C(1:nument,:))*squeeze(R(1,1,jj,:));
            RightFront=(PHTelem(curElem).C(1:nument,:))*squeeze(R(numPts,1,jj,:));
            LeftBack=(PHTelem(curElem).C(1:nument,:))*squeeze(R(1,numPts,jj,:));
            RightBack=(PHTelem(curElem).C(1:nument,:))*squeeze(R(numPts,numPts,jj,:));
            
            coord_up_left(jj,:) = UpLeft'*cpts;
            coord_down_left(jj,:) = DownLeft'*cpts;
            coord_up_right(jj,:) = UpRight'*cpts;
            coord_down_right(jj,:) = DownRight'*cpts;
            coord_up_front(jj,:) = UpFront'*cpts;
            coord_down_front(jj,:) = DownFront'*cpts;
            coord_up_back(jj,:) = UpBack'*cpts;
            coord_down_back(jj,:) = DownBack'*cpts;
            coord_left_front(jj,:) = LeftFront'*cpts;
            coord_right_front(jj,:) = RightFront'*cpts;
            coord_left_back(jj,:) = LeftBack'*cpts;
            coord_right_back(jj,:) = RightBack'*cpts;
            
        end
        coords = [coord_up_left;coord_down_left;coord_up_right;coord_down_right;coord_up_front;coord_down_front;
                  coord_up_back;coord_down_back;coord_left_front;coord_right_front;coord_left_back;coord_right_back];
        
        %check if the coords are inside the marked region
        for iCoord = 1:size(coords,1)
            xCoord = coords(iCoord,1);
            yCoord = coords(iCoord,2);
            zCoord = coords(iCoord,3);
            if (xmin <= xCoord) && (xCoord <= xmax) && (ymin <= yCoord) && (yCoord <= ymax)&& (zmin <= zCoord) && (zCoord <= zmax)
                octupleRef(i) = 1;
                break;
            end
            
        end
        
    end
    
end

