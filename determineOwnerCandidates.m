% Vincente Pericoli
% UC Davis
% 1 November 2015

function [ possibleOwners ]  = determineOwnerCandidates ...
                  ( points, elementLabels, elemConnect, nodesCoords, is3D )
%Given mesh information and desired point locations, determines a small set 
%of "candidate" elements which could possibly own that point (for each
%individual desired point). This is a computationally fast and efficient
%algorithm, especially for a very large set of points (say, greater than 
%10,000)
%
%Partially based on:
%   Silva, G. H., Le Riche, R., Molimard, J., & Vautrin, A. (2009). 
%   "Exact and efficient interpolation using finite elements shape 
%   functions." European Journal of Computational Mechanics/Revue 
%   Européenne de Mécanique Numérique, 18(3-4), 307-331.
%
%Specifically, the use of a "virtual mesh", as described in the paper, is
%implemented here. This provides a significant improvement in efficiency.
%For example, my own tests have shown that, for a dense and complicated 2D
%mesh, the virtual mesh provided a 77% improvement in computation time
%(reducing the time from 352s to 80s for a 10,000 point search!). The 
%efficiency is also mathematically proven in the above paper. Please see 
%the paper for more info and details on the algorithm.
%
%However, only the "preliminary" checks (virtual mesh and bounding box) are
%used here. The "Scalar Cross Product" test from the above paper (2D), and 
%other known/standard tests based on Delaunay triangulation using Qhull 
%(3D), are only guaranteed to work if all elements are a convex hull. 
%Serendipity elements (2D) or 20 node bricks (3D) will NOT necessarily be 
%convex hulls. Therefore, this code performs only rudimentary checks to 
%cull the list of possible owners. The exact ownership can then be
%identified with an inverse mapping from the physical coordinates into the 
%parent coordinates.
%
%
%Inputs --
%   points:        (i x 3) matrix of desired point locations. Each row i
%                  corresponds to a point, while columns are the x,y,z 
%                  coordinates
%
%   elementLabels: the element numbers corresponding to elemConnect
%
%   elemConnect:   (e x n) matrix where each row is an element, and each
%                  column is the nodal connectivity of that element. Row e
%                  corresponds to element elementLabels(e)
%
%   nodesCoords:   (n x 3) matrix where each row is a node, and each column
%                  is the nodal coordinates. Row n corresponds to node n.
%
%   is3D:          Boolean flag, designating the dimensionality of the
%                  input mesh. True (or 1) indicates a 3D mesh, while 
%                  False (or 0) indicates a 2D mesh
%
%Outputs --
%   possibleOwners: (i x 1) cell, where the i-th index corresponds to a 
%                   list of the elementLabels of the possible owners for 
%                   point i
%

% define dimensionality 
if is3D
    ndim = 3;
else
    ndim = 2;
end

%
% determine size of problem
%

%number of points to find the owner of
[npts,~] = size(points);

%number of elements, and number of nodes per element
[numel, ~] = size(elemConnect);

%preallocation
possibleOwners = cell(npts,1);

%
% Build a Virtual Mesh
%

% crude approximation of number of virtual elements, see Silva et al.(2009)
numve = numel/(sqrt(numel/1.25/npts) + 1);

% number of divisions for each dimension. this assumes every dimension has
% equal numbers of virual elements
ndiv = round(numve/ndim);

% initialize virtual mesh
if is3D
    VirtualMesh = cell(ndiv,ndiv,ndiv);
else
    VirtualMesh = cell(ndiv,ndiv);
end

% determine VM boundary
originVM = min(nodesCoords(:,1:ndim));
dxi = ( max(nodesCoords(:,1:ndim)) - originVM )./ndiv;

% add FE mesh into VM
if is3D
    for e = 1:numel
        %for each FE in the mesh, determine the corresponding VM indices
        imin = floor( ( min(nodesCoords(elemConnect(e,:),1:ndim)) ...
                       - originVM )./dxi ) + 1;
        imax = floor( ( max(nodesCoords(elemConnect(e,:),1:ndim)) ...
                       - originVM )./dxi ) + 1;
        
        %imax has the potential to exceed the dimensions of the problem
        imax(imax > ndiv) = ndiv;
        
        %add this element into those indices
        for xi = imin(1):imax(1)
            for yi = imin(2):imax(2)
                for zi = imin(3):imax(3)
                    VirtualMesh{xi,yi,zi}(end+1) = elementLabels(e);
                end
            end
        end
    end
else
    for e = 1:numel
        %for each FE in the mesh, determine the corresponding VM indices
        imin = floor( ( min(nodesCoords(elemConnect(e,:),1:ndim)) ...
                        - originVM )./dxi ) + 1;
        imax = floor( ( max(nodesCoords(elemConnect(e,:),1:ndim)) ...
                        - originVM )./dxi ) + 1;
        
        %imax has the potential to exceed the dimensions of the problem
        imax(imax > ndiv) = ndiv;
                    
        %add this element into those indices
        for xi = imin(1):imax(1)
            for yi = imin(2):imax(2)
                VirtualMesh{xi,yi}(end+1) = elementLabels(e);
            end
        end
    end
end

% for all requested physical points, obtain the possible owners via
% boundingBoxTest
if is3D
    for p = 1:npts
        %point index in VM is:
        pi = floor( ( points(p,1:ndim) - originVM )./dxi ) + 1;
        
        %pi has the potential to exceed the dimensions of the problem
        pi(pi > ndiv) = ndiv;
        
        %elements in that corresponding VM are
        ptElemList = VirtualMesh{pi(1),pi(2),pi(3)};
        
        %for thsoe elements, perform bounding box test
        possibleOwners(p) = boundingBoxTest ...
                                ( points(p,:), ptElemList, ...
                                  elemConnect(ptElemList,:), nodesCoords );
    end
else
    for p = 1:npts
        %point index in VM is:
        pi = floor( ( points(p,1:ndim) - originVM )./dxi ) + 1;
        
        %pi has the potential to exceed the dimensions of the problem
        pi(pi > ndiv) = ndiv;
        
        %elements in that corresponding VM are
        ptElemList = VirtualMesh{pi(1),pi(2)};

        %for thsoe elements, perform bounding box test. this returns cells
        possibleOwners(p) = boundingBoxTest ...
                                ( points(p,:), ptElemList, ...
                                  elemConnect(ptElemList,:), nodesCoords );
    end
end

end

function [ possibleOwners ]  = boundingBoxTest ...
                        ( points, elementLabels, elemConnect, nodesCoords )
%Perform bounding box test to cull list of possible owners to few. This was
%written as the original code (i.e. before the virtual mesh was written).
%Therefore, it is written to be general-purpose, and is not dependent on
%the virtual mesh whatsoever.
%

%number of points to find the owner of
[npts,~] = size(points);

%number of elements, and number of nodes per element
[numel, nnpe] = size(elemConnect);


% ASSEMBLE BOUNDING BOXES
BoundingBox_min = zeros(numel,3);
BoundingBox_max = zeros(numel,3);
tempCoordList   = zeros(1,nnpe);

for e = 1:numel
    %for all elements
	for i = 1:3
        %for all coordinate directions
        for n = 1:nnpe
            %for all nodes in element
            node = elemConnect(e,n);
            tempCoordList(n) = nodesCoords(node,i);
        end
        
        BoundingBox_min(e,i) = min(tempCoordList);
        BoundingBox_max(e,i) = max(tempCoordList);
	end
end


% TEST BOUNDING BOXES, CULL POSSIBLE OWNERSHIP
possibleOwners = cell(npts,1);
 
for p = 1:npts
    %for all points to find ownership of
    possibleOwners{p} = []; %define list
    for e = 1:numel
        %for all elements
        if all(BoundingBox_min(e,:) <= points(p,:)) && ...
           all(points(p,:) <= BoundingBox_max(e,:))
            %then the this element is a possible owner
            possibleOwners{p}(end+1) = elementLabels(e); %append to list
        end
    end
end

end
