% Vincente Pericoli
% UC Davis
% 4 November 2015

function [ owners, xi_ ] = inverseIsoMapping ...
          ( points, elementLabels, elemConnect, nodesCoords, ndim, LMTYPE )
%Given a set of physical points in the space, determine the FE elements who
%own those physical points, and perform an inverse isoparametric map into
%the parent space.
%
%Inputs --
%   points:        (i x 3) matrix of desired physical point locations. Each
%                  row i corresponds to a point, while columns are the
%                  x,y,z coordinates
%
%   elementLabels: the element numbers corresponding to elemConnect
%
%   elemConnect:   (e x n) matrix where each row is an element, and each
%                  column is the nodal connectivity of that element. Row e
%                  corresponds to element elementLabels(e)
%
%   nodesCoords:   (n x 3) matrix where each row is a node, and each column
%                  is the nodal coordinates. Row n corresponds to node n
%                  (i.e. nodesCoords has all nodes in the physical space)
%
%   ndim:          number of dimensions of the input mesh
%
%   LMTYPE:        string name of the (standard) element type, e.g. 'QUAD4'
%
%Outputs --
%   owners:       (i x 1) cell, where the i-th index corresponds to point
%                 points(i,:). Each cell location is a list of element
%                 labels which own that point (more than 1 element could
%                 own a point)
%
%   xi_           (i x 1) cell, where the i-th index corresponds to point
%                 points(i,:). Each cell location is a matrix of parent
%                 coordinates. For example, xi{i}(1,:) corresponds to the
%                 parent coordinates of the owner element owners{i}(1)
%

% check args
assert((ndim == 2) || (ndim == 3));

% determine number of points requested
[npts,~] = size(points);

% preallocate owners 
owners = cell(npts,1); % initialize with empty cell array

% preallocate xi_ vector (parent coordinates for points)
xi_ = cell(npts,1); % initialize with empty cell array

% for each point, cull possible element ownership to a small subset
possibleOwners = determineOwnerCandidates ...
                 ( points, elementLabels, elemConnect, nodesCoords, ndim );

% define tolerance boundaries using an absolute tolerance. the use of an
% absolute tolerance is acceptable because the exact boundary is 1 or -1
upper =  1 + 1e-4;
lower = -1 - 1e-4;
             
% for each point and its corresponding possible owners, perform an inverse
% isoparametric mapping
for p = 1:npts
    %for all points to inverse map
    
    %number of possible owner elements for this point
    nele = length(possibleOwners{p});
    
    for e = 1:nele
        %for all possible owners
        
        %identify element index
        eind = find( elementLabels == possibleOwners{p}(e) ,1);
        
        %therefore, the coordinates of the nodes for this element are:
        %nodesCoords(elemConnect(eind,:),:)
        
        %attempt inverse mapping using optimization algorithm
        [e_xi_, eflag] = inverseIsoObjFun( points(p,:), ...
                              nodesCoords(elemConnect(eind,:),:), LMTYPE );
        
        %check to see if the inverse mapping was successful. use single
        %precision "buffer", because ABAQUS coords are single precision
        if (eflag > 0) && all(e_xi_ <= upper) && all(e_xi_ >= lower)
            %then the inverse mapping was successful
            %and, this is the owner element of the point!
            
            %since a physical point may be shared between elements (like at
            %nodal locations), we want to save to a cell structure.
            owners{p}(end+1) = elementLabels(eind);
            xi_{p}(end+1,:)  = e_xi_;
            
            %don't break--we must loop through all possible owners because
            %the point may be shared.
        end
    end
end

return;
end