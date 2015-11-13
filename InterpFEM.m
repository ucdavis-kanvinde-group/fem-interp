% Vincente Pericoli
% UC Davis
% 4 November 2015
%
% Interpolate element nodal data to other locations in the physical space
%


function [ InterpData, owners ] = InterpFEM ...
           ( points, elementLabels, elemConnect, nodesCoords, LmNodalData )
%For a set of physical coordinates arbitrarily located in the space,
%interpolate element nodal data to those physical coordinates. This is 
%accomplished by first determining which FE element owns each physical 
%point, mapping those physical points into the parent space, and then using 
%the isoparametric shape function relations.
%
%If you are using ABAQUS and need help obtaining LmNodalData, see:
%       https://github.com/ucdavis-kanvinde-group/abaqus-odb-tools
%
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
%   LmNodalData:   (h x k x e) rank-3 matrix of element nodal data which 
%                  you wish to interpolate to the physical points. 
%                  (h,:,:) corresponds to the h-th frame number (history).
%                  (:,:,e) corresponds to element elementLabels(e).
%                  (:,k,:) corresponds to the k-th node in the element, per
%                          ABAQUS node numbering convention
%
%Outputs --
%   InterpData:    (i x 1) cell where the i-th index is (in general) a 
%                  rank-2 matrix of the interpolated data. For example, the
%                  matrix InterpData{i}(h,e) corresponds to the
%                  interpolated data at the h-th frame number for the
%                  element elementLabels(owners{i}(e)). If there is only
%                  one element, then it will be a column vector.
%
%   owners:        (i x 1) cell, where the i-th index corresponds to point
%                  points(i,:). Each cell location is a list of element
%                  labels which own that point (more than 1 element could
%                  own a point)
%


%
% Determine some problem info
%

% use global variable for element type
global LMTYPE

% determine number of points requested
[npts,~] = size(points);

% determine if problem is 3D or 2D
is3D = any(nodesCoords(:,3) ~= 0);

% number of elements, number of nodes per element
[nframe, nnpe, ~] = size(LmNodalData);

% determine the element type
if     (nnpe == 4) && ~is3D
    LMTYPE = 'QUAD4';
elseif (nnpe == 8) && ~is3D
    LMTYPE = 'QUAD8';
elseif (nnpe == 8) && is3D
    LMTYPE = 'BRK8';
elseif (nnpe == 20) && is3D
    LMTYPE = 'BRK20';
else
    error('I haven''t been programmed yet');
end

%
% Begin Interpolation Algorithms
%

% for each physical point, determine which elements "own" that point.
% furthermore, map the physical point into parent space for those owners.
[owners, xi_] = inverseIsoMapping ...
                 ( points, elementLabels, elemConnect, nodesCoords, is3D );

% now that the physical coords are mapped into the parent space, use shape
% functions to interpolate. Use a cell for InterpData because the physical
% point may have multiple results (i.e. if it is shared between elements)
InterpData = cell(npts,1);

for p = 1:npts
    %for all requested points
    for e = 1:length(owners{p})
        %for all element owners of this point
        
        %find owner element index
        eindex = find(elementLabels == owners{p}(e),1);
        
        %evaluate the shape functions in the parent space
        [N, ~, ~] = elementBasisFns( xi_{p}(e,:), LMTYPE );
        
        for h = 1:nframe
            %for all frame values (output history)     
            %use shape functions to evalutate data at point
            %note that data is (1 x n) and N is (n x 1)
            InterpData{p}(h,e) = LmNodalData(h,:,eindex) * N;
        end
    end
end

end

