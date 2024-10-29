function [updated_pts_pos,updated_pts_nei] = refinement_v3(lightup,pts_pos,pts_nei,curr_dx,curr_dy,tot_pts)
%{
Grid Refinement Function

Inputs
lightup: [1 x n_lightup] Array of indexes of cells from coarse grid that
need to be refined
pts_pos: [tot_pts x 2] Array of coordinate positions of each point from
coarse grid (x_coord,y_coord)
pts_nei: [tot_pts x 5] Array of the neighbor indecies of each point from
the coarse grid (North_nei,East_nei,South_nei,West_nei,level)
curr_dx: Value of dx from coarse grid
curr_dy: Value of dy from coarse grid
x_pos: [1 x nx] Array of x coordinate positions of coarse grid
y_pos: [1 x ny] Array of y coordinate positions of coarse grid

Outputs
updated_pts_pos: [tot_pts+n_newnodes x 2] Array of original
coordinate positions of coarse grid with the new fine grid nodes appended
to the end
updated_pts_nei: [tot_pts+n_newnodes x 5] Array of the neighbor indecies
of each point updated to account for finer nodes added
(North_nei,East_nei,South_nei,West_nei,level)
%}
% close all;
%tot_pts = size(pts_pos,1); % Total number of coarse points
n_lightup = length(lightup); % Total number of nodes that need refinement
newpoint_index = 1; % Index counting the number of added finer nodes
hasneigh = zeros(4,1); % Array that records whether the point has a neighbor in each direction or is on the boundary

% Creating finer nodes in the N,E,S,and W directions of each lit up node
for index = 1:n_lightup % Goes through each lit up node
    for dir = 1:4 % Goes through each direction
        if pts_nei(lightup(index),dir) % If the node has a neighbor in this direction
            newpointcoords(newpoint_index,:) = mean([pts_pos(lightup(index),:);pts_pos(pts_nei(lightup(index),dir),:)],1); % New fine point is avg of two points in that direction
            newpoint_index = newpoint_index + 1;
            hasneigh(dir) = 1; % Record that neighbor existed
        else
            hasneigh(dir) = 0; % Record that no neighbor existed
        end
    end
    if sum(hasneigh)==4 % If node has 4 neighbors (center node)
        % Creat new fine points in NE,SE,SW,and NW directions(corner nodes)
        newpointcoords(newpoint_index,:) = pts_pos(lightup(index),:) + [curr_dx/2,curr_dy/2];
        newpointcoords(newpoint_index+1,:) = pts_pos(lightup(index),:) + [-curr_dx/2,curr_dy/2];
        newpointcoords(newpoint_index+2,:) = pts_pos(lightup(index),:) + [curr_dx/2,-curr_dy/2];
        newpointcoords(newpoint_index+3,:) = pts_pos(lightup(index),:) + [-curr_dx/2,-curr_dy/2];
        newpoint_index = newpoint_index + 4;
    elseif sum(hasneigh)==3 % If node has 3 neighbors (edge node)
        boundarydir = find(pts_nei(lightup(index),1:4)==0); % Determine the direction towards the boundary
        if boundarydir == 1 % If boundary is to the North
            newpointcoords(newpoint_index,:) = pts_pos(lightup(index),:) + [curr_dx/2,-curr_dy/2];
            newpointcoords(newpoint_index+1,:) = pts_pos(lightup(index),:) + [-curr_dx/2,-curr_dy/2];
            newpoint_index = newpoint_index + 2;
        elseif boundarydir == 2 % If boundary is to the East
            newpointcoords(newpoint_index,:) = pts_pos(lightup(index),:) + [-curr_dx/2,curr_dy/2];
            newpointcoords(newpoint_index+1,:) = pts_pos(lightup(index),:) + [-curr_dx/2,-curr_dy/2];
            newpoint_index = newpoint_index + 2;
        elseif boundarydir == 3 % If boundary is to the South
            newpointcoords(newpoint_index,:) = pts_pos(lightup(index),:) + [curr_dx/2,curr_dy/2];
            newpointcoords(newpoint_index+1,:) = pts_pos(lightup(index),:) + [-curr_dx/2,curr_dy/2];
            newpoint_index = newpoint_index + 2;
        else % If boundary is to the West
            newpointcoords(newpoint_index,:) = pts_pos(lightup(index),:) + [curr_dx/2,curr_dy/2];
            newpointcoords(newpoint_index+1,:) = pts_pos(lightup(index),:) + [curr_dx/2,-curr_dy/2];
            newpoint_index = newpoint_index + 2;
        end
    elseif sum(hasneigh)==2 % If node has 2 neighbors (corner node)
        if pts_nei(lightup(index),1)&&pts_nei(lightup(index),2) % If bottom-left corner
            newpointcoords(newpoint_index,:) = pts_pos(lightup(index),:) + [curr_dx/2,curr_dy/2];
            newpoint_index = newpoint_index + 1;
        elseif pts_nei(lightup(index),2)&&pts_nei(lightup(index),3) % If top-left corner
            newpointcoords(newpoint_index,:) = pts_pos(lightup(index),:) + [curr_dx/2,-curr_dy/2];
            newpoint_index = newpoint_index + 1;
        elseif pts_nei(lightup(index),3)&&pts_nei(lightup(index),4) % If top-right corner
            newpointcoords(newpoint_index,:) = pts_pos(lightup(index),:) + [-curr_dx/2,-curr_dy/2];
            newpoint_index = newpoint_index + 1;
        else % If bottom-right corner
            newpointcoords(newpoint_index,:) = pts_pos(lightup(index),:) + [-curr_dx/2,curr_dy/2];
            newpoint_index = newpoint_index + 1;
        end
    end
end
% Elminate repeated new finer nodes
newpointcoords = unique(newpointcoords,'rows');
% Collect all fine nodes
finepointcoords = cat(1,pts_pos(tot_pts+1:end,:),newpointcoords);
% Eliminate fine node repeats due to computational inaccuracies
elim_finepts = zeros(size(finepointcoords,1));
for index = 1:size(finepointcoords,1) % For each new node
    % Subtract new node coords. from all new nodes and search for very
    % small differences
    smartsearch = abs(finepointcoords - finepointcoords(index,:));
    condition = find(smartsearch(:,1)<curr_dx/4&smartsearch(:,2)<curr_dy/4); % Record repeated nodes (small differences)
    if size(condition,1)>1 % If more than 1 new node has similar coords.
        elim_finepts(1:size(condition,1)-1,index) = condition(2:end); % Record indecies of repeated nodes to be deleted
    end
end
elim_finepts = unique(nonzeros(elim_finepts)); % Omit zeros and find unique indecies to be deleted
finepointcoords(elim_finepts,:) = []; % Delete repeats due to computational inaccuracies
uniq_finepoints = finepointcoords; % Organize and check for uniqueness once more

updated_pts_pos = cat(1,pts_pos(1:tot_pts,:),uniq_finepoints); % Append new finer node coords. to coarse grid coords.
updated_pts_nei = cat(1,pts_nei(1:tot_pts,:),zeros(size(uniq_finepoints,1),5)); % Append empty neighbor array to be filled with fine node neighbors

% Update neighbor indecies of lit up nodes
for index = 1:n_lightup % Goes through each lit up node
    for dir = 1:4 % Goes through each direction
        if updated_pts_nei(lightup(index),dir) % If the node has a neighbor in this direction
            findnewnode = mean([pts_pos(lightup(index),:);pts_pos(updated_pts_nei(lightup(index),dir),:)],1); % Get coords. of fine node that would be inbetween
            smartsearch = abs(updated_pts_pos - findnewnode); % Find index of node accounting for computational inaccuracies
            updated_pts_nei(lightup(index),dir) = find(smartsearch(:,1)<curr_dx/4&smartsearch(:,2)<curr_dy/4); % Find index of fine node with these coords
            updated_pts_nei(lightup(index),5) = 1; % Fill columns 5 of lit up nodes with level value of 1 (finer grid due to new neighbors)
        end
    end
end
% Fill in neighbor indecies of new finer nodes
for index = tot_pts+1:size(updated_pts_nei,1) % Goes through all fine nodes(starts at end of coarse section)
    if updated_pts_nei(index,1)==0 % If north neighbor index hasn't been filled yet
        pot_neighbor = updated_pts_pos(index,:) + [0,curr_dy/2]; % Calculate coords. of potential neighbor
        smartsearch = abs(updated_pts_pos - pot_neighbor); % Find index of neighbor accounting for computational inaccuracies
        if sum(find(smartsearch(:,1)<curr_dx/4&smartsearch(:,2)<curr_dy/4))==0 % If no north neighbor exists
            % This fine node is on the edge of the fine mesh, keep 0 to denote boundary
        else % If north neighbor does exist
            updated_pts_nei(index,1) = find(smartsearch(:,1)<curr_dx/4&smartsearch(:,2)<curr_dy/4); % Find neighbor index and record
        end
    end
    if updated_pts_nei(index,2)==0 % If east neighbor index hasn't been filled yet
        pot_neighbor = updated_pts_pos(index,:) + [curr_dx/2,0]; % Calculate coords. of potential neighbor
        smartsearch = abs(updated_pts_pos - pot_neighbor); % Find index of neighbor accounting for computational inaccuracies
        if sum(find(smartsearch(:,1)<curr_dx/4&smartsearch(:,2)<curr_dy/4))==0 % If no east neighbor exists
            % This fine node is on the edge of the fine mesh, keep 0 to denote boundary
        else % If east neighbor does exist
            updated_pts_nei(index,2) = find(smartsearch(:,1)<curr_dx/4&smartsearch(:,2)<curr_dy/4); % Find neighbor index and record
        end
    end
    if updated_pts_nei(index,3)==0 % If south neighbor index hasn't been filled yet
        pot_neighbor = updated_pts_pos(index,:) + [0,-curr_dy/2]; % Calculate coords. of potential neighbor
        smartsearch = abs(updated_pts_pos - pot_neighbor); % Find index of neighbor accounting for computational inaccuracies
        if sum(find(smartsearch(:,1)<curr_dx/4&smartsearch(:,2)<curr_dy/4))==0 % If no south neighbor exists
            % This fine node is on the edge of the fine mesh, keep 0 to denote boundary
        else % If south neighbor does exist
            updated_pts_nei(index,3) = find(smartsearch(:,1)<curr_dx/4&smartsearch(:,2)<curr_dy/4); % Find neighbor index and record
        end
    end
    if updated_pts_nei(index,4)==0 % If west neighbor index hasn't been filled yet
        pot_neighbor = updated_pts_pos(index,:) + [-curr_dx/2,0]; % Calculate coords. of potential neighbor
        smartsearch = abs(updated_pts_pos - pot_neighbor); % Find index of neighbor accounting for computational inaccuracies
        if sum(find(smartsearch(:,1)<curr_dx/4&smartsearch(:,2)<curr_dy/4))==0 % If no west neighbor exists
            % This fine node is on the edge of the fine mesh, keep 0 to denote boundary
        else % If west neighbor does exist
            updated_pts_nei(index,4) = find(smartsearch(:,1)<curr_dx/4&smartsearch(:,2)<curr_dy/4); % Find neighbor index and record
        end
    end
end
% Fill column 5 of fine nodes with level value of 1 (finer grid)
updated_pts_nei(tot_pts+1:end,5) = 1;

% % Plot figure to show completeness
% figure(3)
% scatter(pts_pos(1:tot_pts,1),pts_pos(1:tot_pts,2),[],'k') % Black circles to mark coarse grid points
% prev_refined = find(pts_nei(1:tot_pts,5)==1);
% grid on
% hold on
% scatter(pts_pos(prev_refined,1),pts_pos(prev_refined,2),[],'k','filled') % Filled black circles to mark previously refined nodes
% scatter(pts_pos(lightup,1),pts_pos(lightup,2),[],'r','filled') % Filled red circles to mark new lit up nodes
% scatter(uniq_finepoints(:,1),uniq_finepoints(:,2),[],'b','filled') % Blue circles to mark fine nodes
