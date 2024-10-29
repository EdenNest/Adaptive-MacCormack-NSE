function index = slnn(node,direction,n,ptsnei,ptspos,d) % same_level_nth_neighbor
    % node = index of current node(s)
    % direction : 1 = top, 2 = right , 3 = bottom ; 4 = left
    % n = nth neighbor to that direction to left or right 
    % d shuold be either dx or dy based on the direction (BE CAREFUL)
    
    current_pos = ptspos(node,:);
    level = ptsnei(node,5);
    distance = squeeze(n*d./2.^level);
    % dx = n*dx/2^level;
    % disp('size(current_pos)');size(current_pos)
    % distance;
    switch direction
        case 1                     % use dy
            % size(current_pos);
            % size([zeros(size(distance)),distance])
            new_pos = current_pos + [zeros(size(distance)),distance];
        case 2     % use dx
            new_pos = current_pos + [distance,zeros(size(distance))];
        case 3                    % use dy
            new_pos = current_pos + [zeros(size(distance)),-distance];
        case 4    % use dx
            new_pos = current_pos + [-distance,zeros(size(distance))];
    end

   
    new_pos = squeeze(new_pos);
    if isscalar(node)  
         [~,index] = min(sum((ptspos-new_pos).^2,2));
    else
        index = zeros(size(node));
        for i = 1 : length(node)
            [~,index(i)] = min(sum((ptspos-new_pos(i,:)).^2,2));
        end
    end

end
