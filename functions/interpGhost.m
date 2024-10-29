function output = interpGhost(index,direction,field,ptspos,distance)
F = scatteredInterpolant(ptspos(:,1),ptspos(:,2),field);
current_pos = ptspos(index,:);
switch direction

    case 1                     % use dy
        new_pos = current_pos + [0,distance];
    case 2     % use dx
        new_pos = current_pos + [distance,0];
    case 3                    % use dy
        new_pos = current_pos + [0,-distance];
    case 4    % use dx
        new_pos = current_pos + [-distance,0];
end

output = F(new_pos);
end