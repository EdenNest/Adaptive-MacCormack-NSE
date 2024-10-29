function CorrectNeighbors(ptsnei,ptspos)

% Finds all the zero neighbors in ptsnei, then finds the positions of these
% zero neighbors, then checks if there's a point in this position (maybe a
% newly added point), if yes, replaces the zero of ptsnei with the index 
% of this point that it just found with the help of FindIt

dy = 8e-6/79; dx = 1e-5/74;

[index,direction] = find(~ptsnei(:,1:4)); % find the zeros in ptsnei



for i=1:length(index)
  
  switch direction(i)   % find the position of the zero neighbor
      case 1
          k =  [0,dy/2];
      case 2
          k =  [dx/2,0];
      case 3
          k =  [0,-dy/2];
      case 4
          k =  [-dx/2,0];
  end

    nei_pos = ptspos(index , : ) + k ;

    
  % corrects the position of the zero neighbor
ptsnei(index(i),direction(i)) = FindIt(ptspos ,  nei_pos(1) , nei_pos(2));

end


end