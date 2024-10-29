function index = FindIt(ptspos,x_coor,y_coor)
% if there exists a node on (x_coor,y_coor) it will give you the index, if
% not it will give you zero
    
    A = sum(abs(ptspos - [x_coor,y_coor]),2) ; 
    
    index = find(A<1e-8);
    
        if isempty(index)
            index = 0;
        end

end


