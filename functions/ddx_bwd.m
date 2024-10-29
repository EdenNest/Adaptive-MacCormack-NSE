function FirstDerivative=ddx_bwd(ptsnei,ptspos,f,dx)
% [ 4 * 6000]



npts = length(ptsnei);
[a,b]   = size(f);
if a<b
    f=f'; %6000 x 4
end
FirstDerivative = zeros(size(f));


for i=1:npts
    lvl = ptsnei(i,5);
    
    for j = 1:size(f,2)
        ff = f(:,j);
        if ptsnei(i,4) == 0 % left border
            % check if is a border on a refined zone therefore doing cetral
            if ismember(4,refinedborder(i,ptsnei,ptspos))
                % ghost = interpGhost(i,4,f,ptspos,dx); 
                ghost = averageghost(i,ff,4,ptspos);
                FirstDerivative(i,j) = (ff(ptsnei(i,2))-ghost) * 2^lvl / dx / 2; 
                
                % k = slnn(i,2,2,ptsnei,ptspos,dx);
                % FirstDerivative(i,j) = (+ff(k)-4*ff(ptsnei(i,2))+3*ff(i))* 2^lvl / dx / 2; 
               
            else
                FirstDerivative(i,j) = (ff(ptsnei(i,2))-ff(i)) * 2^lvl / dx;
            end
        else
            FirstDerivative(i,j) = (-ff(ptsnei(i,4))+ff(i)) * 2^lvl  / dx;
        end
    end
end



if a<b
    FirstDerivative=FirstDerivative'; %6000 x 4
end


end
