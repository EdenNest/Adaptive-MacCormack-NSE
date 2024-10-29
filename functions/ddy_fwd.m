function FirstDerivative = ddy_fwd(ptsnei,ptspos,f,dy)


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

        if ptsnei(i,1) == 0 % top border
            if ismember(1,refinedborder(i,ptsnei,ptspos))
                % ghost = interpGhost(i,1,f,ptspos,dy);
                ghost = averageghost(i,ff,1,ptspos);
                FirstDerivative(i,j) = (-ff(ptsnei(i,3))+ghost) * 2^lvl / dy / 2;
                % k = slnn(i,3,2,ptsnei,ptspos,dy);
                % FirstDerivative(i,j) = -(-ff(k)+4*ff(ptsnei(i,3))-3*ff(i))* 2^lvl / dy / 2; 
               
            else

                FirstDerivative(i,j) = (-ff(ptsnei(i,3))+ff(i)) * 2^lvl  / dy;
            end
        else
            FirstDerivative(i,j) = (ff(ptsnei(i,1))-ff(i)) * 2^lvl / dy;

        end
    end




end


if a<b
    FirstDerivative=FirstDerivative';
end

end



%
% % FIRST ORDER
% [nx,ny]     = size(f);
%   A=diag(-1*ones(1,ny)) +...
%         diag(1*ones(1,ny-1),1) +...
%         diag(0*ones(1,ny-2),2);
%     A(end,ny-2:end)=[1,-4,3];
%
%
% A=sparse(A);
% dfdy=transpose((1/dy*A*f'));