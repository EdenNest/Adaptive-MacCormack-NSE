
function FirstDerivative=ddy_bwd(ptsnei,ptspos,f,dy)

% f is a column matrix having nodes as rows

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

        if ptsnei(i,3) == 0 % bottom border
            if ismember(3,refinedborder(i,ptsnei,ptspos))
                % ghost = interpGhost(i,3,f,ptspos,dy);
                ghost = averageghost(i,ff,3,ptspos);
                FirstDerivative(i,j) = (ff(ptsnei(i,1))-ghost) * 2^lvl / dy / 2; 
           % k = slnn(i,1,2,ptsnei,ptspos,dy);
           %      FirstDerivative(i,j) = (-ff(k)+4*ff(ptsnei(i,1))-3*ff(i))* 2^lvl / dy / 2; 
               
            else
            FirstDerivative(i,j) = (ff(ptsnei(i,1))-ff(i)) * 2^lvl / dy;
            end
        else
            FirstDerivative(i,j) = (-ff(ptsnei(i,3))+ff(i)) * 2^lvl  / dy;
        end
    end



end

    if a<b
        FirstDerivative=FirstDerivative'; %6000 x 4
    end

%
% [nx,ny]=size(f);
% % first order backward
% A=diag(1*ones(1,ny)) +...
%     diag(-1*ones(1,ny-1),-1) +...
%     diag(0*ones(1,ny-2),-2);
% A(1,1:2)=[-1,1];
% A=sparse(A);
% f=double(f);
%
% % dfdy=permute((1/dy*A*permute(f,[2,1])),[2,1]);
%  dfdy=transpose((1/dy*A*f'));
%
% second order backward
% A=diag(3*ones(1,ny)) +...
%     diag(-4*ones(1,ny-1),-1) +...
%     diag(1*ones(1,ny-2),-2);
% A(1,1:3)=[-3,4,-1];
% A(2,1:4)=[0,-3,4,-1];
% A=A/2;

end


