function FirstDerivative = ddy_central(ptsnei,ptspos,f,dy)

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
            % secondtopnei = ptsnei(ptsnei(i,1),1);
            secondtopnei = slnn(i,1,2,ptsnei,ptspos,dy);
            FirstDerivative(i,j) =...
                (-3*ff(i) +4*ff(ptsnei(i,1)) -ff(secondtopnei)) * 2^lvl / dy / 2;

        elseif ptsnei(i,1) == 0 % top border
            % secondbottomnei = ptsnei(ptsnei(i,3),3);
            secondbottomnei = slnn(i,3,2,ptsnei,ptspos,dy);
            FirstDerivative(i,j) =...
                (3*ff(i) -4*ff(ptsnei(i,3)) + ff(secondbottomnei)) *  2^lvl / dy / 2;

        else
            FirstDerivative(i,j) = ...
                (ff(ptsnei(i,1))-ff(ptsnei(i,3))) *  2^lvl  / dy / 2;
        end


    end

end


if a<b
    FirstDerivative=FirstDerivative'; %6000 x 4
end
end

%% mine
%
% [nx,ny]=size(f);
%
% A=diag(0*ones(1,ny)) +...
%     diag(1*ones(1,ny-1),1) +...
%     diag(-1*ones(1,ny-1),-1);
% A(1,1:3)=[-3,4,-1];
% A(end,ny-2:end)=[1,-4,3];
% % A(1,end)=-1;
% % A(end,1)=1;
% A=sparse(A);
% A=A/2;
%
% % dfdy=permute((1/dy*A*permute(f,[2,1])),[2,1]);
% dfdy=transpose((1/dy*A*f'));
