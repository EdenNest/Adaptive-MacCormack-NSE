function FirstDerivative=ddx_fwd(ptsnei,ptspos,f,dx)



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
            if ptsnei(i,2) == 0 % right border
                if ismember(2,refinedborder(i,ptsnei,ptspos))
                % ghost = interpGhost(i,2,f,ptspos,dx);
                ghost = averageghost(i,ff,2,ptspos);
                FirstDerivative(i,j) = (-ff(ptsnei(i,4))+ghost) * 2^lvl / dx / 2; 
                % k = slnn(i,4,2,ptsnei,ptspos,dx);
                % FirstDerivative(i,j) = -(-ff(k)+4*ff(ptsnei(i,4))-3*ff(i))* 2^lvl / dx / 2; 
               
                else
                FirstDerivative(i,j) = (-ff(ptsnei(i,4))+ff(i)) * 2^lvl  / dx;
                end
            else
                FirstDerivative(i,j) = (ff(ptsnei(i,2))-ff(i)) * 2^lvl  / dx;
          
            end
         
        end

      
end

        if a<b
            FirstDerivative=FirstDerivative'; 
        end

    % 
    % % first order forward
    % 
    % 
    % [nx,ny]=size(f);
    % A=diag(-1*ones(1,nx)) +...
    %     diag(1*ones(1,nx-1),1) +...
    %     diag(0*ones(1,nx-2),2);
    % 
    % % A(end,nx-1:end)=[-1,1];
    % A(end,nx-2:end)=[1,-4,3];
    % 
    % FirstDerivative=1/dx*A*f;
    % 
    % 
    % second order forward
    
        % A=diag(-3*ones(1,nx)) +...
        %     diag(4*ones(1,nx-1),1) +...
        %     diag(-1*ones(1,nx-2),2);
        % A(end,nx-2:end)=[1,-4,3];
        % A(end-1,nx-3:end)=[1,-4,3,0];
        % A=A/2;
        % A=sparse(A);


end
