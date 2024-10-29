function FirstDerivative=ddx_central(ptsnei,ptspos,f,dx)

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

        if ptsnei(i,4) == 0 % left border
          
            secondrightnei = slnn(i,2,2,ptsnei,ptspos,dx);

            FirstDerivative(i,j) = ...
                (-3*ff(i) +4*ff(ptsnei(i,2)) - ff(secondrightnei)) * 2^lvl / dx / 2;

        elseif ptsnei(i,2) == 0
            secondleftnei = slnn(i,4,2,ptsnei,ptspos,dx);

            FirstDerivative(i,j) =  ...
                (3*ff(i) -4*ff(ptsnei(i,4)) + ff(secondleftnei) ) * 2^lvl  / dx / 2;

        else
            FirstDerivative(i,j) = ...
                (ff(ptsnei(i,2))-ff(ptsnei(i,4))) *2^lvl / dx / 2;
        end
    end


end

if a<b
    FirstDerivative=FirstDerivative'; %6000 x 4
end

end
