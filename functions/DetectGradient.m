function indices = DetectGradient(ptsnei,ptspos,v_n,u_n,T_n,percent,dx,dy)
   
  
gg = (1-percent/100);
 
 

indices = do(ptsnei,ptspos,v_n,u_n,T_n,percent,dx,dy);

if isempty(indices)
    percent = percent + 5;
    indices = do(ptsnei,ptspos,v_n,u_n,T_n,percent,dx,dy);
    disp(['increased percent to ' num2str(percent)])
end



    function indices = do(ptsnei,ptspos,v_n,u_n,T_n,percent,dx,dy)

 indices = [];

     all = [ddy_central(ptsnei,ptspos,T_n,dy) ddy_central(ptsnei,ptspos,v_n,dy) ddy_central(ptsnei,ptspos,T_n,dy) ...
        ddx_central(ptsnei,ptspos,u_n,dx) ddx_central(ptsnei,ptspos,v_n,dx) ddx_central(ptsnei,ptspos,T_n,dx)];

    for i = 1:size(all,2)
        phi = all(:,i);
        % threshold = gg * max( phi(ptsnei(:,6)==0 ));
        threshold = gg * max( phi );
        indices = cat(1,indices,find(phi>threshold));
    end
   
    remove = find(ptsnei(:,5)==1);  %  last level of refined
    indices(ismember(indices, remove)) = [];

    % remove2 = find(ptspos(:,2)<3*dy);    % no botton border
    % indices(ismember(indices, remove2)) = [];

    indices = unique(indices);

    end
    

end
% 
% Error Estimation
% Error estimation is crucial for identifying cells that require 
% refinement. There are several techniques for error estimation:
% 
% Residual-based Error Estimation: This method calculates the residual
% of the solution within each cell. Cells with high residuals are 
% candidates for refinement.

% Gradient-based Error Estimation: 
% This method involves computing the
% gradient (or higher derivatives) of the solution. Cells where the 
% gradient is large indicate regions with high solution variation and are
% marked for refinement.

% Solution-based Error Estimation: 
% Directly estimate the error in the 
% solution based on some predefined criteria, such as the difference
% between solutions on different levels of refinement.
% Adjoint-based Error Estimation: Use adjoint equations to estimate the 
% error in a specific functional of interest (e.g., drag in fluid dynamics).


% 2. Refinement Criteria
% Based on the error estimation, cells are selected for refinement. 

% Common criteria include:
% 
% Thresholding: Cells with an error estimate above a certain threshold are
% refined. The threshold can be static or dynamically adjusted based on the
% overall error distribution.

% Relative Error: Refine cells where the error relative to the average or 
% maximum error in the domain is significant.

% Percentage of Cells: Refine the top percentage of cells with the highest
% error estimates
