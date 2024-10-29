
function plotRefinement(fignum , current_iteration ,Number_of_refined,newptsnei, newptspos)
figure(fignum)
 title(['Refinement of interation ' num2str(current_iteration)])
 k=length(newptsnei);

                      for i = 1:Number_of_refined
                     
                          plusplot(k-Number_of_refined+i,newptsnei,newptspos);
                        
                      end
                      drawnow
                      axis equal tight;hold on;
end

% 
% function plotRefinement(fignum , n ,newptsnei, newptspos)
% figure(fignum)
%  title(['Refinement of interation ' num2str(n)])
% 
%                       for i = 1:length(newptspos)
% 
%                           plusplot(i,newptsnei,newptspos);
% 
%                       end
% 
%                       axis equal tight;
% end
