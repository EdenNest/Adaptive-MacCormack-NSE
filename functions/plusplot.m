function plusplot(m,ptsnei,ptspos)
    indices = ptsnei(m,1:4);
    indices(indices==0) = m;
    
    A = ptspos(indices,:);
    x = A(:,1); y=A(:,2);
    line([x(1) x(2);x(3) x(4)],[y(1) y(2);y(3) y(4)],'Color', 'k');%,'lineWidth',2.5);

end
