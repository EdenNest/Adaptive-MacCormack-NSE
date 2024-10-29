function  output = intpRefined(newptspos,oldptspos,vectorfield)
 % recieves the old and new positions and a field to interpolate to (for
 % example u v or T
 
    F = scatteredInterpolant(oldptspos(:,1),oldptspos(:,2),vectorfield);
    output = F(newptspos(:,1),newptspos(:,2));
end


