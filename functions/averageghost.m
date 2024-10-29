function output = averageghost(node,field,direction,ptspos)
dy = 8e-6/79;
dx=1e-5/74;
pos = ptspos(node,:);
switch direction
    case 1
        pos1 = pos+[+dx,dy]/2;
        pos2 = pos+[-dx,dy]/2;
    case 2
        pos1 = pos+[dx,dy]/2;
        pos2 = pos+[dx,-dy]/2;
    case 3
        pos1 = pos+[dx,-dy]/2;
        pos2 = pos+[-dx,-dy]/2;
    case 4
        pos1 = pos+[-dx,dy]/2;
        pos2 = pos+[-dx,-dy]/2;


end
% direction
% pos1
% pos2
% node
% FindIt(ptspos,pos1(1),pos1(2))
% FindIt(ptspos,pos2(1),pos2(2))
value1 = field(FindIt(ptspos,pos1(1),pos1(2)));
value2 = field(FindIt(ptspos,pos2(1),pos2(2)));

output = 1/2*(value1+value2);
end
