function output = refinedborder(index,ptsnei,ptspos)
    L = 1e-5;
    H = 8e-6;
    err = 5e-8;
% outputs the sides which this point is a border for the refined zone
% (not for actualy boundaries)

output = find(~ptsnei(index,1:4));

if ptspos(index,2)>H-err; output(output==1)=[];end
if ptspos(index,1)>L-err; output(output==2)=[];end
if ptspos(index,2)<err; output(output==3)=[];end
if ptspos(index,1)<err; output(output==4)=[];end

end












