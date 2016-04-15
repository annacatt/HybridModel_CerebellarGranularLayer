function [Indices] = pointsincircle(xydata,radius,hk)
dist2 = (hk(1) - xydata(1,:)).^2   +   (hk(2) - xydata(2,:)).^2; 
Indices= find(dist2<radius*radius & dist2>0);
end