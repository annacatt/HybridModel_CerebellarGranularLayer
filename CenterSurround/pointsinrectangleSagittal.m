% function [Indices] = pointsinrectangleSagittal2(xydata,height,hk)
function [Indices] = pointsinrectangleSagittal2_new(xydata,PFaxiswidth,SAGaxiswidth,hk)
% Function to select vertices in a rectangle
% with edges SAGaxiswidth, PFaxiswidth
% xydata(2,:) and hk(2) are coodinates along the Sagittal axis

%Indices= find(xydata(1,:)<hk(1)+height & xydata(1,:)>hk(1)-height);
% Indices= find(xydata(1,:)<hk(1)+height & xydata(1,:)>hk(1)-height & xydata(2,:)<hk(2)+0.65 & xydata(2,:)>hk(2)-0.65);
Indices= find(xydata(1,:)<hk(1)+PFaxiswidth/2 & xydata(1,:)>hk(1)-PFaxiswidth/2 & xydata(2,:)<hk(2)+SAGaxiswidth/2 & xydata(2,:)>hk(2)-SAGaxiswidth/2);