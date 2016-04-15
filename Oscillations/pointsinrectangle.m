function [Indices] = pointsinrectangle(xydata,SAGaxiswidth,hk)
% Function to select vertices along a stripe aligned along the PF axis
% and wide SAGaxiswidth
% xydata(2,:) and hk(2) are coodinates along the Sagittal axis
Indices= find(xydata(2,:)<hk(2)+SAGaxiswidth/2 & xydata(2,:)>hk(2)-SAGaxiswidth/2);
end

% function [Indices] = pointsinrectangle(xydata,height,hk)
% Indices= find(xydata(1,:)<hk(1)+height & xydata(1,:)>hk(1)-height);
% end