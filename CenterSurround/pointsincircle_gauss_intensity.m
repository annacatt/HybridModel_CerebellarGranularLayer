% In base a quello che si vuole ottenere, bisogna giocare con probab=exp(-(dist2(Ind)./radius).^2); es. probab=exp(-100*(dist2(Ind)./radius).^2)
% pointsincircle_gauss_A(coordinates2',1/500*33/2,[0.5 0.5],1)

function [Indices,Intensity] = pointsincircle_gauss_intensity(xydata,radius,hk,K)%,g_sigma)
 dist2 = sqrt((hk(1) - xydata(1,:)).^2   +   (hk(2) - xydata(2,:)).^2); 
    Ind=find(dist2<radius);
    probab=exp(-K*(sqrt(dist2(Ind))./radius));
%     probab = (radius-dist2(Ind))./radius;
%     figure(10)
%     plot(dist2(Ind),probab,'o')
%     bin=binornd(1,probab);
    Indices=Ind;
    Intensity = probab;
%     plot(xydata(1,:),xydata(2,:),'o')
%     hold on
%     plot3(xydata(1,Indices),xydata(2,Indices),Intensity,'ro')
%     xrange(1) = range(xydata(1,Indices))*500;
%     xrange(2) = range(xydata(2,Indices))*500;
end

% 
% function [Indices,count] = pointsincircle_gauss_A(xydata,radius,hk)%,g_sigma)
%  dist2 = sqrt((hk(1) - xydata(1,:)).^2   +   (hk(2) - xydata(2,:)).^2); 
%     Indices=find(dist2<radius*radius);
% %     probab=exp(-(dist2(Ind)./radius^2).^2);
%    probab=exp(-1/4);
%     bin=binornd(1,probab);
%     count=sum(bin);
%     Indices=Ind(bin==1);
% end