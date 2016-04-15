function [Indices,Intensity] = pointsincircle_flattensphere(xydata,radius,centre,divergence)%,g_sigma)
% Calc the prob of selection based on the sphere volume above and below the
% point in the circle
% divergence is the number of nodes that should be contacted by this MF
% this is equal to the divergence of the MF->grc network map which was set
% to 53
dist2 = (centre(1) - xydata(1,:)).^2   +   (centre(2) - xydata(2,:)).^2; 
Ind=find(dist2<radius*radius); % indexes of nodes in the circle
% Calc the height of the sphere dome above the point
% dist2 is the dist from centre
n_dist = sqrt(dist2(Ind))./radius; % dist normalized by max dist, i.e. the radius of the circle
% we are looking for the sin(x) where x is the azimut angle
% dist is equal to cos(x)
% the azimut angle is arccos(d)
dome_height = sin(acos(n_dist)); % this should be <= 1
% max(dome_height)
probab=dome_height;
% The prob of selecting on node is now equl to the the nuber of grcs that
% should be above and below it in the sphere.
% Since each node represents n_grcs the probabilitiy of selecting it should
% be moltiplied by the n_grcs and divided by the total number of grcs in
% the sphere. Since these rescaling affctes all nodes in the circle it
% becomes irrelevant for the calculation of the prob.
% Then the effect of this MF to each node
% should be divided by n_grcs. Since the strength of the Mf->node is now
% arbitrary, this rescaling is meaningless.
% The probability must be normalized by the sum over all nodes
probab = probab ./ sum(probab);
% % Setup for multiple trials 
% % Generate n=divergence trials
% trials = rand([divergence,1]);
% [P,X] = meshgrid(cumsum(probab),trials); % probab is in P colns 
% ends = [];
% for i = 1:size(P,1)
%     p = P(i,:);
%     idx = find(p<trials(i));
%     if ~isempty(idx)
%         ends = [ends idx(end)];
%     end
% end

Indices = Ind;% (ends);
% The total signal intensity is equal to the number of grc activated by
% this MF
Intensity = probab * divergence;
% plot3(xydata(1,Ind),xydata(2,Ind),probab,'o')
end

