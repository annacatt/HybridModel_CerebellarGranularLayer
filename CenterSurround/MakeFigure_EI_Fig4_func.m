function [E,I,D,X,Y]=MakeFigure_EI_Fig4_func(grc_grid_res)


nn2 = grc_grid_res^2;

SaveFigures=1;
show_max = 1;
show_movie = 1;
show_traces = 1;

results_fn = ['data2DMultispecies' num2str(nn2) '_grc_grid_res_' mat2str(grc_grid_res) '.mat'];
WG = load(results_fn);
results_fn = ['data2DMultispecies_Blockinh' num2str(nn2) '_grc_grid_res_' mat2str(grc_grid_res) '.mat'];
WOG = load(results_fn);


nn1=eval(['WG.results',num2str(nn2),'.nn1']);
coordinates1=eval(['WG.results',num2str(nn2),'.coordinates1']);
triangles1=eval(['WG.results',num2str(nn2),'.triangles1']);

nn2=eval(['WG.results',num2str(nn2),'.nn2']);
coordinates2=eval(['WG.results',num2str(nn2),'.coordinates2']);
x__=eval(['WG.results',num2str(nn2),'.x__']);
y__=eval(['WG.results',num2str(nn2),'.y__']);

% Find nodes near the spot
% stim_centre = eval(['WG.results',num2str(nn2),'.stim_centre']);
stim_centre =[ones(10,1)*.5 ones(10,1)*.5]  + rand(10,2).*0.09;
x__idx = find(x__>min(stim_centre(:,1))-0.1 & x__<max(stim_centre(:,1))+0.05);
y__idx = find(y__>min(stim_centre(:,2))-0.1 & y__<max(stim_centre(:,2))+0.05);
[X_idx,Y_idx] = meshgrid(x__idx,y__idx);
% W=eval(['results',num2str(nn2),'.W']);
% W=W(:,3*nn1+1:3*nn1+nn2);
% max(max(W))
eval(['Time = WG.results',num2str(nn2),'.t;'])
eval(['tTGr_in = WG.results',num2str(nn2),'.tTGr_in;'])

% x__=x__*500;
% y__=y__*500;

coords1 = coordinates1;
coords2 = coordinates2;
[X,Y] = meshgrid(x__,y__);
X = X*1e3;
Y = Y*1e3;
xmin = min(min(X(X_idx,Y_idx)));
xmax = max(max(X(X_idx,Y_idx)));
ymin = min(min(Y(X_idx,Y_idx)));
ymax = max(max(Y(X_idx,Y_idx)));

W_grc_idxs = 3*nn1+1:3*nn1+nn2;
W_goc_idxs = 1:nn1;

eval(['WG_W_results = WG.results',num2str(nn2),'.W;']);
eval(['WOG_W_results = WOG.results',num2str(nn2),'.W;']);
eval(['WG_t_results = WG.results',num2str(nn2),'.t;']);
eval(['WOG_t_results = WOG.results',num2str(nn2),'.t;']);


if show_max
%     figure()
%     set(gcf,'Units','normalized','OuterPosition',[0,0,1,1])
    [ma,mx] = max(max(WG_W_results(:,W_grc_idxs),[],2)) % mx is the max frame
    th = 0.7 * ma; % 70% of max in max frame
    max_frame = WG_W_results(mx,3*nn1+1:3*nn1+nn2);
    above_th = find(max_frame>th);
    stim_width(1) = range(coords2(1,above_th));
    stim_width(2) = range(coords2(2,above_th));
    disp(['Spot diameter is ' sprintf('%3.0f',stim_width*1000) ' um'])
%     pause
%     trisurf(triangles2,coords2(1,:),coords2(2,:),W(mx,3*nn1+1:3*nn1+nn2),'edgecolor','k','facecolor','interp','facealpha',0.5);
%     view(2)
end


if show_traces
    figure(1)
    subplot(2,2,1)
    plot(Time/10,WG_W_results(:,W_goc_idxs))
    ylabel('GoC')
    subplot(2,2,3)
    plot(Time/10,WG_W_results(:,W_grc_idxs),'.')
    ylabel('GrC')
    subplot(2,2,2)
    plot(Time/10,WOG_W_results(:,W_goc_idxs))
    ylabel('GoC')
    subplot(2,2,4)
    plot(Time/10,WOG_W_results(:,W_grc_idxs))
    ylabel('GrC Inh blocked')
    
end


if show_movie
    fontsize = 40;
    linewidth = 2
    figure(2)
    set(gcf,'Units','normalized','OuterPosition',[0,0,1,1],'Color','w')
%     subplot(2,2,1)
    [ma,mx] = max(max(WG_W_results(:,3*nn1+1:3*nn1+nn2),[],2)); % mx is the max frame

    i = mx; % get the point max activation in the withGOC sim
    max_act = mat2str((WG_t_results(i)/10-tTGr_in))
    disp(['Max activation of GLN withGOC at t=' max_act ' ms'])
    E = WG_W_results(i,3*nn1+1:3*nn1+nn2);
%     plot3(coordinates2(1,:),coordinates2(2,:),E,'.k')
    
    E = reshape(E,length(y__),length(x__));%length(coords2(1,:)),length(coords2(2,:)));
%     subplot(2,2,2)
%     surf(X,Y,E,'edgecolor','k','facecolor','interp','facealpha',0.5);
    S_E = surf(X(X_idx,Y_idx),Y(X_idx,Y_idx),E(X_idx,Y_idx),'edgecolor','k','facecolor','interp');
    shading interp
% end
%     hold on
%     stem3(WG.results1600.coordinates1(1,:)*500,WG.results1600.coordinates1(2,:)*500,WG.results1600.W(i,1:nn1),'fill','Marker','none','LineWidth',1.5)
%     zlim([-0.1 3])
    xlim([xmin xmax])
    ylim([ymin ymax])
    caxis([0 1])
    set(gca,'FontSize',fontsize,'LineWidth',linewidth)
    t=title(['E peak t=', max_act]);
    xlabel(['Transverse axis ( ' '\mu' 'm)'])
    zlabel(['Normalized Memb. Pot.'])
    ylabel(['Sagittal axis (' '\mu' 'm)'])
    set(t,'FontSize',fontsize)
    hold off
    fig_name = ['CS_Fig4A_' mat2str(grc_grid_res) '.png']
    if SaveFigures
        export_fig(fig_name)
        close(2)
    end
    figure(3)
    set(gcf,'Units','normalized','OuterPosition',[0,0,1,1],'Color','w')
    subplot(2,1,1)
    E2_shift = 4;
    i = find(WG_t_results/10 >= WG_t_results (mx)/10+E2_shift,1,'first');
    disp(['E2 peak of GLN withGOC at t=' mat2str((WG_t_results(i)/10-tTGr_in)) ' ms'])
    g = reshape(WG_W_results(i,3*nn1+1:3*nn1+nn2),length(y__),length(x__));%length(coords2(1,:)),length(coords2(2,:)));
    
    surf(X(X_idx,Y_idx),Y(X_idx,Y_idx),g(X_idx,Y_idx),'edgecolor','k','facecolor','interp');
%     surf(X,Y,g,'edgecolor','k','facecolor','interp','facealpha',0.5);
    %     trisurf(WOG.results1600.triangles2,WG.results1600.coordinates2(1,:)*500,WOG.results1600.coordinates2(2,:)*500,WG.results1600.W(i,3*nn1+1:3*nn1+nn2),'edgecolor','k','facecolor','interp','facealpha',0.5);
%     hold on
%     trisurf(WOG.results1600.triangles2,WG.results1600.coordinates2(1,:)*500,WOG.results1600.coordinates2(2,:)*500,WOG.results1600.W(i,3*nn1+1:3*nn1+nn2),'edgecolor','k','facecolor','interp','facealpha',0.5);
    shading interp
%     stem3(WOG.results1600.coordinates1(1,:)*500,WOG.results1600.coordinates1(2,:)*500,WOG.results1600.W(i,1:nn1),'fill','Marker','none','LineWidth',1.5)
    zlim([-0.1 1])
    xlim([xmin xmax])
    ylim([ymin ymax])
    caxis([0 1])
    view(-37.5,80)
    set(gca,'FontSize',fontsize,'LineWidth',linewidth,'ZTick',[0 1])
    t=title(['E_2 peak t=',num2str((WG_t_results(i)/10-tTGr_in))]);
%     xlabel(['Transverse axis (' ' \mu' 'm)'])
%     ylabel(['Sagittal axis (' ' \mu' 'm)'])
    set(t,'FontSize',fontsize)
    hold off

    subplot(2,1,2)
    g = reshape(WOG_W_results(i,3*nn1+1:3*nn1+nn2),length(y__),length(x__));%length(coords2(1,:)),length(coords2(2,:)));
    surf(X(X_idx,Y_idx),Y(X_idx,Y_idx),g(X_idx,Y_idx),'edgecolor','k','facecolor','interp'); %,'facealpha',0.02);
%     hold on
%     surf(X(X_idx,Y_idx),Y(X_idx,Y_idx),E(X_idx,Y_idx),'edgecolor','k','facecolor','interp');
%     surf(X,Y,g,'edgecolor','k','facecolor','interp','facealpha',0.5);
%     trisurf(WOG.results1600.triangles2,WG.results1600.coordinates2(1,:)*500,WOG.results1600.coordinates2(2,:)*500,WOG.results1600.W(i,3*nn1+1:3*nn1+nn2),'edgecolor','k','facecolor','interp','facealpha',0.5);
    shading interp
%     stem3(WOG.results1600.coordinates1(1,:)*500,WOG.results1600.coordinates1(2,:)*500,WOG.results1600.W(i,1:nn1),'fill','Marker','none','LineWidth',1.5)
    zlim([-0.1 1.5])
    xlim([xmin xmax])
    ylim([ymin ymax])
    caxis([0 1])
    view(-37.5,80)
    set(gca,'FontSize',fontsize,'LineWidth',linewidth,'ZTick',[0 1.2])
    t=title(['E_{2ib}: Inhibition blocked']);
    xlabel(['Transverse axis (' ' \mu' 'm)'])
    ylabel(['Sagittal axis (' '\mu' 'm)'])
    set(t,'FontSize',fontsize)
    hold off
    fig_name = ['CS_Fig4B_' mat2str(grc_grid_res) '.png']
    if SaveFigures
        export_fig(fig_name)
        close(3)
    end

    
    % calc WG-WOG at time of E2 peak to get the net inhibition I
    figure(4)
    set(gcf,'Units','normalized','OuterPosition',[0,0,1,1],'Color','w')
%     subplot(2,2,3)
    I = reshape(WOG_W_results(i,3*nn1+1:3*nn1+nn2)-WG_W_results(i,3*nn1+1:3*nn1+nn2),length(y__),length(x__));
    surf(X(X_idx,Y_idx),Y(X_idx,Y_idx),I(X_idx,Y_idx)-I(1,1),'edgecolor','k','facecolor','interp');
%     hold on
%     surf(X(X_idx,Y_idx),Y(X_idx,Y_idx),E(X_idx,Y_idx),'edgecolor','k','facecolor','interp');
%     surf(X,Y,I,'edgecolor','k','facecolor','interp','facealpha',0.5);
%     trisurf(WG.results1600.triangles2,WG.results1600.coordinates2(1,:)*500,WG.results1600.coordinates2(2,:)*500,WOG.results1600.W(i,3*nn1+1:3*nn1+nn2)-WG.results1600.W(i,3*nn1+1:3*nn1+nn2),'edgecolor','k','facecolor','interp','facealpha',0.5);
    %     shading interp
%     hold on
%     stem3(WG.results1600.coordinates1(1,:)*500,WG.results1600.coordinates1(2,:)*500,WG.results1600.W(i,1:nn1)-WOG.results1600.W(i,1:nn1),'fill','Marker','none','LineWidth',1.5)
    zlim([-0.1 1])
    xlim([xmin xmax])
    ylim([ymin ymax])
    caxis([0 1])
    set(gca,'FontSize',fontsize,'LineWidth',linewidth)
    t=title(['-I=E_{2ib}-E_2']);
    xlabel(['Transverse axis (' ' \mu' 'm)'])
    zlabel(['Normalized Memb. Pot.'])
    ylabel(['Sagittal axis (' '\mu' 'm)'])
    set(t,'FontSize',fontsize)
%     set(gca,'FontSize',fontsize,'LineWidth',linewidth)
    hold off
    fig_name = ['CS_Fig4C_' mat2str(grc_grid_res) '.png']
    if SaveFigures
        export_fig(fig_name)
        close(4)
    end

    % calc E-I to get the center/surround or mexican hat sa in Solinas 2010
    figure(5)
    set(gcf,'Units','normalized','OuterPosition',[0,0,1,1],'Color','w')
%     subplot(2,2,4)
    D = (E-I);
    D_S = surf(X(X_idx,Y_idx),Y(X_idx,Y_idx),D(X_idx,Y_idx) - D(1,1),'edgecolor','k','facecolor','interp');
    
    surf(X,Y,(E-I),'edgecolor','k','facecolor','interp');%,'facealpha',0.5); %ho aggiunto questa riga
%     colorbar
    caxis([-0.2 0.4])

%     trisurf(WG.results1600.triangles2,WG.results1600.coordinates2(1,:)*500,WG.results1600.coordinates2(2,:)*500,E-I,'edgecolor','k','facecolor','interp','facealpha',0.5);

%     hold on
%     stem3(WG.results1600.coordinates1(1,:)*500,WG.results1600.coordinates1(2,:)*500,WG.results1600.W(i,1:nn1)-WOG.results1600.W(i,1:nn1),'fill','Marker','none','LineWidth',1.5)
    zlim([-1 1])
    xlim([xmin xmax])
    ylim([ymin ymax])
    view(-47.50,36)
%     hold on
%     plot3(X(:,1),Y(:,1),D(:,33),'k','LineWidth',4)
%     plot3(X(1,:),Y(1,:),D(33,:),'r','LineWidth',4)
%     
    %     view(-125,84)
%         shading interp
    set(gca,'FontSize',fontsize,'LineWidth',linewidth);
%     set(gca,'DataAspectRatio',[1 1 0.01] ,'ZLim',[-0.5 1])
    t=title(['E-I']);
    xlabel(['Transverse axis (' ' \mu' 'm)'])
    zlabel(['Normalized Memb. Pot.'])
    ylabel(['Sagittal axis (' '\mu' 'm)'])
    set(t,'FontSize',fontsize)
    hold off
    fig_name = ['CS_Fig4D_' mat2str(grc_grid_res) '.png']
    if SaveFigures
        export_fig(fig_name)
        close(5)
    end

%     figure(6)
%     D_SS = D(33,:);
%     D_TT = D(:,33);
%     
%     plot(0:64,D_SS,'k')
%     plot(0:64,D_TT,'r')
    
end
