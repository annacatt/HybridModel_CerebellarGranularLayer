close all
clear all
%
SaveFigures=1;
show_max = 1;
show_movie = 1;

kdata(1).val=['23232']; % # of GrC nodes
for i=1:length(kdata)
    file_name = ['load WholeNetwork_dynamics_',num2str(kdata(i).val),'.mat']; % mat file to load
    eval(file_name);
end

% Output file name
ofn = 'DynamicsWholeNetwork_Oscillations'; % It would be best if this parameter was embedded at simulation time in the data file!!!!!

nn1=eval(['results',kdata(1).val,'.nn1']);
coordinates1=eval(['results',kdata(1).val,'.coordinates1']);
triangles1=eval(['results',kdata(1).val,'.triangles1']);

nn2=eval(['results',kdata(1).val,'.nn2']);
coordinates2=eval(['results',kdata(1).val,'.coordinates2']);
x__=eval(['results',kdata(1).val,'.x__']);
y__=eval(['results',kdata(1).val,'.y__']);

W=eval(['results',kdata(1).val,'.W']);
eval(['Time = results',kdata(1).val,'.t;'])

x__=x__*500;
y__=y__*500;

coords1 = coordinates1 * 500;
coords2 = coordinates2 * 500;
[X,Y] = meshgrid(x__,y__);

if show_max
    [ma,mx] = max(max(W(:,3*nn1+1:3*nn1+nn2),[],2)); % mx is the max frame
    th = 0.7 * ma; % 70% of max in max frame
    max_frame = W(mx,3*nn1+1:3*nn1+nn2);
    above_th = find(max_frame>th);
end


vidObj = VideoWriter([ofn '.avi']); % Prepare the new file
n_frames = length(Time); % Get the number of frames
simed_time = Time(end)/10; % Get the total simulated time; time in ms divide by 10 for numerical reasons
slow_reate = 10; % Slow down the visualization rate: the visualization is 10 times slower than real time
frame_rate = n_frames/simed_time*1000/slow_reate; % Calc the frame rate multiply by 1000 ad simed_time is in ms

vidObj.FrameRate = frame_rate; % Set the final duration of the muvie to 10 sec so it is 10 times slower than real time

open(vidObj);

if show_movie
    figure()
    set(gcf,'Units','normalized','OuterPosition',[0,0,1,1],'Color','w')

    for i=1:2:length(Time)
        g = reshape(W(i,3*nn1+1:3*nn1+nn2),length(y__),length(x__));
        surf(X,Y,g,'edgecolor','k','facecolor','interp','facealpha',0.7);
        hold on
        stem3(coords1(1,:),coords1(2,:),W(i,1:nn1),'fill','Marker','none','LineWidth',1.5)
        zlim([-0.1 3])
        caxis([0 1])
        t=title(['t=',num2str(Time(i)/10)]);
        set(gca,'FontSize',60)
        set(t,'FontSize',60)
        axis tight
        set(gca,'nextplot','replacechildren');

        hold off
        if SaveFigures == 1
            if any(Time(i)/10==[6 8 10 12 13 16 20 22 24])
                if any(Time(i)/10==[8 10 13 16 20 24])
                    set(gca,'ZTickLabel',[])
                end
                export_fig(['CenterSurroundWhole_SS',num2str(Time(i)/10) '_ms.pdf'])
            end
        end
        axis([0 500 0 1500 -0.5 3]);
        set(gca,'DataAspectRatio',[1 1 0.01],'ZLim',[-0.5 3])
        view(-35.5000, 44);

        % Write each frame to the file.
        currFrame = getframe(gcf);
        writeVideo(vidObj,currFrame);

        pause(0.1)
    end
end

disp(['Movie duration is ' mat2str(vidObj.Duration) ' sec for a simulated time of ' mat2str(simed_time) ' ms']) % Set the fina duration of the muvie to 10 sec so it is 10 times slower than real time

close(vidObj); % Close the file