function MultiscaleCS_large_A_func(grc_grid_res,inh_block)

global nn1 nn2 d_Cont d_Disc InElecDisc InElecCont T gsyn_fromCont gsyn_fromDisc InChimFromCont
global InChimFromDisc Wsyn IGr_Mossy IGo_Mossy InElecCont_sparse InElecCont_n 
global InElecDisc_sparse InElecDisc_n InChimFromCont_nfrac InChimFromDisc_nfrac InChimFromCont_sparse 
global InChimFromDisc_sparse tT_end tTGr_in tTGr_in_end
global tTGo_in tTGo_in_end

rebuild = 1;
tic
if rebuild == 1
    SaveData=1;
    
    d_Cont=0;%0.005;%5;%0.005;
    d_Disc=0;%0.005;%0.05;
    gsyn_fromCont=0.8;
    if inh_block == 1
        gsyn_fromDisc=0.1; % It was 0.05 but in the old CS I had 26 Goc concentrated in a small area so the grcs were much more inhibited
    else
        gsyn_fromDisc= inh_block;
    end
    % load grid1_2
    
    
    load grid1rect_CS_200x200_40Goc.mat;
    [InElecDisc_sparse,InElecDisc_n] = Conv_struc2matix(InElecDisc);
    
    % %griglia regolare GrCs:
    nn2_rad= grc_grid_res; %65; %numero Gr per lato
    % Total of 1600 nodes each for 10 grcs
    lung_lato = 0.2; % 1 == 500 um -> 0.2 = 100um
    lung_lato2 =0.2;
    nn2=nn2_rad^2;
    h=lung_lato/(nn2_rad-1);
    x__=linspace(0.4,0.4+lung_lato,nn2_rad);
    y__=linspace(0.4,0.4+lung_lato2,nn2_rad);
    [x,y]=meshgrid(x__,y__);
    siz=size(x,2);
    x_=reshape(x,siz^2,1);
    y_=reshape(y,siz^2,1);
    coordinates2=[x_ y_]; %coordinate Gr
%     
%     plot(coordinates2(:,1),coordinates2(:,2),'ok')
%     hold on
%     plot(coordinates1(:,1),coordinates1(:,2),'or')
%     pause()
    %
    for k=1:nn2
        InElecCont(k).val=pointsincircle(coordinates2',h+0.001,coordinates2(k,:)); %input dal continuo sul discreto
    end
    plot(coordinates2(:,1),coordinates2(:,2),'ok')
    [InElecCont_sparse,InElecCont_n] = Conv_struc2matix(InElecCont);
    
    % Input configuration
    
    % GrCs spot activation
    IGr_Mossy=zeros(nn2,1);
%     mf_grc_div = 53;
    mf_grc_radius = 0.090; % max grc dend length 30 um
%     num_of_active_mfs = 10;
    stim_centre =[ones(10,1)*.5 ones(10,1)*.5]  + rand(10,2).*0.09;
    mf_grc_input_strength = 0.4;
    spot_decay = 0.4;
    % Each active mf contacts mf_grc_div nodes, it can contact more than
    % once the same node as 1 node represents more than one grc
    % We assume that a certain number of mf (num_of_active_mfs) is
    % transmitting the stimulus and then pull together their targets.
    % Get the list of nodes contacted by the 10 mfs.
    [nn2_center,Intensity_mf_grc] = pointsincircle_gauss_intensity(coordinates2',mf_grc_radius,[0.5 0.5],spot_decay);
    Intensity_mf_grc = (Intensity_mf_grc - min(Intensity_mf_grc));
    Intensity_mf_grc = Intensity_mf_grc ./ max(Intensity_mf_grc);
    Intensity_mf_grc = Intensity_mf_grc + 0.05 * rand(size(Intensity_mf_grc));
%     nn2_center = [];
%     for mf_i = 1:num_of_active_mfs
%         nn2_center = [nn2_center pointsincircle_flattensphere(coordinates2',mf_radius,stim_centre(mf_i,:),mf_grc_div)]; % Input dall MF sul continuo
%     end
    % Find repetitions and count them
%     [nn2_center_a,nn2_center_b]=hist(nn2_center,unique(nn2_center)); % a(i) is the number of repetition of b(i)
%     nn2_center_b = round(nn2_center_b);
figure(1)
subplot(2,1,1)
    plot(coordinates2(:,1),coordinates2(:,2),'.k')
    hold on
% %     plot(coordinates2(nn2_center_b,1),coordinates2(nn2_center_b,2),'ro')
% %     % Normalizec the mf stim to 13 grcs in one node
% %     nn2_center_a = nn2_center_a ./ max(nn2_center_a) .* 2;
% %     IGr_Mossy(nn2_center_b)=nn2_center_a * mf_grc_input_strength;
    plot3(coordinates2(nn2_center,1),coordinates2(nn2_center,2),Intensity_mf_grc,'ro')
    view(-20,30)
%     title('Mf->GrC input')
    IGr_Mossy(nn2_center)= Intensity_mf_grc .*  mf_grc_input_strength;

    tTGr_in = [10*ones(nn2,1) 100e3*ones(nn2,1)]; %la correnta inizia dopo
    tTGr_in_end = [15*ones(nn2,1) 110e3*ones(nn2,1)]; %la correnta finisce dopo
%     tT_in = [5*ones(nn1,1) 100e3*ones(nn1,1)]; %la correnta inizia dopo
%     tT_end = [15*ones(nn1,1) 110e3*ones(nn1,1)]; %la correnta finisce dopo
    
    % GoCs spot activation
    IGo_Mossy=zeros(nn1,1);
    mf_goc_input_strength = 0.1;
    mf_goc_radius = mf_grc_radius+0.1; % This has to consider the radius of the GoC basolateral dendrites. It was set to 100 um in Solinas2010
    [nn1_center,Intensity_mf_goc] = pointsincircle_gauss_intensity(coordinates1',mf_goc_radius,[0.5 0.5],spot_decay);
    figure(1)
    subplot(2,1,2)
    plot(coordinates1(:,1),coordinates1(:,2),'.k')
    hold on
    plot3(coordinates1(nn1_center,1),coordinates1(nn1_center,2),Intensity_mf_goc,'ro')
    view(-20,30)
    title('Mf->GoC input')
    IGo_Mossy(nn1_center)= mf_goc_input_strength; %Intensity_mf_goc .* mf_goc_input_strength;
    tTGo_in = 10*ones(nn1,1); % for now I put MF->GoC to 0
    tTGo_in_end = 15*ones(nn1,1); % for now I put MF->GoC to 0
    
    v0=spalloc(3*nn1,1,nn1);
    omega0=spalloc(3*nn2,1,10); %corrispettivo di v0 nel continuo
    %
    T=0.5*[ones(1,nn1)';ones(1,nn2)'];
    W0=[v0;omega0];
    %
    vIsyn=-0.2;%reversal potential Golgi cells
    omegaEsyn=0.93;%reversal potential granular cells
    vIsyn=vIsyn*ones(nn2,1);
    omegaEsyn=omegaEsyn*ones(nn1,1);
    Wsyn=[vIsyn;omegaEsyn];
    %
    %Per calcolare i due input
    %input chimico del discreto sul continuo
    heig=2/10; % 1/10 on either side along the Sag axis
    for k=1:nn1
        InChimFromCont(k).val=pointsinrectangle(coordinates2',heig,coordinates1(k,:)); %input dal continuo sul discreto
    end
    for k=1:nn1
        for i=1:length(InChimFromCont(k).val)
            if InChimFromCont(k).val(i)==k
                InChimFromCont(k).val(i)=[];
                break
            end
        end
    end
    % create sparse matrix for InChimFromCont
    [InChimFromCont_sparse,InChimFromCont_nfrac] = Conv_struc2matix(InChimFromCont);
    InChimFromCont_nfrac = 4./InChimFromCont_nfrac;
    
    % fatt=nn2/23547;
    % heigSag=1/2.8/2;
    PFaxiswidth=0.18; % width of rectancle along the PF axis as 1 corresponds to 500 um 1/2.8=180um
    PFaxiswidth = PFaxiswidth * 2; % double it so that from the ref of 0.5 mm it goes to the 1 short edge lenth of the rectangle
    SAGaxiswidth= 0.65; % width of rectangle along the sagittal axis
    SAGaxiswidth= SAGaxiswidth * 2; % double it so that from the ref of 0.5 mm it goes to the 1 short edge lenth of the rectangle
    
    for k=1:nn2
        %     InChimFromDisc(k).val=pointsinrectangleSagittal2(coordinates1',heigSag,coordinates2(k,:)); %input dal discreto sul continuo
        InChimFromDisc(k).val=pointsinrectangleSagittal2_new(coordinates1',PFaxiswidth,SAGaxiswidth,coordinates2(k,:)); %input dal discreto sul continuo
        if length(InChimFromDisc(k).val)>40%/fatt
            X=randperm(length(InChimFromDisc(k).val),40); %round(40/fatt)
            InChimFromDisc(k).val=InChimFromDisc(k).val(X);
        end
    end
    [InChimFromDisc_sparse, InChimFromDisc_nfrac] = Conv_struc2matix(InChimFromDisc);
    InChimFromDisc_nfrac = 4./InChimFromDisc_nfrac;

    save CS_2.mat;
else
    load CS_2.mat;
    
    d_Cont=0.005;%5;%0.005;
    d_Disc=0.005;%0.05;
    gsyn_fromCont=0.8;
    gsyn_fromDisc=0.05; % It was 0.05 but in the old CS I had 26 Goc concentrated in a small area so the grcs were much more inhibited

    IGr_Mossy=zeros(nn2,1);
    mf_grc_div = 53;
    mf_grc_radius = 0.080; % max grc dend length 30 um
%     num_of_active_mfs = 10;
    stim_centre =[ones(10,1)*.5 ones(10,1)*.5]  + rand(10,2).*0;%0.09;
    mf_grc_input_strength = 0.5; %0.15;
    % Each active mf contacts mf_grc_div nodes, it can contact more than
    % once the same node as 1 node represents more than one grc
    % We assume that a certain number of mf (num_of_active_mfs) is
    % transmitting the stimulus and then pull together their targets.
    % Get the list of nodes contacted by the 10 mfs.
    [nn2_center,Intensity_mf_grc] = pointsincircle_gauss_intensity(coordinates2',mf_grc_radius,[0.5 0.5],80);
    plot3(coordinates2(nn2_center,1),coordinates2(nn2_center,2),Intensity_mf_grc,'ro')
    above_th = find(Intensity_mf_grc>0.1);
    stim_width(1) = range(coordinates2(nn2_center(above_th),1));
    stim_width(2) = range(coordinates2(nn2_center(above_th),2));
    stim_width

%     nn2_center = [];
%     for mf_i = 1:num_of_active_mfs
%         nn2_center = [nn2_center pointsincircle_flattensphere(coordinates2',mf_radius,stim_centre(mf_i,:),mf_grc_div)]; % Input dall MF sul continuo
%     end
    % Find repetitions and count them
%     [nn2_center_a,nn2_center_b]=hist(nn2_center,unique(nn2_center)); % a(i) is the number of repetition of b(i)
%     nn2_center_b = round(nn2_center_b);
% figure(1)
% subplot(2,1,1)
%     plot(coordinates2(:,1),coordinates2(:,2),'.k')
%     hold on
% %     plot(coordinates2(nn2_center_b,1),coordinates2(nn2_center_b,2),'ro')
% %     % Normalizec the mf stim to 13 grcs in one node
% %     nn2_center_a = nn2_center_a ./ max(nn2_center_a) .* 2;
% %     IGr_Mossy(nn2_center_b)=nn2_center_a * mf_grc_input_strength;
%     plot3(coordinates2(nn2_center,1),coordinates2(nn2_center,2),Intensity_mf_grc,'ro')
%     view(-20,30)
%     title('Mf->GrC input')
    IGr_Mossy(nn2_center)= Intensity_mf_grc .*  mf_grc_input_strength;
 
end
%
coordinates1=coordinates1';
coordinates2=coordinates2';

% print_sizes()

setup_time=toc
tic
tspan=[0:0.1:50]; %tspan=[0:5:1000];
% f_multiscaleNewDEFMossyGolgi(0,W0);
[t,W]=ode45('f_multiscaleNewDEFMossyGolgi',tspan,W0);
tempo=toc

W(:,1:nn1)=3*W(:,1:nn1);

t=round(t*(10^13))/(10^12);
if SaveData==1
    eval(['results',num2str(nn2),'.coordinates1=coordinates1;']);
    eval(['results',num2str(nn2),'.triangles1=triangles1;']);
    eval(['results',num2str(nn2),'.nn1=nn1;']);
    eval(['results',num2str(nn2),'.coordinates2=coordinates2;']);
%     eval(['results',num2str(nn2),'.triangles2=triangles2;']);
    eval(['results',num2str(nn2),'.nn2=nn2;']);
    eval(['results',num2str(nn2),'.stim_centre=stim_centre;']);
    eval(['results',num2str(nn2),'.W=W;']);
    eval(['results',num2str(nn2),'.t=t;']);
    eval(['results',num2str(nn2),'.x__=x__;']);
    eval(['results',num2str(nn2),'.y__=y__;']);
    eval(['results',num2str(nn2),'.d_Cont=d_Cont;']);
    eval(['results',num2str(nn2),'.d_Disc=d_Disc;']);
    eval(['results',num2str(nn2),'.tTGr_in=tTGr_in(1);']);
    eval(['results',num2str(nn2),'.gsyn_fromCont=gsyn_fromCont;']);
    eval(['results',num2str(nn2),'.gsyn_fromDisc=gsyn_fromDisc;']);
    eval(['results',num2str(nn2),'.tempo=tempo;']);
    if inh_block ==1
        results_fn = ['data2DMultispecies'];
    else
        results_fn = ['data2DMultispecies_Blockinh'];
    end
    eval(['save ' results_fn num2str(nn2) '_grc_grid_res_' mat2str(grc_grid_res) '.mat results',num2str(nn2)]);
end

% % Make movie
% % MakeVideo_A_regular(results_fn,num2str(nn2),[results_fn '_control']);
% eval(['control = results',num2str(nn2) ';']);
% 
% 
% % Set Goc->Grc to 0 i.e. Block inhibition
% gsyn_fromDisc=0;
% % re-run simulation
% tic
% tspan=[0:1:50]; %tspan=[0:5:1000];
% % f_multiscaleNewDEFMossyGolgi(0,W0);
% [t,W]=ode45('f_multiscaleNewDEFMossyGolgi',tspan,W0);
% tempo=toc
% 
% W(:,1:nn1)=3*W(:,1:nn1);
% 
% t=round(t*(10^13))/(10^12);
% if SaveData==1
%     eval(['results',num2str(nn2),'.coordinates1=coordinates1;']);
%     eval(['results',num2str(nn2),'.triangles1=triangles1;']);
%     eval(['results',num2str(nn2),'.nn1=nn1;']);
%     eval(['results',num2str(nn2),'.coordinates2=coordinates2;']);
% %     eval(['results',num2str(nn2),'.triangles2=triangles2;']);
%     eval(['results',num2str(nn2),'.nn2=nn2;']);
%     eval(['results',num2str(nn2),'.stim_centre=stim_centre;']);
%     eval(['results',num2str(nn2),'.W=W;']);
%     eval(['results',num2str(nn2),'.t=t;']);
%     eval(['results',num2str(nn2),'.x__=x__;']);
%     eval(['results',num2str(nn2),'.y__=y__;']);
%     eval(['results',num2str(nn2),'.d_Cont=d_Cont;']);
%     eval(['results',num2str(nn2),'.d_Disc=d_Disc;']);
%     eval(['results',num2str(nn2),'.gsyn_fromCont=gsyn_fromCont;']);
%     eval(['results',num2str(nn2),'.gsyn_fromDisc=gsyn_fromDisc;']);
%     eval(['results',num2str(nn2),'.tempo=tempo;']);
%     results_fn = ['data2DMultispecies_Blockinh'];
%     eval(['save ' results_fn num2str(nn2) '.mat results',num2str(nn2)]);
% end
% eval(['inhib_block = results',num2str(nn2) ';']);
% 
% % Make movie
% % MakeVideo_A_regular(results_fn,num2str(nn2),[results_fn '_blockinh']);
% 
% % Calc differnece of simulation results
% % MakeFigure_EI
% 
