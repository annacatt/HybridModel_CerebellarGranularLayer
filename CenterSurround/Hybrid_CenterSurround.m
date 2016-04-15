function Hybrid_CenterSurround(grc_grid_res,inh_block)

% global nn1 nn2 d_Cont d_Disc InElecDisc InElecCont T gsyn_fromCont gsyn_fromDisc InChimFromCont
% global InChimFromDisc Wsyn IGr_Mossy IGo_Mossy InElecCont_sparse InElecCont_n 
% global InElecDisc_sparse InElecDisc_n InChimFromCont_nfrac InChimFromDisc_nfrac InChimFromCont_sparse 
% global InChimFromDisc_sparse tT_end tTGr_in tTGr_in_end
% global tTGo_in tTGo_in_end

% GrCs: granule cell nodes
% GoCs: Golgi cells
% MFTs: Mossy Fiber Terminals

d_Cont=0; % d in Eq.(11)
d_Disc=0; % delta in Eq.(12)
gsyn_fromCont=0.8; % gamma_syn in Eq.(11)
gsyn_fromDisc=0.05; % It was 0.05 but in the old CS I had 26 Goc concentrated in a small area so the grcs were much more inhibited


rebuild = 0; %0 if CS.mat has been already saved; 1 otherwise

if rebuild == 1
    SaveData=1;
    
    
    if inh_block == 1
        gsyn_fromDisc=0.1; % g_syn in Eq.(12)
    else
        gsyn_fromDisc= inh_block;
    end
    
    % Loading GoC population grid
    load grid1rect_CS_200x200_40Goc.mat; %200um x 200um
    
    % Building GrC population grid
    nn2_rad= grc_grid_res; % Number of GrC nodes per edge
    edge = 0.2; % 1 == 500 um -> 0.2 = 100um
    nn2=nn2_rad^2; % Total number of GrC nodes
    h=edge/(nn2_rad-1);
    x__=linspace(0.4,0.4+edge,nn2_rad);
    y__=linspace(0.4,0.4+edge,nn2_rad);
    [x,y]=meshgrid(x__,y__);
    siz=size(x,2);
    x_=reshape(x,siz^2,1);
    y_=reshape(y,siz^2,1);
    coordinates2=[x_ y_]; % Coordinates of GrC nodes 
    
    
    % External input configuration (from MFTs)  
    
    % Each active Mossy Fiber contacts mf_grc_div nodes: it can contact more than
    % once the same node as one node represents more than one grc.
    % We assume that a certain number of mf (num_of_active_mfs) is
    % transmitting the stimulus and then pull together their targets.
    IGr_Mossy=zeros(nn2,1);     % GrCs spot activation
    mf_grc_radius = 0.090; % max grc dend length 30 um
    stim_centre =[ones(10,1)*.5 ones(10,1)*.5]  + rand(10,2).*0.09;
    mf_grc_input_strength = 0.4;
    spot_decay = 0.4;
    [nn2_center,Intensity_mf_grc] = pointsincircle_gauss_intensity(coordinates2',mf_grc_radius,[0.5 0.5],spot_decay);
    Intensity_mf_grc = (Intensity_mf_grc - min(Intensity_mf_grc));
    Intensity_mf_grc = Intensity_mf_grc ./ max(Intensity_mf_grc);
    Intensity_mf_grc = Intensity_mf_grc + 0.05 * rand(size(Intensity_mf_grc));
    IGr_Mossy(nn2_center)= Intensity_mf_grc .*  mf_grc_input_strength;  % Get the list of nodes contacted by the 10 mfs.

    % GoCs spot activation
    IGo_Mossy=zeros(nn1,1); % GoC nodes spot activation
    mf_goc_input_strength = 0.1;
    mf_goc_radius = mf_grc_radius+0.1; % This has to consider the radius of the GoC basolateral dendrites. It was set to 100 um in Solinas2010
    [nn1_center,Intensity_mf_goc] = pointsincircle_gauss_intensity(coordinates1',mf_goc_radius,[0.5 0.5],spot_decay);
    IGo_Mossy(nn1_center)= mf_goc_input_strength; 
    
    %%%%%%%%%%%%%%%%
    % Connectivity %
    %%%%%%%%%%%%%%%%

    % 1) InChimFromCont = Phi in Eq.(11)
    heig=2/10;
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
    
    % 2) InChimFromDisc = Psi in Eq.(12)
    PFaxiswidth=0.18; % width of rectancle along the PF axis as 1 corresponds to 500 um 1/2.8=180um
    PFaxiswidth = PFaxiswidth * 2; % double it so that from the ref of 0.5 mm it goes to the 1 short edge lenth of the rectangle
    SAGaxiswidth= 0.65; % width of rectangle along the sagittal axis
    SAGaxiswidth= SAGaxiswidth * 2; % double it so that from the ref of 0.5 mm it goes to the 1 short edge lenth of the rectangle
    
    for k=1:nn2
        %     InChimFromDisc(k).val=pointsinrectangleSagittal2(coordinates1',heigSag,coordinates2(k,:)); %input dal discreto sul continuo
        InChimFromDisc(k).val=pointsinrectangleSagittal(coordinates1',PFaxiswidth,SAGaxiswidth,coordinates2(k,:)); %input dal discreto sul continuo
        if length(InChimFromDisc(k).val)>40%/fatt
            X=randperm(length(InChimFromDisc(k).val),40); %round(40/fatt)
            InChimFromDisc(k).val=InChimFromDisc(k).val(X);
        end
    end
    [InChimFromDisc_sparse, InChimFromDisc_nfrac] = Conv_struc2matix(InChimFromDisc);
    InChimFromDisc_nfrac = 4./InChimFromDisc_nfrac;
    
    % The points 3) and 4) describe gap junctions between GoCs and ephaptic
    % coupling within GrC population. 
    
    % 3) InElecDisc = phi in Eq.(11), i.e. gap junctions among GoCs
    [InElecDisc_sparse,InElecDisc_n] = Conv_struc2matix(InElecDisc);
 
    % 4) InElecDisc = psi in Eq.(12), i.e. ephaptic coupling in GrC population
    for k=1:nn2
        InElecCont(k).val=pointsincircle(coordinates2',h+0.001,coordinates2(k,:)); %input dal continuo sul discreto
    end
    [InElecCont_sparse,InElecCont_n] = Conv_struc2matix(InElecCont);
    
    
    % Defining parameters
    T=0.5*[ones(1,nn1)';ones(1,nn2)']; % Threhsold v_T
    vIsyn=-0.2;  % Reversal potential GoCs
    omegaEsyn=0.93;%reversal potential granular cells
    vIsyn=vIsyn*ones(nn2,1);
    omegaEsyn=omegaEsyn*ones(nn1,1); % Reversal potential GrC nodes
    Wsyn=[vIsyn;omegaEsyn];

    tTGr_in = [10*ones(nn2,1) 100e3*ones(nn2,1)]; % Time at which input MFTs-->GrC population starts
    tTGr_in_end = [15*ones(nn2,1) 110e3*ones(nn2,1)]; % Time at which input MFTs-->GrC population ends 
    tTGo_in = 10*ones(nn1,1); % Time at which input MFTs-->GoCs starts
    tTGo_in_end = 15*ones(nn1,1); % Time at which input MFTs-->GoCs ends
    
    %%%%%%%%%%%%%%%%%%
    % Initialization %
    %%%%%%%%%%%%%%%%%%
    
    v0=spalloc(3*nn1,1,nn1); % GoC population initial datum
    omega0=spalloc(3*nn2,1,10); % GrC population initial datum 
    W0=[v0;omega0];
    

    save CS.mat;
else
    load CS.mat;
end

coordinates1=coordinates1';
coordinates2=coordinates2';

tic
tspan=[0:0.1:50]; 
[t,W]=ode45(@(t,w)fun_centersurround(t,w,nn1,nn2,d_Cont,d_Disc,InElecDisc,InElecCont,T,gsyn_fromCont,gsyn_fromDisc,InChimFromCont,...
InChimFromDisc,Wsyn,IGr_Mossy,IGo_Mossy,InElecCont_sparse,InElecCont_n,... 
InElecDisc_sparse,InElecDisc_n,InChimFromCont_nfrac,InChimFromDisc_nfrac,InChimFromCont_sparse,... 
InChimFromDisc_sparse,tT_end,tTGr_in,tTGr_in_end,tTGo_in,tTGo_in_end),tspan,W0);
time=toc



W(:,1:nn1)=3*W(:,1:nn1);

t=round(t*(10^13))/(10^12);
if SaveData==1
    eval(['results',num2str(nn2),'.coordinates1=coordinates1;']);
    eval(['results',num2str(nn2),'.triangles1=triangles1;']);
    eval(['results',num2str(nn2),'.nn1=nn1;']);
    eval(['results',num2str(nn2),'.coordinates2=coordinates2;']);
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
    eval(['results',num2str(nn2),'.time=time;']);
    if inh_block ==1
        results_fn = ['data2DMultispecies'];
    else
        results_fn = ['data2DMultispecies_Blockinh'];
    end
    eval(['save ' results_fn num2str(nn2) '_grc_grid_res_' mat2str(grc_grid_res) '.mat results',num2str(nn2)]);
end

