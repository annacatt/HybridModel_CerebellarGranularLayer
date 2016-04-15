% Please cite **** TO BE ADDED ****

close all
clear all

% GrCs: granule cell nodes
% GoCs: Golgi cells
% MFTs: Mossy Fiber Terminals

rebuild = 1; %0 if WholeNetwork_data.mat has been already saved; 1 otherwise

% Defining coefficients of coupling terms
d_Disc=0.005; % d in Eq.(11)
d_Cont=0.005; % delta in Eq.(12)
gsyn_fromCont=0.8; % gamma_syn in Eq.(11)
gsyn_fromDisc=0.1; % g_syn in Eq.(12)


if rebuild == 1
    SaveData=1;

    rng(1,'twister'); 
    
    % Loading GoC population grid
    load grid1rect.mat; %GoC coordinates produced by triangular mesh generator BBTR, 
                        %described in (Barbera and Berrone,262 2008)
    
    % Building GrC population grid
    nn2_rad=88; % Number of GrC nodes per edge
    short_edge=1; % Length short edge
    long_edge=3; % Length long edge
    nn2=long_edge*nn2_rad^2;  % Total number of GrC nodes
    h=short_edge/nn2_rad;
    x__=linspace(0,1-h,nn2_rad);
    y__=linspace(0,long_edge-h,long_edge*nn2_rad);
    [x,y]=meshgrid(x__,y__);
    siz=size(x,2);
    x_=reshape(x,3*siz^2,1);
    y_=reshape(y,3*siz^2,1);
    coordinates2=[x_ y_]; % Coordinates of GrC nodes 
    

    % Cell densities [GrC GoC MFT]:
    celldens.mm3 = [4e6 9523 314717]; % #/mm3
    celldens.mm2 = celldens.mm3 / 10; % #/mm2
    
    % Network Area
    NNarea = short_edge * long_edge /4; % mm2
    n_MFT = celldens.mm2(3) * NNarea; % 23.000 MFTs so one MF each 13 GrC = 53/4
    
    % External input configuration (from MFTs)
    mf_grc_div = 53;
    mf_grc_radius = 0.030;
    num=randi([1 nn2],1,round(nn2*3/100));  % Select at random 3% of MFTs that input to GrC nodes
    MFT_Intensity=zeros(nn2,1);
    h=short_edge/nn2_rad;
    for i=1:length(num)
        [also, Intensity] = pointsincircle_flattensphere(coordinates2',mf_grc_radius,coordinates2(num(i),:),mf_grc_div);
        MFT_Intensity(also) = MFT_Intensity(also) + Intensity';
    end
    IGr_Mossy = MFT_Intensity ./ max(MFT_Intensity) .* 0.2;  % GrCs spot activation

    
    % Select at random 3% of MFTs that input to GoCs
    IGo_Mossy=zeros(nn1,1);
    numGo=round(linspace(1,nn1,nn1*3/100));
    IGo_Mossy(numGo)=0.1;
    
    tTGr_in=0*ones(nn2,1); % Time at wich MFTs start input to GrCs
    tTGo_in=10*ones(nn1,1); % Time at wich MFTs start input to GoCs
    

    %%%%%%%%%%%%%%%%
    % Connectivity %
    %%%%%%%%%%%%%%%%
    
    % 1) InChimFromCont = Phi in Eq.(11)
    heig=2/10;
    for k=1:nn1
        InChimFromCont(k).val=pointsinrectangle(coordinates2',heig,coordinates1(k,:)); 
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
        InChimFromDisc(k).val=pointsinrectangleSagittal(coordinates1',PFaxiswidth,SAGaxiswidth,coordinates2(k,:)); %input dal discreto sul continuo
        if length(InChimFromDisc(k).val)>40
            X=randperm(length(InChimFromDisc(k).val),40);
            InChimFromDisc(k).val=InChimFromDisc(k).val(X);
        end
    end
    [InChimFromDisc_sparse, InChimFromDisc_nfrac] = Conv_struc2matix(InChimFromDisc);
    InChimFromDisc_nfrac = 4./InChimFromDisc_nfrac;
    
    % The points 3) and 4) describe gap junctions between GoCs and ephaptic
    % coupling within GrC population. 
    
    % 3) InElecDisc = phi in Eq.(11), i.e. gap junctions among GoCs
    [InElecDisc_sparse,InElecDisc_n] = Conv_struc2matix(InElecDisc); %InElecDisc describes connections within the GoC population

    % 4) InElecDisc = psi in Eq.(12), i.e. ephaptic coupling in GrC population
    for k=1:nn2
        InElecCont(k).val=pointsincircle(coordinates2',h+1e-4,coordinates2(k,:)); 
    end
    [InElecCont_sparse,InElecCont_n] = Conv_struc2matix(InElecCont); %InElecCont describes connections within the GrC population
    
    
    % Defining parameters
    T=0.5*[ones(1,nn1)';ones(1,nn2)'];  % Threhsold v_T
    vIsyn=-0.2; % Reversal potential GoCs
    vIsyn=vIsyn*ones(nn2,1);
    omegaEsyn=0.93; % Reversal potential GrC nodes
    omegaEsyn=omegaEsyn*ones(nn1,1);
    Wsyn=[vIsyn;omegaEsyn]; % Reversal potential whole network
    
    %%%%%%%%%%%%%%%%%%
    % Initialization %
    %%%%%%%%%%%%%%%%%%
    
    v0=spalloc(3*nn1,1,nn1); % GoC population initial datum
    omega0=spalloc(3*nn2,1,10); % GrC population initial datum 
    W0=[v0;omega0];
    
    save WholeNetwork_data.mat;
else
    load WholeNetwork_data.mat;
end

%%%%%%%%%%%%%%%%%%%%%%%%%
% Integrating the model %
%%%%%%%%%%%%%%%%%%%%%%%%%

tic
tspan=0:100; %1000 = 1 second simulation
[t,W]=ode45(@(t,w)fun_WholeNetwork(t,w,d_Cont,d_Disc,nn1,nn2,T,gsyn_fromCont,...
    gsyn_fromDisc,Wsyn,tTGr_in,IGr_Mossy,IGo_Mossy,tTGo_in,InElecCont_sparse,...
    InElecCont_n,InElecDisc_sparse,InElecDisc_n,InChimFromCont_nfrac,...
    InChimFromDisc_nfrac,InChimFromCont_sparse,InChimFromDisc_sparse),tspan,W0);
time=toc


coordinates1=coordinates1';
coordinates2=coordinates2';
W(:,1:nn1)=3*W(:,1:nn1);

t=round(t*(10^13))/(10^12);

if SaveData==1
    eval(['results',num2str(nn2),'.coordinates1=coordinates1;']);
    eval(['results',num2str(nn2),'.triangles1=triangles1;']);
    eval(['results',num2str(nn2),'.nn1=nn1;']);
    eval(['results',num2str(nn2),'.coordinates2=coordinates2;']);
    eval(['results',num2str(nn2),'.nn2=nn2;']);
    eval(['results',num2str(nn2),'.W=W;']);
    eval(['results',num2str(nn2),'.t=t;']);
    eval(['results',num2str(nn2),'.x__=x__;']);
    eval(['results',num2str(nn2),'.y__=y__;']);
    eval(['results',num2str(nn2),'.d_Cont=d_Cont;']);
    eval(['results',num2str(nn2),'.d_Disc=d_Disc;']);
    eval(['results',num2str(nn2),'.gsyn_fromCont=gsyn_fromCont;']);
    eval(['results',num2str(nn2),'.gsyn_fromDisc=gsyn_fromDisc;']);
    eval(['results',num2str(nn2),'.time=time;']);
    eval(['save WholeNetwork_dynamics_',num2str(nn2),'.mat results',num2str(nn2)]);
end