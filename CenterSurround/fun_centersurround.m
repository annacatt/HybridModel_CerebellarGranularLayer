function wdot = fun_centersurround(t,w,nn1,nn2,d_Cont,d_Disc,InElecDisc,InElecCont,T,gsyn_fromCont,gsyn_fromDisc,InChimFromCont,...
    InChimFromDisc,Wsyn,IGr_Mossy,IGo_Mossy,InElecCont_sparse,InElecCont_n,...
    InElecDisc_sparse,InElecDisc_n,InChimFromCont_nfrac,InChimFromDisc_nfrac,InChimFromCont_sparse,...
    InChimFromDisc_sparse,tT_end,tTGr_in,tTGr_in_end,tTGo_in,tTGo_in_end)

a=0.25; b=.001; g=.003; % FitzHugh-Nagumo parameters

% Rates at which synapse turns on and off, respectively
alpha=0.9; % in Eq.(5)
beta=0.1; 
alpha_i=0.9/50; % in Eq.(8)
beta_i=0.1/50;

% Heaviside functions
ditGo_in=t*ones(nn1,1)-tTGo_in;
htGo_in=1./(1+exp(-1000*ditGo_in)); % In the GoC population
ditGr_in=t*ones(nn2,1)-tTGr_in;
htGr_in=1./(1+exp(-1000*ditGr_in)); % In the GrC population

% Variables
Wdisc_e=w(1:nn1,1); % Membrane potential (v) of GoCs
Wcont_e=w(3*nn1+1:3*nn1+nn2,1); % Membrane potential (omega) of GrC nodes
W_e=[Wdisc_e;Wcont_e];
s=w(2*nn1+1:3*nn1,1); % Variable s of GoCs
sigma=w(3*nn1+2*nn2+1:3*nn1+3*nn2,1); % Variable sigma of GrC nodes

% Electrical coupling terms
s_e = -InElecDisc_n' .* Wdisc_e + (Wdisc_e' * InElecDisc_sparse)';
differ_disc=d_Disc*s_e; % Gap junctions among GoCs
sCont_e = -InElecCont_n' .* Wcont_e + (Wcont_e' * InElecCont_sparse)'; %(InElecCont incorporates dirichlet boundary conditions)
differ_cont=d_Cont*sCont_e; % Ephaptic coupling in GrC population 

% Chimical coupling terms
s_fromCont_c = InChimFromCont_nfrac' .* (sigma' * InChimFromCont_sparse')';
I_fromCont=-gsyn_fromCont.*s_fromCont_c.*(w(1:nn1)-Wsyn(nn2+1:nn1+nn2)); % Phi in Eq.(11)
s_fromDisc_c = InChimFromDisc_nfrac' .* (s' * InChimFromDisc_sparse')'; 
I_fromDisc=-gsyn_fromDisc.*s_fromDisc_c.*(w(3*nn1+1:3*nn1+nn2)-Wsyn(1:nn2)); % Psi in Eq. (12)

di=W_e-T;
h=1./(1+exp(-1000*di)); % Heaviside function in Eqs.(5)+(8)

% Eqs.(5)
wdot(1:nn1)=-(w(1:nn1)).*(a-(w(1:nn1))).*(1-(w(1:nn1)))-w(nn1+1:2*nn1)+differ_disc+I_fromCont+IGo_Mossy(1:nn1,1).*htGo_in;
wdot(nn1+1:2*nn1)=b*w(1:nn1)-g*w(nn1+1:2*nn1);
wdot(2*nn1+1:3*nn1)=alpha*(1-w(2*nn1+1:3*nn1)).*h(1:nn1)-beta*w(2*nn1+1:3*nn1);
% Eqs.(8) % boundary conditions (homogeneous Dirichlet) are included in the topology matrices
wdot(3*nn1+1:3*nn1+nn2)=-(w(3*nn1+1:3*nn1+nn2)).*(a-(w(3*nn1+1:3*nn1+nn2))).*(1-(w(3*nn1+1:3*nn1+nn2)))-w(3*nn1+nn2+1:3*nn1+2*nn2)+differ_cont+I_fromDisc+IGr_Mossy(1:nn2,1).*htGr_in;
wdot(3*nn1+nn2+1:3*nn1+2*nn2)=b*w(3*nn1+1:3*nn1+nn2)-g*w(3*nn1+nn2+1:3*nn1+2*nn2);
wdot(3*nn1+2*nn2+1:3*nn1+3*nn2)=alpha_i*(1-w(3*nn1+2*nn2+1:3*nn1+3*nn2)).*h(nn1+1:end)-beta_i*w(3*nn1+2*nn2+1:3*nn1+3*nn2);
wdot=wdot';
