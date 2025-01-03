% Creates parameter files from choice of parameter set specified in Opt
% and saves it in .../DATA/Params.mat

function create_params(Opt,filename)
%----------
% NETWORKS
%----------
% Sim
Sim.t_Start=-1;
Sim.t_End=1;
Sim.dt_step=0.0001;
Stimulus.option='on'; % 'off'
Stimulus.input='Const'; % 'off'
%----------------
% DEFAULT STIMULI
%----------------
% STIMULUS
% numel(feat) is the number of stimuli delivered in the same trial
scnt=0;
% TASTE (generic)
scnt=scnt+1;
tau_som=[0.05 0.2];
feat(scnt).name='som'; % somatosensory stimulus, same for all tastes
feat(scnt).interval=[0 Sim.t_End]; % stimulus interval
gain=0.2;
feat(scnt).gain=gain;
feat(scnt).profile=@(t)(1/(tau_som(2)-tau_som(1)))*(exp(-t/tau_som(2))-exp(-t/tau_som(1))); %@(x)1; % @(x)double_exp
% feat(scnt).profile=@(t)1; %@(x)1; % @(x)double_exp
feat(scnt).selectivity='overlap'; %'one'
feat(scnt).connectivity=0.5; % fraction of neurons targeted
% TASTE (specific)
scnt=scnt+1;
%     tau_th=[0.5 1]; % rise and decay time (s)
feat(scnt).name='US'; % unconditioned stimulus (taste)
feat(scnt).interval=[0 Sim.t_End]; % stimulus interval
gain=0.1;
feat(scnt).gain=gain;
%     feat(scnt).profile=@(t)(1/(tau_th(2)-tau_th(1)))*(exp(-t/tau_th(2))-exp(-t/tau_th(1))); %@(x)1; % @(x)double_exp
feat(scnt).profile=@(t)t; %@(x)1; % @(x)double_exp
feat(scnt).selectivity='one'; %'one';% 
feat(scnt).connectivity=0.5; % fraction of neurons targeted inside cluster
% CUE
tau_cue=[0.2 1];
scnt=scnt+1;
feat(scnt).name='CS'; % conditioned stimulus (cue)
feat(scnt).interval=[-0.5 Sim.t_End]; % stimulus interval
gain=0.2;
feat(scnt).gain=gain;
feat(scnt).profile=@(t)(1/(tau_cue(2)-tau_cue(1)))*(exp(-t/tau_cue(2))-exp(-t/tau_cue(1))); %@(x)1; % @(x)double_exp
%     feat(scnt).profile=@(t)gain; %@(x)1; % @(x)double_exp
feat(scnt).selectivity='inh'; %'one';% 'overlap'=target pops in overlap
%     c=zeros(1,p+1); c(2)=1;
%     feat(scnt).overlap=c; % if .selectivity='overlap', array c of length p+1 (# of clusters,background), c(i)=0 if i-th overlap
%                          % population not targeted, 0<c(i)<1 fraction of
%                          % i-th population targeted. c(i+1)=background pop
feat(scnt).connectivity=0.5; % fraction of neurons targeted inside cluster
% GAUSSIAN CUE up and down
scnt=scnt+1;
feat(scnt).name='CSgauss'; % conditioned stimulus (cue)
feat(scnt).interval=[-0.5 Sim.t_End]; % stimulus interval
gain=0.1; % SD of gaussian
feat(scnt).gain=gain;
% feat(scnt).profile=@(t)(1/(tau_cue(2)-tau_cue(1)))*(exp(-t/tau_cue(2))-exp(-t/tau_cue(1))); %@(x)1; % @(x)double_exp
    feat(scnt).profile=@(t)gain; %@(x)1; % @(x)double_exp
feat(scnt).selectivity='exc'; %'one';% 'overlap'=target pops in overlap
%     c=zeros(1,p+1); c(2)=1;
%     feat(scnt).overlap=c; % if .selectivity='overlap', array c of length p+1 (# of clusters,background), c(i)=0 if i-th overlap
%                          % population not targeted, 0<c(i)<1 fraction of
%                          % i-th population targeted. c(i+1)=background pop
feat(scnt).connectivity=0.50; % fraction of neurons targeted inside cluster%
%
% MIXED US (for stats sims)
scnt=scnt+1;
feat(scnt).name='USdown'; % conditioned stimulus (cue)
feat(scnt).interval=[-0.5 Sim.t_End]; % stimulus interval
gain=-0.2;
feat(scnt).gain=gain;
feat(scnt).profile=@(x)1; % @(x)double_exp
%     feat(scnt).profile=@(t)gain; %@(x)1; % @(x)double_exp
feat(scnt).selectivity='exc'; %'one';% 'overlap'=target pops in overlap
%     c=zeros(1,p+1); c(2)=1;
%     feat(scnt).overlap=c; % if .selectivity='overlap', array c of length p+1 (# of clusters,background), c(i)=0 if i-th overlap
%                          % population not targeted, 0<c(i)<1 fraction of
%                          % i-th population targeted. c(i+1)=background pop
feat(scnt).connectivity=0.25; % fraction of neurons targeted inside cluster
%
scnt=scnt+1;
feat(scnt).name='USup'; % conditioned stimulus (cue)
feat(scnt).interval=[-0.5 Sim.t_End]; % stimulus interval
gain=0.2;
feat(scnt).gain=gain;
feat(scnt).profile=@(x)1; % @(x)double_exp
%     feat(scnt).profile=@(t)gain; %@(x)1; % @(x)double_exp
feat(scnt).selectivity='exc'; %'one';% 'overlap'=target pops in overlap
%     c=zeros(1,p+1); c(2)=1;
%     feat(scnt).overlap=c; % if .selectivity='overlap', array c of length p+1 (# of clusters,background), c(i)=0 if i-th overlap
%                          % population not targeted, 0<c(i)<1 fraction of
%                          % i-th population targeted. c(i+1)=background pop
feat(scnt).connectivity=0.25; % fraction of neurons targeted inside cluster%
%
if strcmp(Opt,'Small2_stats_bernstein')
    % Creates parameter for network with E but not I clusters
    %------------------
    % NETWORK OPTIONS
    %------------------
    % Network.clust='hom'; % homogeneous cluster size
    Network.syn='SingleExp';
    Network.clust='hom'; % heterogeneous cluster size
    Network.clust_std=0.0; % het cluster size: cluster size is sampled from a gaussian with SD=Ncluster*Network.std
    N=800; % network size
    N_e = N*4/5; % E neurons
    N_i = round((N/5)*[0.33,0.33,0.34]); % L2/3 in V1 from Kim et al (Osten lab) Cell 2017
    ni_e = 5; %
    ni_i = 10*[1,1,1]; %
    %--------------------
    % NEURON PARAMETERS
    %--------------------
    tau_i = .02*[1,1,1];     % 0.01
    tausyn_e=0.004;
    tausyn_i=tausyn_e*ones(1,numel(N_i));
    
    Scale=(5000/N)^(1/2);
    %------------------
    % SYNAPTIC WEIGHTS
    %------------------
    % syn weights are drawn from a gaussian distribution with std delta and
    % mean Jab
    % NOTE:
    % - large gPE gives the usual EI network limit! gPE is a
    % "flat direction, one can vary it by orders of magnitude!
    % - g can vary a lot too. large g has selective activity
    % with coupled E-VIP, which is wrong! It should be E-SST!
    % - TO DO: write PSC and Vm correctly for all cell types!
    Jee = Scale*0.012;
    Jii = Scale^2*0.055*[1,1,1;1,1,1;1,1,1];
    Jie = Scale*0.035*[1;1;1];
    Jei = Scale^2*0.06*[1,1,1];      
    delta = 0.; % SD of synaptic weights distribution, eg: ~Jee*(1+delta*randn(N_e))
    Network.delta = delta; % SD of synaptic weights: Jee*(1+delta*randn(N_e)) Larger delta -> clusters look more async in the high state
    Network.deltaEI = delta; % SD of synaptic weights: Jee*(1+delta*randn(N_e)) Larger delta -> clusters look more async in the high state
    Network.deltaIE = delta; % SD of synaptic weights: Jee*(1+delta*randn(N_e)) Larger delta -> clusters look more async in the high state
    Network.deltaII = delta; % SD of synaptic weights: Jee*(1+delta*randn(N_e)) Larger delta -> clusters look more async in the high state
    %------------------
    % CLUSTER PARAMETERS
    %------------------
    %                 Jplus =10.4; % intra-cluster potentiation factor for E clusters
    %                 Jplus =10.4; % intra-cluster potentiation factor for E clusters
    Jplus =8; % intra-cluster potentiation factor for E clusters
    bgr=0.65; % fraction of background excitatory neurons (unclustered)
    Ncluster=100; % average # neurons per cluster
    p = round(N_e*(1-bgr)/Ncluster); % # of clusters
    f = (1-bgr)/p;
    %------------------
    % THRESHOLDS
    %------------------
%     theta_e=6.00918; % exc threshold potential
%     theta_i=[29.0197,46.4366,7.55815]; % inh threshold potentials
    theta_e=152.199; % exc threshold potential
    theta_i=289.975*ones(1,3); % inh threshold potentials
    % reset potentials
    He = 0;%
    Hi = zeros(1,3);%
    %------------------
    % EXTERNAL CURRENT
    %------------------
    % default external currents
    ni_ext = 7; % 7;
    ni_ext_e = 7;
    ni_ext_i = 7*[1,1,1];
    %------------------
    % CONNECTIVITY PARAMETERS
    %------------------
    p_ee = 0.2; % 1600;         % 1600
    p_ie = 0.5*[1;1;1]; % 1600;         % 1600
    p_ii = 0.5*[1,1,1;
        1,1,1;
        1,1,1]; % 400;         % 400
    p_ei = 0.5*[1,1,1]; % 400;         % 400
    p_ext = 0.2; % 1600;
    %------------------
    % EXTERNAL BIAS
    %------------------
    % - low bias: oscillating; high bias: asynchronous
    % higher J_ext -> higher Jplus to get same dynamics
    Jie_ext=2*Scale*0.0915*[1,1,1];% normalized to ni_ext=2
    Jee_ext=2*Scale*0.1027; %
    nn = 5;                 % number of populations
elseif  strcmp(Opt,'Small2_stats_V1')
    % network options
    Network.syn='SingleExp';
    Network.clust='hom';%'het'; % cluster: homogeneous='hom', heterogeneous='het'
    Network.clust_syn='';
    Network.clust_std=0.0; % het clusters: take a gaussian distribution with variance=mean*Network.clust_var
    N=800;
    N_e = N*4/5; % 16000
%     N_i = N/5*[0.5,0.25,0.25]; %4000
%     N_i = N/5*[0.25,0.25,0.5]; % L2/3 in V1 from Kim et al (Osten lab) Cell 2017 
    N_i = (N/5)*[0.4,0.4,0.2]; % L2/3 in V1 from Kim et al (Osten lab) Cell 2017 
    Scale=(5000/N)^(1/2);
    % global spontaneous firing rates (neeed to fix thresholds)
    ni_e = 5; % 3.6;   %6.6   % AB97: 3.0
    ni_i = 10*[1,0.5,0.5]; % 5.2;   %8.2   % AB97: 4.2
    % external currents
%     ni_ext = 7;
    ni_ext_e = 5;
    ni_ext_i = 10*[1,0,0];

    tau_arp = .005;  % .002 (Unita': sec)
    tau_i = .02*[1,1,1];     % 0.01
    tau_e = .02;%15;     % 0.02

    delta = 0.0; % 0.;    
    x = 1.;
    % e=pyr, i=PV,SST,VIP
    % take values from dipoppa, scanziani, 2020 biorxiv 
%     Jee = Scale*0.05; 
%     Jei = Scale*(0.06)*[1,1,0];    
%     Jii = Scale*(0.055)*[1,1,0;
%                          0,0,1 ;
%                          2,2,0];
%     Jie = Scale*(0.035/6)*[1;1;4];
% 

% these values are ok
%     Jee = Scale*0.02; 
%     Jei = Scale*(0.15)*[1.1,2,0];    
%     Jii = Scale*(0.13)*[1,2.2,0;
%                          0,0,2 ;
%                          3,2.1,0];
%     Jie = Scale*(0.035/2)*[2;2.1;6];
g=3;
Jee = g*Scale*0.03;
Jei = g*Scale*(0.4).*[0.625,0.625,2.5].*[1,1,0];
Jii = g*Scale*(0.15)*[0.625,0.625,2.5].*[1,1,0;
                    0,0,  0.2 ;
                    1,1,0];
Jie = g*Scale*(0.03)*[1;2;1];

%     % these weights below with Jplus=1 and give spontaneous rates: 13.2443  13.29  20.3668  2.50531  45.7198
%     Jee = Scale*0.05; 
%     Jei = Scale*(0.15)*[1,2,0];    
%     Jii = Scale*(0.15)*[1,2,0;
%                          0,0,2 ;
%                          2,4,0];
%     Jie = Scale*(0.035/2)*[2;2;6];
%     Jei=Jei.*(1+[-0.1351,    0.0876,   -0.0011]);
%     Jie=Jie.*(1+[0.0645;0.0937;0.0620]);
%     Jii=Jii.*(1+[0.2375,   -0.2134,    0.1139; 0.1308,    0.4544,    0.3877; -0.1031,    0.0742,   -0.3850]);
%     Jei=Jei+0.1*[-0.1351,    0.0876,   -0.0011];
%     Jie=Jie+0.1*[0.0645;0.0937;0.0620];
%     Jii=Jii+[0.1*0.2375,   -0.2134,    0.1139; 0.1308,    0.4544,    0.3877; -0.1031,    0.0742,   -0.3850];
    
    Jplus = 1.;
    %
%     theta_e=50.; %22.8385;%21.7124; 
%     theta_i=50.*[1,1,1] ;%theta_i=11.1897;%11.0887;  27.4304, (I) 58.6757,70.9415,14.0214,
    theta_e=27.4304; %22.8385;%21.7124; 
    theta_i=[58.6757,70.9415,14.0214] ;%theta_i=11.1897;%11.0887;  
    Jplus = 5.;
%     gam = 0.5;%1.125; % 0.5;  %1.3
    gam = 0.5;%1.125; % 0.5;  %1.3

    He = 0;%0.5;
    Hi = 0*[1,1,1];%0.5;

    Qpiu = 0.6;  %0.6
    Qmeno = 0.3;  %0.3
    rho = 2.75;    %(5)
    Jzero =11.5;  %1.

%     % connection probability (from Wang and Huang, 2017)
%     Cee = N_e*0.2; % 1600;         % 1600
%     Cie = N_e*0.5*[1;1;1]; % 1600;         % 1600
%     Cii = 0.5*N_i.*[1,.85,0;
%                     0,0,0.5;
%                     0.2,0.5,0]; % 400;         % 400
%     Cei = 0.2*N_i.*[1,1,0.]; % 400;         % 400
%     Cext = (N_e)*0.2; % 1600; 
        % connection probability
    pee_matrix = 0.2; % 1600;         % 1600
    pie_matrix = 0.5*[1;1;1]; % 1600;         % 1600
    pii_matrix = 0.5*[1,1,0;
                    0,0,1;
                    1,1,0]; % 400;         % 400
    pei_matrix = 0.2*[1,1,0.]; % 400;         % 400
    pext_matrix = 0.2; % 1600; 

    % external currents
    Jie_ext=2*Scale*0.0915*[1,1,1];% normalized to ni_ext=2
    Jee_ext=2*Scale*0.1027; % 
    
    bgr=0.65; % N_e fraction of background neurons
    p = round(N_e*(1-bgr)/100); % keep 100 neurons per cluster
    f = (1-bgr)/p;       % 0.09
%     p_actual = 4;
%     nstories = 10;    % see monte.c
    nn = 5;                 % number of populations
    nfocus = 0;             % Number of populations in focus
    tausyn_e=0.002;
    tausyn_i=tausyn_e;

elseif  strcmp(Opt,'Small2InhDens')
    % network options
    Network.syn='SingleExp';
    Network.clust='hom';%'het'; % cluster: homogeneous='hom', heterogeneous='het'
    Network.clust_syn='';
    Network.clust_std=0.0; % het clusters: take a gaussian distribution with variance=mean*Network.clust_var
    N=800;
    N_e = N*4/5; % 16000
    N_i = (N/5)*[0.4,0.4,0.2]; % L2/3 in V1 from Kim et al (Osten lab) Cell 2017 
    Scale=(5000/N)^(1/2);
    % global spontaneous firing rates (neeed to fix thresholds)
    ni_e = 5; % 3.6;   %6.6   % AB97: 3.0
    ni_i = 10*[1,0.5,0.5]; % 5.2;   %8.2   % AB97: 4.2
    % external currents
    ni_ext_e = 5;
    ni_ext_i = 10*[1,1,1];

    tau_arp = .005;  % .002 (Unita': sec)
    tau_i = .02*[1,1,1];     % 0.01
    tau_e = .02;%15;     % 0.02

    delta = 0.0; % 0.;    
    x = 1.;
gE=2; % these values are great
gI=6; % these values are great
% Jee = g*Scale*0.04;
Jee = gE*Scale*0.1;
Jei = gI*Scale*(0.4).*[0.625,0.625,2.5].*[1,1,0];
Jii = gI*Scale*(0.15)*[0.625,0.625,2.5].*[1,1,0;
                    0,0,  0.2 ;
                    0,1,0];
Jie = gE*Scale*(0.03)*[1;2;1];

% g=8;
% Jee = g*Scale*0.06;
% Jei = g*Scale*(0.4).*[0.625,0.625,2.5].*[6,1,0];
% Jii = g*Scale*(0.15)*[0.625,0.625,2.5].*[6,1,0;
%                     0,0,  0.2 ;
%                     6,1,0];
% Jie = g*Scale*(0.05)*[2;2;1];


    
    Jplus = 1.;
    %
    theta_e=27.4304; %22.8385;%21.7124; 
    theta_i=[58.6757,70.9415,14.0214] ;%theta_i=11.1897;%11.0887;  
    Jplus = 5.;
    gam = 0.5;%1.125; % 0.5;  %1.3

    He = 0;%0.5;
    Hi = 0*[1,1,1];%0.5;

    Qpiu = 0.6;  %0.6
    Qmeno = 0.3;  %0.3
    rho = 2.75;    %(5)
    Jzero =11.5;  %1.

        % connection probability
    pee_matrix = 0.2; % 1600;         % 1600
    pie_matrix = 0.5*[1;1;1]; % 1600;         % 1600
    pii_matrix = 0.5*[1,1,0;
                    0,0,1;
                    1,1,0]; % 400;         % 400
    pei_matrix = 0.2*[1,1,0.]; % 400;         % 400
    pext_matrix = 0.2; % 1600; 

    % external currents
    Jie_ext=2*Scale*0.0915*[1,1,1];% normalized to ni_ext=2
    Jee_ext=2*Scale*0.1027; % 
    
    bgr=0.65; % N_e fraction of background neurons
    p = round(N_e*(1-bgr)/100); % keep 100 neurons per cluster
    f = (1-bgr)/p;       % 0.09
    nn = 5;                 % number of populations
    nfocus = 0;             % Number of populations in focus
    tausyn_e=0.002;
    tausyn_i=tausyn_e;

elseif  strcmp(Opt,'Small2Reversal')
    % network options
    Network.syn='SingleExp';
    Network.clust='hom';%'het'; % cluster: homogeneous='hom', heterogeneous='het'
    Network.clust_syn='';
    Network.clust_std=0.0; % het clusters: take a gaussian distribution with variance=mean*Network.clust_var
    N=800;
    N_e = N*4/5; % 16000
%     N_i = N/5*[0.5,0.25,0.25]; %4000
%     N_i = N/5*[0.25,0.25,0.5]; % L2/3 in V1 from Kim et al (Osten lab) Cell 2017 
    N_i = (N/5)*[0.4,0.4,0.2]; % L2/3 in V1 from Kim et al (Osten lab) Cell 2017 
    Scale=(5000/N)^(1/2);
    % global spontaneous firing rates (neeed to fix thresholds)
    ni_e = 5; % 3.6;   %6.6   % AB97: 3.0
    ni_i = 10*[1,0.5,0.5]; % 5.2;   %8.2   % AB97: 4.2
    % external currents
%     ni_ext = 7;
    ni_ext_e = 5;
    ni_ext_i = 10*[1,1,1];

    tau_arp = .005;  % .002 (Unita': sec)
    tau_i = .02*[1,1,1];     % 0.01
    tau_e = .02;%15;     % 0.02

    delta = 0.0; % 0.;    
    x = 1.;
    % e=pyr, i=PV,SST,VIP
    % take values from dipoppa, scanziani, 2020 biorxiv 
%     Jee = Scale*0.05; 
%     Jei = Scale*(0.06)*[1,1,0];    
%     Jii = Scale*(0.055)*[1,1,0;
%                          0,0,1 ;
%                          2,2,0];
%     Jie = Scale*(0.035/6)*[1;1;4];
% 
% gE=6; gI=9;
% Jee = gE*Scale*0.065;
% Jei = gI*Scale*(0.4).*[0.625,0.625,2.5].*[1,1,0];
% Jii = gI*Scale*(0.15)*[0.625,0.625,2.5].*[1,1,0;
%                     0,0,  0.4 ;
%                     1,0.5,0];
% Jie = gE*Scale*(0.03)*[1;3;1];


% % try higher Jee for better ISN and SST reversal:
% % MSV =  C w_SV[(w_EE -d_E)(w_PP+d_P)-w_EP*w_PE]; where d=1/phi'
% % problem: w_EP*w_PE is much larger than w_PP*w_EE
% % other option: if PV rate goes down in the selective state (because E->SST->PV -> d_P goes up
% % and takes over dynamically! 
% gE=6; gI=9; 
% gEPV=15; % larger gEPV stabilizes the circuit but removes the SOM reversal
% gSST=2; 
% gVIP=2; % does not affect much the fixed points, larger values go against SST reversal, very small values reduce multistable interval
% Jee = gE*Scale*0.25;
% Jei = gI*Scale*(0.4).*[0.625,0.625,2.5].*[0.5*gEPV,gSST,0];
% Jii = gI*Scale*(0.15)*[0.625,0.625,2.5].*[gEPV,5*gSST,0;
%                     0,0,  0.2*gVIP ;
%                     1.,gSST,0];
% Jie = gE*Scale*(0.03)*[gEPV;gSST*5;gVIP];

% gPE controls the SST reversal: it's fixed to preserve w_EE*w_PP/w_EP*w_PE
% fixed! larger gPE leads to more ISN. But we need to keep PV decreasing
% its activity in the selective state by increasing a bit J_PS
% SST reversal works in the window beyond J_c were PV goes down as SST
% rises tracking E
% increasing J_VE increases VIP without affecting the rest
%----------
% these values have ISN in high and low state and SST reversal KEEP THEM
%----------
% g=4; gPE=2;
% Jee = g*2*Scale*0.045*gPE;
% Jei = g*3*Scale*(0.4).*[0.625,0.625,2.5].*[gPE,1,0];
% Jii = g*3*Scale*(0.15)*[0.625,0.625,2.5].*[gPE,1.,0;
%                     0,0,  0.4 ;
%                     1,1,0];
% Jie = g*2*Scale*(0.03)*[gPE;5;2.5];
g=4; gPE=2;
Jee = g*2*Scale*0.045*gPE;
Jei = g*3*Scale*(0.4).*[0.625,0.625,2.5].*[gPE,1,0];
Jii = g*3*Scale*(0.15)*[0.625,0.625,2.5].*[gPE,1.,0;
                    0,0,  0.4 ;
                    1,1,0];
Jie = g*2*Scale*(0.03)*[gPE;5;2.5];



    Jplus = 1.;
    %
%     theta_e=50.; %22.8385;%21.7124; 
%     theta_i=50.*[1,1,1] ;%theta_i=11.1897;%11.0887;  27.4304, (I) 58.6757,70.9415,14.0214,
    theta_e=6.00918;%27.4304; %22.8385;%21.7124; x, (I) 
    theta_i=[29.0197,46.4366,7.55815] ;%theta_i=11.1897;%11.0887;  
%     gam = 0.5;%1.125; % 0.5;  %1.3
    gam = 0.5;%1.125; % 0.5;  %1.3

    He = 0;%0.5;
    Hi = 0*[1,1,1];%0.5;

    Qpiu = 0.6;  %0.6
    Qmeno = 0.3;  %0.3
    rho = 2.75;    %(5)
    Jzero =11.5;  %1.

%     % connection probability (from Wang and Huang, 2017)
%     Cee = N_e*0.2; % 1600;         % 1600
%     Cie = N_e*0.5*[1;1;1]; % 1600;         % 1600
%     Cii = 0.5*N_i.*[1,.85,0;
%                     0,0,0.5;
%                     0.2,0.5,0]; % 400;         % 400
%     Cei = 0.2*N_i.*[1,1,0.]; % 400;         % 400
%     Cext = (N_e)*0.2; % 1600; 
        % connection probability
    pee_matrix = 0.2; % 1600;         % 1600
    pie_matrix = 0.5*[1;1;1]; % 1600;         % 1600
    pii_matrix = 0.5*[1,1,0;
                    0,0,1;
                    1,1,0]; % 400;         % 400
    pei_matrix = 0.2*[1,1,0.]; % 400;         % 400
    pext_matrix = 0.2; % 1600; 

    % external currents
%     Jie_ext=2*Scale*0.0915*[1,1,1];% normalized to ni_ext=2
%     Jee_ext=2*Scale*0.1027; % 
    Jie_ext=2*Scale*0.1*[1,1,1];% normalized to ni_ext=2
    Jee_ext=2*Scale*0.1; % 
    
    bgr=0.65; % N_e fraction of background neurons
    p = round(N_e*(1-bgr)/100); % keep 100 neurons per cluster
    f = (1-bgr)/p;       % 0.09
%     p_actual = 4;
%     nstories = 10;    % see monte.c
    nn = 5;                 % number of populations
    nfocus = 0;             % Number of populations in focus
    tausyn_e=0.002;
    tausyn_i=tausyn_e;

elseif  strcmp(Opt,'Small2_V1_HiRates')
    % network options
    Network.syn='SingleExp';
    Network.clust='hom';%'het'; % cluster: homogeneous='hom', heterogeneous='het'
    Network.clust_syn='';
    Network.clust_std=0.0; % het clusters: take a gaussian distribution with variance=mean*Network.clust_var
    N=800;
    N_e = N*4/5; % 16000
    N_i = (N/5)*[0.4,0.4,0.2]; % L2/3 in V1 from Kim et al (Osten lab) Cell 2017 
    Scale=(5000/N)^(1/2);
    % global spontaneous firing rates (neeed to fix thresholds)
    ni_e = 30; % 3.6;   %6.6   % AB97: 3.0
%     ni_i = 50*[1,0.5,0.5]; % 5.2;   %8.2   % AB97: 4.2
    ni_i = [50,30,20]; % 5.2;   %8.2   % AB97: 4.2
    % external currents
    ni_ext_e = 30;
    ni_ext_i = [50,30,20];

    tau_arp = .005;  % .002 (Unita': sec)
    tau_i = .02*[1,1,1];     % 0.01
    tau_e = .02;%15;     % 0.02

    delta = 0.0; % 0.;    
    x = 1.;
g=2; % these values are great
Jee = g*Scale*0.02;
Jei = g*Scale*(0.4).*[0.625,0.625,2.5].*[1,1,0];
Jii = g*Scale*(0.15)*[0.625,0.625,2.5].*[1,1,0;
                    0,0,  0.2 ;
                    0,1,0];
Jie = g*Scale*(0.03)*[1;2;1];

% g=8;
% Jee = g*Scale*0.06;
% Jei = g*Scale*(0.4).*[0.625,0.625,2.5].*[6,1,0];
% Jii = g*Scale*(0.15)*[0.625,0.625,2.5].*[6,1,0;
%                     0,0,  0.2 ;
%                     6,1,0];
% Jie = g*Scale*(0.05)*[2;2;1];


    
    Jplus = 1.;
    %
    theta_e=27.4304; %22.8385;%21.7124; 
    theta_i=[58.6757,70.9415,14.0214] ;%theta_i=11.1897;%11.0887;  
    Jplus = 5.;
    gam = 0.5;%1.125; % 0.5;  %1.3

    He = 0;%0.5;
    Hi = 0*[1,1,1];%0.5;

    Qpiu = 0.6;  %0.6
    Qmeno = 0.3;  %0.3
    rho = 2.75;    %(5)
    Jzero =11.5;  %1.

        % connection probability
    pee_matrix = 0.2; % 1600;         % 1600
    pie_matrix = 0.5*[1;1;1]; % 1600;         % 1600
    pii_matrix = 0.5*[1,1,0;
                    0,0,1;
                    1,1,0]; % 400;         % 400
    pei_matrix = 0.2*[1,1,0.]; % 400;         % 400
    pext_matrix = 0.2; % 1600; 

    % external currents
    Jie_ext=2*Scale*0.0915*[1,1,1];% normalized to ni_ext=2
    Jee_ext=2*Scale*0.1027; % 
    
    bgr=0.65; % N_e fraction of background neurons
    p = round(N_e*(1-bgr)/100); % keep 100 neurons per cluster
    f = (1-bgr)/p;       % 0.09
    nn = 5;                 % number of populations
    nfocus = 0;             % Number of populations in focus
    tausyn_e=0.002;
    tausyn_i=tausyn_e;

elseif  strcmp(Opt,'Small2_stats_V1_works')
    % network options
    Network.syn='SingleExp';
    Network.clust='hom';%'het'; % cluster: homogeneous='hom', heterogeneous='het'
    Network.clust_syn='';
    Network.clust_std=0.0; % het clusters: take a gaussian distribution with variance=mean*Network.clust_var
    N=800;
    N_e = N*4/5; % 16000
%     N_i = N/5*[0.5,0.25,0.25]; %4000
    N_i = N/5*[0.25,0.25,0.5]; % L2/3 in V1 from Kim et al (Osten lab) Cell 2017
    Scale=(5000/N)^(1/2);
    % global spontaneous firing rates (neeed to fix thresholds)
    ni_e = 5; % 3.6;   %6.6   % AB97: 3.0
    ni_i = 10*[1,0.5,0.5]; % 5.2;   %8.2   % AB97: 4.2
    % external currents
%     ni_ext = 7;
    ni_ext_e = 5;
    ni_ext_i = 10*[1,0,0];

    tau_arp = .005;  % .002 (Unita': sec)
    tau_i = .02*[1,1,1];     % 0.01
    tau_e = .02;%15;     % 0.02

    delta = 0.0; % 0.;    
    x = 1.;
    % e=pyr, i=PV,SST,VIP
    % take values from dipoppa, scanziani, 2020 biorxiv 
%     Jee = Scale*0.05; 
%     Jei = Scale*(0.06)*[1,1,0];    
%     Jii = Scale*(0.055)*[1,1,0;
%                          0,0,1 ;
%                          2,2,0];
%     Jie = Scale*(0.035/6)*[1;1;4];
% 

% these values are ok
%     Jee = Scale*0.02; 
%     Jei = Scale*(0.15)*[1.1,2,0];    
%     Jii = Scale*(0.13)*[1,2.2,0;
%                          0,0,2 ;
%                          3,2.1,0];
%     Jie = Scale*(0.035/2)*[2;2.1;6];
g=3;
Jee = g*Scale*0.03;
Jei = g*Scale*(0.4)*[1,1,0];
Jii = g*Scale*(0.15)*[1,1,0;
                    0,0,  0.2 ;
                    1,1,0];
Jie = g*Scale*(0.03)*[1;2;1];

%     % these weights below with Jplus=1 and give spontaneous rates: 13.2443  13.29  20.3668  2.50531  45.7198
%     Jee = Scale*0.05; 
%     Jei = Scale*(0.15)*[1,2,0];    
%     Jii = Scale*(0.15)*[1,2,0;
%                          0,0,2 ;
%                          2,4,0];
%     Jie = Scale*(0.035/2)*[2;2;6];
%     Jei=Jei.*(1+[-0.1351,    0.0876,   -0.0011]);
%     Jie=Jie.*(1+[0.0645;0.0937;0.0620]);
%     Jii=Jii.*(1+[0.2375,   -0.2134,    0.1139; 0.1308,    0.4544,    0.3877; -0.1031,    0.0742,   -0.3850]);
%     Jei=Jei+0.1*[-0.1351,    0.0876,   -0.0011];
%     Jie=Jie+0.1*[0.0645;0.0937;0.0620];
%     Jii=Jii+[0.1*0.2375,   -0.2134,    0.1139; 0.1308,    0.4544,    0.3877; -0.1031,    0.0742,   -0.3850];
    
    Jplus = 1.;
    %
%     theta_e=50.; %22.8385;%21.7124; 
%     theta_i=50.*[1,1,1] ;%theta_i=11.1897;%11.0887;  
    theta_e=32.0135; %22.8385;%21.7124; 
    theta_i=[60.1,70.6848,16.4979] ;%theta_i=11.1897;%11.0887;  
    Jplus = 5.;
%     gam = 0.5;%1.125; % 0.5;  %1.3
    gam = 0.5;%1.125; % 0.5;  %1.3

    He = 0;%0.5;
    Hi = 0*[1,1,1];%0.5;

    Qpiu = 0.6;  %0.6
    Qmeno = 0.3;  %0.3
    rho = 2.75;    %(5)
    Jzero =11.5;  %1.

%     % connection probability (from Wang and Huang, 2017)
%     Cee = N_e*0.2; % 1600;         % 1600
%     Cie = N_e*0.5*[1;1;1]; % 1600;         % 1600
%     Cii = 0.5*N_i.*[1,.85,0;
%                     0,0,0.5;
%                     0.2,0.5,0]; % 400;         % 400
%     Cei = 0.2*N_i.*[1,1,0.]; % 400;         % 400
%     Cext = (N_e)*0.2; % 1600; 
        % connection probability
    pee_matrix = 0.2; % 1600;         % 1600
    pie_matrix = 0.5*[1;1;1]; % 1600;         % 1600
    pii_matrix = 0.5*[1,1,0;
                    0,0,1;
                    1,1,0]; % 400;         % 400
    pei_matrix = 0.2*[1,1,0.]; % 400;         % 400
    pext_matrix = 0.2; % 1600; 

    % external currents
    Jie_ext=2*Scale*0.0915*[1,0,0.];% normalized to ni_ext=2
    Jee_ext=2*Scale*0.1027; % 
    
    bgr=0.65; % N_e fraction of background neurons
    p = round(N_e*(1-bgr)/100); % keep 100 neurons per cluster
    f = (1-bgr)/p;       % 0.09
%     p_actual = 4;
%     nstories = 10;    % see monte.c
    nn = 5;                 % number of populations
    nfocus = 0;             % Number of populations in focus
    tausyn_e=0.002;
    tausyn_i=tausyn_e;

elseif  strcmp(Opt,'Small2_stats_V1_works_lowSST')
    % network options
    Network.syn='SingleExp';
    Network.clust='hom';%'het'; % cluster: homogeneous='hom', heterogeneous='het'
    Network.clust_syn='';
    Network.clust_std=0.0; % het clusters: take a gaussian distribution with variance=mean*Network.clust_var
    N=800;
    N_e = N*4/5; % 16000
%     N_i = N/5*[0.5,0.25,0.25]; %4000
    N_i = N/5*[0.25,0.25,0.5]; % L2/3 in V1 from Kim et al (Osten lab) Cell 2017
    Scale=(5000/N)^(1/2);
    % global spontaneous firing rates (neeed to fix thresholds)
    ni_e = 5; % 3.6;   %6.6   % AB97: 3.0
    ni_i = 10*[1,0.5,0.5]; % 5.2;   %8.2   % AB97: 4.2
    % external currents
%     ni_ext = 7;
    ni_ext_e = 5;
    ni_ext_i = 10*[1,0,0];

    tau_arp = .005;  % .002 (Unita': sec)
    tau_i = .02*[1,1,1];     % 0.01
    tau_e = .02;%15;     % 0.02

    delta = 0.0; % 0.;    
    x = 1.;
    % e=pyr, i=PV,SST,VIP
    % take values from dipoppa, scanziani, 2020 biorxiv 

    g=3;
    Jee = g*Scale*0.03;
    Jei = g*Scale*(0.4)*[1.2,0.3,0];
    Jii = g*Scale*(0.15)*[1,1,0;
                        0,0,  0.2 ;
                        1,1,0];
    Jie = g*Scale*(0.03)*[1.2;1.5;1];

    Jplus = 1.;
    %
    theta_e=32.0135; %22.8385;%21.7124; 
    theta_i=[60.1,70.6848,16.4979] ;%theta_i=11.1897;%11.0887;  
    Jplus = 5.;
%     gam = 0.5;%1.125; % 0.5;  %1.3
    gam = 0.5;%1.125; % 0.5;  %1.3

    He = 0;%0.5;
    Hi = 0*[1,1,1];%0.5;

    Qpiu = 0.6;  %0.6
    Qmeno = 0.3;  %0.3
    rho = 2.75;    %(5)
    Jzero =11.5;  %1.

    pee_matrix = 0.2; % 1600;         % 1600
    pie_matrix = 0.5*[1;1;1]; % 1600;         % 1600
    pii_matrix = 0.5*[1,1,0;
                    0,0,1;
                    1,1,0]; % 400;         % 400
    pei_matrix = 0.2*[1,1,0.]; % 400;         % 400
    pext_matrix = 0.2; % 1600; 

    % external currents
    Jie_ext=2*Scale*0.0915*[1,0,0.];% normalized to ni_ext=2
    Jee_ext=2*Scale*0.1027; % 
    
    bgr=0.65; % N_e fraction of background neurons
    p = round(N_e*(1-bgr)/100); % keep 100 neurons per cluster
    f = (1-bgr)/p;       % 0.09
%     p_actual = 4;
%     nstories = 10;    % see monte.c
    nn = 5;                 % number of populations
    nfocus = 0;             % Number of populations in focus
    tausyn_e=0.002;
    tausyn_i=tausyn_e;    
elseif strcmp(Opt,'Multifinal')
    % network options
    Network.syn='SingleExp';
    Network.clust='het'; % cluster: homogeneous='hom', heterogeneous='het'
    Network.clust_syn='';
    Network.clust_std=0.01; % het clusters: take a gaussian distribution with variance=mean*Network.clust_var
    N=2000;
    N_e = N*4/5; % 16000
    N_i = N/5; %4000
    Scale=(5000/N)^(1/2);
    % global spontaneous firing rates (neeed to fix thresholds)
    ni_e = 5; % 3.6;   %6.6   % AB97: 3.0
    ni_i = 7; % 5.2;   %8.2   % AB97: 4.2
%     ni_ext = 7;
    ni_ext_e = 7;
    ni_ext_i = 7;

    tau_arp = .005;  % .002 (Unita': sec)
    tau_i = .02;     % 0.01
    tau_e = .02;%15;     % 0.02

    delta = 0.01; % SD of synaptic weights: Jee*(1+delta*randn(N_e))
    x = 1.;         % SE X!=1: Frazione di connessioni eccit. interne al modulo

    Jee = Scale*0.015; 
    Jii = (5/2)^(1/2)*Scale*0.06;
    Jie = Scale*0.02;
    Jei = (5/2)^(1/2)*Scale*0.045;    
    theta_e=66.5052; 
    theta_i=32.0911;  
%     theta_e=66.5052; 
%     theta_i=32.0911;  
    Jplus = 10;
    Jminus = 1.0;
    gam = 0.5;%1.125; % 0.5;  %1.3

    He = 0;%0.5;
    Hi = 0;%0.5;

    Qpiu = 0.6;  %0.6
    Qmeno = 0.3;  %0.3
    rho = 2.75;    %(5)
    Jzero =5;  %1.

    pee_matrix = N_e*0.2; % 1600;         % 1600
    pie_matrix = N_e*0.5; % 1600;         % 1600
    pii_matrix = N_i*0.5; % 400;         % 400
    pei_matrix = N_i*0.5; % 400;         % 400
    pext_matrix = (N_e)*0.2; % 1600; 
%     Jie_ext=Scale*1.025/(Cext*ni_ext*tau_i);% normalized to ni_ext=2
%     Jee_ext=Scale*2*1.15/(Cext*ni_ext*tau_e); % 
    Jie_ext=0.8*Scale*0.0915;% normalized to ni_ext=2
    Jee_ext=0.8*Scale*0.1027; % 
    
    bgr=0.1; % fraction of background neurons
%     p = round(30/Scale);            % 3
    p = round(N_e*(1-bgr)/100); % keep 100 neurons per cluster
    f = (1-bgr)/p;       % 0.09
    p_actual = 4;
    nstories = 10;    % see monte.c
    nn = 4;                 % number of populations
    nfocus = 0;             % Number of populations in focuselse
    tausyn_e=0.004;
    tausyn_i=0.004;
end% STIMULI

Jminus = 1-gam*f*(Jplus-1);
if Jminus<0
    error('ERROR in create_params: Jminus<0. Adjust gam or f\n');
end
tau_1=0.001;
%------------------------------------------------------------------------
% VARIABLES & PARAMETERS: ------------------------------------------------
%------------------------------------------------------------------------/
% PLOTS
Sim.Plotf=0;
Sim.ne_raster=20;
Sim.ni_raster=0;
Sim.plot_length=Sim.t_End-Sim.t_Start;
% indices of ensemble units to store
exc=randperm(N_e);
inh=N_e+1:N_e+N_i;
inh=inh(randperm(numel(inh)));
Sim.ind_p=[exc(1:Sim.ne_raster) inh(1:Sim.ni_raster)]; % ensemble
Sim.weights_save='off';
%
% INITIAL CONDITIONS
%     initial='fixed'; % fix same initial conditions for all trials and events
initial_ext='fixed'; % same external current for all trials and events

%
Stimulus.feat=feat;
Network.stim=Stimulus;

% default directories where to save
sim_folder=fullfile('DATA','SIM',sprintf('%s',Opt)); % ADD PARAMETERS
data_folder=fullfile('DATA','Results',sprintf('%s',Opt));
if ~exist(data_folder,'dir')
    mkdir(data_folder);
end
if ~exist(sim_folder,'dir')
    mkdir(sim_folder);
end 
Sim.save_data=data_folder;
Sim.save_sim=sim_folder;
%
extra=''; % extra string to append on filename

save(filename,'ni_e','ni_i','ni_ext_e','ni_ext_i','tau_arp','tau_i','tau_e','theta_e',...
    'theta_i','delta','f','x','Jee','Jii','Jie','Jei','Jee_ext','Jie_ext',...
    'Jplus','Jminus','gam','He','Hi','N','N_e',...
    'N_i','pee_matrix','pie_matrix','pei_matrix','pii_matrix','pext_matrix','p','nn',...
    'nfocus','Opt','Sim','Network','Stimulus','tausyn_e','tausyn_i','tau_1','initial_ext','extra');