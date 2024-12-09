function create_params(paramsfile,Opt)
%----------
% NETWORKS
%----------
% Start and end of trial (in units of seconds)
Sim.t_Start=-0.2;
Sim.t_End=0.2;
Sim.dt_step=0.001; % integration step (s)

%----------------
% DEFAULT STIMULI
%----------------
% STIMULUS
Stimulus.option='on'; % 'off'
Stimulus.input='Const'; % constant external current
%------------------
% TIME CONSTANTS
%------------------
tau_arp = .005;  % refractory period
tau_i = .02;     % inh membrane time
tau_e = .02;	 % exc membrane time
tausyn_e=0.005;  % exc synaptic time
tausyn_i=0.005;  % inh synaptic time

if strcmp(Opt,'EI')
    %------------------
    % NETWORK OPTIONS
    %------------------
    Network.syn='SingleExp';
    Network.clusters={'EE','EI','IE','II'}; % all weights are clustered
    % Network.clust='hom'; % homogeneous EE cluster size
    Network.clust='het'; % heterogeneous EE cluster size
    Network.clust_std=0.2; % het cluster size: cluster size is sampled from a gaussian with SD=Ncluster*Network.std
    Network.clustEI='hom'; % EI clusters: homogeneous='hom', heterogeneous='het'
    Network.clustIE='hom'; % IE clusters: homogeneous='hom', heterogeneous='het'
    Network.clustII='hom'; % II clusters: homogeneous='hom', heterogeneous='het'
    % Network.clust_syn='';
    N=2000; % network size
    N_e = N*4/5; % exc neurons
    N_i = N/5; % inh neurons
    Scale=(1000/N)^(1/2);
    %------------------
    % SYNAPTIC WEIGHTS
    %------------------
    % syn weights are drawn from a gaussian distribution with std delta and
    % mean Jab
    Jee = Scale*0.02;
    Jii = 1.*Scale*0.12; %     Jii = Scale*0.06;
    Jie = 2*Scale*0.010; %     Jei 3 Scale*0.045;
    Jei = 3.*Scale*0.02;
    %------------------
    % CLUSTER PARAMETERS
    %------------------
    delta=0.2;% SD of synaptic weights distribution, eg: ~Jee*(1+delta*randn(N_e))
    Network.delta = delta; % SD of synaptic weights: Jee*(1+delta*randn(N_e)) Larger delta -> clusters look more async in the high state
    Network.deltaEI = delta; % SD of synaptic weights: Jee*(1+delta*randn(N_e)) Larger delta -> clusters look more async in the high state
    Network.deltaIE = delta; % SD of synaptic weights: Jee*(1+delta*randn(N_e)) Larger delta -> clusters look more async in the high state
    Network.deltaII = delta; % SD of synaptic weights: Jee*(1+delta*randn(N_e)) Larger delta -> clusters look more async in the high state
    Jplus = 14; % EE intra-cluster potentiation factor
    Network.factorEI = 10; % EI intra-cluster potentiation factor
    Network.factorIE = 8; % IE intra-cluster potentiation factor
    Network.factorII = 5; % II intra-cluster potentiation factor
    %
    bgr=0.1; % fraction of background excitatory neurons (unclustered)
    Network.bgrI=0.1; % fraction of background neurons
    Ncluster=80; % average # neurons per cluster
    p = round(N_e*(1-bgr)/Ncluster); % # of clusters
    f = (1-bgr)/p; % fraction of E neurons per cluster
    Network.fI = (1-Network.bgrI)/p;       % fraction of I neurons per cluster
    %------------------
    % THRESHOLDS
    %------------------
    theta_e=1.42824;% exc threshold potential
    theta_i=0.74342;% inh threshold potentials
    % reset potentials
    He = 0;%
    Hi = 0;%
    % EXTERNAL CURRENT
    % default external currents
    ni_ext_e = 7; % 7;
    ni_ext_i = 7; % 7;
%     N_ext=N_e;
    Jie_ext=0.8*Scale*0.0915;% external input synaptic strengths
    Jee_ext=0.8*Scale*0.1027; %
    %------------------
    % CONNECTIVITY PARAMETERS
    %------------------
    p_ee=0.2;
    p_ie=0.5;
    p_ii=0.5;
    p_ei=0.5;
    p_ext = 0.2; % 1600;
    nn = 4;                 % number of populations
    %                 gam=1/(2-f*(p+1));%
elseif strcmp(Opt,'E')
    % Creates parameter for network with E but not I clusters
    %------------------
    % NETWORK OPTIONS
    %------------------
    % Network.clust='hom'; % homogeneous cluster size
    Network.syn='SingleExp';
    Network.clust='hom'; % heterogeneous cluster size
    Network.clust_std=0.01; % het cluster size: cluster size is sampled from a gaussian with SD=Ncluster*Network.std
    N=2000; % network size
    N_e = N*4/5; % E neurons
    N_i = N/5; % I neurons
    Scale=(5000/N)^(1/2);
    ni_e=5; ni_i=7;
    %------------------
    % SYNAPTIC WEIGHTS
    %------------------
    % syn weights are drawn from a gaussian distribution with std delta and
    % mean Jab
    Jee = Scale*0.015; % mean Jee weight
    Jii = (5/2)^(1/2)*Scale*0.06;
    Jie = Scale*0.02;
    Jei = (5/2)^(1/2)*Scale*0.045;
    delta = 0.1; % SD of synaptic weights distribution, eg: ~Jee*(1+delta*randn(N_e))
    Network.delta = delta; % SD of synaptic weights: Jee*(1+delta*randn(N_e)) Larger delta -> clusters look more async in the high state
    Network.deltaEI = delta; % SD of synaptic weights: Jee*(1+delta*randn(N_e)) Larger delta -> clusters look more async in the high state
    Network.deltaIE = delta; % SD of synaptic weights: Jee*(1+delta*randn(N_e)) Larger delta -> clusters look more async in the high state
    Network.deltaII = delta; % SD of synaptic weights: Jee*(1+delta*randn(N_e)) Larger delta -> clusters look more async in the high state
    %------------------
    % CLUSTER PARAMETERS
    %------------------
    Jplus = 9; % intra-cluster potentiation factor for E clusters
    bgr=0.1; % fraction of background excitatory neurons (unclustered)
    Ncluster=100; % average # neurons per cluster
    p = round(N_e*(1-bgr)/Ncluster); % # of clusters
    f = (1-bgr)/p;
    %------------------
    % THRESHOLDS
    %------------------
    theta_e=4; % exc threshold potential
    theta_i=4; % inh threshold potentials
    % reset potentials
    He = 0;%
    Hi = 0;%
    % EXTERNAL CURRENT
    % default external currents
    ni_ext_e = 7; % 7;
    ni_ext_i = 7; % 7;
%     N_ext=N_e;
    Jie_ext=0.8*Scale*0.0915;% external input synaptic strengths
    Jee_ext=0.8*Scale*0.1027; %
    %------------------
    % CONNECTIVITY PARAMETERS
    %------------------
    p_ee=0.2;
    p_ie=0.5;
    p_ii=0.5;
    p_ei=0.5;
    p_ext = 0.2; % 1600;
    nn = 4;                 % number of populations
    %                 gam=1/(2-f*(p+1));%
elseif strcmp(Opt,'Small2_stats_cosyne')
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
    N_ext=N_e;
    Jie_ext=2*Scale*0.0915*[1,1,1];% normalized to ni_ext=2
    Jee_ext=2*Scale*0.1027; %
    nn = 5;                 % number of populations
elseif strcmp(Opt,'InhTypes')
    % Creates parameter for network with E but not I clusters
    %------------------
    % NETWORK OPTIONS
    %------------------
    % Network.clust='hom'; % homogeneous cluster size
    Network.syn='SingleExp';
    Network.clust='hom'; % heterogeneous cluster size
    Network.clust_std=0.0; % het cluster size: cluster size is sampled from a gaussian with SD=Ncluster*Network.std
    N=1000; % network size
    N_e = N*4/5; % E neurons
    N_i = (N/5)*[0.4,0.4,0.2]; % L2/3 in V1 from Kim et al (Osten lab) Cell 2017
    ni_e = 5; %
    ni_i = 10*[1.5,1,1]; %
    %--------------------
    % NEURON PARAMETERS
    %--------------------
    tau_i = .02*[1,1,1];     % 0.01
    tausyn_e=0.002;
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
    g=2.; gPE=2.5;
    Jee = g*0.4*Scale*0.045*gPE;
    Jei = g*3*Scale*(0.3).*[0.625,0.625,2.5].*[gPE,1,0];
    Jii = g*3*Scale*(0.15)*[0.625,0.625,2.5].*[2*gPE,1.,0;
        0,0,  0.4 ;
        1,1,0];
    Jie = g*3*Scale*(0.02)*[gPE;5 ;0.5];
    delta = 0.0; % SD of synaptic weights distribution, eg: ~Jee*(1+delta*randn(N_e))
    Network.delta = delta; % SD of synaptic weights: Jee*(1+delta*randn(N_e)) Larger delta -> clusters look more async in the high state
    Network.deltaEI = delta; % SD of synaptic weights: Jee*(1+delta*randn(N_e)) Larger delta -> clusters look more async in the high state
    Network.deltaIE = delta; % SD of synaptic weights: Jee*(1+delta*randn(N_e)) Larger delta -> clusters look more async in the high state
    Network.deltaII = delta; % SD of synaptic weights: Jee*(1+delta*randn(N_e)) Larger delta -> clusters look more async in the high state
    %------------------
    % CLUSTER PARAMETERS
    %------------------
    %                 Jplus =10.4; % intra-cluster potentiation factor for E clusters
    %                 Jplus =10.4; % intra-cluster potentiation factor for E clusters
    Jplus =10.4; % intra-cluster potentiation factor for E clusters
    bgr=0.5; % fraction of background excitatory neurons (unclustered)
    Ncluster=100; % average # neurons per cluster
    p = round(N_e*(1-bgr)/Ncluster); % # of clusters
    f = (1-bgr)/p;
    %------------------
    % THRESHOLDS
    %------------------
    theta_e=6.00918; % exc threshold potential
    theta_i=[29.0197,46.4366,7.55815]; % inh threshold potentials
    % reset potentials
    He = 0;%
    Hi = zeros(1,3);%
    %------------------
    % EXTERNAL CURRENT
    %------------------
    % default external currents
    ni_ext = 7; % 7;
    ni_ext_e = 5;
    ni_ext_i = 10*[2,1,1];
    %------------------
    % CONNECTIVITY PARAMETERS
    %------------------
    p_ee = 0.2; % 1600;         % 1600
    p_ie = 0.5*[1;1;1]; % 1600;         % 1600
    p_ii = 0.5*[1,1,0;
        0,0,1;
        1,1,0]; % 400;         % 400
    p_ei = 0.2*[1,1,0.]; % 400;         % 400
    p_ext = 0.2; % 1600;
    %------------------
    % EXTERNAL BIAS
    %------------------
    % - low bias: oscillating; high bias: asynchronous
    % higher J_ext -> higher Jplus to get same dynamics
    N_ext=N_e;
    Jie_ext=1.5.*Scale*[0.5,0.5,0.5];% normalized to ni_ext=2
    Jee_ext=1.2*Scale; %
    nn = 5;                 % number of populations
elseif strcmp(Opt,'Marcel0')
    Stimulus.input='Poisson'; % constant external current
    % Creates parameter for network with E but not I clusters
    %------------------
    % NETWORK OPTIONS
    %------------------
    Network.syn='SingleExp';
    % Network.clust='hom'; % homogeneous cluster size
    Network.clust='hom'; % heterogeneous cluster size
    Network.clust_std=0.0; % het cluster size: cluster size is sampled from a gaussian with SD=Ncluster*Network.std
    N=5000; % network size
    N_e = N*4/5; % E neurons
    N_i = (N/5)*[0.5,0.25,0.25]; 
    ni_e = 2.5; %
    ni_i = [12.5,2,0.2]; %
    %--------------------
    % NEURON PARAMETERS
    %--------------------
    tau_arp = .001;  % refractory period
    tau_i = [0.02,0.02,0.02];     % 0.01
    tausyn_e=0.003;
    tausyn_i=[0.004,0.005,0.005];     %*ones(1,numel(N_i));
    
%     Scale=(1/N)^(1/2);
    wEE=0.35;
    %------------------
    % SYNAPTIC WEIGHTS
    %------------------
    % syn weights are drawn from a gaussian distribution with std delta and
    Jee = wEE;
    Jei = wEE*[12,4,0];
    Jii = wEE.*[9,4,0;
                0,0,1.5;
                0.2,6,0];
    Jie = wEE*[4;0.2;0.2];
    delta = 0.0; % SD of synaptic weights distribution, eg: ~Jee*(1+delta*randn(N_e))
    Network.delta = delta; % SD of synaptic weights: Jee*(1+delta*randn(N_e)) Larger delta -> clusters look more async in the high state
    Network.deltaEI = delta; % SD of synaptic weights: Jee*(1+delta*randn(N_e)) Larger delta -> clusters look more async in the high state
    Network.deltaIE = delta; % SD of synaptic weights: Jee*(1+delta*randn(N_e)) Larger delta -> clusters look more async in the high state
    Network.deltaII = delta; % SD of synaptic weights: Jee*(1+delta*randn(N_e)) Larger delta -> clusters look more async in the high state
    %------------------
    % CLUSTER PARAMETERS
    %------------------
    Jplus =1; % intra-cluster potentiation factor for E clusters
    bgr=0.25; % fraction of background excitatory neurons (unclustered)
    Ncluster=200; % average # neurons per cluster
    p = round(N_e*(1-bgr)/Ncluster); % # of clusters
    f = (1-bgr)/p;
    %------------------
    % THRESHOLDS
    %------------------
    Network.ThresholdOpt=0; % use thresholds below
    theta_e=10; % exc threshold potential
    theta_i=10*[1,1,1]; % inh threshold potentials
    % reset potentials
    He = 0;%
    Hi = zeros(1,3);%
    %------------------
    % EXTERNAL CURRENT
    %------------------
    % default external currents
    ni_ext = 5; % 7;
    ni_ext_e = 5;
    ni_ext_i = 5*[1,1,1];
    %------------------
    % CONNECTIVITY PARAMETERS
    %------------------
    p_ee = 0.1; % 1600;         % 1600
    p_ie = 0.1*[1;1;1]; % 1600;         % 1600
    p_ii = 0.1*[1,1,0;
        0,0,1;
        1,1,0]; % 400;         % 400
    p_ei = 0.1*[1,1,0.]; % 400;         % 400
    p_ext = 1; % 1600;
    %------------------
    % EXTERNAL BIAS
    %------------------
    % - low bias: oscillating; high bias: asynchronous
    % higher J_ext -> higher Jplus to get same dynamics
    N_ext=1000;
    Jie_ext=wEE*[ 0.7, 0.25, 0.27];% normalized to ni_ext=2
    Jee_ext=wEE; %
    nn = 5;                 % number of populations
    %                 gam=0.5;%
elseif ~isempty(strfind(Opt,'MarcelPost'))
%     strcmp(Opt,'MarcelPost_disinh')
%     Stimulus.input='Poisson'; % constant external current
    % Creates parameter for network with E but not I clusters
    %------------------
    % NETWORK OPTIONS
    %------------------
    Network.syn='SingleExp';
    % Network.clust='hom'; % homogeneous cluster size
    Network.clust='hom'; % heterogeneous cluster size
    Network.clust_std=0.0; % het cluster size: cluster size is sampled from a gaussian with SD=Ncluster*Network.std
    N=5000; % network size
    N_e = N*4/5; % E neurons
    N_i = (N/5)*[0.5,0.25,0.25]; 
    ni_e = 2.5; %
    ni_i = [12.5,2,0.2]; %
    %--------------------
    % NEURON PARAMETERS
    %--------------------
    tau_arp = .001;  % refractory period
    tau_i = [0.02,0.02,0.02];     % 0.01
    tausyn_e=0.003;
    tausyn_i=[0.004,0.005,0.005];     %*ones(1,numel(N_i));
    
%     Scale=(1/N)^(1/2);
    wEE=0.35;
    %------------------
    % SYNAPTIC WEIGHTS
    %------------------
    % syn weights are drawn from a gaussian distribution with std delta and
    Jee = wEE;
    Jei = wEE*[12,4,0];
    Jii = wEE.*[9,4,0;
                0,0,1.5;
                0.2,6,0];
    Jie = wEE*[4;0.2;0.2];
    delta = 0.0; % SD of synaptic weights distribution, eg: ~Jee*(1+delta*randn(N_e))
    Network.delta = delta; % SD of synaptic weights: Jee*(1+delta*randn(N_e)) Larger delta -> clusters look more async in the high state
    Network.deltaEI = delta; % SD of synaptic weights: Jee*(1+delta*randn(N_e)) Larger delta -> clusters look more async in the high state
    Network.deltaIE = delta; % SD of synaptic weights: Jee*(1+delta*randn(N_e)) Larger delta -> clusters look more async in the high state
    Network.deltaII = delta; % SD of synaptic weights: Jee*(1+delta*randn(N_e)) Larger delta -> clusters look more async in the high state
    %------------------
    % CLUSTER PARAMETERS
    %------------------
    Jplus =1; % intra-cluster potentiation factor for E clusters
    bgr=0.25; % fraction of background excitatory neurons (unclustered)
    Ncluster=200; % average # neurons per cluster
    p = round(N_e*(1-bgr)/Ncluster); % # of clusters
    f = (1-bgr)/p;
    %------------------
    % THRESHOLDS
    %------------------
    Network.ThresholdOpt=0; % use thresholds below
    theta_e=10; % exc threshold potential
    theta_i=10*[1,1,1]; % inh threshold potentials
    % reset potentials
    He = 0;%
    Hi = zeros(1,3);%
    %------------------
    % EXTERNAL CURRENT
    %------------------
    % default external currents
    ni_ext = 5; % 7;
    ni_ext_e = 5;
    ni_ext_i = 5*[1,1,1];
    %------------------
    % CONNECTIVITY PARAMETERS
    %------------------
    p_ee = 0.1; % 1600;         % 1600
    p_ie = 0.1*[1;1;1]; % 1600;         % 1600
    p_ii = 0.1*[1,1,0;
        0,0,1;
        1,1,0]; % 400;         % 400
    p_ei = 0.1*[1,1,0.]; % 400;         % 400
    p_ext = 1; % 1600;
    %------------------
    % EXTERNAL BIAS
    %------------------
    % - low bias: oscillating; high bias: asynchronous
    % higher J_ext -> higher Jplus to get same dynamics
    N_ext=1000;
    Jie_ext=wEE*[ 0.7, 0.25, 0.27];% normalized to ni_ext=2
    Jee_ext=wEE; %
    nn = 5;                 % number of populations
    %                 gam=0.5;%
elseif strcmp(Opt,'Marcel2')
    % Creates parameter for network with E but not I clusters
    %------------------
    % NETWORK OPTIONS
    %------------------
    Network.syn='SingleExp';
    % Network.clust='hom'; % homogeneous cluster size
    Network.clust='hom'; % heterogeneous cluster size
    Network.clust_std=0.0; % het cluster size: cluster size is sampled from a gaussian with SD=Ncluster*Network.std
    N=2000; % network size
    N_e = N*4/5; % E neurons
    N_i = (N/5)*[0.5,0.25,0.25]; % L2/3 in V1 from Kim et al (Osten lab) Cell 2017
    ni_e = 7; %
    ni_i = 7*[1,1,1]; %
    %--------------------
    % NEURON PARAMETERS
    %--------------------
    tau_i = .02*[1,1,1];     % 0.01
    tausyn_e=0.004;
    tausyn_i=tausyn_e*ones(1,numel(N_i));
    
    Scale=10*(1/N)^(1/2);
    %------------------
    % SYNAPTIC WEIGHTS
    %------------------
    % syn weights are drawn from a gaussian distribution with std delta and
    Jee = 1.1*Scale;
    Jei = Scale*[5.0,5.0,0];
    Jii = Scale.*[6.7,1.5,0;
        0,0,  0.55 ;
        1.3,3,0];
    Jie = 1.4*Scale*[1;1;1];
    
    % jee0=1.1/sqrt(N) # in mV
    % 	jse=1.4/sqrt(N)
    % 	jve=1.4/sqrt(N)
    % 	jpe=1.4/sqrt(N)
    %
    % 	jep0=5.0/sqrt(N)
    % 	jpp=6.7/sqrt(N)
    % 	jvp=1.3/sqrt(N)
    %
    % 	jes0=5.0/sqrt(N)
    % 	jps0=1.5/sqrt(N)
    % 	jvs=3.0/sqrt(N)
    % 	jsv=0.55/sqrt(N)
    delta = 0.0; % SD of synaptic weights distribution, eg: ~Jee*(1+delta*randn(N_e))
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
    bgr=0.1; % fraction of background excitatory neurons (unclustered)
    Ncluster=100; % average # neurons per cluster
    p = round(N_e*(1-bgr)/Ncluster); % # of clusters
    f = (1-bgr)/p;
    %------------------
    % THRESHOLDS
    %------------------
    theta_e=3.9; % exc threshold potential
    theta_i=4*[1,1,1]; % inh threshold potentials
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
%     p_ee = 0.1; % 1600;         % 1600
%     p_ie = 0.1*[1;1;1]; % 1600;         % 1600
%     p_ii = 0.1*[1,1,0;
%         0,0,1;
%         1,1,0]; % 400;         % 400
%     p_ei = 0.1*[1,1,0.]; % 400;         % 400
%     p_ext = 0.2; % 1600;
    p_ee = 0.2; % 1600;         % 1600
    p_ie = 0.5*[1;1;1]; % 1600;         % 1600
    p_ii = 0.5*[1,1,0;
        0,0,1;
        1,1,0]; % 400;         % 400
    p_ei = 0.5*[1,1,0.]; % 400;         % 400
    p_ext = 0.2; % 1600;
    %------------------
    % EXTERNAL BIAS
    %------------------
    % - low bias: oscillating; high bias: asynchronous
    % higher J_ext -> higher Jplus to get same dynamics
    N_ext=N*4/5;
    Jie_ext=1000*Scale*[ 0.0052, 0.0020, 0.0028];% normalized to ni_ext=2
    Jee_ext=1000*0.0058*Scale; %
    nn = 5;                 % number of populations
    %                 gam=0.5;%
end
%------------------
% EXTERNAL BIAS
%------------------
% external input parameters, eg: external current given by mu_e_ext=Cext*Jee_ext*ni_ext
mu_E0=(N_ext).*p_ext.*Jee_ext.*ni_ext_e;
mu_I0=(N_ext).*p_ext.*Jie_ext.*ni_ext_i;
%  ext current
fprintf('External baseline: I_{E0}: %0.03gmV ',mu_E0);
for i=1:numel(N_i);  fprintf(' I_{I%d}: %0.03gmV',i,mu_I0(i)); end; fprintf('\n');
Mu=[mu_E0*(ones(N_e,1)+(0./2)*(2*rand([N_e,1])-1))];
for i=1:numel(N_i)
    Mu=[Mu; mu_I0(i)*(ones(N_i(i),1)+(0./2)*(2*rand([N_i(i),1])-1))];     % bias
end
%             gam=1/(2-f*(p+1));%
gam=0.5;%
x = 1.; % for MFT
tau_1=0.001;
Jminus = 1.-gam*f*(Jplus-1.); %Jminus=(2-f*Q-f*Jplus)/(2-f*(Q+1))=[2-f(Q+Jplus)]/(2-f*(Q+1))
%------------------
% CONNECTIVITY PARAMETERS
%------------------
scnt=0;
% TASTE (specific stimulus)
scnt=scnt+1;
feat(scnt).name='US'; % unconditioned stimulus (taste)
feat(scnt).interval=[0 Sim.t_End]; % stimulus interval
gain=0.1; % stimulus value at 1 s
feat(scnt).gain=gain;
feat(scnt).profile=@(t)t; % time course of stimulus, eg a linear ramp
feat(scnt).selectivity='mixed'; % random half of clusters are selective for each stimulus
feat(scnt).selective=rand(1,p)<0.5; % US selective clusters
feat(scnt).connectivity=0.5; % fraction of selective neurons within a selective cluster
% ANTICIPATORY CUE
scnt=scnt+1;
feat(scnt).name='CSgauss'; % conditioned stimulus (cue) with "quenched" noise
feat(scnt).interval=[min(-0.5,Sim.t_Start), Sim.t_End]; % stimulus interval
gain=0.1; % SD of quenched noise across neurons
feat(scnt).gain=gain;
%             tau_cue=[0.5,1]; % rise and decay time of double exp cue time course
%             feat(scnt).profile=@(t)(1/(tau_cue(2)-tau_cue(1)))*(exp(-t/tau_cue(2))-exp(-t/tau_cue(1))); % double exp profile time course
feat(scnt).profile=@(t)1;
feat(scnt).selectivity='exc'; % targets exc neurons only
feat(scnt).selective=ones(1,p); % CS targets all clusters
feat(scnt).connectivity=0.50; % fraction of neurons targeted within each cluster
%
Stimulus.feat=feat;

%------------------------------------------------------------------------
% PLOT PARAMETERS: ------------------------------------------------
%------------------------------------------------------------------------/
% PLOTS
Sim.Plotf=0;
Sim.plot_length=Sim.t_End-Sim.t_Start; % length of plot intervals
% indices of ensemble units to store
exc=randperm(N_e);
inh=N_e+randperm(N_i(1));
Sim.ind_p=[exc(1) inh(1)]; % choosing neuron index for membrane potential plot (one E and one I)
Sim.weights_save='off'; % save weight matrix: 'Yes'
extra='';
save(paramsfile,'ni_e','ni_i','ni_ext_i','ni_ext_e','tau_arp','tau_i','tau_e','theta_e',...
    'theta_i','delta','f','Jee','Jii','Jie','Jei','Jee_ext','Jie_ext',...
    'Jplus','Jminus','gam','He','Hi','N_e','nn',...
    'N_i','N_ext','p_ee','p_ie','p_ei','p_ii','p_ext','p','Sim','Network','Stimulus',...
    'tausyn_e','tausyn_i','extra','Mu','paramsfile');
fprintf('Network parameters saved in %s\n',paramsfile);
end