% Master script for MFT with Inh cell types
% Luca Mazzucato 2020
clear all

%---------------------------------------------
% options for running script
% first time: run RUNMFT=1, RUN=1, PLOTMFT=1.
%---------------------------------------------
RUNMFT=1; % Run MFT with scan on J+, find all fixed points, check stability.
RUN=1; % Run effective MFT to calculate potential wells
PLOTMFT=1; % Plot MFT bifurcation diagrams, potential wells, gain modulation

%-----------------------
% choose network
%-----------------------
% these are parameters for a network 'Small2InhDens' WITHOUT top-down reversal
% Opt='Small2InhDens';%'Small2_stats_V1';%
% Jplus=2.7; Jstop=3.5; Jzero=2; JplusLow=2.2;
% these are parameter for a network 'Small2Reversal' WITH top-down reversal
Opt='Small2Reversal';% LARGE REVERSAL
% Opt='Small2_stats_V1';% NO REVERSAL
% Opt='Small2_stats_V1_works';% SMALL REVERSAL
% Opt='Small2_stats_V1_works_lowSST';% SMALL REVERSAL
%------------------------------------------
% choose parameters for J+ scan and potential wells calculation
%------------------------------------------
if strcmp(Opt,'Small2Reversal')
    Jzero=2; Jstop=2.8; % values of J+ where MFT fixed point scan starts and stops
    Jplus=2.5; JplusLow=2.2; % large (beyond bifurcation) and small (before bifurcation) values of J+ for calculating linear response etc
    % Jplus is the value chosen for calculating effective MFT leading to potential wells, gain, etc.
elseif strcmp(Opt,'Small2_stats_V1')
    Jzero=2.8; Jstop=4; % values of J+ where MFT fixed point scan starts and stops
    Jplus=3.4; JplusLow=3; % large (beyond bifurcation) and small (before bifurcation) values of J+ for calculating linear response etc
    % Jplus is the value chosen for calculating effective MFT leading to potential wells, gain, etc.
elseif strcmp(Opt,'Small2_stats_V1_works')
    Jzero=3.4; Jstop=4.2; % values of J+ where MFT fixed point scan starts and stops
    Jplus=3.9; JplusLow=3.5; % large (beyond bifurcation) and small (before bifurcation) values of J+ for calculating linear response etc
    % Jplus is the value chosen for calculating effective MFT leading to potential wells, gain, etc.
elseif strcmp(Opt,'Small2_stats_V1_works_lowSST')
    Jzero=3; Jstop=4; % values of J+ where MFT fixed point scan starts and stops
    Jplus=3.9; JplusLow=3.5; % large (beyond bifurcation) and small (before bifurcation) values of J+ for calculating linear response etc
    % Jplus is the value chosen for calculating effective MFT leading to potential wells, gain, etc.
end
%---------------------
% CUE OPTIONS
%---------------------
% cue_mode=1 runs MFT with the same network, but different external
% perturbations referred to as 'cue'
%
cue_mode=1; % 0,1 for no cue or yes cue
cue=struct('cue_option',[],'cue_stat',[],'cue_value',[]);
% below, only active if cue_mode=1
cue.cue_option={'VIP'}; % populations affected by the cue, among: 'Pyr','SST','PV','VIP'
cue.cue_stat={'mean'}; % for each population above, choose type of cue: 
% 'gaussian'=increase quenched variance of afferent currents as in Mazzucato et al., 2019
% 'mean'=change mean afferent current
cueset={[0,0.1]}; % for each cue_option, write the set of values in each MFT RUN
stimset=[0]; % sensory stimulus, currently disabled
nbins_grid=[50 50]; % for the effective MFT calculation, # of bins in firing rate space to make a square grid for clusters 1 and 2 [NOTE: large values increase computation time quadratically!!]

% pack options into a structure:
extra=cell(1,numel(cueset{1}));% extra=''; % 
if cue_mode
    for i_c=1:numel(cueset{1})
        for i=1:numel(cue.cue_option)
            extra{i_c}=[extra{i_c} '_' cue.cue_option{i} '[' cue.cue_stat{i}(1) '' num2str(cueset{i}(i_c)) ']'];
        end
    end
end
tmpv=struct('Jplus',Jplus,'Opt',Opt,'nbins_grid',nbins_grid,'cue_mode',cue_mode,'cueset',{cueset},'cue',{cue},'stimset',stimset,'extra',{extra},'RUN',RUN);

%% run scripts
MFT_run;

