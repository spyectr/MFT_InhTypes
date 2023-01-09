
%
% Luca Mazzucato January 2023
clear all
RUNMFT=1;
RUNWELLS=1;
PLOTMFT=1;



% CHOOSE OPTIONS FOR MFT
% 'Opt' is the network name. All networks with their parameters are stored
% in the file 'create_params.m'
Opt='Marcel2';%'Marcel';%'Small2ReversalLowSST';%'Small2Reversal';%'Small2InhDens';%'Small2_stats_V1';%
Jzero=7; % lowest intra-cluster J+ value to start MFT scan
Jstop=8; % highest intra-cluster J+ value to stop MFT scan
Jplus=7.9; % calculates linear response in the regime with 2 stable states beyond the pitchfork bifurcation: 
            % one with 1 cluster active and the rest inactive (selective activity), and one where all clusters have the same rate
JplusLow=7.2; % calculates linear response in the regime with only 1 stable states, before the pitchfork bifurcation


% CUE OPTIONS
% note here we can either use the option 'cue_mode=0' as in the script
% MFT_selective.m; alternatively, we can set up cue_mode=1, but choose a
% value of the cue equal to zero with cueset={[0]}. This is useful because
% the results will be used by the script MFT_energylandscape when comparing
% a nonzero with a zero cue.

%
cue_mode=1; % 0,1 for no cue or yes cue
cue=struct('cue_option',[],'cue_stat',[],'cue_value',[]);
% below, only active if cue_mode=1
cue.cue_option={'Pyr'};%{'VIP'};%{'SST'};%{'PV'};%
% % Bernstein plots v1
% cue.cue_stat={'noise'}; %{'gaussian'};% 'mean', , 'noise'
% cueset={[0,0.1,0.2,0.3]}; % for each cue_option, write the set of values in each MFT RUN
% Bernstein plots v2
cue.cue_stat={'gaussian'}; %{'gaussian'};% 'mean', , 'noise'
cueset={[0]}; % for each cue_option, write the set of values in each MFT RUN
%%

if RUNMFT
    tic
    options_main=struct('a','sel','Opt',Opt,'Jstop',Jstop,'Jzero',Jzero,'cue_mode',cue_mode,'S',0,'stimset',0,'screen',0);
    aux.MFTGainModulation_main_spu(options_main); %  this function runs the MFT scan from Jzero to Jstop
    pltMFT.MFTGainModulation_plot(options_main); % this function plots the MFT results showing the attractor firing rate as a function of J+
    options_main.a='all';%'sel';%'spont';%
    % second, find all the multistable solutions with different sets of
    % active clusters
    aux.MFTGainModulation_main_spu(options_main);
    pltMFT.MFTGainModulation_plot(options_main);
    toc
end
% plot fixed points scanned above
if PLOTMFT
    options_main=struct('a','sel','Opt',Opt,'Jstop',Jstop,'Jzero',Jzero,'cue_mode',cue_mode,'S',0,'stimset',0,'screen',0);
    pltMFT.MFTGainModulation_plot(options_main);
    options_main.a='all';%'sel';%'spont';%
    pltMFT.MFTGainModulation_plot(options_main);
end
