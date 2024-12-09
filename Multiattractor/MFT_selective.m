% Luca Mazzucato January 2023

clear all
RUNMFT=1; % run and plot MFT
PLOTMFT=1; % plot MFT results

% CHOOSE OPTIONS FOR MFT
% 'Opt' is the network name. All networks with their parameters are stored
% in the file 'create_params.m'
Opt='Marcel2';%'Marcel';%'Small2ReversalLowSST';%'Small2Reversal';%'Small2InhDens';%'Small2_stats_V1';%
Jzero=7; % lowest intra-cluster J+ value to start MFT scan
Jstop=10; % highest intra-cluster J+ value to stop MFT scan

cue_mode=0; % 0,1 for no cue or yes cue
%%
options_main=struct('a','sel','Opt',Opt,'Jstop',Jstop,'Jzero',Jzero,'cue_mode',cue_mode,'S',0,'stimset',0,'screen',0);
if RUNMFT
    tic
    options_main.a='sel';
    aux.MFTGainModulation_main_spu(options_main); %  this function runs the MFT scan from Jzero to Jstop
    pltMFT.MFTGainModulation_plot(options_main); % this function plots the MFT results showing the attractor firing rate as a function of J+
    % second, find all the multistable solutions with different sets of
    % active clusters
    options_main.a='all';%'sel';%'spont';%
    aux.MFTGainModulation_main_spu(options_main);
    pltMFT.MFTGainModulation_plot(options_main);
    toc
end
% plot fixed points scanned above
if PLOTMFT
    options_main=struct('a','sel','Opt',Opt,'Jstop',Jstop,'Jzero',Jzero,'cue_mode',cue_mode,'S',0,'stimset',0,'screen',0);
    pltMFT.MFTGainModulation_plot(options_main);
    options_main=struct('a','all','Opt',Opt,'Jstop',Jstop,'Jzero',Jzero,'cue_mode',cue_mode,'S',0,'stimset',0,'screen',0);
    pltMFT.MFTGainModulation_plot(options_main);
end

