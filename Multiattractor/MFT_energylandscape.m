
%
% Luca Mazzucato January 2023
clear all
RUNMFT=0;
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
%
cue_mode=1; % 0,1 for no cue or yes cue
cue=struct('cue_option',[],'cue_stat',[],'cue_value',[]);
% below, only active if cue_mode=1
cue.cue_option={'Pyr'};%{'VIP'};%{'SST'};%{'PV'};%
% % Bernstein plots v1
% cue.cue_stat={'noise'}; %{'gaussian'};% 'mean', , 'noise'
% cueset={[0,0.1,0.2,0.3]}; % for each cue_option, write the set of values in each MFT RUN
% Bernstein plots v2
cue.cue_stat={'mean'}; %{'gaussian'};% 'mean', , 
cueset={[0]}; % for each cue_option, write the set of values in each MFT RUN



% cueset={[0]}; % for each cue_option, write the set of values in each MFT RUN

stimset=[0];%,0.02,0.04];

nbins_grid=[40 40];


% cue mode
extra=cell(1,numel(cueset{1}));
if cue_mode
    for i_c=1:numel(cueset{1})
        for i=1:numel(cue.cue_option)
            extra{i_c}=[extra{i_c} '_' cue.cue_option{i} '[' cue.cue_stat{i}(1) '' num2str(cueset{i}(i_c)) ']'];
        end
    end
end
tmpv=struct('Jplus',Jplus,'Opt',Opt,'nbins_grid',nbins_grid,'cue_mode',cue_mode,'cueset',{cueset},'cue',{cue},'stimset',stimset,'extra',{extra},'RUN',RUNWELLS);

%%
MFT_run;