% transfer function of sel cluster as a function of stimulus,
% with and without cue.
% 3 pop: sel, nonsel, inh
% for each value of stim, find MFT solution for the 3 pop
% plot all of them
% repeat with and without cue.
%
% EXAMPLE: for Small2_stats, choose cues 0.01:0.01:0.05, set Jplus=11.4.
% first run MFT_main_cue with options: with gaussian
% -xSmall2_stats -af11.5 -Jzero11 -Cg0.01
% ...
% -xSmall2_stats -af11.5 -Jzero11 -Cg0.05
% with inh cue
% -xSmall2_stats -af11.5 -Jzero11 -Ci0.01
% then run this script.
%
% Luca Mazzucato December 2016
clear all
RUNMFT=1;
RUN=1;
PLOTMFT=1;

% Opt='Small2InhDens';%'Small2_stats_V1';%'Small2_InhStab';%
% Jplus=4.8; Jstop=5; Jzero=4; JplusLow=4.1;
% Opt='Small2_V1_HiRates';%'Small2InhDens';%'Small2_stats_V1';%'Small2_InhStab';%
% Jplus=7; Jstop=8; Jzero=5; JplusLow=5.5;
% Opt='Small2Reversal';%'Small2_stats_V1';%
% Jplus=2.6; Jstop=3; Jzero=2; JplusLow=2.2;
% Opt='Small2_stats_bernstein';%'Small2InhDens';%'Small2_stats_V1';%
Opt='Small2_stats_cosyne';%'Small2InhDens';%'Small2_stats_V1';%
Jplus=10.3; Jstop=10.8; Jzero=9.8; JplusLow=9.8;
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
cue.cue_stat={'gaussian'}; %{'gaussian'};% 'mean', , 'noise'
cueset={[0,0.1]}; % for each cue_option, write the set of values in each MFT RUN



% cueset={[0]}; % for each cue_option, write the set of values in each MFT RUN

stimset=[0];%,0.02,0.04];

nbins_grid=[80 80];


% cue mode
extra=cell(1,numel(cueset{1}));
if cue_mode
    for i_c=1:numel(cueset{1})
        for i=1:numel(cue.cue_option)
            extra{i_c}=[extra{i_c} '_' cue.cue_option{i} '[' cue.cue_stat{i}(1) '' num2str(cueset{i}(i_c)) ']'];
        end
    end
end
tmpv=struct('Jplus',Jplus,'Opt',Opt,'nbins_grid',nbins_grid,'cue_mode',cue_mode,'cueset',{cueset},'cue',{cue},'stimset',stimset,'extra',{extra},'RUN',RUN);

%%
MFT_run;