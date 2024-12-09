% Demo script for running spiking network simulations and analyses
% 
% by Luca Mazzucato 2023
%
% ----------------------------------------
% This script run simulations of LIF networks with excitatory (E) and inhibitory (I) spiking neurons.
%
%% CHOOSE CLUSTERED ARCHITECTURE
%------------------------
% Opt='EI';%
% Opt='Marcel0';%
% Opt='MarcelPost_disinh';%
% Opt='MarcelPost_inh';%
Opt='Marcel2';%
% Opt='E';%
%------------------------
% LOAD PARAMETERS
%------------------------
paramsfile='params.mat'; % file where all network parameters are saved
create_params(paramsfile,Opt);

stimuli={}; % ongoing activity
% stimuli={'US'}; % stimulus evoked-activity targeting selective clusters
% stimuli={'CSgauss'}; % anticipatory cue speeds up network dynamics
% stimuli={'US','CSgauss'}; % anticipatory cue preceeds stimulu delivery
savedir=fullfile('results',Opt); if ~exist(savedir,'dir'); mkdir(savedir); end % setup directory for saving HMM data
save(paramsfile,'stimuli','savedir','-append');

%% RUN SIMULATION
ntrials=1; % number of trials
file_sim=fullfile(savedir,'results.mat');  % file where simulation results are saved
%---------------------------
% GENERATE SYNAPTIC WEIGHTS
%---------------------------
% J = N x N matrix of synaptic weights
% params = structure containing all network parameters
% [J, params]=auxi.fun_SynWeights(paramsfile,Opt);
[J, params]=auxi.fun_SynWeights_Marcel(paramsfile,Opt);
% end
[stimulus_save, params]=auxi.fun_stim(params); % STIMULUS
%------------------------
% SIMULATION
%%
% ------------------------
tic
firings=cell(1,ntrials); % cell array with all spike times in each trial
PlotData=cell(1,ntrials); % cell array with data for plotting
% parfor iTrial=1:ntrials % uncomment this line if you have a multi-core  machine with 4 or more cores
for iTrial=1:ntrials
    ParamsRun=params;
    ParamsRun.Ext=stimulus_save.Ext;
    ParamsRun.Stimulus=stimulus_save.Stimulus;
    ParamsRun.J=J;
    fprintf('--- Start SIM ...\n');
    if strcmp(ParamsRun.Stimulus.input,'Poisson')
        [firings{iTrial}, PlotData{iTrial}]=auxi.fun_LIF_SIM_Poisson(ParamsRun);
        
    else
    [firings{iTrial}, PlotData{iTrial}]=auxi.fun_LIF_SIM(ParamsRun);
    end
end
% SAVE results
save(file_sim,'params','firings','PlotData','stimulus_save');
fprintf('\nDone. Simulation saved in %s\n',file_sim);
toc
%%
%------------------------
% PLOT EVENTS
%------------------------
iTrial=1; % pick which trial to plot
dataload=load(file_sim); % load simulation results
data=dataload.PlotData{iTrial}; % membrane potentials
firings=dataload.firings{iTrial}; % spikes
Params=dataload.params; % parameters
Params.Ext=dataload.stimulus_save.Ext; % external currents
Params.savedir=savedir;
plotclass.fun_PlotTrial(data,firings,Params);
