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
Opt='Small2Reversal';%'Small2InhDens';%'Small2_stats_V1';%

%------------------------------------------
% choose parameters for J+ scan and potential wells calculation
%------------------------------------------
Jzero=2; Jstop=3; % values of J+ where MFT fixed point scan starts and stops
Jplus=2.5; JplusLow=2.2; % large (beyond bifurcation) and small (before bifurcation) values of J+ for calculating linear response etc
% Jplus is the value chosen for calculating effective MFT leading to potential wells, gain, etc.

%-------------------------------------------
% OPTIONS FOR CHANGING CELL TYPE DENSITIES
%-------------------------------------------
cue_mode=0; % 0,1 for no cue or yes cue

% below: change PV/SST ratio
% InhPops={'PV','SST','VIP'};
% pset=[0.4,0.5,0.6]; % PV/(PV+SST)
pset=[0.5,0.6];%,0.5,0.6]; % PV/(PV+SST)
indInh=1;
%Inh density
InhDensity=cell(1,numel(pset));
InhD=[0.4,0.4,0.2]; % baseline fractions
for ip=1:numel(pset)
    PVpSST=1-InhD(3); 
    InhDtemp=[pset(ip)*PVpSST,(1-pset(ip))*PVpSST,InhD(3)];
    InhDensity{ip}=InhDtemp;
end
Type='PVbySST'; 
extra0=cell(1,numel(pset)); for ip=1:numel(pset); extra0{ip}=num2str(pset(ip)); end

% % below: change one density, and keep the sum fixed
% InhPops={'PV','SST','VIP'};
% pset=[0,1];
% indInh=2;
% %Inh density
% InhDensity=cell(1,numel(pset));
% InhD=[5,5,2.5]; InhD=InhD/sum(InhD); % baseline
% for ip=1:numel(pset)
%     InhDtemp=InhD; InhDtemp(indInh)=InhDtemp(indInh)*(1+pset(ip)); 
% %     InhDtemp=InhDtemp/sum(InhDtemp);
%     InhDensity{ip}=InhDtemp;
% end
% Type=['rho' InhPops{indInh}];
% extra0=cell(1,numel(pset)); for ip=1:numel(pset); extra0{ip}=num2str(pset(ip)); end

% % below: layer distribution
% Type='L';
% % InhPops={'PV','SST','VIP'};
% % pset=[0.4,0.5,0.6]; % PV/(PV+SST)
% %Inh density
% InhDensity=cell(1,4);
% InhD=[0.4,0.4,0.2]; % baseline fractions
% PV=[0.01,0.18,0.33,0.33,0.15];
% SST=[0.04,0.15,0.22,0.33,0.26];
% VIP=[0.13,0.43,0.26,0.1,0.08];
% for ip=1:numel(InhDensity)
%     InhDtemp=InhD.*[PV(ip+1),SST(ip+1),VIP(ip+1)];InhDtemp=InhDtemp/sum(InhDtemp);
%     InhDensity{ip}=InhDtemp;
% end
% extra0=cell(1,numel(pset)); for ip=1:numel(pset); extra0{ip}=num2str(InhDensity(ip)); end

stimset=[0]; % sensory stimulus, currently disabled
nbins_grid=[50 50]; % for the effective MFT calculation, # of bins in firing rate space to make a square grid for clusters 1 and 2 [NOTE: large values increase computation time quadratically!!]

% pack options into a structure:
extra=cell(1,numel(InhDensity));% extra=''; % 
for i_c=1:numel(InhDensity)
    %         for i=1:numel(cue.cue_option)
    %             extra{i_c}=[extra{i_c} '_' cue.cue_option{i} '[' cue.cue_stat{i}(1) '' num2str(cueset{i}(i_c)) ']'];
    %         end
    extra{i_c}=[Type '' extra0{i_c} ''];
end
tmpv=struct('Jplus',Jplus,'Opt',Opt,'nbins_grid',nbins_grid,'InhDensity',{InhDensity},'cue_mode',cue_mode,'stimset',stimset,'extra',{extra},'RUN',RUN);



%%

if RUNMFT
    for i_c=1:numel(InhDensity)
        options_main=struct('a',{'sel'},'Opt',Opt,'Jstop',Jstop,'Jzero',Jzero,'extra',extra{i_c},'cue_mode',cue_mode,...
            'InhDensity',InhDensity{i_c},'S',0,'stimset',[0],'screen',0);
        tic
        options_main.a='sel';%'sel';%'spont';%
        MFTGainModulation_main(options_main);
        aux.MFTGainModulation_plot(options_main);
        options_main.a='all';%'sel';%'spont';%
        MFTGainModulation_main(options_main);
        aux.MFTGainModulation_plot(options_main);
        toc
    end
    % end
end
% plot fixed points scanned above
if PLOTMFT
    for i_c=1:numel(InhDensity)
        options_main=struct('a',{'sel'},'Opt',Opt,'Jstop',Jstop,'Jzero',Jzero,'extra',extra{i_c},'cue_mode',cue_mode,...
            'InhDensity',InhDensity{i_c},'S',0,'stimset',[0],'screen',0);        options_main.a='sel';%'sel';%'spont';%
        aux.MFTGainModulation_plot(options_main);
        options_main.a='all';%'sel';%'spont';%
        aux.MFTGainModulation_plot(options_main);
    end
end


if numel(InhDensity)>1
    tmpv1=tmpv;
    tmpv1.JplusLow=JplusLow;
    % plot changes in firing rates for each perturbation, linear response
    aux.MFT_fun_rate_changes(tmpv1);
    %% run or load results
    if RUN
        options=repmat(struct('Opt',Opt,'extra',[],...
            'Jplus',Jplus,'cue_mode',cue_mode,'nbins_grid',nbins_grid),numel(InhDensity),numel(stimset));
        for i_c=1:numel(InhDensity)
%             for i=1:numel(cue.cue_option); cue.cue_value{i}=cueset{i}(i_c); end
            for i_J=1:numel(stimset)
%                 options(i_c,i_J).cue=cue;%;cueset(i_c);
                options(i_c,i_J).InhDensity=InhDensity{i_c};%;cueset(i_c);
                options(i_c,i_J).stimset=stimset(i_J);
                options(i_c,i_J).extra=extra{i_c};
            end
        end
        Thres=[]; % compute thresholds on first run
        for i_c=1:numel(InhDensity)
            for i_J=1:numel(stimset)
                tic
                % calculate effective MFT and potential wells for the first 2
                % excitatory clusters, integrating out all other
                % populations
                Thres=aux.MFT_fun_amitmascaro_run_Stim_quench(options(i_c,i_J),Thres);
                toc
            end
        end
    end
    %%
    % plot each energy potential separately
    if PLOTMFT
        aux.MFT_fun_amitmascaro_plot_Stim(tmpv);
    end
end