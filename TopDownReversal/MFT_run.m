%% 
% scan fixed points while varying J+
if RUNMFT
    for i_c=1:numel(cueset{1})
        for i=1:numel(cue.cue_option); cue.cue_value{i}=cueset{i}(i_c); end
        options_main=struct('a',{'spont'},'Opt',Opt,'Jstop',Jstop,'Jzero',Jzero,'extra',extra{i_c},'cue_mode',cue_mode,'cue',cue,'S',0,'stimset',[0],'screen',0);
        tic
        options_main.a='sel';%'sel';%'spont';% persistent selective activity
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
    for i_c=1:numel(cueset{1})
        for i=1:numel(cue.cue_option); cue.cue_value{i}=cueset{i}(i_c); end
        options_main=struct('a',{'spont'},'Opt',Opt,'Jstop',Jstop,'Jzero',Jzero,'extra',extra{i_c},'cue_mode',cue_mode,'cue',cue,'S',0,'stimset',[0],'screen',0);
        options_main.a='sel';%'sel';%'spont';%
        aux.MFTGainModulation_plot(options_main);
        options_main.a='all';%'sel';%'spont';%
        aux.MFTGainModulation_plot(options_main);
    end
end

%%
if cue_mode
    tmpv1=tmpv;
    tmpv1.JplusLow=JplusLow;
    % plot changes in firing rates for each perturbation, linear response
    aux.MFT_fun_rate_changes(tmpv1);
    %% run or load results
    if RUN
        options=repmat(struct('Opt',Opt,'extra',[],...
            'Jplus',Jplus,'cue_mode',cue_mode,'nbins_grid',nbins_grid),numel(cueset{1}),numel(stimset));
        for i_c=1:numel(cueset{1})
            for i=1:numel(cue.cue_option); cue.cue_value{i}=cueset{i}(i_c); end
            for i_J=1:numel(stimset)
                options(i_c,i_J).cue=cue;%;cueset(i_c);
                options(i_c,i_J).stimset=stimset(i_J);
                options(i_c,i_J).extra=extra{i_c};
            end
        end
        Thres=[]; % compute thresholds on first run
        for i_c=1:numel(cueset{1})
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