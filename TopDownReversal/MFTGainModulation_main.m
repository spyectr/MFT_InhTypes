%        %***************************************************
% 	  *                                                        *
% 	  *              MFT main script for finding fixed points  *
% 	  *                                                        *
% 	  ***************************************************%
% 
%---------------------------------------------------------------
% The program saves and loads running parameters in DATA/params.mat and overwrites
% thresholds and other default parameters when prompted.
%

function MFTGainModulation_main(options_main)

aux.MFT_preamble

fprintf('# of clusters: %d',p);

temp_cue_mode=flg.cue_mode;

Plotf=1;        % default: plot

Files=auxMFT.File_Info_MFT(nn,p,flg);
if exist(Files.Log,'file')
    delete(Files.Log);
end
diary(Files.Log);

% keep flg.stim_mode=0 for Fix_Threshold
temp_stim_mode=flg.stim_mode;
flg.stim_mode=0;

nn_temp = nn;
%%    ------------------------------------------------------------------------
%     STARTING MEAN FIELD PROGRAM: -------------------------------------------
%     SEARCHING FOR THRESHOLDS.    -------------------------------------------
%     ------------------------------------------------------------------------%

fprintf('\n\n\n    *******************      Mean Field Program      ********************\n\n\n\n');
fprintf('--- Searching For Thresholds... \n');
nn=4; % 2 populations while searching for thresholds - nn_temp holds the actual # of pops

%   % Frequencies we are looking for stable spontaneous states:
Ni=zeros(1,nn);
Ni(1) = ni_e;   % 3.0
Ni(2:4) = ni_i;   % 4.2
fprintf('    ... for ni_e=%g   ni_i=',ni_e);for i=1:nn-1; fprintf('%g,',ni_i(i)); end
fprintf('   ni_ext_e=%g   ni_ext_i=%g\n',ni_ext_e,ni_ext_i(1));
%-----------------------------------------------------------------
% FIX THRESHOLDS in such a way wanted spontaneous point is stable with nn populations.
% (If you want ni_ext=ni_e set second parameter =1)
%-----------------------------------------------------------------%

params=load(flg.paramsfile);
params.paramsfile=flg.paramsfile;
[theta,flg]=aux.fun_Fix_Threshold_InhType(Ni,0,params,flg);
% update parameters file with current thresholds
if (flg.return_value == 0) ; error('\n>>> Error in routine Fix_Threshold called by main()\n\n'); end
fprintf('\n    Done. Thresholds found: (E) %g, (I) ',theta(1)); for i=2:nn; fprintf('%g,',theta(i)); end; fprintf('\n\n');

% Now reinitialize parameters with new thresholds:
Params=aux.In_Param_4pop_InhType(theta,flg.paramsfile,flg);
pltMFT.fun_plot_weights(Params);
% file names are initialized inside In_Param_4pop_InhType
% Params.A
% Params.B
% Research fixed spontaneous point for consistency...:
% flg.screen=1; DEBUG=1; syn_dyn=0;
fprintf('\nCheck stability of state used for finding thresholds...\n');
[Ni,flg]=auxMFT.mnewt(flg.NTRIAL,Ni,Params,flg);
% (if state is instable, exit to program): 
[res, Real_part, Imag_part, EigV, flg]=auxMFT.Test_Stability(Ni,Params,flg);
if (~res)  
    fprintf('>>> State used to calculate thresholds is unstable,.. looking for stable state\n');
    [~,unind]=max(Real_part);
    unstabledir=EigV(:,unind)';
    [Ni,flg]=auxMFT.goAB_simul(Ni+10*unstabledir,0,Params,flg); % 0 := no max number of steps.
    [res, Real_part, Imag_part, EigV,flg]=auxMFT.Test_Stability(Ni,Params,flg);
    if (res) 
        fprintf('>>> The closest stable state is\n');
        fprintf('    %f\n',Ni);
        ni_e=Ni(1); ni_i=Ni(2:end);
        fprintf('>>> Save this state as ni_e, ni_i starting point for analysis...\n');
        save(flg.paramsfile,'ni_e','ni_i','-append');
    else
        fprintf('>>> Error: still unstable! QUIT\n');
        DB_EXIT(flg);
    end
else
    fprintf('Stable\n');
end
% recast original flg.stim_mode value
flg.stim_mode=temp_stim_mode;
flg.cue_mode=temp_cue_mode;
NumFig=1;
nn=nn_temp;
%%
% -----------------------------------------------------------------------
% RATE_VS_J mode: -------------------------------------------------------
% -----------------------------------------------------------------------
% this part calculates the selective persistent activity in Amit-Brunel fashion to
% find the bifurcation point
cond=any(bitand(flg.flags,flg.RATE_VS_J)) & flg.spur_mode==0 & flg.stim_mode==0;
if cond%(any(bitand(flg.flags,flg.RATE_VS_J)) && flg.spur_mode==0) && flg.stim_mode==0
    fprintf('\n--- Starting RATE_VS_J mode...\n');
    if flg.nivsj_mode==1
        fprintf('\n--- selective activity...\n');
        nn=6; % 1 exc selective, 1 exc non-selective, 1 exc bg, 3 inh
        save(flg.paramsfile,'nn','-append');
    else
        fprintf('\n--- spontaneous activity...\n');
    end
    
    if flg.nivsj_mode==0 && flg.all_mode==1
        % if all mode, compute spontaneous activity for later use
        nn=5;
        save(flg.paramsfile,'nn','-append');
    end
    Jplus=Jzero;
    save(flg.paramsfile,'Jplus','-append');
    Params=auxMFT.Learn(Params,flg);
    % plot synaptic weights
    pltMFT.fun_plot_weights(Params);

    %     % Starting point:
    Ni=zeros(1,nn);
    for i=1:nn-3
        Ni(i) = ni_e;
    end
    Ni(nn-2:nn) = ni_i;
% flg.screen=1; DEBUG=0; %syn_dyn=0;
    if (flg.screen); fprintf('--- Initial point:\n\n'); for i=1:nn; fprintf('    Ni[%d]: %f\n',i,Ni(i)); end; fprintf('\n'); end
    fprintf('\nFirst round of goAB_simul\n');
    [Ni,flg]=auxMFT.goAB_simul(Ni,10000,Params,flg);
    fprintf('First round of mnewt\n');
    [Ni,flg]=auxMFT.mnewt(flg.NTRIAL,Ni,Params,flg);
    if (flg.return_value == 0); DB_EXIT(flg);  end
    fprintf('Running RATE_VS_Jfun...\n');
    flg=aux.RATE_VS_Jfun(Ni,Jstop,Params,flg.STEP,Opt,flg);
    if (flg.return_value == 0); fprintf('\n>>> Error in routine RATE_VS_Jfun() called by main()\n'); end
end
% -----------------------------------------------------------------------
% ALL mode: -------------------------------------------------------
% -----------------------------------------------------------------------
% this part calculates the other fixed points generating the full attractor
% landscape
if (any(bitand(flg.flags,flg.RATE_VS_J)) && flg.all_mode==1)% && flg.stim_mode==0
    flg.nivsj_mode=1;
    flg.contr=0;
    if ~isempty(strfind(Opt,'Small2'))
        flg.contr=1;
    end
    fprintf('\n--- Starting All_mode...\n');
    nn=p+4;
    save(flg.paramsfile,'nn','Jplus','-append');
    Params=auxMFT.Learn(Params,flg);
    rates_old=zeros(1,nn);
    % save blank initial conditions in file
    save(Params.Files.rates_old,'rates_old');
    n_hiMax=20;
    n_hiMax=min([n_hiMax,p]);
    % store all stable states
    StableStates=[];     UnstableStates=[]; 
    StableStates(n_hiMax).Jplus=[];    StableStates(n_hiMax).rates=[];    StableStates(n_hiMax).Eig=[];
    UnstableStates(n_hiMax).Jplus=[];    UnstableStates(n_hiMax).rates=[];    UnstableStates(n_hiMax).Eig=[];
    Ni=zeros(1,nn);
    Ni(1:nn-3) = ni_e;
    Ni(nn-2:nn) = ni_i;   
    fprintf('First round of goAB_simul\n');
    [Ni,flg]=auxMFT.goAB_simul(Ni,20000,Params,flg);
    fprintf('First round of mnewt\n');
    [Ni,flg]=auxMFT.mnewt(flg.NTRIAL,Ni,Params,flg);
    for n_hi= 1:n_hiMax
        fprintf(' Running with %d hi populations...\n',n_hi);
        [StableStates(n_hi), UnstableStates(n_hi),flg]=aux.RATE_VS_Jfun_ALL(n_hi,Ni,Jstop,Params,Opt,flg);
        if (flg.return_value == 0)
            fprintf('\n>>> Error in routine RATE_VS_Jfun_ALL() called by main()\n');
        end
    end
    save(Params.Files.stable,'StableStates','UnstableStates');
end




% % -----------------------------------------------------------------------
% % STIM mode: -------------------------------------------------------
% % -----------------------------------------------------------------------
% % Set larger input current on selective populations
% % Default: use GC selectivity stats
% 
% if any(bitand(flg.flags,flg.RATE_VS_J)) && flg.stim_mode==1
%     % PRINT
%     fprintf('\n--- Starting STIMULUS mode...\n');    
%     SelPop=NumStim;%max(sum(selective,2)); % number of selective clusters (higher Ext current)
%     % height_stim=0.001;%box_height(flg.paramsfile);
%     fprintf('    %d stimulated clusters.',SelPop);
%     Network.SelPop=SelPop; % called by Learn>In_Parameters
%     Jplus=Jzero;
%     StimHeightSet=0;%[0.32 0.34];%[];
%     for i=1:numel(StimHeightSet)
%         StimHeight=StimHeightSet(i);
%         fprintf('    Stimulus gain: %d percent\n',round(StimHeight*100));
%         nip0=0; % initial # of ungrouped pops
%         % 
%         % STIMULUS GAIN: 
%         % 0 = no stimulus gain; 
%         % 0.1 = 10% gain, etc;
%         % [] = use gain from simulations through height_box function
%         % 
%         save(flg.paramsfile,'Network','Jplus','nip0','StimHeight','-append');    
%         % FIND FIXED POINTS WITH GROUPED POPS
%         flg.contr=0;
%         fprintf('\n--- Find fixed points using grouped pops...\n\n');
%         Params=Learn(Params,flg);
%         [Params,flg]=STIM_FIND(Params,Jstop,flg.STEP,flg);
% %         Params=STIM_FIND_reverse(Params,Jstop,flg.STEP);
%         flg.return_value=1;
%         % CHECK STABILITY WITH ALL UNGROUPED POPULATIONS
%         Params=STIM_Stability(Params,Jstop,flg.STEP);  
% %         Params=STIM_Stability_reverse(Params,Jstop,flg.STEP);  
%         PlotStableStatesSTIM(p,Params,Jstop,NumFig);
%     end
% end
% 
% 


% EXIT

if (flg.return_value == 0)   
    fprintf('\n--- Exit to program\n\n');
end

% load chirp 
% sound(y,Fs)
diary('off');


return

