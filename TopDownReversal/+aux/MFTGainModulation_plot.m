function MFTGainModulation_plot(options_main)

aux.MFT_preamble;
Plotf=1;        % default: plot
NumFig=1;
if flg.nivsj_mode==1 && flg.all_mode==0
    fprintf('\n--- selective activity...\n');
    nn=6; % 1 exc selective, 1 exc non-selective, 1 exc bg, 3 inh
    save(flg.paramsfile,'nn','-append');
else
    fprintf('\n--- spontaneous activity...\n');
end

if flg.nivsj_mode==0 && flg.all_mode==1
    % if all mode, compute spontaneous activity for later use
    nn=p+3;
    save(flg.paramsfile,'nn','-append');
end
if (any(bitand(flg.flags,flg.RATE_VS_J)) && flg.all_mode==1)% && flg.stim_mode==0
    flg.nivsj_mode=1;
    flg.contr=0;
    if ~isempty(strfind(Opt,'Small2'))
        flg.contr=1;
    end
    fprintf('\n--- Starting All_mode...\n');
    nn=p+4;
    save(flg.paramsfile,'nn','Jplus','-append');
end
%%
%-----------------------
% PLOT
%-----------------------
if Plotf
    if flg.stim_mode==0 && flg.all_mode==0
        if (flg.spur_mode==0) && (flg.all_mode==0)
            pltMFT.PlotFixedPoints(nn,flg.paramsfile,Jstop,NumFig,flg);
        elseif flg.spur_mode==1
            % plot spurious states
            pltMFT.PlotAllFixedPoints(p,Jstop,NumFig,flg);
            NumFig=NumFig+1;
%             PlotClusterRates(p,Params,Jstop,NumFig)
        end
    elseif flg.all_mode==1 || flg.stim_mode==0
        % plot all stable states
        pltMFT.PlotStableStates(p,Jstop,NumFig,flg);
        NumFig=NumFig+1;
%         pltMFT.PlotStableStates_semilogy(p,Jstop,NumFig,flg);
    elseif flg.stim_mode==1
        pltMFT.PlotStableStatesSTIM(p,Params,Jstop,NumFig,flg);
    end
end