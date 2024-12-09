function OUT=MFT_RATE_eff(ni_start,Params,nu_grid,flg)
% Params=Params0;

% global return_value  NTRIAL TOLX TOLF screen contr overlap  nivsj_mode stim_mode syn_dyn nfocus mnewt_flag cue_mode

OUT=struct('rates_eff',[],'nu_out',[],'nu_in',[]);

flg.return_value=1;
flg.screen=0;
flg.mnewt_flag=0;
flg.nfocus=0;
flg.NTRIAL=1000;
flg.TOLX=1.e-6;
flg.TOLF=1.e-6;

ind_mft=Params.ind_mft;
ind_focus=Params.ind_focus; % populations in focus
fieldNames={'fieldNames','A','B','Mu_ext','Sigma_ext','H','Tau','Tausyn','Theta','tau_arp','BS','Files'};
tau_arp=Params.tau_arp; 
BS=Params.BS;
A=Params.A(ind_mft,ind_mft);
B=Params.B(ind_mft,ind_mft);
Mu_ext=Params.Mu_ext(ind_mft);
Sigma_ext=Params.Sigma_ext(ind_mft);
H=Params.H(ind_mft);
Tau=Params.Tau(ind_mft);
Tausyn=Params.Tausyn(ind_mft);
Theta=Params.Theta(ind_mft);
Files=Params.Files;
Params_eff=utils.v2struct(fieldNames);
if flg.cue_mode
    if any(strcmp(Params.Network.cue.cue_stat,'gaussian'))
        tempquench_pops=find(Params.quench_pops==ind_mft); % 
        Params_eff.Mu_extZ=cell(1,numel(tempquench_pops));
        for i=1:numel(tempquench_pops)
            Params_eff.Mu_extZ{i}=Params.Mu_extZ{tempquench_pops(i)};
            Params_eff.quench_pops(i)=tempquench_pops(i)-numel(ind_focus);
        end
    end
%     Params_eff.Network.cue_value=Params.Network.cue_value;
end
Params_eff.Network=Params.Network;
n=numel(ind_mft);
flag=0;
stop_control = 1;
nbins_grid=[length(nu_grid(1).focus),length(nu_grid(2).focus)];
rates_eff=NaN(n,nbins_grid(1),nbins_grid(2));
nu_out=NaN(numel(ind_focus),nbins_grid(1),nbins_grid(2));
nu_in=NaN(numel(ind_focus),nbins_grid(1),nbins_grid(2));
eigenvalues.Real=NaN(n,nbins_grid(1),nbins_grid(2));
eigenvalues.Imag=NaN(n,nbins_grid(1),nbins_grid(2));
% eigenvectors=[]; % dim: 1st=pop #; 2nd=eigenvector #; 3rd=Jplus step
for i1=1:nbins_grid(1)
    for i2=1:nbins_grid(2)
        nu_focus=[nu_grid(1).focus(i1) nu_grid(2).focus(i2)];
        nu_in(:,i1,i2)=nu_focus;
        if flg.screen; fprintf('--- focus: nu_1=%0.03g, nu_2=%0.03g\n',nu_focus(1),nu_focus(2)); end;
        % contribution from pops in focus (shifting Mu_ext and
        % Sigma_ext for mft pops)
        mu_eff=Params.A(ind_mft,ind_focus)*nu_focus';
        sigma_eff=Params.B(ind_mft,ind_focus)*nu_focus';
        Params_eff.Mu_ext=Params.Mu_ext(ind_mft)+mu_eff';
        Params_eff.Sigma_ext=Params.Sigma_ext(ind_mft)+sigma_eff';
        % initial conditions
        if i1==1 && i2==1
            ni=ni_start;
        elseif i1>1 && i2==1
            ni=squeeze(rates_eff(1:n,i1-1,i2))';
        elseif i2>1
            ni=squeeze(rates_eff(1:n,i1,i2-1))';
        end    
%         ni=ni_start;
    % Help mnewt search for low rate fixed points (usually not needed though; goAB_simul could be used as well)

        if flg.stim_mode
            [ni,flg]=auxMFT.goAB_simul(ni,1000,Params_eff,flg);
        end%     if (flg.return_value == 0)  BYE

        nitemp=ni;
        [ni,flg]=auxMFT.mnewt(flg.NTRIAL,ni,Params_eff,flg);
        if (flg.return_value == 0)
            if flg.screen
                fprintf('MFT_RATE_eff: mnewt \n'); % debug
                fprintf('\n>>> Error in routine mnewt() called by MFT_RATE_eff()->run goAB_eff then try mnewt again\n');
            end
            flg.return_value=1;
            [ni,flg]=auxMFT.goAB_simul(nitemp,10000,Params_eff,flg);
            [ni,flg]=auxMFT.mnewt(flg.NTRIAL,ni,Params_eff,flg);
%             continue;
        end
        % handle unstable points:
        % note: the stability analysis doesn't seem to work when the
        % perturbation is inh_gaussian
        [res, Real_part, Imag_part, ~,flg]=auxMFT.Test_Stability(ni,Params_eff,flg);
        % store to save
        eigenvalues.Real(:,i1,i2)=Real_part;
        eigenvalues.Imag(:,i1,i2)=Imag_part;
    %     eigenvectors=cat(3,eigenvectors,EigV);
        %
        if ((~res) && (stop_control == 1)) 
            fprintf('--- Unstable state occurred \n');
            fprintf('\n    Do Dynamics Simulation.\n');
            ni =ni + .001;  % Little perturbation
            [ni,flg]=auxMFT.goAB_simul(ni,0,Params_eff,flg); % 0 := no max number of steps.
        end
        rates_eff(1:n,i1,i2)=ni;
        if flg.screen
            fprintf('--- frozen:');
            for i=1:numel(ind_mft)
                fprintf(' nu(%d)=%0.03g',ind_mft(i),ni(i));
            end
            fprintf('\n');
        end
        % compute nu_out with using (nu_1,nu_2,ni)
        ni_in=[nu_focus,ni];
        [risp_out,flg]=auxMFT.RISP(ni_in,Params,flg);
        nu_out(1:numel(ind_focus),i1,i2)=risp_out(ind_focus);
    end 
end

% Summary of parameters:
fprintf('\n--- Done.\n');
% Summary(paramsfile);

OUT.rates_eff=rates_eff;
OUT.nu_out=nu_out;
OUT.nu_in=nu_in;

return;

%---------------
% SAVE TO FILE
%---------------
% file to store fixed points
file_name=[Files.file_amitmascaro '.mat'];
% file to store eigenvalues and eigenvectors
% fixed point rates
AllParameters=load(Files.paramsfile);
save(file_name,'nu_grid','rates_eff','eigenvalues','AllParameters');
% all eigenvalues




