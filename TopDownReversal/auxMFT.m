classdef auxMFT
    methods(Static)
        
        
        %----------------------------------------------------------------
        %----------------------------------------------------------------
        
        % Determine increment on Jplus, depending on network parameters, in order
        % to optimize behavior around critical points
        % INPUT: Jplus = current value of Jplus
        %          Opt = network option ('AB97','GLC',...)
        % OUTPUT: incremental step on Jplus
        
        function step=StepJplus(Jplus,Opt,STEP,flg)
            
            % global stim_mode
            
            %     % Refine search:
            % step = STEP;
            if strcmp(Opt,'Small2_stats_V1')
                if ((Jplus >= 3.) && (Jplus < 4))
                    %     if ((Jplus >= 1.95) && (Jplus < 2.25))
                    step = 0.01;
                    %     elseif Jplus>2.5
                    %         step=0.1;
                else
                    step = STEP;
                end
            else
                step=0.1;
            end
            
        end
        
        
        
        %----------------------------------------------------------------
        %----------------------------------------------------------------
        
        %-------------
        % FILE INFO
        %-------------
        
        function Files=File_Info_MFT(nn,p,flg)
            
            paramsfile=flg.paramsfile;
            
            extra='';
            load(paramsfile,'Network','Opt','Jplus','extra');
            
            if flg.nivsj_mode==0
                if flg.spur_mode==0
                    kind='spont';
                elseif flg.spur_mode==1
                    kind='spur';
                end
            elseif flg.nivsj_mode==1
                if flg.spur_mode==0
                    kind='sel';
                end
            end
            if flg.all_mode==1
                kind='all';
            end
            %
            % if flg.nivsj_mode==0 && flg.spur_mode==0
            %     kind='spont';
            % elseif flg.nivsj_mode==1 && flg.spur_mode==0
            %     kind='sel';
            % elseif flg.nivsj_mode==0 && flg.spur_mode==1
            %     kind='spur';
            % end
            if flg.stim_mode==1
                load(paramsfile,'StimHeight');
                Temp=[];
                if exist('StimHeight','var')
                    Temp=sprintf('%0.02g',StimHeight);
                    Temp=strrep(Temp,'.','_');
                end
                if flg.nivsj_mode==0
                    kind=['stim_spont_niExt' Temp];
                elseif flg.nivsj_mode==1
                    kind=['stim_sel' Temp];
                end
            end
            
            currentfolder=fullfile('DATA','Results',sprintf('%s',Opt),'MFT');
            paperfolder=fullfile(currentfolder,'figs');
            % if strcmp(Network.syn,'DoubleExp') || strcmp(Network.syn,'SingleExp')
            %     currentfolder=fullfile(currentfolder,'BS');
            %     extra='_BS';
            % end
            if ~exist(currentfolder,'dir')
                mkdir(currentfolder);
            end
            if ~exist(paperfolder,'dir')
                mkdir(paperfolder);
            end
            
            
            % File for eigenvalues dump:
            file_eig=fullfile(currentfolder,sprintf('%s%s_eig_[%d]pop_[%d]_%s.mat',Opt,extra,nn,p,kind));
            % File to store firing rates for full dynamics
            ff = fullfile(currentfolder,sprintf('%s%s_dynamics_[%d]_[%d]_%s.mat',Opt,extra,nn,p,kind));
            % File with dynamics sample
            file_log = fullfile(currentfolder,sprintf('%s%s_dynamics_log_[%d]_[%d]_%s.mat',Opt,extra,nn,p,kind));
            % file to store fixed points
            file_name=fullfile(currentfolder,sprintf('%s%s_nivsj_[%d]_[%d]_%s.mat',Opt,extra,nn,p,kind));
            % File to plot:
            % fixed point pic
            FPpic=fullfile(paperfolder,sprintf('%s%s_FixedPoints_[%d]_[%d]_%s.pdf',Opt,extra,nn,p,kind));
            % eigenvalue pic
            eigpic=fullfile(paperfolder,sprintf('%s%s_eig_[%d]_[%d]_%s.pdf',Opt,extra,nn,p,kind));
            % file to save stable states
            filesave=fullfile(currentfolder,sprintf('%s%s_FP_data[p%d]_%s.mat',Opt,extra,p,kind));
            % stable fixed point pic
            FPstablepic=fullfile(paperfolder,sprintf('%s%s_StableConfig[%d]%s.pdf',Opt,extra,p,kind));
            Allstablepic=fullfile(paperfolder,sprintf('%s%s_AllStable[%d]%s.pdf',Opt,extra,p,kind));
            % unstable fixed point pic
            FPunstablepic=fullfile(paperfolder,sprintf('%s%s_UnstableConfig[%d]%s.pdf',Opt,extra,p,kind));
            % active clusters vs Jplus, heat map pic
            multistable=fullfile(paperfolder,sprintf('%s%s_ActiveClusters[%d]%s.pdf',Opt,extra,p,kind));
            % LOG FILE
            log_file=(fullfile(currentfolder,sprintf('Log_%s%s_[p%d]%s.txt',Opt,extra,p,kind)));
            % file for effective trans fun mode
            file_m=fullfile(currentfolder,sprintf('%s%s_eff_vs_mu_[%d]_[%d]_%s.mat',Opt,extra,nn,p,kind));
            file_f=fullfile(currentfolder,sprintf('%s%s_eff_func_[%d]_[%d]_%s.mat',Opt,extra,nn,p,kind));
            file_bif=fullfile(currentfolder,sprintf('%s%s_eff_bif_[%d]_[%d]_%s.mat',Opt,extra,nn,p,kind));
            % amit-mascaro
            file_amitmascaro=fullfile(paperfolder,sprintf('%s%s_[Jp%g][p%d]_%s',Opt,extra,Jplus,p,kind));
            file_phi_eff=fullfile(paperfolder,sprintf('%s%s_phi_eff_[Jp%g][p%d]_%s',Opt,extra,Jplus,p,kind));
            
            % fig for paper
            % extra1=[];
            % % cue mode
            % if flg.cue_mode
            %     for i=1:numel(Network.cue.cue_option)
            %         extra1=[extra1 '_' Network.cue.cue_option{i} '' Network.cue.cue_stat{i}(1) '_'];
            %     end
            % end
            paperenergy=fullfile(paperfolder,sprintf('%s%s[p%d]energy',Opt,extra,p));
            paperfig=fullfile(paperfolder,sprintf('%s%s_PaperFig[p%d]_%s',Opt,extra,p,kind));
            paperstates=fullfile(paperfolder,sprintf('%s%s_PaperStates[p%d]_%s',Opt,extra,p,kind));
            
            Files.paramsfile=paramsfile;
            Files.file_eig=file_eig;
            Files.ff=ff;
            Files.file_log=file_log;
            Files.FP=file_name;
            Files.FPpic=FPpic;
            Files.FPstablepic=FPstablepic;
            Files.Allstablepic=Allstablepic;
            Files.FPunstablepic=FPunstablepic;
            Files.multistable=multistable;
            Files.eigpic=eigpic;
            Files.stable=filesave;
            Files.Log=log_file;
            Files.m=file_m;
            Files.f=file_f;
            Files.bif=file_bif;
            Files.rates_old=fullfile('DATA','rates_old.mat');
            Files.file_amitmascaro=file_amitmascaro;
            Files.file_phi_eff=file_phi_eff;
            
            % paper
            Files.paperfig=paperfig;
            Files.paperstates=paperstates;
            Files.paperenergy=paperenergy;
            
        end
        %----------------------------------------------------------------
        %----------------------------------------------------------------
        
        
        function result=CheckRates(ni,n_hi)
            
            result=0;
            TolRate=1e-4;
            p=numel(ni)-4;
            MaxRate=max(ni(1:p));
            MaxPops=sum(abs(ni(1:p)-MaxRate)<TolRate); % # of pops with high rate among first n_hi
            if MaxPops==n_hi
                result=1;
            end
        end
        
        
        
        %----------------------------------------------------------------
        %----------------------------------------------------------------
        
        function [QUENCH,quench_pops,nonquench_pops,GaussLimit]=fun_gausscheck(Params,n)
            
            QUENCH=0; quench_pops=[]; nonquench_pops=[];
            if any(strcmp(Params.Network,'cue'))
                if any(strcmp(Params.Network.cue.cue_stat,'gaussian')) && n>4
                    QUENCH=1;
                    quench_pops=Params.quench_pops;
                    nonquench_pops=setxor(1:n,quench_pops);
                end
            end
            GaussLimit=[10*ones(1,n-3),5*ones(1,3)]; % integration limit for eff transf function with gaussian spread, inh is smaller than exc for numerical reasons
            % if any(strcmp(fieldnames(Params),'Network'))
            %     if any(strcmp(fieldnames(Params.Network),'cue'))
            %     if any(strcmp(fieldnames(Params.Network),'cue_stat'))
            % %         if strcmp(Params.Network.cuemode,'gaussian') && n>2
            %         if ~isempty(strfind(Params.Network.cue_stat,'gaussian')) && n>4
            %             % if n=2 override gaussian cue (only used for threshold
            %             % computations)
            %             QUENCH=1;
            %             quench_pops=Params.quench_pops;
            %             nonquench_pops=setxor(1:n,quench_pops);
            %         end
            %     end
            % end
            % if any(strcmp(fieldnames(Params),'Network'))
            %     if any(strcmp(fieldnames(Params.Network),'cuemode'))
            % %         if strcmp(Params.Network.cuemode,'gaussian') && n>2
            %         if ~isempty(strfind(Params.Network.cuemode,'inh_gaussian')) && Params.Network.cue_value>0
            %             % if n=2 override gaussian cue (only used for threshold
            %             % computations)
            %             QUENCH=1;
            %             quench_pops=Params.quench_pops;
            %             nonquench_pops=setxor(1:n,quench_pops);
            %             GaussLimit=5; % integration limit for eff transf function with gaussian spread
            %         end
            %     end
            % end
        end
        %----------------------------------------------------------------
        %----------------------------------------------------------------
        
        
        % Mu
        %
        % from GLC C_overlap_thesis
        % Luca Mazzucato February 2014
        
        
        
        function [Mu,flg]=MU(ni,Params,flg)
            
            
            nn=numel(ni);
            A=Params.A;
            Mu_ext=Params.Mu_ext;
            Mu=zeros(nn,1);
            
            for i=1:nn
                Mu(i) = Mu_ext(i);
                for j=1:nn
                    Mu(i)=Mu(i)+ A(i,j)*ni(j);
                end
                %     if (Mu(i) == 0)
                %         fprintf('\n Allocation failure in routine MU_C\n');
                %         return_value = 0;
                %     end
            end
            
        end
        
        %----------------------------------------------------------------
        %----------------------------------------------------------------
        
        function out1=nerf(w)
            %
            a1=  -1.26551223;
            a2=  1.00002368;
            a3=  .37409196;
            a4=  .09678418;
            a5=  -.18628806;
            a6=  .27886807;
            a7=  -1.13520398;
            a8=  1.48851587;
            a9=  -.82215223;
            a10= .17087277;
            
            z = abs(w);
            t = 1./(1. + 0.5 * z);
            at=a1+t.*(a2+t.*(a3+t.*(a4+t.*(a5+t.*(a6+t.*(a7+t.*(a8+t.*(a9+t.*a10))))))));
            ef=t.*exp(at);
            out1 = 2.*exp(w.*w)-ef;
            if w<0
                out1 = ef;
            end
            %
            % %
            % % note that using either 1) nerf or 2) exp(w^2)(1+erf(w)) for w<0 screws up the transfer function
        end
        
        
        %----------------------------------------------------------------
        %----------------------------------------------------------------
        
        % Sigma
        %
        % from GLC C_overlap_thesis
        % Luca Mazzucato February 2014
        
        
        
        function [Sigma2,flg]=SIGMA2(ni,Params,flg)
            
            % global return_value
            
            % nn=numel(ni);
            % A=Params.A;
            B=Params.B;
            Sigma_ext=Params.Sigma_ext;
            % H=Params.H;
            % Tau=Params.Tau;
            % Theta=Params.Theta;
            
            % Sigma2=zeros(nn,1);
            
            Sigma2=B*ni'+Sigma_ext';
            
            if any(Sigma2<0.0)
                
                fprintf('>>> Error: Sigma[%d]<0\n',find(Sigma2<0));
                %         fprintf('\n\n>>> Error: Setting Sigma[%d] to 0.1*rand\n',i);
                %         Sigma2(i)=0.1*rand(1);
                flg.return_value = 0;
                return;
            end
            % end
        end
        
        
        
        
        
        
        %----------------------------------------------------------------
        %----------------------------------------------------------------
        
        % transfer function
        % INPUT ni   = firing rate row vector
        %       Params = input parameters
        %       Theta (optional): if absent, thresholds are set from Params.Theta
        %                         if present, thresholds are set equal to Theta
        % OUTPUT Risp=firing rate vector;
        
        %
        % from GLC C_overlap_thesis
        % Luca Mazzucato February 2014
        % Luca Mazzucato May 2020
        
        function [Risp,flg]=RISP(ni,Params,flg)
            
            warning('off','all');
            
            % global return_value
            
            n=numel(ni); % # of active populations
            
            % option for quenched noise
            QUENCH=0;
            if flg.cue_mode
                [QUENCH,quench_pops,nonquench_pops,GaussLimit]=auxMFT.fun_gausscheck(Params,n); % check which pops should have gaussian-smeared transfer function
            end
            
            % parameters
            H=Params.H;
            Tau=Params.Tau;
            tau_arp=Params.tau_arp;
            Theta=Params.Theta;
            
            Risp=zeros(n,1);
            % k = 0.9;
            % Brunel Sergi
            Tausyn=Params.Tausyn;
            BS=zeros(n,1);
            if Params.BS==1
                a=1.032626576115609; %-zeta(1/2)/sqrt(2); %
                BS=a*sqrt(Tausyn./Tau);
            end
            
            [Mu,flg]=auxMFT.MU(ni,Params,flg);         % calcola le Mu(ni)
            
            [Sigma,flg]=auxMFT.SIGMA2(ni,Params,flg);      % calcola le Sigma quadro(ni)
            if(flg.return_value == 0)
                fprintf('\n>>> Error in routine SIGMA() called by RISP()\n');
                return;
            end
            if any(Sigma<0)
                fprintf('\n SIGMA<0!!\n ');
            end
            if QUENCH==0
                % usual transfer function
                for i=1:n
                    sigma = sqrt(Sigma(i));
                    thres = (Theta(i)-Mu(i))/sigma+BS(i);
                    %thres = (Theta(i)+(k*ni(i))-Mu(i))/sigma;% Calcium effect
                    reset = (H(i)-Mu(i))/sigma+BS(i);
                    % NOTE: using Qsimp or quad, the integral gives the wrong values using nerf
                    % for w<0, so we need to use Maurizio Mattia's trick to get the correct
                    % numerical result. However, it takes 20x w.r.t. the mex file.
                    %     integ = Qsimp(@auxPhiExp, reset, thres);
                    %     integ = quad(@auxPhiExp, reset, thres);
                    %     integ=1.772453851*integ;
                    % USE MEX FILE TO EVALUATE INTEGRAL (already includes sqrt(pi) factor)
                    % NOTE: the MEX file integrates nerf both for w>=0 and w<0, giving in
                    % both cases the correct result without Mattia's trick!
                    integ=aux.IntAuxPhi_vec(reset,thres);
                    if (integ == -1)
                        fprintf('\n>>> Error in routine qsimp() called by RISP()\n');
                        flg.return_value = 0;
                        return;
                    end
                    
                    inv_risp = (tau_arp + Tau(i)*integ);
                    Risp(i) = 1./inv_risp;
                end
            elseif QUENCH==1
                % quenched gaussian noise on clustered populations only
                for i=quench_pops
                    sigma = sqrt(Sigma(i));
                    funMu_extZ=Params.Mu_extZ{i};
                    MuZ=@(z)(Mu(i)+funMu_extZ(z)); % input current with gaussian quenched noise z
                    thres =@(z) (Theta(i)-MuZ(z))/sigma+BS(i);
                    %thres = (Theta(i)+(k*ni(i))-Mu(i))/sigma;% Calcium effect
                    reset =@(z) (H(i)-MuZ(z))/sigma+BS(i);
                    integfun=@(z)(exp(-z.^2/2)/sqrt(2*pi))./(tau_arp + Tau(i)*aux.IntAuxPhi_vec(reset(z),thres(z)));
                    integ=integral(integfun,-GaussLimit(i),GaussLimit(i));
                    if (integ < 0)
                        fprintf('\n>>> Error in integral of IntAuxPhi_vec called by RISP_DYN()\n');
                        flg.return_value = 0;
                        return;
                    end
                    Risp(i)=integ;
                end
                % background and inhibitory populations without quenched noise
                for i=nonquench_pops
                    sigma = sqrt(Sigma(i));
                    thres = (Theta(i)-Mu(i))/sigma+BS(i);
                    %thres = (Theta(i)+(k*ni(i))-Mu(i))/sigma;% Calcium effect
                    reset = (H(i)-Mu(i))/sigma+BS(i);
                    integ=aux.IntAuxPhi_vec(reset,thres);
                    if (integ == -1)
                        fprintf('\n>>> Error in routine qsimp() called by RISP()\n');
                        flg.return_value = 0;
                        return;
                    end
                    inv_risp = (tau_arp + Tau(i)*integ);
                    Risp(i) = 1./inv_risp;
                end
            end
        end
        
        
        %----------------------------------------------------------------
        %----------------------------------------------------------------
        
        function [Risp_dyn, Mu_dyn, Sigma_dyn, flg]=RISP_DYN(ni,Mu_dyn,Sigma_dyn,Contrast,Params,delta_t,dyn_value,flg)
            
            
            n=numel(ni);
            QUENCH=0;
            if flg.cue_mode
                [QUENCH,quench_pops,nonquench_pops,GaussLimit]=auxMFT.fun_gausscheck(Params,n); % check which pops should have gaussian-smeared transfer function
            end
            
            H=Params.H;
            Tau=Params.Tau;
            tau_arp=Params.tau_arp;
            Theta=Params.Theta;
            Risp_dyn=zeros(1,n);
            % Brunel Sergi
            Tausyn=Params.Tausyn;
            BS=zeros(n,1);
            if Params.BS==1
                a=1.032626576115609; %-zeta(1/2)/sqrt(2); %
                BS=a*sqrt(Tausyn./Tau);
            end
            %
            if dyn_value==1
                % first step: initial values
                [Mu_dyn,flg]=auxMFT.MU(ni,Params,flg); % compute Mu(i)
                [Sigma_dyn,flg]=auxMFT.SIGMA2(ni,Params,flg);        % compute Sigma(i) (= Sigma^2(ni))
            end
            
            [Mu_dyn,flg]=auxMFT.MU_DYN(ni,Mu_dyn,Contrast,Params,delta_t,flg);
            [Sigma_dyn,flg]=auxMFT.SIGMA_DYN(ni,Sigma_dyn,Params,delta_t,flg);
            if (flg.return_value == 0)
                fprintf('\n>>> Error in routine SIGMA_DYN() called by RISP_DYN()\n');
                return;
            end
            if any(Sigma_dyn<0)
                fprintf('\n SIGMA<0!!\n ');
            end
            
            if QUENCH==0
                for i=1:n
                    sigma = sqrt(Sigma_dyn(i));
                    reset=((H(i)-Mu_dyn(i))/sigma)+BS(i);
                    thres=((Theta(i)-Mu_dyn(i))/sigma)+BS(i);
                    integ = aux.IntAuxPhi_vec(reset , thres);
                    if (integ == -1)
                        fprintf('\n>>> Error in routine IntAuxPhi_vec called by RISP_DYN()\n');
                        flg.return_value = 0;
                        return;
                    end
                    inv_risp = (tau_arp + Tau(i)*integ);
                    Risp_dyn(i) = 1./inv_risp;
                end
            elseif QUENCH==1
                % quenched gaussian noise
                for i=quench_pops
                    sigma = sqrt(Sigma_dyn(i));
                    funMu_extZ=Params.Mu_extZ{i};
                    MuZ=@(z)(Mu_dyn(i)+funMu_extZ(z)); % input current with gaussian quenched noise z
                    thres =@(z) (Theta(i)-MuZ(z))/sigma+BS(i);
                    %thres = (Theta(i)+(k*ni(i))-Mu(i))/sigma;% Calcium effect
                    reset =@(z) (H(i)-MuZ(z))/sigma+BS(i);
                    integfun=@(z)(exp(-z.^2/2)/sqrt(2*pi))./(tau_arp + Tau(i)*aux.IntAuxPhi_vec(reset(z),thres(z)));
                    integ=integral(integfun,-GaussLimit(i),GaussLimit(i));
                    if (integ < 0)
                        fprintf('\n>>> Error in integral of IntAuxPhi_vec called by RISP_DYN()\n');
                        flg.return_value = 0;
                        return;
                    end
                    Risp_dyn(i)=integ;
                end
                % background and inhibitory populations without quenched noise
                for i=nonquench_pops
                    sigma = sqrt(Sigma_dyn(i));
                    thres = (Theta(i)-Mu_dyn(i))/sigma+BS(i);
                    %thres = (Theta(i)+(k*ni(i))-Mu(i))/sigma;% Calcium effect
                    reset = (H(i)-Mu_dyn(i))/sigma+BS(i);
                    integ=aux.IntAuxPhi_vec(reset,thres);
                    if (integ == -1)
                        fprintf('\n>>> Error in routine IntAuxPhi_vec() called by RISP_DYN()\n');
                        flg.return_value = 0;
                        return;
                    end
                    inv_risp = (tau_arp + Tau(i)*integ);
                    Risp_dyn(i) = 1./inv_risp;
                end
            end
            
        end
        %-------------------------------------------------------------
        
        function [Sigma_dyn,flg]=SIGMA_DYN(ni,Sigma_dyn,Params,delta_t,flg)         % compute Sigma_dyn^2[1,...,4]
            
            % global flg.return_value flg.syn_dyn
            
            n=numel(ni);
            [Sigma,flg]=auxMFT.SIGMA2(ni,Params,flg);        % compute Sigma(i) (= Sigma^2(ni))
            if flg.syn_dyn
                taus=Params.Tausyn;
            else
                taus=Params.Tau;
            end
            for i=1:n
                Sigma_dyn(i) = (1.-2*delta_t/taus(i))*Sigma_dyn(i)+(2*delta_t/taus(i))*Sigma(i);
            end
            % if (flg.return_value == 0)
            %     fprintf('\n>>> Error in routine SIGMA2() called by SIGMA_DYN()\n');
            %     return;
            % end
            %
        end
        %-------------------------------------------------------------
        
        function [Mu_dyn,flg]=MU_DYN(ni,Mu_dyn,Contrast,Params,delta_t,flg)            % compute Mu_dyn[1,...,4]
            
            % global flg.return_value flg.syn_dyn
            
            n=numel(ni);
            Mu_ext=Params.Mu_ext;
            [Mu,flg]=auxMFT.MU(ni,Params,flg); % compute Mu(i)
            if flg.syn_dyn
                taus=Params.Tausyn;
            else
                taus=Params.Tau;
            end
            for i=1:n
                Mu_dyn(i) = (1.-delta_t/taus(i))*Mu_dyn(i)+(delta_t/taus(i))*(Mu(i)+Contrast(i)*Mu_ext(i));
            end
            % if (flg.return_value == 0)
            %     fprintf('\n>>> Error in routine MU() called by MU_DYN()');
            %     return;
            % end
        end
        
        
        
        
        
        %----------------------------------------------------------------
        %----------------------------------------------------------------
        
        % transfer function for quenched noise computation as a function of quench
        % gaussian variable z
        %
        % INPUT ni   = firing rate row vector
        %       Params = input parameters
        %       Theta (optional): if absent, thresholds are set from Params.Theta
        %                         if present, thresholds are set equal to Theta
        %       i_pop    = integer corresponding to quenched current population: 1<=i_pop<=n-1
        % OUTPUT Risp=function of z;
        
        %
        % from GLC C_overlap_thesis
        % Luca Mazzucato February 2014
        % Luca Mazzucato May 2020
        
        function [Risp,flg]=RISP_fun_quenched(ni,Params,i_pop,flg)
            
            warning('off','all');
            
            % first two arguments are fixed
            
            % ni=[3.6; 5.2];
            % global return_value
            
            % option for quenched noise
            % option for quenched noise
            QUENCH=0;
            n=numel(ni);
            if flg.cue_mode
                [QUENCH,quench_pops,nonquench_pops,GaussLimit]=auxMFT.fun_gausscheck(Params,n); % check which pops should have gaussian-smeared transfer function
            end
            
            % parameters
            H=Params.H;
            Tau=Params.Tau;
            tau_arp=Params.tau_arp;
            Theta=Params.Theta;
            
            nn=numel(ni); % # of active populations
            % Risp=zeros(nn,1);
            % k = 0.9;
            % Brunel Sergi
            Tausyn=Params.Tausyn;
            BS=zeros(nn,1);
            if Params.BS==1
                a=1.032626576115609; %-zeta(1/2)/sqrt(2); %
                BS=a*sqrt(Tausyn./Tau);
            end
            
            [Mu,flg]=auxMFT.MU(ni,Params,flg);         % calcola le Mu(ni)
            if(flg.return_value == 0)
                fprintf('\n>>> Error in routine MU() called by RISP()\n');
                return;
            end
            [Sigma,flg]=auxMFT.SIGMA2(ni,Params,flg);      % calcola le Sigma quadro(ni)
            if(flg.return_value == 0)
                fprintf('\n>>> Error in routine SIGMA() called by RISP()\n');
                return;
            end
            if any(Sigma<0)
                fprintf('\n SIGMA<0!!\n ');
            end
            if QUENCH==1
                i=i_pop;
                % quenched gaussian noise on clustered populations only
                sigma = sqrt(Sigma(i));
                funMu_extZ=Params.Mu_extZ{i};
                MuZ=@(z)(Mu(i)+funMu_extZ(z)); % input current with gaussian quenched noise z
                thres =@(z) (Theta(i)-MuZ(z))/sigma+BS(i);
                %thres = (Theta(i)+(k*ni(i))-Mu(i))/sigma;% Calcium effect
                reset =@(z) (H(i)-MuZ(z))/sigma+BS(i);
                Risp=@(z)1./(tau_arp + Tau(i)*aux.IntAuxPhi_vec(reset(z),thres(z)));
                % background and inhibitory populations without quenched noise
            else
                fprintf('\n error in RISP_fun_quenched, inconsistent pop index i \n');
                flg.return_value = 0;
                return;
            end
        end
        %----------------------------------------------------------------
        %----------------------------------------------------------------
        
        function Summary(paramsfile)
            
            % global overlap
            load(paramsfile);
            
            flag = 0;
            fprintf('\n--- Summary of Parameters:\n\n');
            fprintf('    Total # of Pop.s   : %d   # of clust.s: %d    Coding: %g\n',nn,p,f);
            fprintf('    Synapses           : J+ %0.03g,  J- %0.03g,  Jee %0.03g,   Jei:[ ',Jplus,Jminus,Jee);
            fprintf('%0.03g ',Jei);
            fprintf(']\n                         Jii: [');
            fprintf('%0.03g ',Jii);
            fprintf(']\n                         Jie: [');
            fprintf('%0.03g ',Jie);
            fprintf(']   Jo %g\n',Jzero);
            fprintf('    Learning Parameters: gamma %g   rho %g   delta^2: %g\n',gam,ro,delta);
            fprintf('    Time Costants      : tau_arp %g   tau_e %g  tau_i %g\n',tau_arp,tau_e,tau_i);
            fprintf('    Spont. Act.        : ni_e %g     ni_i   %g\n',ni_e);
            fprintf(' %0.03g ',ni_i);
            fprintf('\n                         ni_ext_e %g     ni_ext_i ',ni_ext_e);
            fprintf('%0.03g ',ni_ext_i);
            fprintf('\n    Thresholds         : theta_e=%0.03g   theta_i=[',theta_e);
            fprintf('%0.03g ',theta_i);
            fprintf(']\n');
            
            %if      (overlap && flag == 0)  flag = 1;
            % if overlap && (flag == 0)
            %     fprintf('\n--- bg_coding: %g\n',d(nn-3));
            % else
            if ((1.-f*p)<0.)
                error('\n>>> ERROR: bg_coding < 0 (%g). Exit to program \n\n',1.-f*p);
            else
                fprintf('\n--- bg_coding: %g\n',1.-f*p);
            end
            % end
            
            % if      (overlap == 1)
            %     fprintf('\n--- PATTERNS IN OVERLAP (General Case)\n\n');
            % elseif (overlap == 2)
            %     fprintf('\n--- Low Loading Overlap Case\n\n');
            % elseif (overlap == 3)
            %     fprintf('\n--- 3 Patterns in OVERLAP \n\n');
            % else
            %     fprintf('\n--- NO OVERLAP \n\n');
            % end
            
        end
        %----------------------------------------------------------------
        %----------------------------------------------------------------
        
        % /**********************************************************\
        % *                  (analytical) usrfun		          *
        % \**********************************************************/
        
        
        function [Alpha, Bet, flg]=usrfun(ni,Params,flg)    % routine di supporto a mnewt()
            
            % global return_value
            
            n=numel(ni);
            Bet=zeros(1,n);
            Alpha=zeros(n);
            
            % option for quenched noise
            QUENCH=0;
            if flg.cue_mode
                [QUENCH,quench_pops,nonquench_pops,GaussLimit]=auxMFT.fun_gausscheck(Params,n); % check which pops should have gaussian-smeared transfer function
            end
            
            [Mu,flg]=auxMFT.MU(ni,Params,flg);
            [Sigma,flg]=auxMFT.SIGMA2(ni,Params,flg);
            if (flg.return_value == 0)
                fprintf('\n>>> Error in routine SIGMA2() called by usrfun()\n');
                return;
            end
            [Risp,flg]=auxMFT.RISP(ni,Params,flg);       % computing Mu, Sigma and transfer function
            %
            if (flg.return_value == 0)
                fprintf('\n>>> Error in routine RISP() called by usrfun()\n');
                return;
            end
            % parameters
            A=Params.A;
            B=Params.B;
            % Sigma_ext=Params.Sigma_ext;
            H=Params.H;
            Tau=Params.Tau;
            Theta=Params.Theta;
            % Brunel Sergi
            Tausyn=Params.Tausyn;
            BS=zeros(n,1);
            if Params.BS==1
                a=1.032626576115609; %-zeta(1/2)/sqrt(2); %
                BS=a*sqrt(Tausyn./Tau);
            end
            sigma = sqrt(Sigma);
            
            
            if QUENCH==0
                % usual transfer function
                X = zeros(n);
                Y = zeros(n);
                for i=1:n
                    aa = (H(i)-Mu(i))/sigma(i)+BS(i);
                    bb = (Theta(i)-Mu(i))/sigma(i)+BS(i);
                    
                    for j=1:n
                        X(i,j) = A(i,j) + (0.5*(aa-BS(i))*B(i,j)/sigma(i));
                        Y(i,j) = A(i,j) + (0.5*(bb-BS(i))*B(i,j)/sigma(i));
                        
                        %         // d(phi(i))/d(Ni(j)):
                        Alpha(i,j) =  ((Tau(i)*Risp(i)*Risp(i))/sigma(i))*(auxMFT.phi(bb)*Y(i,j) - auxMFT.phi(aa)*X(i,j));
                    end
                    Alpha(i,i) =Alpha(i,i)- 1.;
                    Bet(i)       = -Risp(i)+ni(i);
                end
            elseif QUENCH==1
                %------------------
                % transfer function with quenched noise on first n-2 populations
                for i=quench_pops
                    funMu_extZ=Params.Mu_extZ{i};
                    MuZ=@(z)(Mu(i)+funMu_extZ(z)); % input current with gaussian quenched noise z
                    aaZ =@(z) (H(i)-MuZ(z))/sigma(i)+BS(i);
                    bbZ =@(z) (Theta(i)-MuZ(z))/sigma(i)+BS(i);
                    [Risp_funZ,flg]=auxMFT.RISP_fun_quenched(ni,Params,i,flg);
                    phi_b=@(z)auxMFT.phi(bbZ(z));
                    phi_a=@(z)auxMFT.phi(aaZ(z));
                    for j=1:n
                        X =@(z)( A(i,j) + (0.5*(aaZ(z)-BS(i))*B(i,j)/sigma(i)) );
                        Y =@(z)( A(i,j) + (0.5*(bbZ(z)-BS(i))*B(i,j)/sigma(i)) );
                        %         // d(phi(i))/d(Ni(j)):
                        integrand_Alpha=@(z)(exp(-z.^2/2)/sqrt(2*pi)).*((Tau(i).*Risp_funZ(z).*Risp_funZ(z))/sigma(i)).*(phi_b(z).*Y(z) - phi_a(z).*X(z));
                        Alpha(i,j) = integral(integrand_Alpha,-GaussLimit(i),GaussLimit(i))  ; % normalized gaussian integral with infty replaced by 5 (sufficient for total prob to integrate to 1)
                    end
                    Alpha(i,i) =Alpha(i,i)- 1.;
                    Bet(i)       = -Risp(i)+ni(i);
                end
                %-----------
                % background (n-1) and inhibitory (n) populations not quenched
                if (flg.return_value == 0)
                    fprintf('\n>>> Error in routine RISP() called by usrfun()\n');
                    return;
                end
                for i=nonquench_pops
                    % usual transfer function
                    aa = (H(i)-Mu(i))/sigma(i)+BS(i);
                    bb = (Theta(i)-Mu(i))/sigma(i)+BS(i);
                    for j=1:n
                        X = A(i,j) + (0.5*(aa-BS(i))*B(i,j)/sigma(i));
                        Y = A(i,j) + (0.5*(bb-BS(i))*B(i,j)/sigma(i));
                        %         // d(phi(i))/d(Ni(j)):
                        Alpha(i,j) =  ((Tau(i)*Risp(i)*Risp(i))/sigma(i))*(auxMFT.phi(bb)*Y - auxMFT.phi(aa)*X);
                    end
                    Alpha(i,i) =Alpha(i,i)- 1.;
                    Bet(i)       = -Risp(i)+ni(i);
                end
            end
            
            
            % % /* Debug code:
            % fprintf(' \nMatrix:\n\n');
            % for i=1:n
            %     for j=1:n
            %         fprintf(' %g  ', Alpha(i,j));
            %     end
            % fprintf('\n');
            % end
            % fprintf('\n');
            
        end
        
        %-----------------------------------------------
        % AUXILIARY FUNCTIONS
        %-----------------------------------------------
        
        %-------------------------------------------------------------
        
        function out=phi(w)
            
            out=1.772453851*auxMFT.nerf(w);
            
            
        end
        
        
        %----------------------------------------------------------------
        %----------------------------------------------------------------
        
        %/***************************************************
        % *                                                 *
        % *                Test_Stability                   *
        % * USAGE: [res Real Imag]=Test_Stability(ni,Params)*
        % * res= 1 if fixed point is stable, 0 if it's   *
        % * unstable
        % * Real, Imag are real and imaginary parts of eigenvalues
        % *                                                 *
        % ***************************************************/
        
        function [res, Real_part, Imag_part, EigV, flg, varargout]=Test_Stability(ni,Params,flg)
            
            % global return_value nfocus screen mnewt_flag syn_dyn
            
            n=numel(ni);
            
            % initialization
            Imag_part = zeros(1,n);
            Real_part = zeros(1,n);
            matrix = zeros(n);
            phiprime = NaN(n,1);
            if (flg.screen)
                fprintf('    Testing stability of fixed point using eigenvalues...\n\n');
            end
            linresp=NaN(n);
            
            % % Creating Stability Matrix:
            [matrix,flg,phiprime,linresp]=auxMFT.Stability_Matrix(ni,matrix,Params,flg);
            if (flg.screen); fprintf('>>> Stability Matrix: \n'); for i=1:n; for j=1:n; fprintf(' %g',matrix(i,j));  end; fprintf('\n'); end; end
            if (flg.return_value == 0); fprintf('>>> Error in routine Stability_Matrix called by Test_Stability\n'); end
            varargout{1}=matrix;
            varargout{2}=phiprime;
            varargout{3}=linresp;
            
            [EigV lambda]=eig(matrix);
            Real_part=real(diag(lambda));
            Imag_part=imag(diag(lambda));
            if (flg.return_value == 0)
                fprintf('\n>>> Error in routine hqr() called by Test_Stability\n');
                flg.return_value = 1; %% Ignore and continue.
                res= -1;
            end
            
            %  % Testing stability:
            count = 0;
            if (flg.screen); fprintf('    Eigenvalues:\n'); end
            for i=1:n
                if (flg.screen); fprintf('    Real_part[%d] %9.6g        Imag_part[%d] %g\n',i,Real_part(i),i,Imag_part(i)); end
                if (Real_part(i) < 0.0)
                    count=count+1;
                end
            end
            if (count == n)
                if (flg.screen || flg.mnewt_flag == 1); fprintf('    Fixed point is stable\n\n'); end
                res=1;
            else
                if (flg.screen || flg.mnewt_flag == 1); fprintf('>>> Fixed point is NOT stable\n\n'); end
                res=0;
            end
            return;
            
        end
        
        
        %----------------------------------------------------------------
        %----------------------------------------------------------------
        
        %/*****************************************************
        % *                                                   *
        % *                  Stability_Matrix                 *
        % *                                                   *
        % * Compute stability matrix analytically             *
        % *                                                   *
        % *****************************************************/
        
        function [m,flg,varargout]=Stability_Matrix(ni,m,Params,flg)
            
            % global return_value DEBUG syn_dyn
            
            paramsfile=Params.Files.paramsfile;
            load(paramsfile,'tau_e'); % needed below
            n=numel(ni);
            
            filemm = 'Data/stabil.mat'; % saving current stability matrix
            QUENCH=0;
            if flg.cue_mode
                [QUENCH,quench_pops,nonquench_pops,GaussLimit]=auxMFT.fun_gausscheck(Params,n); % check which pops should have gaussian-smeared transfer function
            end
            
            % parameters
            A=Params.A;
            B=Params.B;
            H=Params.H;
            Tau=Params.Tau;
            Theta=Params.Theta;
            % computing Mu, Sigma and transfer function
            [Mu,flg]=auxMFT.MU(ni,Params,flg);
            [Sigma2,flg]=auxMFT.SIGMA2(ni,Params,flg);
            [Risp,flg]=auxMFT.RISP(ni,Params,flg);
            % Brunel Sergi
            Tausyn=Params.Tausyn;
            BS=zeros(n,1);
            if Params.BS==1
                a=1.032626576115609; %-zeta(1/2)/sqrt(2); %
                BS=a*sqrt(Tausyn./Tau);
            end
            %
            if (flg.return_value == 0);  fprintf('\nError in routine RISP() called by Stability_Matrix()\n');
                return;
            end
            
            sigma = sqrt(Sigma2);
            aa = (H-Mu)./sigma+BS;
            bb = (Theta-Mu)./sigma+BS;
            phiprime=zeros(n,1);
            if QUENCH==0
                for i=1:n
                    if bb(i)>26; bbi=26; else bbi=bb(i); end % avoid Inf setting bbi to largest value giving finite phi(bbi)
                    phiprime(i) = ((Tau(i)*Risp(i)*Risp(i))/sigma(i))*(auxMFT.phi(bbi) - auxMFT.phi(aa(i)));
                    for j=1:n
                        %
                        %         X(i,j) = A(i,j) + (0.5*(aa(i)-BS(i))*B(i,j)/sigma(i));
                        %         Y(i,j) = A(i,j) + (0.5*(bb(i)-BS(i))*B(i,j)/sigma(i));
                        
                        % debug: ------------------------------------------------------------
                        %         //fprintf('\n   delta[%d][%d]: %g  \n', i,j,(((Tau(i)*Risp(i)*Risp(i))/(2*Sigma(i)))*B(i,j)*(bb*phi(bb) - aa*phi(aa))));
                        % -------------------------------------------------------------------
                        
                        % matrix d(phi(i))/d(Ni[j]):
                        %         m(i,j) = ((Tau(i)*Risp(i)*Risp(i))/sigma(i))*(phi(bb(i))*    Y(i,j) - phi(aa(i))*X(i,j));
                        
                        % EXACT STABILITY MATRIX:
                        % AMIT-MASCARO 99 eq. (48)
                        m(i,j) = phiprime(i)*A(i,j);
                        
                    end
                    m(i,i) =m(i,i) - 1.;
                end
            elseif QUENCH==1
                %------------------
                % transfer function with quenched noise
                for i=quench_pops
                    funMu_extZ=Params.Mu_extZ{i};
                    MuZ=@(z)(Mu(i)+funMu_extZ(z)); % input current with gaussian quenched noise z
                    aaZ =@(z) (H(i)-MuZ(z))/sigma(i)+BS(i);
                    bbZ =@(z) (Theta(i)-MuZ(z))/sigma(i)+BS(i);
                    [Risp_funZ,flg]=auxMFT.RISP_fun_quenched(ni,Params,i,flg);
                    phi_b=@(z)auxMFT.phi(bbZ(z));
                    phi_a=@(z)auxMFT.phi(aaZ(z));
                    integrand_Alpha=@(z)(exp(-z.^2/2)/sqrt(2*pi)).*(Tau(i).*Risp_funZ(z).*Risp_funZ(z)/sigma(i)).*(phi_b(z) - phi_a(z));
                    phiprime(i) = integral(integrand_Alpha,-GaussLimit(i),GaussLimit(i))  ; % normalized gaussian integral with infty replaced by 5 (sufficient for total prob to integrate to 1)
                    for j=1:n
                        m(i,j)= phiprime(i)*A(i,j);
                    end
                    m(i,i) =m(i,i) - 1.;
                end
                
                %-----------
                for i=nonquench_pops
                    % usual transfer function
                    if bb(i)>26; bbi=26; else bbi=bb(i); end
                    phiprime(i) = ((Tau(i)*Risp(i)*Risp(i))/sigma(i))*(auxMFT.phi(bbi) - auxMFT.phi(aa(i)));
                    for j=1:n
                        m(i,j) = phiprime(i)*A(i,j);
                    end
                    m(i,i) =m(i,i) - 1.;
                end
            end
            % fprintf('m(1,1): %g -- m(1,2)=%g \n',m(1,1),m(1,2));
            if flg.syn_dyn
                taus=Tausyn;
            else
                taus=Tau;
            end
            for i=1:n
                for j=1:n
                    m(i,j) =m(i,j)* (tau_e/taus(i)); %% eigenvalues in units of tau_e
                end
            end
            
            if flg.DEBUG
                % /* Debug code: */
                fprintf(' \n--- Matrix:\n\n');
                for i=1:n
                    fprintf('   ');
                    for j=1:n
                        fprintf(' %g ', m(i,j));
                    end
                    fprintf('\n');
                end
                fprintf('\n');
            end %% flg.DEBUG
            
            linresp=inv(diag(1./phiprime)-A);
            varargout{1}=phiprime;
            varargout{2}=linresp;
            
            % return
            
        end
        
        %-----------------------------------------------
        %-----------------------------------------------
        
        % /*************************************************
        %  *                                               *
        %  *                    Learn                      *
        %  *                                               *
        %  *                                               *
        %  *************************************************/
        
        function Params_out=Learn(Params,flg)
            
            % global screen overlap nivsj_mode spur_mode stim_mode all_mode
            paramsfile=Params.Files.paramsfile;
            load(paramsfile);
            flag = 0;
            % % AB97 rule:
            % Jminus = (2-f*(p+Jplus))/(2-f*(p+1));
            % % B99 rule:
            % Jminus
            % if Jminus==1
            Jminus = 1.-gam*f*(Jplus-1.);
            % end
            % % OVL rule:
            % if (flg.overlap)
            %     corr = gam*f*(Jplus-Jzero);
            %     Jminus = Jzero-corr;
            % end
            % if strcmp(Opt,'DOIRON_good')
            %     Jminus=1;
            % end
            if Jminus<0
                error('ERROR in create_params: Jminus<0. Adjust gam or f\n');
            end
            
            % update Jminus in paramsfile called by the following functions
            save(paramsfile,'Jminus','-append');
            
            if (flg.screen)
                Summary(paramsfile);
                %     sleep(sleep_value);
            else
                %     if ((flag == 0) && flg.overlap)
                %         fprintf('\n--- Initializing Parameters ...\n');
                %         flag = 1;
                %     end
            end
            
            if (nn == 4)
                Params_out=aux.In_Param_4pop_InhType([theta_e theta_i],paramsfile,flg);
            else
                Params_out=aux.In_Parameters_InhType([theta_e theta_i],paramsfile,flg);  %% NO Overlap
            end
        end
        
        %-----------------------------------------------
        %-----------------------------------------------
        
        
        
        
        function [ni,flg]=mnewt(ntrial,ni,Params,flg)
            
            
            
            if flg.screen
                fprintf('    Computing fixed point by mnewt()...\n\n');
            end
            
            n=numel(ni);
            indx=zeros(1,n);
            
            for k=1:ntrial
                if flg.nfocus == 0
                    [Alpha, Bet, flg]=auxMFT.usrfun(ni,Params,flg);
                    if (flg.return_value == 0)
                        fprintf('\n>>> Error in routine usrfun() called by mnewt()\n');
                        return;
                    end
                    
                end
                % diagnostic  code
                for i=1:n
                    if (Bet(i) >= 1000)
                        fprintf('\n>>> Error in routine mnewt(): |Phi[%d] - Ni[%d]| is too big (>1000)\n',i,i);
                        flg.return_value = 0;
                        return;
                    end
                end
                % end diagnostic
                
                %  /*------------------------- SCREEN SECTION ----------------------------*/
                %  /*----------------- HAVE WE REACHED A FIXED POINT ? -------------------*/
                %  /*--------------------------- FIRST CASE  -----------------------------*/
                errf=0.0;
                errf =sum(abs(Bet));
                
                if (errf <= flg.TOLF)
                    %          If flg.screen flag is enabled show the output:
                    if (flg.screen  || flg.mnewt_flag)
                        fprintf('    Fixed point by mnewt(): \n');
                        fprintf('    ( mnewt() tolerance on Transfer Function: %g )\n\n',flg.TOLF);
                        auxMFT.Info_Newt(ni,Params,flg);
                    end
                    
                    return;
                end
                
                
                %-----------------------------------------
                % Matrix inversion using Matlab functions
                %-----------------------------------------
                % A*x=b-> x=A\b
                Bet=Alpha\Bet';
                Bet=Bet';
                %-----------------------------------------
                
                % /*------------------------- SCREEN SECTION ----------------------------*/
                % /*----------------- HAVE WE REACHED A FIXED POINT ? -------------------*/
                % /*-------------------------- SECOND CASE  -----------------------------*/
                errx=0.0;
                errx =sum(abs(Bet));
                ni =ni + Bet;
                
                if (errx <= flg.TOLX)
                    %        If flg.screen flag is enabled show the output:
                    if (flg.screen  || flg.mnewt_flag)
                        fprintf('    Fixed point by mnewt(): \n');
                        fprintf('    ( mnewt() tolerance on Fixed Point: %g )\n\n',flg.TOLX);
                        auxMFT.Info_Newt(ni,Params,flg);
                    end
                    
                    return;
                end
            end
            fprintf('\n\n>>> Too many steps (>%d) in routine mnewt\n',ntrial);
            flg.return_value = 0;
            return;
            
        end
        
        %-------------------------------------------------------------
        
        
        % /**********************************************************\
        %  *                  	 auxMFT.Info_Newt                        *
        %  *                                                        *
        %  *               To put informations to flg.screen            *
        % \**********************************************************/
        
        
        function Info_Newt(ni,Params,flg)
            
            % global overlap DETAILS DEBUG
            
            n=numel(ni);
            bg=n-1;
            in=n;
            paramsfile=Params.Files.paramsfile;
            load(paramsfile,'p_actual');
            % computing Mu, Sigma
            [Mu,flg]=auxMFT.MU(ni,Params,flg);
            [Sigma,flg]=auxMFT.SIGMA2(ni,Params,flg);
            
            fprintf('        [Ni:]      [Mu:]        [Sigma:]   \n');
            for i=1:n
                fprintf('    [%d] %3.5f    %4.5f    %4.5f\n',i,ni(i),Mu(i),sqrt(Sigma(i)));
            end
            fprintf('\n');
            
            
            if flg.DETAILS
                
                fprintf('    Contributions:\n\n');
                for i=1:n
                    fprintf('    Mu[%d]: ',i);
                    for j=1:nn+1
                        fprintf(' %1.5f ',Mu_perc(i,j));
                    end
                    fprintf('\n');
                end
                fprintf('\n');
                
                if flg.DEBUG
                    %  same percentage as currents
                    for i=1:n
                        fprintf('    Sigma[%d]: ',i);
                        for j=1:nn+1
                            fprintf(' %1.5f ',Mu_perc(i,j));
                        end
                        fprintf('\n');
                    end
                end
            end
            
        end
        
        %-------------------------------------------------------------
        
        % /*------------------------------------*
        % *       Active_bits(int i,int p)     *
        % *                                    *
        % * Returns the number of active bits  *
        % * of population 'i' (i.e. the number *
        % * of patterns to which pop. 'i' is   *
        % * selective).                        *
        % * p = total number of patterns.      *
        % *------------------------------------*/
        
        function n=Active_bits(i,p)
            
            n=0;
            k=1;
            for j=0:p-1
                
                if (i && k)
                    n=n+1;
                end
                k = bitsll(k,1);
            end
            
        end
        
        %---------------------------------------------------------------------
        %---------------------------------------------------------------------
        
        
        
        % /**********************************************\
        % *                                            *
        % *                goAB_simul                  *
        % * *
        % * Simulate dynamics of Amit-Brunel. Saves population activities in nn files (ni(i).*
        % * dyn), 'simul' in dir Data.            *
        % * Parameter 'int kmax' is the largest number of steps allowed: if goAB_simul is called with kmax = 0 it stops at the fixed point, otherwise it continues indefinitely.                  *
        % * Every tot steps verifies with mnewt() if it's close to a fixed point and verifies stability. The fixed point is not modified by mnewt and the simulations restarts from it.
        % * Parameter 'int contr' allows contrast (default=1, disallowed). To allow it, select -c option *
        % *    GLC (2000), LM (2020)                                         *
        % \**********************************************/
        
        function [ni_out,flg]=goAB_simul(ni,kmax,Params,flg)    % Amit-Brunel dynamics simulation
            
            % global return_value contr screen NTRIAL TOLF TOLX syn_dyn
            
            % load parameters
            paramsfile=Params.Files.paramsfile;
            load(paramsfile);
            
            nn=numel(ni);
            if nn<=50
                MAXPOP=nn;
            else
                MAXPOP=50;
            end
            CONTRAST=5e-1;
            TOTAL_COUNT=5000;
            STEP_DYN_SIMUL=1e-4;  % 1e-4 Integration step in sec
            k = 1; % loop count
            tstart = 1;
            tSTOP = 501;
            tot_count = TOTAL_COUNT;  % To check for fixed points
            
            % Allocate and Initialize:
            Contrast = zeros(1,nn);
            Ni_temp = zeros(1,nn);
            delta_t = STEP_DYN_SIMUL;
            
            ff=Params.Files.ff;
            file_log=Params.Files.file_log;
            % contains a sample of firing rates during a run of the firing rate dynamics.';
            
            % report parameters
            fprintf('--- goAB_simul() parameters:\n');
            fprintf('    delta_t       : %g ms\n',delta_t*1000);
            if (kmax ~= 0)
                fprintf('    kmax          : %d\n',kmax);
            end
            fprintf('    nn            : %d\n',nn);
            fprintf('    f             : %g\n',f);
            fprintf('    p             : %d\n',p);
            fprintf('    J+            : %g;  J-  :%g\n',Jplus,Jminus);
            if (flg.contr == 1)
                fprintf('    Contrast      : %g\n',CONTRAST);
            end
            fprintf('    starting point:');
            for j=1:nn
                fprintf(' %g ',ni(j));
            end
            fprintf('\n');
            
            % store variables for saving to file
            time_steps=[];
            rates=[];
            sample_rates=[];
            sample_time=[];
            ni_out=[];
            
            % Starting Dynamics:
            %---------------------------------------------------------------------
            while(1)
                
                % Handle contrasts:
                if (flg.contr == 1)
                    Contrast=zeros(1,nn);
                    if ( k >= tstart) && (k<= tSTOP)
                        Contrast(1) = CONTRAST;
                    end
                end
                %-------------
                % START DYNAMICS
                %-------------
                if ( k == 1 )
                    dyn_value = 1;
                    % RISP_DYN does not use Mu_dyn and Sigma_dyn but iself generates initial
                    % values
                    Mu_dyn=[];
                    Sigma_dyn=[];
                else
                    % RISP_DYN accepts values at previous dynamical step
                    dyn_value = 0;
                end
                % dynamics:
                [Risp_dyn Mu_dyn Sigma_dyn,flg]=auxMFT.RISP_DYN(ni,Mu_dyn,Sigma_dyn,Contrast,Params,delta_t,dyn_value,flg);
                %     Risp_dyn=RISP_DYN(ni,Contrast,Params,delta_t);      % (Mu,Sigma)n --> ... --> (Risp_dyn)n
                if (flg.return_value == 0)
                    fprintf('\n>>> Error in routine RISP_DYN() called by goAB_simul()\n');
                    %         ni_out=ni;
                    return;
                end
                %-------------
                % END DYNAMICS
                %-------------
                
                
                % Have we reached a fixed point (for TOTAL_COUNT times at least)
                count = 0;
                for j=1:nn
                    %         if ( Risp_dyn(j)-ni(j)==0. )
                    if ( abs(Risp_dyn(j)-ni(j))<flg.TOLF )
                        count=count+1;
                    end
                end
                if ( count == nn )
                    tot_count=tot_count-1;
                else
                    tot_count = TOTAL_COUNT;
                end
                
                % If yes, stop:
                if ( tot_count == 0 )
                    fprintf('\n--- Reached Stable Fixed Point by goAB_simul():\n');
                    fprintf('   ');
                    for j=1:nn
                        fprintf(' %g ',ni(j));
                    end
                    fprintf('   delta_t: %g ms\n',delta_t*1000);
                    %%%%%%%%%%%% SAVE TO FILE
                    ni_out=ni;
                    save(ff,'time_steps','rates');
                    save(file_log,'sample_rates','sample_time');
                    return;
                end
                
                ni = Risp_dyn;
                % STORE FOR SAVING
                % every 10 k cycles, take a snapshot of dynamics
                if (rem(k,10)== 0)
                    time_steps=[time_steps delta_t*k];
                    rates=[rates; ni(1:MAXPOP)];
                end
                
                % SCREEN SECTION:
                if (rem(k,100)== 0) && flg.screen
                    for j=1:MAXPOP-4
                        if (Contrast(j) > 0.0)
                            fprintf('(+)%f ',ni(j));
                        elseif (Contrast(j) < 0.0)
                            fprintf('(-)%f ',ni(j));
                        else
                            fprintf(' %f ',ni(j));
                        end
                    end
                    fprintf(' %f(bg)  %f(PV) %f(SST) %f(VIP) ',ni(nn-3),ni(nn-2),ni(nn-1),ni(nn));
                    fprintf('  time: %gms\n',k*delta_t*1000);
                end %if (rem(k,10)== 0)
                
                % Do mnewt investigation:
                if (rem(k,2000) == 0 ) % That's every 100 ms if delta_t = 0.1 ms.
                    %         fprintf('\n\n --- goAB_simul-> Mnewt investigation...\n');
                    
                    % store to print to file .file_log temporary frequencies:
                    sample_rates=[sample_rates; ni];
                    sample_time=[sample_time k*delta_t];
                    
                    %         Ni_temp = ni;
                    %         Ni_temp=mnewt(flg.NTRIAL,Ni_temp,1e-9,1e-9,Params);
                    if (flg.return_value == 0)
                        % if error occurred skip mnewt:
                        flg.return_value = 1;
                        fprintf('\n>>> Error in routine mnewt() called by goAB_simul()\n');
                        fprintf('\n--- Mnewt skipped. Dynamics again...  (J+: %g)\n',Jplus);
                        % if(delay) sleep(2);
                    else
                        %             [res,~,~,~]=Test_Stability(Ni_temp,Params);
                        %             if res
                        %                 fprintf('   mnewt-fixed point is stable - back to dynamics...\n');
                        %             else
                        %                 fprintf('   mnewt-fixed point is unstable - back to dynamics...\n');
                        %             end%sleep(2);
                        %             fprintf('\n--- Dynamics again...  (J+: %g) -- total_count=%d\n\n',Jplus,tot_count);
                        for j=1:MAXPOP-4
                            fprintf(' %0.04g ',ni(j));
                        end
                        fprintf(' %f(bg)  %f(PV) %f(SST) %f(VIP) ',ni(nn-3),ni(nn-2),ni(nn-1),ni(nn));
                        fprintf(' [count:%d]',count);
                        fprintf('>time:%gms\n',k*delta_t*1000);
                    end
                end
                k=k+1;
                
                % Have we a maximum number of steps?
                if (kmax ~= 0)  &&  (k >= kmax)
                    fprintf('\n');
                    ni_out=ni;
                    save(ff,'time_steps','rates');
                    save(file_log,'sample_rates','sample_time');
                    return;
                end
                
            end
            %%%%%%%%%%%% SAVE TO FILE
            ni_out=ni;
            save(ff,'time_steps','rates');
            save(file_log,'sample_rates','sample_time');
            
            return;
        end
        
        
        
        
    end
end