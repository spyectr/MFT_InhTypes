
function [res_stable res_unstable flg]=RATE_VS_Jfun_ALL(n_hi,ni,Jstop,Params,Opt,flg)
% ni=Ni


% Files
paramsfile=Params.Files.paramsfile;
file_eig=Params.Files.file_eig;
% file_single_eig=Params.Files.file_single_eig;
if (flg.screen)
    fprintf('Loading %s',paramsfile);
end

load(paramsfile);
nn=numel(ni);
if 1%nn<50
    MAX_DATA=nn;
else
    MAX_DATA=50;
end
% previous version:
% flg.STEP=0.1;

step = flg.STEP;
% k=10;
flag=0;
% stop_control = 1;

% variables to be saved to file
% STABLE
rates=[];
Jplus_store=[];
Eig=struct('Real',[],'Imag',[],'vector',[],'matrix',[],'phiprime',[],'linresp',[]); % dim: 1st=pop #; 2nd=eigenvector #; 3rd=Jplus step
% UNSTABLE
rates_un=[];
Jplus_store_un=[];
Eig_un=struct('Real',[],'Imag',[],'vector',[],'matrix',[],'phiprime',[],'linresp',[]); % dim: 1st=pop #; 2nd=eigenvector #; 3rd=Jplus step

% Starting point:
Jplus = Jzero;
% create JplusArray
JplusArray=[];
while Jplus<=Jstop
    JplusArray=[JplusArray Jplus];
    step=auxMFT.StepJplus(Jplus,Opt,flg.STEP,flg);
    Jplus =Jplus + step;
end
load(Params.Files.rates_old,'rates_old','J_old');
if n_hi==p
    % load spontaneous activity
    filespont=auxMFT.File_Info_MFT(5,p,flg);
    SpontData=load(filespont.FP,'Jplus_store','rates');
    rates_old=SpontData.rates;
    J_old=SpontData.Jplus_store;
end
rates_old_store=NaN(numel(JplusArray),nn); % store rates for next n_hi cycle
cnt=0;
for Jplus=JplusArray
    skip_mode=0; flg.return_value=1;
    cnt=cnt+1;
    % update Jplus in paramsfile to be loaded by function Learn
    fprintf('Jplus=%g\n',Jplus);
    save(paramsfile,'Jplus','-append');
    Params=auxMFT.Learn(Params,flg);
    % Investigate Selective Activity?
    % if ni_hi==p and Jplus gives unstable spont activity, end
    %     if n_hi==p && Jplus>J_old(end)
    %         continue;
    %     end
    if (flg.nivsj_mode == 1)
        % NO OVERLAP CASE:
        if ~isempty(strfind(Opt,'Small2'))
            if all(all(rates_old==0))
                JKick=1.1;
                Nkick=50;
%                 if flg.cue_mode==1
%                     Nkick=40;
%                     acue=abs(Network.cue_value(1));
%                     if acue<=0.01
%                         JKick=10.25;
%                     elseif acue>0.01 && acue<=0.02
%                         JKick=10.4;
%                     elseif acue>0.02 && acue<=0.03
%                         JKick=10.55;
%                     elseif acue>0.03 && acue<0.04
%                         JKick=10.7;
%                     elseif acue>=0.04 && acue<0.045
%                         JKick=10.9;
%                     elseif acue>=0.045 && acue<0.05
%                         JKick=11.05;
%                     elseif acue>=0.05
%                         JKick=11.1;
%                     end
%                 end
%                 if (flag < 5)  % this means do it only 3 times (because it's time consuming)
                    if (Jplus>JKick) && n_hi<p
                        ni(1:n_hi) = Nkick; % ni(nn)=10.;
                        ni(n_hi+1:end-4) = 2; % ni(nn)=10.;
                        ni(end-3) = 2; % ni(nn)=10.;
                        ni(end-2:end) = 10; % ni(nn)=10.;
                        ni0=ni;
                        fprintf('RATE_VS_Jfun: ni(1)=%0.03g Hz kick starts at Jplus>%0.03g \n',ni(1),JKick); % debug
                        [ni,flg]=auxMFT.goAB_simul(ni,5000,Params,flg);
                        if isempty(ni)
                            ni=ni0;
                        end
                        % %                             flag=flag+1;
                    end
%                 end
            elseif any(rates_old(cnt,:))
                if flg.stim_mode==0 || (flg.stim_mode==1 && n_hi>1)
                    fprintf('Loading rates from rates_old...\n');
                    ni=rates_old(cnt,:);
                    ni(1:n_hi)=ni(1)/2;
                    if n_hi==p % use spontaneous
                        ni=zeros(1,p+4);
                        ni(1:n_hi)=rates_old(cnt,1);
                        ni(end-3:end)=rates_old(cnt,end-3:end);
                    end
                    %                         if n_hi==4 && Jplus<4.4
                    %                             ni(1:n_hi)=ni(1:n_hi)-3; % more active clusters -> slightly lower rates
                    %                         else
                    %                             ni(1:n_hi)=ni(1:n_hi)-5; % more active clusters -> slightly lower rates
                    %                         end
                elseif flg.stim_mode==1 && n_hi==1
                    if Jplus>JKick
                        fprintf('\n Extra kick ni(1)+20Hz...\n')
                        % give extra kick to first pop
                        ni(1)=ni(1)+30;
                    end
                end
                if n_hi<p
                    [ni,flg]=auxMFT.goAB_simul(ni,20000,Params,flg);
                end
            elseif all(isnan(rates_old(cnt,:)))
                % if previous n_hi didn't have active clusters -> skip
                skip_mode=1;
            end
        end
    end % if(flg.nivsj_mode == 1)
    if ~skip_mode
        %     /*
        %     % Help mnewt search for low rate fixed points (usually not needed though; goAB_simul could be used as well)
        %     goAB(ni,600);
        %     if (flg.return_value == 0)  BYE
        %     */
        if (flg.nivsj_mode==0) && (flg.stim_mode==1)
            fprintf('help out with dynamics...');
            [ni,flg]=auxMFT.goAB_simul(ni,5000,Params,flg);
        end
        % Search fixed point:
        fprintf('RATE_VS_Jfun_ALL: mnewt \n'); % debug
        NTRIAL=5000;
        ni0=ni; % store values of ni kicked
        [ni,flg]=auxMFT.mnewt(NTRIAL,ni0,Params,flg);
        if (flg.return_value == 0)
            fprintf('\n>>> Error in routine mnewt() called by RATE_VS_Jfun_ALL()\n');
            fprintf('>>> Run goAB_simul instead\n');
            flg.return_value =1;
            [ni,flg]=auxMFT.goAB_simul(ni0,20000,Params,flg);
            if isempty(ni)
                ni=ni0;
                continue; % skip this step and move to next Jplus value
            else
                fprintf('>>> Recheck goAB_simul output with mnewt\n');
                [ni,flg]=auxMFT.mnewt(NTRIAL,ni,Params,flg); % double check new fixed point
            end
        end
        
        
        % handle unstable points:
        %     fprintf('RATE_VS_Jfun: Test_Stability \n'); % debug
        [res, Real_part, Imag_part, EigV,flg,matrix,phiprime,linresp]=auxMFT.Test_Stability(ni,Params,flg);        
        if (~res) %&& (stop_control == 1))
            fprintf('Unstable fixed point.\n');
            % save unstable state
            Jplus_store_un=[Jplus_store_un Jplus];
            rates_un=[rates_un; ni];
            Eig_un.Real=[Eig_un.Real Real_part];
            Eig_un.Imag=[Eig_un.Imag Imag_part];
            Eig_un.vector=cat(3,Eig_un.vector,EigV);
            Eig_un.matrix=cat(3,Eig_un.matrix,matrix);
            Eig_un.phiprime=cat(2,Eig_un.phiprime,phiprime);
            Eig_un.linresp=cat(3,Eig_un.linresp,linresp);
            %         % run dynamics anyway, see where it goes
            %         fprintf('\n\n--- Unstable state occurred at J+ = %g\n',Jplus);
            %         for j=1:nn
            %             fprintf(' %g ',ni(j));
            %         end
            %         fprintf('\n');
            %         % each time unstable point occurs, stop and ask for dynamics simulation:
            %         ni=Cp_ni;
            %         fprintf('\n    Do Dynamics Simulation...\n');
            %         ni =ni + .001;  % Little perturbation
            %         ni=goAB_simul(ni,0,Params); % 0 := no max number of steps.
            %         ni_Test=ni;
            %         fprintf('\n    Done dynamics sim... checking stability. \n');
            %         [res, Real_part, Imag_part, EigV]=Test_Stability(ni_Test,Params);
            %         fprintf('\n Stability:%d\n',res);
        end % end if ~res stop_control
        fprintf('-> ');
        for j=1:nn
            fprintf(' %g ',ni(j));
        end
        fprintf('\n');
        result=auxMFT.CheckRates(ni,n_hi);
        fprintf('    fixed point:');
        for j=1:nn
            fprintf(' %g ',ni(j));
        end
        fprintf('\n');
        % store stable states only
        if result && res
            % store to save fixed points to file 'nivsj.mat':
            Jplus_store=[Jplus_store Jplus];
            rates=[rates; ni(1:MAX_DATA)];
            % store to save
            Eig.Real=[Eig.Real Real_part];
            Eig.Imag=[Eig.Imag Imag_part];
            Eig.vector=cat(3,Eig.vector,EigV);
            Eig.matrix=cat(3,Eig.matrix,matrix);
            Eig.phiprime=cat(2,Eig.phiprime,phiprime);
            Eig.linresp=cat(3,Eig.linresp,linresp);
            fprintf('\n --- Fixed point with %d Active Clusters ---> stored...\n',n_hi);
            % store rates in rates_old for next n_hi cycle
            % check if fixed point has n_hi pops at high rate
            rates_old_store(cnt,:)=ni;
        elseif ~result
            fprintf('\n --- Fixed point has wrong number of active clusters ---> rejected... \n ');
        end
    end
    % % % % % %     if Jplus==13 && flg.DEBUG
    % % % % % %         data(cnt).Params=Params;
    % % % % % %         data(cnt).rates=ni;
    % % % % % %     end
end % while(Jplus<=Jstop)
% stable
res_stable.Jplus=Jplus_store;
res_stable.rates=rates;
res_stable.Eig=Eig;
% unstable
res_unstable.Jplus=Jplus_store_un;
res_unstable.rates=rates_un;
res_unstable.Eig=Eig_un;
% save rates_old as initial guesses for next n_hi cycle
rates_old=rates_old_store;
J_old=Jplus_store;
save(Params.Files.rates_old,'rates_old','J_old');
% Summary of parameters:
fprintf('\n\n--- Done.\n');
auxMFT.Summary(paramsfile);

% save flg.DEBUG
if flg.DEBUG
    AllParameters=load(paramsfile);
    save('debug_values.mat','AllParameters','data');
end
%
% %---------------
% % SAVE TO FILE
% %---------------
% % file to store fixed points
% file_name=Params.Files.FP;
% % file to store eigenvalues and eigenvectors
% file_eig=Params.Files.file_eig;
% % fixed point rates
% AllParameters=load(paramsfile);
% save(file_name,'Jplus_store','rates','AllParameters');
% % all eigenvalues
% save(file_eig,'Jplus_store','eigenvalues','eigenvectors','AllParameters');
%

%----------------
% AUX FUNCTIONS
%----------------
%----------------
% check candidate fixed point has n_hi pops at high rate in the first
% 1:n_hi pops
% 
% function result=auxMFT.CheckRates(ni,n_hi)
% 
% result=0;
% TolRate=1e-4;
% p=numel(ni)-4;
% MaxRate=max(ni(1:p));
% MaxPops=sum(abs(ni(1:p)-MaxRate)<TolRate); % # of pops with high rate among first n_hi
% if MaxPops==n_hi
%     result=1;
% end


return;
