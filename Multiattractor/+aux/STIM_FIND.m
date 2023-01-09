% Compute spurious states starting from nn=3 and increasing population
% # until nn=p+2.
% 1) For nn=3, compute the stable spontaneous state (n_+(3),nbg(3),nI(3)), existing for all
% Jplus. Check it is a stable fixed point (FP)
% 2) Ungroup new selective population n1 and set initial point n1(4)=n_+(4)=n+(3),
% nbg(4)=nbg(3), nI(4)=nI(3). Check this it is an unstable FP, compute
% eig and find Eig.vec v(4) corresponding to positive eigenvalue.
% 3) Perturb n along v(4) with negative n1(4) and positive n+(4) and
% evolve dynamics until new stable FP reached (spurious state).
% 4) iterate procedure until spurious state is stable for large enough
% nn.

function [Params,flg]=STIM_FIND(Params,Jstop,STEP,flg)

% global return_value NTRIAL TOLX flg.TOLF TOLSPU screen contr overlap nivsj_mode spur_mode BLA syn_dyn stim_modes

% parameter file for spurious
paramsfile=Params.Files.paramsfile;
load(paramsfile);

if ~any(strcmp(fieldnames(Network),'SelPop'))
    fprintf('\n Number of selective pops not specified. Setting: 1 sel pop...\n');
    SelPop=1;
end
nn=4; % start with stim,non-stim,bg,inib only
flg.nivsj_mode=0;
%----------------------------------------------
% STARTING POINT, spontaneous activity for nn=4
%----------------------------------------------
Jplus = Jzero;
nip0=0; % set temporarily to zero;
save(paramsfile,'Jplus','nn','Network','-append');
Params=auxMFT.Learn(Params,flg);
% % % paramsfile=Params.Files.paramsfile;
% Starting point for nn=3
for i=1:nn-1  
    Ni(i) = ni_e;
end
Ni(nn) = ni_i;
if (flg.screen)
    fprintf('--- Initial point for nn=%d:\n\n',nn);
    for i=1:nn  
        fprintf('    Ni[%d]: %f\n',i,Ni(i));
    end
    fprintf('\n');
end
% % check it's a fixed point
% Ni=mnewt(flg.NTRIAL,Ni,flg.TOLX,flg.TOLF,Params);
% save nn=3 spontaneous activity rates in ->spur files
if (flg.screen)
    fprintf('--- STIM_FIND->RATE_VS_Jfun for nn=4 \n');
end

aux.RATE_VS_Jfun(Ni,Jstop,Params,STEP,Opt,flg);
if (flg.return_value == 0)
    fprintf('\n--- Error in routine RATE_VS_Jfun() called by main()\n');
end
fprintf('\n--- STIM_FIND->RATE_VS_Jfun() computed nn=4 spont. activity\n');

% Starting point:
Jplus = Jzero;
% create JplusArray
JplusArray=[];
while Jplus<=Jstop
    JplusArray=[JplusArray Jplus];
    step=auxMFT.StepJplus(Jplus,Opt,STEP,flg);
    Jplus =Jplus + step;
end

% guess for location of first critical point given value of ni_ext for
% stimulated pops
x=[1 1.1 1.15 1.2 1.3 1.4]; % ni_ext values 
J1=[4.1 4.5 4.7 4.9 5.3 5.7]; % (from selective activity with all cluster stimulated)
P=polyfit(x,J1,2);
JKickFun=@(t)P(1)*t.^2+P(2)*t+P(3);
if exist('StimHeight','var')
    if isempty(StimHeight)
        height_stim=box_height(paramsfile);
    elseif ~isempty(StimHeight)
        height_stim=StimHeight;%box_height(paramsfile);
    end
end
JKick=JKickFun(1+StimHeight)-0.1;
fprintf('\n JKick starts at %g...\n',JKick);

% figure(1); clf; h=[]; 
% h(1)=plot(x,J1,'o','color','r','markersize',10); 
% hold on;
% h(2)=plot(x,Y,'color','r','linewidth',1.5);

%-----------------------
% START SEARCH for nn>4
%-----------------------
% at each loop iteration, ungroup one more selective population
SelPop=Network.SelPop;
nip0=0;
for nn=5:SelPop+3%p+2
    nip0=nip0+1;
    % Starting points:
    load(Params.Files.FP,'rates','Jplus_store');
    Jplus_old=Jplus_store; %use nn-1 array of 
    rates_old=rates; % load nn-1 spurious fixed points
    Jplus = Jzero;
    save(paramsfile,'Jplus','nn','nip0','-append');
    Params=auxMFT.Learn(Params,flg);
    fprintf('\n--- Learned parameters for nn=%d\n',nn);
    % store variables to save
    Jplus_store=[];
    J_high=[]; % store Jplus values below 2nd critical point (stable symmetric activity)
    J_unstable=[]; % store Jplus values above 2nd critical point (stable spurious states)
    rates=[]; % save all stable fixed points
    rates_high=[]; % save new fixed points only
    rates_unstable=[];
    Eig.Real=[];
    Eig.Imag=[];
    Eig.vec=[]; % dim: 1st=pop #; 2nd=eigenvector #; 3rd=Jplus step
    % file to store fixed points
    file_name=Params.Files.FP;
    IndMatch=[];
    NewGuess=0;
    for Jplus=JplusArray
        % update Jplus in paramsfile to be loaded by function Learn
        save(paramsfile,'Jplus','-append');
        Params=auxMFT.Learn(Params,flg);
        % info on flg.screen
        if (flg.screen) 
            fprintf('\n    (Jplus = %g)\n\n',Jplus);
        else
            fprintf('--- Jplus = %g\n',Jplus);
        end
        if NewGuess
            % use previous Jplus step as initial condition
            Ni_start=rates_high(end,:);
            Ni_start(1:nip0)=Ni_start(nip0)+2;
            Ni_start(nip0+1)=Ni_start(nip0+1)-1;
        else
            % Does Jplus appear in Jplus_old fixed points for nn-1?
            % - if it did not appear, pick as new starting point cnt+1 from
            % rates_old
            % - if it appeared, create candidate fixed point from
            % rates_old(cnt,:)
            IndMatch=find(Jplus==Jplus_old);
    %         fprintf('size rates_old: [%d %d]',size(rates_old,1),size(rates_old,2));
    %         Ni_old(1)=rates_old(cnt,1); % new ungrouped selective pop
            if isempty(IndMatch)
                IndMatch=find(Jplus_old>Jplus,1); 
                fprintf('\n --- Jplus=%g missing from nn=%d data, use Jplus=%g\n',Jplus,nn-1,Jplus_old(IndMatch));
            end
            % Set rate initial conditions at current Jplus value 
            Ni_start = rates_old(IndMatch,:);
            Ni_start=[Ni_start(1) Ni_start]; % add ungrouped stimulated pop in 1st col
            if nn==5
                if Jplus>JKick-0.05
                    % kick it
                    Ni_start(1)=30;
                    if Ni_start(2)>3
                        Ni_start(2)=Ni_start(2)-3;
                    end
                    fprintf('\n Kick...\n');
                end
            else % if nn>5
                if (flg.screen) 
                    fprintf('Rates loaded from nn=%d\n',nn-1);
                    for i=1:nn  
                        fprintf('    Ni[%d]: %f\n',i,Ni_start(i));
                    end
                    fprintf('\n');
                end
                % KICK:
                % 1) if all stim pops have same rate (high), bring down grouped
                % ones
                % 2) if ungrouped ones are lower already, leave it
                roundedData = round(Ni_start*1e4)/1e4; % round to 1e-6
                UniqueRates=unique(roundedData(1:nip0+1));
                NumRates=numel(UniqueRates); % number of different firing rates
                if numel(UniqueRates)==1
                    fprintf('\n Bring down grouped pops; up ungrouped pops\n');
    %                     Ni_start(1:nip0)=Ni_start(1:nip0);
                    Ni_start(nip0+1)=Ni_start(nip0+2);
                else
                    HiRate=sum(roundedData==UniqueRates(2)); % # of high rates
                    LowRate=sum(roundedData==UniqueRates(1)); % # of low rates
                    if LowRate>HiRate
                        fprintf('\n LowRate>HiRate: problem, initial guess 2 steps earlier...\n');
                        % if # of low rates is larger than high rates, something
                        % went wrong at the previous step -> get initial points
                        % from two steps earlier
                        Ni_start = rates_old(IndMatch-1,:);
                        Ni_start = [Ni_start(1) Ni_start];                    
                        Ni_start(nip0+1)=Ni_start(nip0+2);
                    end
                end
                if Jplus>JKick
                    Ni_start(1:nip0)=Ni_start(1:nip0)+20;
                end
            end
        end
        if (flg.screen)
            fprintf('New initial conditions\n');
            for i=1:nn  
                fprintf('    Ni[%d]: %f\n',i,Ni_start(i));
            end
            fprintf('\n');
        end
        % give it a little kick
        Cycles=5000;
        if StimHeight>0.1
            Cycles=10000;
        end
        [Ni_start,flg]=auxMFT.goAB_simul(Ni_start,Cycles,Params,flg);
        % check it is still a fixed point
        [Ni_check,flg]=auxMFT.mnewt(flg.NTRIAL,Ni_start,flg.TOLX,flg.TOLF,Params,flg);
        MoreDynamics=0;
        if (flg.return_value == 0)
            fprintf('\n--- Error in routine mnewt() called by STIM_FIND\n');
            flg.return_value=1;
            MoreDynamics=1;
        end
        if (flg.screen)
            fprintf('Rates found by mnewt\n');
            for i=1:nn  
                fprintf('    Ni[%d]: %f\n',i,Ni_check(i));
            end
            fprintf('\n');
        end
        errx=sum(abs(Ni_start-Ni_check));
%             if ( errx==0. )  
        if (errx>flg.TOLSPU ) || MoreDynamics  
            fprintf('\n --- Using nn=%d fixed point as initial guess for nn=%d',nn-1,nn);
            fprintf('\n --- Look for new fixed point using dynamics\n');
%             flg.return_value=0;
            if MoreDynamics
                Ni_check=Ni_start; % because mnewt screwed it up
                Ni_start=[Ni_start(1) Ni_start]; % because mnewt screwed it up
            end
            [Ni_check,flg]=auxMFT.goAB_simul(Ni_check,0,Params,flg);
%             continue;
            fprintf('\n --- Found new fixed point \n');
        end
        % Does new fixed point have higher rate for ungrouped vs. grouped
        % pops?
        CheckHiRate=2;
        roundedData = round(Ni_check*1e4)/1e4; % round to 1e-6
        UniqueRates=unique(roundedData(1:nip0+1));
        NumRates=numel(UniqueRates); % number of different firing rates
        LowRate=sum(roundedData==UniqueRates(1)); % # of low rates
        NewGuess=0;
        if NumRates==1
            fprintf('   Same rate for all stim pops: drop fixed point\n');
        elseif NumRates>2
            fprintf('   Fixed point has %d different rates for stim pops: weird, drop it.\n',NumRates);
        elseif NumRates==2 && LowRate==1
            % if # of high rates is nip0 and low rates =1, keep it, if not, drop it (already
            % stored at an earlier step)
            fprintf('   %d ungrouped pop > grouped pops\n',nip0);
            fprintf('   new fixed point with %d active clusters\n',nip0);
            % check fixed point stability
            [res, Real_part, Imag_part, EigV,flg]=auxMFT.Test_Stability(Ni_check,Params,flg);
    %         % store to save
            Eig.Real=[Eig.Real Real_part];
            Eig.Imag=[Eig.Imag Imag_part];
            Eig.vec=cat(3,Eig.vec,EigV);
            %----------------------------
            % IS FIXED POINT STABLE?
            %----------------------------
            if ~any(Real_part>0) 
                fprintf('   Stored.\n');
                % If NumRates==5, we have a new fixed point
                J_high=[J_high Jplus]; % Jplus values below 2nd critical point
                rates_high=[rates_high; Ni_check];
                Eig.Real=[Eig.Real Real_part];
                Eig.Imag=[Eig.Imag Imag_part];            
                Eig.vec=cat(3,Eig.vec,EigV);
                NewGuess=1;
            %----------------------------
            % IS FIXED POINT UNSTABLE? (above second critical point)
            %----------------------------
            elseif any(Real_part>0)
                fprintf('   fixed point with %d active clusters is unstable\n',nip0);
                J_unstable=[J_unstable Jplus]; % Jplus values below 2nd critical point        
                rates_unstable=[rates_unstable Ni_check]; % Jplus values below 2nd critical point        
            end
            % store variables
        elseif NumRates==2 && LowRate>1
            fprintf('\n Already stored fixed point earlier...\n');
        end
        rates=[rates; Ni_check];
        Jplus_store=[Jplus_store Jplus];
        
    end % end of Jplus
    %---------------
    % SAVE TO FILE
    %---------------
    % fixed point rates
    if (flg.screen)
        fprintf('--- saving (nn=%d) to %s.\n',nn,file_name);
    end
    AllParameters=load(paramsfile);
    save(file_name,'J_high','rates_high','rates','J_unstable',...
        'rates_unstable','Jplus_store','AllParameters','Eig');
    
    if (flg.return_value == 0)
        fprintf('---Error: STIM_FIND for nn=%d...',nn);
        return;
    end
end % end of nn
fprintf('--- STIM_FIND ended! Congratulations...\n');
return;

