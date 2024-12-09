function Thres=MFT_fun_amitmascaro_run_Stim_quench(OPTIONS,Thres)

% OPTIONS=options(i_c,i_J);

utils.v2struct(OPTIONS);

fprintf('Requires all network attractors from MFT_main\n');
fprintf(' --- running Amit-Mascaro for J+=%0.03g\n',Jplus)
% MFT parameters
ntrial=20000;
tolx=1e-6;
tolf=1e-6;
kmax=2000;


% MFT_preamble;
options_main=OPTIONS;
aux.MFT_preamble;

params=load(paramsfile);
Network=params.Network;
Network.cue_mode=cue_mode;
if cue_mode
    flg.cue_mode=1;
else
    flg.cue_mode=0;
end
params.paramsfile=paramsfile;
params.name=name;
params.Opt=Opt;
flg.paramsfile=paramsfile;
stim_mode=0;
flg.nivsj_mode=0; flg.spur_mode=0; flg.all_mode=1; flg.DEBUG=0;
Jplus=OPTIONS.Jplus;
save(paramsfile,'Jplus','-append');
% bla
%% load stable states
StableStates=[]; % dim: [#states, pops]
p=params.p;
cnt=0; 

% all other fixed points
flg.all_mode=1;
Files=auxMFT.File_Info_MFT(p+4,p,flg);
fileload=Files.stable;
AllData=load(fileload);
% find states
Ind=find(cell2mat(arrayfun(@(x)~isempty(x(:).Jplus),AllData.StableStates,'uniformoutput',false))); % states with Ind active clusters
% Select attractors for current Jplus
for n=1:numel(Ind)
    indJ=find(abs(AllData.StableStates(Ind(n)).Jplus-Jplus)<1e-3);
%     [~,indJ]=min(abs(AllData.StableStates(Ind(n)).Jplus-Jplus));
    if ~isempty(indJ)
        cnt=cnt+1;
        StableStates(cnt,:)=AllData.StableStates(Ind(n)).rates(indJ,:);
    end
end
flg.all_mode=0;

%% thresholds
ni=[params.ni_e params.ni_i];
% compute thresholds on first run
if isempty(Thres)
    [Thres,flg]=aux.fun_Fix_Threshold_InhType(ni,0,params,flg);
end

% STIM MODE
Network.cue_mode=OPTIONS.cue_mode;
if exist('stimset','var')
    stim_mode=1;
    nivsj_mode=1;
    Network.SelPop=1;
    StimHeight=stimset;
    save(paramsfile,'Network','StimHeight','-append');
end
Files=auxMFT.File_Info_MFT(p+4,p,flg);
filesave=[Files.file_amitmascaro '' extra '_results.mat'];
nn=p+4;
save(paramsfile,'nn','-append');
%% MFT
% focus=pops in focus
% mft=population to be integrated out (running mft)
ind_focus=1:2; % indices of populations in focus
ind_mft=setxor(1:p+4,ind_focus);
nu_grid=repmat(struct('focus',[]),1,numel(ind_focus));
% create (nu_1,nu_2) grid
nu_min=min(min(StableStates(:,ind_focus)));
nu_max=max(max(StableStates(:,ind_focus)));
nu_min=nu_min*0.1; nu_max=nu_max*1.5;
for i=1:numel(ind_focus)
    nu_grid(i).focus=linspace(nu_min,nu_max,nbins_grid(i));
end

% initial conditions: for smallest value of nu_1,nu_2, get closest fixed
% point.
firstpt=[nu_grid(1).focus(1) nu_grid(2).focus(1)]; % smallest value of nu_1,nu_2 (bottom left in heat map plot)
[~,ind0]=min(sum((repmat(firstpt,size(StableStates,1),1)-StableStates(:,ind_focus)).^2,2)); % closest FP values for mft pops
ni_start(1:numel(ind_mft))=StableStates(ind0,ind_mft); % initial condition

% reduce parameters
            % save extra options in paramsfile here
            % if cue on, put nonzero sigma_ext
%             stim_mode=0;
%             stim_mode=1;
Params0=aux.In_Parameters_InhType(Thres,paramsfile,flg);
Params0.ind_mft=ind_mft;
Params0.ind_focus=ind_focus;
Params0.Network=Network;
tic
OUT=aux.MFT_RATE_eff(ni_start,Params0,nu_grid,flg);
toc
OUT.nu_grid=nu_grid;
OUT.StableStates=StableStates;    
% rates_eff with dim [p,nbins_grid(1),nbins_grid(2)] contains the value of the integrated out nu_3...nu_p=2
% given the fixed values nu_grid(1).focus, nu_grid(2).focus
save(filesave,'OUT','params');
fprintf('results saved in %s\n',filesave);

