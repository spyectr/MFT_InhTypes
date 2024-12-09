%/********************************************************
% *                                                      *
% *       Initializations in nn populations case         *
% *       ------  Important: nn >= 3(*) --------         *
% *                                                      *
% * (*)It works also with nn=3 for in this case          *
% *    nn_sel=0 and therefore loops on selective         *
% *    populations are skipped.                          *
% *    See code for a better understanding (see also     *
% *    routine In_Param_3pop() for a match).             *
% *                                                      *
% ********************************************************/

function Params=In_Parameters_Stim1(theta,paramsfile,flg)

% global flg.DEBUG  flg.flg.flg.nivsj_mode spur_mode stim_mode all_mode

% load parameters
load(paramsfile);
% connectivity
Cee = N_e*pee_matrix; % 1600;         % 1600
Cie = N_e*pie_matrix; % 1600;         % 1600
Cii = N_i.*pii_matrix; % 400;         % 400
Cei = N_i.*pei_matrix; % 400;         % 400
Cext = (N_e)*pext_matrix; % 1600;

% set thresholds
theta_e=theta(1); % exc
theta_i=theta(2:end); % inh cell types

Params=[]; % output structure 

% POPULATIONS
% minimal number of populations is nn=5: exc (OTHER_SEL), exc (BACK_EXC), inh PV, inh SST, inh VIP 
OTHER_SEL=nn-4;
BACK_EXC =nn-3;
INIB     =nn-2:nn;

% % external currents
% if ~exist('ni_ext_e','var') && ~exist('ni_ext_i','var')
%     ni_ext_e=ni_ext;
%     ni_ext_i=ni_ext;
% end

% nn_sel = number of populations selective to a single pattern, which we do not
% want to group with all the other p-nn_sel
nn_sel = nn-5;
%------------------------------------------------------------------
% Synapses from sel. neurons to bg neurons.
% AB97: Jminus.  
%------------------------------------------------------------------*/
Jminus = 1.-gam*f*(Jplus-1.);
if Jminus<0
    error('ERROR in create_params: Jminus<0. Adjust gam or f\n');
end
Jsel_to_bg = Jminus; %% 1.;  % B99
Jcor = 0.;%%-0.04*Jplus; % 1.;  % B99
% % TO AVOID PROBLEMS:
x = 1.; % in AB97 local exc currents have x while the external current 1-x

% higher intra-cluster connectivity in DOIRON
HiIntraCluster=0;
% if any(strcmp(fieldnames(Network),'Ree')) && strcmp(Opt,'DOIRON')
%     HiIntraCluster=1;
%     pee=Cee/N_e;
%     REE=Network.Ree;
%     denom=REE*(N_e/p-1)+(N_e/p)*(p-1);
%     peeout=pee*(N_e-1)/denom; % inter-cluster connectivity
%     peein=REE*peeout;  % intra-cluster connectivity
%     Cee=peeout*N_e;
% end

%/*-------------------------------------------------------
% *                                                      *
% *                Matrix A initialization               *
% *                                                      *
% -------------------------------------------------------*/

%/*-----------------------------------------------
% connections TO splitted selective populations:
% -----------------------------------------------*/
for i=1:nn_sel
%     % sel -> itself:
    A(i,i) = tau_e*f*Cee*x*Jplus*Jee;
    if HiIntraCluster
        A(i,i) = A(i,i)*peein/peeout;
    end
    for j=1:nn
        %       % sel(j) -> sel(i):
        if (j~=i) && (j<=nn_sel)  
          A(i,j) = tau_e*f*Cee*x*Jminus*Jee;
        end


        %       % all other sel -> sel:
        if (j==OTHER_SEL)  
          A(i,j) = tau_e*f*(p-nn_sel)*Cee*x*Jminus*Jee;
        end

        %       % bg -> sel:
        if (j==BACK_EXC)   
%           A(i,j) = tau_e*(1.-p*f)*Cee*x*Jminus*Jee;
          A(i,j) = tau_e*(1.-p*f)*Cee*x*Jsel_to_bg*Jee;
        end

        %       % inib -> sel:
        if any(j==INIB)       
          A(i,j) = -tau_e*Cei(j==INIB)*Jei(j==INIB); 
        end
    end
end



% /*--------------------------------------------------------------
% connections TO other p-nn_sel COLLECTED selective populations:
% --------------------------------------------------------------*/
i=OTHER_SEL;
for j=1:nn_sel
    % sel -> all other sel collected: 
    A(i,j) = tau_e*f*Cee*x*Jminus*Jee;
end
% all other sel collected -> itself:
Temp=((p-nn_sel-1)*Jminus+Jplus);
A(i,i)        = tau_e*x*f*Cee*Temp*Jee;
if HiIntraCluster
    Temp1=((p-nn_sel-1)*Jminus+(peein/peeout)*Jplus);
    A(i,i) = A(i,i)*Temp1/Temp;
end
% bg  -> all other sel collected:
A(i,BACK_EXC) = tau_e*(1.-f*p)*Cee*x*Jsel_to_bg*Jee;
% inib  -> all other sel collected:
A(i,INIB)     = -tau_e*Cei.*Jei;    

% /*--------------------------------------------------------------
% connections TO bg (=eccit. nonselective) population:
% --------------------------------------------------------------*/
i=BACK_EXC;
for j=1:nn_sel
    % sel ->  bg: 
    A(i,j) = tau_e*f*Cee*x*Jsel_to_bg*Jee;
end
% all other sel (collected) ->  bg:
A(i,OTHER_SEL) = tau_e*f*(p-nn_sel)*Cee*x*Jsel_to_bg*Jee;
% bg  ->  itself:
A(i,i)         = tau_e*(1-p*f)*Cee*x*Jee*(1+Jcor);
% inib  ->  bg:
A(i,INIB)      = -tau_e*Cei.*Jei;

% /*--------------------------------------------------------------
% connections TO inhibitory population:
% --------------------------------------------------------------*/
i=INIB;
for j=1:nn_sel
    % sel -> inib:
    A(i,j) = f*x*tau_i'.*Cie.*Jie;
end
% all other sel -> inib:

A(i,OTHER_SEL)= f*(p-nn_sel)*x*tau_i'.*Cie.*Jie;
% bg -> inib:
A(i,BACK_EXC)  = (1-p*f)*x*tau_i'.*Cie.*Jie;
% inib -> itself:
A(i,i)         = -repmat(tau_i',1,3).*Cii.*Jii; 

if flg.DEBUG
  printf('Matrix A:\n');
  for i=1:nn
    for j=1:nn  
        printf(' %f ',A(i,j)); 
    end
    printf('\n');
  end
end

% /*-------------------------------------------------------
% *                                                      *
% *                Matrix B initialization               *
% *                                                      *
% -------------------------------------------------------*/

% /*-----------------------------------------------
% connections TO splitted selective populations:
% -----------------------------------------------*/
for i=1:nn_sel
    % sel -> itself:
    B(i,i) = Jee*Jee*tau_e*f*Cee*x*Jplus*Jplus;
    if HiIntraCluster
        B(i,i) = B(i,i)*(peein/peeout);
    end
    for j=1:nn
        % sel(j) -> sel(i):
        if (j~=i && j<=nn_sel)  
          B(i,j) = Jee*Jee*tau_e*f*Cee*x*Jminus*Jminus;
        end
        % all other sel -> sel:
        if (j==OTHER_SEL)  
          B(i,j) = Jee*Jee*tau_e*f*(p-nn_sel)*Cee*x*Jminus*Jminus;
        end
        % bg -> sel:
        if (j==BACK_EXC)   
          B(i,j) = Jee*Jee*tau_e*(1-p*f)*Cee*x*Jminus*Jminus;
        end
        % inib -> sel:
        if any(j==INIB)       
          B(i,j) = tau_e*Cei(j==INIB)*Jei(j==INIB)*Jei(j==INIB);
        end
    end
end


% /*--------------------------------------------------------------
% connections TO other p-nn_sel COLLECTED selective populations:
% --------------------------------------------------------------*/
i=OTHER_SEL;
for j=1:nn_sel
    % sel -> all other sel collected: 
    B(i,j) = Jee*Jee*tau_e*f*Cee*x*Jminus*Jminus;
end
% all other sel collected -> itself:
Temp=((p-nn_sel-1)*Jminus*Jminus+Jplus*Jplus);
B(i,i)        = tau_e*x*f*Cee*Temp*Jee*Jee;
if HiIntraCluster
    Temp1=((p-nn_sel-1)*Jminus*Jminus+(peein/peeout)*Jplus*Jplus);
    B(i,i) = B(i,i)*Temp1/Temp;
end
% bg  -> all other sel collected: 
B(i,BACK_EXC) = Jee*Jee*tau_e*(1.-f*p)*Cee*x*Jminus*Jminus;
% inib  -> all other sel collected:
B(i,INIB)     = tau_e*Cei.*Jei.*Jei;

%   /*--------------------------------------------------------------
%     connections TO bg (=eccit. nonselective) population:
%     --------------------------------------------------------------*/
i=BACK_EXC;
for j=1:nn_sel
    % sel ->  bg: 
    B(i,j) = tau_e*f*Cee*x*Jminus*Jee*Jminus*Jee;
end
% all other sel (collected) ->  bg:
B(i,OTHER_SEL) = tau_e*f*(p-nn_sel)*Cee*x*Jminus*Jee*Jminus*Jee;
% bg  ->  itself:
B(i,i)         = tau_e*(1-p*f)*Cee*x*Jee*Jee*(1+Jcor*Jcor);
% inib  ->  bg:
B(i,INIB)      = tau_e*Cei.*Jei.*Jei;


% /*--------------------------------------------------------------
% connections TO inhibitory  population:
% --------------------------------------------------------------*/
i=INIB;
for j=1:nn_sel
    % sel -> inib:
    B(i,j) = f*x*tau_i'.*Cie.*Jie.*Jie; 
end
% all other sel -> inib: 
B(i,OTHER_SEL) = f*(p-nn_sel)*x*tau_i'.*Cie.*Jie.*Jie;
% bg -> inib:
B(i,BACK_EXC)  = (1-p*f)*x*tau_i'.*Cie.*Jie.*Jie;
% inib -> itself:
B(i,i)         = Cii*repmat(tau_i',1,3).*Jii.*Jii;
if (delta ~= 0) 
    for i=1:nn
        for j=1:nn 
            B(i,j) =B(i,j) *( 1.+delta);
        end
    end
end

if flg.DEBUG
  printf('Matrix B:\n');
  for i=1:nn
    for j=1:nn 
        printf(' %f ',B(i,j)); 
    end
    printf('\n');
  end
end 

%/*-------------------------------------------------------
% *                                                      *
% *         All Other Parameters Initialization          *
% *                                                      *
% -------------------------------------------------------*/


% external currents: mean
Mu_ext=zeros(1,nn);
for i=1:nn-3
    Mu_ext(i)  = Cext*tau_e*ni_ext_e*Jee_ext;
end
Mu_ext(INIB) = Cext*tau_i.*ni_ext_i.*Jie_ext;

% external currents: variance
Sigma_ext=zeros(1,nn);
if strcmp(Network.stim,'Poisson') 
    for i=1:nn-3
        Sigma_ext(i)  = (1.+delta)*tau_e*Cext*ni_ext_e*Jee_ext^2;
    end
    Sigma_ext(INIB) = (1.+delta)*Cext*tau_i.*ni_ext_i.*Jie_ext.^2;
elseif strcmp(Opt,'DOIRON') 
    Sigma_ext(1:nn-3)=(delta)^2*tau_e*Cext*ni_ext_e*Jee_ext/(1.15); % see normalization of Jee_ext in createparamsair.m
    Sigma_ext(INIB) = (delta/2)^2*Cext*tau_i.*ni_ext_i.*Jie_ext/1.025;
end


%%%%%%%%%%%%%%%%%%%%%%%
% BELOW IS STILL WITHOUT INH TYPES
%%%%%%%%%%%%%%%%%%%%%%%

%----------
% STIMULUS
%----------
if flg.stim_mode
    % # of selective populations
    SelPop=Network.SelPop;
    indfeat=find(cellfun(@(x)~isempty(x),strfind({Stimulus.feat(:).name},'US')));
    p_stim=Stimulus.feat(indfeat).connectivity; % stimulus connectivity

    %-----------------
    % EXTERNAL (with stimulus)
    %-----------------
    %
    % higher kick (box stimulus equivalent to thalamic)
    if exist('StimHeight','var')
        height_stim=StimHeight;%box_height(paramsfile);
    end
    % CURRENT
    % stimulated
    Mu_ext(1:SelPop)  = Cext*ni_ext_e*tau_e*Jee_ext*(1+p_stim*height_stim);
    % VARIANCE (zero unless 'Poisson')
    if strcmp(Network.stim,'Poisson') 
        Sigma_ext(1:nn_sel+1)  = (1.+delta)*Cext*ni_ext_e*tau_e*Jee_ext^2*(1+p_stim*height_stim);
        Sigma_ext(nn_sel+2:nn-1)  = (1.+delta)*Cext*ni_ext_e*tau_e*Jee_ext^2;
        Sigma_ext(INIB) = (1.+delta)*Cext*ni_ext_i*tau_i*Jie_ext^2;
%     elseif strcmp(Opt,'DOIRON') 
%         Sigma_ext(1:nn-1)=(delta)^2*Cext*ni_ext*tau_e*Jee_ext/(1.15); % see normalization of Jee_ext in createparamsair.m
%         Sigma_ext(INIB) = (delta/2)^2*Cext*ni_ext*tau_i*Jie_ext/1.025;
    end
    

end

%-----------
% CUE
%-----------
% contributes to variance of Sigma_ext to inh neurons
if any(strcmp(fieldnames(Network),'cuemode'))
    cue_value=Network.cue_value;
    switch Network.cuemode
        case 'inh_mean'
            Mu_ext(end)=Mu_ext(end)*(1+cue_value);
        case 'exc_mean'
            Mu_ext(1:end-2)=Mu_ext(1:end-2)*(1+cue_value);
        case 'exc_gaussian'
            %--------------
            % GAUSSIAN CUE
            %--------------
            % function to be integrated later for quenched noise computation with
            % gaussian envelope
            Mu_extZ  = @(Z)Cext*ni_ext_e*tau_e*Jee_ext*cue_value*Z.*heaviside(1+cue_value*Z);
            Params.Mu_extZ=Mu_extZ;
            if nn>=3
                Params.quench_pops=1:nn-2;
            end
        case 'inh_gaussian'
            %--------------
            % GAUSSIAN CUE
            %--------------
            % function to be integrated later for quenched noise computation with
            % gaussian envelope
            Mu_extZ  = @(Z)Cext*ni_ext_i*tau_i*Jie_ext*cue_value*Z.*heaviside(1+cue_value*Z);
            Params.Mu_extZ=Mu_extZ;
            if nn>=3
                Params.quench_pops=nn;
            end
        case 'Chance'
            for i=1:nn-1
                Sigma_ext(i)  = (1.+delta)*Cext*ni_ext_e*tau_e*Jee_ext^2*cue_value;
            end
            Sigma_ext(INIB) = (1.+delta)*Cext*ni_ext_i*tau_i*Jie_ext^2*cue_value;     
    end
%         pcue=Network.pcue;
%         mu_cue=0; 
%         sigma_cue=0;
%         if numel(cue_value)~=numel(pcue)
%             error('ERROR in cue_mode: elements in cue_value different from pcue...\n');
%         end
%         switch Network.cuemode
%             case 'mixedcue'
%                 for icue=1:numel(cue_value)
%                     % current
%                     mu_cue=mu_cue+tau_e*Cext*ni_ext_e*Jee_ext*pcue(icue)*cue_value(icue);
%                     sigma_cue=sigma_cue+tau_e*Cext*(ni_ext_e*Jee_ext)^2*pcue(icue)*(2*ni_ext_e*cue_value(icue)+cue_value(icue)^2);
%                 end
%             case 'gaussian'
%                 sigma_cue=sigma_cue+tau_e*pcue*(Cext*ni_ext_e*Jee_ext*cue_value)^2;
%             case 'inh'   
%         end
%         Sigma_ext(1:end-1)=Sigma_ext(1:end-1)+sigma_cue;
%         Mu_ext(1:end-1)=Mu_ext(1:end-1)+mu_cue;
%     end
end

% if any(strcmp(fieldnames(Network),'cuemode'))
%     cue_value=Network.cue_value;
%     switch Network.cuemode
%         case 'inh'
%             Mu_ext(end)=Mu_ext(end)*(1+cue_value);
%         case 'up_exc'
%             Mu_ext(1:end-2)=Mu_ext(1:end-2)*(1+cue_value);
%         case 'down_exc'
%             Mu_ext(1:end-2)=Mu_ext(1:end-2)*(1+cue_value);
%         case 'gaussian'
%             %--------------
%             % GAUSSIAN CUE
%             %--------------
%             % function to be integrated later for quenched noise computation with
%             % gaussian envelope
%             Mu_extZ  = @(Z)Cext*ni_ext_e*tau_e*Jee_ext*cue_value*Z;
%             Params.Mu_extZ=Mu_extZ;
%             if nn>=3
%                 Params.quench_pops=1:nn-2;
%             end
%         case 'Chance'
%             for i=1:nn-1
%                 Sigma_ext(i)  = (1.+delta)*Cext*ni_ext_e*tau_e*Jee_ext^2*cue_value;
%             end
%             Sigma_ext(INIB) = (1.+delta)*Cext*ni_ext_i*tau_i*Jie_ext^2*cue_value;     
%     end
% %         pcue=Network.pcue;
% %         mu_cue=0; 
% %         sigma_cue=0;
% %         if numel(cue_value)~=numel(pcue)
% %             error('ERROR in cue_mode: elements in cue_value different from pcue...\n');
% %         end
% %         switch Network.cuemode
% %             case 'mixedcue'
% %                 for icue=1:numel(cue_value)
% %                     % current
% %                     mu_cue=mu_cue+tau_e*Cext*ni_ext_e*Jee_ext*pcue(icue)*cue_value(icue);
% %                     sigma_cue=sigma_cue+tau_e*Cext*(ni_ext_e*Jee_ext)^2*pcue(icue)*(2*ni_ext_e*cue_value(icue)+cue_value(icue)^2);
% %                 end
% %             case 'gaussian'
% %                 sigma_cue=sigma_cue+tau_e*pcue*(Cext*ni_ext_e*Jee_ext*cue_value)^2;
% %             case 'inh'   
% %         end
% %         Sigma_ext(1:end-1)=Sigma_ext(1:end-1)+sigma_cue;
% %         Mu_ext(1:end-1)=Mu_ext(1:end-1)+mu_cue;
% %     end
% end



% reset
for i=1:nn-3
    H(i) = He*theta_e*Jee;
end
H(INIB) = Hi.*theta_i*Jee;

% thresholds
for i=1:nn-3
    Theta(i) = theta_e*Jee;
end
Theta(INIB) = theta_i*Jee;

% membrane time
for i=1:nn-3
    Tau(i) = tau_e;
end
Tau(INIB) = tau_i;

% synaptic time
for i=1:nn-3
    Tausyn(i) = tausyn_e;
end
Tausyn(INIB) = tausyn_i;
if strcmp(Network.syn,'DoubleExp')
    Tausyn=Tausyn+tau_1;
end

% % IMPOSE RANDOM CONNECTIVITY BETWEEN CLUSTERS
% % connection probability between clusters
% if strcmp(Network.clust_syn,'random')
%     if abs(p+2-nn)>1e-6
%         error('\n # Random connectivity: of clusters p must be equal to n+2\n');
%     end
%     ProbConn=Network.clust_conn;
%     Bsmall=~(rand(p)<ProbConn); % 1 where connectivity matrix should be 0
%     for i=1:p
%         A(i,Bsmall(i,:))=0;
%         B(i,Bsmall(i,:))=0;
%     end
% end
    
    




Params.A=A;
Params.B=B;
Params.Mu_ext=Mu_ext; 
Params.Sigma_ext=Sigma_ext; 
Params.H=H' ;
Params.Tau=Tau';
Params.Tausyn=Tausyn';
Params.Theta=Theta';
Params.tau_arp=tau_arp;
% vector d (dimensions) is defined in _espo only if overlap>0; for the
% moment let's define it as a zeros:
Params.d=zeros(1,nn);
% synaptic dynamics
if strcmp(Network.syn,'DoubleExp') || strcmp(Network.syn,'SingleExp')
    Params.BS=1; % activates Brunel-Sergi thresholds
else
    Params.BS=0; % usual transfer function
end



% FILE INFO
Files=auxMFT.File_Info_MFT(nn,p,flg);
Params.Files=Files;
Params.Network=Network;

  