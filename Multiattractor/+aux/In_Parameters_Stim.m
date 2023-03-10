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

function Params=In_Parameters(theta,paramsfile)

global DEBUG  nivsj_mode spur_mode stim_mode all_mode

% load parameters
load(paramsfile);
% connectivity
Cee = N_e*p_ee; % 1600;         % 1600
Cie = N_e*p_ie; % 1600;         % 1600
Cii = N_i.*p_ii; % 400;         % 400
Cei = N_i.*p_ei; % 400;         % 400
Cext = (N_ext)*p_ext; % 1600;

% set thresholds
theta_e=theta(1);
theta_i=theta(2);

OTHER_SEL=nn-2;
BACK_ECC =nn-1;
INIB     =nn;

% nn_sel = number of populations selective to a single pattern, which we do not
% want to group with all the other p-nn_sel
nn_sel = nn-3;
%------------------------------------------------------------------
% Synapses from sel. neurons to bg neurons.
% AB97: Jminus.  
%------------------------------------------------------------------*/
Jminus = 1.-gam*f*(Jplus-1.);
if Jminus<0
    error('ERROR in create_params: Jminus<0. Adjust gam or f\n');
end
Jsel_to_bg = Jminus; %% 1.;  % B99
%------------------------------------------------------------------
% Corrections to synapses from  bg neurons to itself.
% AB97: 0.  
%------------------------------------------------------------------*/
Jcor = 0.;%%-0.04*Jplus; % 1.;  % B99
% % TO AVOID PROBLEMS:
x = 1.; % in AB97 local exc currents have x while the external current 1-x

% higher intra-cluster connectivity in DOIRON
HiIntraCluster=0;
if any(strcmp(fieldnames(Network),'Ree')) && strcmp(Opt,'DOIRON')
    HiIntraCluster=1;
    pee=Cee/N_e;
    REE=Network.Ree;
    denom=REE*(N_e/p-1)+(N_e/p)*(p-1);
    peeout=pee*(N_e-1)/denom; % inter-cluster connectivity
    peein=REE*peeout;  % intra-cluster connectivity
    Cee=peeout*N_e;
end

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
        if (j==BACK_ECC)   
%           A(i,j) = tau_e*(1.-p*f)*Cee*x*Jminus*Jee;
          A(i,j) = tau_e*(1.-p*f)*Cee*x*Jsel_to_bg*Jee;
        end

        %       % inib -> sel:
        if (j==INIB)       
          A(i,j) = -tau_e*Cei*Jei; 
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
A(i,BACK_ECC) = tau_e*(1.-f*p)*Cee*x*Jsel_to_bg*Jee;
% inib  -> all other sel collected:
A(i,INIB)     = -tau_e*Cei*Jei;   

% /*--------------------------------------------------------------
% connections TO bg (=eccit. nonselective) population:
% --------------------------------------------------------------*/
i=BACK_ECC;
for j=1:nn_sel
    % sel ->  bg: 
    A(i,j) = tau_e*f*Cee*x*Jsel_to_bg*Jee;
end
% all other sel (collected) ->  bg:
A(i,OTHER_SEL) = tau_e*f*(p-nn_sel)*Cee*x*Jsel_to_bg*Jee;
% bg  ->  itself:
A(i,i)         = tau_e*(1-p*f)*Cee*x*Jee*(1+Jcor);
% inib  ->  bg:
A(i,INIB)      = -tau_e*Cei*Jei;

% /*--------------------------------------------------------------
% connections TO inhibitory population:
% --------------------------------------------------------------*/
i=INIB;
for j=1:nn_sel
    % sel -> inib:
    A(i,j) = tau_i*f*Cie*x*Jie;
end
% all other sel -> inib:
A(i,OTHER_SEL)= tau_i*f*(p-nn_sel)*Cie*x*Jie;
% bg -> inib:
A(i,BACK_ECC)  = tau_i*(1-p*f)*Cie*x*Jie;
% inib -> itself:
A(i,i)         = -tau_i*Cii*Jii; 

if DEBUG
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
        if (j==BACK_ECC)   
          B(i,j) = Jee*Jee*tau_e*(1-p*f)*Cee*x*Jminus*Jminus;
        end
        % inib -> sel:
        if (j==INIB)       
          B(i,j) = tau_e*Cei*Jei*Jei;
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
B(i,BACK_ECC) = Jee*Jee*tau_e*(1.-f*p)*Cee*x*Jminus*Jminus;
% inib  -> all other sel collected:
B(i,INIB)     = tau_e*Cei*Jei*Jei;

%   /*--------------------------------------------------------------
%     connections TO bg (=eccit. nonselective) population:
%     --------------------------------------------------------------*/
i=BACK_ECC;
for j=1:nn_sel
    % sel ->  bg: 
    B(i,j) = tau_e*f*Cee*x*Jminus*Jee*Jminus*Jee;
end
% all other sel (collected) ->  bg:
B(i,OTHER_SEL) = tau_e*f*(p-nn_sel)*Cee*x*Jminus*Jee*Jminus*Jee;
% bg  ->  itself:
B(i,i)         = tau_e*(1-p*f)*Cee*x*Jee*Jee*(1+Jcor*Jcor);
% inib  ->  bg:
B(i,INIB)      = tau_e*Cei*Jei*Jei;


% /*--------------------------------------------------------------
% connections TO inhibitory  population:
% --------------------------------------------------------------*/
i=INIB;
for j=1:nn_sel
    % sel -> inib:
    B(i,j) = tau_i*f*Cie*x*Jie*Jie; 
end
% all other sel -> inib: 
B(i,OTHER_SEL) = tau_i*f*(p-nn_sel)*Cie*x*Jie*Jie;
% bg -> inib:
B(i,BACK_ECC)  = tau_i*(1-p*f)*Cie*x*Jie*Jie;
% inib -> itself:
B(i,i)         = tau_i*Cii*Jii*Jii;
if (delta ~= 0) 
    for i=1:nn
        for j=1:nn 
            B(i,j) =B(i,j) *( 1.+delta);
        end
    end
end

if DEBUG
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
for i=1:nn-1
    Mu_ext(i)  = Cext*ni_ext*tau_e*Jee_ext;
end
Mu_ext(INIB) = Cext*ni_ext*tau_i*Jie_ext;

% external currents: variance
Sigma_ext=zeros(1,nn);
if strcmp(Network.stim,'Poisson') 
    for i=1:nn-1
        Sigma_ext(i)  = (1.+delta)*Cext*ni_ext*tau_e*Jee_ext^2;
    end
    Sigma_ext(INIB) = (1.+delta)*Cext*ni_ext*tau_i*Jie_ext^2;
elseif strcmp(Opt,'DOIRON') 
    Sigma_ext(1:nn-1)=(delta)^2*Cext*ni_ext*tau_e*Jee_ext/(1.15); % see normalization of Jee_ext in createparamsair.m
    Sigma_ext(INIB) = (delta/2)^2*Cext*ni_ext*tau_i*Jie_ext/1.025;
end

%----------
% STIMULUS
%----------
if stim_mode
    % # of selective populations
    SelPop=Network.SelPop;
    p_stim=0.5; % stimulus connectivity
    if SelPop>p
        error('# of selective pops > # of clusters!');
    end
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
    Mu_ext(1:SelPop)  = Mu_ext(1:SelPop)*(1+p_stim*height_stim);
end

%-----------
% MIXED CUE
%-----------
% contributes to variance of Sigma_ext to exc neurons
% cue mode; does not change mean ext current
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
%     if ~strcmp(Network.cuemode,'none')
%         cue_value=Network.cue_value;
%         pcue=Network.pcue;
%         mu_cue=0; 
%         sigma_cue=0;
%         if numel(cue_value)~=numel(pcue)
%             error('ERROR in cue_mode: elements in cue_value different from pcue...\n');
%         end
%         if strcmp(Network.cuemode,'mixedcue')
%             for icue=1:numel(cue_value)
%                 % current
%                 mu_cue=mu_cue+tau_e*Cext*ni_ext*Jee_ext*pcue(icue)*cue_value(icue);
%                 sigma_cue=sigma_cue+tau_e*Cext*(ni_ext*Jee_ext)^2*pcue(icue)*(2*ni_ext*cue_value(icue)+cue_value(icue)^2);
%             end
%         elseif strcmp(Network.cuemode,'gaussian')
%             sigma_cue=sigma_cue+tau_e*pcue*(Cext*ni_ext*Jee_ext*cue_value)^2;
%         end
%             
%         Sigma_ext(1:end-1)=Sigma_ext(1:end-1)+sigma_cue;
%         Mu_ext(1:end-1)=Mu_ext(1:end-1)+mu_cue;
%     end
% end



% reset
for i=1:nn-1
    H(i) = He*theta_e*Jee;
end
H(INIB) = Hi*theta_i*Jee;

% thresholds
for i=1:nn-1
    Theta(i) = theta_e;
end
Theta(INIB) = theta_i;

% membrane time
for i=1:nn-1
    Tau(i) = tau_e;
end
Tau(INIB) = tau_i;

% synaptic time
for i=1:nn-1
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
    
    




Params=[]; % output structure 
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
Files=File_Info_MFT(nn,p,paramsfile);
Params.Files=Files;

  