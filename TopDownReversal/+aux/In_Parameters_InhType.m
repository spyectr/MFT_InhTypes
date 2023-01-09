%/********************************************************
% *                                                      *
% *       Initializations in nn populations case         *
% *       ------  Important: nn >= 3(*) --------         *
% *                                                      *
% * (*)It works also with nn=4 for in this case          *
% *    nn_sel=0 and therefore loops on selective         *
% *    populations are skipped.                          *
% *    See code for a better understanding.              *
% *    LM 2020                                           *
% ********************************************************/

function Params=In_Parameters_InhType(theta,paramsfile,flg)

% global DEBUG  nivsj_mode spur_mode stim_mode all_mode

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
Jsel_to_bg = Jminus; %% 1.;  %from Amit&Brunel, B99
Jcor = 0.; % from Amit&Brunel, B99
% % TO AVOID PROBLEMS:
x = 1.; % in Amit&Brunel 97 local exc currents have x while the external current 1-x

% % higher intra-cluster connectivity in DOIRON
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

%----------
% STIMULUS
%----------

%%%%%%%%%%%%%%%%%%%%%%%
% BELOW IS STILL WITHOUT INH TYPES
%%%%%%%%%%%%%%%%%%%%%%%

if flg.stim_mode
    A=zeros(nn); B=zeros(nn);
    % # of selective populations which are stimulated
    SelPop=Network.SelPop;
    % distribution of grouped populations
    % (ni_{+0}(1:nn_sel),ni_+,ni_-,ni_bg,ni_I)
    %-------------
    % A (currents)
    %-------------
    if nn<p+2
        nn_sel=nip0; % # of ungrouped stimulated populations
        nn_other_sel=SelPop-nn_sel; % # of group stimulated populations
        OTHER_SEL=nn_sel+1;
        NON_STIM=nn-2;
        nn_stim=p-SelPop;
    elseif nn==p+2
        nn_sel=SelPop-1;
        nn_other_sel=1;
        OTHER_SEL=SelPop;
        NON_STIM=SelPop+1:nn-2;
        nn_stim=1;
    end
    BACK_EXC =nn-1;
    INIB     =nn;
    %/*-----------------------------------------------
    % connections TO splitted selective populations:
    % -----------------------------------------------*/
    for i=1:nn_sel
    %     % sel -> itself:
        A(i,i) = tau_e*f*Cee*x*Jplus*Jee;
        if HiIntraCluster
            A(i,i) = A(i,i)*(peein/peeout);
        end
        for j=1:nn
            %       % sel(j) -> sel(i):
            if (j~=i) && (j<=nn_sel)  
              A(i,j) = tau_e*f*Cee*x*Jminus*Jee;
            end
            %       % all other sel -> sel:
            if (j==OTHER_SEL)  
              A(i,j) = tau_e*f*(nn_other_sel)*Cee*x*Jminus*Jee;
            end
            %       % non stimulated -> sel:
            if any(j==NON_STIM)  
              A(i,j) = tau_e*f*(nn_stim)*Cee*x*Jminus*Jee;
            end
            %       % bg -> sel:
            if (j==BACK_EXC)   
              A(i,j) = tau_e*(1.-p*f)*Cee*x*Jminus*Jee;
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
        % sel -> all other stimulated collected: 
        A(i,j) = tau_e*f*Cee*x*Jminus*Jee;
    end
    % all other stimulated collected -> itself:
    Temp=((nn_other_sel-1)*Jminus+Jplus);
    A(i,i)        = tau_e*x*f*Cee*Temp*Jee;
    if HiIntraCluster
        Temp1=((nn_other_sel-1)*Jminus+(peein/peeout)*Jplus);
        A(i,i) = A(i,i)*Temp1/Temp;
    end
    % non-stimulated pops -> all other stimulated collected: 
    A(i,NON_STIM) = tau_e*Cee*x*f*(nn_stim*Jminus)*Jee;
    % bg  -> all other sel collected: 
    A(i,BACK_EXC) = tau_e*(1.-f*p)*Cee*x*Jminus*Jee;
    % inib  -> all other sel collected:
    A(i,INIB)     = -tau_e*Cei*Jei;
    
    % /*--------------------------------------------------------------
    % connections TO p-SelPop NON-STIM COLLECTED selective populations:
    % --------------------------------------------------------------*/
    for i=NON_STIM;
        for j=1:nn_sel
            % sel -> all other stimulated collected: 
            A(i,j) = tau_e*f*Cee*x*Jminus*Jee;
        end
        % all other stimulated collected -> NON-STIM:
        A(i,OTHER_SEL)        = tau_e*Cee*x*f*((nn_other_sel)*Jminus)*Jee;
        % non-stimulated pops -> all other stimulated collected: 
        Temp=((nn_stim-1)*Jminus+Jplus);
        A(i,i)        = tau_e*x*f*Cee*Temp*Jee;
        if HiIntraCluster
            Temp1=((nn_stim-1)*Jminus+(peein/peeout)*Jplus);
            A(i,i) = A(i,i)*Temp1/Temp;
        end
        for j=NON_STIM
            if (j~=i)
              A(i,j) = tau_e*f*Cee*x*Jminus*Jee;
            end
        end
        % bg  -> all other sel collected: 
        A(i,BACK_EXC) = tau_e*(1.-f*p)*Cee*x*Jminus*Jee;
        % inib  -> all other sel collected:
        A(i,INIB)     = -tau_e*Cei*Jei;
    end
    
    % /*--------------------------------------------------------------
    % connections TO bg (=eccit. nonselective) population:
    % --------------------------------------------------------------*/
    i=BACK_EXC;
    for j=1:nn_sel
        % sel ->  bg: 
        A(i,j) = tau_e*f*Cee*x*Jsel_to_bg*Jee;
    end
    % all other sel (collected) ->  bg:
    A(i,OTHER_SEL) = tau_e*f*(nn_other_sel)*Cee*x*Jsel_to_bg*Jee;
    % all other NON-STIM (collected) ->  bg:
    A(i,NON_STIM) = tau_e*f*(nn_stim)*Cee*x*Jsel_to_bg*Jee;
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
    A(i,OTHER_SEL)= tau_i*f*(nn_other_sel)*Cie*x*Jie;
    % all other NON-STIM (collected) ->  bg:
    A(i,NON_STIM) = tau_e*f*(nn_stim)*Cee*x*Jsel_to_bg*Jee;
    % bg -> inib:
    A(i,BACK_EXC)  = tau_i*(1-p*f)*Cie*x*Jie;
    % inib -> itself:
    A(i,i)         = -tau_i*Cii*Jii; 

    if flg.DEBUG
      printf('Matrix A:\n');
      for i=1:nn
        for j=1:nn  
            printf(' %f ',A(i,j)); 
        end
        printf('\n');
      end
    end
    
    %-------------
    % B (variances)
    %-------------
    %/*-----------------------------------------------
    % connections TO splitted selective populations:
    % -----------------------------------------------*/
    for i=1:nn_sel
    %     % sel -> itself:
        B(i,i) = tau_e*f*Cee*x*(Jplus*Jee)^2;
        if HiIntraCluster
            B(i,i) = B(i,i)*(peein/peeout);
        end
        for j=1:nn
            %       % sel(j) -> sel(i):
            if (j~=i) && (j<=nn_sel)  
              B(i,j) = tau_e*f*Cee*x*(Jminus*Jee)^2;
            end
            %       % all other sel -> sel:
            if (j==OTHER_SEL)  
              B(i,j) = tau_e*f*(nn_other_sel)*Cee*x*(Jminus*Jee)^2;
            end
            %       % non stimulated -> sel:
            if any(j==NON_STIM)  
              B(i,j) = tau_e*f*(nn_stim)*Cee*x*(Jminus*Jee)^2;
            end
            %       % bg -> sel:
            if (j==BACK_EXC)   
              B(i,j) = tau_e*(1.-p*f)*Cee*x*(Jminus*Jee)^2;
            end
            %       % inib -> sel:
            if (j==INIB)       
              B(i,j) = tau_e*Cei*(Jei)^2; 
            end
        end
    end
    
    % /*--------------------------------------------------------------
    % connections TO other p-nn_sel COLLECTED selective populations:
    % --------------------------------------------------------------*/
    i=OTHER_SEL;
    for j=1:nn_sel
        % sel -> all other stimulated collected: 
        B(i,j) = tau_e*f*Cee*x*(Jminus*Jee)^2;
    end
    % all other stimulated collected -> itself:
    B(i,i)        = tau_e*Cee*x*f*((nn_other_sel-1)*Jminus^2+Jplus^2)*Jee^2;
    Temp=((nn_other_sel-1)*Jminus^2+Jplus^2);
    B(i,i)        = tau_e*x*f*Cee*Temp*Jee^2;
    if HiIntraCluster
        Temp1=((nn_other_sel-1)*Jminus^2+(peein/peeout)*Jplus^2);
        B(i,i) = B(i,i)*Temp1/Temp;
    end
    % non-stimulated pops -> all other stimulated collected: 
    B(i,NON_STIM) = tau_e*Cee*x*f*(nn_stim)*(Jminus*Jee)^2;
    % bg  -> all other sel collected: 
    B(i,BACK_EXC) = tau_e*(1.-f*p)*Cee*x*(Jminus*Jee)^2;
    % inib  -> all other sel collected:
    B(i,INIB)     = tau_e*Cei*(Jei)^2;
    
    % /*--------------------------------------------------------------
    % connections TO p-SelPop NON-STIM COLLECTED selective populations:
    % --------------------------------------------------------------*/
    for i=NON_STIM;
        for j=1:nn_sel
            % sel -> all other stimulated collected: 
            B(i,j) = tau_e*f*Cee*x*(Jminus*Jee)^2;
        end
        % all other stimulated collected -> NON-STIM:
        B(i,OTHER_SEL)        = tau_e*Cee*x*f*((nn_other_sel)*Jminus^2)*Jee^2;
        % non-stimulated pops -> all other stimulated collected: 
        Temp=((nn_stim-1)*Jminus^2+Jplus^2);
        B(i,i)        = tau_e*x*f*Cee*Temp*Jee^2;
        if HiIntraCluster
            Temp1=((nn_stim-1)*Jminus^2+(peein/peeout)*Jplus^2);
            B(i,i) = B(i,i)*Temp1/Temp;
        end
        for j=NON_STIM
            if (j~=i)
              B(i,j) = tau_e*f*Cee*x*(Jminus*Jee)^2;
            end
        end
        % bg  -> all other sel collected: 
        B(i,BACK_EXC) = tau_e*(1.-f*p)*Cee*x*(Jminus*Jee)^2;
        % inib  -> all other sel collected:
        B(i,INIB)     = tau_e*Cei*(Jei)^2;
    end
    
    % /*--------------------------------------------------------------
    % connections TO bg (=eccit. nonselective) population:
    % --------------------------------------------------------------*/
    i=BACK_EXC;
    for j=1:nn_sel
        % sel ->  bg: 
        B(i,j) = tau_e*f*Cee*x*(Jsel_to_bg*Jee)^2;
    end
    % all other sel (collected) ->  bg:
    B(i,OTHER_SEL) = tau_e*f*(nn_other_sel)*Cee*x*(Jsel_to_bg*Jee)^2;
    % all other NON-STIM (collected) ->  bg:
    B(i,NON_STIM) = tau_e*f*(nn_stim)*Cee*x*(Jsel_to_bg*Jee)^2;
    % bg  ->  itself:
    B(i,i)         = tau_e*(1-p*f)*Cee*x*Jee^2*(1+Jcor^2);
    % inib  ->  bg:
    B(i,INIB)      = tau_e*Cei*(Jei)^2;

    % /*--------------------------------------------------------------
    % connections TO inhibitory population:
    % --------------------------------------------------------------*/
    i=INIB;
    for j=1:nn_sel
        % sel -> inib:
        B(i,j) = tau_i*f*Cie*x*(Jie)^2;
    end
    % all other sel -> inib:
    B(i,OTHER_SEL)= tau_i*f*(nn_other_sel)*Cie*x*(Jie)^2;
    % all other NON-STIM (collected) ->  bg:
    B(i,NON_STIM) = tau_e*f*(nn_stim)*Cee*x*(Jsel_to_bg*Jee)^2;
    % bg -> inib:
    B(i,BACK_EXC)  = tau_i*(1-p*f)*Cie*x*(Jie)^2;
    % inib -> itself:
    B(i,i)         = tau_i*Cii*(Jii)^2; 
    
    %-----------------
    % EXTERNAL (with stimulus)
    %-----------------
    %
    % higher kick (box stimulus equivalent to thalamic)
    if exist('StimHeight','var')
        if isempty(StimHeight)
            height_stim=box_height(paramsfile);
        elseif ~isempty(StimHeight)
            height_stim=StimHeight;%box_height(paramsfile);
        end
    end
    if strcmp(Opt,'AB97')
        height_stim=0;
    end
    % CURRENT
    Mu_ext=zeros(1,nn);
    % stimulated
    Mu_ext(1:nn_sel+1)  = Cext*ni_ext_e*tau_e*Jee_ext*(1+height_stim);
    Mu_ext(nn_sel+2:nn-1)  = Cext*ni_ext_e*tau_e*Jee_ext;
    Mu_ext(INIB) = Cext*ni_ext_i*tau_i*Jie_ext;

    % VARIANCE (zero unless 'Poisson')
    Sigma_ext=zeros(1,nn);
    if strcmp(Network.stim,'Poisson') 
        Sigma_ext(1:nn_sel+1)  = (1.+delta)*Cext*ni_ext_e*tau_e*Jee_ext^2*(1+height_stim);
        Sigma_ext(nn_sel+2:nn-1)  = (1.+delta)*Cext*ni_ext_e*tau_e*Jee_ext^2;
        Sigma_ext(INIB) = (1.+delta)*Cext*ni_ext_i*tau_i*Jie_ext^2;
%     elseif strcmp(Opt,'DOIRON') 
%         Sigma_ext(1:nn-1)=(delta)^2*Cext*ni_ext*tau_e*Jee_ext/(1.15); % see normalization of Jee_ext in createparamsair.m
%         Sigma_ext(INIB) = (delta/2)^2*Cext*ni_ext*tau_i*Jie_ext/1.025;
    end

end



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




%-----------
% CUE
%-----------
% contributes to variance of Sigma_ext to inh neurons
cue_pops=cell(1,nn);
[cue_pops{1:nn-3}] = deal('Pyr');
[cue_pops{nn-2:nn}] = deal('PV','SST','VIP');
if flg.cue_mode
    cue_option=Network.cue.cue_option;
    cue_stat=Network.cue.cue_stat;
    cue_value=Network.cue.cue_value;
    Mu_extZ=cell(1,nn);
    for i=1:numel(cue_option)
        indpop=find(strcmp(cue_pops,cue_option{i}));
        if strcmp(cue_stat{i},'mean')
            Mu_ext(indpop)=Mu_ext(indpop)*(1+cue_value{i});
            % if SST and VIP have zero baseline current take the PV baseline current and replace (1+cue_value{i}) with just cue_value{i}
            indinh=find(strcmp({'PV','SST','VIP'},cue_option{i}));
            if ~isempty(indinh)
                if ni_ext_i(indinh)==0 && any(ni_ext_i(1))
                    Mu_ext(indpop)=Mu_ext(nn-2)*cue_value{i};
                elseif ni_ext_i(indinh)==0 && ni_ext_i(1)==0
                    fprintf('Error in InParameters_InhType in cue_mode: external current vanishes\n'); flg.return_value=0; 
                end
            end

        elseif strcmp(cue_stat{i},'gaussian')
            Params.quench_pops=indpop;
            if strcmp(cue_option{i},'Pyr')
                Mu_extZ0  = @(Z)Cext*ni_ext_e*Tau(indpop(1))*Jee_ext*cue_value{i}*Z.*heaviside(1+cue_value{i}*Z);
            else
                indinh=find(strcmp({'PV','SST','VIP'},cue_option{i}));
                Mu_extZ0  = @(Z)Cext*ni_ext_i(indinh)*tau_i(indinh)*Jie_ext(indinh)*cue_value{i}*Z;%.*heaviside(1+cue_value{i}*Z);
                % if SST and VIP have zero baseline current take the PV baseline current and replace (1+cue_value{i}) with just cue_value{i}
                if ~isempty(indinh)
                    if ni_ext_i(indinh)==0 && any(ni_ext_i(1))
                        Mu_extZ0  = @(Z)Cext*ni_ext_i(1)*tau_i(indinh)*Jie_ext(indinh)*cue_value{i}*Z;%.*heaviside(1+cue_value{i}*Z);
                    elseif ni_ext_i(indinh)==0 && ni_ext_i(1)==0
                        fprintf('Error in InParameters_InhType in cue_mode: external current vanishes\n'); flg.return_value=0;
                    end
                end
            end
            Mu_extZ{indpop(:)}=deal(Mu_extZ0);
        elseif strcmp(cue_stat{i},'noise')
            if strcmp(cue_option{i},'Pyr')
                Sigma_ext(indpop)  = (1.+delta)*tau_e*Cext*ni_ext_e*Jee_ext^2*cue_value{i};
            else
                indinh=find(strcmp({'PV','SST','VIP'},cue_option{i}));
                Sigma_ext(indpop) = (1.+delta)*Cext*tau_i(indinh).*ni_ext_i(indinh).*Jie_ext(indinh).^2*cue_value{i};
            end            
        end
    end
    Params.Mu_extZ=Mu_extZ;
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

  