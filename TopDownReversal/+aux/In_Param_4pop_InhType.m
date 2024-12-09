% initializing parameters for Fix_Threshold
% loads parameters from paramsfile except thresholds
% INPUT theta=thresholds
% OUTPUT Params structure

function Params=In_Param_4pop_InhType(theta,paramsfile,flg)

% global  nivsj_mode spur_mode

%Jie=gi*Jii;
%Jei=ge*Jee;

n=numel(theta);
% load parameters
load(paramsfile);

% connectivity

    Cee = N_e*pee_matrix; % 1600;         % 1600
    Cie = N_e*pie_matrix; % 1600;         % 1600
    Cii = N_i.*pii_matrix; % 400;         % 400
    Cei = N_i.*pei_matrix; % 400;         % 400
    Cext = (N_e)*pext_matrix; % 1600; 



% external currents

% set thresholds
    Theta(1) = theta(1)*Jee;
    Theta(2:n) = theta(2:n)*Jee;  

% TO AVOID PROBLEMS:
x = 1.;

A(1,1) = tau_e*Cee*x*Jee;
A(1,2:n) = -tau_e*Cei.*Jei;
A(2:n,1) = x*tau_i'.*Cie.*Jie;
A(2:n,2:n) = -repmat(tau_i',1,3).*Cii.*Jii;

B(1,1) = tau_e*Cee*x*Jee^2; 
B(1,2:n) = tau_e*Cei.*Jei.^2;
B(2:n,1) = x*tau_i'.*Cie.*Jie.^2;
B(2:n,2:n) = repmat(tau_i',1,3).*Cii.*Jii.^2;

if (delta ~= 0) 
    for i=1:n
        for j=1:n B(i,j) = B(i,j) *(1.+delta); end
    end
end

  % if (x!=1)  Cext = C*(1-x);

Mu_ext(1) = Cext*ni_ext_e*tau_e*Jee_ext;
Mu_ext(2:n) = Cext*ni_ext_i.*tau_i.*Jie_ext;

Sigma_ext(1) = (1.+delta)*Cext*ni_ext_e*tau_e*Jee_ext^2;
Sigma_ext(2:n) = (1.+delta)*Cext*ni_ext_i.*tau_i.*Jie_ext.^2;
if strcmp(Network.stim.input,'Const')
    Sigma_ext=zeros(1,n);
end

H(1) = He*theta(1)*Jee;
H(2:n) = Hi.*theta(2:n)*Jee;

Tau(1) = tau_e;
Tau(2:n) = tau_i;

Tausyn(1) = tausyn_e;
Tausyn(2:n) = tausyn_i; 

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

