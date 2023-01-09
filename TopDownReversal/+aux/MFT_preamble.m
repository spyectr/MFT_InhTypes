%----------------------
% PREAMBLE
%----------------------
% flags keep track of options and are passed to all functions where they
% are updated and passed back to main program
flg=struct('DEFAULT',[0 0 0 0 0 0 0 1],'ERROR',[0 0 0 0 0 0 1 0],'RATE_VS_J',[0 0 0 0 1 0 0 0],'DYNAMICS',[0 0 1 0 0 0 0 0],'flags',[0 0 0 0 0 0 0 0],...
    'DEBUG',0,...
    'return_value',1,...% Exit code
    'dyn_mode',[],...% To select dynamics
    'syn_dyn',1,... % if 1, use synaptic dynamics (La Camera, Neural Comp. 2004)
    'contr',0,...% To enable contrast
    'screen',0,...% To put informations to screen
    'mnewt_flag',1,...% To enable screen section in routine mnewt()
    'NTRIAL',5000,... % # of steps for dynamics in goAB_simul()
    'TOLX',1.e-9,...% tolerance for convergence
    'TOLF',1.e-9,...% tolerance for convergence
    'nivsj_mode',0,...% To enact selective mode in RATE_VS_J
    'spur_mode',0,...% to enact spurious mode
    'stim_mode',0,...% to enact stimulus
    'cue_mode',0,...% to enact perturbation
    'all_mode',0,...% to find all stable states
    'nfocus',0,...
    'STEP',0.05,...
    'paramsfile',[]); % file where all current parameters are saved and loaded by all functions

%----------------------
% FLAGS
%----------------------

NumStim=1;      % default # of stimuli
% flg.nfocus=0;
% ---------------------------------
% Parsing Command Line Arguments: 
%---------------------------------
%
% FILES
% File to save and update current value of parameters in params.mat, to be
% called by all functions
name = utils.getComputerName();
paramsfile=fullfile('DATA',sprintf('params%sMFT.mat',name));
flg.paramsfile=paramsfile;
paramStr = ['create_params'];
paramFun = str2func(paramStr);
Opt=options_main.Opt; % which network to run
% load network parameters
paramFun(Opt,flg.paramsfile);
load(flg.paramsfile);
flg.flags =bitor(flg.flags,flg.RATE_VS_J);
flg.mnewt_flag = 0;

flg.nivsj_mode = 0;% 'Enact' Spontaneous mode:
if any(strcmp(fieldnames(options_main),'a'))
    switch options_main.a
        case 'sel'; flg.nivsj_mode = 1;% 'Enact' Selective mode:
        case 'spurious'; flg.spur_mode = 1;% Enact spurious mode
        case 'all'; flg.all_mode = 1;% Enact spurious mode
        otherwise
    end
end

%----------
% STIMULUS MODE
%----------
if any(strcmp(fieldnames(options_main),'S'))
    if  options_main.S>0
        flg.stim_mode=1;
        NumStim=options_main.S;
    end
end
%----------
% CUE MODE
%----------
if options_main.cue_mode
    flg.cue_mode=1;
    Network.cue=options_main.cue;
%     Network.pcue=1; % fraction of neurons targeted by cue

end
% OPTIONAL ARGUMENT FOR SPECIAL CONTROLS
Jzero=1; Jstop=5;
if any(strcmp(fieldnames(options_main),'screen')); flg.screen = options_main.screen; end
if any(strcmp(fieldnames(options_main),'p')); p=options_main.p; end % number of patterns:
if any(strcmp(fieldnames(options_main),'n')); nn=options_main.n; end % number of populations:
if any(strcmp(fieldnames(options_main),'Jzero')); Jzero=options_main.Jzero; end % first J+ value
if any(strcmp(fieldnames(options_main),'Jstop')); Jstop=options_main.Jstop; end % last J+ value
if ~exist('ro','var'); ro=1; end

% CELL DENSITY
if any(strcmp(fieldnames(options_main),'InhDensity'))
    InhDensity=options_main.InhDensity; 
    N_i=(N/5)*InhDensity;
end
if any(strcmp(fieldnames(options_main),'extra')); extra=options_main.extra; end % first J+ value

save(flg.paramsfile,'ni_e','ni_i','ni_ext_i','ni_ext_e','tau_arp','tau_i','tau_e','theta_e','theta_i','delta','f',...
    'x','Jee','Jii','Jie','Jei','Jee_ext','Jie_ext','Jplus','Jminus','gam','He','Hi','Qpiu','Qmeno','ro','Jzero','N_e',...
    'N_i','pee_matrix','pie_matrix','pei_matrix','pii_matrix','pext_matrix','p','nn','nfocus','Opt',...
    'Sim','Network','Stimulus','tausyn_e','tausyn_i','tau_1','NumStim','extra');
