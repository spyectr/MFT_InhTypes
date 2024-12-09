% Finds Exc and Inh thresholds given firing rates in spontaneous global state
% INPUT theta (old)
% OUTPUT theta (new)
% Overwrites new theta's to current paramsfile
%
% from GLC C_overlap_thesis
% Luca Mazzucato February 2014
% Luca Mazzucato May 2020


function theta_out=fun_Fix_Weights_InhType(ni,control,params)

TOL_THRES=1e-9;
TOL_POINT=1e-9;
global DEBUG return_value screen 

% load parameters, including 
% load(paramsfile);
utils.v2struct(params);
if DEBUG
    disp('::: threshold.c in debugging mode \n');
end

% ------------------------------------------------------------------
%                       Starting conditions:                         
% ------------------------------------------------------------------

% To avoid problems:
% ni_e = 0.0;
% ni_i = 0.0;
n=numel(ni);
ni_e=ni(1);
ni_i=ni(2:n);
%
ni_e_wanted = ni(1);
ni_i_wanted = ni(2:n)';
if (control == 1)
    ni_ext_e = ni_e_wanted;
    ni_ext_i = ni_i_wanted;
end
% flags
flag_e = 1;   % also 0: it is (almost) the same
flag_i = zeros(1,n-1);

params=aux.In_Param_4pop_InhType([theta_e, theta_i],paramsfile);
% params.A
% params.B
% Risp=auxMFT.RISP(ni,params);
Mu=auxMFT.MU(ni,params); % RUN FUNCTION
Sigma2=auxMFT.SIGMA2(ni,params); % RUN FUNCTION
if return_value == 0 
    fprintf('\n>>> Error in MU_C or SIGMA_C called by routine Fix_Threshold\n');
    return; % exit from function
end

%   Initial thresholds:
theta_e = Mu(1)+3*sqrt(Sigma2(1));
theta_i = max(Mu(2:n)',zeros(size(Mu(2:n)')))+3*sqrt(Sigma2(2));

%   Initial steps:
step_e = sqrt(Sigma2(1))/10.;
step_i = sqrt(Sigma2(2:n))/10.;

if DEBUG 
    fprintf('::: Starting Fix_Threshold...\n');
    fprintf('::: Initial parameters: Theta     [e] %5.3f  [i] %5.3f\n',theta_e,theta_i(1));
    fprintf('::: Initial parameters: ni       [e] %5.3f  [i] %5.3f\n',ni_e,ni_i(1));
    fprintf('                        steps    [e] %5.3f  [i] %5.3f\n\n',step_e,step_i(1));
    fprintf('\n');
end

%------------------------------------------------------------------
%                                 Loop:                         
%------------------------------------------------------------------
k=1; % counting loops
while 1
    if k>20000 
        fprintf('\n>>> Error: too many steps (>%d) in routine fix_threshold \n',k-1);
        return_value = 0;
        return; % exit from function
    end

    previous_flag_e = flag_e;
    previous_flag_i = flag_i;
    params=aux.In_Param_4pop_InhType([theta_e, theta_i],paramsfile);
    
    
    
    Risp=auxMFT.RISP(ni,params);
    if return_value == 0 
        fprintf('\n>>> Error in RISP() called by routine Fix_Threshold()\n');
        return;
    end
    %
    ni_e = Risp(1);
    ni_i = Risp(2:n);

    if (ni_e-ni_e_wanted) > 0.0  
      flag_e = 1;
      theta_e = theta_e+step_e;
    else 
      flag_e = 0;
      theta_e =theta_e- step_e;
    end
    for i=1:n-1
        if (ni_i(i)-ni_i_wanted(i)) > 0.0
            flag_i(i) = 1;
            theta_i(i) = theta_i(i) + step_i(i);
        else
            flag_i(i) = 0;
            theta_i(i) = theta_i(i)- step_i(i);
        end
    end
    % 
    if ((abs(ni_e-ni_e_wanted) < TOL_POINT) && (sum(abs(ni_i-ni_i_wanted)) < TOL_POINT)) || (step_e< TOL_THRES && sum(step_i)< TOL_THRES) 
        if screen
            fprintf('Threshold fixed by routine Fix_Threshold:\n');
            fprintf('Theta_e = %g   Theta_i = ',theta_e);
            for i=1:n-1; fprintf('%g,',theta_i(i)); end
            fprintf('with parameters:\n');
            fprintf('Ni_Exc: %g    Ni_Inh: ',ni_e);
            for i=1:n-1; fprintf('%g,',ni_i(i)); end
            fprintf('Ni_ext_e: %g   Ni_ext_i: %g\n\n',ni_ext_e,ni_ext_i(1));
            fprintf('and tolerances:\n');
            fprintf('On thresholds: %g    On fixed point: %g \n\n', TOL_THRES,TOL_POINT);
        end
        theta_out =[theta_e, theta_i]; % output
%         if any(strcmp(fieldnames(params),'paramsfile'))
            save(paramsfile,'theta_e','theta_i','-append'); % overwrite new threholds to current paramsfile
%         end
        return; % exit function
    end 

    if DEBUG && rem(k,100)==0
        fprintf('step_e=%0.03g,step_i=',step_e);
        for i=1:n-1; fprintf('%0.03g,',step_i(i));end
        fprintf(' theta_e=%0.03g,theta_i=',theta_e);
        for i=1:n-1; fprintf('%0.03g,',theta_i(i)); end
        fprintf('ni_e=%0.03g,ni_i=',ni_e);
        for i=1:n-1; fprintf('%0.03g,',ni_i(i));end
        fprintf('\n');
    end

    if (flag_e - previous_flag_e ~= 0)   
        step_e = step_e/2;          % do change!
    end
    for i=1:n-1
        if (flag_i(i) - previous_flag_i(i) ~= 0)   
            step_i(i) = step_i(i)/2;
        end
    end
	
    k=k+1;
end

%







