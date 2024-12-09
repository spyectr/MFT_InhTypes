%   /***************************************************
%    *                                                 *
%    *                  RATE_VS_Jfun                   *
%    *                                                 *
%    * Plots activities as function of J+/J.         *
%    * It stops when it finds an unstable point, or if Jplus reaches the Jstop *
%    * The starting point for activities is the spontaneous state with 4 pops found with Fix_Threshold with Jplus=1
%    * The routine has 2 modalities: selective (option -as) or spontaneous (option -a)
%    * Jstop is set with options -a<Jstop> or -as<Jstop>, otherwise it uses default value
%    *                                                 *
%    ***************************************************/

function flg=RATE_VS_Jfun(ni,Jstop,Params,STEP,Opt,flg)

% ni=Ni;

% global return_value NTRIAL TOLX TOLF screen contr overlap  nivsj_mode spur_mode stim_mode syn_dyn cue_mode

% Files
paramsfile=Params.Files.paramsfile;
file_eig=Params.Files.file_eig;
if (flg.screen); fprintf('Loading %s',paramsfile); end

load(paramsfile);
nn=numel(ni);
if 1%nn<50
    MAX_DATA=nn;
else
    MAX_DATA=50;
end
% previous version:
% STEP=0.1;

step = STEP;
% k=10;
flag=0;
stop_control = 1;

% variables to be saved to file
rates=[]; 
Jplus_store=[];
eigenvalues.Real=[];
eigenvalues.Imag=[];
eigenvectors=[]; % dim: 1st=pop #; 2nd=eigenvector #; 3rd=Jplus step

% Starting point:
Jplus = Jzero;
% create JplusArray
JplusArray=[];
while Jplus<=Jstop
    JplusArray=[JplusArray Jplus];
    step=auxMFT.StepJplus(Jplus,Opt,STEP,flg);
    Jplus =Jplus + step;
end
for Jplus=JplusArray
        % update Jplus in paramsfile to be loaded by function Learn
    save(paramsfile,'Jplus','-append');
    Params=auxMFT.Learn(Params,flg);
    % Investigate Selective Activity
    if (flg.nivsj_mode == 1)
        % NO OVERLAP CASE:
%         if (~overlap) % the kick is customized for each network by hand
            if ~isempty(strfind(Opt,'Small2'))
                JKick=1.1;
                Nkick=30;
%                 if flg.cue_mode==1    
%                     Nkick=40;
%                     acue=abs(Network.cue_value(1));
%                     if acue<=0.01
%                         JKick=10.22;
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
%                         JKick=11.25;
%                     end
%                 end
%                 if (flag < 3)% this means do it only 3 times (because it's time consuming)
                    if (Jplus>JKick)
                        ni(1) = Nkick; % ni(nn)=10.;
                        ni(2:end-4) = 2; % ni(nn)=10.;
                        ni(end-3) = 2; % ni(nn)=10.;
                        ni(end-2:end) = 11; % ni(nn)=10.;
                        ni0=ni;
                        fprintf('RATE_VS_Jfun: ni(1)=%0.03g Hz kick starts at Jplus>%0.03g \n',ni(1),JKick); % debug
                        [ni,flg]=auxMFT.goAB_simul(ni,5000,Params,flg);
                        if isempty(ni)
                            ni=ni0;
                        end
                        flag=flag+1;
                    end
%                 end
            end
%             if ~isempty(strfind(Opt,'Multifinal'))
% %                 JKick=15.84;
%                 JKick=7.85;
%                 Nkick=25;
%                 if (flag < 5)  % this means do it only 3 times (because it's time consuming)
%                     if (Jplus>JKick)
%                         ni(1) = Nkick; % ni(nn)=10.;
%                         ni(2:end-4) = 5; % ni(nn)=10.;
%                         ni(end-3) = 2; % ni(nn)=10.;
%                         ni(end-2:end) = 7; % ni(nn)=10.;
%                         ni0=ni;
%                         fprintf('RATE_VS_Jfun: ni(1)=%0.03g Hz kick starts at Jplus>%0.03g \n',ni(1),JKick); % debug
%                         ni=auxMFT.goAB_simul(ni,20000,Params);
%                         if isempty(ni)
%                             ni=ni0;
%                         end
%                         flag=flag+1;
%                     end
%                 end
%             end
%         end
    else
    end % if(flg.nivsj_mode == 1)

%     /*
    % Help mnewt search for low rate fixed points (usually not needed though; goAB_simul could be used as well)
    
    if flg.stim_mode
        [ni,flg]=auxMFT.goAB_simul(ni,2000,Params,flg);
    end%     if (flg.return_value == 0)  BYE
%     */

    % Search fixed point:
    fprintf('RATE_VS_Jfun: mnewt \n'); % debug
    ni0=ni;
    [ni,flg]=auxMFT.mnewt(flg.NTRIAL,ni0,Params,flg);
    if (flg.return_value == 0)
        fprintf('\n>>> Error in routine mnewt() called by RATE_VS_Jfun(): rerun goAB_simul\n');
        flg.return_value=1;
        [ni,flg]=auxMFT.goAB_simul(ni0,40000,Params,flg);
        if isempty(ni)
            ni=ni0;
            continue;
        else
%             ni=mnewt(flg.NTRIAL,ni0,flg.TOLX,flg.TOLF,Params);
        end
    end
    % store to save fixed points to file 'nivsj.mat':
    Jplus_store=[Jplus_store Jplus];
    rates=[rates; ni(1:MAX_DATA)];
    
    % info on flg.screen
    if (flg.screen); fprintf('\n    (Jplus = %g)\n',Jplus); else; fprintf('--- Jplus = %g\n',Jplus); end
    % handle unstable points:
%         ni=auxMFT.goAB_simul(ni0,40000,Params);
    [res, Real_part, Imag_part, EigV, flg,matrix]=auxMFT.Test_Stability(ni,Params,flg);
    % store to save
    eigenvalues.Real=[eigenvalues.Real Real_part];
    eigenvalues.Imag=[eigenvalues.Imag Imag_part];
    eigenvectors=cat(3,eigenvectors,EigV);
    %
    if ((~res) && (stop_control == 1)) && (flg.nivsj_mode == 1)
        fprintf('\n--- Unstable state occurred at J+ = %g\n',Jplus);
        
        % each time unstable point occurs, stop and ask for dynamics simulation:
        %         if (stop_control == 1)
        fprintf('\n    Do Dynamics Simulation.\n');
        ni =ni + 0.1;  % Little perturbation
        [ni,flg]=auxMFT.goAB_simul(ni,10000,Params,flg); % 0 := no max number of steps.
        [res, Real_part, Imag_part, EigV,flg]=auxMFT.Test_Stability(ni,Params,flg);
        [~,unind]=max(Real_part); % find largest positive eig
        unstabledir=EigV(:,unind)'; % corresponding eigenvector
        [ni,flg]=auxMFT.goAB_simul(ni+0.2*unstabledir,10000,Params,flg); % 0 := no max number of steps.
        [ni,flg]=auxMFT.mnewt(flg.NTRIAL,ni0,Params,flg);
        [res, Real_part, Imag_part, EigV, flg,matrix]=auxMFT.Test_Stability(ni,Params,flg);
%         res
%         Real_part
%         EigV
%         bla
    end % end if ~res stop_control
    if res
        fprintf('    stable fixed point:');
    else
        fprintf('    unstable fixed point:');
    end
    for j=1:nn; fprintf(' %g ',ni(j)); end
    fprintf('\n');
end % while(Jplus<=Jstop)

% Summary of parameters:
fprintf('\n--- Done.\n');
auxMFT.Summary(paramsfile);

%---------------
% SAVE TO FILE
%---------------
% file to store fixed points
file_name=Params.Files.FP;
% file to store eigenvalues and eigenvectors
file_eig=Params.Files.file_eig;
% fixed point rates
AllParameters=load(paramsfile);
save(file_name,'Jplus_store','rates','AllParameters');
% all eigenvalues
save(file_eig,'Jplus_store','eigenvalues','eigenvectors','AllParameters');



return;
