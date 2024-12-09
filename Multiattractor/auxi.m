classdef auxi
    methods(Static)
        
        
        %% Creates parameter for network with E and I clusters
        % from ClustersOption='EI'
        % Luca Mazzucato December 2020
        

        
        %% function [J, params]=SynWeights(params)
        %
        % OUTPUT
        %       J                =synaptic matrix
        %       params.popsize      =array of dim # pops, how many exc neurons in each
        %                         pop
        %       params.clustermatrix=matrix of dimension # pops x # clusters, each row shows
        %                   which clusters belong to that pop (1s and 0s)
        %
        % Luca Mazzucato December 2020
        function [J, params]=fun_SynWeights(paramsfile,Opt)
            
            % LOAD PARAMETERS
            params=load(paramsfile);
            utils.v2struct(params);
            Network=params.Network;
            utils.v2struct(Network);
            Sim=params.Sim;
            %-----------------------
            % PARAMETERS VALUES
            %-----------------------
            numfig=1;
            %             Next=N_e; % external units
            % CLUSTERS
            Q=p; % number of clusters
            %-----------------------
            % SYNAPTIC WEIGHTS
            %-----------------------
            % WEIGHTS
            % depression
            %             Jminus = (1-Jplus*f)/(1-f) ;
            params.Jminus = Jminus;
            %
            jee=Jee;
            jee_out=Jminus*Jee; % intra-cluster potentiation
            jee_in=Jplus*Jee; % inter-cluster depression
            jei=-Jei;
            jie=Jie;
            jii=-Jii;
            % connection probability
            p_eeout=p_ee;
            p_eein=p_ee;
            p_eiout=p_ei;
            p_eiin=p_ei;
            p_ieout=p_ie;
            p_iein=p_ie;
            p_iiout=p_ii;
            p_iiin=p_ii;
            %
            fprintf('  --- Jplus=%0.03g, Jminus=%0.03g\n',Jplus,Jminus);
            %----------------------------
            % SYNAPTIC MATRIX
            %----------------------------
            % generate a distribution of synaptic weights with mean J and
            % variance delta^2 J^2
            p_eeout=p_ee;
            p_eein=p_ee;
            % check #clusters and coding level are consistent
            NcUnits=round(f*N_e);    %  number of Exc units per cluster
            fprintf('  --- Synaptic weights: %d units/cluster \n',NcUnits);
            Numbg=round(N_e*(1-f*p)); % number of background (i.e. non-selective) Exc units
            fprintf('  --- fraction of bg Exc units: %0.03g',Numbg/N_e);
            jee_in=jee_in*ones(1,Q);
            switch Network.clust
                case 'hom'
                    popsize=repmat(NcUnits,Q,1); % size of each Exc cluster
                case 'het'
                    Nc=[];
                    clust_std=Network.clust_std;
                    while (sum(Nc)-(N_e-Numbg))~=0 || any(Nc<0)
                        Nc=round(NcUnits+(NcUnits*clust_std)*randn(Q,1));
                    end
                    popsize=Nc; % array of cluster sizes
                    if any(sum(popsize)-(N_e-Numbg))
                        fprintf('\n---ERROR: Heterogeneous clusters: Problem with cluster numbers\n');
                    end
                    fprintf('\n  --- het cluster sizes->rescale Jplus in each cluster. JEE+:');
                    jee_in=jee_in*mean(popsize)./popsize';
                    disp(jee_in);
            end
            cusumNcE=[0 cumsum(popsize)'];
            % background units (if present), if not, override in next line
            JEE=(jee*(ones(N_e)+delta*randn(N_e,N_e))).*(rand([N_e,N_e])<p_eeout);
            
            % weights involving inhibitory neurons (including InhTypes)
            JEI=[]; JIE=[]; JII=[];
            for i=1:numel(N_i)
                JEI=[JEI, (jei(i)*(ones(N_e,N_i(i))+deltaEI*randn(N_e,N_i(i)))).*(rand([N_e,N_i(i)])<p_ei(i))];
                JIE=[JIE; (jie(i)*(ones(N_i(i),N_e)+deltaIE*randn(N_i(i),N_e))).*(rand([N_i(i),N_e])<p_ie(i))];
                JIItemp=[];
                for j=1:numel(N_i)
                    JIItemp=[JIItemp,(jii(i,j)*(ones(N_i(i),N_i(j))+delta*randn(N_i(i),N_i(j)))).*(rand([N_i(i),N_i(j)])<p_ii(i,j))];
                end
                JII=[JII;JIItemp];
            end
            
            if strcmp(Network.clust,'het') || strcmp(Network.clust,'hom')
                % clustered units: inter-cluster weights
                JEE(1:cusumNcE(Q+1),1:cusumNcE(Q+1))=...
                    (jee_out*(ones(cusumNcE(Q+1))+delta*randn(cusumNcE(Q+1),...
                    cusumNcE(Q+1)))).*(rand([cusumNcE(Q+1),cusumNcE(Q+1)])<p_eeout); % inter-cluster weights
                for clu=2:Q+1 % intra-cluster higher weights
                    JEE(1+cusumNcE(clu-1):cusumNcE(clu),1+cusumNcE(clu-1):cusumNcE(clu))=...
                        (jee_in(clu-1)*(ones(popsize(clu-1))+delta*randn(popsize(clu-1),popsize(clu-1)))).*...
                        (rand([popsize(clu-1),popsize(clu-1)])<p_eein);
                end
            end
            clustermatrix=eye(Q);
            
            if strcmp(Opt,'MarcelPost_disinh')
                % use weights at the end of training
                filename=fullfile('marcel_parameters','marcel.mat');
                dataload=load(filename,'data');
                names={dataload.data(:).names};
                indEpoch=find(dataload.data(10).disinh==dataload.data(1).endTraining); % find epoch of end training
                jee_in=dataload.data(find(strcmp(names,'Intraassembly'))).disinh(indEpoch);
                jee_in=dataload.data(find(strcmp(names,'Intraassembly'))).disinh(indEpoch);
                bla
            end
            
            if strcmp(Opt,'EI')
                % INHIBITORY CLUSTERS
                if any(strcmp(clusters,'EI'))
                    %     JminusEI = 1.-gam*fI*(JplusEI-1.);
                    JplusEI = 1/(1/p+(1-1/p)/factorEI);
                    JminusEI = JplusEI/factorEI;
                    params.JminusEI = JminusEI;
                    jei_out=-JminusEI*Jei; % intra-cluster
                    jei_in=-JplusEI*Jei; % inter-cluster
                    fprintf('JplusEI=%0.03g, JminusEI=%0.03g\n',JplusEI,JminusEI);
                end
                if any(strcmp(clusters,'IE'))
                    %     JminusIE = 1.-gam*fI*(JplusIE-1.);
                    JplusIE = 1/(1/p+(1-1/p)/factorIE);
                    JminusIE = JplusIE/factorIE;
                    params.JminusIE = JminusIE;
                    jie_out=JminusIE*Jie; % intra-cluster
                    jie_in=JplusIE*Jie; % inter-cluster
                    fprintf('JplusIE=%0.03g, JminusIE=%0.03g\n',JplusIE,JminusIE);
                end
                if any(strcmp(clusters,'II'))
                    JplusII=factorII;
                    JminusII = 1.-gam*fI*(JplusII-1.);
                    params.JminusII = JminusII;
                    jii_out=-JminusII*Jii; % intra-cluster
                    jii_in=-JplusII*Jii; % inter-cluster
                    fprintf('JplusII=%0.03g, JminusII=%0.03g\n',JplusII,JminusII);
                end
                %-------------
                % EI CLUSTERS
                %-------------
                % check #clusters and coding level are consistent
                if any([strcmp(clusters,'EI'), strcmp(clusters,'EI'), strcmp(clusters,'II')])
                    NcUnits=round(fI*N_i);    %  number of Exc units per cluster
                    fprintf('  --- Synaptic weights: %d units/cluster \n',NcUnits);
                    Numbg=round(N_i*(1-fI*p)); % number of background (i.e. non-selective) Exc units
                    fprintf('  --- fraction of bg Inh units: %0.03g',Numbg/N_i);
                    switch Network.clustEI
                        case 'hom'
                            popsizeI=repmat(NcUnits,Q,1);
                        case 'het'
                            Nc=[];
                            clust_std=Network.clust_std;
                            while (sum(Nc)-(N_i-Numbg))~=0 || any(Nc<0)
                                Nc=round(NcUnits+(NcUnits*clust_std)*randn(Q,1));
                            end
                            popsizeI=Nc; % array of cluster sizes
                            if any(sum(popsizeI)-(N_i-Numbg))
                                fprintf('\n---ERROR: Heterogeneous clusters: Problem with cluster numbers\n');
                            end
                    end
                    cusumNcI=[0 cumsum(popsizeI)'];
                    
                    %-------------
                    % EI weights
                    %-------------
                    if any(strcmp(Network.clusters,'EI'))
                        % background units (if present), if not, override in next line
                        if strcmp(Network.clustEI,'het') || strcmp(Network.clustEI,'hom')
                            % clustered units: inter-cluster weights
                            JEI(1:cusumNcE(Q+1),1:cusumNcI(Q+1))=...
                                (jei_out*(ones(cusumNcE(Q+1),cusumNcI(Q+1))+deltaEI*randn(cusumNcE(Q+1),...
                                cusumNcI(Q+1)))).*(rand([cusumNcE(Q+1),cusumNcI(Q+1)])<p_eiout); % inter-cluster weights
                            for clu=2:Q+1 % intra-cluster higher weights
                                JEI(1+cusumNcE(clu-1):cusumNcE(clu),1+cusumNcI(clu-1):cusumNcI(clu))=...
                                    (jei_in*(ones(popsize(clu-1),popsizeI(clu-1))+deltaEI*randn(popsize(clu-1),popsizeI(clu-1)))).*...
                                    (rand([popsize(clu-1),popsizeI(clu-1)])<p_eiin);
                            end
                        end
                    end
                    
                    %-------------
                    % IE weights
                    %-------------
                    if any(strcmp(Network.clusters,'IE'))
                        % background units (if present), if not, override in next line
                        if strcmp(Network.clustIE,'het') || strcmp(Network.clustIE,'hom')
                            % clustered units: inter-cluster weights
                            JIE(1:cusumNcI(Q+1),1:cusumNcE(Q+1))=...
                                (jie_out*(ones(cusumNcI(Q+1),cusumNcE(Q+1))+deltaIE*randn(cusumNcI(Q+1),...
                                cusumNcE(Q+1)))).*(rand([cusumNcI(Q+1),cusumNcE(Q+1)])<p_ieout); % inter-cluster weights
                            for clu=2:Q+1 % intra-cluster higher weights
                                JIE(1+cusumNcI(clu-1):cusumNcI(clu),1+cusumNcE(clu-1):cusumNcE(clu))=...
                                    (jie_in*(ones(popsizeI(clu-1),popsize(clu-1))+deltaIE*randn(popsizeI(clu-1),popsize(clu-1)))).*...
                                    (rand([popsizeI(clu-1),popsize(clu-1)])<p_iein);
                            end
                        end
                    end
                    
                    %-------------
                    % II weights
                    %-------------
                    if any(strcmp(Network.clusters,'II'))
                        % background units (if present), if not, override in next line
                        if strcmp(Network.clustII,'het') || strcmp(Network.clustII,'hom')
                            % clustered units: inter-cluster weights
                            JII(1:cusumNcI(Q+1),1:cusumNcI(Q+1))=...
                                (jii_out*(ones(cusumNcI(Q+1),cusumNcI(Q+1))+deltaII*randn(cusumNcI(Q+1),...
                                cusumNcI(Q+1)))).*(rand([cusumNcI(Q+1),cusumNcI(Q+1)])<p_iiout); % inter-cluster weights
                            for clu=2:Q+1 % intra-cluster higher weights
                                JII(1+cusumNcI(clu-1):cusumNcI(clu),1+cusumNcI(clu-1):cusumNcI(clu))=...
                                    (jii_in*(ones(popsizeI(clu-1),popsizeI(clu-1))+deltaII*randn(popsizeI(clu-1),popsizeI(clu-1)))).*...
                                    (rand([popsizeI(clu-1),popsizeI(clu-1)])<p_iiin);
                            end
                        end
                    end
                    params.popsizeI=popsizeI;
                end
            end
            
            JEI(JEI>0)=0;
            JIE(JIE<0)=0;
            JII(JII>0)=0;
            JEE(JEE<0)=0;
            J=[JEE JEI; JIE JII];
            J=J-diag(diag(J)); % eliminate self-couplings
            fprintf('  --- New synaptic weights set: done...\n');
            InhNames={'PV','SST','VIP'};
            if numel(jie)>1
                fprintf('      Overall: Jee=%g\n',jee);
                for i=1:numel(N_i)
                    fprintf('      J[%s,E]=%g -- J[E,%s]=%g \n',InhNames{i},jie(i),InhNames{i},jei(i));
                end
                for i=1:numel(N_i); for j=1:numel(N_i)
                        fprintf('      J[%s,%s]=%g  ',InhNames{i},InhNames{j},jii(i,j));
                    end; fprintf('\n');
                end
            else
                fprintf('      Overall: Jee=%g -- Jie=%g -- Jei=%g -- Jii=%g \n',jee,jie,jei,jii);
            end
            fprintf('      Var[J]=(Jx%0.03g)^2\n',delta);
            
            params.popsize=popsize;
            params.clustermatrix=clustermatrix;
            
            if strcmp(Network.clust,'het')
                for i=1:Q
                    a(i)=sum(popsize(clustermatrix(:,i)>0));
                end
                fprintf('  --- clusters size: mean=%0.03g neurons/cluster, sd/mean=%0.03g\n',mean(a),std(a)/mean(a));
            end
            fprintf('  --- fraction of bg Exc units: %0.03g\n',Numbg/N_e);
            
%             %--------------
%             % POISSON INPUT
%             %--------------
%             % generate weights for external neurons
%             if strcmp(Stimulus.input,'Poisson')
%                 JExt=[];
%                 JExt=[JExt; Jee_ext*(rand([N_e,N_ext])<p_ext)];
%                 for i=1:numel(Jie_ext)
%                     JExt=[JExt; Jie_ext(i)*(rand([N_i(i),N_ext])<p_ext)];
%                 end
%                 save(fullfile(params.savedir,'Ext_weights.mat'),'JExt');
%             end
            
            
            %-------------------
            % PLOT weight matrix
            %-------------------
            figure(1); clf;
            subplot(2,1,1);
            colormap(utils.redblue); %xy=J; fun_colormapLim;
            imagesc(J);
            utils.figset(gca,'neurons','neurons','weights',10);
            colorbar;
            subplot(2,1,2);
            lambda=eig(J);
            plot(real(lambda),imag(lambda),'.');
            utils.figset(gca,'Re(\lambda)','Im(\lambda)','eig(weights)',10);
            saveas(gcf,fullfile(params.savedir,'weights.pdf'),'pdf');
            
            %             % spike thresholds for each population
            if any(strcmp(fieldnames(Network),'ThresholdOpt'))
                if Network.ThresholdOpt
                    flg.return_value=1; flg.screen=1;
                    Ni=[params.ni_e,params.ni_i];
                    [theta,~]=auxMFT.fun_Fix_Threshold_InhType(Ni,params,flg);
                    params.theta_e=theta(1); params.theta_i=theta(2:end);
                    fprintf('Spike thresholds Theta calculated for each population and stored in params\n');
                end
            end
%             params.theta_e=3.9; params.theta_i=[4,4,4];
            Theta=[params.theta_e, params.theta_i];
            fprintf('Spike thresholds: theta_e=%0.03g,theta_i=',Theta(1)');
            for i=1:numel(theta_i); fprintf('%0.03g,',Theta(i+1)); end
            params.Theta=Theta;
        end
        
        
        
        %% function [J, params]=SynWeights(params)
        %
        % OUTPUT
        %       J                =synaptic matrix
        %       params.popsize      =array of dim # pops, how many exc neurons in each
        %                         pop
        %       params.clustermatrix=matrix of dimension # pops x # clusters, each row shows
        %                   which clusters belong to that pop (1s and 0s)
        %
        % Luca Mazzucato December 2020
        function [J, params]=fun_SynWeights_Marcel(paramsfile,Opt)
            
            % LOAD PARAMETERS
            params=load(paramsfile);
            utils.v2struct(params);
            Network=params.Network;
            utils.v2struct(Network);
            Sim=params.Sim;
            %-----------------------
            % PARAMETERS VALUES
            %-----------------------
            numfig=1;
            %             Next=N_e; % external units
            % CLUSTERS
            Q=p; % number of clusters
            %-----------------------
            % SYNAPTIC WEIGHTS
            %-----------------------
            % WEIGHTS
            % depression
            %             Jminus = (1-Jplus*f)/(1-f) ;
            params.Jminus = Jminus;
            %
            jee=Jee;
            jee_out=Jminus*Jee; % intra-cluster potentiation
            jee_in=Jplus*Jee; % inter-cluster depression
            jei=-Jei;
            jii=-Jii;
            jie=Jie;
            %
            jei_in=jei;
            jei_out=jei;
            jee_bg_bg=jee;
            jee_in_bg=jee;
            jee_bg_in=jee;
            % connection probability
            p_eeout=p_ee;
            p_eein=p_ee;
            p_eiout=p_ei;
            p_eiin=p_ei;
            p_ieout=p_ie;
            p_iein=p_ie;
            p_iiout=p_ii;
            p_iiin=p_ii;
            %
            fprintf('  --- Jplus=%0.03g, Jminus=%0.03g\n',Jplus,Jminus);
            %----------------------------
            % SYNAPTIC MATRIX
            %----------------------------
            % generate a distribution of synaptic weights with mean J and
            % variance delta^2 J^2
            p_eeout=p_ee;
            p_eein=p_ee;
            
            
            if strcmp(Opt,'MarcelPost_disinh')
                % use weights at the end of training
                filename=fullfile('marcel_parameters','marcel.mat');
                dataload=load(filename,'data');
                names={dataload.data(:).names};
                indEpoch=find(dataload.data(10).disinh==dataload.data(1).endTraining); % find epoch of end training
%                 indEpoch=numel(dataload.data(10).disinh); % find epoch of end training
                jee_in=dataload.data(find(strcmp(names,'Intraassembly'))).disinh(indEpoch);
                jee_out=dataload.data(find(strcmp(names,'Interassembly'))).disinh(indEpoch);
                jee_bg_in=dataload.data(find(strcmp(names,'Assembly to bkg'))).disinh(indEpoch);
                jee_in_bg=dataload.data(find(strcmp(names,'Bkg to assembly'))).disinh(indEpoch);
                jee_bg_bg=dataload.data(find(strcmp(names,'Bkg to bkg'))).disinh(indEpoch);
                jei_in=-[dataload.data(find(strcmp(names,'PV to assembly'))).disinh(indEpoch),dataload.data(find(strcmp(names,'SST to assembly'))).disinh(indEpoch),0];
                jei_out=-[dataload.data(find(strcmp(names,'PV to bkg'))).disinh(indEpoch),dataload.data(find(strcmp(names,'SST to bkg'))).disinh(indEpoch),0];
            elseif strcmp(Opt,'MarcelPost_inh')
                % use weights at the end of training
                filename=fullfile('marcel_parameters','marcel.mat');
                dataload=load(filename,'data');
                names={dataload.data(:).names};
                indEpoch=find(dataload.data(10).inh==dataload.data(1).endTraining); % find epoch of end training
                jee_in=dataload.data(find(strcmp(names,'Intraassembly'))).inh(indEpoch);
                jee_out=dataload.data(find(strcmp(names,'Interassembly'))).inh(indEpoch);
                jee_bg_in=dataload.data(find(strcmp(names,'Assembly to bkg'))).inh(indEpoch);
                jee_in_bg=dataload.data(find(strcmp(names,'Bkg to assembly'))).inh(indEpoch);
                jee_bg_bg=dataload.data(find(strcmp(names,'Bkg to bkg'))).inh(indEpoch);
                jei_in=-[dataload.data(find(strcmp(names,'PV to assembly'))).inh(indEpoch),dataload.data(find(strcmp(names,'SST to assembly'))).inh(indEpoch),0];
                jei_out=-[dataload.data(find(strcmp(names,'PV to bkg'))).inh(indEpoch),dataload.data(find(strcmp(names,'SST to bkg'))).inh(indEpoch),0];
            elseif strcmp(Opt,'Marcel0')
                % use weights at the end of training
                filename=fullfile('marcel_parameters','marcel.mat');
                dataload=load(filename,'data');
                names={dataload.data(:).names};
                indEpoch=1; % start of training
                jee_in=dataload.data(find(strcmp(names,'Intraassembly'))).inh(indEpoch);
                jee_out=dataload.data(find(strcmp(names,'Interassembly'))).inh(indEpoch);
                jee_bg_in=dataload.data(find(strcmp(names,'Assembly to bkg'))).inh(indEpoch);
                jee_in_bg=dataload.data(find(strcmp(names,'Bkg to assembly'))).inh(indEpoch);
                jee_bg_bg=dataload.data(find(strcmp(names,'Bkg to bkg'))).inh(indEpoch);
                jei_in=-[dataload.data(find(strcmp(names,'PV to assembly'))).inh(indEpoch),dataload.data(find(strcmp(names,'SST to assembly'))).inh(indEpoch),0];
                jei_out=-[dataload.data(find(strcmp(names,'PV to bkg'))).inh(indEpoch),dataload.data(find(strcmp(names,'SST to bkg'))).inh(indEpoch),0];
            end            
            
            % check #clusters and coding level are consistent
            NcUnits=round(f*N_e);    %  number of Exc units per cluster
            fprintf('  --- Synaptic weights: %d units/cluster \n',NcUnits);
            Numbg=round(N_e*(1-f*p)); % number of background (i.e. non-selective) Exc units
            fprintf('  --- fraction of bg Exc units: %0.03g',Numbg/N_e);
%             jee_in=jee_in*ones(1,Q);
            switch Network.clust
                case 'hom'
                    popsize=repmat(NcUnits,Q,1); % size of each Exc cluster
                case 'het'
                    Nc=[];
                    clust_std=Network.clust_std;
                    while (sum(Nc)-(N_e-Numbg))~=0 || any(Nc<0)
                        Nc=round(NcUnits+(NcUnits*clust_std)*randn(Q,1));
                    end
                    popsize=Nc; % array of cluster sizes
                    if any(sum(popsize)-(N_e-Numbg))
                        fprintf('\n---ERROR: Heterogeneous clusters: Problem with cluster numbers\n');
                    end
                    fprintf('\n  --- het cluster sizes->rescale Jplus in each cluster. JEE+:');
%                     jee_in=jee_in*mean(popsize)./popsize';
                    disp(jee_in);
            end
            cusumNcE=[0 cumsum(popsize)'];
            % background units (if present), if not, override in next line
            JEE=(jee*(ones(N_e)+delta*randn(N_e,N_e))).*(rand([N_e,N_e])<p_eeout);
            
            % weights involving inhibitory neurons (including InhTypes)
            JEI=[]; JEIbg=[]; JIE=[]; JII=[];
            for i=1:numel(N_i)
                % to assemblies
                JEI=[JEI, (jei_in(i)*(ones(cusumNcE(Q+1),N_i(i))+deltaEI*randn(cusumNcE(Q+1),N_i(i)))).*(rand([cusumNcE(Q+1),N_i(i)])<p_ei(i))];
                % to bg
                JEIbg=[JEIbg, (jei_out(i)*(ones(N_e-cusumNcE(Q+1),N_i(i))+deltaEI*randn(N_e-cusumNcE(Q+1),N_i(i)))).*(rand([N_e-cusumNcE(Q+1),N_i(i)])<p_ei(i))];
                % all E to I
                JIE=[JIE; (jie(i)*(ones(N_i(i),N_e)+deltaIE*randn(N_i(i),N_e))).*(rand([N_i(i),N_e])<p_ie(i))];
                JIItemp=[];
                for j=1:numel(N_i)
                    JIItemp=[JIItemp,(jii(i,j)*(ones(N_i(i),N_i(j))+delta*randn(N_i(i),N_i(j)))).*(rand([N_i(i),N_i(j)])<p_ii(i,j))];
                end
                JII=[JII;JIItemp];
            end
            JEI=[JEI;JEIbg];
            
            if strcmp(Network.clust,'het') || strcmp(Network.clust,'hom')
                % clustered units: inter-cluster weights
                JEE(1:cusumNcE(Q+1),1:cusumNcE(Q+1))=...
                    (jee_out*(ones(cusumNcE(Q+1))+delta*randn(cusumNcE(Q+1),...
                    cusumNcE(Q+1)))).*(rand([cusumNcE(Q+1),cusumNcE(Q+1)])<p_eeout); % inter-cluster weights
                for clu=2:Q+1 % intra-cluster higher weights
                    JEE(1+cusumNcE(clu-1):cusumNcE(clu),1+cusumNcE(clu-1):cusumNcE(clu))=...
                        (1.*jee_in*(ones(popsize(clu-1))+delta*randn(popsize(clu-1),popsize(clu-1)))).*...
                        (rand([popsize(clu-1),popsize(clu-1)])<p_eein);
                end
                % bg to bg
                JEE(cusumNcE(Q+1)+1:N_e,cusumNcE(Q+1)+1:N_e)=...
                    (jee_bg_bg*(ones(N_e-cusumNcE(Q+1))+delta*randn(N_e-cusumNcE(Q+1),...
                    N_e-cusumNcE(Q+1)))).*(rand([N_e-cusumNcE(Q+1),N_e-cusumNcE(Q+1)])<p_eeout); % inter-cluster weights
                % bg to clusters
                JEE(1:cusumNcE(Q+1),cusumNcE(Q+1)+1:N_e)=...
                    (jee_in_bg*(ones(cusumNcE(Q+1),N_e-cusumNcE(Q+1))+delta*randn(cusumNcE(Q+1),...
                    N_e-cusumNcE(Q+1)))).*(rand([cusumNcE(Q+1),N_e-cusumNcE(Q+1)])<p_eeout); % inter-cluster weights
                % clusters to bg
                JEE(cusumNcE(Q+1)+1:N_e,1:cusumNcE(Q+1))=...
                    (jee_bg_in*(ones(N_e-cusumNcE(Q+1),cusumNcE(Q+1))+delta*randn(N_e-cusumNcE(Q+1),...
                    cusumNcE(Q+1)))).*(rand([N_e-cusumNcE(Q+1),cusumNcE(Q+1)])<p_eeout); % inter-cluster weights
            end
            clustermatrix=eye(Q);
            

           
            
            JEI(JEI>0)=0;
            JIE(JIE<0)=0;
            JII(JII>0)=0;
            JEE(JEE<0)=0;
            J=[JEE JEI; JIE JII];
            J=J-diag(diag(J)); % eliminate self-couplings
            fprintf('  --- New synaptic weights set: done...\n');
            InhNames={'PV','SST','VIP'};
            if numel(jie)>1
                fprintf('      Overall: Jee=%g\n',jee);
                for i=1:numel(N_i)
                    fprintf('      J[%s,E]=%g -- J[Eclust,%s]=%g \n',InhNames{i},jie(i),InhNames{i},jei_in(i));
                end
                for i=1:numel(N_i); for j=1:numel(N_i)
                        fprintf('      J[%s,%s]=%g  ',InhNames{i},InhNames{j},jii(i,j));
                    end; fprintf('\n');
                end
            else
                fprintf('      Overall: Jee=%g -- Jie=%g -- Jei=%g -- Jii=%g \n',jee,jie,jei,jii);
            end
            fprintf('      Var[J]=(Jx%0.03g)^2\n',delta);
            fprintf('      J+=%0.03g\n',mean(jee_in));
            fprintf('      J-(clust,clust)=%0.03g\n',mean(jee_out));
            
            params.popsize=popsize;
            params.clustermatrix=clustermatrix;
            
            if strcmp(Network.clust,'het')
                for i=1:Q
                    a(i)=sum(popsize(clustermatrix(:,i)>0));
                end
                fprintf('  --- clusters size: mean=%0.03g neurons/cluster, sd/mean=%0.03g\n',mean(a),std(a)/mean(a));
            end
            fprintf('  --- fraction of bg Exc units: %0.03g\n',Numbg/N_e);
            
%             %--------------
%             % POISSON INPUT
%             %--------------
%             % generate weights for external neurons
%             if strcmp(Stimulus.input,'Poisson')
%                 JExt=[];
%                 JExt=[JExt; Jee_ext*(rand([N_e,N_ext])<p_ext)];
%                 for i=1:numel(Jie_ext)
%                     JExt=[JExt; Jie_ext(i)*(rand([N_i(i),N_ext])<p_ext)];
%                 end
%                 save(fullfile(params.savedir,'Ext_weights.mat'),'JExt');
%             end
            
            
            %-------------------
            % PLOT weight matrix
            %-------------------
            figure(1); clf;
            subplot(2,1,1);
            colormap(utils.redblue); %xy=J; fun_colormapLim;
            imagesc(J);
            utils.figset(gca,'neurons','neurons','weights',10);
            colorbar;
            subplot(2,1,2);
            lambda=eig(J);
            plot(real(lambda),imag(lambda),'.');
            utils.figset(gca,'Re(\lambda)','Im(\lambda)','eig(weights)',10);
            saveas(gcf,fullfile(params.savedir,'weights.pdf'),'pdf');
            
            %             % spike thresholds for each population
            if any(strcmp(fieldnames(Network),'ThresholdOpt'))
                if Network.ThresholdOpt
                    flg.return_value=1; flg.screen=1;
                    Ni=[params.ni_e,params.ni_i];
                    [theta,~]=auxMFT.fun_Fix_Threshold_InhType(Ni,params,flg);
                    params.theta_e=theta(1); params.theta_i=theta(2:end);
                    fprintf('Spike thresholds Theta calculated for each population and stored in params\n');
                end
            end
%             params.theta_e=3.9; params.theta_i=[4,4,4];
            Theta=[params.theta_e, params.theta_i];
            fprintf('Spike thresholds: theta_e=%0.03g,theta_i=',Theta(1)');
            for i=1:numel(theta_i); fprintf('%0.03g,',Theta(i+1)); end
            params.Theta=Theta;
        end
        
        
                
        %% fun_stim generate stimulus profile and selectivity indices
        
        function [stimulus_save, params]=fun_stim(params)
            
            % unpack vars
            p=params.p;
            Stimulus=params.Stimulus;
            stimuli=params.stimuli;
            Sim=params.Sim;
            
            % for each event, create external current and  stim properties in .Ext
            % and add clusters selective to the event .Stimulus.feat(n).StimClust
            stimulus_save=struct('Ext',[],'Stimulus',[]);
            
            % select stimuli
            temp_Stimulus=struct('input',Stimulus.input);
            indfeat=zeros(1,numel(stimuli));
            
            for ev=1:numel(stimuli)
                fprintf('Stimulus %s',stimuli{ev});
                % match current stimuli to features in Stimulus
                indfeat(ev)=find(cell2mat(arrayfun(@(x)strcmp(stimuli{ev},x(:).name),...
                    Stimulus.feat(:),'uniformoutput',false)));
            end
            fprintf('\n');
            if ~isempty(indfeat)
                temp_Stimulus.feat(1:numel(indfeat))=Stimulus.feat(indfeat);
                
                % define stimulus selectivity: which clusters each stimulus
                % targets
                for n=1:numel(indfeat)
                    sclust=[];
                    
                    switch Stimulus.feat(n).selectivity
                        case {'mixed','allsel'}
                            % this is the option for the US stimulus
                            pstim=0.5; % probability that a cluster is selective to a stimulus
                            selective=rand(1,p)<pstim;
                        case 'exc'
                            % this is the option for the perturbation: CSgauss stimulus (all neurons
                            % are selective)
                            selective=ones(1,p);
                    end
                    temp_Stimulus.feat(n).selective=selective;
                    
                    
                    if ~isempty(temp_Stimulus.feat(n).selective)
                        sclust=find(temp_Stimulus.feat(n).selective(1,:));
                    end
                    temp_Stimulus.feat(n).StimClust=sclust;
                end
            end
            Stimulus=temp_Stimulus;
            Ext=struct('Mu',[]);
            
            % LOAD PARAMETERS
            fieldNames={'Sim','Network','p','popsize','clustermatrix','N_e','N_i','tau_e','tau_i','fieldNames'};
            utils.v2struct(params,fieldNames);
            cusumNcE=[0 cumsum(popsize)'];
            Tseq=Sim.t_Start:Sim.dt_step:Sim.t_End;
            
            if ~isempty(stimuli)
                feat=Stimulus.feat;
                nstim=numel(feat); % number of stimuli in current trials
                stim=repmat(struct('profile',[],'ind',[],'interval',[]),1,nstim);
                temp_ind=repmat(struct('ind',[]),1,nstim); % stores indices for mixed cue (see below)
                for n=1:nstim
                    % stimulus interval
                    interv=feat(n).interval;
                    Gain=feat(n).gain;
                    if ~isempty(strfind(feat(n).name,'gauss'))
                        Gain=1; % with gaussian stim set profile to peak at 1, then multiply each profile by gaussian with SD feat(n).gain for each neuron in feat(n).gauss
                    end
                    Profile=feat(n).profile;
                    Profile=@(t)Profile(t-interv(1));
                    MaxThInput=max(abs(Profile(Tseq(Tseq>interv(1) & Tseq<interv(2)))));
                    Profile=@(t)Gain*Profile(t)/MaxThInput;
                    stim(n).profile=@(t)Profile(t); % fraction increase above baseline
                    % selective neurons
                    StimClust=Stimulus.feat(n).StimClust; % clusters activated by current stimulus
                    % units selective to stimulus
                    ind=[]; % indices of stim sel units
                    switch feat(n).selectivity
                        case 'mixed'
                            for c=StimClust
                                pop_ind=find(clustermatrix(:,c));
                                for k=1:numel(pop_ind)
                                    ind=[ind cusumNcE(pop_ind(k))+1:cusumNcE(pop_ind(k)+1)]; % stim selective units
                                end
                            end
                        case 'exc'
                            ind=1:N_e;
                        otherwise
                            ind=1:N_e;
                    end
                    % sparsify
                    a=randperm(numel(ind));
                    temp_ind(n).ind=ind;
                    ind=ind(a(1:round(feat(n).connectivity*numel(ind))));
                    % gaussian stimulus, draw from randn
                    if ~isempty(strfind(feat(n).name,'gauss'))
                        stim(n).gauss=feat(n).gain*randn(numel(ind),1);
                    end
                    %
                    stim(n).ind=ind;
                    stim(n).interval=interv;
                    stim(n).name=feat(n).name;
                    stim(n).StimClust=StimClust;
                    stim(n).selectivity=feat(n).selectivity;
                    
                end
                Ext.stim=stim;
            end
            Ext.Mu=params.Mu;
            stimulus_save.Ext=Ext;
            stimulus_save.Stimulus=temp_Stimulus;
            
        end
        
        
        %%
        % SIM of one trials given parameters
        %
        % Luca Mazzucato March 2014
        
        % SET OPTIONS
        % ParamsRun = structure containing parameters for simulation
        
        
        %%
        % SIM of one trials given parameters
        %
        % Luca Mazzucato March 2014
        
        % SET OPTIONS
        % ParamsRun = structure containing parameters for simulation
        
        
        function [all_firings, PlotData]=fun_LIF_SIM(ParamsRun)
            
            
            Theta=ParamsRun.Theta; Sim=ParamsRun.Sim; stimuli=ParamsRun.stimuli;%Stimulus=ParamsRun.Stimulus;
            Ext=ParamsRun.Ext; J=ParamsRun.J; N_e=ParamsRun.N_e; N_i=ParamsRun.N_i; p=ParamsRun.p; He=ParamsRun.He;
            Hi=ParamsRun.Hi; tau_e=ParamsRun.tau_e; tau_i=ParamsRun.tau_i; tausyn_e=ParamsRun.tausyn_e;
            tausyn_i=ParamsRun.tausyn_i; tau_arp=ParamsRun.tau_arp;
            %
            all_firings=[];
            dt=Sim.dt_step;            % time step (s)
            Tseq=Sim.t_Start:dt:Sim.t_End;
            
            %--------------------
            % PARAMETERS
            %--------------------
            % CELL
            VEreset=He*Theta(1);
            VIreset=Hi.*Theta(2:end);
            %
            %----------
            % STIMULUS
            %----------
            % add stimuli on top of baseline: for each stimulus provide
            %              - profile (perc. increase on baseline current)
            %              - index of selective neurons
            % BASELINE EXTERNAL CURRENT
            mu=Ext.Mu; % mu is an (N_e+N_i)-th array
            if ~isempty(stimuli)
                stim=Ext.stim;
                nstim=numel(stim); % number of stimuli in current trials
            end
            %----------------
            % SYNAPTIC FILTER
            %----------------
            Tau.tausyn_e=tausyn_e; % exc synaptic time (fall time)
            Tau.tausyn_i=tausyn_i; % exc synaptic time (fall time)
            F=synaptic_trace(Tau,dt,N_e,N_i); % traces for recurrent connections
            %--------------------------
            % SIMULATION
            %--------------------------
            % preallocate memory for stored variable firings_tmp
            % INITIAL CONDITIONS: random
            v=[(Theta(1)-VEreset)/2*ones(N_e,1)+(Theta(1)-VEreset)/2*(2*rand(N_e,1)-1)];   % Initial values of v
            VTh=[Theta(1)*ones(N_e,1)]; % THRESHOLD VECTOR
            c=[VEreset*ones(N_e,1)]; % reset potentials
            tau=[tau_e*ones(N_e,1)]; % membrane time
            for i=1:numel(N_i)
                v=[v; (Theta(1+i)-VIreset(i))/2*ones(N_i(i),1)+(Theta(1+i)-VIreset(i))/2*(2*rand(N_i(i),1)-1)];
                VTh=[VTh; Theta(i+1)*ones(N_i(i),1)];
                c=[c;  VIreset(i)*ones(N_i(i),1)];
                tau=[tau;       tau_i(i)*ones(N_i(i),1)];
            end
            %
            firings=zeros(10*numel(Tseq),2);
            firings_cnt=0;
            tic
            %--------------------
            % PLOT
            %--------------------
            PlotData=[];
            PlotData.Ne_plot=N_e; % number of exc neuron to plot
            PlotData.Ni_plot=(N_i); % number of inh neurons to plot
            cusumpop=1+[0 cumsum([N_e,N_i])];
            ind_plot=cusumpop(1:end-1); % indices of neurons to plot, one per population
            if ~isempty(stimuli)
                indcue=find(cellfun(@(x)~isempty(x),strfind({stim(:).name},'CS')));
                if ~isempty(indcue)
                    ind_plot(1)=stim(indcue).ind(1);
                end
            end
            nplot=numel(ind_plot); % number of neurons to plot (membrane potential plot)
            vi=0; % running index for vplot
            PlotData.vplot = zeros(nplot,round(Sim.plot_length/dt)); % store membrane potential for plots; rows=neurons, cols=time steps;
            PlotData.iEplot = zeros(2,round(Sim.plot_length/dt)); % store EPSC for plots; rows=neurons, cols=time steps;
            PlotData.iExtplot = zeros(2,round(Sim.plot_length/dt)); % store IPSC for plots; rows=neurons, cols=time steps;
            PlotData.iIplot = zeros(2,round(Sim.plot_length/dt),numel(N_i)); % store IPSC for plots; rows=neurons, cols=time steps; 3rd: inh cell types
            PlotData.p=p;
            PlotData.VTh=VTh;
            PlotData.tau=tau;
            PlotData.ind_plot=ind_plot;
            %----------------------------
            % RUN
            %----------------------------
            refr=zeros(size(mu,1),1);       % neurons in refractory state
            for t=1:numel(Tseq)         % siMulation of 1000 ms
                fired=find(v>VTh); % indices of spikes
                Isyn=zeros(N_e+sum(N_i),1);
                % spikes
                if ~isempty(fired)
                    v(fired)=c(fired);
                    refr(fired)=tau_arp;
                end
                % recurrent synaptic current
                F=syn_evolve(F,fired);
                % integrate
                muRun=mu;
                if ~isempty(stimuli)
                    for n=1:nstim
                        if Tseq(t)>=stim(n).interval(1) && Tseq(t)<=stim(n).interval(2)
                            if strcmp(stim(n).name,'CSgauss')
                                muRun(stim(n).ind)=muRun(stim(n).ind)+stim(n).profile(Tseq(t))*mu(stim(n).ind).*stim(n).gauss;
                            else
                                muRun(stim(n).ind)=muRun(stim(n).ind)+stim(n).profile(Tseq(t))*mu(stim(n).ind);
                            end
                        end
                    end
                end
                Isyn=Isyn+J*F.f;
                v=v-v*dt./tau+muRun(:)*dt+Isyn*dt; % units: [V]=[tau*mu]=[tau*Isyn]
                % neurons in refractory state
                refr=max(-0.001,refr-dt);
                v(refr>0)=c(refr>0);
                % store spikes
                if ~isempty(fired)
                    % if firings_tmp has no more space, preallocate more memory
                    if firings_cnt+numel(fired)>size(firings,1)
                        firings=[firings; zeros(10*numel(Tseq),2)];
                    end
                    firings(firings_cnt+1:firings_cnt+numel(fired),1:2)=[Tseq(t)+0*fired, fired];
                    firings_cnt=firings_cnt+numel(fired);
                end
                % store values for plotting, only last Sim.plot_length interval
                if Tseq(t)>Sim.t_End-Sim.plot_length
                    vi=vi+1;
                    % membrane potential
                    PlotData.vplot(1:nplot,vi)=v(ind_plot);
                    % input currents
                    PlotData.iEplot(1:nplot,vi)=J(ind_plot,1:N_e)*F.f(1:N_e);
                    for i_i=1:numel(N_i)
                        PlotData.iIplot(1:nplot,vi,i_i)=J(ind_plot,cusumpop(1+i_i):cusumpop(2+i_i)-1)*F.f(cusumpop(1+i_i):cusumpop(2+i_i)-1);
                    end
                    PlotData.iExtplot(1:nplot,vi)=muRun(ind_plot,1);
                end
            end
            fprintf('--- End of trial...\n');
            toc
            %---------------------------------------------------------------------------
            if ~any(any(firings))
                fprintf('\n --- NO SPIKES GENERATED... \n');
            else
                % find last spike in firings
                IndexEnd=find(firings(:,2)==0,1)-1;
                if isempty(IndexEnd)
                    IndexEnd=size(firings,1);
                end
                all_firings=firings(1:IndexEnd,[1 2]);
            end
            
            
            function F=synaptic_trace(Tau,dt,N_e,N_i)
                
                F=struct();
                tau_sE=Tau.tausyn_e; % exc synaptic time (fall time)
                tau_sI=Tau.tausyn_i; % inh synaptic time (fall time)
                fexp=repmat(exp(-dt/tau_sE),N_e,1); % Multiplicative step (fp)
                fSpike=repmat((1/tau_sE),N_e,1); % add to fp with a spike
                for i_n=1:numel(N_i)
                    fexp=[fexp; repmat(exp(-dt/tau_sI(i_n)),N_i(i_n),1)];
                    fSpike=[fSpike; repmat((1/tau_sI(i_n)),N_i(i_n),1)];
                end
                f=zeros(N_e+sum(N_i),1);
                F.fexp=fexp;
                F.fSpike=fSpike;
                F.f=f;
                
            end
            
            function F=syn_evolve(F,fired)
                
                
                % update synaptic filter
                F.f=F.fexp.*F.f;
                if ~isempty(fired)
                    % update of synaptic filter with spikes
                    F.f(fired)=F.f(fired)+F.fSpike(fired);
                end
            end
            
            
        end        
        
        function [all_firings, PlotData]=fun_LIF_SIM_Poisson(ParamsRun)
            
            
            Theta=ParamsRun.Theta; Sim=ParamsRun.Sim; stimuli=ParamsRun.stimuli;%Stimulus=ParamsRun.Stimulus;
            Ext=ParamsRun.Ext; J=ParamsRun.J; N_e=ParamsRun.N_e; N_i=ParamsRun.N_i; p=ParamsRun.p; He=ParamsRun.He;
            Hi=ParamsRun.Hi; tau_e=ParamsRun.tau_e; tau_i=ParamsRun.tau_i; tausyn_e=ParamsRun.tausyn_e;
            tausyn_i=ParamsRun.tausyn_i; tau_arp=ParamsRun.tau_arp; 
            Stimulus=ParamsRun.Stimulus; N_ext=ParamsRun.N_ext;
            ni_ext_i=ParamsRun.ni_ext_i; ni_ext_e=ParamsRun.ni_ext_e;
            Jee_ext=ParamsRun.Jee_ext; Jie_ext=ParamsRun.Jie_ext;
            % 
            all_firings=[];
            dt=Sim.dt_step;            % time step (s)
            Tseq=Sim.t_Start:dt:Sim.t_End;
            
            %--------------------
            % PARAMETERS
            %--------------------
            % CELL
            VEreset=He*Theta(1);
            VIreset=Hi.*Theta(2:end);
            %
            %----------
            % STIMULUS
            %----------
            % add stimuli on top of baseline: for each stimulus provide
            %              - profile (perc. increase on baseline current)
            %              - index of selective neurons
            % BASELINE EXTERNAL CURRENT
            mu=Ext.Mu; % mu is an (N_e+N_i)-th array
            if ~isempty(stimuli)
                stim=Ext.stim;
                nstim=numel(stim); % number of stimuli in current trials
            end
            if strcmp(Stimulus.input,'Poisson')
%                 if p_ext==1 
%                     JExt=[Jee_ext*(rand([N_e,Next])<pext); Jie_ext*(rand([N_i,Next])<pext)]; % external synaptic matrix
                    ext_freq_e=N_ext*ni_ext_e; % mean external poisson rate
                    ext_freq_i=N_ext*ni_ext_i; % mean external poisson rate
                    JExt=Jee_ext*(ones(N_e,1)); for i=1:numel(N_i); JExt=[JExt; Jie_ext(i)*(ones(N_i(i),1))]; end % external synaptic matrix
%                 end
                lambda=ext_freq_e*ones(N_e,1);%+0*(ext_freq_e/4)*randn(N_ext,1); % mean firing rate (Hz)
                for i=1:numel(N_i); lambda=[lambda; ext_freq_i(i)*ones(N_i(i),1)]; end%+0*(ext_freq_e/4)*randn(N_ext,1); % mean firing rate (Hz)
                ns=round(ext_freq_e*(Sim.t_End-Sim.t_Start)*5);             % number of spikes to be generated (10x more than needed)                
                isi1=-repmat((1./lambda),1,ns).*log(rand(size(JExt,1),ns));   %generation of exponential distribution of ISI
                spk_times=Sim.t_Start+cumsum(isi1,2); % input spike times
                avecheck=sum((spk_times>=Sim.t_Start & spk_times<Sim.t_End),2)/(Sim.t_End-Sim.t_Start);
                fprintf('ave ext Poisson input neuron 1=%0.03g Hz, neuron 2=%0.03g Hz',avecheck(1),avecheck(2));
            end
            
            %----------------
            % SYNAPTIC FILTER
            %----------------
            Tau.tausyn_e=tausyn_e; % exc synaptic time (fall time)
            Tau.tausyn_i=tausyn_i; % exc synaptic time (fall time)
            F=synaptic_trace(Tau,dt,N_e,N_i); % traces for recurrent connections
            if strcmp(Stimulus.input,'Poisson')
                FExt=synaptic_trace(Tau,dt,N_e,N_i); % traces for external input
            end
            %--------------------------
            % SIMULATION
            %--------------------------
            % preallocate memory for stored variable firings_tmp
            % INITIAL CONDITIONS: random
            v=[(Theta(1)-VEreset)/2*ones(N_e,1)+(Theta(1)-VEreset)/2*(2*rand(N_e,1)-1)];   % Initial values of v
            VTh=[Theta(1)*ones(N_e,1)]; % THRESHOLD VECTOR
            c=[VEreset*ones(N_e,1)]; % reset potentials
            tau=[tau_e*ones(N_e,1)]; % membrane time
            for i=1:numel(N_i)
                v=[v; (Theta(1+i)-VIreset(i))/2*ones(N_i(i),1)+(Theta(1+i)-VIreset(i))/2*(2*rand(N_i(i),1)-1)];
                VTh=[VTh; Theta(i+1)*ones(N_i(i),1)];
                c=[c;  VIreset(i)*ones(N_i(i),1)];
                tau=[tau;       tau_i(i)*ones(N_i(i),1)];
            end
            %
            firings=zeros(10*numel(Tseq),2);
            firings_cnt=0;
            tic
            %--------------------
            % PLOT
            %--------------------
            PlotData=[];
            PlotData.Ne_plot=N_e; % number of exc neuron to plot
            PlotData.Ni_plot=(N_i); % number of inh neurons to plot
            cusumpop=1+[0 cumsum([N_e,N_i])];
            ind_plot=cusumpop(1:end-1); % indices of neurons to plot, one per population
            if ~isempty(stimuli)
                indcue=find(cellfun(@(x)~isempty(x),strfind({stim(:).name},'CS')));
                if ~isempty(indcue)
                    ind_plot(1)=stim(indcue).ind(1);
                end
            end
            nplot=numel(ind_plot); % number of neurons to plot (membrane potential plot)
            vi=0; % running index for vplot
            PlotData.vplot = zeros(nplot,round(Sim.plot_length/dt)); % store membrane potential for plots; rows=neurons, cols=time steps;
            PlotData.iEplot = zeros(2,round(Sim.plot_length/dt)); % store EPSC for plots; rows=neurons, cols=time steps;
            PlotData.iExtplot = zeros(2,round(Sim.plot_length/dt)); % store IPSC for plots; rows=neurons, cols=time steps;
            PlotData.iIplot = zeros(2,round(Sim.plot_length/dt),numel(N_i)); % store IPSC for plots; rows=neurons, cols=time steps; 3rd: inh cell types
            PlotData.p=p;
            PlotData.VTh=VTh;
            PlotData.tau=tau;
            PlotData.ind_plot=ind_plot;
            %----------------------------
            % RUN
            %----------------------------
            refr=zeros(size(mu,1),1);       % neurons in refractory state
            for t=1:numel(Tseq)         % siMulation of 1000 ms
                fired=find(v>VTh); % indices of spikes
                Isyn=zeros(N_e+sum(N_i),1);
                % spikes
                if ~isempty(fired)
                    v(fired)=c(fired);
                    refr(fired)=tau_arp;
                end
                % recurrent synaptic current
                muRun=mu;
                F=syn_evolve(F,fired);
                if strcmp(Stimulus.input,'Poisson')
                    ext_fired=(sum(spk_times>=Tseq(t) & spk_times<Tseq(t)+dt,2)); % spikes from all external units
                    FExt=syn_evolve_ext(FExt,ext_fired);
                    Isyn=Isyn+JExt.*FExt.f;
                    muRun=muRun*0; % constant external bias
                end
                % integrate
                if ~isempty(stimuli)
                    for n=1:nstim
                        if Tseq(t)>=stim(n).interval(1) && Tseq(t)<=stim(n).interval(2)
                            if strcmp(stim(n).name,'CSgauss')
                                muRun(stim(n).ind)=muRun(stim(n).ind)+stim(n).profile(Tseq(t))*mu(stim(n).ind).*stim(n).gauss;
                            else
                                muRun(stim(n).ind)=muRun(stim(n).ind)+stim(n).profile(Tseq(t))*mu(stim(n).ind);
                            end
                        end
                    end
                end
                Isyn=Isyn+J*F.f;
                v=v-v*dt./tau+muRun(:)*dt+Isyn*dt; % units: [V]=[tau*mu]=[tau*Isyn]
                % neurons in refractory state
                refr=max(-0.001,refr-dt);
                v(refr>0)=c(refr>0);
                % store spikes
                if ~isempty(fired)
                    % if firings_tmp has no more space, preallocate more memory
                    if firings_cnt+numel(fired)>size(firings,1)
                        firings=[firings; zeros(10*numel(Tseq),2)];
                    end
                    firings(firings_cnt+1:firings_cnt+numel(fired),1:2)=[Tseq(t)+0*fired, fired];
                    firings_cnt=firings_cnt+numel(fired);
                end
                % store values for plotting, only last Sim.plot_length interval
                if Tseq(t)>Sim.t_End-Sim.plot_length
                    vi=vi+1;
                    % membrane potential
                    PlotData.vplot(1:nplot,vi)=v(ind_plot);
                    % input currents
                    PlotData.iEplot(1:nplot,vi)=J(ind_plot,1:N_e)*F.f(1:N_e);
                    for i_i=1:numel(N_i)
                        PlotData.iIplot(1:nplot,vi,i_i)=J(ind_plot,cusumpop(1+i_i):cusumpop(2+i_i)-1)*F.f(cusumpop(1+i_i):cusumpop(2+i_i)-1);
                    end
                    PlotData.iExtplot(1:nplot,vi)=muRun(ind_plot,1);
                end
            end
            fprintf('--- End of trial...\n');
            toc
            %---------------------------------------------------------------------------
            if ~any(any(firings))
                fprintf('\n --- NO SPIKES GENERATED... \n');
            else
                % find last spike in firings
                IndexEnd=find(firings(:,2)==0,1)-1;
                if isempty(IndexEnd)
                    IndexEnd=size(firings,1);
                end
                all_firings=firings(1:IndexEnd,[1 2]);
            end
            
            
            function F=synaptic_trace(Tau,dt,N_e,N_i)
                
                F=struct();
                tau_sE=Tau.tausyn_e; % exc synaptic time (fall time)
                tau_sI=Tau.tausyn_i; % inh synaptic time (fall time)
                fexp=repmat(exp(-dt/tau_sE),N_e,1); % Multiplicative step (fp)
                fSpike=repmat((1/tau_sE),N_e,1); % add to fp with a spike
                for i_n=1:numel(N_i)
                    fexp=[fexp; repmat(exp(-dt/tau_sI(i_n)),N_i(i_n),1)];
                    fSpike=[fSpike; repmat((1/tau_sI(i_n)),N_i(i_n),1)];
                end
                f=zeros(N_e+sum(N_i),1);
                F.fexp=fexp;
                F.fSpike=fSpike;
                F.f=f;
                
            end
            
            function F=syn_evolve(F,fired)
                
                
                % update synaptic filter
                F.f=F.fexp.*F.f;
                if ~isempty(fired)
                    % update of synaptic filter with spikes
                    F.f(fired)=F.f(fired)+F.fSpike(fired);
                end
            end
            
            function F=syn_evolve_ext(F,fired)
                
                
                % update synaptic filter
                F.f=F.fexp.*F.f;
                if any(fired)
                    % update of synaptic filter with spikes
                    F.f(fired>0)=F.f(fired>0)+fired(fired>0).*F.fSpike(fired>0);
                end
            end            
        end                
        
        
        
    end
end