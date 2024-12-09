function MFT_fun_rate_changes(tmpv1)

options_main=tmpv1;
Jcases={'Low','Hi'};
Jvalue.Hi=options_main.Jplus;
Jvalue.Low=options_main.JplusLow;


utils.v2struct(options_main);
ncueset=1;
if cue_mode
    if iscell(cueset)
        ncueset=numel(cueset{1});
    else
        ncueset=numel(cueset);
        cueset={cueset};
    end
elseif any(strcmp(fieldnames(options_main),'InhDensity'))
    ncueset=numel(InhDensity);
    cueset=repmat({1:numel(InhDensity)},1,ncueset);
end

% global cue_mode all_mode
% all_mode=1;
% if strcmp(Network.cuemode,'none')
%     cue_mode=0;
% else
%     cue_mode=1;
% end
aa1=struct('rate',[],'Clust',[],'cue',[],'stabmat',[],'phiprime',[],'linresp',[],'A',[]);
data=struct('Hi',aa1,'Low',aa1);
% plot separately
for i_c=1:ncueset
    %     cue=cueset(i_c);
    if cue_mode
        for i=1:numel(cue.cue_option); cue.cue_value{i}=cueset{i}(i_c); options_main.cue.cue_value{i}=cueset{i}(i_c); end
    elseif any(strcmp(fieldnames(options_main),'InhDensity'))
        options_main.InhDensity=tmpv1.InhDensity{i_c};
    end
    options_main.extra=tmpv1.extra{i_c};
    
    aux.MFT_preamble;
    for iJ=1:numel(Jcases)
        Jplus=Jvalue.(Jcases{iJ});
        % high firing rate
        save(paramsfile,'Jplus','-append');
        params=load(paramsfile);
        %         flg.stim_mode=1;
        flg.all_mode=0;
        %         flg.nivsj_mode=1;
        % preamble
    %     StimHeight=stimset(i_J);
    %     save(paramsfile,'StimHeight','-append');
        % load stable states
        StableStates=[]; % dim: [#states, pops]
    %     p=params.p;
        cnt=0;

        % all other fixed points
        flg.all_mode=1;
        Files=auxMFT.File_Info_MFT(p+4,p,flg);
        fileload=Files.stable;
        AllData=load(fileload);
        % find states
        Ind=find(cell2mat(arrayfun(@(x)~isempty(x(:).Jplus),AllData.StableStates,'uniformoutput',false))); % states with Ind active clusters
        % Select attractors for current Jplus
        for n=Ind
            indJ=find(abs(AllData.StableStates(n).Jplus-Jplus)<1e-3);
            %     [~,indJ]=min(abs(AllData.StableStates(Ind(n)).Jplus-Jplus));
            if ~isempty(indJ)
                cnt=cnt+1;
                StableStates=AllData.StableStates(n).rates(indJ,:);
                data.(Jcases{iJ}).rate=[data.(Jcases{iJ}).rate; StableStates];
                data.(Jcases{iJ}).Clust=[data.(Jcases{iJ}).Clust; n];
                data.(Jcases{iJ}).cue=[data.(Jcases{iJ}).cue; cueset{1}(i_c)];
                % stability matrix
                matrix=squeeze(AllData.StableStates(n).Eig.matrix(:,:,indJ));
                % compute Exc subset
                excm=matrix(1:end-3,1:end-3);
                [V,D]=eig(excm);
                data.(Jcases{iJ}).stabmat=cat(3,data.(Jcases{iJ}).stabmat, matrix);
                % phiprime
                phiprime=squeeze(AllData.StableStates(n).Eig.phiprime(:,indJ));
                data.(Jcases{iJ}).phiprime=cat(2,data.(Jcases{iJ}).phiprime, phiprime);
                % linear response matrix
                linresp=squeeze(AllData.StableStates(n).Eig.linresp(:,:,indJ));
                data.(Jcases{iJ}).linresp=cat(3,data.(Jcases{iJ}).linresp, linresp);
                factortau=repmat([ones(p+1,1)*params.tausyn_e;ones(3,1)*params.tausyn_i]/params.tau_e,1,p+4);
                factorphi=repmat(1./phiprime,1,p+4);
                Aphip=factortau.*matrix+eye(p+4);  
                A=factorphi.*Aphip;
%                 linresp1=inv(diag(1./phiprime)-A);
                data.(Jcases{iJ}).A=cat(3,data.(Jcases{iJ}).A, A);
            end
        end


        flg.all_mode=0;
    end
end
    
%% PLOT Firing rate changes
nn=p+4;
Legends=cell(nn,1);
Legends{1}='sel. pop';
for i=2:nn-4; Legends{i}=sprintf('pop %d',i); end
Legends(nn-3:nn+1)={'bg (Exc)','PV','SST','VIP','mean Inh'};

inh_frac=params.N_i/sum(params.N_i); % fraction of inh types
for iJ=1:numel(Jcases)
    figure(iJ); clf;
    rate=data.(Jcases{iJ}).rate; Clust=data.(Jcases{iJ}).Clust; cueval=data.(Jcases{iJ}).cue;
    cnt=0; YLIM=[min(rate(:)),max(rate(:))];
    for ic=unique(Clust)'
        if sum(Clust==ic)>1
            cnt=cnt+1;
            subplot(1,numel(unique(Clust)),cnt);
            y=rate(Clust==ic,:);
            x=cueval(Clust==ic)';
            h=plot(x,y,'linewidth',2); hold on;
            % plot sum of all Inh pops
            inhtot=inh_frac*rate(Clust==ic,end-2:end)';
            h=[h;plot(x,inhtot,':','linewidth',2)]; hold on;
            tit=sprintf('J+=%0.03g,ActiveCl=%d',Jvalue.(Jcases{iJ}),ic);
            utils.figset(gca,'pert [%]','Rate [spks/s]',tit,20);
            legend(h,Legends,'location','northeast','fontsize',8);
            ylim(YLIM);
        end
    end
    filesave=Files.paperenergy;
    saveas(gcf,[filesave '[Jp' sprintf('%0.03g',Jvalue.(Jcases{iJ})) ']_RateChange.pdf'],'pdf');
end

%% PLOT Firing rate changes
nn=p+4;
Legends=cell(nn,1);
Legends{1}='sel. pop';
for i=2:nn-4; Legends{i}=sprintf('pop %d',i); end
Legends(nn-3:nn+2)={'bg (Exc)','PV','SST','VIP','mean_I','mean_E'};

inh_frac=params.N_i/sum(params.N_i); % fraction of inh types
exc_frac=[ones(1,p)*params.f*params.N_e,(1-params.p*params.f)*params.N_e]/params.N_e;
for iJ=1:numel(Jcases)
    figure(iJ); clf;
    rate=data.(Jcases{iJ}).rate; Clust=data.(Jcases{iJ}).Clust; cueval=data.(Jcases{iJ}).cue;
    cnt=0; YLIM=[min(rate(:)),max(rate(:))];
    for ic=unique(Clust)'
        if sum(Clust==ic)>1
            cnt=cnt+1;
            subplot(1,numel(unique(Clust)),cnt);
            y=rate(Clust==ic,:);
            % plot sum of all Inh pops
            inhtot=inh_frac*rate(Clust==ic,end-2:end)';
            exctot=exc_frac*rate(Clust==ic,1:p+1)';
            x=cueval(Clust==ic)';
            zerocue=(x==0); if ~any(zerocue); zerocue(1)=1; end
            yref=repmat(y(zerocue,:),numel(zerocue),1);
            inhtot=inhtot-repmat(inhtot(zerocue),1,numel(zerocue));
            exctot=exctot-repmat(exctot(zerocue),1,numel(zerocue));
            y=y-yref;
            h=plot(x,y,'linewidth',2); hold on;
            h=[h;plot(x,inhtot,':','linewidth',2)]; hold on;
            h=[h;plot(x,exctot,'--','linewidth',2)]; hold on;
            tit=sprintf('J+=%0.03g,ActiveCl=%d',Jvalue.(Jcases{iJ}),ic);
            utils.figset(gca,'pert [%]','Rate [spks/s]',tit,20);
            legend(h,Legends,'location','northwest','fontsize',8);
%             ylim(YLIM);
        end
    end
    filesave=Files.paperenergy;
    saveas(gcf,[filesave '[Jp' sprintf('%0.03g',Jvalue.(Jcases{iJ})) ']_DeltaRateChange.pdf'],'pdf');
end


%% PLOT changes in input current

Legends=cell(nn,1);
Legends{1}='sel. pop';
for i=2:nn-4; Legends{i}=sprintf('pop %d',i); end
Legends(nn-3:nn)={'bg (Exc)','PV','SST','VIP'};
inh_frac=params.N_i/sum(params.N_i); % fraction of inh types
for iJ=2%1:numel(Jcases)
    figure(iJ); clf;
    rate=data.(Jcases{iJ}).rate; Atot=data.(Jcases{iJ}).A; Clust=data.(Jcases{iJ}).Clust; cueval=data.(Jcases{iJ}).cue;
    cnt=0; YLIM=[min(rate(:)),max(rate(:))];
    for ic=unique(Clust)'
        if sum(Clust==ic)>1
            cnt=cnt+1;
            subplot(1,numel(unique(Clust)),cnt);
            A=(Atot(:,:,Clust==ic));
            y=rate(Clust==ic,:);
            x=cueval(Clust==ic)';
            inhCurr=NaN(numel(x),nn);
            excCurr=NaN(numel(x),nn);
            for icue=1:numel(x)
                inhCurr(icue,1:nn)=squeeze(A(:,end-2:end,icue))*y(icue,end-2:end)';
                excCurr(icue,1:nn)=squeeze(A(:,1:end-3,icue))*y(icue,1:end-3)';
            end
            if any(x==0)
                yrefinh=repmat(inhCurr(x==0,:),numel(x),1);
                inhCurr=inhCurr-repmat(inhCurr(x==0,:),numel(x),1);
                excCurr=excCurr-repmat(excCurr(x==0,:),numel(x),1);
            else
                yrefinh=repmat(inhCurr(1,:),numel(x),1);
                inhCurr=inhCurr-repmat(inhCurr(1,:),numel(x),1);
                excCurr=excCurr-repmat(excCurr(1,:),numel(x),1);
            end            
            
            h1=plot(x,inhCurr','-','linewidth',2); hold on;
            plot(0,0,'k'); hold on;
            h2=plot(x,excCurr',':','linewidth',2); hold on;
            % plot sum of all Inh pops
%             inhtot=inh_frac*rate(Clust==ic,end-2:end)';
%             h=[h;plot(x,inhtot,':','linewidth',2)]; hold on;
            tit=sprintf('-=I,:=E;J+=%0.03g,ActiveCl=%d',Jvalue.(Jcases{iJ}),ic);
            utils.figset(gca,'pert [%]','Current [mV]',tit,10);
            if cnt==1; legend([h1;h2],[Legends;Legends],'location','northeast','fontsize',8); end
%             ylim(YLIM);
        end
    end
    filesave=Files.paperenergy;
    saveas(gcf,[filesave '[Jp' sprintf('%0.03g',Jvalue.(Jcases{iJ})) ']_DeltaCurrChange.pdf'],'pdf');
end
    
%% Linear response

poplabel={'E(bg)','PV','SST','VIP'};
% Jtot=Params.A;
% n=size(Jtot,1);
for i=1:params.p
    poplabel=[sprintf('E%d',i),poplabel];
end

for iJ=1:numel(Jcases)
    figure(iJ); clf;
    colormap(utils.redblue);
    linresp=data.(Jcases{iJ}).linresp; Clust=data.(Jcases{iJ}).Clust; cueval=data.(Jcases{iJ}).cue;
    totplot=size(linresp,3);
    cnt=0; %YLIM=[min(rate(:)),max(rate(:))];
    for i=1:size(linresp,3)%ic=unique(Clust)'
%         for icue=unique(cueval)'
            cnt=cnt+1;
            subplot(numel(unique(Clust)),numel(unique(cueval)),cnt);
            x=squeeze(linresp(:,:,cnt));
            imagesc(x);
            for i=1:size(x,1)
                for j=1:size(x,2)
                    text(j-0.25,i,sprintf('%0.01g',x(i,j)),'fontsize',8); hold on;
                end
            end
            jmax=max(abs(x(:)));
            set(gca,'CLim',[-jmax,jmax]);
            colorbar;
            set(gca,'Xtick',1:size(x,1), 'XTickLabel',poplabel,'Ytick',1:size(x,2), 'YTickLabel', poplabel);
            tit=sprintf('J+=%0.03g,ActCl=%d,cue=%0.02g',Jvalue.(Jcases{iJ}),Clust(cnt),cueval(cnt));
            title(tit);
%         end
    end
            p=utils.mtit('rows=out,cols=in');
    filesave=Files.paperenergy;
    saveas(gcf,[filesave '[Jp' sprintf('%0.03g',Jvalue.(Jcases{iJ})) ']_LinResp.pdf'],'pdf');
end

    
