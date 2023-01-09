classdef pltMFT
    methods(Static)
        
        
        %----------------------------------------------------------------
        %----------------------------------------------------------------
        
        
        function fun_plot_weights(Params)
            
            load(Params.Files.paramsfile);
            
            Jtot=[Jee -Jei; Jie, -Jii];
            
            figure(1); clf;
            col1=colormap(utils.redblue);
            poplabel={'Exc','PV','SST','VIP'};
            imagesc(Jtot); hold on;
            for i=1:size(Jtot,1)
                for j=1:size(Jtot,2)
                    text(j-0.2,i,sprintf('%0.02g',Jtot(i,j)),'fontsize',20); hold on;
                end
            end
            jmax=max(abs(Jtot(:)));
            set(gca,'CLim',[-jmax,jmax]);
            colorbar;
            set(gca,'Xtick',1:size(Jtot,1), 'XTickLabel',poplabel,'Ytick',1:size(Jtot,2), 'YTickLabel', poplabel);
            title('rows=post,cols=pre');
            filename=fullfile([Params.Files.paperfig 'weights.pdf']);
            saveas(gcf,filename,'pdf');
            
            
            figure(2); clf;
            col1=colormap(utils.redblue);
            poplabel={'E(bg)','PV','SST','VIP'};
            Jtot=Params.A;
            n=size(Jtot,1);
            for i=1:n-4
                poplabel=[sprintf('E%d',i),poplabel];
            end
            imagesc(Jtot); hold on;
            for i=1:size(Jtot,1)
                for j=1:size(Jtot,2)
                    text(j-0.2,i,sprintf('%0.02g',Jtot(i,j)),'fontsize',20); hold on;
                end
            end
            jmax=max(abs(Jtot(:)));
            set(gca,'CLim',[-jmax,jmax]);
            colorbar;
            set(gca,'Xtick',1:size(Jtot,1), 'XTickLabel',poplabel,'Ytick',1:size(Jtot,2), 'YTickLabel', poplabel);
            title('rows=post,cols=pre');
            filename=fullfile([Params.Files.paperfig '[' num2str(n) ']A-matrix.pdf']);
            saveas(gcf,filename,'pdf');
        end
        
        
        %----------------------------------------------------------------
        %----------------------------------------------------------------
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % function PlotStableStates
        %
        % overlays all stable states for ungrouped populations nn=p+2
        
        function PlotAllFixedPoints(p,Jstop,NumFig,flg)
            
            % global flg.nivsj_mode flg.flg.spur_mode
            % preamble
            colors=lines;
            % colors=colors(randperm(size(colors,1)),:);
            %------------------------
            % FILE INFO
            %------------------------
            if flg.nivsj_mode==0 && flg.spur_mode==0
                kind='spont';
            elseif flg.nivsj_mode==1
                kind='sel';
            elseif flg.nivsj_mode==0 && flg.spur_mode==1
                kind='spur';
            end
            
            
            %------------------------
            % PLOT FIXED POINTS
            %------------------------
            %
            % Load data file
            
            paramsfile  =flg.paramsfile;
            load(paramsfile);
            nn=p+4;
            Files=auxMFT.File_Info_MFT(nn,p,flg);
            filesave=Files.stable;
            fprintf('%s\n',filesave);
            load(filesave);
            
            figure(NumFig); clf; h=[]; Legends=[];
            % subplot(211)
            % OVERLAY STABLE STATES
            % Each stable network configuration has its own color
            % Inh dashed
            cnt=0;
            % Only keep configurations stable after second critical point
            % Endpoint=zeros(1,numel(FP_data.spur));
            PopInd=find(cell2mat(arrayfun(@(x)~isempty(x(:).Jplus),FP_data.spur,'uniformoutput',false)));
            % for nn=PopInd
            %     EndPoint(nn)=max(FP_data.spur(nn).Jplus);
            % end
            % NOTATIONS: nn=total number of populations;
            %            n=p+3-nn= number of populations firing at hi rate (legend)
            for nn=PopInd
                %     if (EndPoint(nn)==max(EndPoint))
                %         if max(FP_data.spur(nn).Jplus)==EndPoint
                cnt=cnt+1;
                Jplus=FP_data.spur(nn).Jplus;
                rates=FP_data.spur(nn).rates;
                for i=1:p
                    %                 plot(Jplus, rates(:,i),'diamond','color',colors(nn,:),'markersize',8);
                    plot(Jplus, rates(:,i),'diamond','color',[0 0 0]+(nn-2)/(p+5),'markersize',8);
                    hold on
                end
                plot(Jplus,rates(:,p+1),'sq','color','b'); % background Exc
                hold on
                plot(Jplus,rates(:,p+2),'.','color','r'); % PV
                hold on
                plot(Jplus,rates(:,p+3),':','color','r'); % SST
                hold on
                plot(Jplus,rates(:,p+4),'--','color','r'); % VIP
                hold on
                % legend construction
                h(cnt)=plot(Jplus(1),rates(1,1),'diamond','color',[0 0 0]+(nn-2)/(p+5),'markersize',8); % for legend purpose
                hold on
                %             Legends{cnt}=sprintf('n=%d',13-nn);
                Legends{cnt}=sprintf('n=%d',p+5-nn);
                % overlay saddle points
                if any(strcmp(fieldnames(FP_data.spur(nn)),'unstable'))
                    if ~isempty(FP_data.spur(nn).unstable)
                        Junstable=FP_data.spur(nn).unstable.Jplus;
                        rates_unstable=FP_data.spur(nn).unstable.rates;
                        for i=1:p
                            %                     plot(Junstable, rates_unstable(:,i),'x','color',colors(nn,:));
                            plot(Junstable, rates_unstable(:,i),'x','color',[0 0 0]+(nn-2)/(p+5));
                            hold on
                        end
                    end
                end
                %         end
                %   end
            end
            % legend
            h(cnt+1)=plot(Jplus(1),rates(1,p+1),'sq','color',colors(nn,:)); % background Exc
            Legends{cnt+1}='bg (Exc)';
            hold on
            h(cnt+2)=plot(Jplus(1),rates(1,p+2),'.','color',colors(nn,:)); % Inh
            Legends{cnt+2}='PV';
            hold on
            h(cnt+3)=plot(Jplus(1),rates(1,p+3),':','color',colors(nn,:)); % Inh
            Legends{cnt+3}='SST';
            hold on
            h(cnt+4)=plot(Jplus(1),rates(1,p+4),'--','color',colors(nn,:)); % Inh
            Legends{cnt+4}='VIP';
            hold on
            h(cnt+5)=plot(0, 0,'x','color','k');
            Legends{cnt+5}='saddle';
            hold off
            % if strcmp(kind,'spur') && ~isempty('J_below')
            %     h(nn+1)=line([J_below(end) J_below(end)],[min(min(rates)) max(max(rates))],'linestyle','--','linewidth',1.2,'color','k');
            %     Legends{nn+1}='J"';
            % end
            hh=h(1);
            Leg=Legends(1);
            Leg{2}='...';
            hh(2)=0;
            % Leg=[Leg Legends(p:p+3)];
            % hh=[hh h(p:p+3)];
            % legend(hh,Leg,'location','northeastoutside');
            title(sprintf('n=#populations at hi rate, p=%d',p));
            xlabel('J+ (mV/s)');
            ylabel('Firing rate (spks/s)');
            fntsz=15;
            xlim([Jzero Jstop]);
            % utils.figset(gca,xlab,ylab,'',fntsz);
            FPpic=Params.Files.FPstablepic;
            saveas(gcf,FPpic,'pdf')
            %
            % subplot(212)
            % STATES BASIN OF ATTRACTION
            % plot real part of eigenvalues
            
            
            
            return
            nn=size(rates,2);
            %
            h(1)=plot(Jplus_store, rates(:,1),'color',colors(1,:),'linestyle',':','linewidth',1.1);
            hold on
            for i=2:nn
                h(i)=plot(Jplus_store, rates(:,i),'color',colors(i,:));
                hold on
            end
            % legend
            Legends=cell(nn,1);
            Legends{1}='sel. pop';
            for i=2:nn-2
                Legends{i}=sprintf('pop %d',i);
            end
            Legends{nn-1}='bg (Exc)';
            Legends{nn}='Inh';
            % if strcmp(kind,'spur') && ~isempty('J_below')
            %     h(nn+1)=line([J_below(end) J_below(end)],[min(min(rates)) max(max(rates))],'linestyle','--','linewidth',1.2,'color','k');
            %     Legends{nn+1}='J"';
            % end
            legend(h,Legends,'location','northeastoutside');
            title(sprintf('%s. activity: Fixed points (nn=%d, p=%d)',kind,nn,p));
            xlabel('J+ (mV/s)');
            ylabel('Firing rate (spks/s)');
            % xlim([Jplus_store(1) Jplus_store(end)]);
            FPpic=Params.Files.FPpic;
            saveas(gcf,FPpic,'pdf')
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % PLOT EIGENVALUES
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            file_eig=Params.Files.file_eig;
            load(file_eig);
            % Real
            figure(2); clf; h=[];
            subplot(211)
            h(1)=plot(Jplus_store, eigenvalues.Real(1,:),'color',colors(1,:),'linestyle',':','linewidth',1.1);
            hold on
            for i=2:nn
                h(i)=plot(Jplus_store, eigenvalues.Real(i,:),'color',colors(i,:));
                hold on
            end
            hold on
            line([Jplus_store(1) Jplus_store(end)],[0 0],'linestyle','--','linewidth',1.2,'color','k');
            % if strcmp(kind,'spur') && ~isempty('J_below')
            %     h(nn+1)=line([J_below(end) J_below(end)],[min(min(rates)) max(max(rates))],'linestyle','--','linewidth',1.2,'color','k');
            %     Legends{nn+1}='J"';
            % end
            hold off
            legend(h,Legends,'location','northeastoutside');
            xlabel('J+ (mV/s)');
            ylabel('Re(eigenvalues)');
            title(sprintf('%s. activity: Stab. matrix eigenvalues (nn=%d, p=%d)',kind,nn,p));
            ylim([min(min(eigenvalues.Real)) 0.2]);
            % xlim([Jplus_store(1) Jplus_store(end)]);
            % eigpic=Params.Files.eigpic;
            % saveas(gcf,eigpic,'pdf')
            % Imag
            subplot(212)
            h(1)=plot(Jplus_store, eigenvalues.Imag(1,:),'color',colors(1,:),'linestyle',':','linewidth',1.1);
            hold on
            for i=1:nn
                h(i)=plot(Jplus_store, eigenvalues.Imag(i,:),'color',colors(i,:));
                hold on
            end
            hold on
            line([Jplus_store(1) Jplus_store(end)],[0 0],'linestyle','--','linewidth',1.2,'color','k');
            % if strcmp(kind,'spur') && ~isempty('J_below')
            %     h(nn+1)=line([J_below(end) J_below(end)],[min(min(rates)) max(max(rates))],'linestyle','--','linewidth',1.2,'color','k');
            %     Legends{nn+1}='J"';
            % end
            hold off
            legend(h,Legends,'location','northeastoutside');
            xlabel('J+ (mV/s)');
            ylabel('Im(eigenvalues)');
            % xlim([Jplus_store(1) Jplus_store(end)]);
            title(sprintf('%s. activity: Stab. matrix eigenvalues (nn=%d, p=%d)',kind,nn,p));
            eigpic=Params.Files.eigpic;
            saveas(gcf,eigpic,'pdf')
            
            
        end
        
        %----------------------------------------------------------------
        %----------------------------------------------------------------
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % function PlotFixedPoints(nn,p,Opt)
        %
        % INPUT: nn = number of populations
        %        Params = structure with current parameters and file info
        %
        % OUTPUT figure(1) = fixed points of selective activity as J+ increases
        %        figure(2) = Re(eigenvalues) of stability matrix as J+ increases
        
        function PlotFixedPoints(nn,paramsfile,Jstop,NumFig,flg)
            
            % global flg.nivsj_mode spur_mode stim_mode
            % preamble
            colors=lines;
            
            %------------------------
            % FILE INFO
            %------------------------
            if flg.nivsj_mode==0 && flg.spur_mode==0
                kind='spont';
            elseif flg.nivsj_mode==1
                kind='sel';
            elseif flg.nivsj_mode==0 && flg.spur_mode==1
                kind='spur';
            end
            if nn==5
                J_below=[];
            end
            
            
            
            %------------------------
            % PLOT FIXED POINTS
            %------------------------
            %
            fntsz=15;
            % file to store fixed points
            load(paramsfile);
            Files=auxMFT.File_Info_MFT(nn,p,flg);
            file_name=Files.FP;
            fprintf('%s\n',file_name);
            load(file_name);
            nn=size(rates,2);
            %
            %
            file_eig=Files.file_eig;
            load(file_eig);
            %
            % stable states
            ind=find(~any(eigenvalues.Real>0));
            figure(NumFig); clf; NumFig=NumFig+1; h=[]; Legends=[];
            
            h(1)=plot(Jplus_store(ind), rates(ind,1),'color',colors(1,:),'linestyle',':','linewidth',2);
            hold on
            for i=2:nn
                h(i)=plot(Jplus_store(ind), rates(ind,i),'color',colors(i,:),'linewidth',2);
                hold on
            end
            
            % legend
            Legends=cell(nn,1);
            Legends{1}='sel. pop';
            for i=2:nn-4
                Legends{i}=sprintf('pop %d',i);
            end
            Legends{nn-3}='bg (Exc)';
            Legends{nn-2}='PV';
            Legends{nn-1}='SST';
            Legends{nn}='VIP';
            % if strcmp(kind,'spur') && ~isempty('J_below')
            %     h(nn+1)=line([J_below(end) J_below(end)],[min(min(rates)) max(max(rates))],'linestyle','--','linewidth',1.2,'color','k');
            %     Legends{nn+1}='J"';
            % end
            legend(h,Legends,'location','northeastoutside');
            tt=sprintf('%s. activity: Fixed points (nn=%d, p=%d)',kind,nn,p);
            xlab='J+ (mV/s)';
            ylab='Firing rate (spks/s)';
            fntsz=15;
            utils.figset(gca,xlab,ylab,tt,fntsz);
            xlim([Jzero Jstop]);
            FPpic=Files.FPpic;
            fprintf('FixedPoint plot saved in %s\n',FPpic);
            saveas(gcf,FPpic,'pdf')
            
            return
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % PLOT EIGENVALUES
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            file_eig=Files.file_eig;
            load(file_eig);
            % Real
            figure(NumFig); clf; h=[];
            subplot(211)
            h(1)=plot(Jplus_store, eigenvalues.Real(1,:),'color',colors(1,:),'linestyle',':','linewidth',3);
            hold on
            for i=2:nn
                h(i)=plot(Jplus_store, eigenvalues.Real(i,:),'color',colors(i,:),'linewidth',3);
                hold on
            end
            hold on
            xlab='J+ (mV/s)';
            ylab='Re(eigenvalues)';
            tt=sprintf('%s. activity: Stab. matrix eigenvalues (nn=%d, p=%d)',kind,nn,p);
            utils.figset(gca,xlab,ylab,tt,fntsz);
            line([Jplus_store(1) Jplus_store(end)],[0 0],'linestyle','--','linewidth',3,'color','k');
            % if strcmp(kind,'spur') && ~isempty('J_below')
            %     h(nn+1)=line([J_below(end) J_below(end)],[min(min(rates)) max(max(rates))],'linestyle','--','linewidth',1.2,'color','k');
            %     Legends{nn+1}='J"';
            % end
            hold off
            legend(h,Legends,'location','northeastoutside');
            
            ylim([min(min(eigenvalues.Real)) 0.2]);
            % xlim([Jplus_store(1) Jplus_store(end)]);
            % eigpic=Files.eigpic;
            % saveas(gcf,eigpic,'pdf')
            % Imag
            subplot(212)
            h=[];
            h(1)=plot(Jplus_store, eigenvalues.Imag(1,:),'color',colors(1,:),'linestyle',':','linewidth',3);
            hold on
            for i=1:nn
                h(i)=plot(Jplus_store, eigenvalues.Imag(i,:),'color',colors(i,:),'linewidth',2);
                hold on
            end
            hold on
            line([Jplus_store(1) Jplus_store(end)],[0 0],'linestyle','--','linewidth',3,'color','k');
            % if strcmp(kind,'spur') && ~isempty('J_below')
            %     h(nn+1)=line([J_below(end) J_below(end)],[min(min(rates)) max(max(rates))],'linestyle','--','linewidth',1.2,'color','k');
            %     Legends{nn+1}='J"';
            % end
            hold off
            legend(h,Legends,'location','northeastoutside');
            xlab='J+ (mV/s)';
            ylab='Im(eigenvalues)';
            tt=sprintf('%s. activity: Stab. matrix eigenvalues (nn=%d, p=%d)',kind,nn,p);
            utils.figset(gca,xlab,ylab,tt,fntsz);
            xlim([Jzero Jstop]);
            eigpic=Files.eigpic;
            saveas(gcf,eigpic,'pdf')
            
        end
        
        %----------------------------------------------------------------
        %----------------------------------------------------------------
        
        
        % function PlotStableStates(p,Params,Jstop,NumFig)
        
        function PlotStableStates(p,Jstop,NumFig,flg)
            
            % global flg.nivsj_mode spur_mode all_mode flg.stim_mode
            
            plotstyle='diamond';%'-';%
            temp_all=flg.all_mode;
            %------
            % LOAD
            %------
            
            paramsfile=flg.paramsfile;
            load(paramsfile);
            nn=p+4;
            Files=auxMFT.File_Info_MFT(nn,p,flg);
            AllData=load(Files.stable);
            
            % set spont mode to retrieve filename
            flg.nivsj_mode=0; flg.spur_mode=0; flg.all_mode=1;
            % nn=3;
            % Files=auxMFT.File_Info_MFT(nn,p,paramsfile);
            % SpontData=load(Files.FP);
            flg.nivsj_mode=1;
            flg.all_mode=temp_all;
            
            % DATA
            % critical points
            Ind=find(cell2mat(arrayfun(@(x)~isempty(x(:).Jplus),AllData.StableStates,'uniformoutput',false)));
            % remove spontaneous activity
            Ind(Ind==p)=[];
            %
            N=numel(Ind);
            JCrit=[];
            YLimRate=100;
            cnt_crit=0;
            for n=Ind
                cnt_crit=cnt_crit+1;
                Jplus=AllData.StableStates(n).Jplus;
                if n==1
                    % max rate for plot
                    YLimRate=max(AllData.StableStates(n).rates(Jplus<Jstop,1));
                end
                JCrit(cnt_crit)=Jplus(1);
            end
            
            % PLOT
            figure(NumFig); clf; h=[]; NumFig=NumFig+1;
            cnt_crit=0;
            % spont
            for n=Ind
                cnt_crit=cnt_crit+1;
                % stable states
                Jplus=AllData.StableStates(n).Jplus;
                rates=AllData.StableStates(n).rates;
                plot(Jplus, rates(:,1),plotstyle,'color',[0 0 0]+(2*(N-cnt_crit))/(2*N+8),'markersize',12);
                hold on;
                % critical points
                line([JCrit(cnt_crit) JCrit(cnt_crit)],[0 YLimRate],'linewidth',1,'linestyle','-.','color','r');
                hold on;
            end
            %     XLim=5.5;
            XLim=Jstop;
            % spontaneous activity
            rates=AllData.StableStates(p).rates;
            Jplus=AllData.StableStates(p).Jplus;
            
            % % Jplus=SpontData.Jplus_store;
            % % rates=SpontData.rates;
            % end of spont activity
            if ~isempty(rates)
                indEnd=find(diff(rates(:,1))>10*mean(diff(rates(:,1))));
                if isempty(indEnd)
                    indEnd=numel(Jplus);
                end
                if ~isempty(indEnd)
                    line([Jplus(indEnd) Jplus(indEnd)],[0 YLimRate],'linewidth',1.5,'linestyle','-.','color','b');
                    hold on;
                    plot(Jplus(1:indEnd), rates(1:indEnd,1),plotstyle,'color','b','markersize',8);
                    hold on;
                    line([Jplus(indEnd) Jplus(indEnd)],[0 YLimRate],'linewidth',1.5,'linestyle','-.','color','b');
                end
            end
            xlab='J+';
            ylab='Firing rate [spks/s]';
            fntsz=15;
            xlim([Jzero XLim]);
            % ylim([0 80]);
            utils.figset(gca,xlab,ylab,'',fntsz);
            
            % SAVE PLOT
            FPpic=Files.Allstablepic;
            saveas(gcf,FPpic,'pdf')
            
            
            %-------------------------------------
            % 2nd PLOT: NUMBER OF ACTIVE CLUSTERS
            %-------------------------------------
            Y=[0:numel(JCrit) numel(JCrit)];
            X=[0 JCrit Jstop];
            XLimStable=[Jzero Jstop];
            figure(NumFig); clf; h=[];
            for n=1:numel(JCrit)
                % critical points
                line([JCrit(n) JCrit(n)],[0 YLimRate],'linewidth',1.5,'linestyle','-.','color','r');
                hold on;
            end
            stairs(X,Y,'k','linewidth',1.5);
            xlim(XLimStable);
            ylim([0 max(Y)+0.5]);
            xlab='J_+';
            ylab='# of active clusters';
            utils.figset(gca,xlab,ylab,'',fntsz);
            % SAVE
            FPpicpdf=[Files.paperstates '.pdf'];
            saveas(gcf,FPpicpdf,'pdf')
            fprintf('\n Figure saved in %s\n',FPpicpdf);
            % FPpicai=[Files.paperstates '.ai'];
            % set(gcf, 'renderer', 'painters');
            % print( gcf, '-painters', FPpicai, '-dill')
            % fprintf('\n Figure saved in %s\n',FPpicai);
            
        end
        
        %----------------------------------------------------------------
        %----------------------------------------------------------------
        
        
        
        function AllData=PlotStableStates_semilogy(p,Jstop,NumFig,flg)
            
            % global nivsj_mode spur_mode all_mode flg.stim_mode
            
            %------
            % LOAD
            %------
            paramsfile=flg.paramsfile;
            load(paramsfile);
            nn=p+4;
            Files=auxMFT.File_Info_MFT(nn,p,flg);
            AllData=load(Files.stable);
            LINE=1.5;
            mrksz=18;
            flg.nivsj_mode=1; flg.spur_mode=0; flg.all_mode=1;
            
            % DATA
            % critical points
            Ind=find(cell2mat(arrayfun(@(x)~isempty(x(:).Jplus),AllData.StableStates,'uniformoutput',false)));
            Ind(Ind==p)=[];
            N=numel(Ind);
            JCrit=[];
            YLimRate=100;
            cnt_crit=0;
            for n=Ind
                cnt_crit=cnt_crit+1;
                Jplus=AllData.StableStates(n).Jplus;
                if n==1
                    % max rate for plot
                    YLimRate=max([YLimRate; AllData.StableStates(n).rates(Jplus<Jstop,1)]);
                end
                JCrit(cnt_crit)=Jplus(1);
            end
            %%
            % PLOT
            figure(NumFig); clf; h=[]; NumFig=NumFig+1;
            cnt_crit=0;
            % spont
            for n=Ind
                cnt_crit=cnt_crit+1;
                % ACTIVE CLUSTERS
                Jplus=AllData.StableStates(n).Jplus;
                rates=AllData.StableStates(n).rates;
                if n==1
                    YLim=max(max(rates))+1;
                end
                % % % % % % % %     % FIX BY HAND
                % % % % % % % %     if n==4
                % % % % % % % %         Jplus=[4.3 Jplus];
                % % % % % % % %         AddRate=[rates_old(2,1)-3.7 rates_old(2,2:end)];
                % % % % % % % %         rates=[AddRate; rates];
                % % % % % % % %     end
                %     plot(Jplus, rates(:,1),'color',[0 0 0]+(2*(N-cnt_crit))/(2*N+10));
                %     hold on;
                semilogy(Jplus, rates(:,1),'-','color',[0 0 0.3]+(2*(N-cnt_crit))/(2*N+10),'linewidth',LINE);
                hold on;
                % Add a point at onset
                semilogy(Jplus(1), rates(1,1),'.','color',[0 0 0.3]+(2*(N-cnt_crit))/(2*N+10),'markersize',mrksz);
                hold on;
                %         % critical points: vertical lines
                %         line([JCrit(cnt_crit) JCrit(cnt_crit)],[0 YLimRate],'linewidth',1.5,'linestyle','-.','color','r');
                %         hold on;
                if 0 % no inactive clusters
                    % INACTIVE CLUSTERS
                    semilogy(Jplus, rates(:,p),'linestyle','--','color',[0 0 0.3]+(2*(N-cnt_crit))/(2*N+10),'linewidth',LINE);
                    hold on;
                    % Add a point at onset
                    semilogy(Jplus(1), rates(1,p),'.','color',[0 0 0.3]+(2*(N-cnt_crit))/(2*N+10),'markersize',mrksz);
                    hold on;
                end
                % NUMBER of active clusters
                %     h=text(Jplus(1)-0.05,rates(1,1)-0.6,['\bf' num2str(n)],'color',...
                %         [0 0 0]+(2*(N-cnt_crit))/(2*N+10),'fontsize',15);
                h=text(Jplus(end),rates(end,1),['\bf' num2str(n)],'color',...
                    [0 0 0]+(2*(N-cnt_crit))/(2*N+10),'fontsize',15);
                %     semilogy(Jplus(1),rates(1,1),'diamond','color',[0 0 0]+(2*(N-cnt_crit))/(2*N+10),'markersize',10);
                hold on;
                % INHIBITORY UNITS
                % PV
                semilogy(Jplus, rates(:,p+2),'-','color',[0.3 0 0]+(2*(N-cnt_crit))/(2*N+10),'linewidth',LINE);
                hold on;
                % Add a point at onset
                semilogy(Jplus(1), rates(1,p+2),'.','color',[0.3 0 0]+(2*(N-cnt_crit))/(2*N+10),'markersize',mrksz);
                hold on;
                % SST
                semilogy(Jplus, rates(:,p+3),':','color',[0.3 0 0]+(2*(N-cnt_crit))/(2*N+10),'linewidth',LINE);
                hold on;
                % Add a point at onset
                semilogy(Jplus(1), rates(1,p+3),'.','color',[0.3 0 0]+(2*(N-cnt_crit))/(2*N+10),'markersize',mrksz);
                hold on;
                % VIP
                semilogy(Jplus, rates(:,p+4),'--','color',[0.3 0 0]+(2*(N-cnt_crit))/(2*N+10),'linewidth',LINE);
                hold on;
                % Add a point at onset
                semilogy(Jplus(1), rates(1,p+4),'.','color',[0.3 0 0]+(2*(N-cnt_crit))/(2*N+10),'markersize',mrksz);
                hold on;
                if 0 % No bg units
                    % BACKGROUND UNITS
                    semilogy(Jplus, rates(:,p+1),'-.','color',[0 0 0.3]+(2*(N-cnt_crit))/(2*N+10),'linewidth',LINE);
                    hold on;
                    % Add a point at onset
                    semilogy(Jplus(1), rates(1,p+1),'.','color',[0 0 0.3]+(2*(N-cnt_crit))/(2*N+10),'markersize',mrksz);
                    hold on;
                    rates_old=rates;
                end
            end
            % if flg.stim_mode==0
            % %     XLim=5.5;
            % elseif flg.stim_mode==1
            %     XLim=Jplus(end);
            % end
            XLim=Jstop;
            
            %%
            %----------------------
            % SPONTANEOUS ACTIVITY
            %----------------------
            % spontaneous activity
            rates=AllData.StableStates(p).rates;
            Jplus=AllData.StableStates(p).Jplus;
            % end of spont activity
            
            if ~isempty(rates)
                indEnd=numel(Jplus);
                if flg.stim_mode==0
                    %         indEnd=find(diff(rates(:,1))>10*mean(diff(rates(:,1))));
                elseif flg.stim_mode==1
                    indEnd=numel(Jplus);
                end
                %     line([Jplus(indEnd) Jplus(indEnd)],[0 YLimRate],'linewidth',1.5,'linestyle','-.','color','b');
                semilogy(Jplus(1:indEnd), rates(1:indEnd,1),'-','color',[0.6 0.4 0.2],'linewidth',LINE);
                % Add a point at onset
                semilogy(Jplus(indEnd), rates(indEnd,1),'.','color',[0.6 0.4 0.2],'markersize',mrksz);
                hold on;
                %     line([Jplus(indEnd) Jplus(indEnd)],[0 YLimRate],'linewidth',1.5,'linestyle','-.','color','b');
            end
            xlab='J+';
            ylab='Firing rate [spks/s]';
            fntsz=22;
            xlim([Jzero XLim]);
            % ylim([0 YLim]);
            ylim([5 30]);
            utils.figset(gca,xlab,ylab,'',fntsz);
            
            %-------------------
            % UNSTABLE STATES
            %-------------------
            
            
            
            
            
            % SAVE PLOT
            FPpic=[Files.paperfig 'LinesSemilogy.pdf'];
            % FPpicai=[Params.Files.paperfig '.ai'];
            saveas(gcf,FPpic,'pdf')
            % set(gcf, 'renderer', 'painters');
            % print( gcf, '-painters', FPpicai, '-dill')
            fprintf('\n Figure saved in %s\n',FPpic);
            
            
        end
        
        
        
        %----------------------------------------------------------------
        %----------------------------------------------------------------
        
        % function PlotStableStates(p,Params,Jstop,NumFig)
        
        function PlotStableStatesSTIM(p,Params,Jstop,NumFig,flg)
            
            % global flg.flg.nivsj_mode spur_mode all_mode flg.flg.stim_mode
            
            % just to be sure
            flg.all_mode=0;
            flg.spur_mode=0;
            TEXT_ON=1;
            
            paramsfile=Params.Files.paramsfile;
            load(paramsfile);
            SelPop=Network.SelPop;
            AllData=load(Params.Files.stable);
            
            % DATA
            % critical points
            %
            % indices with stable states
            % nn=4-> spontaneous activity, where all 1:SelPop pops have same value
            % nn=5... stable states with nn-4 active clusters
            Ind=find(cell2mat(arrayfun(@(x)~isempty(x(:).Jplus),AllData.FP_data(:),'uniformoutput',false)));
            N=numel(Ind);
            JCrit=[];
            cnt_crit=0;
            YRates=[]; % for ylim
            XLim=[];
            for n=Ind'
                cnt_crit=cnt_crit+1;
                Jplus=AllData.FP_data(n).Jplus;
                JCrit(cnt_crit)=Jplus(1);
                rates=AllData.FP_data(n).rates;
                XLim=[XLim Jplus];
                YRates=[YRates max(rates(:,1))];
            end
            XLim=max(XLim);
            YLim=max(YRates);
            % PLOT
            figure(NumFig); clf; h=[]; NumFig=NumFig+1;
            cnt_crit=0;
            % spont
            Legend=[];
            for n=Ind'
                cnt_crit=cnt_crit+1;
                Legend{cnt_crit}=sprintf('%d clust.',n);
                Cols=(2*(N-cnt_crit+1))/(3*N);
                % stable states
                Jplus=AllData.FP_data(n).Jplus;
                rates=AllData.FP_data(n).rates;
                if TEXT_ON
                    text(Jplus',rates(:,1),repmat(num2str(n),numel(Jplus),1),'fontsize',15);
                else
                    h(cnt_crit)=plot(Jplus, rates(:,1),'diamond','color',[0 0 0]+Cols,'markersize',8);%
                end
                hold on;
                % NON-STIMULATED POPULATIONS
                if 1
                    % stim-low rates
                    if n~=SelPop
                        plot(Jplus, rates(:,SelPop),'diamond','color',[Cols Cols 1],'markersize',8);%
                        hold on;
                    end
                    % % %         % non-stim
                    % % %         plot(Jplus, rates(:,end-2),'diamond','color',[Cols 1 Cols],'markersize',6);%
                    % % %         hold on;
                    % inh
                    plot(Jplus, rates(:,end),'o','color',[1 Cols Cols],'markersize',6);%
                    hold on;
                end
                % critical points
                line([JCrit(cnt_crit) JCrit(cnt_crit)],[0 YLim],'linewidth',0.5,'linestyle','-.','color','r');
            end
            % EXTRA VALUES FOR LEGEND
            cnt_crit=cnt_crit+1;
            h(cnt_crit)=plot(Jplus(end), rates(end,SelPop),'diamond','color',[Cols Cols 1],'markersize',8);%
            hold on;
            % % % cnt_crit=cnt_crit+1;
            % % % h(cnt_crit)=plot(Jplus(end), rates(end,end-2),'diamond','color',[Cols 1 Cols],'markersize',6);%
            % % % hold on;
            cnt_crit=cnt_crit+1;
            h(cnt_crit)=plot(Jplus(end), rates(end,end),'o','color',[1 Cols Cols],'markersize',6);%
            Legend=[Legend 'stim.(low)'  'Inh'];
            if ~TEXT_ON
                legend(h,Legend,'location','NorthEastOutside');
            end
            % % % % % spontaneous activity
            % % % % Jplus=AllData.FP_data(4).Jplus;
            % % % % rates=AllData.FP_data(4).rates;
            % end of spont activity
            
            % % % % % % if ~isempty(rates)
            % % % % % %     indEnd=numel(Jplus);
            % % % % % % %     line([Jplus(indEnd) Jplus(indEnd)],[0 YLimRate],'linewidth',1.5,'linestyle','-.','color','k');
            % % % % % % %     hold on;
            % % % % % %     plot(Jplus, rates(:,1),'diamond','color','k','markersize',8);
            % % % % % % %     line([Jplus(indEnd) Jplus(indEnd)],[0 YLimRate],'linewidth',1.5,'linestyle','-.','color','k');
            % % % % % % end
            xlab='J+';
            ylab='Firing rate [spks/s]';
            fntsz=15;
            xlim([Jzero XLim]);
            ylim([0 YLim+5]);
            utils.figset(gca,xlab,ylab,'',fntsz);
            
            
            
            
            
            % SAVE PLOT
            FPpic=Params.Files.Allstablepic;
            saveas(gcf,FPpic,'pdf')
            FPpicai=[Params.Files.Allstablepic(1:end-3) 'ai'];
            set(gcf, 'renderer', 'painters');
            print( gcf, '-painters', FPpicai, '-dill')
            fprintf('\n Figure saved in %s\n',FPpicai);
            
            return
            %-------------------------------------
            % 2nd PLOT: NUMBER OF ACTIVE CLUSTERS
            %-------------------------------------
            Y=[0:numel(JCrit) numel(JCrit)];
            X=[0 JCrit Jstop];
            XLimStable=[Jzero Jstop];
            figure(NumFig); clf; h=[];
            for n=1:numel(JCrit)
                % critical points
                line([JCrit(n) JCrit(n)],[0 YLimRate],'linewidth',1.5,'linestyle','-.','color','r');
                hold on;
            end
            stairs(X,Y,'k','linewidth',1.5);
            xlim(XLimStable);
            ylim([0 max(Y)+0.5]);
            xlab='J_+';
            ylab='# of active clusters';
            utils.figset(gca,xlab,ylab,'',fntsz);
            % SAVE
            FPpicpdf=[Params.Files.paperstates '.pdf'];
            saveas(gcf,FPpicpdf,'pdf')
            FPpicai=[Params.Files.paperstates '.ai'];
            set(gcf, 'renderer', 'painters');
            print( gcf, '-painters', FPpicai, '-dill')
            fprintf('\n Figure saved in %s\n',FPpicai);
            
        end
        
        
        %----------------------------------------------------------------
        %----------------------------------------------------------------
        
        %----------------------------------------------------------------
        %----------------------------------------------------------------
        
        
    end
end