function MFT_fun_amitmascaro_plot_Stim(tmpv)
options_main=tmpv;


% a=load(paramsfile);



v2struct(options_main);
ncueset=1;
if cue_mode
    if iscell(cueset)
        ncueset=numel(cueset{1});
    else
        ncueset=numel(cueset);
        cueset={cueset};
    end
end
values=repmat(struct('X0',[],'Y0',[],'force',[]),ncueset,numel(stimset));

% global cue_mode all_mode
% all_mode=1;
% if strcmp(Network.cuemode,'none')
%     cue_mode=0;
% else
%     cue_mode=1;
% end

% plot separately
for i_c=1:ncueset
%     cue=cueset(i_c);
    for i=1:numel(cue.cue_option); cue.cue_value{i}=cueset{i}(i_c); options_main.cue.cue_value{i}=cueset{i}(i_c); end
    options_main.extra=tmpv.extra{i_c};
    for i_J=1:numel(stimset)
%         extra='';%['[' cuemode 'cue' num2str(cueset(i_c)) ']']; % filename for plot
        aux.MFT_preamble;
        Jplus=options_main.Jplus;
        save(paramsfile,'Jplus','-append');
%         flg.stim_mode=1;
        flg.all_mode=0;
%         flg.nivsj_mode=1;
        % preamble
        StimHeight=stimset(i_J);
        save(paramsfile,'StimHeight','-append');
        
%         p=params.p;
        nn=6;
        Files=aux.File_Info_MFT(nn,p,flg);
        filesave=[Files.file_amitmascaro '' extra '_results.mat'];
        load(filesave);
        v2struct(OUT);
        nbins_grid=[length(nu_grid(1).focus),length(nu_grid(2).focus)];
        %% PLOT
        numfig=1;
        x1=nu_grid(1).focus; %x1=x1-diff(x1(1:2));
        x2=nu_grid(2).focus; %x2=x2-diff(x2(1:2));
        % ENERGY & FORCE
        % find values of nu_out=phi(nu_in) at each point
        force=OUT.nu_out-OUT.nu_in; % dim= [2D components,x1,x2]
        % energy0=atan(squeeze(sum(temp_diff.^2,1)));
        energy0=sqrt(squeeze(sum(force.^2,1)));
        
        % trajectory:
        % 1) compute trajectory twice running down both i1 and i2
        % 2) find point of max distance from origin, join the two trajectories
        % there
        scan=struct();
        scan(1).path_ind=NaN(2,nbins_grid(1));
        scan(1).path_ind(1,:)=1:nbins_grid(1); % rows=x1 and x2 components; columns= indices along trajectory
        scan(2).path_ind=NaN(2,nbins_grid(2));
        scan(2).path_ind(2,:)=1:nbins_grid(2); % rows=x1 and x2 components; columns= indices along trajectory
        for i1=1:nbins_grid(1)
            [val,indv]=min(energy0(i1,:));
            % while loop below used to avoid getting a minimum at the
            % origin
            cnt=0;
            while indv<=3 && i1==1
                [val,indv]=min(energy0(i1,indv+1:end));
                cnt=cnt+1; if cnt>=nbins_grid(1)-1; break; end
            end
            indv=indv+cnt;
            scan(2).landscape(i1)=val;
            scan(2).path_ind(1,i1)=indv;
        end
        for i2=1:nbins_grid(2)
            [val,indv]=min(energy0(:,i2));
            % while loop below used to avoid getting a minimum at the
            % origin
            cnt=0;
            while indv<=3 && i2==1
                [val,indv]=min(energy0(indv+1:end,i2));
                cnt=cnt+1; if cnt>=nbins_grid(2)-1; break; end
            end
            indv=indv+cnt;
            scan(1).landscape(i2)=val;
            scan(1).path_ind(2,i2)=indv;
        end
        % switch when path intersects diagonal
        [~,ind1]=min(abs(scan(1).path_ind(1,:)-scan(1).path_ind(2,:))); % switch index on nu1
        [~,ind2]=min(abs(scan(2).path_ind(1,:)-ind1)); % switch index on nu2
        % merge two half trajectories on i1 index
        path_ind=[scan(1).path_ind(:,1:ind1) fliplr(scan(2).path_ind(:,1:ind2-1))];
        landscape=[scan(1).landscape(1:ind1) fliplr(scan(2).landscape(1:ind2-1))];
        npoints=size(path_ind,2)-1;
        % trajectory direction at each step
        dx=diff([x1(path_ind(1,:)); x2(path_ind(2,:))],1,2); % rows=x1 and x2 components; columns= indices along trajectory ODD
        dx_norm=dx./repmat(sqrt(sum(dx.^2,1)),2,1); % normalized trajectory vector
        deltax=sqrt(sum(dx.^2,1)); % trajectory increment EVEN
        path_length=cumsum([0 sum(deltax,1)]); % trajectory distance between midbins
        path_length=mean([path_length(1:npoints); path_length(2:npoints+1)],1);
        
        % projection of force onto trajectory
        ind_force=sub2ind(nbins_grid,path_ind(1,:),path_ind(2,:));
        f2=squeeze(force(1,:,:));
        f1=squeeze(force(2,:,:));
        force_comp=[f1(ind_force); f2(ind_force)];
        % average consecutive points to get force in midbin
        force_comp=[mean([force_comp(1,1:npoints); force_comp(1,2:npoints+1)],1);...
            mean([force_comp(2,1:npoints); force_comp(2,2:npoints+1)],1)];
        force_path=dot(force_comp,dx_norm); % force projected on trajectory
        % energy
        energy_path=NaN(1,npoints);
        for i=1:npoints
            energy_path(i)=force_path(1:i)*deltax(1:i)';
        end
        %     return
        
        
        
        
        %--------
        % PLOT 1: Amit-Mascaro potential
        %--------
        figure(numfig); numfig=numfig+1; clf; hold on;
%         energy=log(0.5+energy0);
        energy=log(0.5+energy0);
        Z=energy;
        % % set max and min to interquartile
        % te=energy(:);
        % te=min(te)+(max(te)-min(te))*linspace(1/100,1,100);
        % intq=prctile(te,[1,100]);
        % Z(find(Z<=intq(1)))=intq(1);
        % Z(find(Z>=intq(2)))=intq(2);
        Cmin=min(min(Z));
        Cmax=max(max(Z));
        imagesc(x1,x2,Z); axis xy;
        colormap jet
        caxis([Cmin, Cmax]);
        % equal energy contour
        % contour(x1,x2,energy);
        % force
        % subsample x1 and x2 at 1Hz
        [~,ind11]=min(abs(x1-10));
        [~,ind22]=min(abs(x1-12));
        binstep=ind22-ind11+1;
        x1_arrow=x1(1:binstep:end);
        x2_arrow=x2(1:binstep:end);
        f1_arrow=f1(:,1:binstep:end);
        f1_arrow=f1_arrow(1:binstep:end,:);
        f2_arrow=f2(:,1:binstep:end);
        f2_arrow=f2_arrow(1:binstep:end,:);
        [X1,X2] = meshgrid(x1_arrow,x2_arrow);
        h=quiver(X1,X2,f1_arrow,f2_arrow);
        set(h,'color','k','linewidth',1.5,'MaxHeadSize',0.4);
        
        % % trajectory
        plot(x1(path_ind(1,:)),x2(path_ind(2,:)),'w','linewidth',2);
        
        if 0
            % fixed points
            for i_fp=1:size(StableStates,1)
                a=StableStates(i_fp,ind_focus);
                if (a(1)>=x1(1) && a(1)<=x1(end)) && (a(2)>=x2(1) && a(2)<=x2(end))
                    plot(a(1),a(2),'x','color','r','markersize',10,'linewidth',3);
                end
                % symmetric fixed point
                if (a(2)>=x1(1) && a(2)<=x1(end)) && (a(1)>=x2(1) && a(1)<=x2(end))
                    plot(a(2),a(1),'x','color','r','markersize',10,'linewidth',3);
                end
            end
        end
        % colormap gray;
        % colormap(1-colormap);
        % xlim([nu_grid(2).focus(1) nu_grid(2).focus(end)]);
        % ylim([nu_grid(1).focus(1) nu_grid(1).focus(end)]);
        t=colorbar; get(t,'ylabel');
        set(get(t,'ylabel'),'String', 'Energy = log\Sigma_i(\nu_i^{in}-\nu_i^{out} )^2');
        hold off
        xlab='\nu_2 [spks/s]';
        ylab='\nu_1 [spks/s]';
        tt=sprintf('J+=%g',Jplus);
        fntsz=15;
        figset(gca,xlab,ylab,tt,fntsz);
        
        filepdf=[Files.file_amitmascaro '' extra 'energy2D.pdf'];
        saveas(gcf,filepdf,'pdf');
        
        %%
        %--------
        % PLOT 2: potential along 1D trajectory
        %--------
        
        binwidth=4;
        [path_length, force_path]=gaussfilt(path_length,force_path,2*binwidth,0);
        [path_length, energy_path]=gaussfilt(path_length,energy_path,binwidth,0);
%         path_length=flipud(path_length);
        %         energy_path=flipud(energy_path);
        % Y0=landscape;
        Y0=energy_path-min(energy_path);
        X0=path_length;
        
        % % % chop height to saddle point value on left and right
        % xmin0=findpeaks(-Y0); % local minima
        % Ymax=max(Y0(xmin0.loc(1):xmin0.loc(end)));
        % [~,indhalf0]=min(abs(Y0(1:xmin0.loc(1))-Ymax*1.5));
        % [~,indhalf1]=min(abs(Y0(xmin0.loc(end):end)-Ymax*1.5));
        % Y0=Y0(indhalf0:xmin0.loc(end)+indhalf1-1);
        % X0=X0(indhalf0:xmin0.loc(end)+indhalf1-1);
        %
% % % % % % % % %         % % shift x-axis
% % % % % % % % %         if numel(stimset)==1 && StimHeight==0
% % % % % % % % %             [indhalf, Ybarrier]=peakseek(-abs(force_path)); % middle peak
% % % % % % % % %             X0=X0-X0(indhalf(ceil(numel(indhalf)/2)));
% % % % % % % % %             %             indhalf=findpeaks(-abs(force_path));
% % % % % % % % %             %             X0=X0-X0(indhalf.loc(ceil(numel(indhalf.loc)/2)));
% % % % % % % % %         end
% % % % % % % % %         % plot
% % % % % % % % %         figure(numfig); numfig=numfig+1; clf; hold on;
% % % % % % % % %         x1=X0-X0(1); y1=-force_path+x1;
% % % % % % % % %         h(1)=plot(x1,y1,'b','linewidth',2);
% % % % % % % % %         hold on;
% % % % % % % % %         line([x1(1) x1(end)],[x1(1) x1(end)],'linestyle','--','color','k');
% % % % % % % % %         %         legend(h,'force');
% % % % % % % % %         figset(gca,'\nu [spk/s]','\phi(\nu) [spk/s]','Transfer ftc.',15);
% % % % % % % % %         filepdf=[Files.file_amitmascaro '' extra 'phi.pdf'];
% % % % % % % % %         saveas(gcf,filepdf,'pdf');


        figure(numfig); numfig=numfig+1; clf; hold on;
        indx=sub2ind(size(nu_in),ones(1,numel(path_ind(2,:))),path_ind(1,:),path_ind(2,:));
        indy=sub2ind(size(nu_out),ones(1,numel(path_ind(2,:))),path_ind(1,:),path_ind(2,:));
        x1=squeeze(nu_in(indx)); y1=squeeze(nu_out(indy));
        h(1)=plot(x1,y1,'b','linewidth',2);
        hold on;
        line([x1(1) x1(end)],[x1(1) x1(end)],'linestyle','--','color','k');
        %         legend(h,'force');
        figset(gca,'\nu [spk/s]','\phi(\nu) [spk/s]','Transfer ftc.',15);
        filepdf=[Files.file_amitmascaro '' extra 'phi.pdf'];
        saveas(gcf,filepdf,'pdf');



        figure(numfig); numfig=numfig+1; clf; hold on;
        % plot
        subplot(2,1,1); % force
        h(1)=plot(X0,force_path,'b','linewidth',2);
        hold on;
        line([X0(1) X0(end)],[0 0],'linestyle','--','color','k');
        %         legend(h,'force');
        figset(gca,'Firing rate [spk/s]','Force=\nu-\phi(\nu) [spk/s]','Force',15);
        
        subplot(2,1,2); % energy
        h=plot(X0,Y0,'k');
        set(h,'linewidth',2);
            ylim([min(Y0) max(Y0)*1.1]);
        % plot(X0,force_path);
        figset(gca,'Firing rate [spk/s]','Energy = \int[\nu-\phi(\nu)] ','Potential',15);
        % save
        filepdf=[Files.file_amitmascaro '' extra 'energy1D.pdf'];
        saveas(gcf,filepdf,'pdf');
        
        % store values
        values(i_c,i_J).X0=X0;
        values(i_c,i_J).Y0=Y0;
        values(i_c,i_J).force=-force_path;
        
    end
end

if numel(values)==1
    return;
end
%% overlay plots
%--------
% PLOT overlay potentials and transfer function for different parameters
%--------
if 0%any(cueset==0.05)
    cueset(end)=[];
end

figure(numfig); numfig=numfig+1; clf; hold on;
nplot=max(ncueset*numel(stimset),10);
colors=distinguishable_colors(nplot);
TypePlot={'Y0','force'}; TitlePlot={'Energy = \int[\nu-\phi(\nu)] ','\phi(\nu)-\nu'};
for i_plot=1:numel(TypePlot)
    subplot(2,1,i_plot);
    legs=[];
    plot_cnt=0;
    i_col=0;
    for i_J=1:numel(stimset)
        ColorFactor=1:ncueset;
        if numel(stimset)==1
            Color=[1 0 1];
        elseif ncueset==1
            Color=[0 1 1];
            ColorFactor=1:numel(stimset);
        else
            Color=colors(i_J,:);
        end
        ColorFactor=ColorFactor/max(ColorFactor);
        for i_c=1:ncueset
            plot_cnt=plot_cnt+1;
            X0=values(i_c,i_J).X0;
            Y0=values(i_c,i_J).(TypePlot{i_plot});
            h(plot_cnt)=plot(X0,Y0);
            hold on;
            if numel(stimset)==1 || ncueset==1
                i_col=i_col+1;
            else
                i_col=i_c;
            end
            set(h(plot_cnt),'linewidth',2,'color',Color*ColorFactor(i_col));
            legs{plot_cnt}='[Pert;stim]=[';
            for iopt=1:numel(cue.cue_option)
                legs{plot_cnt}=[legs{plot_cnt} '' cue.cue_option{iopt} '' num2str(cueset{iopt}(i_c))];
            end
            legs{plot_cnt}=[legs{plot_cnt} ';' num2str(stimset(i_J)) ']'];
        end
        figset(gca,'Firing rate [spk/s]',TitlePlot{i_plot},'',12);
        if i_plot==2
            legend(h,legs,'location','northwest');
        end
    end
end
% save
filepdf=[Files.paperenergy 'potentials1D.pdf'];
saveas(gcf,filepdf,'pdf');

%%
%------------------------------------------------
% PROBABILITY OF SEL CLUSTER ACTIVATION, NON SEL CLUSTER INACTIVATION
%------------------------------------------------
% PLOT probabilities
%--------
% activation of selective cluster: jump the barrier from high to low
flg.stim_mode=1;
res=struct();
if flg.stim_mode
    %     logprob_sel2nonsel=NaN(ncueset,numel(stimset)); % prob of transition sel2nonsel clusters (low2high)
    GainTF=NaN(ncueset,numel(stimset)); % gain of transfer function
    Delta_nonsel2sel=NaN(ncueset,numel(stimset)); % prob of transition sel2nonsel clusters (low2high)
    Delta_sel2nonsel=NaN(ncueset,numel(stimset)); % prob of transition sel2nonsel clusters (low2high)
    
    for i_J=1:numel(stimset)
        for i_c=1:ncueset
            % barrier height
%             plot_cnt=plot_cnt+1;
            X0=values(i_c,i_J).X0;
            Y0=values(i_c,i_J).Y0;
            force_path=-values(i_c,i_J).force;
            % local minima
            [indmin, Ymin]=peakseek(-Y0);
%             if numel(indmin)<2% try with the force
%                 [indbarrier1, ~]=peakseek(-abs(force_path)); % middle peak
%                 Ymin=-Y0(indmin)'; indmin=indbarrier1;
%             end
            % NEW 06012020
            if numel(indmin)<2% try with the force
                indbarrier1=find(sign(force_path(2:end).*force_path(1:end-1))<0); % where force changes sign
%                 [indbarrier1, ~]=peakseek(-abs(force_path)); % middle peak
                if numel(indbarrier1)==3
                    Ymin=-Y0(indmin)'; indmin=indbarrier1;
                end
            end

            
            if numel(indmin)>=2
                %                 [Ymin_high,ind_high]=max(-Ymin);
                %                 [Ymin_low,ind_low]=min(-Ymin);
                ind_high=1;
                ind_low=numel(indmin);
                Ymin_high=Y0(indmin(ind_high));
                Ymin_low=Y0(indmin(ind_low));
                % barrier
                if ind_high<ind_low
                    indlow2high=indmin(ind_high):indmin(ind_low);
                else
                    indlow2high=indmin(ind_low):indmin(ind_high);
                end
%                 [indbarrier, Ybarrier]=peakseek(force_path(indlow2high)); % middle peak
                [indbarrier, Ybarrier]=peakseek(Y0(indlow2high)); % middle peak
                
                
                Ypeak=NaN;
                if ~isempty(Ybarrier)
                    Ypeak=max(Ybarrier); % in case of ties
                end
                Eplus=Ypeak-Ymin_high;
                Eminus=Ypeak-Ymin_low;
                Delta_nonsel2sel(i_c,i_J)=Eplus;
                Delta_sel2nonsel(i_c,i_J)=Eminus;
                % gain calculation
                [indforce, forcepeak]=peakseek(abs(force_path)); % middle peak
                x1=X0-min(X0); y1=force_path;
                GainTF(i_c,i_J)=(diff(y1(indforce)))/(abs(diff(x1(indforce))));
            end
        end
    end
    res=struct('Gain',GainTF,'Delta',Delta_nonsel2sel,'cueset',cueset);
    
    
    if numel(stimset)==1 && StimHeight==0
        % CUE ONLY MODE
        figure(numfig); numfig=numfig+1; clf; hold on;
        plot(1:ncueset,Delta_nonsel2sel,'-','color','k','linewidth',3); hold on;
        for i=1:numel(Delta_nonsel2sel)
            h(i)=plot(i,Delta_nonsel2sel(i),'.','color',[1 0 1]*ColorFactor(i),'markersize',50); hold on;
        end
        legend(legs);
        figset(gca,'Pert','\Delta [Hz^2]','barrier height',25);
        % save
        filepdf=[Files.file_amitmascaro '' extra 'BarrierHeight_ongoing.pdf'];
        saveas(gcf,filepdf,'pdf');
    
        % CUE ONLY MODE
        figure(numfig); numfig=numfig+1; clf; hold on;
        plot(1:ncueset,GainTF,'-','color','k','linewidth',3); hold on;
        for i=1:numel(GainTF)
            h(i)=plot(i,GainTF(i),'.','color',[1 0 1]*ColorFactor(i),'markersize',50); hold on;
        end
        figset(gca,'Pert','Gain','',25);
        legend(legs,'fontsize',10);
        % save
        filepdf=[Files.paperenergy '' extra 'Gain_ongoing.pdf'];
        saveas(gcf,filepdf,'pdf');
    
%         % CUE ONLY MODE
%         figure(numfig); numfig=numfig+1; clf; hold on;
% %         plot(GainTF,'-','color','k','linewidth',3);
%         for i=1:numel(GainTF)
%             h(i)=plot(GainTF(i),Delta_nonsel2sel(i),'.','color',[1 0 1]*ColorFactor(i),'markersize',50); hold on;
%         end
%         % linear fit
%         Test=LinearFit(GainTF,Delta_nonsel2sel);
%         plot(GainTF,Test.a+GainTF*Test.b,'linewidth',2,'color','k');% hold on;
%         figset(gca,'Gain','\Delta [Hz^2]',...
%         sprintf('fit:y=%0.03g+%0.03g*x,p=%0.03g,R2=%0.03g',Test.a,Test.b,Test.pLinear,Test.Rsquare),20);
%         % save
%         filepdf=[Files.file_amitmascaro '' extra 'GainBarrier.pdf'];
%         saveas(gcf,filepdf,'pdf');
    end
    try
        if abs(StimHeight)>0
            legs=[];
            % if there's a stimulus, draw time course of probability
            % p_sel-p_nonsel
            %         Delta=Delta_sel2nonsel(1,1);
            %         p_s2ns=exp(-Delta_sel2nonsel/Delta);
            %         p_ns2s=exp(-Delta_nonsel2sel/Delta);
            %         f_sel_m_nonsel=@(t,p0_s2ns,p0_ns2s)(-1+2*(p0_ns2s/(p0_ns2s+p0_s2ns)+...
            %             (0.5-p0_ns2s/(p0_ns2s+p0_s2ns)).*exp(-t*(p0_ns2s+p0_s2ns))));
            % plot
            figure(numfig); numfig=numfig+1; clf; hold on; h=[];
            %         t=0:0.01:1;
            ColorFactor=1:numel(stimset);
            ColorFactor=ColorFactor/numel(stimset);
            LineStyle={'-','--'};
            for i_c=1:ncueset
                h(i_c)=plot(stimset,Delta_nonsel2sel(i_c,:),'color','k','linewidth',3,'linestyle',LineStyle{i_c});
                plot(stimset,Delta_sel2nonsel(i_c,:),'color','c','linewidth',3,'linestyle',LineStyle{i_c});
%                 legs{i_c}=['\sigma_{ext}=' num2str(cueset(i_c))];
                legs{i_c}='[\sigma_{ext}=[';
                for iopt=1:numel(cue.cue_option)
                    legs{i_c}=[legs{i_c} '' cue.cue_option{iopt} '' num2str(cueset{iopt}(i_c))];
                end
                legs{i_c}=[legs{i_c}  ']'];
                if stimset(i_J)>0
                    Color=colors(i_J,:);
                    for i_J=1:numel(stimset)
                        for i=1:numel(stimset)
                            plot(stimset(i_J),Delta_nonsel2sel(i_c,i_J),'.','color',[0 1 1].*ColorFactor(i_J),'markersize',50);
                            plot(stimset(i_J),Delta_sel2nonsel(i_c,i_J),'.','color',[0 1 1].*ColorFactor(i_J),'markersize',50);
                        end
                        %                     if ~isnan(Delta_nonsel2sel(i_c,i_J))
                        %                         h(plot_cnt)=plot(t,f_sel_m_nonsel(t,p_s2ns(i_c,i_J),p_ns2s(i_c,i_J)));
                        %                         set(h(plot_cnt),'color',Color*ColorFactor(i_c),'linewidth',2);
                        %                     end
                    end
                end
            end
            figset(gca,'stim (% baseline)','\Delta [Hz^2]','barrier height',25);
            legend(h,legs,'location','northeast','fontsize',10);
            filepdf=[Files.paperenergy '' extra 'BarrierHeight_evoked.pdf'];
            saveas(gcf,filepdf,'pdf');
        end
    catch; end
end


% filesave=fullfile('figs',[cuemode '_' Opt '.mat']);
save([Files.paperenergy '.mat'],'res');
