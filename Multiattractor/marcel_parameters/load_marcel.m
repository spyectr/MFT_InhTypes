% load marcel's python files
filename=fullfile('marcel_parameters','weight_disinh.h5');
info = h5info(filename);
trace_datasets = info.Datasets;
dataset_names = {trace_datasets.Name};
data=struct('disinh',cell(1,numel(info.Datasets)),'inh',cell(1,numel(info.Datasets)),'names',dataset_names);
ngroups=length(dataset_names);
for K = 1 : ngroups
    data(K).disinh=double(h5read(filename,strcat('/',dataset_names{K})));
end
figure(1); clf;
[nplot1,nplot2]=nplotxy(ngroups);
for i=1:ngroups
    subplot(nplot1,nplot2,i)
    plot(data(ngroups).disinh,data(i).disinh)
    figset(gca,dataset_names{ngroups},'',dataset_names{i},8);
end
filename=fullfile('marcel_parameters','weight_disinh.pdf');
saveas(gcf,filename,'pdf');

% load marcel's python files
filename=fullfile('marcel_parameters','weight_inh.h5');
info = h5info(filename);
trace_datasets = info.Datasets;
dataset_names = {trace_datasets.Name};
ngroups=length(dataset_names);
for K = 1 : ngroups
    data(K).inh=double(h5read(filename,strcat('/',dataset_names{K})));
end
figure(2); clf;
[nplot1,nplot2]=nplotxy(ngroups);
for i=1:ngroups
    subplot(nplot1,nplot2,i)
    plot(data(ngroups).inh,data(i).inh)
    figset(gca,dataset_names{ngroups},'',dataset_names{i},8);
end
filename=fullfile('marcel_parameters','weight_inh.pdf');
saveas(gcf,filename,'pdf');

data(1).endTraining=360;

filename=fullfile('marcel_parameters','marcel.mat');
save(filename,'data');