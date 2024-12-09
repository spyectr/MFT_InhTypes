clear all
% addpath(genpath('/srv/data/luca/neuro/packages/Luca copy'));
name=getComputerName();
recipient='luca.mazzucato@gmail.com';
subject=['MATLAB - ' name];
% myCluster = parcluster('local');
% NumWorkers=myCluster.NumWorkers;
% parpool('local',min(NumWorkers-3,30));
try
    %
    MFTGainModulation_AllInh;
    MFTGainModulation_VIP;
    MFTGainModulation_PV;
    MFTGainModulation_SST;
    MFTGainModulation_Pyr;
    name=getComputerName();
    recipient='luca.mazzucato@gmail.com';
    subject=['MATLAB - ' name];
    maintext='Finished successfully';
    aux.matlabmail(recipient,maintext,subject);
catch theErrorInfo
    maintext={'Error. ',theErrorInfo.message};
    name=getComputerName();
    recipient='luca.mazzucato@gmail.com';
    subject=['MATLAB - ' name];
    aux.matlabmail(recipient,maintext,subject);
    disp(theErrorInfo.message);
end
delete(gcp);
