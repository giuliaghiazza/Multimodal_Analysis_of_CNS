Settings.inputFile = 'DPF_MmAB0001.DAT';
Settings.nChannels = 4096;
Settings.nCurvesTot = 20;

curves = TDNIRS_dataLoad(Settings);
size(curves)
plot(curves(:,1:5))
