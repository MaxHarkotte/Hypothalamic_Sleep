% Prepare intracranial recordings for further analysis 
% Author: Max Harkotte (maximilian.harkotte@gmail.com)
% last updated: June 2024
clear 
close all
clc 

%% Recording file 
rec_ID_nlx             = '2024-05-21_10-28-00';               

%% Paths
% local_connect = 'Z:\';

serv_connect = '/gpfs01/born/animal/';
dirProject    = 'Max/05_Downscaling/';
dirAnalysis   = '03_Analysis/01_Pilot/';
dirDat        = '02_Raw_Data/01_Pilot/01_Ephys/';
dirPipelines  = '03_Analysis/00_Pipelines/';

% Add paths to necessary functions for reading binary and header info and
% to fieldtrip
addpath(strcat(serv_connect, dirProject, dirPipelines)); % change dir to local if needed
addpath(strcat(serv_connect, dirProject, dirPipelines, '00_Resources/Neuralynx_Import_MEX'));
addpath(strcat(serv_connect, dirProject, '03_Analysis/00_Pipelines/00_Resources/fieldtrip-20210709/fieldtrip-20210709')); % change dir to local if needed

%% References
% Reference table
opts                               = detectImportOptions(strcat(serv_connect, dirProject, dirDat, 'Reference.txt')); % change dir to local if needed
opts.DataLines                     = [1,Inf ]; 
opts.VariableTypes(1, 1:end)       = {'char'};
opts.Delimiter                     = {'\t'};
reference_raw                      = readmatrix(strcat(serv_connect, dirProject, dirDat, 'Reference.txt'),opts); % change dir to local if needed
reference_raw                      = reference_raw(:,1:end);
reference                          = cell2table(reference_raw(2:end,:));
reference.Properties.VariableNames = reference_raw(1,:);
reference.Properties.RowNames      = table2cell(reference(:,1));
clear reference_raw opts;

% General recording info
rec_info   = reference(rec_ID_nlx,:); 
dirRec     = strcat(serv_connect, dirProject, dirDat, rec_ID_nlx); % change dir to local if needed
clear rec_node;

% list all existing channels in the recording 
cd(dirRec)

%% Read in data

% Rereferencing during read in
cfgp            = [];
cfgp.dataset    = strcat(serv_connect, dirProject, dirDat, rec_ID_nlx);
cfgp.channel    = {'Probe_17', 'Probe_11'};
cfgp.reref      = 'yes';
cfgp.refchannel = {'Probe_11'};
cfgp.refmethod  = 'avg';
datp            = ft_preprocessing(cfgp);

%% Prepare mulit unit signal with spike extraction
% MUA filtered
fs                   = datp.fsample; 
MUA_order            = 4;             % order of butterworth filter 
MUA_fcutlow          = 300;           % low cut frequency in Hz
MUA_fcuthigh         = 3000;          % high cut frequency in Hz
[d,c]                = butter(MUA_order,[MUA_fcutlow,MUA_fcuthigh]/(fs/2),'bandpass');
data_MUA             = datp;
data_MUA.trial{1,1}  = filtfilt(d,c,datp.trial{1,1}(2,:));

% Spike extraction
threshold = -40; 
sptimes= find(data_MUA.trial{1,1} < threshold);

%% Prepare LFP for event detection 
% Downsample to 1000 Hz 
cfgr            = [];
cfgr.resamplefs = 1000;
LFP_downsample = ft_resampledata(cfgr, datp);

% Save matlab struct
cd(strcat(serv_connect, dirProject, dirAnalysis, '02_Data/01_Downsampled_data/', rec_ID_nlx)) % change dir to local if needed

save('LFP_downsampled.mat', 'LFP_downsample');
