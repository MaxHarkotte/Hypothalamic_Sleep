% Comparison of spindles in the cortex and hypothalamus
% Author: Max Harkotte (maximilian.harkotte@gmail.com)
% Date: June 2024
clear
close all
clc 

%%
id_nlx = '2024-05-21_10-28-00';

%% Paths
local_connect = 'Z:\';
dirProject    = 'Max\05_Downscaling\';
dirDat        = '03_Analysis\01_Pilot\02_Data\01_Downsampled_data\';
dirHypno      = '03_Analysis\01_Pilot\02_Data\02_Sleep_scorings\';

% Add paths to necessary functions for event detection and to fieldtrip
addpath(strcat(local_connect, dirProject, '03_Analysis\00_Pipelines\00_Resources\MATLAB\event_detector\')); % change dir to local if needed

%% Read in recordings
current_dir = strcat(local_connect, dirProject, dirDat, id_nlx);
cd(current_dir)

% EEG
tmp_rec_EEG = load('EEG_EMG_downsampled.mat');
rec_EEG = tmp_rec_EEG.rec_downsample;

% LFP
tmp_rec_LFP = load('LFP_downsampled.mat');
rec_LFP = tmp_rec_LFP.LFP_downsample;

% Concatenate
cfg = [];
data_merged =  ft_appenddata(cfg, rec_EEG, rec_LFP);

clear tmp_rec_EEG tmp_rec_LFP rec_EEG rec_LFP;

%% Filter in spindle band 
fs                   = data_merged.fsample; 
Spi_order            = 3;           % order of butterworth filter 
Spi_fcutlow          = 10;          % low cut frequency in Hz
Spi_fcuthigh         = 16;          % high cut frequency in Hz
[d,c]                = butter(Spi_order,[Spi_fcutlow,Spi_fcuthigh]/(fs/2),'bandpass');
data_spi             = data_merged;

data_spi.trial{1,1}(1,:)  = filtfilt(d,c,data_merged.trial{1,1}(1,:));
data_spi.trial{1,1}(6,:)  = filtfilt(d,c,data_merged.trial{1,1}(6,:));

%% Read in detections
current_dir = strcat(local_connect, dirProject, '03_Analysis\01_Pilot\02_Data\03_Sleep_oscillation_detections\',  id_nlx);
cd(current_dir)

detections = load('Detections.mat', 'detection_output');
detections = detections.detection_output; 
% Channel 1: EEG frontal 
% Channel 2: EEG parietal
% Channel 3: Probe_17 LFP













%% Plotting
% Plot signal timelocked to cortical spindles
for iSpindle = 1:size(detections.spi.events{1,1}, 2)
    subplot(2,2,1)
    plot(data_merged.trial{1,1}(1, detections.spi.events{1,1}(1,iSpindle)-2000:detections.spi.events{1,1}(2,iSpindle)+2000))
    xline(2000, 'Color', 'r')
    xline(2000 + detections.spi.events{1,1}(2,iSpindle)-detections.spi.events{1,1}(1,iSpindle), 'Color', 'r')
    ylim([-400, 400])
    title('Frontal EEG')
    
    subplot(2,2,3)
    plot(data_merged.trial{1,1}(6, detections.spi.events{1,1}(1,iSpindle)-500:detections.spi.events{1,1}(2,iSpindle)+500))
    ylim([-400, 400])
    title('Lateral Hypothalamus')

    subplot(2,2,2)
    plot(data_spi.trial{1,1}(1, detections.spi.events{1,1}(1,iSpindle)-2000:detections.spi.events{1,1}(2,iSpindle)+2000))
    xline(2000, 'Color', 'r')
    xline(2000 + detections.spi.events{1,1}(2,iSpindle)-detections.spi.events{1,1}(1,iSpindle), 'Color', 'r')
    ylim([-400, 400])
    title('Filtered: Frontal EEG')
    
    subplot(2,2,4)
    plot(data_spi.trial{1,1}(6, detections.spi.events{1,1}(1,iSpindle)-500:detections.spi.events{1,1}(2,iSpindle)+500))
    ylim([-400, 400])
    title('Filtered: Lateral Hypothalamus')
    waitforbuttonpress; 

end

% Plot signal timelocked to detections from lateral hypothalamus
for iSpindle = 1:size(detections.spi.events{3,1}, 2)
    subplot(2,2,1)
    plot(data_merged.trial{1,1}(6, detections.spi.events{3,1}(1,iSpindle)-2000:detections.spi.events{3,1}(2,iSpindle)+2000))
    xline(2000, 'Color', 'r')
    xline(2000 + detections.spi.events{3,1}(2,iSpindle)-detections.spi.events{3,1}(1,iSpindle), 'Color', 'r')
    ylim([-400, 400])
    title('Lateral Hypothalamus')
    
    subplot(2,2,3)
    plot(data_merged.trial{1,1}(1, detections.spi.events{3,1}(1,iSpindle)-500:detections.spi.events{3,1}(2,iSpindle)+500))
    ylim([-400, 400])
    title('Frontal EEG')

    subplot(2,2,2)
    plot(data_spi.trial{1,1}(6, detections.spi.events{3,1}(1,iSpindle)-2000:detections.spi.events{3,1}(2,iSpindle)+2000))
    xline(2000, 'Color', 'r')
    xline(2000 + detections.spi.events{3,1}(2,iSpindle)-detections.spi.events{3,1}(1,iSpindle), 'Color', 'r')
    ylim([-400, 400])
    title('Filtered: Lateral Hypothalamus')
    
    subplot(2,2,4)
    plot(data_spi.trial{1,1}(1, detections.spi.events{3,1}(1,iSpindle)-500:detections.spi.events{3,1}(2,iSpindle)+500))
    ylim([-400, 400])
    title('Filtered: Frontal EEG')
    waitforbuttonpress; 
end
