% Visualization spindles detected from preprocessed data
% Author: Max Harkotte (maximilian.harkotte@gmail.com)
% Date: January 2025
clear
close all
clc 

% Takes in preprocessed data from HYDO_Preprocessing.m and event detections
% from HYDO01_Event_detection.m

% Requirements: 

%% Paths
script_path = which('HYDO01_Spindle_vis.m');
script_path = strrep(char(script_path), '\', '/');
root = strsplit(char(script_path),'src/');
file_server_path = 'Z:/'; % cluster: '/gpfs01/born/animal/'; 

% Paths to toolboxes and functions
addpath(strcat(file_server_path, 'Hypothalamic_Sleep/fieldtrip/fieldtrip-20240722')); 
ft_defaults

clear script_path; 
%% Recording information
reference = readtable(strcat(file_server_path, 'Hypothalamic_Sleep/data/raw/reference.xlsx'));

%% Read in recording
selection = {'2024-05-21_10-28-00'};   
iRec = 1;
    
% Recording information
rec_info         = reference(strcmp(reference.file, selection(iRec)), :);
recording_length = rec_info.rec_length*60; % in seconds

% read in data and combine in one ft structure
cd(char(strcat(file_server_path, 'Hypothalamic_Sleep/data/processed/', rec_info.id, '/', rec_info.file, '/')))

cfg         = [];
cfg.dataset = strcat(char(rec_info.file), '_EEGEMG_250Hz.edf');
cfg.channel = {'EEG_frontal', 'EEG_parietal'};
tmp_EEGEMG  = ft_preprocessing(cfg);

cfg         = [];
cfg.dataset = strcat(char(rec_info.file), '_LFP_avg_250Hz.edf');
tmp_LFP     = ft_preprocessing(cfg);

cfg         = [];
rec         =  ft_appenddata(cfg, tmp_EEGEMG, tmp_LFP);

clear cfg tmp_LFP tmp_EEGEMG;

% Cut to correct length
fs             = rec.fsample; % in Hz
cfg            = [];
cfg.begsample  = 1;
cfg.endsample  = recording_length*fs;
rec            = ft_redefinetrial(cfg, rec);
rec.sampleinfo = [1,recording_length*fs];

clear cfg;

%% Read in detections
cd(char(strcat(char(root(1)), 'results/detections/', rec_info.id, '/', rec_info.file, '/')))
load('Detections.mat', 'detection_output');

%% Select channel for visualization
chan    = 'EEG_frontal'; % change accordingly
chan_id = find(ismember(detection_output.info.channel,chan));

% filter in spindle band 
fs               = rec.fsample; % in Hz
data_raw         = rec.trial{1,1};
data_raw         = data_raw(chan_id,:);
spindle_order    = 6;             
spindle_band     = [10 16];

[SpiFilterHigh1,SpiFilterHigh2] = butter(spindle_order ,2*spindle_band(1)/fs,'high');
[SpiFilterLow1,SpiFilterLow2] = butter(spindle_order,2*spindle_band(2)/fs,'low');

SpiBand = filtfilt(SpiFilterHigh1,SpiFilterHigh2,data_raw');
recFilt_Spi = filtfilt(SpiFilterLow1,SpiFilterLow2,SpiBand);
        
% Smoothed Hilbert transformation for threshold and plotting 
recHil_Spi      = smooth(abs(hilbert(recFilt_Spi)),0.1 * fs);

% Divide recording in segments around spindle
spi_bounds    = detection_output.spi.events{chan_id,1};
window_bounds(1,:) = spi_bounds(1,:)-2.5*fs; % start of segment
window_bounds(2,:) = spi_bounds(1,:)+3*fs; % end of segment

for iSpi = 1:size(window_bounds, 2)
    figure; 
    subplot(2,1,1)
    plot(recFilt_Spi(window_bounds(1,iSpi):window_bounds(2,iSpi)))
    hold on;
    plot(recHil_Spi(window_bounds(1,iSpi):window_bounds(2,iSpi)))
    yline(detection_output.spi.thr(1, chan_id), 'r')
    yline(detection_output.spi.thr(2, chan_id), 'r')
    yline(detection_output.spi.thr(3, chan_id), 'r')
    xline(2.5*fs);
    xline(2.5*fs + (spi_bounds(2,iSpi) - spi_bounds(1,iSpi)));

    subplot(2,1,2)
    plot(data_raw(window_bounds(1,iSpi):window_bounds(2,iSpi)))
    xline(2.5*fs);
    xline(2.5*fs + (spi_bounds(2,iSpi) - spi_bounds(1,iSpi)));
    waitforbuttonpress;
    close all
end

