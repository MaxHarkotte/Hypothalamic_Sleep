% Offline detection of SOs and spindles from surface EEG signal
% Author: Max Harkotte (maximilian.harkotte@gmail.com)
% uses detection toolbox by Niels Niethard and Jens Klinzing (2022)
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

% Exclude Channels
cfg = [];
cfg.channel = {'EEG_frontal', 'EEG_parietal', 'Probe_17'}';
data_merged = ft_selectdata(cfg, data_merged);

% recording infos
fs                      = data_merged.fsample; % in Hz
recording_length        = size(data_merged.trial{1,1}, 2); % in seconds
data_merged.sampleinfo    = [1,recording_length];


% Hypnogram 
current_dir = strcat(local_connect, dirProject, dirHypno, id_nlx);
cd(current_dir)

tmp_hypno  = load('EEG_frontal.mat', 'SlStNew'); 
hypno = double(tmp_hypno.SlStNew.codes(1:end,1));

clear tmp_hypno 

%% Event detection
% Scoring params
cfg                        = [];
cfg.name                   = 'Hypothalamus Animal 1';
cfg.scoring                = hypno;
cfg.scoring_epoch_length   = 10; % scoring epoch changed to 2 seconds 
cfg.code_NREM              = [2 4];
cfg.code_REM               = 3;
cfg.code_WAKE              = 1;
% cfg.artfctdef              = artifact_snips;
cfg.artfctpad			   = 0;	
cfg.spectrum               = 1;					
cfg.invertdata             = 0;

% SO detection params
cfg.slo		               = 0;					

% Theta burst detection params
cfg.spi					   = 1;		
cfg.spi_dur_min			   = [0.2 0.2];		
cfg.spi_dur_max			   = [2.5 2.5];
cfg.spi_thr(1,1)		   = 2.5;
cfg.spi_thr(2,1)		   = 2.5;
cfg.spi_thr(3,1)		   = 2.5;
cfg.spi_thr_chan		   = [];
cfg.spi_freq			   = [4 10];
cfg.spi_peakdist_max	   = 0.25;
cfg.spi_filt_ord	       = 6;
cfg.spi_indiv			   = 0;
%cfg.spi_indiv_win		   = 2;
%cfg.spi_indiv_chan		   = [];

% Ripples 
cfg.rip                    = 0;

% Run detection
detection_output = detectEvents_V2(cfg, data_merged);

%% Save detections

if ~exist(strcat(local_connect, dirProject, '03_Analysis\01_Pilot\02_Data\03_Sleep_oscillation_detections\',  id_nlx), 'dir')
    mkdir(strcat(local_connect, dirProject, '03_Analysis\01_Pilot\02_Data\03_Sleep_oscillation_detections\',  id_nlx))
end

cd(strcat(local_connect, dirProject, '03_Analysis\01_Pilot\02_Data\03_Sleep_oscillation_detections\',  id_nlx))

% Plot detections 
plotDetectedEvents(detection_output, 'Detections_theta.fig')

% Save detections
save('Detections_theta','detection_output','-v7.3')

