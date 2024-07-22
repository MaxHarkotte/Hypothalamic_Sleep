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

clear tmp_rec_EEG tmp_rec_LFP 
%% Read in detections
current_dir = strcat(local_connect, dirProject, '03_Analysis\01_Pilot\02_Data\03_Sleep_oscillation_detections\',  id_nlx);
cd(current_dir)

detections = load('Detections_theta.mat', 'detection_output');
detections = detections.detection_output; 
% Channel 1: EEG frontal 
% Channel 2: EEG parietal
% Channel 3: Probe_17 LFP

%% Define trials -/+ 3 seconds around negative peak of OFFLINE detected SOs
tmp_trials  = detections.spi.events{3,1}; % select Spi start and end from first channel = frontal EEG (3rd channel is LFP)
trl_begin   = tmp_trials(1,:)' - 2*rec_LFP.fsample;
trl_end     = trl_begin + 6*rec_LFP.fsample;
offline_trl = horzcat(trl_begin, trl_end, zeros(size(trl_begin,1),1));

clear tmp_trials trl_begin trl_end; 
%% Select trials of offline detected SOs for all channels 
cfg = [];
cfg.trl = offline_trl;
EEG_theta_segments = ft_redefinetrial(cfg, rec_EEG);
LFP_theta_segments = ft_redefinetrial(cfg, rec_LFP);
%% Calculate time frequency representation (TFR) for spindle band 
% EEG 
cfg              = [];
cfg.output       = 'pow';
cfg.channel      = 'EEG_frontal';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.foi          = 1:0.5:20;                      % analysis 5 to 20 Hz in steps of 0.5 Hz
cfg.t_ftimwin    = [1:6/38:7]./cfg.foi;                   
cfg.toi          = 0:0.05:6;                      % time window "slides" from 0 to 4 sec in steps of 0.05 sec (50 ms)

EEG_TFRhann_theta = ft_freqanalysis(cfg, EEG_theta_segments);

% LFP
cfg              = [];
cfg.output       = 'pow';
cfg.channel      = 'Probe_17';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.foi          = 1:0.5:20;                      % analysis 5 to 20 Hz in steps of 0.5 Hz
cfg.t_ftimwin    = [1:6/38:7]./cfg.foi;                       
cfg.toi          = 0:0.05:6;                      % time window "slides" from 0 to 4 sec in steps of 0.05 sec (50 ms)

LFP_TFRhann_theta = ft_freqanalysis(cfg, LFP_theta_segments);
    
% Baseline normalization of TFR
cfg              = [];
cfg.baseline     = [0.5 1.5];
cfg.baselinetype = 'relative';
EEG_TFRhann_theta_corr = ft_freqbaseline(cfg, EEG_TFRhann_theta); 
LFP_TFRhann_theta_corr = ft_freqbaseline(cfg, LFP_TFRhann_theta); 

%% Plotting
cfg = [];
EEG_theta_grand_avg = ft_freqgrandaverage(cfg, EEG_TFRhann_theta_corr);
LFP_theta_grand_avg = ft_freqgrandaverage(cfg, LFP_TFRhann_theta_corr);

figure;
subplot(2,1,1)
imagesc(EEG_theta_grand_avg.time, EEG_theta_grand_avg.freq, squeeze(EEG_theta_grand_avg.powspctrm(1,:,:)));
axis xy % flip vertically
colormap('jet');
colorbar; 
ylabel('Frequency (Hz)');
title('EEG frontal')

xticks([1 1.5 2 2.5 3 3.5 4 4.5 5]);
xticklabels({'-1', '-0.5', '0', '0.5', '1', '1.5', '2', '2.5', '3'});
xlim([1, 5])

subplot(2,1,2)
imagesc(LFP_theta_grand_avg.time, LFP_theta_grand_avg.freq, squeeze(LFP_theta_grand_avg.powspctrm(1,:,:)));
axis xy % flip vertically
colormap('jet');
colorbar; 
ylabel('Frequency (Hz)');
title('Lateral Hypothalamus LFP')

xticks([1 1.5 2 2.5 3 3.5 4 4.5 5]);
xticklabels({'-1', '-0.5', '0', '0.5', '1', '1.5', '2', '2.5', '3'});
xlim([1, 5])
xlabel('Time (s) relative to Theta burst onset');

