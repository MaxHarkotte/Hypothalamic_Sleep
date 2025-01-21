% Offline detection of SOs and spindles from surface EEG signal
% Author: Max Harkotte (maximilian.harkotte@gmail.com)
% uses detection toolbox by Niels Niethard and Jens Klinzing (2022)
% Date: January 2025
clear
close all
clc 

% Takes in preprocessed data from HYDO_Preprocessing.m and runs 
% event detection (SO, Spi, Ripples) on all channels

% Requirements: 
% - Signal Processing Toolbox

%% Paths
script_path = which('HYDO01_Event_detection.m');
script_path = strrep(char(script_path), '\', '/');
file_server_path = 'Z:/'; % cluster: '/gpfs01/born/animal/'; 

% Paths to toolboxes and functions
root = strsplit(char(script_path),'src/');
addpath(strcat(char(root(1)), 'src/functions/matlab/event_detector/')); 
addpath(strcat(file_server_path, 'Hypothalamic_Sleep/fieldtrip/fieldtrip-20240722')); 
ft_defaults

clear script_path; 
%% Recording information
reference = readtable(strcat(file_server_path, 'Hypothalamic_Sleep/data/raw/reference.xlsx'));

%% Select recordings for detection 
selection = {'2024-05-21_10-28-00'};   

%% Event detection parameters
cfg_det                        = [];
cfg_det.scoring_epoch_length   = 10; 
cfg_det.code_NREM              = [2 4];
cfg_det.code_REM               = 3;
cfg_det.code_WAKE              = 1;
cfg_det.artfctpad			   = 0;	
cfg_det.spectrum               = 0;					
cfg_det.invertdata             = 0;

% SO detection params
cfg_det.slo		               = 0;					
cfg_det.slo_dur_min		       = 0.5;			
cfg_det.slo_dur_max		       = 2.0;			
cfg_det.slo_freq			   = [0.1 4];
cfg_det.slo_filt_ord	       = 3;
cfg_det.slo_rel_thr            = 33; 
cfg_det.slo_dur_max_down       = 0.300; %in s

% Spindle detection params
cfg_det.spi					   = 1;		
cfg_det.spi_dur_min			   = [0.5 0.25];		
cfg_det.spi_dur_max			   = [2.5 2.5];
cfg_det.spi_thr(1,1)		   = 1.5;
cfg_det.spi_thr(2,1)		   = 2;
cfg_det.spi_thr(3,1)		   = 2.5;
cfg_det.spi_thr_chan		   = [];
cfg_det.spi_freq			   = [10 16];
cfg_det.spi_peakdist_max	   = 0.125;
cfg_det.spi_filt_ord	       = 6;
cfg_det.spi_indiv			   = 0;

% Ripples 
cfg_det.rip                    = 0;

%% Read in recordings and run event detection for each channel

for iRec = 1:size(selection,2)
    
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

    % Hypnogram 
    cd(char(strcat(file_server_path, 'Hypothalamic_Sleep/data/scoring/', rec_info.id, '/', rec_info.file, '/')))

    tmp_hypno       = load(strcat(char(rec_info.file), '_EEGEMG_250Hz.mat'), 'SlStNew'); 
    hypno           = double(tmp_hypno.SlStNew.codes(1:recording_length/10,1));
    cfg_det.scoring = hypno;

    clear tmp_hypno hypno;

    % Run event detection detection
    detection_output = detectEvents_V2(cfg_det, rec);

    % Save detections
    if ~exist(char(strcat(char(root(1)), 'results/detections/', rec_info.id, '/', rec_info.file, '/')), 'dir')
        mkdir(char(strcat(char(root(1)), 'results/detections/', rec_info.id, '/', rec_info.file, '/')))
    end
    
    cd(char(strcat(char(root(1)), 'results/detections/', rec_info.id, '/', rec_info.file, '/')))
    save('Detections','detection_output','-v7.3')

end
