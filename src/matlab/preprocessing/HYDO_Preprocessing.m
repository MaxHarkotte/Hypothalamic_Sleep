% Preprocessing pipeline Hypothalamic Sleep Project
% Author: Max Harkotte (maximilian.harkotte@gmail.com)
% last updated: July 2024

% Takes in raw data and outputs downsampled EEG/EMG and downsampled +
% re-referenced LFP

clear 
close all
clc          

%% Paths
script_path = which('HYDO_Preprocessing.m');
file_server_path = 'Z:/'; % cluster: '/gpfs01/born/animal/'; 

% Paths to toolboxes and functions
addpath(strcat(file_server_path, 'Hypothalamic_Sleep/fieldtrip/fieldtrip-20240722')); 
ft_defaults

%% Recording information
reference = readtable(strcat(file_server_path, 'Hypothalamic_Sleep/data/raw/reference.xlsx'));

%% Select recordings for preprocessing 
selection             = {'2024-05-21_10-28-00', '2024-07-22_10-47-00'};    

%% Preprocessing
% Prepare .edf file to write data to 
cfg            = [];
cfg.dataset    = strcat(strcat(file_server_path, 'Hypothalamic_Sleep/data/edfexample/Osas2002.edf')); % change dir to local if needed
cfg.continuous = 'yes';
cfg.channel    = 'EEG Fpz-M2';
ft_dummy_dat   = ft_preprocessing(cfg); 

% preprocess recordings 
for iRec = 1:size(selection,2)

    % Recording information
    rec_info = reference(strcmp(reference.file, selection(iRec)), :);
    elecs = dir(fullfile(char(strcat(file_server_path, 'Hypothalamic_Sleep/data/raw/', rec_info.id, '/', rec_info.file, '/')),  '*.ncs'));
    elecs = {elecs.name}';
    elecs = erase(elecs, '.ncs');

    % Read in EEG and EMG 
    EEGEMG_chans = startsWith(elecs, 'E');
    EEGEMG_chans = elecs(EEGEMG_chans);

    for iChan = 1:size(EEGEMG_chans, 1)
        % read in data 
        cfgp         = [];
        cfgp.dataset = char(strcat(file_server_path, 'Hypothalamic_Sleep/data/raw/', rec_info.id, '/', rec_info.file, '/'));
        cfgp.channel = EEGEMG_chans(iChan);
        datp         = ft_preprocessing(cfgp);
    
        % Downsample to 250 Hz 
        cfgr            = [];
        cfgr.resamplefs = 250;
        datr{iChan}     = ft_resampledata(cfgr, datp);

        clear datp;
    end

    cfg               = [];
    EEGEMG_downsample = ft_appenddata(cfg, datr{:}); 

    clear datr; 

    % Read in LFP
    Probe_chans = startsWith(elecs, 'P');
    Probe_chans = elecs(Probe_chans);

    for iChan = 1:size(Probe_chans, 1)
        % read in data 
        cfgp         = [];
        cfgp.dataset = char(strcat(file_server_path, 'Hypothalamic_Sleep/data/raw/', rec_info.id, '/', rec_info.file, '/'));
        cfgp.channel = Probe_chans(iChan);
        datp         = ft_preprocessing(cfgp);
    
        % Downsample to 250 Hz 
        cfgr            = [];
        cfgr.resamplefs = 250;
        datr{iChan}     = ft_resampledata(cfgr, datp);

        clear datp;
    end

    cfg               = [];
    LFP_downsample = ft_appenddata(cfg, datr{:}); 

    clear datr; 

    % Re-reference channels of Probe 
    cfg = [];
    cfg.channel = 'all'; 
    cfg.reref = 'yes';
    cfg.refmethod = 'avg';
    cfg.refchannel = 'all';
    LFP_reref = ft_preprocessing(cfg, LFP_downsample);

    % Write data as .edf file 
    cfg       = [];
    rec_clean = ft_appenddata(cfg, EEGEMG_downsample, LFP_reref);

    % Change header information
    



end



