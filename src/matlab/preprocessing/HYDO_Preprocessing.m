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

%% Channel map for A1x32-Poly3-10mm-50-177
channel_neighbors = containers.Map;
channel_neighbors('Probe_2') = {'Probe_1', 'Probe_16'};
channel_neighbors('Probe_1') = {'Probe_2', 'Probe_3', 'Probe_18'};
channel_neighbors('Probe_3') = {'Probe_1', 'Probe_4', 'Probe_15'};
channel_neighbors('Probe_4') = {'Probe_3', 'Probe_5', 'Probe_19'};
channel_neighbors('Probe_5') = {'Probe_4', 'Probe_6', 'Probe_14'};
channel_neighbors('Probe_6') = {'Probe_5', 'Probe_7', 'Probe_20'};
channel_neighbors('Probe_7') = {'Probe_6', 'Probe_8', 'Probe_13'};
channel_neighbors('Probe_8') = {'Probe_7', 'Probe_9', 'Probe_21'};
channel_neighbors('Probe_9') = {'Probe_8', 'Probe_10', 'Probe_12'};
channel_neighbors('Probe_10') = {'Probe_9', 'Probe_22'};
channel_neighbors('Probe_17') = {'Probe_16'};
channel_neighbors('Probe_16') = {'Probe_2', 'Probe_17', 'Probe_18', 'Probe_31'};
channel_neighbors('Probe_18') = {'Probe_1', 'Probe_15', 'Probe_16', 'Probe_32'};
channel_neighbors('Probe_15') = {'Probe_3', 'Probe_18', 'Probe_19', 'Probe_30'};
channel_neighbors('Probe_19') = {'Probe_4', 'Probe_14', 'Probe_15', 'Probe_29'};
channel_neighbors('Probe_14') = {'Probe_5', 'Probe_19', 'Probe_20', 'Probe_28'};
channel_neighbors('Probe_20') = {'Probe_6', 'Probe_13', 'Probe_14', 'Probe_27'};
channel_neighbors('Probe_13') = {'Probe_7', 'Probe_20', 'Probe_21', 'Probe_26'};
channel_neighbors('Probe_21') = {'Probe_8', 'Probe_12', 'Probe_13', 'Probe_25'};
channel_neighbors('Probe_12') = {'Probe_9', 'Probe_21', 'Probe_22', 'Probe_24'};
channel_neighbors('Probe_22') = {'Probe_10', 'Probe_11', 'Probe_12', 'Probe_23'};
channel_neighbors('Probe_11') = {'Probe_22'};
channel_neighbors('Probe_31') = {'Probe_16', 'Probe_32'};
channel_neighbors('Probe_32') = {'Probe_18', 'Probe_30', 'Probe_31'};
channel_neighbors('Probe_30') = {'Probe_15', 'Probe_29', 'Probe_32'};
channel_neighbors('Probe_29') = {'Probe_19', 'Probe_28', 'Probe_30'};
channel_neighbors('Probe_28') = {'Probe_14', 'Probe_27', 'Probe_29'};
channel_neighbors('Probe_27') = {'Probe_20', 'Probe_26', 'Probe_28'};
channel_neighbors('Probe_26') = {'Probe_13', 'Probe_25', 'Probe_27'};
channel_neighbors('Probe_25') = {'Probe_21', 'Probe_24', 'Probe_26'};
channel_neighbors('Probe_24') = {'Probe_12', 'Probe_23', 'Probe_25'};
channel_neighbors('Probe_23') = {'Probe_22', 'Probe_24'};

%% Preprocessing
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

    % Re-reference channels of Probe to average of all elecs
    cfg = [];
    cfg.channel = 'all'; 
    cfg.reref = 'yes';
    cfg.refmethod = 'avg';
    cfg.refchannel = 'all';
    LFP_avg = ft_preprocessing(cfg, LFP_downsample);

    % Re-reference channels of Probe locally using bipolar/laplace scheme 
    LFP_bipolar = LFP_downsample;
    exclude_channels = {};
    iExcl = 1;

    for iChan = 1:size(Probe_chans, 1)
        % list possible reference channels
        neighbors = channel_neighbors(char(Probe_chans(iChan)));
        [isPresent, idx] = ismember(neighbors, LFP_downsample.label);
       
        % Re-reference with the average of all existing channels that are
        % 50 microns away from the elec
        if isPresent
            LFP_bipolar.trial{1,1}(iChan,:) = LFP_downsample.trial{1,1}(iChan, :) - mean(LFP_downsample.trial{1,1}(idx(isPresent),:),1);
        else
            exclude_channels(iExcl) = Probe_chans(iChan);
            iExcl = iExcl +1;
        end
    end

    % Exclude channels for bipolar re-referencing, if no neighoring channel
    % was available
    if ~isempty(exclude_channels)
        exclude_channels = strcat('-', exclude_channels);
        cfg = [];
        cfg.channel = [{'all'}, exclude_channels]; 
        LFP_bipolar = ft_selectdata(cfg, LFP_bipolar);
    end

    % Create general header info 
    hdr.Fs = 250; 
    hdr.nSamples = size(EEGEMG_downsample.trial{1,1},2);
    hdr.SamplesPre = 0;
    hdr.nTrials = 1;

    % create channel specific header info 
    EEGEMG_downsample.hdr          = hdr; 
    EEGEMG_downsample.hdr.label    = EEGEMG_downsample.label; 
    EEGEMG_downsample.hdr.nCans    = size(EEGEMG_downsample.label,1);
    EEGEMG_downsample.hdr.chantype = EEGEMG_downsample.label;
    units                          = cell(size(EEGEMG_downsample.label,1), 1);
    units(:)                       = {'uV'};
    EEGEMG_downsample.hdr.chanunit = units;
    hdr_EEGEMG                     = ft_fetch_header(EEGEMG_downsample);
    clear units; 

    LFP_avg.hdr          = hdr; 
    LFP_avg.hdr.label    = LFP_avg.label; 
    LFP_avg.hdr.nCans    = size(LFP_avg.label,1);
    LFP_avg.hdr.chantype = LFP_avg.label;
    units                = cell(size(LFP_avg.label,1), 1);
    units(:)             = {'uV'};
    LFP_avg.hdr.chanunit = units;
    hdr_LFP_avg          = ft_fetch_header(LFP_avg);
    clear units;

    LFP_bipolar.hdr          = hdr; 
    LFP_bipolar.hdr.label    = LFP_bipolar.label; 
    LFP_bipolar.hdr.nCans    = size(LFP_bipolar.label,1);
    LFP_bipolar.hdr.chantype = LFP_bipolar.label;
    units                    = cell(size(LFP_bipolar.label,1), 1);
    units(:)                 = {'uV'};
    LFP_bipolar.hdr.chanunit = units;
    hdr_LFP_bipolar          = ft_fetch_header(LFP_bipolar);
    clear units; 

    % create directory to store downsampled data
    if ~exist(char(strcat(file_server_path, 'Hypothalamic_Sleep/data/processed/', rec_info.id, '/', rec_info.file, '/')), 'dir')
        mkdir(char(strcat(file_server_path, 'Hypothalamic_Sleep/data/processed/', rec_info.id, '/', rec_info.file, '/')))
    end
    
    % create directory to store sleep scorings (done manually later)
    if ~exist(char(strcat(file_server_path, 'Hypothalamic_Sleep/data/scoring/', rec_info.id, '/', rec_info.file, '/')), 'dir')
        mkdir(char(strcat(file_server_path, 'Hypothalamic_Sleep/data/scoring/', rec_info.id, '/', rec_info.file, '/')))
    end
    
    % Write .edf files 
    cd(char(strcat(file_server_path, 'Hypothalamic_Sleep/data/processed/', rec_info.id, '/', rec_info.file, '/')))
    
    % EEG and EMG channels
    ft_write_data(strcat(char(rec_info.file), '_EEGEMG_250Hz.edf') ...
    , EEGEMG_downsample.trial{1}, 'header', hdr_EEGEMG,'chaninx', EEGEMG_downsample.label', 'dataformat', "edf")
    
    % LFP channels rereferenced to average of all probe elecs
    ft_write_data(strcat(char(rec_info.file), '_LFP_avg_250Hz.edf') ...
    , LFP_avg.trial{1}, 'header', hdr_LFP_avg,'chaninx', LFP_avg.label', 'dataformat', "edf")

    % LFP channels rereferenced with laplace/biploar method
    ft_write_data(strcat(char(rec_info.file), '_LFP_bipolar_250Hz.edf') ...
    , LFP_bipolar.trial{1}, 'header', hdr_LFP_bipolar,'chaninx', LFP_bipolar.label', 'dataformat', "edf")

end



