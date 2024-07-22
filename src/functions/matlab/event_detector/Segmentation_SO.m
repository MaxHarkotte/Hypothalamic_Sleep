function [trl, event] = Segmentation_SO(cfg)

% read the header information
hdr           = ft_read_header(cfg.dataset);

% read the events from the data
SOs = load(cfg.SOs);
SOs = SOs.detection_output.slo.neg_peaks{1,1}; % select SO troughs from first channel = frontal EEG left
f = {'type', 'sample','value', 'offset', 'duration'};
tmp_string           = cell(1, size(SOs,1));
tmp_string(:)        = {'Trough'};
f{2,1}               = tmp_string; % add triggers for each neg. peak of an SO
f{2,2}               = num2cell(SOs'); % timepoint of neg. peak
f{2,3}               = repmat({ones(1,1)}, 1, size(SOs,1)); % value = 1 
f{2,4}               = repmat({double.empty}, 1, size(SOs,1));
f{2,5}               = repmat({double.empty}, 1, size(SOs,1));
event                = struct(f{:});

% define trials around the events
trl           = [];
pretrig       = cfg.trialdef.pre  * hdr.Fs; % e.g., 1 sec before trigger
posttrig      = cfg.trialdef.post * hdr.Fs; % e.g., 2 sec after trigger

for i = 1:numel(event)
  offset    = -hdr.nSamplesPre;  % number of samples prior to the trigger
  trlbegin  = event(i).sample - pretrig;
  trlend    = event(i).sample + posttrig;
  newtrl    = [trlbegin trlend offset];
  trl       = [trl; newtrl]; % store in the trl matrix
end

