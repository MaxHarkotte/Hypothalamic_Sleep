function [trl, event] = Segmentation_TTL(cfg)

% read the header information
hdr           = ft_read_header(cfg.dataset);

% read the events from the data
chanindx = find(ismember(hdr.label, ft_channelselection('TTLSO', hdr.label)));
detectflank   = 'up';
threshold     = 4; % TTl surpasses 4 uV
event         = ft_read_event(cfg.dataset, 'chanindx', chanindx, 'detectflank', detectflank, 'threshold', threshold);

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
