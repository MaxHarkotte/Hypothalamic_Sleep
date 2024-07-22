function trl = hyp2trl(hyp, epoch_length, stages)
% Turns a hypnogram into a trial structure as expected by FieldTrip. Can be
% used to extract specific sleep stages from a dataset using
% ft_preprocessing. Each trial in the output will represent an episode of
% consecutive epochs belonging to the set of provided sleep stages.
%
% INPUT VARIABLES:
% hyp				hypnogram (num_epochs x 1); one integer for each epoch
% epoch_length		length of an epoch in samples or seconds (e.g., 30 * sampling_frequency)
%					Output will adhere to the units used here. If you want to forward the trl to Fieldtrip use samples!
% stages			stages to create output trials from, e.g. [2 3 4] for NREM
%
% OUTPUT VARIABLES:
% trl				trials (num_episodes x 3), given as [begin end 0].
%					Each trial spans a time window of consecutive epochs provided in 'stages'.
%					trl wil be returned in the same units as epoch_length.
%
% AUTHOR:
% Jens Klinzing, klinzing@princeton.edu

begs			= strfind(any(hyp==stages,2)',[0 1]); % episode onset
ends			= strfind(any(hyp==stages,2)',[1 0]); % episode offset
begs			= begs+1; % because it always finds the epoch before
if any(hyp(1,1)==stages,2), begs = [1 begs]; end % in case recordings starts with a kept stage
if any(hyp(end,1)==stages,2), ends = [ends length(hyp)]; end % in case recordings ends with a kept stage

trl				= [];
trl = [(begs'-1)*epoch_length+1, ends'*epoch_length, zeros(length(ends), 1)]; 
