function [stage, solid] = isStage(time, hyp, epoch_length, solid_margin, solid_range)
% Returns the sleep stage for a given time point. Also returns whether this
% time point is solid within that sleep stage (= more than solid_margin
% away from another stage) or any other sleep stage given in solid_range.
%
% INPUT VARIABLES:
% data							time (in sec) for which the sleep stage is 
%								requested
% hyp							hypnogram (num_epochs x 1)
% epoch_length					length of each epoch in hyp (in sec)
% solid_margin					optional; margin (in sec) around time point
%								to check for sleep stage changes (+/- solid_margin; e.g., 1.5)
% solid_range					optional (only required of solid_margin is
%								provided); array of ints (n x 1); sleep
%								stages other than the current one
%								permissable for stage to be solid (e.g.,
%								for checking whether a time point is in
%								solid NREM sleep: [2 3 4])
% OUTPUT VARIABLES:
% stage							sleep stage of provided time point (taken from hyp)
% solid							logical; 1 if sleep stage does not change (or 
%								does not change to stage other than solid_range) 
%								within +/- solid_margin seconds
if nargin < 3
	error('At least three ingredients required.')
end







end