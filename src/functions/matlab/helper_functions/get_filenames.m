function filename = get_filenames(path, varargin)
% Returns a list of files in a given folder, a filename starting with a matching string, 
% or the filename of the nth file in that folder.
%
% Use as
% filename = get_filenames(path, match, fullpath)
% eg.: filename = get_filenames('C:\files\', 5)
%      filelist = get_filenames('C:\files\')
%      filelist = get_filenames('C:\files\', 3, 'full')
%	   filelist = get_filenames('C:\files\', 'subject_1', 'full')
%
% INPUT VARIABLES:
% path              String; input folder
% match             string or int (optional)
%                   Every string other than 'full' will be matched with the
%                   beginning of all file names. Only one matching filename
%                   will be returned. If there is more than one match an
%                   error will be thrown. If there is no match a warning
%                   will be thrown and the ouput will be empty. Instead of
%                   a matching string, the number of the desired file can
%                   be provided (when ordered alphabetically).
% fullpath          string (optional)
%                   If 'full' the full path(s) will be returned
%
% .mff-Files (as provided by some EEG recording systems) are treated as
% files, although windows treats them as folders. DS_Store files that Macs
% leave behind are ignored.
%                   
% OUTPUT VARIABLES:
% filename          names or list of names (in a cell array) of requested
%                   files IN THE ONLY PROPER WAY TO SORT FILENAMES, even
%                   if matlab is of a different opinion!
%                   (a1.mat, a5.mat, a12.mat...)
%                   
% AUTHOR:
% Jens Klinzing, jens.klinzing@uni-tuebingen.de

%% SETUP
all             = dir(path);               % Get the data for the current directory
idx_dirs        = [all.isdir];             % Find the index for directories

% Make sure .mff files are treated as files not folders and Mac nonsense is
% treated as such and is not even mentioned.
names = {all.name};
for i = 1:numel(names)
	if length(names{i}) > 4 && strcmpi(names{i}(end-3:end),'.mff')
		idx_dirs(i) = 0;
	elseif strcmpi(names{i},'.DS_Store') || strcmpi(names{i},'._.DS_Store')
		idx_dirs(i) = 1;
	end
end

fileList        = natsort({all(~idx_dirs).name}');  % ...delete folders out of the file list and sort them reasonably

match           = [];                      % default: dont try to match a string
n               = 0;                       % default: return all files in the directory
full            = false;                   % default: dont add path to filename

if nargin < 1 || nargin > 4
	error('Unexpected number of input arguments.');
elseif nargin > 1 
	if any(strcmp(varargin, 'full'))      
		full = true;
		varargin(strcmp(varargin, 'full')) = []; % makes life easier looking for another string
	end
	if any(cellfun(@ischar, varargin))     % look for match argument
		if sum(cellfun(@ischar, varargin)) == 1
			match = varargin{cellfun(@ischar, varargin)};
		else
			error('More than one matching string provided. Cannot compute.')
		end
	end
	if any(cellfun(@isnumeric, varargin))   % look for n argument
		if numel(varargin{cellfun(@isnumeric, varargin)}) ~= 1
			error('More than one numerical value provided. Cannot compute.')
		elseif ~isempty(match)
			error('You cannot provide a number selector n and a match string at the same time.')
		else
			n = varargin{cellfun(@isnumeric, varargin)};
		end
	end
end

%% START
% Are any selectors given?
if n ~= 0
    fileList = fileList{n};
    disp(['Selecting file ' fileList '.'])    
elseif ~isempty(match)
    for iFile = 1:numel(fileList)
        if length(fileList{iFile}) >= length(match)+1 && (strcmp(fileList{iFile}(1:length(match)+1), [match '_']) || strcmp(fileList{iFile}(1:length(match)+1), [match '.']))
            if n == 0 % if there hasnt been any match yet
                n = iFile;
            else
                error('More than one match was found.')
            end
        end
    end
    if n == 0
        warning(['No file beginning with ''' match ''' was found.'])
        filename = [];
        return
    end
    fileList = fileList{n};
    disp(['Selecting file ' fileList '.'])
end

% Filename or full path?
if full
    filename = fullfile(path, fileList);
else
    filename = fileList;
end


