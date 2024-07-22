function varargout = trial_selector(varargin)
% GUI allowing to manually exclude trials from a fieldtrip structure. Shows
% 20 trials at once. Input data is expected to be a standard fieldtrip
% dataset with one channel and the fields 'trials' and 'time'.
%
% Use as:
%  output = trial_selector(cfg)
%  output = trial_selector(cfg, data)
%  trial_selector(cfg, data)
%
% INPUT VARIABLES:
% cfg               structure; configuration with fields:
% .inputfile		string (optional); path to dataset
%					alternatively, data can be provided as a second input
%					argument
% .outputfile		string (optional); path to (potentially existing)
%					output file
% .page				int; on which page to start on (default: 1)
% .latency			int array; time limits for displaying trials
%					(e.g. [-1 4]
% .comment			cell array of strings with as many entries as there are
%					trials (optional); will be shown in the corner of each
%					plot; could for instance be used to show trial
%					conditions
% .selected			array of logicals with as many entries as there are
%					trials; pre-selects those trials (for instance from a
%					previous run. Alternatively, existing selections can be
%					provided as a file using:
% .load_selected 	logical (optional); will try to load existing
%					selections from cfg.outputfile (default: false)
%
% OUTPUT VARIABLES:
% output			array of logicals; is only returned if cfg.outputfile
%					is not specified or empty
%
% COMMENTS:
% All important stuff happens in:
% trial_selector_OpeningFcn		runs once at the beginning
% initializePage				is called for displaying each new page
%
% WANT TO CONTRIBUTE?
% Cool! Here are some suggestions for things to do (just send me a pull
% request on github):
% . currently the function can only handle single-channel input. this could
%   be extended to (lets say) 5 channels that are shown in a butterfly plot
%   with a nice color scheme.
% . the output to command line is implemented an a dodgy way. matlab did
%   not make it easy for me to have that output independent of whether the
%   user clicked the X to close the window or the 'close' button. if you
%   know a more elegant solution, please contribute!
% . there might be further tasks hidden in the comments (search for
%   'TODO')
% . currently two vertical lines are displayed in each plot, one at -.5 and
%	one at 0. this is motivated by my own project. it'd be preferably if
%	those were input parameters (e.g. an array of x-axis values and an
%	array of color definitions of the same length?)
%
% AUTHOR:
% Jens Klinzing, jens.klinzing@uni-tuebingen.de

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
	'gui_Singleton',  gui_Singleton, ...
	'gui_OpeningFcn', @trial_selector_OpeningFcn, ...
	'gui_OutputFcn',  @trial_selector_OutputFcn, ...
	'gui_LayoutFcn',  [] , ...
	'gui_Callback',   []);
if nargin && ischar(varargin{1})
	gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
	[varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
	gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before trial_selector is made visible.
function trial_selector_OpeningFcn(hObject, ~, handles, varargin)

% For debugging
% if numel(varargin) == 0
% 	varargin{1} = load_file('debug_cfg.mat');
% 	varargin{2} = varargin{1}.data;
% 	varargin{1} = rmfield(varargin{1}, 'data');
% end
% hObject.Visible = 'on';

% Make all fields in input available
if numel(varargin) > 0 && isstruct(varargin{1})
	f = fields(varargin{1});
	for iField = 1:numel(f)
		handles.(f{iField}) = varargin{1}.(f{iField});
	end
else
	error('No configuration given as input.')
end

% Set some settings and defaults
handles.color_selected = [1 .7 .7]; % redish
handles.color_unselected = [1 1 1]; % white
if ~isfield(handles, 'page'), handles.page = 1; end
if ~isfield(handles, 'ylim'), handles.ylim = [-150 150]; end
if ~isfield(handles, 'xlim'), handles.xlim = [-2 6]; end
if ~isfield(handles, 'load_selected'), handles.load_selected = false; end
if isfield(handles, 'inputfile') && isempty(handles.inputfile)
	handles = rmfield(handles, 'inputfile');
end
if isfield(handles, 'outputfile') && isempty(handles.outputfile)
	handles = rmfield(handles, 'outputfile');
end
if isfield(handles, 'selected') && isempty(handles.selected)
	handles = rmfield(handles, 'selected');
end

% Load and check data
if numel(varargin) == 2
	if isfield(handles, 'inputfile')
		error('You cannot provide both cfg.inputfile and data as input.')
	else
		handles.data = varargin{2};
	end
else
	temp = load(handles.inputfile);
	names = fieldnames(temp);
	if length(names) ~= 1
		error('Unexpected content in preprocessed data file. File should contain one single structure.')
	end
	handles.data = temp.(names{1});
	clear temp names
end
requiredFields = {'trial', 'time'};
for i = requiredFields
	if ~isfield(handles.data,i)
		error(['Required field missing in data: ' i{1} '.']);
	end
end
if size(handles.data.trial, 1) ~= 1
	error('Data should have only one channel.')
end

% Check other input
if isfield(handles, 'selected')
	if handles.load_selected
		error('You cannot try to load an existing selection file and provide cfg.selected')
	else
		handles.output = handles.selected;
	end
end
if isfield(handles, 'outputfile')
	% Test whether existing selections exist and fit the data
	if handles.load_selected
		if exist(handles.outputfile, 'file')
			temp = load(handles.outputfile);
			names = fieldnames(temp);
			if length(names) ~= 1
				error('Unexpected content in existing trial selection file. File should contain one single structure.')
			end
			handles.output = temp.(names{1});
			clear temp names
			if ~islogical(handles.output) || size(handles.output, 2) ~= 1
				error('Something is wrong with existing trial selections.')
			elseif size(handles.output, 1) ~= size(handles.data.trial,2)
				error('Existing trial selections do not fit the provided data.')
			end
		else
			error('Cannot load existing selections, cfg.outputfile does not exist.')
		end
	end
	% Test whether we will be able to write that file
	[fid,errmsg] = fopen(handles.outputfile, 'a');
	if ~isempty(errmsg)&&strcmp(errmsg,'Permission denied')
		error('You do not have write permission to outputfile (%s).', handles.outputfile);
	elseif ~isempty(errmsg)&&strcmp(errmsg,'No such file or directory')
		error('Outputfile directory does not exist (%s).', handles.outputfile);
	else
		fclose(fid);
	end
end
if isfield(handles, 'comment')
	if isempty(handles.comment)
		handles = rmfield(handles, 'comment');
	elseif ~iscellstr(handles.comment) || length(handles.comment) ~= size(handles.data.trial,2)
		error('Something is wrong with the comments you provided.')
	end
end

% Bookkeeping
handles.ax_handles	= {handles.axes1, handles.axes2, handles.axes3, handles.axes4, handles.axes5, handles.axes6, handles.axes7, handles.axes8, handles.axes9, handles.axes10, handles.axes11, handles.axes12, handles.axes13, handles.axes14, handles.axes15, handles.axes16, handles.axes17, handles.axes18, handles.axes19, handles.axes20}; % TODO: Add all
handles.pl_handles	= cell(numel(handles.ax_handles));
handles.num_axes	= numel(handles.ax_handles);
handles.num_trials	= size(handles.data.trial,2);
handles.num_pages	= ceil(handles.num_trials / handles.num_axes);
if ~isfield(handles, 'output')
	handles.output	= false(handles.num_trials,1); % eventual output
end
if handles.page <= 0 || handles.page > handles.num_pages || mod(handles.page, 1) ~= 0
	error('Invalid cfg.page.')
end

% Bring page selection in order
set(handles.pagetext2, 'String', [' / ' num2str(handles.num_pages)]);

% Make sure all shenanigans are save in the almighty cloud
handles.initialpage = true; % the first page is handled a bit differently than the others..
guidata(hObject, handles);

% Trigger the first page (mostly plotting)
initializePage(hObject, handles)

% Output voodoo in case we'll return result to user
if ~isfield(handles, 'outputfile'), uiwait(handles.figure1); end


% --- Shows a new page (potting, bookkeeping etc.)
function initializePage(hObject, handles)
% Determine what to plot on the current page
handles.ax_trial		= zeros(20,1);
first_trial				= (handles.page-1) * handles.num_axes + 1;
last_trial				= handles.page * handles.num_axes;
if last_trial > handles.num_trials, last_trial = handles.num_trials; end
cnt = 1;
for iTrial = first_trial:last_trial
	handles.ax_trial(cnt)	= iTrial;
	cnt						= cnt + 1;
end

% Update buttons
set(handles.pagesel, 'String', num2str(handles.page));
if handles.page == 1
	handles.previous.Enable		= 'off';
	handles.next.Enable			= 'on';
elseif handles.page == handles.num_pages
	handles.next.Enable			= 'off';
	handles.previous.Enable		= 'on';
else
	handles.next.Enable			= 'on';
	handles.previous.Enable		= 'on';
end

% Do the actual plotting
for iAx = 1:handles.num_axes
	ax		= handles.ax_handles{iAx}; % current axis
	tr		= handles.ax_trial(iAx);   % which trial to plot
	
	% Create line plots
	if tr ~= 0 % if plot is not empty
		if handles.initialpage  % if it is the first page we need to create figure handles
			[~,l1]	= min(abs(handles.data.time{tr} - handles.xlim(1))); % lower time limit
			[~,l2]	= min(abs(handles.data.time{tr} - handles.xlim(2))); % upper time limit
			handles.pl_handles{iAx} = plot(ax, handles.data.time{tr}(l1:l2), handles.data.trial{tr}(l1:l2));
			hold(ax, 'on')
			pl2				= plot(ax, [0 0], handles.ylim); % vertical line at 0
			pl3				= plot(ax, [-.5 -.5], handles.ylim); % vertical line at -.5
		else
			[~,l1]	= min(abs(handles.data.time{tr} - handles.xlim(1)));  % lower time limit
			[~,l2]	= min(abs(handles.data.time{tr} - handles.xlim(2)));  % upper time limit
			handles.pl_handles{iAx}.XData = handles.data.time{tr}(l1:l2); % faster than calling plot() again
			handles.pl_handles{iAx}.YData = handles.data.trial{tr}(l1:l2);
		end
	else % in case trial is empty we still need to plot something so that axes are correct
		handles.pl_handles{iAx} = plot(ax, [handles.data.time{1}(l1) handles.data.time{1}(l2)],[0 0]); % take over time limits from first trial
		hold(ax, 'on')
		pl2				= plot(ax, [0 0], handles.ylim); % vertical line at 0
		pl3				= plot(ax, [-.5 -.5], handles.ylim); % vertical line at -.5
	end
	
	handles.pl_handles{iAx}.LineWidth = 1.0;
	handles.pl_handles{iAx}.Color	= 'k';
	set(handles.pl_handles{iAx},'HitTest','off') % otherwise those plots will register a button down

	if exist('pl2', 'var')
		pl2.LineWidth	= .5;
		pl2.Color		= 'r';
		set(pl2,'HitTest','off') %... and not our axes
	end
	if exist('pl3', 'var')
		pl3.LineWidth	= .5;
		pl3.Color		= [.6 .6 .6];
		set(pl3,'HitTest','off')
	end
	
	% Behavior when clicked
	if tr ~= 0
		set(ax,'ButtonDownFcn',{@toggleAxes, iAx}) % for the axes, please call our own button-down function
	else
		set(ax,'ButtonDownFcn','')
	end
	
	% Add the annotations
	if isfield(handles, 'comment') && ~isempty(handles.comment)
		current_annotag = ['anno' num2str(iAx)];
		if handles.initialpage
			if tr ~= 0
				annotation('textbox',ax.Position,'String', handles.comment{tr},'FitBoxToText','on', 'EdgeColor','none', 'Tag', current_annotag); % + [.001 -.001 0 0]
			else
				annotation('textbox',ax.Position,'String', [],'FitBoxToText','on', 'EdgeColor','none', 'Tag', current_annotag); % + [.001 -.001 0 0]
			end
		else
			an = findall(gcf,'Tag',current_annotag);
			if tr ~= 0
				an.String = handles.comment{tr};
			else
				an.String = [];
			end
		end
	end
	
	% Set some axes properties
	ax.YLim			= handles.ylim;
	ax.YTickLabel	= '';
	if tr ~= 0
		if handles.output(handles.ax_trial(iAx))
			handles.ax_handles{iAx}.Color = handles.color_selected;
		else
			handles.ax_handles{iAx}.Color = handles.color_unselected;
		end
	end
	if tr ~= 0 && iAx > numel(handles.ax_trial) - 5 && iAx <= numel(handles.ax_trial)
		ax.XTick		= handles.data.time{tr}(l1)-rem(handles.data.time{tr}(l1),2):2:handles.data.time{tr}(l2)-rem(handles.data.time{tr}(l2),2); % :P
	elseif iAx > numel(handles.ax_trial) - 5 && iAx <= numel(handles.ax_trial)
		ax.XTick		= handles.data.time{1}(l1)-rem(handles.data.time{1}(l1),2):2:handles.data.time{1}(l2)-rem(handles.data.time{1}(l2),2); % :P
	else
		ax.XTick		= [];
	end
	
	hold(ax, 'off')
end
handles.initialpage = false;
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = trial_selector_OutputFcn(hObject, ~, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure. We need to check
% whether outputfile is provided because a) nargout is not available within
% a GUI for some reason (nargout is always 0) and b) because this output
% function is called right away if uiwait is not enabled (absolute
% nonsense, thanks matlab).
% If you know of a better solution to the whole output disaster, please let
% me know!
if ~isfield(handles, 'outputfile')
	varargout{1} = handles.output;
	delete(hObject);
end


% --- Selects/unselects plots.
function toggleAxes(hObject, ~, iAx)
handles = guidata(hObject); % handing over handles as input parameter somehow did not work
if ~handles.output(handles.ax_trial(iAx))
	handles.ax_handles{iAx}.Color			= handles.color_selected;
	handles.output(handles.ax_trial(iAx))	= true;
else
	handles.ax_handles{iAx}.Color			= handles.color_unselected;
	handles.output(handles.ax_trial(iAx))	= false;
end
guidata(hObject, handles);

% --- Executes on button press in previous.
function previous_Callback(hObject, ~, handles)
% hObject    handle to previous (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.page > 1
	handles.page = handles.page-1;
	guidata(hObject, handles);
	initializePage(hObject, handles);
end


% --- Executes on button press in next.
function next_Callback(hObject, ~, handles)
% hObject    handle to next (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.page < handles.num_pages
	handles.page = handles.page+1;
	guidata(hObject, handles);
	initializePage(hObject, handles);
end


% --- Executes when number is entered into page field.
function pagesel_Callback(hObject, ~, handles)
% hObject    handle to pagesel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pagesel as text
%        str2double(get(hObject,'String')) returns contents of pagesel as a double
in = str2double(get(hObject,'String'));
if in < 1 || in > handles.num_pages || mod(in,1) ~= 0
	set(hObject,'String', num2str(handles.page));
	set(handles.pageinvalid, 'Visible', 'on');
	initializePage(hObject, handles);
else
	handles.page = in;
	guidata(hObject, handles);
	set(handles.pageinvalid, 'Visible', 'off');
	initializePage(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function pagesel_CreateFcn(hObject, ~, ~)
% hObject    handle to pagesel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
	set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in close.
function close_Callback(~, ~, ~)
% hObject    handle to close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(gcf); % this leads to a call to figure1_CloseRequestFcn, which closereq doesn't


% --- Executes when user attempts to close figure1.
% Only executed if X is pressed.
function figure1_CloseRequestFcn(hObject, ~, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% This is taken from the internet. If we resume here, the output function
% will still have access to the handles, otherwise its empty. Not really
% sure why otherwise we delete the object, but that code is never reached
% anyways..
if isequal(get(hObject,'waitstatus'),'waiting')
	uiresume(hObject);
	guidata(hObject,handles);
else
	if ~isfield(handles, 'outputfile')|| isempty(handles.outputfile)
		warning('Something went wrong, no outputfile provided.')
	else
		out = handles.output;
		disp(['Saving output to ' handles.outputfile])
		save(handles.outputfile, 'out');
	end
	delete(hObject); % closes the figure
end


% --- Executes on key press with focus on figure1 or any of its controls.
function figure1_WindowKeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
if strcmp(eventdata.Key, 'leftarrow')
	previous_Callback(hObject, [], handles);
elseif strcmp(eventdata.Key, 'rightarrow')
	next_Callback(hObject, [], handles);
end

