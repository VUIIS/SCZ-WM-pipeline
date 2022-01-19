function scz_entrypoint(varargin)
% This function serves as the entrypoint to the matlab part of the
% pipeline. Its purpose is to parse the command line arguments, then call
% the main function that actually does the work.

% If there are multiple fmris, store them as a cell when inputing file
% names to function (i.e. use {})

%% Just quit, if requested - needed for Singularity build
if numel(varargin)==1 && strcmp(varargin{1},'quit') && isdeployed
	disp('Exiting as requested')
	exit
end

%% Parse the inputs and parameters

% Matlab's input parser is very convenient. We will add all arguments as
% "optional", providing default values when appropriate.
P = inputParser;

addOptional(P,'t1_niigz','')
addOptional(P,'fmri_niigz','')

% When processing runs on XNAT, we generally have the project, subject,
% session, and scan labels from XNAT available in case we want them. Often
% the only need for these is to label the QA PDF.
addOptional(P,'xnat_project','UNK_PROJ');
addOptional(P,'xnat_subject','UNK_SUBJ');
addOptional(P,'xnat_session','UNK_SESS');


% where to store the outputs.
addOptional(P,'out_dir','/OUTPUTS');
addOptional(P,'in_dir','/INPUTS');

% Parse
parse(P,varargin{:});

% Display the command line parameters - very helpful for running on XNAT,
% as this will show up in the outlog.
disp(P.Results)


%% Run the actual pipeline
scz_main(P.Results);


%% Exit
% But only if we're running the compiled executable. If we're testing in a
% matlab session, exiting is not helpful.
if isdeployed
	exit
end

