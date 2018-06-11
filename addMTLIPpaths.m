addpath(fullfile(pwd, 'code'))
addpath(fullfile(pwd, 'code', 'dependencies'))
addpath(fullfile(pwd, 'analysis_scripts'))

filepath = fileparts(which(mfilename));

% this is where your data live
dataPath = fullfile(filepath, 'data'); % EDIT THIS LINE %%%%%%%%%%%%%%%%%%
if isempty(dataPath)
    error('dataPath must be specified');
end
% set matlab preference so that other functions can access this
setpref('mtlipglm', 'dataPath', dataPath)

% set neuroGLM path
neuroGLMpath = 'D:\Dropbox\MatlabCode\Repos\neuroGLM\'; % EDIT THIS LINE %%%%%%%%%%%%%%%%%%
if isempty(neuroGLMpath)
    error('neuroGLM - must be downloaded from https://github.com/jcbyts/neuroGLM');
end
addpath(neuroGLMpath)
if exist('neuroGLM') ~= 2
    error('neuroGLM path invalid');
end

% --- setup directory structure
if isdir(dataPath)
    fit_dir = fullfile(dataPath, 'main_fits');
    if ~isdir(fit_dir)
        mkdir(fit_dir)
    end
    
    fit_dir = fullfile(dataPath, 'lip_trunc_fits');
    if ~isdir(fit_dir)
        mkdir(fit_dir)
    end
end