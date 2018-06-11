function stim = getStim(exname, dataPath)
% get stimulus data for experiment 'exname', in folder <dataPath>/stim
% stim = getStim(exname, dataPath)
if nargin < 2
    dataPath = getpref('mtlipglm', 'dataPath');
end
stim = load(fullfile(dataPath, 'stim', [exname '_stim.mat']));