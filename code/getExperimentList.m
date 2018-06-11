function experiments = getExperimentList()
% get experiment names

dataPath = getpref('mtlipglm', 'dataPath');

fl = dir(fullfile(dataPath, 'stim', '*_stim.mat'));

experiments = arrayfun(@(x) strrep(x.name, '_stim.mat', ''), fl, 'uni', false);