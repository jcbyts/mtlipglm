function nw = getAllNeurons(dataPath)
% nw = getAllNeurons(tags)
    
experiments = getExperimentList();
nExperiments = numel(experiments);

for kEx = 1:nExperiments
    exname = experiments{kEx};
    if kEx == 1
        nw = getNeurons(exname, dataPath);
        continue
    end
    
    nw = [nw getNeurons(exname, dataPath)];
end

if isempty(nw)
    error('No neurons loaded! Make sure data files are in [%s]', dataPath);
end