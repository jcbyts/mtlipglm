% This script executes plotting code for figure 3 of Yates et al. (2017)

% path to the data base directory. dataPath should have subdirectories
% neurons, stim, main_fits, lip_trunc_fits
dataPath = getpref('mtlipglm', 'dataPath');

%% Load up LIP fits

% The raw data have been fit with the required models and analyzed. To
% regenerate fits, use cmdFitMT.m
% Note: The fits in cmdFitMT.m are generated with a fixed hyperparameter 
% so the exact numbers may subtly differ from the published text

S = loadUpFits(dataPath, 'MT');

%% Plot Model PSTH on top of data PSTH
% The three models compared in figure 5 are the "Stimulus-to-LIP" model,
% the "MT-to-LIP" model and "MT-to-LIP (full choice)"

modelNames=arrayfun(@(x) x.name, S(1).model, 'UniformOutput', false);

% The three models have different names in the code:
% Stimulus-to-LIP = Poisson
% MT-to-LIP       = MTsim
% MT-to-LIP (full choice) = MTsimChoice
% Compare them to the data in each plot
modelIxs={[find(strcmp(modelNames, 'data')) find(strcmp(modelNames, 'Poisson'))]};


for iModel=1:numel(modelIxs)
    figure(iModel); clf
    
    modelIx  = modelIxs{iModel};
    nModels  = numel(modelIx);
    nNeurons = numel(S);
    
    plotFields={'psthCoh', 'ptaRaw'};
    nPlotFields = numel(plotFields);
    
    for f=1:nPlotFields
        subplot(1,nPlotFields, f)
        
        field=plotFields{f};
        
        
        % --- Change Color Map
        if any(strfind(field, 'pta'))
            cmap=hot(12);
        else
            cmap=cbrewer('jake', 'rdbu', 8);
        end
        
        for kM=1:nModels
            
            kModel=modelIx(kM);
            
            % get normalized population response
            cohpData = reshape(cell2mat(arrayfun(@(x) x.model(kModel).(field)/max(x.model(1).(field)(:)), S, 'UniformOutput', false)), [size(S(1).model(kModel).(field)) nNeurons]);
            popcoh   = nanmean(cohpData,3);
            
            % set time axis
            if any(strfind(field, 'pta'))
                timex=S(1).model(1).ptaTime;
                xdim = [0 2000];
            else
                timex=S(1).model(1).psthTime(:);
                xdim = [-500 1500];
            end
            
            % change marker for data and model
            if kModel==1
                marker='.';
            else
                marker='-';
            end
            
            % plot
            for kP=1:size(popcoh,2)
                if size(timex,2)>1
                    txx=timex(:,kP);
                else
                    txx=timex;
                end
                
                plot(txx, popcoh(:,kP), marker, 'Color', cmap(kP,:), 'MarkerSize', 2); hold on
            end
            
            % title is model name
            title(S(1).model(kModel).name)
            xlim(xdim)
            
        end
        
    end
    
end
