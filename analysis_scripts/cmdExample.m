dataPath = getpref('mtlipglm', 'dataPath');

% get list of experiments to choose from
experimentList = getExperimentList();
disp(experimentList(:))


%% Load up an experiment
exname=experimentList{1};

mstruct=mtlipglm(exname);
mstruct.buildTrialStruct('IncludeContrast', true, 'MotionEpoch', false, 'MTsim', false);

%% plot sample trials
figure(1); clf
subplot(212)
mstruct.plotTrial

subplot(211)
mstruct.neurons(1).plotTaskGeometry

%% Plot choice sorted PSTH for each neuron
figure(3); clf
sp = ceil(sqrt(numel(mstruct.neurons)));

for kNeuron = 1:numel(mstruct.neurons)
    subplot(sp,sp,kNeuron)
    neuronName = mstruct.neurons(kNeuron).getName(0);
    mstruct.plotChoPSTH(neuronName, 'motionon', 'smoothing', 100, 'window', [-100 1.3e3])
    title(sprintf('%d: %s', kNeuron, neuronName))
    xlim([-100 1.2e3])
end


%% Plot PTA for each neuron
figure(3); clf
sp = ceil(sqrt(numel(mstruct.neurons)));

for kNeuron = 1:numel(mstruct.neurons)
    subplot(sp,sp,kNeuron)
    neuronName = mstruct.neurons(kNeuron).getName(0);
    mstruct.plotPTA(neuronName)
    title(sprintf('%d: %s', kNeuron, neuronName))
%     xlim([-100 1.2e3])
end


%% Try fitting an example neuron
% This script demonstrates how to use neuroGLM to fit an example neuron
% with different covariates

kNeuron = 3;  % neuron id
binSize = 10;  % in ms

% --- initialize model
% This sets up a neuroGLM object.
g=glmspike('example');
g.setDesignOptions('binSize', binSize)
binfun = g.binfun;
g.fitNeuron=mstruct.neurons(kNeuron).getName(0);

%% Explore basis function options
% basisFactory can create a number of temporal basis fcuntions
help basisFactory
nBasisFunctions = 10;
figure(4); clf

subplot(1,3,1)
basisType = 'nonlinearly scaled cosine';
bs=basisFactory.makeSmoothTemporalBasis(basisType, [100 500], nBasisFunctions, g.binfun, 1500);
bs.normalize
bs.plot
xlabel('Time (bins)')
title(basisType)

subplot(1,3,2)
basisType = 'raised cosine';
bs=basisFactory.makeSmoothTemporalBasis(basisType, 1500, nBasisFunctions, g.binfun);
bs.normalize
bs.plot
xlabel('Time (bins)')
title(basisType)

subplot(1,3,3)
basisType = 'boxcar';
bs=basisFactory.makeSmoothTemporalBasis(basisType, 1500, nBasisFunctions, g.binfun);
bs.plot
xlabel('Time (bins)')
title(basisType)

%% Add covariates

nBasisFunctions = 15;
% --- Add motion stimulus 
basisType = 'nonlinearly scaled cosine';
bs=basisFactory.makeSmoothTemporalBasis(basisType, [100 1000], nBasisFunctions, g.binfun, 100);
bs.normalize

% add motion onset/offset as boxcar
% g.addCovariateBoxcar(mstruct.trial, 'motion', 'motionon', 'motionoff', 'contrast', bs);

% or add it as an impulse
g.addCovariateTiming(mstruct.trial, 'motion', 'motionon', 'contrast', bs, 0);


% add pulses as boxcar
basisType = 'nonlinearly scaled cosine';
bs=basisFactory.makeSmoothTemporalBasis(basisType, [100 1000], nBasisFunctions, g.binfun, 1e3);
bs.normalize

help neuroGLM/addCovariateBoxcar
g.addCovariateBoxcar(mstruct.trial, 'direction', 'pulseon', 'pulseoff', 'pulses', 'Motion Pulses', bs)

% % or add pulses as impluses
% stimHandle=@(trial) sparse(binfun(trial.pulseon), 1, trial.pulses, binfun(trial.duration),1);
% g.addCovariate(trial, 'direction', 'direction', stimHandle, bs)

% add choice terms
covariateDuration=3e3;
bs=basisFactory.makeSmoothTemporalBasis('raised cosine', covariateDuration, ceil(covariateDuration/100), g.binfun);
bs.normalize

% acausal, aligned to saccade
offset=-2500;
g.addCovariateTiming(mstruct.trial, 'saccade1', 'saccade', 'Sac T1', bs, offset, @(trial) (trial.choice == 1))
g.addCovariateTiming(mstruct.trial, 'saccade2', 'saccade', 'Sac T2', bs, offset, @(trial) (trial.choice == 0))

% % causal, aligned to motionon
% offset=0;
% g.addCovariateTiming(mstruct.trial, 'saccade1', 'motionon', 'Sac T1', bs, offset, @(trial) (trial.choice == 1))
% g.addCovariateTiming(mstruct.trial, 'saccade2', 'motionon', 'Sac T2', bs, offset, @(trial) (trial.choice == 0))

% add targets
bs=basisFactory.makeSmoothTemporalBasis('raised cosine', 3000, 30, g.binfun);
bs.normalize
g.addCovariateTiming(mstruct.trial, 'Targets', 'targson', 'Targets', bs)

%% Use recorded spikes as covariates
% % add history filter

% bs=basisFactory.basisFactory('history', g.binfun, 200);
% bs.plot
% g.addCovariateSpiketrain(mstruct.trial, 'history', g.fitNeuron, 'Post Spike Filter', bs)

% % add coupling filters (causal)
% neuronNames = arrayfun(@(x) x.getName(0), mstruct.neurons, 'UniformOutput', false);
% coupledNeurons = setdiff(neuronNames, g.fitNeuron);
% for kCouple = 1:numel(coupledNeurons)
%     g.addCovariateSpiketrain(mstruct.trial, coupledNeurons{kCouple}, coupledNeurons{kCouple}, ['Coupling: '  coupledNeurons{kCouple}], bs)
% end

%% compile design matrix
% pick trials to use for fitting
% convolve added covariates with their respective basis functions to build
% the design matrix
idx=1:numel(mstruct.trial); % trials to use

g.compileDesignMatrix(mstruct.trial, idx)

%% Fit and evaluate model
% try ridge regression
rho = 1; % ridge parameter (play around with this to see how it effects fitting. Use cross validation on the training set to set it properly)

g.addBiasColumn('right'); % augment with column of ones

% fit with k-fold cross validation
kFolds = 5;
g.fitCV(mstruct.trial, kFolds, rho)

% Evaluate model performance

S = mstruct.modelComparison(g);


%% plot weights

% Visualize
ws = g.combineWeights(mean([g.modelfit.khat],2));
wvar = g.combineWeights(mean([g.modelfit.SDebars],2));

fig = figure(2913); clf;
subplot(2,3,1)

field='motion';
if isfield(ws, field)
    plot(ws.(field).tr, ws.(field).data, 'k'); hold on
    plot(wvar.(field).tr, ws.(field).data+wvar.(field).data, 'k--')
    plot(wvar.(field).tr, ws.(field).data-wvar.(field).data, 'k--')
    axis tight
end
title(field)

field='contrast';
if isfield(ws, field)
    plot(ws.(field).tr, ws.(field).data, 'k'); hold on
    plot(wvar.(field).tr, ws.(field).data+wvar.(field).data, 'k--')
    plot(wvar.(field).tr, ws.(field).data-wvar.(field).data, 'k--')
    axis tight
end
title(field)

subplot(2,3,2)
field='direction';
if isfield(ws, field)
    plot(ws.(field).tr, ws.(field).data, 'k'); hold on
    plot(wvar.(field).tr, ws.(field).data+wvar.(field).data, 'k--')
    plot(wvar.(field).tr, ws.(field).data-wvar.(field).data, 'k--')
    axis tight
else
    axis off
end
title(field)

subplot(2,3,3)
field='history';
if isfield(ws, field)
    plot(ws.(field).tr, ws.(field).data, 'k'); hold on
    plot(wvar.(field).tr, ws.(field).data+wvar.(field).data, 'k--')
    plot(wvar.(field).tr, ws.(field).data-wvar.(field).data, 'k--')
    axis tight
else
    axis off
end
title(field)

subplot(2,3,4)
field='choice1';
if isfield(ws, field)
plot(ws.(field).tr, ws.(field).data, 'r'); hold on
field='choice2';
plot(ws.(field).tr, ws.(field).data, 'b')
axis tight
else
    axis off
end
title('Choice')

subplot(2,3,5)
field='saccade1';
if isfield(ws, field)
plot(ws.(field).tr, ws.(field).data, 'r'); hold on
plot(wvar.(field).tr, ws.(field).data+wvar.(field).data, 'r--')
plot(wvar.(field).tr, ws.(field).data-wvar.(field).data, 'r--')

field='saccade2';
plot(ws.(field).tr, ws.(field).data, 'b')
plot(wvar.(field).tr, ws.(field).data+wvar.(field).data, 'b--')
plot(wvar.(field).tr, ws.(field).data-wvar.(field).data, 'b--')

axis tight
else
    axis off
end
title('Saccade')

subplot(2,3,6)
field='Targets';
if isfield(ws, field)
    plot(ws.(field).tr, ws.(field).data, 'k'); hold on
    plot(wvar.(field).tr, ws.(field).data+wvar.(field).data, 'k--')
    plot(wvar.(field).tr, ws.(field).data-wvar.(field).data, 'k--')
    axis tight
else
    axis off
end
title(field)


%% plot psth

figure(10); clf
set(gcf, 'DefaultAxesColorOrder', jet(size(S.model(1).psthCoh,2)))
subplot(121)
plot(S.model(1).psthTime, S.model(1).psthCoh)
yd = ylim;
xlabel('Time (ms)')
ylabel('Firing rate')
title('Data')
subplot(122)
plot(S.model(2).psthTime, S.model(2).psthCoh)
ylim(yd)
xlabel('Time (ms)')
ylabel('Firing rate')
title('GLM')
%% plot PTA
figure(11); clf
set(gcf, 'DefaultAxesColorOrder', hot(12))
subplot(121)
plot(S.model(1).ptaTime, S.model(1).ptaRaw)
xlabel('Time (ms)')
ylabel('Firing rate')
axis tight
title('Data')
subplot(122)
plot(S.model(2).ptaTime, S.model(2).ptaRaw)
xlabel('Time (ms)')
ylabel('Firing rate')
title('GLM')
axis tight


