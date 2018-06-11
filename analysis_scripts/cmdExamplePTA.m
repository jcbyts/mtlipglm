dataPath = getpref('mtlipglm', 'dataPath');

experiments = getExperimentList();
%% Load up an experiment
exname=experiments{1}; % this session has example LIP units with strong pulse-triggered responses


neurons = getNeurons(exname, dataPath);
stim    = getStim(exname, dataPath); % load stimulus struct

disp(numel(neurons))

%% Pick a neuron and plot PSTH and PTA
kNeuron = 3;

spikeTimes  = neurons(kNeuron).spikeTimes;
motionOnset = [stim.timing(:).motionon] + [stim.timing(:).plxstart];

binSize = 1/60; % frame rate (in seconds)
window  = [-.5 2.5]; % 100 ms before motion onset to 1.5 seconds after (in seconds)

% get spike count aligned to motion onset
[spcnt, bins]  = binSpTimes(spikeTimes, motionOnset, window, binSize);

% smooth spike count with a boxcar filter
sm = 3; % size of boxcar (number of bins)
tmp = filter(boxcar(sm)/sm, 1, spcnt');
tmp = flipud(tmp);
tmp = filter(boxcar(sm)/sm, 1, tmp);
tmp = flipud(tmp);

spcnt = tmp';

% conver to spike rate
sprate = spcnt/binSize;

goodTrials  = stim.goodtrial & ~any(isnan(spcnt), 2);

% --- plot choice-sorted PSTH
figure(1); clf
set(gcf, 'DefaultAxesColorOrder', lines)
subplot(2,1,1)
plot(bins, mean(sprate(goodTrials & stim.targchosen==1, :))); hold on
plot(bins, mean(sprate(goodTrials & stim.targchosen==2, :)));
xlim([-.1 1.2])
xlabel('Time from Motion Onset')
ylabel('Spike Rate')
title('Choice Sorted PSTH', 'FontWeight', 'normal')
legend({'Choice 1', 'Choice 2'}, 'Location', 'BestOutside')
% --- Get Whitened PTA

% get residuals (subtract off the effect of motion onset)
residuals = bsxfun(@minus, sprate, nanmean(sprate));


% pulse values
pulses = sum(stim.pulses(goodTrials,:,:),3);
% pulse times (in bins relative to spcnt column 1)
pulseTimes = find(binSpTimes(stim.timing(find(stim.goodtrial,1)).pulses - stim.timing(find(stim.goodtrial,1)).motionon, 0, window, binSize));

nTimeLags = 60;
pta = pulseSTA(pulses, residuals(goodTrials,:), pulseTimes, nTimeLags);
ptaBins = (1:nTimeLags) * binSize;
subplot(2,1,2)
set(gcf, 'DefaultAxesColorOrder', hot(12))
plot(ptaBins, pta)

leg_str = regexp(sprintf('Pulse %d,', (1:7)'), ',', 'split');
legend(leg_str{1:7}, 'Location', 'BestOutside')
xlabel('Time from Pulse Onset')
ylabel('Spikes / sec / gabor')
title('PTA', 'FontWeight', 'normal')

%% How does the PTA code work?
% pulseSTA.m is optimized for speed and takes advantage of the fact that
% pulses always occur at the same time (relative to motion onset) to do
% some tricks. Here, I'll compute the PTA a more readable (slower) way. Each
% of these steps are comparable to the equation in the manuscript under
% Methods, Neural Analysis

% get spike times
spikeTimes = neurons(kNeuron).spikeTimes;

% Get the time that each pulse came on. Each trial in the stim struct has
% timing relative to trial start. Add the trial start in ephys (plexon)
% clock to convert to the same clock as the spikes
pulseTimes  = cell2mat(arrayfun(@(x) x.pulses(:)' + x.plxstart, stim.timing, 'UniformOutput', false));
pulseValues = sum(stim.pulses,3);

% index into completed trials
goodTrials = stim.goodtrial & ~any(isnan(pulseTimes), 2);

% Function for binning spikes
binfun = @(t) (t==0) + ceil(t/binSize);

% Find the time the recording stopped
plexonStarts = [stim.timing.plxstart];
plexonStarts(isnan(plexonStarts)) = [];
lastTimeBin = binfun(plexonStarts(end) + 10);

% --- Bin spike times
spikeTimesBin = binfun(spikeTimes); % convert time to bins

% create a vector of binned spike counts (sparse does this for you). Could
% also use histc, or something manual.
binnedSpikes  = sparse(spikeTimesBin, ones(numel(spikeTimesBin),1), ones(numel(spikeTimesBin),1), lastTimeBin, 1);
binnedSpikes  = smooth(full(binnedSpikes), 5); % smooth the binned spikes
binnedSpikes  = binnedSpikes/binSize; % convert to spike rate

% --- Visualize the basic idea:
% cross-correlate the pulse values with the spike rate, but do it for each
% pulse separately, while accounting for their covariance
figure(3); clf
ax = gca();
bar(binnedSpikes, 'FaceColor', .5*[1 1 1])
ax2 = axes('Position', ax.Position);
stem(binfun(pulseTimes(:)), pulseValues(:), 'r', 'MarkerFaceColor', 'r');
xlim(ax, [0 2e3])
xlim(ax2, [0 2e3])
% legend({'pulses', 'spike count'})
xlabel(ax,'Time Bin')
ylabel(ax2,'Pulse Value');
ylabel(ax,'Spike Rate');
set(ax2, 'YAxisLocation', 'left', 'color', 'none', 'box', 'off', 'XTick', [], 'YColor', 'r')
set(ax, 'YAxisLocation', 'right', 'color', 'none', 'box', 'off')

%% --- setup PTA analysis
nTimeLags = 60; % how many bins after a pulse to analyze

% "Unwrapping the convolution"
% The PTA is a linear model that the spike rate of the neuron is equal to
% the convolution of the PTA Kernels with the pulse impulses. Here, we
% create two copies of the identity basis that are sized for the number of
% time lags we're interested in recovering. By convolving these lagged bins
% with the pulse impulses, we can learn the PTA kernels using linear
% regression.

% Hear are the two identity matrices that we use to unwrap the convolution
filterPulse       = diag(ones(nTimeLags, 1));
filterMotionOnset = diag(ones(2*nTimeLags, 1)); % length of motion + pta

% Next we build the X matrix from the manuscript. We convolve each pulse
% onset with the filter from above, which unwraps each time-lag as a column
% of a matrix. We then concatenate all those matrices so that each column
% holds the values of a certain pulse at a certain time-lag


% First, the DCterm for motion onset effects. In pulseSTA, we remove this
% by analyzing the residuals. That is also how we analyze the
% choice-corrected PTA (by running the PTA analysis on the residual spike
% counts after subtracting off the average for each choice accordingly)

% get motion onset times
motionOnsets = binfun(pulseTimes(goodTrials,1));
n = numel(motionOnsets);

% create a vector of delta functions at the times of motion onset
tmp = sparse(motionOnsets, ones(n,1), ones(n,1), lastTimeBin,1);

% create a matrix where each column is a time-lag relative to motion onset
tmp = conv2(full(tmp), filterMotionOnset);
tmp = tmp(1:lastTimeBin,:);
    
% store that matrix in our deisgn matrix, X
Xdesign = tmp;

% loop over pulses and do the same thing for each pulse
for kPulse = 1:7
    
    % get pulse times
    binnedPulseOnsets = binfun(pulseTimes(goodTrials,kPulse));
    n = numel(binnedPulseOnsets);
    % create a vector of the pulse values (at the time the pulse came on)
    tmp = sparse(binnedPulseOnsets, ones(n,1), pulseValues(goodTrials,kPulse), lastTimeBin,1);
    % create a matrix where each column is those values at a particular
    % time lag
    tmp = conv2(full(tmp), filterPulse);
    tmp = tmp(1:lastTimeBin,:);
    
    % concatenate to keep building up the full design matrix
    Xdesign = [Xdesign tmp];
end

% project the design matrix on the binned spike rate (The raw pulse
% triggered response)
XY = Xdesign'*binnedSpikes;

% get the sample covariance matrix for the pulses (used to whiten)
XX = Xdesign'*Xdesign;

% Least squares equation (whitened PTA)
wls = XX\XY;

DC  = wls(1:(2*nTimeLags));
PTA = reshape(wls((2*nTimeLags+1):end), [], 7);

% --- Plot it
figure(2); clf
subplot(121)
set(gcf, 'DefaultAxesColorOrder', lines)
plot((1:(2*nTimeLags))*binSize,DC)
xlabel('Time from motion onset')
ylabel('Spike Rate')
title('Gabor onset response', 'fontweight', 'normal')

subplot(122)
set(gcf, 'DefaultAxesColorOrder', hot(12))
plot( (1:nTimeLags)*binSize, PTA)
xlabel('Time from Pulse Onset')
ylabel('Spikes / Sec / Gabor')
title('Pulse-Triggered Average', 'fontweight', 'normal')

