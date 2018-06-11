function PTA = pulseSTA(pulses, rtrial, ptimes, nkt, addDC)
% PULSE STA computes the whitened pulse-triggered average
% Inputs:
%   pulses [nTrials x nPulses] - signed value of pulses on each trial
%   rtrial [nTrials x nBins]   - binned spike rate (or count)
%   ptimes [1 x nPulses]       - bin that each pulse started at
%   nkt    [1 x 1]             - number of time lags to include
%   addDC  [boolean]           - include DC term (default: false)
% Output:
%   PTA [nkt x nPulses]        - the whitened Pulse-Triggered Average

assert(~any(isnan(pulses(:))), 'There are NaNs in the variable: pulses')
assert(~any(isnan(rtrial(:))), 'There are NaNs in the variable: rtrial')

nTimeBins = size(rtrial,2);
nPulses   = size(pulses,2);
nTrials   = size(rtrial,1);

assert(nkt<nTimeBins, 'Kernel Size must be smaller than the window of spikes')

% Two ways to compute it: the fast way and the readable way
if ~exist('quickFlag', 'var')
    quickFlag = true; % default to using the quick way
end

if ~exist('addDC', 'var')
    addDC = false;
end

% the fast way: pre compute the covariance matrix of the pulse stimulus and
% the covariance matrix of the time lagged pulse onsets. Combine these at
% the end. This takes advantage of the fact that the pulse onset times are
% always the same (relative to motion onset)
filtk=diag(ones(nkt,1));
kFilt=[zeros(ptimes(1)-1, nkt); filtk; zeros(nTimeBins-nkt-ptimes(1)+1, nkt)];

if addDC
    Xs=kFilt;
    XY=sum(rtrial*kFilt)';
else
    Xs = [];
    XY = [];
end

for kPulse=1:nPulses
    kFilt=[zeros(ptimes(kPulse)-1, nkt); filtk; zeros(nTimeBins-nkt-ptimes(kPulse)+1, nkt)];
    Xs=[Xs kFilt];
    XY=[XY; sum(bsxfun(@times, (rtrial*kFilt), pulses(:,kPulse)))'];
end

% Loop over trials and build up pulse covariance matrix
YY=0;
XX=0;
for kTrial=1:nTrials
    if addDC
        thisTrialCov = [1 pulses(kTrial,:)]'*[1 pulses(kTrial,:)];
    else
        thisTrialCov = pulses(kTrial,:)'*pulses(kTrial,:);
    end
    
    XX=XX+thisTrialCov;
    YY = YY + rtrial(kTrial,:)*rtrial(kTrial,:)';
end

% rescale pulse covariance matrix to capture each time bin
XX=(Xs'*Xs).*imresize(XX, nkt, 'nearest');

nY = numel(rtrial);

wls = XX\XY; % linear regression

if addDC
    PTA = (reshape(wls, [], nPulses+1));
    PTA(:,1) = []; % remove bias
else
    PTA = reshape(wls, [], nPulses);
end