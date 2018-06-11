function PTA = pulseSTASingle(pulses, rtrial, ptimes, nkt)
% PTA = pulseSTA(pulses, rtrial, ptimes, nkt)

nb = size(rtrial,2);

nTrial = size(rtrial,1);


if nkt>size(rtrial,2)
    PTA=nan(nkt, 1);
    return
end

tstarts=(0:nTrial-1)*nb;
t=bsxfun(@plus, ptimes(:), tstarts);
t=t(:);
p=reshape(pulses', [], 1);
n=numel(t);
X=sparse(t, ones(n,1), p, nb*nTrial, 1);
Xs=makeStimRows(X, nkt);

X0=ones(nb*nTrial,1);
Xs=[X0 Xs];

y=reshape(rtrial', [], 1);
PTA=(Xs'*Xs)\(Xs'*y);
PTA(1)=[];
PTA=flipud(PTA);