function h=plotWeights(xy, wf, sz, cmin, cmax)
% h=plotWeights(xy, wf, sz, cmin, cmax)
    if nargin < 5
        cmax=nan;
    end
    
    if nargin < 4
        cmin = nan;
    end
    if nargin < 3 || isnan(sz)
        sz = 200;
    end
    
    clrs = getColors(wf,cmin,cmax);

    hold on;
    axis off;
    axis equal;
    set(gcf,'color','w');
    for ii = 1:numel(wf)
        h(ii)=plot(xy(ii,1), xy(ii,2), 'Marker', '.', 'MarkerSize', sz, ...
            'Color', clrs(ii,:), 'LineStyle', 'none');
    end
end
