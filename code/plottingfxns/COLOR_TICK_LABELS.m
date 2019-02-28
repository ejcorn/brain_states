function COLOR_TICK_LABELS(x,y,numClusters,clusterColors)


if numClusters == 5
	ax = gca;
	if ~exist('clusterColors','var')
		clusterColors = {'AB484F','591A23', 'AA709F','527183','7E874B'};
	end
	clusterColors = hex2rgb(clusterColors);
	if x
		for K = 1:5
			ax.XTickLabel{K} = sprintf('\\color[rgb]{%f,%f,%f}%s', ...
        	clusterColors(K,:), ax.XTickLabel{K});
		end
	end
	if y
		for K = 1:5	
			ax.YTickLabel{K} = sprintf('\\color[rgb]{%f,%f,%f}%s', ...
        	clusterColors(K,:), ax.YTickLabel{K});
		end
	end
end