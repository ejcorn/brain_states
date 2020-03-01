function COLOR_TICK_LABELS(x,y,numClusters,clusterColors)


if numClusters == 5 | numClusters == 6
	ax = gca;
	if ~exist('clusterColors','var')
		clusterColors = GET_CLUSTER_COLORS(numClusters);
	end
	clusterColors = hex2rgb(clusterColors);
	if x
		for K = 1:numClusters
			ax.XTickLabel{K} = sprintf('\\color[rgb]{%f,%f,%f}%s', ...
        	clusterColors(K,:), ax.XTickLabel{K});
		end
	end
	if y
		for K = 1:numClusters	
			ax.YTickLabel{K} = sprintf('\\color[rgb]{%f,%f,%f}%s', ...
        	clusterColors(K,:), ax.YTickLabel{K});
		end
	end
end