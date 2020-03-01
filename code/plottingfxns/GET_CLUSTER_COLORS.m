function clusterColors = GET_CLUSTER_COLORS(k)
	% return cluster color hex codes for plotting for a given number of clusters
	% only works right now for 5 and 6, otherwise returns black

if k == 5
    clusterColors = {'AB484F','591A23', 'AA709F','527183','7E874B'};
elseif k == 6
    clusterColors = {'AB484F','591A23', 'D98FCA','945F8A','C6CDF7','7294D4'};
else
	clusterColors = repmat('000000',[1 k]);
end