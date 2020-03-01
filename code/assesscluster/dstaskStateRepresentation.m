% makes Fig S4

addpaths;
load(fullfile(basedir,['data/Demographics',name_root,'.mat']));
load(fullfile(datadir,['TimeSeriesIndicators',name_root,'.mat']));

scanlab = {'RestComb','nBackComb'};
minNumClusters = 2; maxNumClusters = 11; 
clusterRange = minNumClusters:maxNumClusters; 
numK = 1 + maxNumClusters-minNumClusters;
rTR = 120; 
scanInd = repelem([0 1]', rTR*nobs);	% make new scan label for 50-50 rest task composition
subjInd = repelem([1:nobs,1:nobs]',rTR); % make new subject label for 50-50 rest task composition

restAbsentStateCount = zeros(1,numK);
nbackAbsentStateCount = zeros(1,numK);
dwelltime = cell(numK,2);
for C = 1:numK
	numClusters = clusterRange(C);
	load([masterdir,'/analyses/choosing_k/rvn/DownSampleTaskCluster_k',num2str(numClusters),name_root,'.mat']);
	kClusterAssignments = partition;
	
	for S = 1:numel(scanlab)
		dwelltime{C,S} = zeros(nobs,numClusters);
		for N = 1:nobs
			subjStates = kClusterAssignments(scanInd == (S-1) & subjInd == N);
			for K = 1:numClusters
				dwelltime{C,S}(N,K) = sum(subjStates == K);
			end
		end
	end
	restAbsentStateCount(C) = 100* sum(sum(dwelltime{C,1} == 0,2) > 0) / nobs;	% percent of subjects with >0 missing states
	nbackAbsentStateCount(C) = 100* sum(sum(dwelltime{C,2} == 0,2) > 0) / nobs;
end

RNcolors = {'#005C9F','#FF8400'};
RNcolors = hex2rgb(RNcolors);
dsTaskStateRepresentation = [restAbsentStateCount;nbackAbsentStateCount]';

% Fig S4f
f = figure;
bar(dsTaskStateRepresentation,'FaceAlpha',.5);
xticks(1:2:numK);xticklabels(minNumClusters:2:maxNumClusters);
ylabel('% Subjects'); xlabel('k','FontWeight','bold'); title('Absent States');
legend({'Rest','n-back'},'Location','northwest')
prettifyEJC;
f.PaperUnits = 'centimeters';
f.PaperSize = [8 4];
f.PaperPosition = [0 0 8 4];

saveas(f,[masterdir,'/analyses/choosing_k/rvn/DownSampleTaskStateRepresentation',name_root,'.pdf']);

