addpaths;
masterdir = fullfile(basedir,'results',name_root);
savedir = fullfile(masterdir,'analyses','choosing_k');
mkdir(savedir);

load([masterdir,'/clusterAssignments/repkmeansPartitions',distanceMethod,name_root,'.mat'],'nreps');

minNumClusters = 2; maxNumClusters = 11; clusterRange = minNumClusters:maxNumClusters; numK = 1 + maxNumClusters-minNumClusters;
k_offset = minNumClusters - 1;
zRandAllK = zeros(numK,2);
%
for K = 1:numK
	load(fullfile(masterdir,'clusterAssignments',['zrand_k',num2str(clusterRange(K)),name_root,'.mat']),'zRandK');
	zRandAllK(K,:) = zRandK(clusterRange(K),:);
end

barcolors = [53 183 121; 68 1 84] / 255;
f = figure;
errorbar(zRandAllK(:,1),zRandAllK(:,2),'Color',barcolors(2,:));
xticks(1:(maxNumClusters-k_offset)); xticklabels(minNumClusters:maxNumClusters)
ylabel('Z-Scored Rand Index'); xlabel('k','FontWeight','bold');
prettifyEJC;

f.PaperUnits = 'centimeters';
f.PaperSize = [8 4];
f.PaperPosition = [0 0 8 4];
saveas(f,fullfile(savedir,['zRandbyK',num2str(nreps),'reps.pdf']),'pdf');

f = figure;
cv = zRandAllK(:,2)./zRandAllK(:,1)
categorical(clusterRange)
zRandAllK

%col = [53 183 121]/255; %viridis green
bar(1:(maxNumClusters-k_offset),cv,'FaceAlpha',.5,'FaceColor',barcolors(1,:));
xticks(1:(maxNumClusters-k_offset));xticklabels(minNumClusters:maxNumClusters);
ylabel('CV of ZRI'); xlabel('k','FontWeight','bold');
prettifyEJC;

f.PaperUnits = 'centimeters';
f.PaperSize = [8 4];
f.PaperPosition = [0 0 8 4];
saveas(f,fullfile(savedir,['CVbyK',num2str(nreps),'reps',name_root,'.pdf']),'pdf');

disp('zrand done');