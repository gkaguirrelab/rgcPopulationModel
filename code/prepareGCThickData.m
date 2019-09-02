function data = prepareGCThickData()

% Script to take the GC thickness profile generated by retinaTOMEAnalysis
% and turn it into a file suitable for analysis by this repo.

rawDataFile = '/Users/aguirre/Dropbox (Aguirre-Brainard Lab)/AOSO_analysis/OCTExplorerExtendedHorizontalData/LineAnalysisResults.mat';
load(rawDataFile,'meanGCVec','XPos_Degs');

data(1).label = 'temporal';
data(1).angle = 180;
data(1).supportDeg = abs(XPos_Degs(XPos_Degs<=0));
data(1).thickMM = (meanGCVec(XPos_Degs<=0)'./1000)';

data(2).label = 'nasal';
data(2).angle = 0;
data(2).supportDeg = abs(XPos_Degs(XPos_Degs>=0));
data(2).thickMM = (meanGCVec(XPos_Degs>=0)'./1000)';

end