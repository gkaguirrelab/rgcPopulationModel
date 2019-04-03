function cellProportionModel

cardinalMeridianAngles = [0 90 180 270];
cardinalMeridianNames =  {'nasal' 'superior' 'temporal' 'inferior'};
lineColor = {'r','g','b','c'};

% The total number of retinal ganglion cells in the retina is ~1 million.
totalRetinalRGCs = 1e6;

% Regular support for cumulative sum
regularSupportCumSum = 0:0.01:1.2;

% Create a function handle that provides spline interpolated values
splineInterp = @(x,y,xq) ppval(spline(x,y),xq);

% Obtain the average cumulative RGC density across meridians across
% regualar support
regularSupportDeg=0:0.01:90;
loopVar = zeros(size(regularSupportDeg));
for mm = 1:length(cardinalMeridianNames)
    fitRGCDensitySqDegVisual = getSplineFitToRGCDensitySqDegVisual(cardinalMeridianAngles(mm));
    loopVar = loopVar+fitRGCDensitySqDegVisual(regularSupportDeg);
end
meanRGCDensity = loopVar./length(cardinalMeridianNames);
meanRGCCumSum = calcRingCumulative(regularSupportDeg,meanRGCDensity);
maxMeanRGCCumSum = max(meanRGCCumSum);
meanRGCCumSum = meanRGCCumSum./maxMeanRGCCumSum;

% Generate a function that expresses amacrine cell count across regular
% support and with a given profile of RGC density. The amacrine values were
% reported for the average across merdians, and thus will be related to the
% average RGCCumSum
amacrine = cell.amacrine( 0, {'nasal'} );
amacrineCumSum = calcRingCumulative(regularSupportDeg,amacrine.density.fitDegSq.nasal(regularSupportDeg)');
amacrineCumSum = amacrineCumSum./maxMeanRGCCumSum;

figure
plot(regularSupportDeg,meanRGCDensity,'-k');
hold on
plot(regularSupportDeg,amacrine.density.fitDegSq.nasal(regularSupportDeg)','-r');


% Remove the initial run of points for which the meanRGCCumSum is all zero
initialZeroIdx = find(meanRGCCumSum==0);
meanRGCCumSum=meanRGCCumSum(initialZeroIdx(end):end);
amacrineCumSum=amacrineCumSum(initialZeroIdx(end):end);
refucedSupportDeg = regularSupportDeg(initialZeroIdx(end):end);


figure
plot(refucedSupportDeg,meanRGCCumSum,'-k')    
hold on
plot(refucedSupportDeg,amacrineCumSum,'-r')    
hold off

figure
plot(meanRGCCumSum,amacrineCumSum)
plot(amacrineCumSum./meanRGCCumSum)
hold on


% And now a function to return that provides the thickness of the RGC
% portion of the RGC+IPL layer given that total thickness over regular
% support.
amacrineDensityFunc = @(regularSupportDeg, rgcDensity) splineInterp(regularSupportCumSum,splineInterp(meanRGCCumSum,amacrineDensity',regularSupportCumSum),relativeCumulativeCount(regularSupportDeg, rgcDensity, totalRetinalRGCs));
for mm = 2:length(cardinalMeridianNames)
    fitRGCDensitySqDegVisual = getSplineFitToRGCDensitySqDegVisual(cardinalMeridianAngles(mm));
    plot(regularSupportDeg,amacrineDensityFunc(regularSupportDeg, fitRGCDensitySqDegVisual(regularSupportDeg)),'-','Color',lineColor{mm})
end


%     plot(x,cellCountCumSum,'-','Color',lineColor{mm})
%     [fitParams,fVal] = fourParamLogitFit(x,cellCountCumSum);
%     fitParams
%     plot(x2,fourParamLogitFunc(fitParams(1),fitParams(2),fitParams(3),fitParams(4),x2),'.','Color',lineColor{mm})    

end

