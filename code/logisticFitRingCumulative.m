cardinalMeridianNames =  {'nasal' 'superior' 'temporal' 'inferior'};

totalRGC = cell.totalRGC([0 90 180 270],cardinalMeridianNames);
    
% a - Hill's slope; the steepness of the curve
% b - inflection point; the point at which the curvature changes sign
% c - maximum asymptote
% d - asymmetry factor (symmetric around inflection when d = 1)
fourParamLogitFunc = @(a,b,c,d,x) c./((1+(x./b).^a).^d);

x=0:0.01:15;
[ ~, refPointIdx ] = min(abs(x-10));

figure
hold on

lineColor = {'r','g','b','k'};

for mm = 1:length(cardinalMeridianNames)

thicknessCumSum = calcRingCumulative(x,totalRGC.density.fitDegSq.(cardinalMeridianNames{mm})(x)');

%thicknessCumSum = thicknessCumSum ./ thicknessCumSum(refPointIdx);
thicknessCumSum = thicknessCumSum ./ 6e5;

plot(x,thicknessCumSum,'-','Color',lineColor{mm})

% Fit a logistic function
myObj = @(p) sum((thicknessCumSum - fourParamLogitFunc(p(1),p(2),p(3),p(4),x)).^2);
[fitVal, fVal] = fmincon(myObj,[-2,4,1,3],[],[]);

plot(x,fourParamLogitFunc(fitVal(1),fitVal(2),fitVal(3),fitVal(4),x),'.','Color',lineColor{mm})
end


