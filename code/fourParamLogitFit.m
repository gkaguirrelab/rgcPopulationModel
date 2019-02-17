function [fitParams,fVal] = fourParamLogitFit(x,y,x0)

if nargin==2
    x0 = [-2,4,1,3];
end

% a - Hill's slope; the steepness of the curve
% b - inflection point; the point at which the curvature changes sign
% c - maximum asymptote
% d - asymmetry factor (symmetric around inflection when d = 1)
fourParamLogitFunc = @(a,b,c,d) c./((1+(x./b).^a).^d);

myObj = @(p) sum((y - fourParamLogitFunc(p(1),p(2),p(3),p(4))).^2);
[fitParams,fVal] = fmincon(myObj,x0,[],[]);


end

