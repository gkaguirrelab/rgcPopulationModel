function [fitParams,fVal] = fourParamLogitFit(x,y,x0)

if nargin==2
    x0 = [-1,1,1,1];
end

lb = [-3 0 0.25 0.5];
ub = [ 0 40 4 6];

% a - Hill's slope; the steepness of the curve
% b - inflection point; the point at which the curvature changes sign
% c - maximum asymptote
% d - asymmetry factor (symmetric around inflection when d = 1)
fourParamLogitFunc = @(a,b,c,d) c./((1+(x./b).^a).^d);

options = optimoptions(@fmincon,...
    'Display','off');

myObj = @(p) sum((y - fourParamLogitFunc(p(1),p(2),p(3),p(4))).^2);
[fitParams,fVal] = fmincon(myObj,x0,[],[],[],[],lb,ub,[],options);


end

